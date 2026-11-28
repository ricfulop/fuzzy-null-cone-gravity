"""Export bounded L=2 angular weak-slot diagnostic snapshots.

This is a fixed-radius angular diagnostic for the prepared oblate weak forms.
It exports complex equal-m angular averages of the four radial derivative slots
for the potential and lapse weak densities, plus the separate radial Green-flux
slots.  It intentionally performs no radial Galerkin projection, finite K/E/D
assembly, rank test, endpoint-domain choice, seed search, or spectrum.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
import time
from pathlib import Path
from typing import Any, Callable

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import kinematic_blocks as kinematics  # noqa: E402
import weak_forms as weak  # noqa: E402


L_CUTOFF = 2
M_VALUES = (0, 1, 2)
R_VALUES = (sp.Rational(3, 2), sp.Rational(2, 1))
EPSILON_VALUE = sp.Rational(1, 4)
QUADRATURE_ORDERS = (24, 48, 96)
REFINEMENT_RTOL = 2e-10
CONTROL_BATCHES = 3
ANGULAR_CASES = len(M_VALUES) * len(R_VALUES)
SLOTS = ("00", "01", "10", "11")
ZERO_TOL = 1e-11
RESULT_NAME = "angular-weak-slot-checks.json"
PROGRESS_NAME = "angular-weak-slot-progress.jsonl"
CHECKPOINT_NAME = "angular-weak-slot-checkpoint.json"
DENOMINATOR_CACHE: set[tuple[sp.Rational, ...]] = set()


class RegularityError(ValueError):
    """Raised when a fixed-chart density is not a supported rational function."""


def now() -> float:
    return time.monotonic()


def label_tuple(label: tuple[Any, ...]) -> str:
    return "_".join(str(part) for part in label)


def real_labels(labels: list[tuple[Any, ...]], m: int) -> list[str]:
    base = [label_tuple(label) for label in labels]
    if m == 0:
        return base
    return [f"{label}_cos" for label in base] + [f"{label}_sin" for label in base]


def cscalar(value: complex | np.complexfloating) -> list[float]:
    z = complex(value)
    return [float(z.real), float(z.imag)]


def matrix_to_json(matrix: np.ndarray) -> list[list[list[float]]]:
    return [[cscalar(value) for value in row] for row in matrix.tolist()]


def real_matrix_to_json(matrix: np.ndarray) -> list[list[float]]:
    return [[float(value) for value in row] for row in matrix.tolist()]


def emit(payload: dict[str, Any]) -> None:
    print("::orx:: " + json.dumps(payload, sort_keys=True), flush=True)


def stage(name: str, **payload: Any) -> None:
    print(f"STAGE {name} {json.dumps(payload, sort_keys=True)}", flush=True)


def substage(name: str, event: str, start: float | None = None, **payload: Any) -> float:
    elapsed = None if start is None else now() - start
    print(
        f"SUBSTAGE {name} {event} "
        + json.dumps({"elapsed_seconds": elapsed, **payload}, sort_keys=True),
        flush=True,
    )
    return now()


def source_sha256() -> dict[str, str]:
    paths = [
        Path(__file__).resolve(),
        HERE / "weak_forms.py",
        HERE / "kinematic_blocks.py",
        HERE / "potential_hessian.py",
    ]
    return {
        str(path.relative_to(HERE)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in paths
    }


def safe_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(payload, allow_nan=False, indent=2, sort_keys=True) + "\n")
    tmp.replace(path)


def write_progress(output_dir: Path | None, payload: dict[str, Any]) -> None:
    emit(payload)
    if output_dir is not None:
        with (output_dir / PROGRESS_NAME).open("a") as handle:
            handle.write(json.dumps(payload, sort_keys=True) + "\n")


def write_checkpoint(
    output_dir: Path | None,
    started: float,
    stage_name: str,
    payload: dict[str, Any],
) -> None:
    if output_dir is None:
        return
    checkpoint = {
        "status": "running",
        "stage": stage_name,
        "elapsed_seconds": now() - started,
        "payload": payload,
    }
    safe_write_json(output_dir / CHECKPOINT_NAME, checkpoint)


def runtime_plan() -> list[dict[str, Any]]:
    return [
        {"phase": "setup", "of": 1, "unit": "environment"},
        {"phase": "controls", "of": CONTROL_BATCHES, "unit": "batches"},
        {"phase": "angular_export", "of": ANGULAR_CASES, "unit": "snapshots"},
        {"phase": "export", "of": 1, "unit": "checker result"},
    ]


def clean_expr(value: sp.Expr, x: sp.Symbol, cache: dict[str, sp.Expr]) -> sp.Expr:
    key = sp.srepr(value)
    if key in cache:
        return cache[key]
    if value.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise RegularityError("nonexact or nonfinite assembled chart expression")
    try:
        reduced = sp.cancel(sp.expand(value))
    except Exception as exc:
        raise RegularityError("could not cancel assembled chart expression") from exc
    if reduced.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise RegularityError("nonexact or nonfinite cancelled chart expression")
    for power in reduced.atoms(sp.Pow):
        if power.exp.is_Rational and not power.exp.is_Integer:
            raise RegularityError(f"unsupported radical after same-m cancellation: {power}")
    if reduced.has(sp.Abs):
        raise RegularityError(f"unsupported Abs after same-m cancellation: {reduced}")
    numerator, denominator = sp.fraction(reduced)
    try:
        sp.Poly(numerator, x, domain=sp.QQ)
    except Exception as exc:
        raise RegularityError("numerator is not a rational polynomial in the chart") from exc
    denominator = sp.expand(denominator)
    if denominator == 0:
        raise RegularityError("zero denominator after same-m cancellation")
    _check_denominator_roots(denominator, x)
    cache[key] = reduced
    return reduced


def _check_denominator_roots(denominator: sp.Expr, x: sp.Symbol) -> None:
    try:
        poly = sp.Poly(denominator, x, domain=sp.QQ)
    except Exception as exc:  # noqa: BLE001 - turn SymPy coercions into typed chart failures.
        raise RegularityError(
            f"denominator is not a rational polynomial in the chart: {denominator}"
        ) from exc
    if poly.is_zero:
        raise RegularityError("zero denominator polynomial")
    if poly.eval(-1) == 0 or poly.eval(1) == 0:
        raise RegularityError(f"denominator endpoint root on -1<=x<=1: {poly.as_expr()}")
    monic = poly.monic()
    cache_key = tuple(sp.Rational(coefficient) for coefficient in monic.all_coeffs())
    if cache_key in DENOMINATOR_CACHE:
        return
    try:
        interior_roots = monic.count_roots(sp.Rational(-1), sp.Rational(1))
    except Exception as exc:  # noqa: BLE001 - keep Sturm failures typed as regularity failures.
        raise RegularityError(
            f"could not certify denominator roots on -1<x<1: {monic.as_expr()}"
        ) from exc
    if interior_roots != 0:
        raise RegularityError(f"denominator interior root on -1<x<1: {monic.as_expr()}")
    DENOMINATOR_CACHE.add(cache_key)


def split_complex(value: sp.Expr, x: sp.Symbol, cache: dict[str, sp.Expr]) -> tuple[sp.Expr, sp.Expr]:
    expanded = sp.expand(value)
    real = clean_expr(expanded.coeff(sp.I, 0), x, cache)
    imag = clean_expr(expanded.coeff(sp.I, 1), x, cache)
    residual = sp.expand(expanded - real - sp.I * imag)
    if residual != 0:
        residual = clean_expr(residual, x, cache)
        if residual != 0:
            raise RegularityError(f"nonlinear I residual after split: {residual}")
    return real, imag


class FixedChartEvaluator:
    def __init__(self, x: sp.Symbol, substitutions: dict[sp.Symbol, sp.Rational]):
        self.x = x
        self.substitutions = substitutions
        self.cache: dict[str, sp.Expr] = {}
        self.lambda_cache: dict[str, Callable[[np.ndarray], Any]] = {}
        self.nodes: dict[int, np.ndarray] = {}
        self.weights: dict[int, np.ndarray] = {}

    def normalize(self, value: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
        return split_complex(value.subs(self.substitutions), self.x, self.cache)

    def integrate(self, value: sp.Expr, order: int) -> complex:
        real, imag = self.normalize(value)
        return complex(self._integrate_real(real, order), self._integrate_real(imag, order))

    def integrate_matrix(self, matrix: sp.Matrix, order: int) -> np.ndarray:
        output = np.zeros(matrix.shape, dtype=np.complex128)
        for i in range(matrix.rows):
            for j in range(matrix.cols):
                output[i, j] = self.integrate(matrix[i, j], order)
        return output

    def _integrate_real(self, value: sp.Expr, order: int) -> float:
        key = sp.srepr(value)
        if key not in self.lambda_cache:
            self.lambda_cache[key] = sp.lambdify(self.x, value, "numpy")
        if order not in self.nodes:
            nodes, weights = np.polynomial.legendre.leggauss(order)
            self.nodes[order] = nodes
            self.weights[order] = weights
        evaluated = np.asarray(self.lambda_cache[key](self.nodes[order]), dtype=np.float64)
        if evaluated.shape == ():
            evaluated = np.full(order, float(evaluated), dtype=np.float64)
        if not np.all(np.isfinite(evaluated)):
            raise FloatingPointError("nonfinite quadrature integrand value")
        integral = float(np.dot(self.weights[order], evaluated) / 2)
        if not math.isfinite(integral):
            raise FloatingPointError("nonfinite quadrature integral")
        return integral


def realify_rectangular(matrix: np.ndarray, m: int) -> np.ndarray:
    if m < 0:
        raise ValueError("Use nonnegative-m sectors")
    if m == 0:
        if np.max(np.abs(matrix.imag)) > ZERO_TOL:
            raise ValueError(f"m=0 block has imaginary leakage {np.max(np.abs(matrix.imag))}")
        return matrix.real
    real, imag = matrix.real, matrix.imag
    return np.vstack([np.hstack([real, imag]), np.hstack([-imag, real])])


def substitute_jet(jet: dict[str, Any], substitutions: dict[sp.Symbol, sp.Rational]) -> dict[str, Any]:
    return {
        "h": jet["h"].subs(substitutions),
        "nabla": [nabla_a.subs(substitutions) for nabla_a in jet["nabla"]],
        "trace": jet["trace"].subs(substitutions),
        "trace_gradient": [value.subs(substitutions) for value in jet["trace_gradient"]],
        "divergence_lower": [value.subs(substitutions) for value in jet["divergence_lower"]],
    }


def fixed_metric_data(
    data: dict[str, Any],
    substitutions: dict[sp.Symbol, sp.Rational],
) -> dict[str, Any]:
    fixed = dict(data)
    fixed["metric"] = data["metric"].subs(substitutions)
    fixed["inverse_diagonal"] = [
        value.subs(substitutions) for value in data["inverse_diagonal"]
    ]
    fixed["volume_density"] = data["volume_density"].subs(substitutions)
    return fixed


def fixed_density_and_flux_slots(
    data: dict[str, Any],
    substitutions: dict[sp.Symbol, sp.Rational],
) -> tuple[dict[str, dict[str, sp.Matrix]], dict[str, dict[str, sp.Matrix]]]:
    gamma = weak.connection(data)
    tensors = data["spatial_tangents"]
    fixed_data = fixed_metric_data(data, substitutions)
    jets = {
        0: [
            substitute_jet(weak.tensor_jet(tensor, data, 1, 0, gamma), substitutions)
            for tensor in tensors
        ],
        1: [
            substitute_jet(weak.tensor_jet(tensor, data, 0, 1, gamma), substitutions)
            for tensor in tensors
        ],
    }
    harmonics = [
        kinematics.harmonic_amplitudes(ell, data["m"], data["x"])[0]
        for _, ell, _ in data["lapse_labels"]
    ]
    density_potential: dict[str, sp.Matrix] = {}
    density_lapse: dict[str, sp.Matrix] = {}
    flux_potential: dict[str, sp.Matrix] = {}
    flux_lapse: dict[str, sp.Matrix] = {}
    for a in (0, 1):
        for b in (0, 1):
            key = f"{a}{b}"
            density_potential[key] = sp.Matrix(
                [
                    [weak.potential_density(left, right, fixed_data) for right in jets[b]]
                    for left in jets[a]
                ]
            )
            density_lapse[key] = sp.Matrix(
                [
                    [
                        weak.lapse_density(y, int(a == 0), int(a == 1), right, fixed_data)
                        for right in jets[b]
                    ]
                    for y in harmonics
                ]
            )
            flux_potential[key] = sp.Matrix(
                [
                    [weak.potential_boundary_flux(left, right, fixed_data) for right in jets[b]]
                    for left in jets[a]
                ]
            )
            flux_lapse[key] = sp.Matrix(
                [
                    [
                        weak.lapse_boundary_flux(y, int(a == 0), right, fixed_data)
                        for right in jets[b]
                    ]
                    for y in harmonics
                ]
            )
    return (
        {"potential": density_potential, "lapse": density_lapse},
        {"potential": flux_potential, "lapse": flux_lapse},
    )


def polynomial_profile(r: sp.Symbol, index: int, side: str, family: int) -> sp.Expr:
    offset = sp.Rational(index + 1, 7 if side == "left" else 11)
    if side == "left":
        base = r**2 + offset * r + sp.Rational(3, 5)
    elif side == "right":
        base = r**3 - offset * r + sp.Rational(4, 9)
    else:
        base = r**2 - offset * r + sp.Rational(5, 8)
    return base if family == 0 else base + r


def profile_values(
    r: sp.Symbol,
    index: int,
    side: str,
    family: int,
) -> tuple[sp.Expr, sp.Expr]:
    profile = polynomial_profile(r, index, side, family)
    return profile, sp.diff(profile, r)


def dense_spatial_jet(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    side: str,
    family: int,
    gamma: Any,
) -> tuple[dict[str, Any], list[sp.Expr], list[sp.Expr]]:
    r = data["r"]
    values, derivatives = [], []
    field = sp.zeros(3)
    for index, (coefficient, tensor) in enumerate(zip(coeffs, data["spatial_tangents"])):
        profile, derivative = profile_values(r, index, side, family)
        values.append(coefficient * profile)
        derivatives.append(coefficient * derivative)
        field += coefficient * profile * tensor
    return weak.tensor_jet(field, data, 1, 0, gamma), values, derivatives


def dense_lapse_profile(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    family: int,
) -> tuple[sp.Expr, sp.Expr]:
    r = data["r"]
    x = data["x"]
    scalar = 0
    for index, (coefficient, (_, ell, _)) in enumerate(zip(coeffs, data["lapse_labels"])):
        y = kinematics.harmonic_amplitudes(ell, data["m"], x)[0]
        scalar += coefficient * polynomial_profile(r, index, "lapse", family) * y
    return scalar, sp.diff(scalar, r)


def lapse_density_for_scalar(
    scalar: sp.Expr,
    scalar_r: sp.Expr,
    right: dict[str, Any],
    data: dict[str, Any],
) -> sp.Expr:
    gradient = [scalar_r, sp.diff(scalar, data["x"]), sp.I * data["m"] * scalar]
    return -data["volume_density"] * sum(
        data["inverse_diagonal"][j] * weak.chart_conjugate(gradient[j])
        * (right["divergence_lower"][j] - right["trace_gradient"][j])
        for j in range(3)
    ) / 2


def lapse_flux_for_scalar(
    scalar: sp.Expr,
    right: dict[str, Any],
    data: dict[str, Any],
) -> sp.Expr:
    return data["volume_density"] * weak.chart_conjugate(scalar) * data["inverse_diagonal"][0] * (
        right["divergence_lower"][0] - right["trace_gradient"][0]
    ) / 2


def fixed_dense_spatial_jet(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    side: str,
    family: int,
    gamma: Any,
    substitutions: dict[sp.Symbol, sp.Rational],
) -> tuple[dict[str, Any], list[complex], list[complex]]:
    jet, values, derivatives = dense_spatial_jet(data, coeffs, side, family, gamma)
    return (
        substitute_jet(jet, substitutions),
        [complex(sp.N(value.subs(substitutions), 30)) for value in values],
        [complex(sp.N(value.subs(substitutions), 30)) for value in derivatives],
    )


def fixed_dense_lapse_profile(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    family: int,
    substitutions: dict[sp.Symbol, sp.Rational],
) -> tuple[sp.Expr, sp.Expr, list[complex], list[complex]]:
    scalar, scalar_r = dense_lapse_profile(data, coeffs, family)
    values: list[complex] = []
    derivatives: list[complex] = []
    for index, coefficient in enumerate(coeffs):
        profile, derivative = profile_values(data["r"], index, "lapse", family)
        values.append(complex(sp.N((coefficient * profile).subs(substitutions), 30)))
        derivatives.append(complex(sp.N((coefficient * derivative).subs(substitutions), 30)))
    return scalar.subs(substitutions), scalar_r.subs(substitutions), values, derivatives


def combine_slot_matrices(
    matrices: dict[str, np.ndarray],
    left_values: list[complex],
    left_derivatives: list[complex],
    right_values: list[complex],
    right_derivatives: list[complex],
) -> complex:
    left = [np.asarray(left_values), np.asarray(left_derivatives)]
    right = [np.asarray(right_values), np.asarray(right_derivatives)]
    total = 0j
    for a in (0, 1):
        for b in (0, 1):
            total += (left[a].conjugate().T @ matrices[f"{a}{b}"] @ right[b]).item()
    return complex(total)


def combine_lapse_matrices(
    matrices: dict[str, np.ndarray],
    left_values: list[complex],
    left_derivatives: list[complex],
    right_values: list[complex],
    right_derivatives: list[complex],
) -> complex:
    left = [np.asarray(left_values), np.asarray(left_derivatives)]
    right = [np.asarray(right_values), np.asarray(right_derivatives)]
    total = 0j
    for a in (0, 1):
        for b in (0, 1):
            total += (left[a].conjugate().T @ matrices[f"{a}{b}"] @ right[b]).item()
    return complex(total)


def evaluate_slot_family(
    evaluator: FixedChartEvaluator,
    matrices: dict[str, sp.Matrix],
) -> dict[str, dict[str, Any]]:
    evaluated: dict[str, dict[str, Any]] = {}
    for key in SLOTS:
        stage("integrate_slot", slot=key, rows=matrices[key].rows, cols=matrices[key].cols)
        orders = {order: evaluator.integrate_matrix(matrices[key], order) for order in QUADRATURE_ORDERS}
        err24_48 = float(np.max(np.abs(orders[48] - orders[24])))
        err48_96 = float(np.max(np.abs(orders[96] - orders[48])))
        scale = np.maximum(1.0, np.maximum(np.abs(orders[48]), np.abs(orders[96])))
        scaled = float(np.max(np.abs(orders[96] - orders[48]) / scale))
        ratio = None if err24_48 == 0.0 else err48_96 / err24_48
        if scaled > REFINEMENT_RTOL:
            raise AssertionError(f"quadrature refinement failed for slot {key}: {scaled}")
        evaluated[key] = {
            "order_96": matrix_to_json(orders[96]),
            "max_abs_24_48": err24_48,
            "max_abs_48_96": err48_96,
            "max_abs_48_96_over_24_48": ratio,
            "max_abs_error_nonincreasing": err48_96 <= err24_48,
            "max_scaled_48_96": scaled,
            "all_final_integrands_finite": True,
            "regularity_guard": "post-cancel assembled scalar entries have no denominator roots on -1<=x<=1",
        }
    return evaluated


def raw_order96(evaluated: dict[str, dict[str, Any]]) -> dict[str, np.ndarray]:
    return {
        key: np.asarray(
            [[complex(real, imag) for real, imag in row] for row in payload["order_96"]],
            dtype=np.complex128,
        )
        for key, payload in evaluated.items()
    }


def assert_close(value: complex, expected: complex, label: str) -> float:
    if not all(math.isfinite(part) for z in (complex(value), complex(expected))
               for part in (z.real, z.imag)):
        raise FloatingPointError(f"nonfinite comparison: {label}")
    error = abs(value - expected)
    scale = max(1.0, abs(value), abs(expected))
    if error / scale > REFINEMENT_RTOL:
        raise AssertionError(f"{label} mismatch: {value} vs {expected}")
    return float(error / scale)


def assert_nonzero(value: sp.Expr, substitutions: dict[sp.Symbol, sp.Rational], label: str) -> None:
    if sp.simplify(value.subs(substitutions)) == 0:
        raise AssertionError(f"profile family lost two-slot span: {label}")


def profile_span_determinants(
    data: dict[str, Any],
    substitutions: dict[sp.Symbol, sp.Rational],
) -> dict[str, list[str]]:
    r = data["r"]
    records: dict[str, list[str]] = {}
    for side, count in (
        ("left", len(data["spatial_labels"])),
        ("right", len(data["spatial_labels"])),
        ("lapse", len(data["lapse_labels"])),
    ):
        side_records = []
        for index in range(count):
            a, ap = profile_values(r, index, side, 0)
            b, bp = profile_values(r, index, side, 1)
            determinant = sp.cancel((a * bp - b * ap).subs(substitutions))
            if determinant == 0:
                raise AssertionError(f"profile family lost two-slot span: {side} index {index}")
            side_records.append(str(determinant))
        records[side] = side_records
    return records


def verify_rational_guard_controls(x: sp.Symbol) -> dict[str, Any]:
    removable = sp.Mul(1 - x**2, sp.Pow(1 - x**2, -1, evaluate=False), evaluate=False)
    if clean_expr(removable, x, {}) != 1:
        raise AssertionError("removable endpoint factor was not cancelled")
    smooth = x**3 / (x**2 + 2)
    if sp.cancel(clean_expr(smooth, x, {}) - smooth) != 0:
        raise AssertionError("smooth rational guard control changed value")
    rejects = {
        "endpoint_pole": 1 / (1 - x),
        "interior_pole": 1 / x,
        "nonrational_function": sp.sin(x),
        "nonrational_coefficient": sp.pi * x,
        "unrelated_radical": sp.sqrt(x**2),
        "inexact_coefficient": sp.Float("0.1") * x,
        "nonfinite_nan": sp.nan,
        "nonfinite_infinity": sp.oo,
    }
    for label, value in rejects.items():
        try:
            clean_expr(value, x, {})
        except RegularityError:
            continue
        raise AssertionError(f"regularity guard accepted negative control: {label}")
    return {"positive_controls": 2, "rejected_controls": list(rejects)}


def verify_legendre_norms() -> dict[str, Any]:
    x = sp.Symbol("x", real=True)
    regularity = verify_rational_guard_controls(x)
    evaluator = FixedChartEvaluator(x, {})
    checked = []
    for ell in range(L_CUTOFF + 1):
        for m in range(ell + 1):
            y = sp.assoc_legendre(ell, m, x)
            expected = sp.factorial(ell + m) / ((2 * ell + 1) * sp.factorial(ell - m))
            actual = sp.integrate(y * y, (x, -1, 1)) / 2
            if sp.simplify(actual - expected) != 0:
                raise AssertionError(f"Legendre norm failed for {(ell, m)}")
            numeric = {}
            for order in QUADRATURE_ORDERS:
                value = evaluator.integrate(y * y, order)
                error = assert_close(value, complex(expected), f"Legendre quadrature {(ell, m, order)}")
                numeric[str(order)] = {
                    "value": cscalar(value),
                    "scaled_error": error,
                }
            checked.append({
                "ell": ell,
                "m": m,
                "norm_dx_over_2": str(expected),
                "quadrature": numeric,
            })
    return {"checked": checked, "regularity_guard_controls": regularity}


def verify_manual_round_anchor() -> dict[str, Any]:
    data = kinematics.sector(2, 0)
    r = data["r"]
    x = data["x"]
    F, G, Fp, Gp, f, fp = sp.symbols("F G Fp Gp f fp")
    data0 = dict(data)
    data0["epsilon"] = sp.S.Zero
    data0["metric"] = data["metric"].subs(data["epsilon"], 0)
    data0["inverse_diagonal"] = [value.subs(data["epsilon"], 0) for value in data["inverse_diagonal"]]
    data0["volume_density"] = data["volume_density"].subs(data["epsilon"], 0)
    data0["spatial_tangents"] = [tensor.subs(data["epsilon"], 0) for tensor in data["spatial_tangents"]]
    gamma = weak.connection(data0)
    lambda0 = data0["spatial_tangents"][0]
    left = weak.tensor_jet(lambda0, data0, F, Fp, gamma)
    right = weak.tensor_jet(lambda0, data0, G, Gp, gamma)
    potential = sp.cancel(weak.potential_density(left, right, data0) + 2 * F * G + r * (F * Gp + Fp * G))
    potential_flux = sp.cancel(weak.potential_boundary_flux(left, right, data0) - r * F * G)
    lapse = sp.cancel(weak.lapse_density(sp.S.One, f, fp, right, data0) + r * fp * G)
    lapse_flux = sp.cancel(weak.lapse_boundary_flux(sp.S.One, f, right, data0) - r * f * G)
    residuals = [potential, potential_flux, lapse, lapse_flux]
    if any(residual != 0 for residual in residuals):
        raise AssertionError(f"round manual anchor failed: {residuals}")
    if any(residual.has(x) for residual in residuals):
        raise AssertionError("round manual anchor retained angular dependence")
    return {
        "epsilon": "0",
        "m": 0,
        "lambda_ell0": "passed",
        "potential_density": "-2*F*G - r*(F*Gprime + Fprime*G)",
        "potential_slot_V00": "-2",
        "potential_slot_V01": "-r",
        "potential_slot_V10": "-r",
        "potential_slot_V11": "0",
        "potential_flux": "r*F*G",
        "potential_flux_slot_00": "r",
        "lapse_density": "-r*fprime*G",
        "lapse_slot_10": "-r",
        "lapse_flux": "r*f*G",
        "lapse_flux_slot_00": "r",
    }


def verify_realification_controls() -> dict[str, Any]:
    square = np.array([[1 + 2j, 3 + 5j], [7 - 11j, 13 - 17j]], dtype=np.complex128)
    rectangular = np.array([[1j, 2 - 3j]], dtype=np.complex128)
    square_real = realify_rectangular(square, 1)
    rectangular_real = realify_rectangular(rectangular, 1)
    wrong = np.vstack(
        [np.hstack([rectangular.real, -rectangular.imag]),
         np.hstack([rectangular.imag, rectangular.real])]
    )
    if np.allclose(rectangular_real, wrong):
        raise AssertionError("wrong-sign realification was not rejected")
    rectangular_expected = np.array([[0.0, 2.0, 1.0, -3.0], [-1.0, 3.0, 0.0, 2.0]])
    if not np.array_equal(rectangular_real, rectangular_expected):
        raise AssertionError("manufactured imaginary rectangular anchor failed")

    direct = direct_positive_sine_phi_average(square, 12)
    if not np.allclose(square_real, direct, rtol=0, atol=1e-12):
        raise AssertionError("positive-sine direct real-field phi average failed")
    return {
        "positive_sine_square": "passed",
        "positive_sine_rectangular": "passed",
        "direct_phi_nodes": 12,
        "direct_phi_oracle": (
            "uniform average of sqrt(2)*cos(phi), sqrt(2)*sin(phi) with "
            "Re(H)*f*g + Im(H)*f*partial_phi(g)"
        ),
        "wrong_sign_rejected": True,
    }


def direct_positive_sine_phi_average(matrix: np.ndarray, nodes: int) -> np.ndarray:
    phis = np.arange(nodes, dtype=np.float64) * (2 * math.pi / nodes)
    basis = np.vstack([np.sqrt(2) * np.cos(phis), np.sqrt(2) * np.sin(phis)])
    basis_phi = np.vstack([-np.sqrt(2) * np.sin(phis), np.sqrt(2) * np.cos(phis)])
    real = matrix.real
    imag = matrix.imag
    output = np.zeros((2 * matrix.shape[0], 2 * matrix.shape[1]))
    for left_component in range(2):
        for right_component in range(2):
            for i in range(matrix.shape[0]):
                for j in range(matrix.shape[1]):
                    density = (
                        real[i, j] * basis[left_component] * basis[right_component]
                        + imag[i, j] * basis[left_component] * basis_phi[right_component]
                    )
                    output[
                        left_component * matrix.shape[0] + i,
                        right_component * matrix.shape[1] + j,
                    ] = float(np.mean(density))
    return output


def chart_real(value: sp.Expr) -> sp.Expr:
    return (value + weak.chart_conjugate(value)) / 2


def chart_imag(value: sp.Expr) -> sp.Expr:
    return (value - weak.chart_conjugate(value)) / (2 * sp.I)


def realify_phased_expr(value: sp.Expr, phase: sp.Expr, component: int) -> sp.Expr:
    phased = value * phase
    projected = chart_real(phased) if component == 0 else chart_imag(phased)
    return sp.sqrt(2) * projected


def realify_jet(jet: dict[str, Any], phase: sp.Expr, component: int) -> dict[str, Any]:
    return {
        "h": jet["h"].applyfunc(lambda value: realify_phased_expr(value, phase, component)),
        "nabla": [
            nabla_a.applyfunc(lambda value: realify_phased_expr(value, phase, component))
            for nabla_a in jet["nabla"]
        ],
        "trace": realify_phased_expr(jet["trace"], phase, component),
        "trace_gradient": [
            realify_phased_expr(value, phase, component) for value in jet["trace_gradient"]
        ],
        "divergence_lower": [
            realify_phased_expr(value, phase, component) for value in jet["divergence_lower"]
        ],
    }


def realify_lapse_gradient(
    scalar: sp.Expr,
    scalar_r: sp.Expr,
    data: dict[str, Any],
    phase: sp.Expr,
    component: int,
) -> list[sp.Expr]:
    return [
        realify_phased_expr(scalar_r, phase, component),
        realify_phased_expr(sp.diff(scalar, data["x"]), phase, component),
        realify_phased_expr(sp.I * data["m"] * scalar, phase, component),
    ]


def lapse_density_for_gradient(
    gradient: list[sp.Expr],
    right: dict[str, Any],
    data: dict[str, Any],
) -> sp.Expr:
    return -data["volume_density"] * sum(
        data["inverse_diagonal"][j] * weak.chart_conjugate(gradient[j])
        * (right["divergence_lower"][j] - right["trace_gradient"][j])
        for j in range(3)
    ) / 2


def max_case_refinement(case: dict[str, Any]) -> float:
    return max(
        slot["max_scaled_48_96"]
        for family in (
            "potential_density",
            "lapse_density",
            "potential_radial_flux",
            "lapse_radial_flux",
        )
        for slot in case[family].values()
    )


def hermiticity_residuals(blocks: dict[str, np.ndarray]) -> dict[str, float]:
    residuals = {
        "B00_minus_B00_dagger": float(np.max(np.abs(blocks["00"] - blocks["00"].conjugate().T))),
        "B11_minus_B11_dagger": float(np.max(np.abs(blocks["11"] - blocks["11"].conjugate().T))),
        "B01_minus_B10_dagger": float(np.max(np.abs(blocks["01"] - blocks["10"].conjugate().T))),
    }
    if max(residuals.values()) > 2e-9:
        raise AssertionError(f"potential Hermiticity failed: {residuals}")
    return residuals


def structural_zero_residuals(
    potential_flux: dict[str, np.ndarray],
    lapse_flux: dict[str, np.ndarray],
) -> dict[str, dict[str, float]]:
    residuals = {
        "potential_radial": {
            key: float(np.max(np.abs(potential_flux[key]))) for key in ("10", "11")
        },
        "lapse_radial": {
            key: float(np.max(np.abs(lapse_flux[key]))) for key in ("10", "11")
        },
    }
    largest = max(value for family in residuals.values() for value in family.values())
    if largest > ZERO_TOL:
        raise AssertionError(f"radial flux structural zero failed: {residuals}")
    return residuals


def realfield_oracle(
    evaluator: FixedChartEvaluator,
    fixed_data: dict[str, Any],
    left_spatial: list[dict[str, Any]],
    right_spatial: list[dict[str, Any]],
    left_lapse: list[tuple[sp.Expr, sp.Expr]],
    expected: dict[str, dict[str, np.ndarray]],
    m: int,
) -> dict[str, Any]:
    if m == 0:
        return {"performed": False, "reason": "m=0 is already real"}

    phases = [sp.I ** k for k in range(4)]
    errors: dict[str, float] = {}
    for name, density, left_items, right_items, expected_matrix in (
        ("potential_density", weak.potential_density, left_spatial, right_spatial,
         expected["potential_density"]),
        ("potential_radial_flux", weak.potential_boundary_flux, left_spatial, right_spatial,
         expected["potential_radial_flux"]),
    ):
        actual = np.zeros((4, 4), dtype=np.complex128)
        for lf, left in enumerate(left_items):
            for rf, right in enumerate(right_items):
                for lc in range(2):
                    for rc in range(2):
                        averaged_density = sum(
                            density(
                                realify_jet(left, phase, lc),
                                realify_jet(right, phase, rc),
                                fixed_data,
                            )
                            for phase in phases
                        ) / len(phases)
                        actual[lc * 2 + lf, rc * 2 + rf] = evaluator.integrate(
                            averaged_density, 96
                        )
        target = realify_rectangular(expected_matrix, m)
        errors[name] = float(np.max(np.abs(actual - target)))

    for name, density, expected_matrix in (
        ("lapse_density", lapse_density_for_scalar, expected["lapse_density"]),
        ("lapse_radial_flux", lapse_flux_for_scalar, expected["lapse_radial_flux"]),
    ):
        actual = np.zeros((4, 4), dtype=np.complex128)
        for lf, (scalar, scalar_r) in enumerate(left_lapse):
            for rf, right in enumerate(right_spatial):
                for lc in range(2):
                    for rc in range(2):
                        averaged_density = 0
                        for phase in phases:
                            left_scalar = realify_phased_expr(scalar, phase, lc)
                            right_jet = realify_jet(right, phase, rc)
                            if name == "lapse_density":
                                value = lapse_density_for_gradient(
                                    realify_lapse_gradient(
                                        scalar, scalar_r, fixed_data, phase, lc
                                    ),
                                    right_jet,
                                    fixed_data,
                                )
                            else:
                                value = density(left_scalar, right_jet, fixed_data)
                            averaged_density += value / len(phases)
                        actual[lc * 2 + lf, rc * 2 + rf] = evaluator.integrate(
                            averaged_density, 96
                        )
        target = realify_rectangular(expected_matrix, m)
        errors[name] = float(np.max(np.abs(actual - target)))

    largest = max(errors.values())
    if largest > 1e-9:
        raise AssertionError(f"actual positive-sine real-field oracle failed: {errors}")
    return {
        "performed": True,
        "phase_average": "exact fourth roots of unity; density averaged before x quadrature",
        "phase_count": 4,
        "equivalent_uniform_phi_nodes": 4 * m,
        "max_abs_errors": errors,
    }


def evaluate_case(m: int, r_value: sp.Rational) -> dict[str, Any]:
    started = now()
    data = kinematics.sector(L_CUTOFF, m)
    x = data["x"]
    substitutions = {data["r"]: r_value, data["epsilon"]: EPSILON_VALUE}
    fixed_data = fixed_metric_data(data, substitutions)
    evaluator = FixedChartEvaluator(x, substitutions)

    tick = substage("build_symbolic_slots", "start", None, m=m, r=str(r_value))
    density_slots, flux_slots = fixed_density_and_flux_slots(data, substitutions)
    substage("build_symbolic_slots", "done", tick, m=m, r=str(r_value))

    tick = substage("cancel_and_quad_density", "start", None, m=m, r=str(r_value))
    potential_density = evaluate_slot_family(evaluator, density_slots["potential"])
    lapse_density = evaluate_slot_family(evaluator, density_slots["lapse"])
    substage("cancel_and_quad_density", "done", tick, m=m, r=str(r_value))

    tick = substage("cancel_and_quad_flux", "start", None, m=m, r=str(r_value))
    potential_flux = evaluate_slot_family(evaluator, flux_slots["potential"])
    lapse_flux = evaluate_slot_family(evaluator, flux_slots["lapse"])
    substage("cancel_and_quad_flux", "done", tick, m=m, r=str(r_value))

    potential_density96 = raw_order96(potential_density)
    lapse_density96 = raw_order96(lapse_density)
    potential_flux96 = raw_order96(potential_flux)
    lapse_flux96 = raw_order96(lapse_flux)
    flux_structural_zeros = structural_zero_residuals(potential_flux96, lapse_flux96)

    nq = len(data["spatial_labels"])
    na = len(data["lapse_labels"])
    profile_spans = profile_span_determinants(data, substitutions)
    left_coeffs = [sp.Rational(i + 1, 13) for i in range(nq)]
    right_coeffs = [sp.Rational(2 * i + 3, 17) for i in range(nq)]
    lapse_coeffs = [sp.Rational(3 * i + 1, 19) for i in range(na)]

    tick = substage("direct_profile_oracle", "start", None, m=m, r=str(r_value))
    gamma = weak.connection(data)
    left_families = [
        fixed_dense_spatial_jet(data, left_coeffs, "left", family, gamma, substitutions)
        for family in range(2)
    ]
    right_families = [
        fixed_dense_spatial_jet(data, right_coeffs, "right", family, gamma, substitutions)
        for family in range(2)
    ]
    lapse_families = [
        fixed_dense_lapse_profile(data, lapse_coeffs, family, substitutions)
        for family in range(2)
    ]

    reconstruction_errors: dict[str, list[float]] = {
        "potential_density_scaled_error": [],
        "potential_flux_scaled_error": [],
        "lapse_density_scaled_error": [],
        "lapse_flux_scaled_error": [],
    }
    family_matrices = {
        "potential_density": np.zeros((2, 2), dtype=np.complex128),
        "potential_radial_flux": np.zeros((2, 2), dtype=np.complex128),
        "lapse_density": np.zeros((2, 2), dtype=np.complex128),
        "lapse_radial_flux": np.zeros((2, 2), dtype=np.complex128),
    }
    for lf, (left_jet, left_v, left_d) in enumerate(left_families):
        for rf, (right_jet, right_v, right_d) in enumerate(right_families):
            direct_potential_density = evaluator.integrate(
                weak.potential_density(left_jet, right_jet, fixed_data), 96
            )
            direct_potential_flux = evaluator.integrate(
                weak.potential_boundary_flux(left_jet, right_jet, fixed_data), 96
            )
            family_matrices["potential_density"][lf, rf] = direct_potential_density
            family_matrices["potential_radial_flux"][lf, rf] = direct_potential_flux
            reconstruction_errors["potential_density_scaled_error"].append(assert_close(
                direct_potential_density,
                combine_slot_matrices(potential_density96, left_v, left_d, right_v, right_d),
                f"potential density reconstruction m={m} r={r_value} {(lf, rf)}",
            ))
            reconstruction_errors["potential_flux_scaled_error"].append(assert_close(
                direct_potential_flux,
                combine_slot_matrices(potential_flux96, left_v, left_d, right_v, right_d),
                f"potential flux reconstruction m={m} r={r_value} {(lf, rf)}",
            ))

    for lf, (lapse_scalar, lapse_scalar_r, lapse_v, lapse_d) in enumerate(lapse_families):
        for rf, (right_jet, right_v, right_d) in enumerate(right_families):
            direct_lapse_density = evaluator.integrate(
                lapse_density_for_scalar(lapse_scalar, lapse_scalar_r, right_jet, fixed_data), 96
            )
            direct_lapse_flux = evaluator.integrate(
                lapse_flux_for_scalar(lapse_scalar, right_jet, fixed_data), 96
            )
            family_matrices["lapse_density"][lf, rf] = direct_lapse_density
            family_matrices["lapse_radial_flux"][lf, rf] = direct_lapse_flux
            reconstruction_errors["lapse_density_scaled_error"].append(assert_close(
                direct_lapse_density,
                combine_lapse_matrices(lapse_density96, lapse_v, lapse_d, right_v, right_d),
                f"lapse density reconstruction m={m} r={r_value} {(lf, rf)}",
            ))
            reconstruction_errors["lapse_flux_scaled_error"].append(assert_close(
                direct_lapse_flux,
                combine_lapse_matrices(lapse_flux96, lapse_v, lapse_d, right_v, right_d),
                f"lapse flux reconstruction m={m} r={r_value} {(lf, rf)}",
            ))
    reconstruction = {
        key: {"max_scaled_error": max(errors), "combinations": len(errors)}
        for key, errors in reconstruction_errors.items()
    }
    substage("direct_profile_oracle", "done", tick, m=m, r=str(r_value))

    tick = substage("hermiticity_and_realify", "start", None, m=m, r=str(r_value))
    hermiticity = hermiticity_residuals(potential_density96)
    real_shapes = {
        "potential": {
            key: list(realify_rectangular(value, m).shape)
            for key, value in potential_density96.items()
        },
        "lapse": {
            key: list(realify_rectangular(value, m).shape)
            for key, value in lapse_density96.items()
        },
    }
    real_matrices = {
        "potential_density": {
            key: real_matrix_to_json(realify_rectangular(value, m))
            for key, value in potential_density96.items()
        },
        "lapse_density": {
            key: real_matrix_to_json(realify_rectangular(value, m))
            for key, value in lapse_density96.items()
        },
        "potential_radial_flux": {
            key: real_matrix_to_json(realify_rectangular(value, m))
            for key, value in potential_flux96.items()
        },
        "lapse_radial_flux": {
            key: real_matrix_to_json(realify_rectangular(value, m))
            for key, value in lapse_flux96.items()
        },
    }
    realfield = realfield_oracle(
        evaluator,
        fixed_data,
        [item[0] for item in left_families],
        [item[0] for item in right_families],
        [(item[0], item[1]) for item in lapse_families],
        family_matrices,
        m,
    )
    substage("hermiticity_and_realify", "done", tick, m=m, r=str(r_value))

    return {
        "m": m,
        "r": str(r_value),
        "epsilon": str(EPSILON_VALUE),
        "spatial_labels": [label_tuple(label) for label in data["spatial_labels"]],
        "lapse_labels": [label_tuple(label) for label in data["lapse_labels"]],
        "real_spatial_labels": real_labels(data["spatial_labels"], m),
        "real_lapse_labels": real_labels(data["lapse_labels"], m),
        "all_pairs": True,
        "quadrature_orders": list(QUADRATURE_ORDERS),
        "all_final_integrands_finite": True,
        "nq": nq,
        "na": na,
        "potential_density": potential_density,
        "lapse_density": lapse_density,
        "potential_radial_flux": potential_flux,
        "lapse_radial_flux": lapse_flux,
        "potential_hermiticity": hermiticity,
        "direct_combined_profile_reconstruction": reconstruction,
        "profile_span_determinants": profile_spans,
        "real_block_shapes": real_shapes,
        "realified_matrices": real_matrices,
        "actual_realfield_oracle": realfield,
        "structural_zero_flux_slots": flux_structural_zeros,
        "elapsed_seconds": now() - started,
    }


def verify_export(output_dir: Path | None) -> dict[str, Any]:
    started = now()
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / PROGRESS_NAME).write_text("")

    write_progress(output_dir, {
        "phase": "setup",
        "step": 1,
        "of": 1,
        "unit": "environment",
        "plan": runtime_plan(),
    })
    write_progress(output_dir, {
        "phase": "controls",
        "step": 0,
        "of": CONTROL_BATCHES,
        "unit": "batches",
        "m": {"completed_control_batches": 0},
    })

    controls: dict[str, Any] = {}
    control_steps: tuple[tuple[str, Callable[[], dict[str, Any]]], ...] = (
        ("legendre_norms", verify_legendre_norms),
        ("round_weak_anchor", verify_manual_round_anchor),
        ("manufactured_realification", verify_realification_controls),
    )
    for index, (name, control) in enumerate(control_steps, start=1):
        tick = substage(name, "start")
        controls[name] = control()
        substage(name, "done", tick)
        write_checkpoint(output_dir, started, f"control_{name}", controls[name])
        write_progress(output_dir, {
            "phase": "controls",
            "step": index,
            "of": CONTROL_BATCHES,
            "unit": "batches",
            "m": {"completed_control_batches": index},
        })
    write_checkpoint(output_dir, started, "controls", controls)

    cases = []
    completed = 0
    max_scaled_refinement_error = 0.0
    write_progress(output_dir, {
        "phase": "angular_export",
        "step": completed,
        "of": ANGULAR_CASES,
        "unit": "snapshots",
        "m": {
            "completed_snapshots": completed,
            "max_scaled_refinement_error": max_scaled_refinement_error,
        },
    })
    for m in M_VALUES:
        for r_value in R_VALUES:
            stage("angular_case", m=m, r=str(r_value))
            case = evaluate_case(m, r_value)
            if output_dir is not None:
                safe_write_json(
                    output_dir / f"angular-snapshot-m{m}-r{str(r_value).replace('/', '_')}.json",
                    case,
                )
            cases.append(case)
            completed += 1
            max_scaled_refinement_error = max(
                max_scaled_refinement_error,
                max_case_refinement(case),
            )
            write_checkpoint(
                output_dir,
                started,
                f"angular_m{m}_r{r_value}",
                {
                    "m": m,
                    "r": str(r_value),
                    "max_scaled_refinement_error": max_scaled_refinement_error,
                    "elapsed_seconds": case["elapsed_seconds"],
                },
            )
            write_progress(output_dir, {
                "phase": "angular_export",
                "step": completed,
                "of": ANGULAR_CASES,
                "unit": "snapshots",
                "m": {
                    "completed_snapshots": completed,
                    "max_scaled_refinement_error": max_scaled_refinement_error,
                },
            })

    result = {
        "status": "passed",
        "scope": "fixed-radius angular weak-slot diagnostic snapshots",
        "kind": "angular_snapshot",
        "diagnostic_snapshots": True,
        "L": L_CUTOFF,
        "m_values": list(M_VALUES),
        "r_values": [str(value) for value in R_VALUES],
        "epsilon": str(EPSILON_VALUE),
        "angular_measure": "dx/2 after equal-m azimuthal averaging",
        "slot_convention": "ab: radial derivative order a on conjugated left profile, b on right",
        "formal_slot_mapping": {"P": "11", "Q": "01", "R": "10", "S": "00"},
        "realification": "positive-sine [[Re(H),Im(H)],[-Im(H),Re(H)]]",
        "quadrature": {
            "orders": list(QUADRATURE_ORDERS),
            "precision": "NumPy float64 Gauss-Legendre",
            "scaled_refinement_threshold_48_vs_96": REFINEMENT_RTOL,
            "rigorous_error_bound_claimed": False,
        },
        "density_slots": {
            "potential": list(SLOTS),
            "lapse": list(SLOTS),
        },
        "radial_flux_slots": {
            "potential": list(SLOTS),
            "lapse": list(SLOTS),
            "structural_zero_left_derivative_slots": ["10", "11"],
            "separate_boundary_bilinear": True,
        },
        "negative_scope": {
            "radial_projection_performed": False,
            "finite_KED_performed": False,
            "operator_spectrum_performed": False,
            "domain_selected": False,
            "full_radial_operator_export": False,
            "finite_assembly": False,
            "radial_projection": False,
            "boundary_domain": False,
            "kernel_preservation": False,
            "spectrum": False,
            "seed": False,
            "nonlinear_theorem": False,
        },
        "controls": controls,
        "cases": cases,
        "max_scaled_refinement_error": max_scaled_refinement_error,
        "source_sha256": source_sha256(),
        "python": platform.python_version(),
        "dependencies": {"numpy": np.__version__, "sympy": sp.__version__},
        "elapsed_seconds": now() - started,
    }
    if output_dir is not None:
        write_progress(output_dir, {
            "phase": "export",
            "step": 0,
            "of": 1,
            "unit": "checker result",
        })
        safe_write_json(output_dir / RESULT_NAME, result)
        write_progress(output_dir, {
            "phase": "export",
            "step": 1,
            "of": 1,
            "unit": "checker result",
            "m": {
                "angular_export_passed": 1,
                "elapsed_seconds": now() - started,
            },
        })
    print("PASS angular weak-slot diagnostic snapshots", flush=True)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    verify_export(args.output_dir)


if __name__ == "__main__":
    main()

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
PROBE_RESULT_NAME = "direct-profile-probe-checks.json"
ARITHMETIC_PROBE_RESULT_NAME = "profile-arithmetic-probe-checks.json"
COMPACT_PROFILE_PROBE_RESULT_NAME = "compact-profile-probe-checks.json"
PROGRESS_NAME = "angular-weak-slot-progress.jsonl"
CHECKPOINT_NAME = "angular-weak-slot-checkpoint.json"
DENOMINATOR_CACHE: set[tuple[sp.Rational, ...]] = set()
ARITHMETIC_PROBE_PAIRS = 12


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


def reset_denominator_cache() -> int:
    size = len(DENOMINATOR_CACHE)
    DENOMINATOR_CACHE.clear()
    return size


def runtime_plan(snapshot_count: int = ANGULAR_CASES) -> list[dict[str, Any]]:
    return [
        {"phase": "setup", "of": 1, "unit": "environment"},
        {"phase": "controls", "of": CONTROL_BATCHES, "unit": "batches"},
        {"phase": "angular_export", "of": snapshot_count, "unit": "snapshots"},
        {"phase": "export", "of": 1, "unit": "checker result"},
    ]


def arithmetic_runtime_plan() -> list[dict[str, Any]]:
    return [
        {"phase": "setup", "of": 1, "unit": "environment"},
        {"phase": "controls", "of": CONTROL_BATCHES, "unit": "batches"},
        {"phase": "arithmetic", "of": ARITHMETIC_PROBE_PAIRS, "unit": "pairs"},
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


def _validate_qq_i_rational(value: sp.Expr, x: sp.Symbol) -> None:
    qq_i = sp.QQ.algebraic_field(sp.I)
    numerator, denominator = sp.fraction(value)
    try:
        sp.Poly(numerator, x, domain=qq_i)
        denominator_poly = sp.Poly(denominator, x, domain=qq_i)
    except Exception as exc:
        raise RegularityError(f"component is not rational over QQ(I)(x): {value}") from exc
    if denominator_poly.is_zero:
        raise RegularityError(f"component has zero denominator: {value}")


def compact_component(value: sp.Expr, x: sp.Symbol) -> sp.Expr:
    if value.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise RegularityError(f"nonexact or nonfinite component: {value}")
    try:
        compacted = sp.cancel(sp.together(value))
    except Exception as exc:
        raise RegularityError(f"could not compact component: {value}") from exc
    if compacted.has(sp.Float, sp.nan, sp.zoo, sp.oo, -sp.oo):
        raise RegularityError(f"nonexact or nonfinite compacted component: {compacted}")
    for power in compacted.atoms(sp.Pow):
        if power.exp.is_Rational and not power.exp.is_Integer:
            raise RegularityError(f"unsupported component radical: {power}")
    if compacted.has(sp.Abs):
        raise RegularityError(f"unsupported component Abs: {compacted}")
    _validate_qq_i_rational(compacted, x)
    return compacted


def compact_matrix(matrix: sp.Matrix, x: sp.Symbol) -> sp.Matrix:
    return matrix.applyfunc(lambda value: compact_component(value, x))


def compact_jet(jet: dict[str, Any], x: sp.Symbol) -> dict[str, Any]:
    return {
        "h": compact_matrix(jet["h"], x),
        "nabla": [compact_matrix(nabla_a, x) for nabla_a in jet["nabla"]],
        "trace": compact_component(jet["trace"], x),
        "trace_gradient": [
            compact_component(value, x) for value in jet["trace_gradient"]
        ],
        "divergence_lower": [
            compact_component(value, x) for value in jet["divergence_lower"]
        ],
    }


def assert_compacted_jet_equal(
    raw: dict[str, Any],
    compacted: dict[str, Any],
    x: sp.Symbol,
    cache: dict[str, sp.Expr],
) -> int:
    assert_jet_shape(raw, "equality_raw")
    assert_jet_shape(compacted, "equality_compact")
    checks = 0
    for raw_value, compact_value in zip(list(raw["h"]), list(compacted["h"])):
        if clean_expr(raw_value - compact_value, x, cache) != 0:
            raise AssertionError("compacted h component changed")
        checks += 1
    for raw_nabla, compact_nabla in zip(raw["nabla"], compacted["nabla"]):
        for raw_value, compact_value in zip(list(raw_nabla), list(compact_nabla)):
            if clean_expr(raw_value - compact_value, x, cache) != 0:
                raise AssertionError("compacted nabla component changed")
            checks += 1
    if clean_expr(raw["trace"] - compacted["trace"], x, cache) != 0:
        raise AssertionError("compacted trace changed")
    checks += 1
    for raw_value, compact_value in zip(
        raw["trace_gradient"], compacted["trace_gradient"]
    ):
        if clean_expr(raw_value - compact_value, x, cache) != 0:
            raise AssertionError("compacted trace gradient changed")
        checks += 1
    for raw_value, compact_value in zip(
        raw["divergence_lower"], compacted["divergence_lower"]
    ):
        if clean_expr(raw_value - compact_value, x, cache) != 0:
            raise AssertionError("compacted divergence changed")
        checks += 1
    return checks


def verify_component_compaction_controls(x: sp.Symbol) -> dict[str, Any]:
    positives = {
        "endpoint_pole_allowed": 1 / (1 - x),
        "removable_endpoint_pole": (1 - x**2) / (1 - x),
        "imaginary_rational": sp.I * x / (x**2 + 2),
    }
    for label, value in positives.items():
        compacted = compact_component(value, x)
        if clean_expr(value - compacted, x, {}) != 0:
            raise AssertionError(f"component compaction changed {label}")
    try:
        clean_expr(positives["endpoint_pole_allowed"], x, {})
    except RegularityError:
        pass
    else:
        raise AssertionError("final scalar guard accepted a component endpoint pole")
    rejects = {
        "unsupported_function": sp.sin(x),
        "unsupported_radical": sp.sqrt(x),
        "inexact_float": sp.Float("0.1") * x,
        "nonfinite_nan": sp.nan,
        "nonfinite_infinity": sp.oo,
    }
    for label, value in rejects.items():
        try:
            compact_component(value, x)
        except RegularityError:
            continue
        raise AssertionError(f"component compaction accepted negative control: {label}")
    return {
        "positive_controls": list(positives),
        "component_pole_rejected_by_final_scalar_guard": True,
        "rejected_controls": list(rejects),
        "endpoint_poles_allowed": True,
        "final_density_regularization_applied_to_components": False,
    }


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


class TimedChartEvaluator(FixedChartEvaluator):
    def __init__(
        self,
        x: sp.Symbol,
        substitutions: dict[sp.Symbol, sp.Rational],
        context: dict[str, Any] | None = None,
    ):
        super().__init__(x, substitutions)
        self.context = context or {}

    def timed_normalize(
        self,
        value: sp.Expr,
        context: dict[str, Any],
    ) -> tuple[sp.Expr, sp.Expr, dict[str, float]]:
        timings: dict[str, float] = {}
        tick = substage("fixed_substitution", "start", None, **context)
        substituted = value.subs(self.substitutions)
        timings["fixed_substitution_seconds"] = now() - tick
        substage("fixed_substitution", "done", tick, **context)

        tick = substage("expand_split", "start", None, **context)
        expanded = sp.expand(substituted)
        timings["expand_split_seconds"] = now() - tick
        substage("expand_split", "done", tick, **context)

        tick = substage("cancel_qq_root_guard", "start", None, **context)
        real = clean_expr(expanded.coeff(sp.I, 0), self.x, self.cache)
        imag = clean_expr(expanded.coeff(sp.I, 1), self.x, self.cache)
        residual = sp.expand(expanded - real - sp.I * imag)
        if residual != 0:
            residual = clean_expr(residual, self.x, self.cache)
            if residual != 0:
                raise RegularityError(f"nonlinear I residual after split: {residual}")
        timings["cancel_qq_root_guard_seconds"] = now() - tick
        substage("cancel_qq_root_guard", "done", tick, **context)
        return real, imag, timings

    def timed_integrate(
        self,
        value: sp.Expr,
        order: int,
    ) -> tuple[complex, dict[str, float], sp.Expr, sp.Expr]:
        context = {**self.context, "order": order}
        real, imag, timings = self.timed_normalize(value, context)
        tick = substage("legendre_nodes", "start", None, order=order, **self.context)
        if order not in self.nodes:
            nodes, weights = np.polynomial.legendre.leggauss(order)
            self.nodes[order] = nodes
            self.weights[order] = weights
        timings["legendre_nodes_seconds"] = now() - tick
        substage("legendre_nodes", "done", tick, order=order, **self.context)

        parts = []
        lambdify_seconds = 0.0
        evaluation_seconds = 0.0
        lambdify_tick = substage("lambdify", "start", None, order=order, **self.context)
        for value_part in (real, imag):
            key = sp.srepr(value_part)
            tick = now()
            if key not in self.lambda_cache:
                self.lambda_cache[key] = sp.lambdify(self.x, value_part, "numpy")
            lambdify_seconds += now() - tick
        substage("lambdify", "done", lambdify_tick, order=order, **self.context)
        quadrature_tick = substage("quadrature", "start", None, order=order, **self.context)
        for value_part in (real, imag):
            key = sp.srepr(value_part)
            tick = now()
            evaluated = np.asarray(self.lambda_cache[key](self.nodes[order]), dtype=np.float64)
            if evaluated.shape == ():
                evaluated = np.full(order, float(evaluated), dtype=np.float64)
            if not np.all(np.isfinite(evaluated)):
                raise FloatingPointError("nonfinite quadrature integrand value")
            integral = float(np.dot(self.weights[order], evaluated) / 2)
            if not math.isfinite(integral):
                raise FloatingPointError("nonfinite quadrature integral")
            parts.append(integral)
            evaluation_seconds += now() - tick
        substage("quadrature", "done", quadrature_tick, order=order, **self.context)
        timings["lambdify_seconds"] = lambdify_seconds
        timings["quadrature_evaluation_seconds"] = evaluation_seconds
        return complex(parts[0], parts[1]), timings, real, imag


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


def fixed_profile_spatial_jets(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    side: str,
    family: int,
    gamma: Any,
    substitutions: dict[sp.Symbol, sp.Rational],
) -> tuple[list[dict[str, Any]], list[complex], list[complex]]:
    r = data["r"]
    jets = []
    values: list[complex] = []
    derivatives: list[complex] = []
    for index, (coefficient, tensor) in enumerate(zip(coeffs, data["spatial_tangents"])):
        profile, derivative = profile_values(r, index, side, family)
        field = coefficient * profile * tensor
        jet = weak.tensor_jet(field, data, 1, 0, gamma)
        jets.append(substitute_jet(jet, substitutions))
        values.append(complex(sp.N((coefficient * profile).subs(substitutions), 30)))
        derivatives.append(complex(sp.N((coefficient * derivative).subs(substitutions), 30)))
    return jets, values, derivatives


def fixed_single_profile_spatial_jet(
    data: dict[str, Any],
    coeff: sp.Expr,
    side: str,
    family: int,
    index: int,
    gamma: Any,
    substitutions: dict[sp.Symbol, sp.Rational],
) -> dict[str, Any]:
    r = data["r"]
    profile, _ = profile_values(r, index, side, family)
    field = coeff * profile * data["spatial_tangents"][index]
    return substitute_jet(weak.tensor_jet(field, data, 1, 0, gamma), substitutions)


def assert_jet_shape(jet: dict[str, Any], label: str) -> None:
    if jet["h"].shape != (3, 3):
        raise AssertionError(f"{label} h has shape {jet['h'].shape}")
    if len(jet["nabla"]) != 3:
        raise AssertionError(f"{label} nabla has {len(jet['nabla'])} entries")
    for index, matrix in enumerate(jet["nabla"]):
        if matrix.shape != (3, 3):
            raise AssertionError(f"{label} nabla[{index}] has shape {matrix.shape}")
    if len(jet["trace_gradient"]) != 3:
        raise AssertionError(
            f"{label} trace_gradient has {len(jet['trace_gradient'])} entries"
        )
    if len(jet["divergence_lower"]) != 3:
        raise AssertionError(
            f"{label} divergence_lower has {len(jet['divergence_lower'])} entries"
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


def fixed_profile_lapse_profiles(
    data: dict[str, Any],
    coeffs: list[sp.Expr],
    family: int,
    substitutions: dict[sp.Symbol, sp.Rational],
) -> tuple[list[tuple[sp.Expr, sp.Expr]], list[complex], list[complex]]:
    r = data["r"]
    x = data["x"]
    profiles = []
    values: list[complex] = []
    derivatives: list[complex] = []
    for index, (coefficient, (_, ell, _)) in enumerate(zip(coeffs, data["lapse_labels"])):
        y = kinematics.harmonic_amplitudes(ell, data["m"], x)[0]
        profile, derivative = profile_values(r, index, "lapse", family)
        scalar = coefficient * profile * y
        profiles.append((scalar.subs(substitutions), sp.diff(scalar, r).subs(substitutions)))
        values.append(complex(sp.N((coefficient * profile).subs(substitutions), 30)))
        derivatives.append(complex(sp.N((coefficient * derivative).subs(substitutions), 30)))
    return profiles, values, derivatives


def compact_lapse_profiles(
    profiles: list[tuple[sp.Expr, sp.Expr]],
    x: sp.Symbol,
    cache: dict[str, sp.Expr],
) -> tuple[list[tuple[sp.Expr, sp.Expr]], int]:
    compacted = []
    checks = 0
    for scalar, scalar_r in profiles:
        compact_scalar = compact_component(scalar, x)
        compact_scalar_r = compact_component(scalar_r, x)
        if clean_expr(scalar - compact_scalar, x, cache) != 0:
            raise AssertionError("compacted lapse scalar changed")
        checks += 1
        if clean_expr(scalar_r - compact_scalar_r, x, cache) != 0:
            raise AssertionError("compacted lapse radial derivative changed")
        checks += 1
        compacted.append((compact_scalar, compact_scalar_r))
    return compacted, checks


def compact_profile_spatial_families(
    families: list[tuple[list[dict[str, Any]], list[complex], list[complex]]],
    x: sp.Symbol,
    side: str,
) -> tuple[list[tuple[list[dict[str, Any]], list[complex], list[complex]]], dict[str, Any]]:
    tick = substage("compact_profile_family", "start", None, side=side)
    equality_cache: dict[str, sp.Expr] = {}
    equality_checks = 0
    compacted = []
    for family, (jets, values, derivatives) in enumerate(families):
        family_tick = substage("compact_profile_family", "start", None, side=side, family=family)
        compact_jets = []
        for index, jet in enumerate(jets):
            compacted_jet = compact_jet(jet, x)
            equality_checks += assert_compacted_jet_equal(jet, compacted_jet, x, equality_cache)
            assert_jet_shape(compacted_jet, f"{side}[{family}][{index}]")
            compact_jets.append(compacted_jet)
        compacted.append((compact_jets, values, derivatives))
        substage(
            "compact_profile_family",
            "done",
            family_tick,
            side=side,
            family=family,
            jets=len(jets),
        )
    elapsed = now() - tick
    substage(
        "compact_profile_family",
        "done",
        tick,
        side=side,
        equality_checks=equality_checks,
        families=len(families),
        elapsed_seconds=elapsed,
    )
    return compacted, {
        "equality_checks": equality_checks,
        "families": len(families),
        "elapsed_seconds": elapsed,
    }


def compact_profile_lapse_families(
    families: list[tuple[list[tuple[sp.Expr, sp.Expr]], list[complex], list[complex]]],
    x: sp.Symbol,
) -> tuple[list[tuple[list[tuple[sp.Expr, sp.Expr]], list[complex], list[complex]]], dict[str, Any]]:
    tick = substage("compact_profile_family", "start", None, side="lapse")
    equality_cache: dict[str, sp.Expr] = {}
    equality_checks = 0
    compacted = []
    for family, (profiles, values, derivatives) in enumerate(families):
        family_tick = substage(
            "compact_profile_family", "start", None, side="lapse", family=family
        )
        compact_profiles, family_checks = compact_lapse_profiles(
            profiles, x, equality_cache
        )
        equality_checks += family_checks
        compacted.append((compact_profiles, values, derivatives))
        substage(
            "compact_profile_family",
            "done",
            family_tick,
            side="lapse",
            family=family,
            profiles=len(profiles),
        )
    elapsed = now() - tick
    substage(
        "compact_profile_family",
        "done",
        tick,
        side="lapse",
        equality_checks=equality_checks,
        families=len(families),
        elapsed_seconds=elapsed,
    )
    return compacted, {
        "equality_checks": equality_checks,
        "families": len(families),
        "elapsed_seconds": elapsed,
    }


def fsum_complex(values: list[complex]) -> complex:
    return complex(
        math.fsum(value.real for value in values),
        math.fsum(value.imag for value in values),
    )


def integrate_spatial_pair_sum(
    evaluator: FixedChartEvaluator,
    fixed_data: dict[str, Any],
    left_jets: list[dict[str, Any]],
    right_jets: list[dict[str, Any]],
    density: Callable[[dict[str, Any], dict[str, Any], dict[str, Any]], sp.Expr],
    block: str,
    lf: int,
    rf: int,
) -> tuple[complex, dict[str, Any]]:
    expected_terms = len(left_jets) * len(right_jets)
    tick = substage("direct_pair_sum", "start", None, block=block, lf=lf, rf=rf,
                    expected_terms=expected_terms)
    row_timings = []
    terms = []
    for left_index, left in enumerate(left_jets):
        row_started = substage("direct_pair_row", "start", None, block=block,
                               lf=lf, rf=rf, left_index=left_index,
                               expected_terms=len(right_jets))
        completed_terms = 0
        for right in right_jets:
            terms.append(evaluator.integrate(density(left, right, fixed_data), 96))
            completed_terms += 1
        row_timings.append({
            "left_index": left_index,
            "completed_terms": completed_terms,
            "elapsed_seconds": now() - row_started,
        })
        substage("direct_pair_row", "done", row_started, block=block, lf=lf, rf=rf,
                 left_index=left_index, completed_terms=completed_terms,
                 cumulative_terms=len(terms), expected_total_terms=expected_terms)
    substage(
        "direct_pair_sum",
        "done",
        tick,
        block=block,
        lf=lf,
        rf=rf,
        terms=len(terms),
    )
    return (
        fsum_complex(terms),
        {"block": block, "lf": lf, "rf": rf, "completed_terms": len(terms),
         "expected_terms": expected_terms, "elapsed_seconds": now() - tick, "rows": row_timings},
    )


def integrate_lapse_pair_sum(
    evaluator: FixedChartEvaluator,
    fixed_data: dict[str, Any],
    lapse_profiles: list[tuple[sp.Expr, sp.Expr]],
    right_jets: list[dict[str, Any]],
    density: Callable[..., sp.Expr],
    block: str,
    lf: int,
    rf: int,
) -> tuple[complex, dict[str, Any]]:
    expected_terms = len(lapse_profiles) * len(right_jets)
    tick = substage("direct_pair_sum", "start", None, block=block, lf=lf, rf=rf,
                    expected_terms=expected_terms)
    row_timings = []
    terms = []
    for left_index, (scalar, scalar_r) in enumerate(lapse_profiles):
        row_started = substage("direct_pair_row", "start", None, block=block,
                               lf=lf, rf=rf, left_index=left_index,
                               expected_terms=len(right_jets))
        completed_terms = 0
        for right in right_jets:
            if block == "lapse_density":
                terms.append(evaluator.integrate(density(scalar, scalar_r, right, fixed_data), 96))
            else:
                terms.append(evaluator.integrate(density(scalar, right, fixed_data), 96))
            completed_terms += 1
        row_timings.append({
            "left_index": left_index,
            "completed_terms": completed_terms,
            "elapsed_seconds": now() - row_started,
        })
        substage("direct_pair_row", "done", row_started, block=block, lf=lf, rf=rf,
                 left_index=left_index, completed_terms=completed_terms,
                 cumulative_terms=len(terms), expected_total_terms=expected_terms)
    substage(
        "direct_pair_sum",
        "done",
        tick,
        block=block,
        lf=lf,
        rf=rf,
        terms=len(terms),
    )
    return (
        fsum_complex(terms),
        {"block": block, "lf": lf, "rf": rf, "completed_terms": len(terms),
         "expected_terms": expected_terms, "elapsed_seconds": now() - tick, "rows": row_timings},
    )


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


def verify_legendre_norms(include_component_compaction: bool = False) -> dict[str, Any]:
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
    result = {
        "checked": checked,
        "regularity_guard_controls": regularity,
    }
    if include_component_compaction:
        result["component_compaction_controls"] = verify_component_compaction_controls(x)
    return result


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


def evaluate_case(
    m: int,
    r_value: sp.Rational,
    compact_actual_profiles: bool = False,
) -> dict[str, Any]:
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
    left_families = []
    for family in range(2):
        build_tick = substage("build_profile_jets", "start", None, side="left", family=family)
        left_families.append(
            fixed_profile_spatial_jets(
                data, left_coeffs, "left", family, gamma, substitutions
            )
        )
        substage("build_profile_jets", "done", build_tick, side="left", family=family)
    right_families = []
    for family in range(2):
        build_tick = substage("build_profile_jets", "start", None, side="right", family=family)
        right_families.append(
            fixed_profile_spatial_jets(
                data, right_coeffs, "right", family, gamma, substitutions
            )
        )
        substage("build_profile_jets", "done", build_tick, side="right", family=family)
    lapse_families = []
    for family in range(2):
        build_tick = substage("build_profile_jets", "start", None, side="lapse", family=family)
        lapse_families.append(
            fixed_profile_lapse_profiles(data, lapse_coeffs, family, substitutions)
        )
        substage("build_profile_jets", "done", build_tick, side="lapse", family=family)

    profile_compaction = {
        "performed": False,
        "m_gt_0_radical_compaction": False,
        "spatial_equality_checks": 0,
        "lapse_equality_checks": 0,
        "total_equality_checks": 0,
        "by_side": {},
    }
    if compact_actual_profiles:
        if m != 0:
            raise AssertionError("compact-profile-probe intentionally supports only m=0")
        compaction_tick = substage(
            "compact_actual_profile_oracle", "start", None, m=m, r=str(r_value)
        )
        left_families, left_compaction = compact_profile_spatial_families(
            left_families, x, "left"
        )
        right_families, right_compaction = compact_profile_spatial_families(
            right_families, x, "right"
        )
        lapse_families, lapse_compaction = compact_profile_lapse_families(
            lapse_families, x
        )
        spatial_checks = (
            left_compaction["equality_checks"] + right_compaction["equality_checks"]
        )
        lapse_checks = lapse_compaction["equality_checks"]
        profile_compaction = {
            "performed": True,
            "m_gt_0_radical_compaction": False,
            "spatial_equality_checks": spatial_checks,
            "lapse_equality_checks": lapse_checks,
            "total_equality_checks": spatial_checks + lapse_checks,
            "by_side": {
                "left": left_compaction,
                "right": right_compaction,
                "lapse": lapse_compaction,
            },
            "component_endpoint_poles_allowed": True,
            "final_scalar_regularity_guard": "existing FixedChartEvaluator/clean_expr path",
        }
        substage(
            "compact_actual_profile_oracle",
            "done",
            compaction_tick,
            m=m,
            r=str(r_value),
            equality_checks=profile_compaction["total_equality_checks"],
        )

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
    direct_term_timings = {
        "potential_density": [],
        "potential_radial_flux": [],
        "lapse_density": [],
        "lapse_radial_flux": [],
    }
    for lf, (left_jets, left_v, left_d) in enumerate(left_families):
        for rf, (right_jets, right_v, right_d) in enumerate(right_families):
            direct_potential_density, timing = integrate_spatial_pair_sum(
                evaluator,
                fixed_data,
                left_jets,
                right_jets,
                weak.potential_density,
                "potential_density",
                lf,
                rf,
            )
            direct_term_timings["potential_density"].append(timing)
            direct_potential_flux, timing = integrate_spatial_pair_sum(
                evaluator,
                fixed_data,
                left_jets,
                right_jets,
                weak.potential_boundary_flux,
                "potential_radial_flux",
                lf,
                rf,
            )
            direct_term_timings["potential_radial_flux"].append(timing)
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

    for lf, (lapse_profiles, lapse_v, lapse_d) in enumerate(lapse_families):
        for rf, (right_jets, right_v, right_d) in enumerate(right_families):
            direct_lapse_density, timing = integrate_lapse_pair_sum(
                evaluator,
                fixed_data,
                lapse_profiles,
                right_jets,
                lapse_density_for_scalar,
                "lapse_density",
                lf,
                rf,
            )
            direct_term_timings["lapse_density"].append(timing)
            direct_lapse_flux, timing = integrate_lapse_pair_sum(
                evaluator,
                fixed_data,
                lapse_profiles,
                right_jets,
                lapse_flux_for_scalar,
                "lapse_radial_flux",
                lf,
                rf,
            )
            direct_term_timings["lapse_radial_flux"].append(timing)
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
    reconstruction["term_timing"] = direct_term_timings
    reconstruction["oracle_method"] = (
        "pairwise integration of compacted independently profile-weighted fields"
        if compact_actual_profiles
        else "pairwise integration of independently profile-weighted fields"
    )
    reconstruction["actual_profile_compaction"] = profile_compaction
    reconstruction["raw_slot_export_unchanged"] = compact_actual_profiles
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
    if m == 0:
        realfield = realfield_oracle(evaluator, fixed_data, [], [], [], family_matrices, m)
    else:
        realfield = realfield_oracle(
            evaluator,
            fixed_data,
            [
                fixed_dense_spatial_jet(
                    data, left_coeffs, "left", family, gamma, substitutions
                )[0]
                for family in range(2)
            ],
            [
                fixed_dense_spatial_jet(
                    data, right_coeffs, "right", family, gamma, substitutions
                )[0]
                for family in range(2)
            ],
            [
                fixed_dense_lapse_profile(data, lapse_coeffs, family, substitutions)[:2]
                for family in range(2)
            ],
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
        "actual_profile_compaction": profile_compaction,
        "structural_zero_flux_slots": flux_structural_zeros,
        "elapsed_seconds": now() - started,
    }


def timed_quadrature_ladder(
    evaluator: TimedChartEvaluator,
    value: sp.Expr,
) -> tuple[dict[str, Any], float, str | None, dict[str, float]]:
    orders: dict[int, complex] = {}
    timings: dict[str, float] = {}
    exact_nonzero = None
    for order in QUADRATURE_ORDERS:
        orders[order], order_timings, real, imag = evaluator.timed_integrate(value, order)
        if exact_nonzero is None:
            if real != 0:
                exact_nonzero = str(real)
            elif imag != 0:
                exact_nonzero = f"I*({imag})"
        timings[str(order)] = sum(order_timings.values())
        for key, elapsed in order_timings.items():
            timings[f"{order}_{key}"] = elapsed
    err24_48 = abs(orders[48] - orders[24])
    err48_96 = abs(orders[96] - orders[48])
    scale = max(1.0, abs(orders[48]), abs(orders[96]))
    scaled = err48_96 / scale
    if scaled > REFINEMENT_RTOL:
        raise AssertionError(f"arithmetic probe quadrature refinement failed: {scaled}")
    payload = {
        "values": {str(order): cscalar(orders[order]) for order in QUADRATURE_ORDERS},
        "max_abs_24_48": float(err24_48),
        "max_abs_48_96": float(err48_96),
        "max_scaled_48_96": float(scaled),
        "timings": timings,
    }
    return payload, float(scaled), exact_nonzero, timings


def compact_profile_jets(
    left_raw: dict[str, Any],
    right_raw: list[dict[str, Any]],
    x: sp.Symbol,
) -> tuple[dict[str, Any], list[dict[str, Any]], int]:
    tick = substage("compact_profile_jets", "start", None)
    left_compact = compact_jet(left_raw, x)
    right_compact = [compact_jet(jet, x) for jet in right_raw]
    equality_cache: dict[str, sp.Expr] = {}
    equality_checks = assert_compacted_jet_equal(left_raw, left_compact, x, equality_cache)
    for raw, compacted in zip(right_raw, right_compact):
        equality_checks += assert_compacted_jet_equal(raw, compacted, x, equality_cache)
    substage(
        "compact_profile_jets",
        "done",
        tick,
        equality_checks=equality_checks,
        jets=1 + len(right_raw),
    )
    return left_compact, right_compact, equality_checks


def evaluate_arithmetic_pair(
    output_dir: Path | None,
    fixed_data: dict[str, Any],
    x: sp.Symbol,
    left_raw: dict[str, Any],
    left_compact: dict[str, Any],
    right_raw: dict[str, Any],
    right_compact: dict[str, Any],
    right_index: int,
    metadata: dict[str, Any],
) -> tuple[dict[str, Any], float]:
    pair_started = now()
    write_progress(output_dir, {
        "phase": "arithmetic",
        "step": right_index,
        "of": ARITHMETIC_PROBE_PAIRS,
        "unit": "pairs",
        "m": {"active_right_index": right_index},
    })

    tick = substage("arithmetic_density", "start", None, arm="raw", right_index=right_index)
    raw_density = weak.potential_density(left_raw, right_raw, fixed_data)
    raw_construct_seconds = now() - tick
    substage(
        "arithmetic_density",
        "done",
        tick,
        arm="raw",
        right_index=right_index,
    )

    tick = substage("arithmetic_density", "start", None, arm="compact", right_index=right_index)
    compact_density = weak.potential_density(left_compact, right_compact, fixed_data)
    compact_construct_seconds = now() - tick
    substage(
        "arithmetic_density",
        "done",
        tick,
        arm="compact",
        right_index=right_index,
    )

    reset_raw = reset_denominator_cache()
    raw_evaluator = TimedChartEvaluator(
        x,
        {},
        {"arm": "raw", "right_index": right_index},
    )
    raw_payload, raw_scaled, raw_nonzero, _ = timed_quadrature_ladder(
        raw_evaluator, raw_density
    )
    reset_compact = reset_denominator_cache()
    compact_evaluator = TimedChartEvaluator(
        x,
        {},
        {"arm": "compact", "right_index": right_index},
    )
    compact_payload, compact_scaled, compact_nonzero, _ = timed_quadrature_ladder(
        compact_evaluator, compact_density
    )

    exact_tick = substage("exact_density_difference", "start", None, right_index=right_index)
    exact_cache: dict[str, sp.Expr] = {}
    exact_difference = clean_expr(raw_density - compact_density, x, exact_cache)
    if exact_difference != 0:
        raise AssertionError(f"raw compact density mismatch at right_index={right_index}")
    substage("exact_density_difference", "done", exact_tick, right_index=right_index)

    for order in QUADRATURE_ORDERS:
        raw_value = complex(*raw_payload["values"][str(order)])
        compact_value = complex(*compact_payload["values"][str(order)])
        assert_close(
            raw_value,
            compact_value,
            f"raw compact arithmetic probe right_index={right_index} order={order}",
        )

    pair = {
        **metadata,
        "right_index": right_index,
        "raw": raw_payload,
        "compact": compact_payload,
        "raw_density_construct_seconds": raw_construct_seconds,
        "compact_density_construct_seconds": compact_construct_seconds,
        "exact_raw_minus_compact_zero": True,
        "raw_nonzero": raw_nonzero is not None,
        "compact_nonzero": compact_nonzero is not None,
        "raw_nonzero_witness": raw_nonzero,
        "compact_nonzero_witness": compact_nonzero,
        "denominator_cache_reset_before_raw": reset_raw,
        "denominator_cache_reset_before_compact": reset_compact,
        "diagnostic_timing_only": True,
        "diagnostic_only": True,
        "snapshot": False,
        "quadrature_orders": list(QUADRATURE_ORDERS),
        "refinement_rtol": REFINEMENT_RTOL,
        "raw_compact_rtol": REFINEMENT_RTOL,
        "angular_measure": "dx/2",
        "elapsed_seconds": now() - pair_started,
    }
    if output_dir is not None:
        safe_write_json(output_dir / f"profile-arithmetic-pair-r{right_index}.json", pair)
    return pair, max(raw_scaled, compact_scaled)


def verify_profile_arithmetic_probe(output_dir: Path | None) -> dict[str, Any]:
    started = now()
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / PROGRESS_NAME).write_text("")

    write_progress(output_dir, {
        "phase": "setup",
        "step": 1,
        "of": 1,
        "unit": "environment",
        "plan": arithmetic_runtime_plan(),
        "mode": "profile-arithmetic-probe",
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
        ("legendre_norms", lambda: verify_legendre_norms(True)),
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

    m = 0
    r_value = sp.Rational(3, 2)
    lf = rf = 0
    left_index = 10
    pairs = []
    max_scaled = 0.0
    nonzero_witnesses = []
    write_progress(output_dir, {
        "phase": "arithmetic",
        "step": 0,
        "of": ARITHMETIC_PROBE_PAIRS,
        "unit": "pairs",
        "m": {"completed_arithmetic_pairs": 0, "max_scaled_refinement_error": max_scaled},
    })
    data = kinematics.sector(L_CUTOFF, m)
    if len(data["spatial_labels"]) != ARITHMETIC_PROBE_PAIRS:
        raise AssertionError(f"m=0 probe expected 12 spatial labels, got {len(data['spatial_labels'])}")
    if len(data["spatial_tangents"]) != ARITHMETIC_PROBE_PAIRS:
        raise AssertionError(
            f"m=0 probe expected 12 spatial tangents, got {len(data['spatial_tangents'])}"
        )
    x = data["x"]
    substitutions = {data["r"]: r_value, data["epsilon"]: EPSILON_VALUE}
    fixed_data = fixed_metric_data(data, substitutions)
    gamma = weak.connection(data)
    left_coeff = sp.Rational(left_index + 1, 13)
    if data["spatial_labels"][left_index] != ("C_polar", 2, 0):
        raise AssertionError(
            f"probe left index changed: {data['spatial_labels'][left_index]}"
        )
    right_coeffs = [sp.Rational(2 * index + 3, 17) for index in range(len(data["spatial_labels"]))]

    build_tick = substage("arithmetic_build_jets", "start", None)
    left_raw = fixed_single_profile_spatial_jet(
        data, left_coeff, "left", lf, left_index, gamma, substitutions
    )
    right_raw = [
        fixed_single_profile_spatial_jet(
            data, coefficient, "right", rf, right_index, gamma, substitutions
        )
        for right_index, coefficient in enumerate(right_coeffs)
    ]
    substage("arithmetic_build_jets", "done", build_tick)
    assert_jet_shape(left_raw, "left_raw")
    for right_index, jet in enumerate(right_raw):
        assert_jet_shape(jet, f"right_raw[{right_index}]")
    left_compact, right_compact, equality_checks = compact_profile_jets(left_raw, right_raw, x)
    assert_jet_shape(left_compact, "left_compact")
    for right_index, jet in enumerate(right_compact):
        assert_jet_shape(jet, f"right_compact[{right_index}]")

    for right_index, right in enumerate(right_raw):
        pair_metadata = {
            "m": m,
            "r": str(r_value),
            "epsilon": str(EPSILON_VALUE),
            "lf": lf,
            "rf": rf,
            "left_index": left_index,
            "left_label": label_tuple(data["spatial_labels"][left_index]),
            "left_coefficient": str(left_coeff),
            "left_profile": str(polynomial_profile(data["r"], left_index, "left", lf)),
            "right_label": label_tuple(data["spatial_labels"][right_index]),
            "right_coefficient": str(right_coeffs[right_index]),
            "right_profile": str(polynomial_profile(data["r"], right_index, "right", rf)),
        }
        pair, scaled = evaluate_arithmetic_pair(
            output_dir,
            fixed_data,
            x,
            left_raw,
            left_compact,
            right,
            right_compact[right_index],
            right_index,
            pair_metadata,
        )
        pairs.append(pair)
        max_scaled = max(max_scaled, scaled)
        if pair["raw_nonzero"]:
            nonzero_witnesses.append(right_index)
        write_checkpoint(output_dir, started, f"arithmetic_right_{right_index}", pair)
        write_progress(output_dir, {
            "phase": "arithmetic",
            "step": right_index + 1,
            "of": ARITHMETIC_PROBE_PAIRS,
            "unit": "pairs",
            "m": {
                "completed_arithmetic_pairs": right_index + 1,
                "max_scaled_refinement_error": max_scaled,
            },
        })
    if not nonzero_witnesses:
        raise AssertionError("profile arithmetic probe saw only zero raw densities")

    result = {
        "status": "passed",
        "mode": "profile-arithmetic-probe",
        "scope": "single-row raw versus component-compacted arithmetic diagnostic",
        "diagnostic_only": True,
        "row_arithmetic_comparison_only": True,
        "full_slots_built": False,
        "production_path_compaction": False,
        "profile_block_completed": False,
        "snapshot_file_written": False,
        "six_snapshot_export_performed": False,
        "production_compaction_integrated": False,
        "full_snapshot_performed": False,
        "angular_snapshot": False,
        "m_values": [m],
        "r_values": [str(r_value)],
        "m": m,
        "r": str(r_value),
        "epsilon": str(EPSILON_VALUE),
        "lf": lf,
        "rf": rf,
        "left_index": left_index,
        "right_terms": ARITHMETIC_PROBE_PAIRS,
        "right_indices": list(range(ARITHMETIC_PROBE_PAIRS)),
        "left_label": label_tuple(data["spatial_labels"][left_index]),
        "right_labels": [label_tuple(label) for label in data["spatial_labels"]],
        "left_coefficient": str(left_coeff),
        "right_coefficients": [str(coefficient) for coefficient in right_coeffs],
        "left_profile": str(polynomial_profile(data["r"], left_index, "left", lf)),
        "right_profiles": [
            str(polynomial_profile(data["r"], right_index, "right", rf))
            for right_index in range(ARITHMETIC_PROBE_PAIRS)
        ],
        "component_equality_checks": equality_checks,
        "component_endpoint_poles_allowed": True,
        "final_scalar_regularity_guard": "existing FixedChartEvaluator/clean_expr path",
        "nonzero_witness_right_indices": nonzero_witnesses,
        "pairs": pairs,
        "max_scaled_refinement_error": max_scaled,
        "controls": controls,
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
        safe_write_json(output_dir / ARITHMETIC_PROBE_RESULT_NAME, result)
        write_progress(output_dir, {
            "phase": "export",
            "step": 1,
            "of": 1,
            "unit": "checker result",
            "m": {
                "profile_arithmetic_probe_passed": 1,
                "elapsed_seconds": now() - started,
            },
        })
    print("PASS profile-arithmetic-probe", flush=True)
    return result


def mode_cases(mode: str) -> list[tuple[int, sp.Rational]]:
    if mode == "full":
        return [(m, r_value) for m in M_VALUES for r_value in R_VALUES]
    if mode in ("direct-profile-probe", "compact-profile-probe"):
        return [(0, sp.Rational(3, 2))]
    raise ValueError(f"unsupported mode: {mode}")


def verify_export(output_dir: Path | None, mode: str = "full") -> dict[str, Any]:
    started = now()
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / PROGRESS_NAME).write_text("")

    cases_to_run = mode_cases(mode)
    if mode == "full":
        result_name = RESULT_NAME
        pass_metric = "angular_export_passed"
    elif mode == "compact-profile-probe":
        result_name = COMPACT_PROFILE_PROBE_RESULT_NAME
        pass_metric = "compact_profile_probe_passed"
    else:
        result_name = PROBE_RESULT_NAME
        pass_metric = "direct_profile_probe_passed"
    write_progress(output_dir, {
        "phase": "setup",
        "step": 1,
        "of": 1,
        "unit": "environment",
        "plan": runtime_plan(len(cases_to_run)),
        "mode": mode,
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
        ("legendre_norms", lambda: verify_legendre_norms(mode == "compact-profile-probe")),
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
        "of": len(cases_to_run),
        "unit": "snapshots",
        "m": {
            "completed_snapshots": completed,
            "max_scaled_refinement_error": max_scaled_refinement_error,
        },
    })
    for m, r_value in cases_to_run:
        stage("angular_case", m=m, r=str(r_value))
        compact_profiles = mode == "compact-profile-probe"
        case = evaluate_case(m, r_value, compact_actual_profiles=compact_profiles)
        if output_dir is not None:
            snapshot_prefix = (
                "compact-profile-snapshot" if compact_profiles else "angular-snapshot"
            )
            safe_write_json(
                output_dir / f"{snapshot_prefix}-m{m}-r{str(r_value).replace('/', '_')}.json",
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
            "of": len(cases_to_run),
            "unit": "snapshots",
            "m": {
                "completed_snapshots": completed,
                "max_scaled_refinement_error": max_scaled_refinement_error,
            },
        })

    result = {
        "status": "passed",
        "scope": (
            "fixed-radius angular weak-slot diagnostic snapshots"
            if mode == "full"
            else "first-case compact actual-profile diagnostic snapshot"
            if mode == "compact-profile-probe"
            else "first-case direct-profile diagnostic snapshot"
        ),
        "kind": "angular_snapshot",
        "mode": mode,
        "diagnostic_snapshots": True,
        "L": L_CUTOFF,
        "m_values": list(dict.fromkeys(m for m, _ in cases_to_run)),
        "r_values": list(dict.fromkeys(str(r_value) for _, r_value in cases_to_run)),
        "first_case_only": mode in ("direct-profile-probe", "compact-profile-probe"),
        "full_mode_unchanged": True,
        "raw_slot_export_unchanged": True,
        "m_gt_0_radical_compaction": False,
        "compact_actual_profile_oracle": mode == "compact-profile-probe",
        "full_six_snapshots_performed": mode == "full",
        "six_snapshot_export_performed": mode == "full",
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
        "compact_profile_spatial_equality_checks": sum(
            case["actual_profile_compaction"]["spatial_equality_checks"]
            for case in cases
        ),
        "compact_profile_lapse_equality_checks": sum(
            case["actual_profile_compaction"]["lapse_equality_checks"]
            for case in cases
        ),
        "compact_profile_total_equality_checks": sum(
            case["actual_profile_compaction"]["total_equality_checks"]
            for case in cases
        ),
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
        safe_write_json(output_dir / result_name, result)
        write_progress(output_dir, {
            "phase": "export",
            "step": 1,
            "of": 1,
            "unit": "checker result",
            "m": {
                pass_metric: 1,
                "elapsed_seconds": now() - started,
            },
        })
    print(f"PASS {mode} angular weak-slot diagnostic snapshots", flush=True)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--mode",
        choices=(
            "full",
            "direct-profile-probe",
            "profile-arithmetic-probe",
            "compact-profile-probe",
        ),
        default="full",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.mode == "profile-arithmetic-probe":
        verify_profile_arithmetic_probe(args.output_dir)
    else:
        verify_export(args.output_dir, args.mode)


if __name__ == "__main__":
    main()

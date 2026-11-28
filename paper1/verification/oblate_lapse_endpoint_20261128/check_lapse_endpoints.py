#!/usr/bin/env python3
"""Check the oblate lapse endpoint Green identities.

Prepared standalone SymPy checker. Do not execute without the separately
bounded runner that owns wall and memory limits.
"""

from __future__ import annotations

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time

import sympy as s


EPSILON = s.Rational(1, 100)
SHELL = (s.Integer(1), s.Integer(2))
REVIEWED_REPORT = (
    "paper1/reports/oblate-common-domain-preparation-2026-11-28/README.md"
)
REVIEWED_REPORT_COMMIT = "ae32314d27e20ac682c5ed1b71b4d7f11bee7133"
REVIEWED_REPORT_SHA256 = (
    "fa659f17f3ba9a8dce1ff0d2201cb5418fa2a10faeb64251ebcc20f690d3bc7d"
)
SOURCE_OF_DERIVATION_COMMIT = "773ccc4e3e8a1a89e06086163042fbf45e171c76"


def emit_raw(event, **payload):
    print(json.dumps({"event": event, **payload}, sort_keys=True), flush=True)


class Progress:
    def __init__(self):
        self._phase_stack = []

    def _numeric_metrics(self, metrics):
        return {
            key: int(value) if isinstance(value, bool) else value
            for key, value in metrics.items()
            if isinstance(value, (bool, int, float))
        }

    def report(self, **payload):
        print(
            "::orx:: " + json.dumps(payload, sort_keys=True, allow_nan=False),
            flush=True,
        )

    @contextmanager
    def stage(self, name, unit="check"):
        start = time.monotonic()
        self.report(phase=name, step=0, of=1, unit=unit)
        emit_raw("STAGE_BEGIN", stage=name)
        self._phase_stack.append(name)
        try:
            yield
        except Exception as exc:
            emit_raw("STAGE_FAIL", stage=name, error=type(exc).__name__)
            raise
        finally:
            self._phase_stack.pop()
        elapsed = time.monotonic() - start
        self.report(
            phase=name,
            step=1,
            of=1,
            unit=unit,
            m={"elapsed_seconds": elapsed},
        )
        emit_raw("STAGE_END", stage=name, seconds=elapsed)

    def metric(self, name, value, **extra):
        phase = self._phase_stack[-1] if self._phase_stack else "unscoped_metric"
        metrics = self._numeric_metrics({name: value, **extra})
        payload = {"phase": phase, "step": 0, "of": 1, "unit": "check"}
        if metrics:
            payload["m"] = metrics
        self.report(**payload)
        emit_raw("PROGRESS", metric=name, value=value, **extra)


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def canonical(expr):
    return s.cancel(s.together(expr))


def assert_zero(expr, label):
    require(canonical(expr) == s.S.Zero, label)


def assert_nonzero(expr, label):
    require(canonical(expr) != s.S.Zero, label)


def christoffel(metric, coords):
    inv = metric.inv()
    n = len(coords)
    return [
        [
            [
                canonical(
                    sum(
                        inv[a, k]
                        * (
                            s.diff(metric[k, c], coords[b])
                            + s.diff(metric[k, b], coords[c])
                            - s.diff(metric[b, c], coords[k])
                        )
                        for k in range(n)
                    )
                    / 2
                )
                for c in range(n)
            ]
            for b in range(n)
        ]
        for a in range(n)
    ]


def ricci_tensor(Gamma, coords):
    n = len(coords)
    ricci = s.zeros(n)
    for i in range(n):
        for j in range(n):
            derivative = sum(
                s.diff(Gamma[k][i][j], coords[k])
                - s.diff(Gamma[k][i][k], coords[j])
                for k in range(n)
            )
            products = sum(
                Gamma[k][i][j] * Gamma[l][k][l]
                - Gamma[k][i][l] * Gamma[l][k][j]
                for k in range(n)
                for l in range(n)
            )
            ricci[i, j] = canonical(derivative + products)
    return ricci


def varied_christoffel(metric, perturbation, coords, Gamma, inv):
    d_inv = -inv * perturbation * inv
    n = len(coords)
    return [
        [
            [
                canonical(
                    sum(
                        d_inv[a, k]
                        * (
                            s.diff(metric[k, c], coords[b])
                            + s.diff(metric[k, b], coords[c])
                            - s.diff(metric[b, c], coords[k])
                        )
                        + inv[a, k]
                        * (
                            s.diff(perturbation[k, c], coords[b])
                            + s.diff(perturbation[k, b], coords[c])
                            - s.diff(perturbation[b, c], coords[k])
                        )
                        for k in range(n)
                    )
                    / 2
                )
                for c in range(n)
            ]
            for b in range(n)
        ]
        for a in range(n)
    ]


def varied_ricci_tensor(Gamma, dGamma, coords):
    n = len(coords)
    d_ricci = s.zeros(n)
    for i in range(n):
        for j in range(n):
            derivative = sum(
                s.diff(dGamma[k][i][j], coords[k])
                - s.diff(dGamma[k][i][k], coords[j])
                for k in range(n)
            )
            products = sum(
                dGamma[k][i][j] * Gamma[l][k][l]
                + Gamma[k][i][j] * dGamma[l][k][l]
                - dGamma[k][i][l] * Gamma[l][k][j]
                - Gamma[k][i][l] * dGamma[l][k][j]
                for k in range(n)
                for l in range(n)
            )
            d_ricci[i, j] = canonical(derivative + products)
    return d_ricci


def hessian(scalar, Gamma, coords):
    n = len(coords)
    return s.Matrix(
        n,
        n,
        lambda i, j: canonical(
            s.diff(scalar, coords[i], coords[j])
            - sum(Gamma[k][i][j] * s.diff(scalar, coords[k]) for k in range(n))
        ),
    )


def trace(inv, tensor):
    n = inv.rows
    return canonical(sum(inv[i, j] * tensor[i, j] for i in range(n) for j in range(n)))


def raise_both(inv, tensor):
    return inv * tensor * inv


def dewit_raw_density(tensor_a, tensor_b, inv, sqrt_det):
    raised_b = raise_both(inv, tensor_b)
    contraction = sum(
        tensor_a[i, j] * raised_b[i, j]
        for i in range(inv.rows)
        for j in range(inv.cols)
    )
    return canonical(sqrt_det * (contraction - trace(inv, tensor_a) * trace(inv, tensor_b)))


def divergence_second_index(P, Gamma, coords):
    n = len(coords)
    return [
        canonical(
            sum(s.diff(P[i, j], coords[j]) for j in range(n))
            + sum(
                Gamma[i][j][k] * P[k, j] + Gamma[j][j][k] * P[i, k]
                for j in range(n) for k in range(n)
            )
        )
        for i in range(n)
    ]


def coordinate_divergence(flux, coords):
    return canonical(sum(s.diff(flux[i], coords[i]) for i in range(len(coords))))


def radial_dirichlet_subs(F, G, r, F1, F2, G1, G2):
    return {
        F: 0,
        s.diff(F, r): F1,
        s.diff(F, r, 2): F2,
        s.diff(F, r, 3): s.Symbol("F3"),
        G: 0,
        s.diff(G, r): G1,
        s.diff(G, r, 2): G2,
        s.diff(G, r, 3): s.Symbol("G3"),
    }


def apply_endpoint_jet(expr, replacements):
    return canonical(expr.subs(replacements, simultaneous=True))


def atomic_json_new(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(payload, indent=2, sort_keys=True).encode("utf-8") + b"\n"
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o644)
    try:
        with os.fdopen(fd, "wb") as handle:
            handle.write(data)
            handle.flush()
            os.fsync(handle.fileno())
        os.link(tmp, path)
    finally:
        try:
            os.unlink(tmp)
        except FileNotFoundError:
            pass


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def read_config(path):
    if path is None:
        return None, None
    data = Path(path).read_bytes()
    return json.loads(data.decode("utf-8")), sha256_bytes(data)


def source_sha256():
    return sha256_bytes(Path(__file__).read_bytes())


def sign_at(expr, r, x, r_value):
    return canonical(expr.subs({r: r_value, x: s.S.Zero}))


def pass_case(cases, name, **payload):
    cases.append({"name": name, "passed": True, **payload})


def run_checks(config_sha256=None, checker_sha256=None):
    progress = Progress()
    controls = []
    cases = []
    r, x, phi = s.symbols("r x phi", real=True)
    F = s.Function("F")(r)
    G = s.Function("G")(r)
    F1, F2, G1, G2 = s.symbols("F1 F2 G1 G2", real=True)

    with progress.stage("metric_and_variation"):
        d = 1 - x * x
        rho = r * r + EPSILON
        Sigma = r * r + EPSILON * x * x
        coords = [r, x, phi]
        metric = s.diag(Sigma / rho, Sigma / d, rho * d)
        inv = metric.inv()
        sqrt_det = Sigma
        Y = s.legendre(2, x)
        U = s.zeros(3)
        U[1, 1] = Y / d
        U[2, 2] = Y * d
        f = F * Y
        g = G * Y
        metric_det_residual = canonical(metric.det() - Sigma**2)
        assert_zero(metric_det_residual, "oblate metric determinant equals Sigma^2")
        pass_case(cases, "metric_determinant_sigma_squared",
                  residual=str(metric_det_residual))
        progress.metric("matrix_dimension", 3)
        progress.metric("epsilon_numerator", int(EPSILON.p))
        progress.metric("epsilon_denominator", int(EPSILON.q))

    with progress.stage("christoffel_and_linearized_ricci"):
        Gamma = christoffel(metric, coords)
        ricci = ricci_tensor(Gamma, coords)
        for i in range(3):
            for j in range(3):
                assert_zero(ricci[i, j], f"background Ricci[{i},{j}]")
        pass_case(cases, "background_ricci_zero", residual="0")
        dGamma = varied_christoffel(metric, U, coords, Gamma, inv)
        dRicci = varied_ricci_tensor(Gamma, dGamma, coords)
        d_inv = -inv * U * inv
        delta_R = canonical(
            sum(inv[i, j] * dRicci[i, j] + d_inv[i, j] * ricci[i, j]
                for i in range(3) for j in range(3))
        )
        controls.append("linearized Ricci scalar constructed from varied Christoffels")

    with progress.stage("scalar_curvature_green_identity"):
        Hf = hessian(f, Gamma, coords)
        trU = trace(inv, U)
        U_up = raise_both(inv, U)
        P = s.Matrix(3, 3, lambda i, j: canonical(U_up[i, j] - inv[i, j] * trU))
        divP = divergence_second_index(P, Gamma, coords)
        flux1 = [
            canonical(sqrt_det * (f * divP[i] - sum(P[i, j] * s.diff(f, coords[j])
                                                       for j in range(3))))
            for i in range(3)
        ]
        density1 = canonical(sqrt_det * f * delta_R - dewit_raw_density(Hf, U, inv, sqrt_det))
        scalar_residual = canonical(density1 - coordinate_divergence(flux1, coords))
        assert_zero(scalar_residual,
                    "scalar curvature Green density equals coordinate flux divergence")
        pass_case(cases, "scalar_curvature_green_density", residual=str(scalar_residual))
        progress.metric("scalar_flux_components", len(flux1))
        controls.append("strong-minus-adjoint density equals radial plus angular flux divergence")

    with progress.stage("hessian_flux_identity"):
        Hg = hessian(g, Gamma, coords)
        lap_g = trace(inv, Hg)
        Hg_up = raise_both(inv, Hg)
        Pg = s.Matrix(3, 3, lambda i, j: canonical(Hg_up[i, j] - inv[i, j] * lap_g))
        divPg = divergence_second_index(Pg, Gamma, coords)
        for component in divPg:
            assert_zero(component, "Hessian P tensor is divergence-free")
        flux2 = [
            canonical(sqrt_det * sum(Pg[i, j] * s.diff(f, coords[j])
                                     for j in range(3)))
            for i in range(3)
        ]
        density2 = dewit_raw_density(Hf, Hg, inv, sqrt_det)
        hessian_residual = canonical(density2 - coordinate_divergence(flux2, coords))
        assert_zero(hessian_residual,
                    "two-Hessian DeWitt density equals coordinate flux divergence")
        pass_case(cases, "two_hessian_density", residual=str(hessian_residual))
        controls.append("two-Hessian density equals radial plus angular flux divergence")

    with progress.stage("endpoint_flux_reduction"):
        endpoint_subs = radial_dirichlet_subs(F, G, r, F1, F2, G1, G2)
        first_radial = apply_endpoint_jet(flux1[0], endpoint_subs)
        hessian_radial = apply_endpoint_jet(flux2[0], endpoint_subs)
        expected_first = canonical(rho * F1 * Y**2 * (1 / Sigma + 1 / rho))
        expected_hessian = canonical(
            -rho**2 * r * (1 / rho + 1 / Sigma) * F1 * G1 * Y**2 / Sigma
        )
        first_radial_residual = canonical(first_radial - expected_first)
        hessian_radial_residual = canonical(hessian_radial - expected_hessian)
        assert_zero(first_radial_residual, "Dirichlet scalar radial flux")
        assert_zero(hessian_radial_residual, "Dirichlet Hessian radial flux")
        pass_case(cases, "scalar_dirichlet_radial_flux",
                  residual=str(first_radial_residual))
        pass_case(cases, "hessian_dirichlet_radial_flux",
                  residual=str(hessian_radial_residual))
        for label, component in [("scalar", flux1[1]), ("hessian", flux2[1])]:
            south = canonical(s.limit(canonical(component), x, -1))
            north = canonical(s.limit(canonical(component), x, 1))
            assert_zero(south, f"{label} south pole flux")
            assert_zero(north, f"{label} north pole flux")
            pass_case(cases, f"{label}_angular_pole_fluxes",
                      south_residual=str(south), north_residual=str(north))
        scalar_clamped = canonical(first_radial.subs(F1, 0))
        hessian_clamped = canonical(hessian_radial.subs({F1: 0, G1: 0}))
        scalar_wrong_sign = canonical(first_radial + expected_first)
        hessian_wrong_factor = canonical(hessian_radial - expected_hessian / 2)
        assert_zero(scalar_clamped, "clamped scalar flux vanishes")
        assert_zero(hessian_clamped, "clamped Hessian flux vanishes")
        assert_nonzero(scalar_wrong_sign, "wrong scalar sign is detected")
        assert_nonzero(hessian_wrong_factor,
                       "wrong Hessian factor is detected")
        pass_case(cases, "clamped_fluxes",
                  scalar_residual=str(scalar_clamped),
                  hessian_residual=str(hessian_clamped))
        pass_case(cases, "sign_and_factor_defect_controls",
                  scalar_wrong_sign=str(scalar_wrong_sign),
                  hessian_wrong_factor=str(hessian_wrong_factor))
        progress.metric("endpoint_count", 2)
        controls.append("Dirichlet radial fluxes match reviewed endpoint predictions")
        controls.append("angular pole fluxes vanish")
        controls.append("clamped endpoint slopes kill both fluxes")
        controls.append("sign and factor defects stay nonzero")

    with progress.stage("endpoint_sign_controls"):
        first_outer = sign_at(expected_first.subs(F1, 1), r, x, SHELL[1])
        first_inner = sign_at(expected_first.subs(F1, 1), r, x, SHELL[0])
        hessian_outer = sign_at(expected_hessian.subs({F1: 1, G1: 1}), r, x, SHELL[1])
        hessian_inner = sign_at(expected_hessian.subs({F1: 1, G1: 1}), r, x, SHELL[0])
        require(first_outer > 0, "outer scalar endpoint must be positive")
        require(first_inner > 0, "coordinate inner scalar endpoint must be positive")
        require(-first_inner < 0, "oriented inner scalar endpoint must be negative")
        require(hessian_outer < 0, "outer Hessian endpoint must be negative")
        require(-hessian_inner > 0, "oriented inner Hessian endpoint must be positive")
        pass_case(cases, "endpoint_orientation_signs",
                  scalar_outer=str(first_outer),
                  scalar_oriented_inner=str(canonical(-first_inner)),
                  hessian_outer=str(hessian_outer),
                  hessian_oriented_inner=str(canonical(-hessian_inner)))
        controls.append("outer-only and inner-only endpoint orientations have expected signs")

    with progress.stage("result_assembly"):
        result = {
            "state": "prepared lapse endpoint identities verified",
            "reviewed_report": REVIEWED_REPORT,
            "reviewed_report_commit": REVIEWED_REPORT_COMMIT,
            "reviewed_report_sha256": REVIEWED_REPORT_SHA256,
            "source_of_derivation_commit": SOURCE_OF_DERIVATION_COMMIT,
            "reviewed_mathematical_definitions": {
                "report": REVIEWED_REPORT,
                "report_commit": REVIEWED_REPORT_COMMIT,
                "report_sha256": REVIEWED_REPORT_SHA256,
                "source_of_derivation_commit": SOURCE_OF_DERIVATION_COMMIT,
            },
            "input_sha256": {
                "config": config_sha256,
                "checker": checker_sha256,
            },
            "runtime": {
                "sys_version": sys.version,
                "python_version": platform.python_version(),
                "python_version_info": {
                    "major": sys.version_info.major,
                    "minor": sys.version_info.minor,
                    "micro": sys.version_info.micro,
                },
                "python_executable": sys.executable,
                "sympy_version": s.__version__,
            },
            "fixed_parameters": {
                "epsilon": str(EPSILON),
                "N": 3,
                "m": 0,
                "ell": 2,
                "radial_shell": [str(SHELL[0]), str(SHELL[1])],
                "Y": str(Y),
            },
            "normalization": {
                "scalar_curvature_identity_prefactor": "1/(8*pi)",
                "hessian_identity_prefactor": "1/(4*pi)",
            },
            "identities": {
                "metric_determinant_matches_sigma_squared": True,
                "linearized_ricci_from_varied_christoffels": True,
                "scalar_curvature_density_divergence": True,
                "two_hessian_density_divergence": True,
                "angular_pole_fluxes_vanish": True,
                "dirichlet_radial_fluxes_match": True,
            },
            "endpoint_fluxes_before_prefactors": {
                "scalar_dirichlet_radial": str(expected_first),
                "hessian_dirichlet_radial": str(expected_hessian),
            },
            "sign_samples": {
                "x": "0",
                "outer_r": str(SHELL[1]),
                "inner_r": str(SHELL[0]),
                "scalar_outer_coordinate": str(first_outer),
                "scalar_inner_coordinate": str(first_inner),
                "hessian_outer_coordinate": str(hessian_outer),
                "hessian_inner_coordinate": str(hessian_inner),
                "orientation_reverses_inner_endpoint": True,
            },
            "case_results": cases,
            "controls": [{"name": name, "passed": True} for name in controls],
            "scope": {
                "full_operator_domain_established": False,
                "evolution_invariance_established": False,
                "physical_null_boundary_selected": False,
                "final_Dirac_classification_established": False,
                "physical_DOF_count_established": False,
            },
        }
        progress.metric("control_count", len(controls))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    progress = Progress()
    with progress.stage("input_config"):
        _, config_sha256 = read_config(args.config)
        checker_sha256 = source_sha256()
        progress.metric("python_major", sys.version_info.major)
        progress.metric("python_minor", sys.version_info.minor)
        if config_sha256 is not None:
            progress.metric("config_bytes_read", 1)
    result = run_checks(config_sha256=config_sha256,
                        checker_sha256=checker_sha256)
    with progress.stage("atomic_final_export"):
        atomic_json_new(args.output, result)
    emit_raw("POST_EXPORT", state="atomic export written", output=str(args.output))


if __name__ == "__main__":
    main()

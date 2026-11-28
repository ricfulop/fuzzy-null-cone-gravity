#!/usr/bin/env python3
"""Prepared check; DO NOT execute without separate user approval.

This verifies geometry, coefficient projections and the full linear-Euler
map against independent identities. It does not provide an on-shell physical
seed or evaluate a nonlinear obstruction.
"""
from pathlib import Path
import hashlib
import importlib.metadata
import json
import platform
import subprocess
import time

import numpy as np
import sympy as s

import null_compatibility as n

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'reports/oblate-nonlinear-null-preparation-2026-11-27'


def assert_zero(values):
    assert all(s.cancel(s.together(v)) == s.S.Zero for v in values)


def lie_metric(metric, vector, coords):
    return s.Matrix(len(coords), len(coords), lambda i, j: sum(
        vector[k]*s.diff(metric[i, j], coords[k])
        + metric[k, j]*s.diff(vector[k], coords[i])
        + metric[i, k]*s.diff(vector[k], coords[j])
        for k in range(len(coords))))


def geometry(coords, epsilon):
    t, r, x, phi = coords
    eta = s.symbols("eta", real=True)
    f = s.Function("f")(t, r)
    d, rho, Sigma = 1-x*x, r*r+epsilon, r*r+epsilon*x*x
    g0 = n.background(coords, epsilon)
    transformed = s.Matrix([t, r, x, phi+eta*f*x])
    J = transformed.jacobian(coords)
    pullback = J.T*g0*J
    gamma = s.Matrix([[Sigma/d+eta*eta*rho*d*f*f, eta*rho*d*f],
                      [eta*rho*d*f, rho*d]])
    fields = n.Fields(-1, 0, Sigma/rho, gamma,
                      s.Matrix([0, eta*x*s.diff(f, r)]),
                      s.Matrix([0, eta*x*s.diff(f, t)]))
    assert_zero(n.reconstruct(fields)-pullback)
    # Separate embedding curvature from the second coordinate-orbit tangent.
    gamma_linear = s.Matrix([[Sigma/d, eta*rho*d*f],
                             [eta*rho*d*f, rho*d]])
    linear_path = n.Fields(-1, 0, Sigma/rho, gamma_linear, fields.V, fields.W)
    embedding_second = n.reconstruct(linear_path).diff(eta, 2).subs(eta, 0)
    defect = pullback.diff(eta, 2).subs(eta, 0)-embedding_second
    expected = s.zeros(4)
    expected[2, 2] = 2*rho*d*f*f
    assert_zero(defect-expected)
    # Independent Lie-derivative construction, including f_t and f_r terms.
    xi = [0, 0, 0, f*x]
    assert_zero(lie_metric(g0, xi, coords)-pullback.diff(eta).subs(eta, 0))
    assert_zero(lie_metric(lie_metric(g0, xi, coords), xi, coords)
                -pullback.diff(eta, 2).subs(eta, 0))
    trace = (d*defect[2, 2]+defect[3, 3]/d)/2
    coefficient = s.cancel(n.average(trace*s.legendre(4, x), x)*9)
    assert_zero([coefficient-s.Rational(8, 35)*rho*f*f])
    # Base-dependent rotation phi -> phi+eta f has unchanged angular metric.
    rotation = s.Matrix([t, r, x, phi+eta*f]).jacobian(coords)
    rotation_fields = n.Fields(-1, 0, Sigma/rho, g0[2:, 2:],
                              s.Matrix([0, eta*s.diff(f, r)]),
                              s.Matrix([0, eta*s.diff(f, t)]))
    assert_zero(n.reconstruct(rotation_fields)-rotation.T*g0*rotation)
    return pullback.diff(eta).subs(eta, 0), {
        "quadrupole_pullback_matches_full_reconstruction": True,
        "changing_tangent_defect_xx": str(expected[2, 2]),
        "trace_j4_coefficient": str(coefficient),
        "time_and_radial_profile_derivatives_retained": True,
        "rotation_control_exactly_retained": True,
        "finite_orbit_leak_is_a_dynamical_obstruction": False,
    }


def harmonic_map(coords):
    _, r, x, _ = coords
    w = s.Function("w")(r)
    q, z = n.coefficient_labels()
    coefficients = {(kind, ell): s.symbols(kind+"_"+str(ell), real=True)
                    for kind, ell in q+z}
    fields = n.coefficient_fields(coefficients, x)
    transformed = n.linear_angular_generator(fields, x, w)
    projected = n.project_coefficients(transformed, x)
    expected = n.expected_discarded_generator(coefficients, w)
    for label, value in projected.items():
        if label[1] > 2:
            assert_zero([value-expected.get(label, s.S.Zero)])
    # Completeness through degree four, independently reconstructing every field.
    reconstructed = n.coefficient_fields(projected, x)
    assert_zero(reconstructed.gamma-transformed.gamma)
    assert_zero(reconstructed.V-transformed.V)
    assert_zero(reconstructed.W-transformed.W)
    rotation = n.linear_angular_generator(fields, x, w, s.S.One)
    assert_zero(rotation.gamma)
    assert_zero(rotation.V)
    assert_zero(rotation.W)
    # A distinct numerical integration of the scalar orbit coefficient.
    integrand = s.lambdify(x, 9*(1-x*x)**2*s.legendre(4, x), "numpy")
    quadrature = []
    for count in [24, 48]:
        nodes, weights = np.polynomial.legendre.leggauss(count)
        value = float(weights@integrand(nodes)/2)
        error = abs(value-float(s.Rational(8, 35)))
        assert error < 1e-12
        quadrature.append({"nodes": count, "coefficient": value, "error": error})
    return {
        "q_labels": q, "z_labels": z,
        "discarded_generator_coefficients": [
            {"label": label, "expression": str(value)}
            for label, value in expected.items()],
        "independent_projection_and_reconstruction_agree": True,
        "rotation_generator_defect_zero": True,
        "scalar_orbit_quadrature": quadrature,
    }


def continuum_controls(coords, epsilon, gauge_perturbation):
    g0 = n.background(coords, epsilon)
    gauge_G = n.linearized_einstein(g0, gauge_perturbation, coords)
    assert_zero(gauge_G)
    # Independent nonzero off-shell control: conformal variation in four dimensions.
    t, r, _, _ = coords
    psi = s.Function("psi")(t, r)
    conformal_G = n.linearized_einstein(g0, 2*psi*g0, coords)
    Gamma = n.christoffel(g0, coords)
    Hess = s.Matrix(4, 4, lambda i, j:
        s.diff(psi, coords[i], coords[j])
        - sum(Gamma[k][i][j]*s.diff(psi, coords[k]) for k in range(4)))
    inv = g0.inv()
    box = sum(inv[i, j]*Hess[i, j] for i in range(4) for j in range(4))
    expected = -2*Hess+2*g0*box
    assert_zero(conformal_G-expected)
    assert any(s.cancel(value) != 0 for value in expected)
    # At the flat background E[test] = -Sigma deltaG^up:test, without extra r².
    q, z = n.coefficient_labels()
    probe = n.coefficient_fields({label: s.S.One for label in q+z}, coords[2])
    assert_zero([n.euler_pair_density(gauge_G, probe, coords, epsilon)])
    # Nonzero density control tests the sign, Sigma, index raising and both tr entries.
    b_test = n.coefficient_fields({("b", 0): s.S.One}, coords[2])
    density = n.euler_pair_density(conformal_G, b_test, coords, epsilon)
    expected_density = -4*(r*r+epsilon)*s.diff(psi, t, r)
    assert_zero([density-expected_density])
    return {
        "pure_continuum_gauge_full_linear_Einstein_zero": True,
        "conformal_off_shell_formula_verified": True,
        "conformal_b_test_euler_density": str(expected_density),
        "conformal_control_is_an_admissible_seed": False,
        "discarded_temporal_components_were_not_projected_away": True,
    }


def verify():
    start = time.monotonic()
    t, r, x, phi, epsilon = s.symbols("t r x phi epsilon", real=True)
    coords = [t, r, x, phi]
    gauge, geometry_checks = geometry(coords, epsilon)
    harmonic_checks = harmonic_map(coords)
    continuum_checks = continuum_controls(coords, epsilon, gauge)
    result = {
        "status": "passed map/control checks",
        "geometry": geometry_checks,
        "harmonic_map": harmonic_checks,
        "continuum_euler_controls": continuum_checks,
        "admissible_non_gauge_seed_supplied": False,
        "quadratic_source_evaluated_on_admissible_non_gauge_seed": False,
        "nonlinear_obstruction_established": False,
        "full_nonlinear_Dirac_classification_established": False,
    }
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/"map-checks.json").write_text(json.dumps(result, indent=2)+"\n")
    sources = [Path(__file__), Path(n.__file__), ROOT/"sections/nonlinear_action.tex"]
    manifest = {
        "base_commit_at_execution": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "python": platform.python_version(),
        "dependencies": {name: importlib.metadata.version(name)
                         for name in ["numpy", "sympy"]},
        "seconds": round(time.monotonic()-start, 3),
        "GPU_jobs": 0, "remote_jobs": 0, "tracked_experiment_launches": 0,
        "sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in sources},
    }
    (OUT/"map-reproducibility.json").write_text(json.dumps(manifest, indent=2)+"\n")
    print("PASS: full pullback, changing tangent, and exact rotation control.")
    print("PASS: eight discarded coefficient channels, independent harmonic projections.")
    print("PASS: full linear Einstein map on pure-gauge and conformal controls.")
    print("No admissible non-gauge seed or nonlinear obstruction has been evaluated.")


if __name__ == "__main__":
    verify()

#!/usr/bin/env python3
"""Constraint bracket of a declared Hamiltonian-first Gaunt candidate.

This is not the Legendre transform of the original metric-coefficient
Galerkin Lagrangian. Canonical metric and momentum tensors are restricted
first; continuum ADM constraints are then evaluated pointwise with complete
Gaunt contractions. Lapse and shift smearings are independent.

At the flat spatial metric, constant Cartesian K_ij=delta_iz delta_jz gives
pi=diag(-1,-1,0), satisfying both vacuum constraints pointwise. The exact
angular projector test is checked independently by integrating the full
matrix Hamiltonian flow of the spatial generator, including the curvature
potential, in an orthonormal spherical frame.

No scattering amplitude, full Dirac stabilization, or quantization is
implemented here. Only Python's standard library is required.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from pathlib import Path
import json
import os
import platform
import sys
import time

from run import (
    K_COMPONENTS, ONE, WEIGHT, X, ZERO, add, as_decimal,
    check_constraint_data, deriv, eval_decimal, gauss_nodes, inner,
    legendre, mul, power, project_spin, report, scale, sub, weighted_inner,
)

ROOT = Path(__file__).resolve().parent
SPINS = (0, 0, 1, 2)
RADIAL_F = mul(power(sub(X, ONE), 2), power(sub(scale(ONE, 2), X), 2))
RADIAL_DF = deriv(RADIAL_F)
RADIAL_DDF = deriv(RADIAL_DF)


def integrate_between(polynomial, lower, upper):
    return sum((coefficient * F(upper ** (n + 1) - lower ** (n + 1), n + 1)
                for n, coefficient in enumerate(polynomial)), F(0))


def frobenius(a, b):
    """Positive tensor pairing; this is not the DeWitt pairing."""
    ar, at, av, ad = a
    br, bt, bv, bd = b
    return (inner(ar, br) + F(1, 2) * inner(at, bt)
            + 2 * weighted_inner(av, bv, 1)
            + 2 * weighted_inner(ad, bd, 2))


def pointwise_frobenius(a, b):
    ar, at, av, ad = a
    br, bt, bv, bd = b
    return add(add(mul(ar, br), scale(mul(at, bt), F(1, 2))),
               add(scale(mul(WEIGHT, mul(av, bv)), 2),
                   scale(mul(power(WEIGHT, 2), mul(ad, bd)), 2)))


def project_tensor(tensor, cutoff):
    return tuple(project_spin(component, cutoff, spin)
                 for component, spin in zip(tensor, SPINS))


def complement_tensor(tensor, cutoff):
    return tuple(sub(a, b) for a, b in zip(tensor, project_tensor(tensor, cutoff)))


def angular_tensors(ell):
    u = legendre(ell)
    du = deriv(u)
    x2 = power(X, 2)
    target = tuple(mul(u, a) for a in K_COMPONENTS)
    radial_gradient = (
        scale(mul(WEIGHT, u), -1),
        mul(add(ONE, x2), u),
        ZERO,
        scale(u, -F(1, 2)),
    )
    tangential_gradient = (
        scale(mul(WEIGHT, add(u, mul(X, du))), 2),
        ZERO,
        add(mul(X, u), mul(x2, du)),
        ZERO,
    )
    return target, radial_gradient, tangential_gradient


def gaunt_bracket(cutoff, ell):
    target, dr, dt = angular_tensors(ell)
    # Before restriction, pi K=0 makes the kinetic bracket vanish pointwise.
    assert pointwise_frobenius(dr, target) == ZERO
    assert pointwise_frobenius(dt, target) == ZERO
    qt = complement_tensor(target, cutoff)
    cr = frobenius(complement_tensor(dr, cutoff), qt)
    ct = frobenius(complement_tensor(dt, cutoff), qt)
    # The integration by parts int r^2 f f'=-int r f^2 is checked below.
    return 2 * (cr - ct), cr, ct


def edge_formula(ell):
    return F(2 * (ell + 1) * (ell + 2) * (2 * (ell + 1) ** 2 - 1),
             (2 * ell + 1) ** 2 * (2 * ell + 3) ** 2)


def exact_case(cutoff, ell, radial_factor):
    value, cr, ct = gaunt_bracket(cutoff, ell)
    expected = edge_formula(ell) if ell >= cutoff - 1 else F(0)
    assert value == expected, (cutoff, ell, value, expected)
    if ell >= cutoff - 1:
        c = F((ell + 1) * (ell + 2), (2 * ell + 1) * (2 * ell + 3))
        assert cr == -c * c / ((ell + 1) * (ell + 2))
        assert ct == -2 * c * c * (ell + 1) / (ell + 2)
        assert value > 0
    # Increasing the CANONICAL space removes this particular witness;
    # it is not an auxiliary numerical-quadrature change.
    enlarged = max(cutoff, ell + 2)
    enlarged_value, _, _ = gaunt_bracket(enlarged, ell)
    assert enlarged_value == 0
    u = legendre(ell)
    # Lie_beta gamma and the potential gradient are retained exactly.
    laplace_u = sub(mul(WEIGHT, deriv(deriv(u))), scale(mul(X, deriv(u)), 2))
    assert laplace_u == scale(u, -ell * (ell + 1))
    assert project_spin(u, cutoff, 0) == u
    assert project_spin(deriv(u), cutoff, 1) == deriv(u)
    assert project_spin(deriv(deriv(u)), cutoff, 2) == deriv(deriv(u))
    return {
        "L": cutoff, "N": cutoff + 1, "ell": ell,
        "radial_tail_pairing": str(cr), "tangential_tail_pairing": str(ct),
        "bracket_divided_by_k2_radial_factor": str(value),
        "full_radially_integrated_bracket": str(value * radial_factor),
        "unit_L2_smearings_bracket_divided_by_k2_radial_factor": str((2 * ell + 1) * value),
        "continuum_bracket": "0", "enlarged_canonical_cutoff": enlarged,
        "enlarged_canonical_bracket": str(enlarged_value),
    }


def matrix_add(a, b):
    return [[a[i][j] + b[i][j] for j in range(3)] for i in range(3)]


def matrix_scale(a, scalar):
    return [[scalar * a[i][j] for j in range(3)] for i in range(3)]


def matrix_transpose(a):
    return [list(row) for row in zip(*a)]


def matrix_product(a, b):
    return [[sum((a[i][k] * b[k][j] for k in range(3)), D(0))
             for j in range(3)] for i in range(3)]


def matrix_pair(a, b):
    return sum((a[i][j] * b[i][j] for i in range(3) for j in range(3)), D(0))


def tensor_matrix(polynomials, x, sine):
    a, b, c, d = (eval_decimal(p, x) for p in polynomials)
    weight = sine * sine
    return [[a, sine * c, D(0)],
            [sine * c, b / 2 + weight * d, D(0)],
            [D(0), D(0), b / 2 - weight * d]]


def numerical_integral(cutoff, ell, angular_order, radial_order, precision):
    """Direct canonical flow, without the exact tail-contraction formula."""
    u_poly = legendre(ell)
    du_poly, ddu_poly = deriv(u_poly), deriv(deriv(u_poly))
    target = project_tensor(angular_tensors(ell)[0], cutoff)
    totals = {name: D(0) for name in (
        "hamiltonian_flow_bracket", "gradient_bracket", "potential_piece",
        "unprojected_continuum_bracket",
    )}
    max_local_potential = D(0)
    # These quadratures are distinct from exact polynomial integration.
    radial_data = []
    for y, wr in gauss_nodes(radial_order, precision):
        radius = (3 + y) / 2
        radial_data.append((radius, wr, eval_decimal(RADIAL_F, radius),
                            eval_decimal(RADIAL_DF, radius),
                            eval_decimal(RADIAL_DDF, radius)))
    for x, wx in gauss_nodes(angular_order, precision):
        sine = (1 - x * x).sqrt()
        u, du, ddu = (eval_decimal(p, x) for p in (u_poly, du_poly, ddu_poly))
        retained_target = tensor_matrix(target, x, sine)
        k_matrix = tensor_matrix(K_COMPONENTS, x, sine)
        pi_matrix = [[k_matrix[i][j] - int(i == j) for j in range(3)] for i in range(3)]
        full_target = matrix_scale(k_matrix, u)
        for radius, wr, f, df, ddf in radial_data:
            # M^i_j = nabla_j beta^i, beta = f(r) u(theta) e_r.
            deformation = [[df * u, -f / radius * sine * du, D(0)],
                           [D(0), f / radius * u, D(0)],
                           [D(0), D(0), f / radius * u]]
            transpose = matrix_transpose(deformation)
            div_beta = (df + 2 * f / radius) * u
            metric_generator = matrix_add(deformation, transpose)
            momentum_gradient = matrix_add(
                matrix_add(matrix_product(deformation, pi_matrix),
                           matrix_product(pi_matrix, transpose)),
                matrix_scale(pi_matrix, -div_beta),
            )
            # C_gamma = -Hess(alpha) + I Delta(alpha), including curvature.
            hessian = [
                [ddf * u, -(df / radius - f / radius ** 2) * sine * du, D(0)],
                [-(df / radius - f / radius ** 2) * sine * du,
                 df / radius * u + f / radius ** 2 * ((1 - x * x) * ddu - x * du), D(0)],
                [D(0), D(0), df / radius * u - f / radius ** 2 * x * du],
            ]
            laplace = sum((hessian[i][i] for i in range(3)), D(0))
            potential_gradient = [[-hessian[i][j] + int(i == j) * laplace
                                   for j in range(3)] for i in range(3)]
            potential = -matrix_pair(metric_generator, potential_gradient)
            max_local_potential = max(max_local_potential, abs(potential))

            def direct_flow(target_matrix):
                qdot = matrix_scale(target_matrix, 2 * f)
                radial_qdot = matrix_scale(target_matrix, 2 * df)
                lie_qdot = matrix_add(
                    matrix_scale(radial_qdot, f * u),
                    matrix_add(matrix_product(transpose, qdot),
                               matrix_product(qdot, deformation)),
                )
                return matrix_pair(pi_matrix, lie_qdot)

            kinetic_flow = direct_flow(retained_target)
            kinetic_gradient = matrix_pair(momentum_gradient,
                                           matrix_scale(retained_target, 2 * f))
            measure = radius ** 2 * wr * wx
            totals["hamiltonian_flow_bracket"] += measure * (kinetic_flow + potential)
            totals["gradient_bracket"] += measure * (kinetic_gradient + potential)
            totals["potential_piece"] += measure * potential
            totals["unprojected_continuum_bracket"] += measure * (direct_flow(full_target) + potential)
    return {**totals, "max_abs_local_potential_piece": max_local_potential}


def numerical_case(cutoff, ell, radial_factor, config):
    with localcontext() as ctx:
        ctx.prec = config["decimal_precision"]
        target = as_decimal(radial_factor * gaunt_bracket(cutoff, ell)[0])
        checks = []
        tolerance = D(config["absolute_tolerance"])
        for angular_order, radial_order in config["quadratures"]:
            values = numerical_integral(cutoff, ell, angular_order, radial_order, ctx.prec)
            error = max(abs(values["hamiltonian_flow_bracket"] - target),
                        abs(values["gradient_bracket"] - target),
                        abs(values["potential_piece"]),
                        abs(values["unprojected_continuum_bracket"]))
            assert error < tolerance, (cutoff, ell, error)
            assert values["max_abs_local_potential_piece"] > 0
            checks.append({
                "angular_nodes": angular_order, "radial_nodes": radial_order,
                **{name: str(value) for name, value in values.items()},
                "max_verification_error": str(error),
            })
        difference = abs(D(checks[-1]["hamiltonian_flow_bracket"])
                         - D(checks[-2]["hamiltonian_flow_bracket"]))
        assert difference < tolerance
        return {
            "L": cutoff, "N": cutoff + 1, "ell": ell, "target": str(target),
            "checks": checks, "quadrature_difference": str(difference),
        }


def main(config):
    start = time.monotonic()
    print(json.dumps({"configuration": config, "python": sys.version,
                      "platform": platform.platform(), "commit": os.getenv("ORX_COMMIT"),
                      "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    report(phase="constraint_data", step=0, of=1, unit="cases")
    data = check_constraint_data()
    data.pop("metric_path")
    data.pop("lapse")
    data["candidate"] = "Canonical phase-space restriction with independent lapse/shift smearings; distinct from the original Lagrangian regulator."
    data["normalization"] = "16 pi G=1; normalized angular average dOmega/(4 pi); k=1; r in [1,2]."
    data["smearings"] = "alpha=f(r) P_ell, beta=f(r) P_ell e_r, f=(r-1)^2(2-r)^2."
    radial_factor = integrate_between(mul(X, power(RADIAL_F, 2)), 1, 2)
    radial_derivative = integrate_between(mul(power(X, 2), mul(RADIAL_F, RADIAL_DF)), 1, 2)
    assert radial_factor == F(1, 420)
    assert radial_derivative == -radial_factor
    data["radial_factor_int_r_f_squared"] = str(radial_factor)
    data["radial_factor_int_r_squared_f_df"] = str(radial_derivative)
    print("CONSTRAINT DATA " + json.dumps(data), flush=True)
    report(step=1, m={"constraint_residual": 0})

    total = sum(L + 1 for L in config["cutoffs_L"])
    exact = []
    report(phase="exact_brackets", step=0, of=total, unit="cases")
    for cutoff in config["cutoffs_L"]:
        for ell in range(cutoff + 1):
            result = exact_case(cutoff, ell, radial_factor)
            exact.append(result)
            value = float(F(result["bracket_divided_by_k2_radial_factor"]))
            if ell >= cutoff - 1:
                report(row="edges", m={"N": cutoff + 1, "ell": ell,
                                       "shell": "ell=L" if ell == cutoff else "ell=L-1",
                                       "mixed_bracket": value})
            if (cutoff, ell) in ((2, 1), (2, 2), (3, 1), (3, 3)):
                report(row="witnesses", m={"case": f"N={cutoff+1}, ell={ell}",
                    "mixed_bracket": value, "continuum_bracket": 0, "enlarged_bracket": 0})
            report(step=len(exact), m={"mixed_bracket": value,
                                      "exact_identity_residual": 0})
            print("EXACT BRACKET " + json.dumps(result), flush=True)

    report(phase="independent_flow", step=0, of=len(config["numerical_cases"]), unit="cases")
    numerical = []
    for done, (cutoff, ell) in enumerate(config["numerical_cases"], 1):
        result = numerical_case(cutoff, ell, radial_factor, config)
        numerical.append(result)
        error = max(D(row["max_verification_error"]) for row in result["checks"])
        report(row="flow_checks", m={"case": f"N={cutoff+1}, ell={ell}",
            "integrated_bracket": float(D(result["target"])),
            "verification_error_text": f"{error:.3E}"})
        report(step=done, m={"max_verification_error": float(error)})
        print("INDEPENDENT FLOW " + json.dumps(result), flush=True)

    report(phase="export", step=0, of=1, unit="reports")
    maximum_error = max(D(check["max_verification_error"])
                        for row in numerical for check in row["checks"])
    nonzero = sum(F(row["bracket_divided_by_k2_radial_factor"]) != 0 for row in exact)
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"),
        "run": os.getenv("ORX_RUN_ID"), "data": data, "exact_brackets": exact,
        "independent_checks": numerical, "nonzero_cases": nonzero,
        "max_verification_error": str(maximum_error),
        "interpretation": "This canonical Gaunt restriction does not inherit the continuum first-class mixed spatial/lapse constraint algebra.",
        "limitations": [
            "The candidate changes the original Lagrangian regulator; no equivalence is assumed.",
            "A nonzero bracket on the primary-secondary constraint surface is not a completed Dirac stabilization or a universal no-go.",
            "Cutoff-shell smearings vary with L. This is not a fixed-physical-momentum large-N scattering limit.",
            "Enlarging the canonical phase space removes this witness only; numerical auxiliary padding is not the same change.",
            "No finite-N one-loop four-graviton amplitude, dynamical corner quantization or physical Euclidean reconstruction is calculated."
        ],
        "status": "all specified exact and independent full-flow checks passed",
        "wall_seconds": time.monotonic() - start,
    }
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / "canonical_mixed_bracket_results.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    report(step=1, m={"max_verification_error": float(maximum_error),
                      "exact_identity_residual": 0, "checks_passed": len(exact) + len(numerical),
                      "nonzero_cases": nonzero, "constraint_residual": 0})
    print("FINAL " + json.dumps({"result": str(path), "exact_cases": len(exact),
        "independent_cases": len(numerical), "nonzero_cases": nonzero,
        "max_verification_error": str(maximum_error),
        "N3_ell1_normalized_bracket": str(edge_formula(1)),
        "N3_ell2_normalized_bracket": str(edge_formula(2)),
        "radial_factor": str(radial_factor), "status": result["status"],
        "wall_seconds": result["wall_seconds"]}), flush=True)

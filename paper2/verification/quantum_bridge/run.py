#!/usr/bin/env python3
"""Tracked local calculation for the common-action quantum research program.

This baseline tests the actual Galerkin metric ansatz at a point satisfying
the continuum vacuum constraints. It is not a one-loop graviton amplitude,
a full Dirac algorithm, or a proof of a physical boundary Hilbert space.

Only Python's standard library is required. Exact rational angular projections
are cross-checked against a separately evaluated finite-amplitude Hamiltonian
using Decimal Gauss--Legendre quadrature and matrix elimination.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from functools import lru_cache
from math import cos, factorial, pi
from pathlib import Path
import json
import os
import platform
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common_regulator_identities import (
    ONE, X, ZERO, add, average, deriv, inner, legendre, mul, power, scale, sub,
)

ROOT = Path(__file__).resolve().parent
WEIGHT = sub(ONE, power(X, 2))
# Orthonormal-frame K_rr, angular trace, K_rtheta/sin(theta), and
# ((K_thetatheta-K_phiphi)/2)/sin(theta)^2, for K_zz=1 in Cartesian axes.
K_COMPONENTS = (power(X, 2), WEIGHT, scale(X, -1), scale(ONE, F(1, 2)))


def report(**data):
    print("::orx:: " + json.dumps(data, allow_nan=False), flush=True)


def as_decimal(value):
    value = F(value)
    return D(value.numerator) / D(value.denominator)


@lru_cache(None)
def basis_poly(ell, spin):
    result = legendre(ell)
    for _ in range(spin):
        result = deriv(result)
    return result


@lru_cache(None)
def weighted_inner(a, b, spin):
    return average(mul(power(WEIGHT, spin), mul(a, b)))


def basis_norm(ell, spin):
    return F(factorial(ell + spin), (2 * ell + 1) * factorial(ell - spin))


@lru_cache(None)
def project_spin(a, L, spin):
    """Project (1-x^2)^(spin/2) a(x); return its polynomial part."""
    result = ZERO
    for ell in range(spin, L + 1):
        basis = basis_poly(ell, spin)
        coeff = weighted_inner(basis, a, spin) / basis_norm(ell, spin)
        result = add(result, scale(basis, coeff))
    return result


def tensor_inner(a, b):
    """<A_ij B^ij - tr(A) tr(B)>; all six tensor entries are accounted for.

    The absent rphi and theta-phi components vanish for these axisymmetric
    data. Their velocity blocks are retained in the theory and decouple at
    this background; this calculation does not truncate them out of the action.
    """
    ar, at, av, ad = a
    br, bt, bv, bd = b
    return (-inner(ar, bt) - inner(at, br) - F(1, 2) * inner(at, bt)
            + 2 * weighted_inner(av, bv, 1)
            + 2 * weighted_inner(ad, bd, 2))


def leakage(L, ell):
    u = legendre(ell)
    return tuple(sub(mul(u, component), project_spin(mul(u, component), L, spin))
                 for component, spin in zip(K_COMPONENTS, (0, 0, 1, 2)))


def check_constraint_data():
    a, b, c, d = K_COMPONENTS
    assert add(a, b) == ONE
    density = add(add(scale(mul(a, b), -2), scale(power(b, 2), -F(1, 2))),
                  add(scale(mul(WEIGHT, power(c, 2)), 2),
                      scale(mul(power(WEIGHT, 2), power(d, 2)), 2)))
    assert density == ZERO
    # For pi^ij/sqrt(gamma)=K^ij-gamma^ij K, the Cartesian components
    # are diag(-1,-1,0), so their derivatives vanish. Independently sum the
    # spherical-coordinate divergence terms. Factors r and r^2*sin(theta)
    # remove coordinate denominators in the r and theta components.
    radial = add(add(sub(ONE, scale(power(X, 2), 2)), add(power(X, 2), ONE)),
                 add(scale(WEIGHT, -2), scale(power(X, 2), -1)))
    polar = add(scale(mul(X, WEIGHT), -1), sub(X, power(X, 3)))
    assert radial == polar == ZERO
    assert tensor_inner(K_COMPONENTS, K_COMPONENTS) == 0
    for spin in (0, 1, 2):
        for ell in range(spin, 15):
            for other in range(spin, 15):
                expected = basis_norm(ell, spin) if ell == other else 0
                assert weighted_inner(basis_poly(ell, spin),
                                      basis_poly(other, spin), spin) == expected
    return {
        "spatial_geometry": "flat Euclidean metric on an annulus; Lambda=0",
        "extrinsic_curvature_cartesian": "K_ij=k delta_iz delta_jz; constant k=1 in normalized results",
        "hamiltonian_constraint_polynomial": "0",
        "momentum_constraint_cartesian": [0, 0, 0],
        "momentum_constraint_spherical_scaled": [str(radial[0]), str(polar[0]), "0"],
        "metric_path": "h_tt=-1-2 epsilon f(r) P_ell(cos(theta)); h_tr=0, h_rr=1, X=r^2, C=0, A=0",
        "lapse": "alpha=sqrt(1+2 epsilon f(r) P_ell); it is derived, not separately bandlimited",
        "normalization": "divide Hamiltonian density by r^2/(4G); k=1, f(r)=1 locally; use compactly supported f in the integrated test",
    }


def exact_case(L):
    leaks = [leakage(L, ell) for ell in range(L + 1)]
    for components in leaks:
        assert add(components[0], components[1]) == ZERO
    hessian = [[-2 * tensor_inner(a, b) for b in leaks] for a in leaks]
    for i, row in enumerate(hessian):
        for j, value in enumerate(row):
            if i != j or i < L - 1:
                assert value == 0
            else:
                assert value < 0
    if L == 2:
        assert hessian[1][1] == -F(4, 15)
        assert hessian[2][2] == -F(24, 175)
    if L >= 3:
        assert hessian[1][1] == 0
    return {
        "L": L, "N": L + 1,
        "basis": "unnormalized Legendre P_ell for ell=0,...,L",
        "hessian": [[str(value) for value in row] for row in hessian],
        "axisymmetric_lapse_block_rank": 2,
        "fixed_ell_1_curvature": str(hessian[1][1]),
        "cutoff_ell_L_curvature": str(hessian[L][L]),
        "unit_L2_cutoff_curvature": str((2 * L + 1) * hessian[L][L]),
        "unit_L2_penultimate_curvature": str((2 * L - 1) * hessian[L - 1][L - 1]),
    }


def eval_decimal(poly, x):
    result = D(0)
    for coeff in reversed(poly):
        result = result * x + as_decimal(coeff)
    return result


def legendre_decimal(n, x):
    previous, value = D(1), x
    if n == 0:
        return previous, D(0)
    for j in range(2, n + 1):
        previous, value = value, ((2 * j - 1) * x * value - (j - 1) * previous) / j
    derivative = n * (x * value - previous) / (x * x - 1)
    return value, derivative


@lru_cache(None)
def gauss_nodes(n, precision):
    """Independent quadrature; initial double root guesses are Newton-refined."""
    with localcontext() as ctx:
        ctx.prec = precision
        tolerance = D(10) ** (-(precision - 8))
        pairs = []
        for i in range(1, n // 2 + 1):
            root = D(str(cos(pi * (i - 0.25) / (n + 0.5))))
            for _ in range(30):
                value, derivative = legendre_decimal(n, root)
                step = value / derivative
                root -= step
                if abs(step) < tolerance:
                    break
            else:
                raise ArithmeticError("Gauss root did not converge")
            _, derivative = legendre_decimal(n, root)
            # Half the usual weight implements the normalized sphere average.
            weight = 1 / ((1 - root * root) * derivative * derivative)
            pairs.extend([(-root, weight), (root, weight)])
        if n % 2:
            _, derivative = legendre_decimal(n, D(0))
            pairs.append((D(0), 1 / (derivative * derivative)))
        pairs.sort()
        assert abs(sum(weight for _, weight in pairs) - 1) < D("1e-45")
        return tuple(pairs)


def solve_many(matrix, columns):
    n = len(matrix)
    rows = [list(row) + [column[i] for column in columns]
            for i, row in enumerate(matrix)]
    width = n + len(columns)
    for pivot in range(n):
        selected = max(range(pivot, n), key=lambda i: abs(rows[i][pivot]))
        rows[pivot], rows[selected] = rows[selected], rows[pivot]
        divisor = rows[pivot][pivot]
        if not divisor:
            raise ArithmeticError("Singular kinetic matrix")
        for j in range(pivot, width):
            rows[pivot][j] /= divisor
        for i in range(n):
            if i == pivot:
                continue
            factor = rows[i][pivot]
            for j in range(pivot, width):
                rows[i][j] -= factor * rows[pivot][j]
    return [[rows[i][n + j] for i in range(n)] for j in range(len(columns))]


def dot(a, b):
    return sum((x * y for x, y in zip(a, b)), D(0))


@lru_cache(None)
def numerical_inputs(L, ell, order, precision):
    nodes = gauss_nodes(order, precision)
    u = [eval_decimal(legendre(ell), x) for x, _ in nodes]
    tables, moments = {}, {}
    for spin in (0, 1, 2):
        basis = [basis_poly(j, spin) for j in range(spin, L + 1)]
        tables[spin] = [[eval_decimal(b, x) for b in basis] for x, _ in nodes]
        indices = (0, 1) if spin == 0 else (spin + 1,)
        moments[spin] = [[as_decimal(weighted_inner(b, K_COMPONENTS[index], spin))
                          for b in basis] for index in indices]
    return nodes, u, tables, moments


def full_hamiltonian(L, ell, epsilon, order, precision):
    nodes, u, tables, moments = numerical_inputs(L, ell, order, precision)
    terms = []
    for spin in (0, 1, 2):
        n = L - spin + 1
        matrix = [[D(0) for _ in range(n)] for _ in range(n)]
        for k, (x, quadrature_weight) in enumerate(nodes):
            lapse_squared = 1 + 2 * epsilon * u[k]
            if lapse_squared <= 0:
                raise ValueError("Metric has left the Lorentzian domain")
            weight = quadrature_weight * (1 - x * x) ** spin / lapse_squared.sqrt()
            values = tables[spin][k]
            for i in range(n):
                for j in range(i + 1):
                    matrix[i][j] += weight * values[i] * values[j]
        for i in range(n):
            for j in range(i):
                matrix[j][i] = matrix[i][j]
        solved = solve_many(matrix, moments[spin])
        if spin == 0:
            a, b = moments[spin]
            inverse_a, inverse_b = solved
            assert abs(dot(a, inverse_b) - dot(b, inverse_a)) < D("1e-45")
            terms.append(-2 * dot(a, inverse_b) - D("0.5") * dot(b, inverse_b))
        else:
            terms.append(2 * dot(moments[spin][0], solved[0]))
    return sum(terms, D(0))


def numerical_case(L, ell, config):
    precision = config["decimal_precision"]
    with localcontext() as ctx:
        ctx.prec = precision
        target_fraction = -2 * tensor_inner(leakage(L, ell), leakage(L, ell))
        target = as_decimal(target_fraction)
        steps = [D(value) for value in config["difference_steps"]]
        order_results = []
        for order in config["quadrature_orders"]:
            zero = full_hamiltonian(L, ell, D(0), order, precision)
            assert abs(zero) < D("1e-40")
            values = []
            for step in steps:
                positive = full_hamiltonian(L, ell, step, order, precision)
                negative = full_hamiltonian(L, ell, -step, order, precision)
                curvature = (positive - 2 * zero + negative) / (step * step)
                values.append(curvature)
            richardson = (4 * values[-1] - values[-2]) / 3
            error = abs(richardson - target)
            assert error < D(config["richardson_absolute_tolerance"]), (L, ell, error)
            order_results.append({
                "quadrature_order": order,
                "central_curvatures": [str(value) for value in values],
                "richardson_curvature": str(richardson),
                "absolute_error": str(error),
                "hamiltonian_at_zero": str(zero),
            })
        convergence = abs(D(order_results[-1]["richardson_curvature"])
                          - D(order_results[-2]["richardson_curvature"]))
        assert convergence < D(config["quadrature_absolute_tolerance"])
        return {"L": L, "N": L + 1, "ell": ell,
                "exact_curvature": str(target_fraction), "steps": [str(s) for s in steps],
                "precision": precision, "checks": order_results,
                "quadrature_difference": str(convergence)}


def main():
    start = time.monotonic()
    config = json.loads((ROOT / "config.json").read_text())
    if config["calculation"] == "complete_quartic_loops":
        from quartic_loop_check import main as quartic_loop_main
        quartic_loop_main(config)
        return
    if config["calculation"] == "reduced_weak_field_vertices":
        from reduced_vertices_check import main as vertices_main
        vertices_main(config)
        return
    if config["calculation"] == "reduced_gaunt_completion":
        from reduced_gaunt_completion import main as reduced_main
        reduced_main(config)
        return
    if config["calculation"] == "continuum_one_loop_all_plus":
        from one_loop_all_plus import main as all_plus_main
        all_plus_main(config)
        return
    if config["calculation"] == "kinetic_completion":
        from kinetic_completion import main as completion_main
        completion_main(config)
        return
    if config["calculation"] == "canonical_mixed_bracket":
        from canonical_mixed_bracket import main as bracket_main
        bracket_main(config)
        return
    if config["calculation"] != "constraint_surface_lapse":
        raise ValueError("Unsupported calculation")
    print(json.dumps({"configuration": config, "python": sys.version,
                      "platform": platform.platform(), "commit": os.getenv("ORX_COMMIT"),
                      "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    report(phase="constraint_data", step=0, of=1, unit="cases")
    data = check_constraint_data()
    report(step=1, m={"constraint_residual": 0})
    print("CONSTRAINT DATA " + json.dumps(data), flush=True)
    exact = []
    report(phase="exact_hessians", step=0, of=len(config["cutoffs_L"]), unit="cutoffs")
    for done, L in enumerate(config["cutoffs_L"], 1):
        result = exact_case(L)
        exact.append(result)
        report(row="cutoffs", m={"N": L + 1, "L": L,
                "boundary_defect": abs(float(F(result["unit_L2_cutoff_curvature"]))),
                "interior_curvature": float(F(result["fixed_ell_1_curvature"])),
                "lapse_block_rank": result["axisymmetric_lapse_block_rank"]})
        report(step=done, m={"boundary_defect": abs(float(F(result["unit_L2_cutoff_curvature"]))),
                             "lapse_block_rank": result["axisymmetric_lapse_block_rank"]})
        print("EXACT HESSIAN " + json.dumps(result), flush=True)
    numerical = []
    report(phase="independent_quadrature", step=0, of=len(config["numerical_cases"]), unit="cases")
    for done, (L, ell) in enumerate(config["numerical_cases"], 1):
        result = numerical_case(L, ell, config)
        numerical.append(result)
        error = float(D(result["checks"][-1]["absolute_error"]))
        report(row="quadrature", m={"case": f"N={L+1}, ell={ell}", "fd_error": error,
                                     "quadrature_difference": float(D(result["quadrature_difference"]))})
        report(step=done, m={"fd_error": error})
        print("NUMERICAL CHECK " + json.dumps(result), flush=True)
    report(phase="export", step=0, of=1, unit="reports")
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"),
        "run": os.getenv("ORX_RUN_ID"), "data": data,
        "exact_hessians": exact, "independent_checks": numerical,
        "interpretation": "A lapse primary/secondary bracket is nonzero on the primary-secondary constraint surface of the Galerkin pullback. The continuum first-class structure cannot be imported unchanged.",
        "limitations": [
            "Full Dirac stabilization, tertiary constraints and physical degree counting remain open.",
            "The curvature is normalized at a local flat-slice phase-space point; it is not a scattering observable.",
            "Cutoff-shell lapse variations change with N. Their nonuniform defect is not a no-go theorem for fixed-physical-scale convergence.",
            "No one-loop four-graviton amplitude, corner quantization or Euclidean reconstruction has been completed by this test."
        ],
        "status": "all specified exact and independent numerical checks passed",
        "wall_seconds": time.monotonic() - start,
    }
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / "lapse_constraint_results.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    maximum_error = max(float(D(row["checks"][-1]["absolute_error"])) for row in numerical)
    report(step=1, m={"max_fd_error": maximum_error, "constraint_residual": 0,
                      "checks_passed": len(exact) + len(numerical) + 1})
    print("FINAL " + json.dumps({"result": str(path), "exact_cutoffs": len(exact),
                                "independent_cases": len(numerical), "max_fd_error": maximum_error,
                                "N3_ell1_curvature": exact[0]["fixed_ell_1_curvature"],
                                "N3_ell2_curvature": exact[0]["cutoff_ell_L_curvature"],
                                "status": result["status"], "wall_seconds": result["wall_seconds"]}), flush=True)


if __name__ == "__main__":
    main()

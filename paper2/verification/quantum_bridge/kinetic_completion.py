#!/usr/bin/env python3
"""Schur completion and angular determinant terms for the ADM kinetic test.

The input matrices are complete Gaunt contractions of the same Einstein
kinetic restriction as the parent calculation. This is not the full
gravitational Hessian: spatial dynamics, non-axisymmetric fluctuation
determinants, constraints/ghosts and boundary terms remain to be derived.
In particular no reported number is a one-loop four-graviton amplitude.
"""

from decimal import Decimal as D, localcontext
from fractions import Fraction as F
from math import comb
from pathlib import Path
import json
import os
import platform
import sys
import time

from run import (
    K_COMPONENTS, ROOT, as_decimal, basis_norm, basis_poly, full_hamiltonian,
    legendre, mul, power, report, tensor_inner, leakage, weighted_inner,
)


def zeros(n, m=None):
    return [[F(0) for _ in range(n if m is None else m)] for _ in range(n)]


def eye(n):
    return [[F(i == j) for j in range(n)] for i in range(n)]


def add(A, B):
    return [[a + b for a, b in zip(ar, br)] for ar, br in zip(A, B)]


def scale(A, value):
    return [[value * a for a in row] for row in A]


def multiply(A, B):
    result = zeros(len(A), len(B[0]))
    for i, row in enumerate(A):
        for k, value in enumerate(row):
            if not value:
                continue
            for j, other in enumerate(B[k]):
                if other:
                    result[i][j] += value * other
    return result


def trace(A):
    return sum((row[i] for i, row in enumerate(A)), F(0))


def block(A, rows, columns):
    return [[A[i][j] for j in columns] for i in rows]


def series_product(A, B):
    result = []
    for n in range(len(A)):
        value = zeros(len(A[0]), len(B[0][0]))
        for k in range(n + 1):
            value = add(value, multiply(A[k], B[n - k]))
        result.append(value)
    return result


def series_inverse(A):
    n = len(A[0])
    assert A[0] == eye(n)
    inverse = [eye(n)]
    for order in range(1, len(A)):
        value = zeros(n)
        for k in range(1, order + 1):
            value = add(value, multiply(A[k], inverse[order - k]))
        inverse.append(scale(value, -1))
    product = series_product(A, inverse)
    assert product[0] == eye(n)
    assert all(value == zeros(n) for value in product[1:])
    return inverse


def logdet_series(A):
    """Derive coefficients by integrating tr(A^-1 dA/depsilon)."""
    inverse = series_inverse(A)
    result = [F(0)]
    for n in range(1, len(A)):
        value = sum(((k + 1) * trace(multiply(inverse[n - 1 - k], A[k + 1]))
                     for k in range(n)), F(0)) / n
        result.append(value)
    if len(A) > 4:
        X1, X2, X3, X4 = A[1:5]
        square = multiply(X1, X1)
        # Independent cyclic fourth-order expansion, including the term
        # that would originate from a sixth-order action vertex.
        cyclic = (trace(X4) - trace(multiply(X1, X3))
                  - F(1, 2) * trace(multiply(X2, X2))
                  + trace(multiply(square, X2))
                  - F(1, 4) * trace(multiply(square, square)))
        assert cyclic == result[4]
    return result


def angular_series(J, ell, spin, order):
    """D0^-1 M(eps), with M_ij=<B_i B_j / sqrt(1+2 eps u)>.

    D0 is the diagonal unweighted Gram matrix. The constant term is I,
    and all coefficients are exact rationals in the unnormalized basis.
    The logarithm has the same background-dependent coefficients as log M.
    """
    basis = [basis_poly(j, spin) for j in range(spin, J + 1)]
    u = legendre(ell)
    result = []
    for n in range(order + 1):
        coefficient = F((-1) ** n * comb(2 * n, n), 2 ** n)
        un = power(u, n)
        matrix = [[coefficient * weighted_inner(a, mul(un, b), spin)
                   / basis_norm(i + spin, spin)
                   for b in basis] for i, a in enumerate(basis)]
        result.append(matrix)
    assert result[0] == eye(len(basis))
    return result


def bilinear(left, A, right):
    return sum((left[i] * value * right[j]
                for i, row in enumerate(A) for j, value in enumerate(row)), F(0))


def kinetic_coefficients(inverses, cutoff):
    """Hamiltonian coefficients at fixed parent canonical momenta."""
    order = len(inverses[0]) - 1
    result = [F(0) for _ in range(order + 1)]
    for spin in (0, 1, 2):
        basis = [basis_poly(j, spin) for j in range(spin, cutoff + 1)]
        indices = (0, 1) if spin == 0 else (spin + 1,)
        moments = [[weighted_inner(b, K_COMPONENTS[index], spin) for b in basis]
                   for index in indices]
        coefficients = [[value / basis_norm(j + spin, spin) for j, value in enumerate(column)]
                        for column in moments]
        for n, inverse in enumerate(inverses[spin]):
            if spin == 0:
                a, b = moments
                ac, bc = coefficients
                assert bilinear(a, inverse, bc) == bilinear(b, inverse, ac)
                result[n] += (-2 * bilinear(a, inverse, bc)
                              - F(1, 2) * bilinear(b, inverse, bc))
            else:
                result[n] += 2 * bilinear(moments[0], inverse, coefficients[0])
    return result


def exact_block_case(L, ell, J, order):
    naive, effective, enlarged = {}, {}, {}
    determinant = []
    for spin in (0, 1, 2):
        A = angular_series(J, ell, spin, order)
        split = L - spin + 1
        retained = range(split)
        discarded = range(split, J - spin + 1)
        PP = [block(a, retained, retained) for a in A]
        PQ = [block(a, retained, discarded) for a in A]
        QP = [block(a, discarded, retained) for a in A]
        QQ = [block(a, discarded, discarded) for a in A]
        correction = series_product(series_product(PQ, series_inverse(QQ)), QP)
        E = [add(pp, scale(c, -1)) for pp, c in zip(PP, correction)]
        inverse_E = series_inverse(E)
        inverse_A = series_inverse(A)
        assert inverse_E == [block(a, retained, retained) for a in inverse_A]
        naive[spin] = series_inverse(PP)
        effective[spin] = inverse_E
        enlarged[spin] = inverse_A
        full_log, low_log = logdet_series(A), logdet_series(E)
        high_log, naive_log = logdet_series(QQ), logdet_series(PP)
        assert all(full_log[n] == low_log[n] + high_log[n] for n in range(order + 1))
        determinant.append({
            "spin_block": spin,
            "full_logdet_coefficients": [str(x) for x in full_log],
            "effective_low_logdet_coefficients": [str(x) for x in low_log],
            "discarded_logdet_coefficients": [str(x) for x in high_log],
            "naive_low_logdet_coefficients": [str(x) for x in naive_log],
            "quartic_error_if_discarded_determinant_omitted": str(high_log[4]),
            "quartic_error_if_mixed_schur_terms_omitted": str(low_log[4] - naive_log[4]),
            "exact_factorization": True,
            "cyclic_fourth_order_formula": True,
        })
    hn = kinetic_coefficients(naive, L)
    he = kinetic_coefficients(effective, L)
    ha = kinetic_coefficients(enlarged, J)
    assert he == ha
    assert hn[0] == hn[1] == he[0] == he[1] == 0
    predicted = -tensor_inner(leakage(L, ell), leakage(L, ell))
    assert hn[2] == predicted
    assert J >= ell + 2 and he[2] == 0
    correction = he[2] - hn[2]
    assert correction == -predicted > 0
    # These controls must be sensitive: dropping determinant blocks must
    # actually change at least one fourth-order coefficient in this case.
    assert any(F(d["quartic_error_if_discarded_determinant_omitted"]) != 0 for d in determinant)
    assert any(F(d["quartic_error_if_mixed_schur_terms_omitted"]) != 0 for d in determinant)
    return {
        "L": L, "N": L + 1, "ell": ell, "auxiliary_J": J,
        "naive_H_coefficients": [str(x) for x in hn],
        "schur_H_coefficients": [str(x) for x in he],
        "enlarged_H_coefficients": [str(x) for x in ha],
        "naive_lapse_curvature": str(2 * hn[2]),
        "completed_lapse_curvature": str(2 * he[2]),
        "required_H_quartic_correction": str(correction),
        "determinants": determinant,
        "scope": "axisymmetric angular kinetic matrices only; no frequency/radial integral, other m loops, ghosts or gravitational measure",
    }


def finite_amplitude(J, ell, amplitude, config):
    with localcontext() as context:
        context.prec = config["decimal_precision"]
        epsilon = D(amplitude)
        values = [full_hamiltonian(J, ell, epsilon, order, context.prec)
                  for order in config["quadrature_orders"]]
        difference = abs(values[-1] - values[-2])
        assert difference < D(config["quadrature_absolute_tolerance"]), (J, amplitude, difference)
        return {"auxiliary_J": J, "ell": ell, "epsilon": amplitude,
                "H_values": [str(v) for v in values],
                "absolute_H": str(abs(values[-1])),
                "quadrature_difference": str(difference),
                "continuum_phase_space_restricted_H": "0"}


def main(config):
    started = time.monotonic()
    print(json.dumps({"configuration": config, "python": sys.version,
                      "platform": platform.platform(), "commit": os.getenv("ORX_COMMIT"),
                      "run": os.getenv("ORX_RUN_ID"), "experiment": os.getenv("ORX_EXPERIMENT")}), flush=True)
    exact = []
    report(phase="exact_completion", step=0, of=len(config["block_cases"]), unit="cases")
    for done, (L, ell, J) in enumerate(config["block_cases"], 1):
        result = exact_block_case(L, ell, J, config["series_order"])
        exact.append(result)
        report(row="completion", m={"case": f"N={L+1}, ell={ell}, J={J}",
                "naive_curvature": float(F(result["naive_lapse_curvature"])),
                "completed_curvature": float(F(result["completed_lapse_curvature"])),
                "quartic_correction": float(F(result["required_H_quartic_correction"]))})
        for entry in result["determinants"]:
            report(row="determinants", m={"case": f"N={L+1}, ell={ell}, s={entry['spin_block']}",
                "omitted_high_error": float(F(entry["quartic_error_if_discarded_determinant_omitted"])),
                "omitted_mixed_error": float(F(entry["quartic_error_if_mixed_schur_terms_omitted"])),
                "determinant_identity_residual": 0})
        report(step=done, m={"completed_curvature": 0, "determinant_identity_residual": 0})
        print("EXACT COMPLETION " + json.dumps(result), flush=True)
    finite = []
    total = len(config["auxiliary_cutoffs"]) * len(config["amplitudes"])
    report(phase="finite_amplitude", step=0, of=total, unit="cases")
    for amplitude in config["amplitudes"]:
        for J in config["auxiliary_cutoffs"]:
            result = finite_amplitude(J, config["finite_amplitude_ell"], amplitude, config)
            finite.append(result)
            report(row="padding", m={"auxiliary_J": J, "epsilon": amplitude,
                                      "absolute_H": float(D(result["absolute_H"]))})
            report(step=len(finite), m={"quadrature_difference": float(D(result["quadrature_difference"]))})
            print("FINITE AMPLITUDE " + json.dumps(result), flush=True)
    report(phase="export", step=0, of=1, unit="reports")
    result = {
        "config": config, "source_commit": os.getenv("ORX_COMMIT"), "run": os.getenv("ORX_RUN_ID"),
        "parent_run": "82f9a0a7-3fd7-432b-9913-27738446ad38",
        "exact_completion": exact, "finite_amplitude_checks": finite,
        "limitations": [
            "This completes only the tested kinetic restriction; the full gravitational constraint algebra is not derived.",
            "The determinant checks are for m=0 angular kinetic blocks, not the full fluctuation operator or a graviton amplitude.",
            "Finite padding is exact only to the stated perturbative order; the finite-amplitude residual is recorded separately.",
            "A dynamical integration of discarded fields also needs spatial derivatives, propagators, boundary data and the actual gauge/constraint measure.",
            "Corner-state derivation and Euclidean--Lorentzian reconstruction for that physical gravitational theory remain open."
        ],
        "status": "all specified exact completion and independent finite-amplitude checks passed",
        "wall_seconds": time.monotonic() - started,
    }
    output = ROOT / "output"
    output.mkdir(exist_ok=True)
    path = output / "kinetic_completion_results.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    report(step=1, m={"completed_curvature": 0, "determinant_identity_residual": 0,
                      "checks_passed": len(exact) + len(finite)})
    print("FINAL " + json.dumps({"result": str(path), "exact_cases": len(exact),
             "finite_amplitude_cases": len(finite), "completed_lapse_curvature": "0 exactly in every stated case",
             "determinant_factorization_residual": "0 exactly through order four",
             "status": result["status"], "wall_seconds": result["wall_seconds"]}), flush=True)

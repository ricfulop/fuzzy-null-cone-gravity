#!/usr/bin/env python3
"""Exact checks for the proposed common-regulator convention.

Uses only Python's standard library and rational polynomial arithmetic.
The pairing is the normalized sphere integral for axisymmetric scalars,
<f,g> = (1/2) integral_{-1}^1 f(x) g(x) dx, x = cos(theta).
P_L retains Legendre degrees 0,...,L; matrix rank is N=L+1.

These check product/evaluation identities and one already-derived conformal
restriction of Paper I's bulk action. They are not a full gravitational vertex,
constraint, path-integral, or scattering calculation.
"""

from fractions import Fraction as F
from functools import lru_cache
from itertools import product
import json


def poly(values):
    values = list(map(F, values))
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values) or (F(0),)


ZERO = poly([0])
ONE = poly([1])
X = poly([0, 1])


def add(a, b):
    return poly((a[i] if i < len(a) else 0) +
                (b[i] if i < len(b) else 0)
                for i in range(max(len(a), len(b))))


def scale(a, c):
    return poly(c * v for v in a)


def sub(a, b):
    return add(a, scale(b, -1))


@lru_cache(None)
def mul(a, b):
    values = [F(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            values[i + j] += x * y
    return poly(values)


def power(a, n):
    result = ONE
    for _ in range(n):
        result = mul(result, a)
    return result


def deriv(a):
    return poly(i * a[i] for i in range(1, len(a)))


def average(a):
    return sum((c / (i + 1) for i, c in enumerate(a) if i % 2 == 0), F(0))


@lru_cache(None)
def inner(a, b):
    return average(mul(a, b))


@lru_cache(None)
def legendre(j):
    if j == 0:
        return ONE
    if j == 1:
        return X
    return scale(sub(scale(mul(X, legendre(j - 1)), 2 * j - 1),
                     scale(legendre(j - 2), j - 1)), F(1, j))


@lru_cache(None)
def project(a, L):
    result = ZERO
    for j in range(L + 1):
        result = add(result, scale(legendre(j), (2 * j + 1) * inner(a, legendre(j))))
    return result


def coefficient(a, j):
    return (2 * j + 1) * inner(a, legendre(j))


def coefficients(a):
    return {str(j): str(coefficient(a, j))
            for j in range(len(a)) if coefficient(a, j)}


def projected_power(a, n, L):
    result = ONE
    for _ in range(n):
        result = project(mul(result, a), L)
    return result


def recoupled_pair(a, b, c, d, intermediate_L):
    """Complete scalar four-harmonic tensor via two three-harmonic tensors."""
    ab, cd = mul(a, b), mul(c, d)
    return sum(((2 * j + 1) * inner(ab, legendre(j)) * inner(cd, legendre(j))
                for j in range(intermediate_L + 1)), F(0))


def main():
    for i, j in product(range(9), repeat=2):
        assert inner(legendre(i), legendre(j)) == (F(1, 2 * i + 1) if i == j else 0)

    # Exhaustive small-band tensor checks, with direct monomial integration
    # independent of the Legendre recoupling expression.
    tensor_checks = []
    for L in (2, 3):
        basis = [legendre(j) for j in range(L + 1)]
        for a, b, c in product(basis, repeat=3):
            assert inner(project(mul(a, b), L), c) == inner(mul(a, b), c)
        for a, b, c, d in product(basis, repeat=4):
            ab, cd = mul(a, b), mul(c, d)
            pab, pcd = project(ab, L), project(cd, L)
            qab, qcd = sub(ab, pab), sub(cd, pcd)
            full = inner(ab, cd)
            assert full == inner(pab, pcd) + inner(qab, qcd)
            assert full == recoupled_pair(a, b, c, d, 2 * L)
            assert full == recoupled_pair(a, c, b, d, 2 * L)
            assert full == recoupled_pair(a, d, b, c, 2 * L)
        tensor_checks.append({"L": L, "N": L + 1,
                              "cubic_tuples": (L + 1) ** 3,
                              "quartic_tuples": (L + 1) ** 4,
                              "exact_identities": "passed"})

    u, L = legendre(2), 2
    u2 = power(u, 2)
    low, high = project(u2, L), sub(u2, project(u2, L))
    full4, projected4, missing4 = inner(u2, u2), inner(low, low), inner(high, high)
    assert full4 == F(3, 35)
    assert projected4 == F(69, 1225)
    assert missing4 == F(36, 1225)
    assert missing4 / full4 == F(12, 35)

    # For the same external mode at L=4, quartic agreement does not imply
    # sextic agreement. The lost l=6 channel contributes to the six-field tensor.
    L6 = 4
    u3 = power(u, 3)
    q_u3 = sub(u3, project(u3, L6))
    full6 = average(power(u, 6))
    projected6 = average(projected_power(u, 6, L6))
    assert average(projected_power(u, 4, L6)) == full4
    assert full6 - projected6 == inner(q_u3, q_u3) == F(324, 77077)
    assert inner(project(u3, 6), project(u3, 6)) == full6

    # A within-band star inverse solves P_L[(1+eps*u)*v]=1.
    # Its order-k coefficient is (-P_L M_u P_L)^k 1; compare with
    # projecting the ordinary pointwise inverse's geometric series.
    for k in range(3):
        assert project(power(u, k), L) == projected_power(u, k, L)
    inverse_difference3 = sub(scale(project(power(u, 3), L), -1),
                              scale(projected_power(u, 3, L), -1))
    assert inverse_difference3 == scale(u, -F(36, 245))

    # Finite bands do not commute with nonlinear field coordinates:
    # exp(P_L log(1+eps*u)) - (1+eps*u) = eps^2 Q_L(u^2)/2 + O(eps^3).
    reparam_difference2 = scale(sub(u2, project(u2, L)), F(1, 2))
    assert reparam_difference2 == scale(legendre(4), F(9, 35))

    # The conformal family g=(1+eps*P_2)*g_Schwarzschild is representable
    # in Paper I at N>=3. Use its already-derived angular Einstein density,
    # (3/2)<|grad A|^2/A>, with common dimensional prefactors suppressed.
    # The comparator is explicitly defined here; Paper II has no unique
    # fully specified projected inverse/action to silently identify it with.
    grad_sq = mul(sub(ONE, power(X, 2)), power(deriv(u), 2))
    conformal = []
    for n in (2, 3, 4):
        exact = F(3, 2) * (-1) ** (n - 2) * inner(grad_sq, power(u, n - 2))
        projected = F(3, 2) * (-1) ** (n - 2) * inner(
            project(grad_sq, L), projected_power(u, n - 2, L))
        if n < 4:
            assert exact == projected
        conformal.append({"order": n, "full_pointwise": str(exact),
                          "specified_projected_comparator": str(projected),
                          "difference_full_minus_projected": str(exact - projected)})
    assert conformal[-1]["full_pointwise"] == "9/35"
    assert conformal[-1]["specified_projected_comparator"] == "531/1225"
    assert F(3, 2) * inner(sub(grad_sq, project(grad_sq, L)), high) == -F(216, 1225)

    result = {
        "arithmetic": "exact rational; no floating-point tolerances",
        "scope": "axisymmetric scalar angular identities and the specified conformal Einstein restriction; not full spin-2 scattering or quantum consistency",
        "tensor_checks": tensor_checks,
        "quartic_example": {"L": L, "N": L + 1, "field": "P_2(cos(theta))",
                            "product_legendre_coefficients": coefficients(u2),
                            "full": str(full4), "projected": str(projected4),
                            "missing": str(missing4), "relative_missing": str(missing4 / full4)},
        "sextic_example": {"L": L6, "N": L6 + 1, "field": "P_2(cos(theta))",
                           "quartic_agrees": True, "full": str(full6),
                           "repeated_projected_product": str(projected6),
                           "missing": str(full6 - projected6),
                           "discarded_cube_coefficients": coefficients(q_u3)},
        "inverse_series_difference_at_order_3": coefficients(inverse_difference3),
        "field_redefinition_difference_at_order_2": coefficients(reparam_difference2),
        "conformal_Einstein_density": conformal,
        "status": "all exact checks passed",
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

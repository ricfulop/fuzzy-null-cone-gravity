"""Exact terminal preservation of the original constraints over QQ(r,T)."""

from contextlib import nullcontext

from sympy import Float, nan, oo, zoo
from sympy.polys.matrices import DomainMatrix

from exact_linear_algebra import FIELD


def require_exact_finite(expr):
    if expr.has(Float, nan, oo, -oo, zoo):
        raise ValueError("Inexact or nonfinite terminal coefficient")


def field_matrix(matrix):
    entries = {}
    for (row, column), expr in matrix.todok().items():
        require_exact_finite(expr)
        value = FIELD.from_sympy(expr)
        if value:
            entries.setdefault(row, {})[column] = value
    return DomainMatrix.from_dod(entries, matrix.shape, FIELD)


def field_derivative(matrix):
    r, t = FIELD.gens
    t_prime = -10 / (100*r**2 + 1)
    entries = {}
    for i, row in matrix.to_dod().items():
        for j, value in row.items():
            derivative = value.diff(r) + t_prime * value.diff(t)
            if derivative:
                entries.setdefault(i, {})[j] = derivative
    return DomainMatrix.from_dod(entries, matrix.shape, FIELD)


def verify_terminal_preservation(constraints, dynamics, controls, embedding, *, stage=nullcontext):
    with stage("preservation_field_conversion"):
        c, k, b, m = map(field_matrix, (constraints, dynamics, controls, embedding))
    with stage("preservation_DC"):
        dc = field_derivative(c)
    with stage("preservation_CK"):
        ck = c.matmul(k)
    with stage("preservation_sum"):
        s = dc.add(ck)
    with stage("preservation_CB"):
        cb = c.matmul(b)
    with stage("preservation_SM"):
        sm = s.matmul(m)
    with stage("preservation_zero_checks"):
        if not cb.is_zero_matrix:
            raise ValueError("remaining controls violate constraint preservation")
        if not sm.is_zero_matrix:
            raise ValueError("constraint subspace is not invariant")

"""Exact right kernel of the accumulated stationary constraints over QQ(r,T)."""

from contextlib import nullcontext

from sympy import Float, Matrix
from sympy.polys.matrices import DomainMatrix

from exact_linear_algebra import FIELD


def expression_matrix(matrix):
    values = [FIELD.to_sympy(value) for value in matrix.to_list_flat()]
    return Matrix(*matrix.shape, values)


def kernel_embedding(constraints, *, stage=nullcontext):
    """Return M, free coordinates, nonzero RREF rows, and pivot columns."""
    with stage("kernel_field_conversion"):
        entries = {}
        for (row, column), expr in constraints.todok().items():
            if expr.has(Float):
                raise ValueError("Inexact coefficient in terminal constraint matrix")
            value = FIELD.from_sympy(expr)
            if value:
                entries.setdefault(row, {})[column] = value
        original = DomainMatrix.from_dod(entries, constraints.shape, FIELD)

    with stage("kernel_constraint_rref"):
        if 0 in original.shape:
            numerator, denominator, pivots = original, FIELD.one, ()
        else:
            numerator, denominator, pivots = original.rref_den(method="CD", keep_domain=False)

    with stage("kernel_denominator_normalization"):
        scale = FIELD.convert(denominator, numerator.domain)
        if not scale:
            raise ValueError("Terminal RREF denominator is zero")
        reduced = numerator.convert_to(FIELD) / scale
        rows, columns = original.shape
        rank = len(pivots)
        if reduced.extract(range(rows), pivots) != DomainMatrix.eye((rows, rank), FIELD):
            raise ValueError("Terminal RREF pivot columns do not have unit pivots")
        if not reduced.extract(range(rank, rows), range(columns)).is_zero_matrix:
            raise ValueError("Terminal RREF has a nonzero trailing row")

    with stage("kernel_build_embedding"):
        free = [column for column in range(columns) if column not in pivots]
        entries = {column: {j: FIELD.one} for j, column in enumerate(free)}
        for i, pivot in enumerate(pivots):
            for j, column in enumerate(free):
                value = reduced[i, column].element
                if value:
                    entries.setdefault(pivot, {})[j] = -value
        embedding = DomainMatrix.from_dod(entries, (columns, len(free)), FIELD)
        if embedding.extract(free, range(len(free))) != DomainMatrix.eye(len(free), FIELD):
            raise ValueError("Terminal embedding free-coordinate rows do not form the identity")

    with stage("kernel_original_constraint_product"):
        if not original.matmul(embedding).is_zero_matrix:
            raise ValueError("Original constraint matrix times terminal embedding is nonzero")

    with stage("kernel_to_expressions"):
        basis = reduced.extract(range(rank), range(columns))
        return expression_matrix(embedding), free, expression_matrix(basis), list(pivots)

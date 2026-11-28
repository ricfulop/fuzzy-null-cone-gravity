"""Optional exact row operations for the fixed-epsilon stationary reducer.

Inputs are SymPy matrices over QQ(r,T), with the reducer's positive symbols.
No derivative, radius specialization, or numerical rank threshold is used.
"""
from sympy import Float, Matrix, QQ, symbols, zeros
from sympy.polys.matrices import DomainMatrix


R, T = symbols("r T", positive=True)
FIELD = QQ.frac_field(R, T)


def _rref_den(matrix):
    """Use the associated polynomial ring; keep the denominator for the caller."""
    entries = {}
    for (row, col), value in matrix.todok().items():
        if value.has(Float):
            raise ValueError("Inexact entry in a QQ(r,T) matrix")
        element = FIELD.from_sympy(value)
        if element:
            entries.setdefault(row, {})[col] = element
    exact = DomainMatrix.from_dod(entries, matrix.shape, FIELD)
    return exact.rref_den(method="CD", keep_domain=False)


def rref_with_side(left, right):
    """Solve left pivots; rows below those pivots retain all state constraints.

    The right block is also row-reduced. Upper rows can differ from control-only
    elimination by multiples of the residual state equations.
    """
    left, right = Matrix(left), Matrix(right)
    if left.rows != right.rows:
        raise ValueError("Control and state blocks must have the same row count")
    if left.rows == 0 or left.cols + right.cols == 0:
        return left, right, []
    numerator, denominator, pivots = _rref_den(left.row_join(right))
    scale = FIELD.convert(denominator, numerator.domain)
    reduced = (numerator.convert_to(FIELD) / scale).to_Matrix()
    return (reduced[:, :left.cols], reduced[:, left.cols:],
            [pivot for pivot in pivots if pivot < left.cols])


def independent_extension(existing, candidates):
    """Select original candidate rows extending the existing span, in order."""
    existing, candidates = Matrix(existing), Matrix(candidates)
    if existing.cols != candidates.cols:
        raise ValueError("All constraint rows must have the same state width")
    if candidates.rows == 0 or candidates.cols == 0:
        return zeros(0, candidates.cols)
    _, _, pivots = _rref_den(existing.col_join(candidates).T)
    selected = [pivot - existing.rows for pivot in pivots if pivot >= existing.rows]
    return candidates.extract(selected, range(candidates.cols))

"""Optional exact row operations for the fixed-epsilon stationary reducer.

Inputs are SymPy matrices over QQ(r,T), with the reducer's positive symbols.
No derivative, radius specialization, or numerical rank threshold is used.
"""
from sympy import Float, Matrix, QQ, symbols, zeros
from sympy.polys.matrices import DomainMatrix


R, T = symbols("r T", positive=True)
FIELD = QQ.frac_field(R, T)


def _field_value(value):
    if value.has(Float):
        raise ValueError("Inexact entry in a QQ(r,T) matrix")
    return FIELD.from_sympy(value)


def _rref_den(matrix):
    """Use the associated polynomial ring; keep the denominator for the caller."""
    entries = {}
    for (row, col), value in matrix.todok().items():
        element = _field_value(value)
        if element:
            entries.setdefault(row, {})[col] = element
    exact = DomainMatrix.from_dod(entries, matrix.shape, FIELD)
    return exact.rref_den(method="CD", keep_domain=False)


def rref_with_side(left, right):
    """Apply the historical first-nonzero control pivots to both blocks over QQ(r,T)."""
    left, right = Matrix(left), Matrix(right)
    if left.rows != right.rows:
        raise ValueError("Control and state blocks must have the same row count")
    rows = [[_field_value(value) for value in control_row + state_row]
            for control_row, state_row in zip(left.tolist(), right.tolist())]
    pivots = []
    for column in range(left.cols):
        pivot_row = len(pivots)
        found = next((i for i in range(pivot_row, left.rows) if rows[i][column]), None)
        if found is None:
            continue
        rows[pivot_row], rows[found] = rows[found], rows[pivot_row]
        pivot = rows[pivot_row][column]
        normalized = [value / pivot for value in rows[pivot_row]]
        rows[pivot_row] = normalized
        for i, row in enumerate(rows):
            factor = row[column]
            if i != pivot_row and factor:
                rows[i] = [value - factor * basis for value, basis in zip(row, normalized)]
        pivots.append(column)
        if len(pivots) == left.rows:
            break

    controls = [FIELD.to_sympy(value) for row in rows for value in row[:left.cols]]
    states = [FIELD.to_sympy(value) for row in rows for value in row[left.cols:]]
    return Matrix(left.rows, left.cols, controls), Matrix(right.rows, right.cols, states), pivots


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

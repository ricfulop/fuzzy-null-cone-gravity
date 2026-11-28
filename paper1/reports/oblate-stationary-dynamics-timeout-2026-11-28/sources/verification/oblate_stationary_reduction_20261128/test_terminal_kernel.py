"""Small exact controls for the terminal kernel's reducer-facing return convention."""

from sympy import Float, Matrix, cancel, eye, zeros

from exact_linear_algebra import R, T


def matrices_equal(left, right):
    return left.shape == right.shape and all(cancel(value) == 0 for value in left - right)


def verify_controls(kernel_embedding):
    # The zero first row requires a swap. The final row repeats the rational constraint.
    a = (T + 1) / (R + 2)
    constraints = Matrix([[0, 0, 0], [R + 1, a, 0], [0, 0, T], [2*(R + 1), 2*a, 0]])
    embedding, free, basis, pivots = kernel_embedding(constraints)
    coeff = a / (R + 1)
    assert pivots == [0, 2] and free == [1]
    assert matrices_equal(basis, Matrix([[1, coeff, 0], [0, 0, 1]]))
    assert matrices_equal(embedding, Matrix([[-coeff], [1], [0]]))
    assert matrices_equal(constraints * embedding, zeros(4, 1))

    for nrows, ncols in ((0, 0), (0, 3), (2, 0), (2, 3)):
        embedding, free, basis, pivots = kernel_embedding(zeros(nrows, ncols))
        assert free == list(range(ncols)) and pivots == []
        assert matrices_equal(embedding, eye(ncols)) and basis.shape == (0, ncols)

    embedding, free, basis, pivots = kernel_embedding(eye(2))
    assert pivots == [0, 1] and free == []
    assert embedding.shape == (2, 0) and matrices_equal(basis, eye(2))

    try:
        kernel_embedding(Matrix([[Float("0.5")]]))
    except ValueError:
        pass
    else:
        raise AssertionError("Terminal kernel accepted an inexact coefficient")

    return {"state": "terminal kernel controls passed",
            "cases": ["dependent_rational_rows", "empty_dimensions", "zero_matrix",
                      "full_column_rank", "float_rejected"]}

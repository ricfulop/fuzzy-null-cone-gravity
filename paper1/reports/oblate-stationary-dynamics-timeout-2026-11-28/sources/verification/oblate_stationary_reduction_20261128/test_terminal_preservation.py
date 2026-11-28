"""Manufactured controls for terminal field identities and pre-preservation encoding."""

from sympy import Float, Matrix, eye, nan, oo, zeros, zoo

from exact_linear_algebra import R, T
from terminal_preservation import field_derivative, field_matrix, verify_terminal_preservation


def must_reject(operation, message):
    try:
        operation()
    except ValueError as exc:
        assert message in str(exc), str(exc)
    else:
        raise AssertionError(f"Expected rejection: {message}")


def verify_controls(snapshot_matrix):
    f = T/R
    df = -T/R**2 - 10/(R*(100*R**2 + 1))
    c = Matrix([[1, f]])
    k = Matrix([[1, f - df], [0, 0]])
    b = Matrix([[-f], [1]])
    m = b.copy()
    cf = field_matrix(c)
    dc = field_derivative(cf)
    assert dc == field_matrix(Matrix([[0, df]]))
    assert field_derivative(field_matrix(eye(2))).is_zero_matrix
    assert dc.add(cf.matmul(field_matrix(k))) == cf
    verify_terminal_preservation(c, k, b, m)

    bad_k = k + Matrix([[0, 1], [0, 0]])
    defect = dc.add(cf.matmul(field_matrix(bad_k))).matmul(field_matrix(m))
    assert defect == field_matrix(eye(1))
    must_reject(lambda: verify_terminal_preservation(c, bad_k, b, m),
                "constraint subspace is not invariant")
    bad_b = b + Matrix([[1], [0]])
    assert cf.matmul(field_matrix(bad_b)) == field_matrix(eye(1))
    must_reject(lambda: verify_terminal_preservation(c, k, bad_b, m),
                "remaining controls violate constraint preservation")

    verify_terminal_preservation(c, k, zeros(2, 0), m)
    verify_terminal_preservation(zeros(0, 2), k, b, eye(2))
    verify_terminal_preservation(zeros(0), zeros(0), zeros(0), zeros(0))
    for index in range(4):
        args = [matrix.copy() for matrix in (c, k, b, m)]
        args[index][0, 0] = Float("0.5")
        must_reject(lambda: verify_terminal_preservation(*args), "Inexact or nonfinite")

    assert snapshot_matrix(c) == {"shape": [1, 2], "entries": [["1", str(f)]]}
    assert snapshot_matrix(zeros(2, 0)) == {"shape": [2, 0], "entries": [[], []]}
    assert snapshot_matrix(zeros(0, 2)) == {"shape": [0, 2], "entries": []}
    for value in (Float("0.5"), nan, oo, -oo, zoo):
        must_reject(lambda: snapshot_matrix(Matrix([[value]])), "Inexact or nonfinite")

    return {"state": "terminal preservation controls passed",
            "cases": ["known_total_derivative", "nonzero_S_equals_C", "unit_SM_defect_rejected",
                      "unit_CB_defect_rejected", "empty_dimensions", "float_in_each_input_rejected",
                      "exact_snapshot_strings_and_shapes", "nonfinite_snapshot_rejected"]}

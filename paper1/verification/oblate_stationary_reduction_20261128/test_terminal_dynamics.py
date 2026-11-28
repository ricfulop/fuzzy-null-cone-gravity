"""Manufactured defects for the all-row checks and the full reconstruction."""

from sympy import Float, Matrix, eye, nan, oo, zeros, zoo

from exact_linear_algebra import R, T
from terminal_dynamics import reconstruct_dynamics, verify_original_rows
from terminal_preservation import field_matrix


def must_reject(operation, expected_message):
    try:
        operation()
    except ValueError as error:
        assert expected_message in str(error), str(error)
    else:
        raise AssertionError(f"Expected rejection: {expected_message}")


def check_original_rows():
    operators = [zeros(23, 7) for _ in range(3)]
    k, b, m = zeros(14), zeros(14, 1), eye(14)
    k[7, 0] = k[7, 7] = -1
    b[8, 0] = 1
    # On retained row 21, v0 + v0' + v0'' vanishes. Its control coefficient is zero.
    for block in operators:
        block[21, 0] = 1
    operators[2][22, 1] = 7
    assert (operators[2] * b[7:, :])[22, 0] == 7
    fields = tuple(map(field_matrix, (k, b, m)))
    verify_original_rows(operators, *fields)

    state_bad = [block.copy() for block in operators]
    state_bad[0][1, 0] = 1  # The X0 equation remains after the X0 column is removed.
    must_reject(lambda: verify_original_rows(state_bad, *fields), "do not vanish")

    control_bad = [block.copy() for block in operators]
    control_bad[2][2, 1] = 1  # Retained lambda1 row outside the eight discovery rows.
    must_reject(lambda: verify_original_rows(control_bad, *fields), "depend on remaining controls")

    for bad in (Float("0.5"), nan, oo, -oo, zoo):
        invalid = [block.copy() for block in operators]
        invalid[1][1, 0] = bad
        must_reject(lambda: verify_original_rows(invalid, *fields), "Inexact or nonfinite")


def check_reconstruction():
    f = T/R
    df = -T/R**2 - 10/(R*(100*R**2 + 1))
    m = Matrix([[1], [f]])
    b = m.copy()
    def reconstruct(k, control):
        return reconstruct_dynamics(*map(field_matrix, (k, control, m)), [0])

    for rate in (0, 2):
        k = Matrix([[rate, 0], [df + rate*f, 0]])
        reduced_k, reduced_b = reconstruct(k, b)
        assert reduced_k == Matrix([[rate]]) and reduced_b == eye(1)
        # Both defects occur outside free row 0; row selection alone would pass.
        must_reject(lambda: reconstruct(k + Matrix([[0, 0], [1, 0]]), b),
                    "K does not reconstruct")
        must_reject(lambda: reconstruct(k, b + Matrix([[0], [1]])),
                    "controls do not reconstruct")
    for index, bad in enumerate((oo, Float("0.5"), zoo)):
        args = [matrix.copy() for matrix in (k, b, m)]
        args[index][0, 0] = bad
        must_reject(lambda: reconstruct_dynamics(*map(field_matrix, args), [0]),
                    "Inexact or nonfinite")


def check_empty_widths():
    k, b = zeros(14), zeros(14, 0)
    for width in (0, 14):
        m = eye(14)[:, :width]
        fields = tuple(map(field_matrix, (k, b, m)))
        verify_original_rows([zeros(23, 7) for _ in range(3)], *fields)
        reduced_k, reduced_b = reconstruct_dynamics(*fields, list(range(width)))
        assert reduced_k.shape == (width, width) and reduced_b.shape == (width, 0)


def verify_controls():
    check_original_rows()
    check_reconstruction()
    check_empty_widths()
    return {"state": "terminal original-row and dynamics controls passed",
            "cases": ["all_22_valid_with_nonzero_discarded_control_response",
                      "X0_state_defect_rejected", "outside_discovery_control_defect_rejected",
                      "variable_embedding_with_valid_nonzero_control", "non_free_K_defect_rejected",
                      "non_free_B_defect_rejected", "empty_widths", "invalid_coefficients_rejected"]}

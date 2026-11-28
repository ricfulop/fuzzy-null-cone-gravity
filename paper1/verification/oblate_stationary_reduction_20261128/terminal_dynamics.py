"""Field identities for the retained Euler rows and reduced radial dynamics."""

from contextlib import nullcontext

from terminal_kernel import expression_matrix
from terminal_preservation import field_matrix, field_derivative


def verify_original_rows(operators, k, b, m, *, stage=nullcontext):
    """Use all 22 retained rows and both acceleration blocks. Row 22 is the response."""
    with stage("original_rows_field_blocks"):
        retained = range(22)
        if len(operators) != 3 or any(op.shape != (23, 7) for op in operators):
            raise ValueError("Expected three gauge-fixed 23 by 7 operator matrices")
        e0, e1, e2 = (field_matrix(op.extract(retained, range(7))) for op in operators)
        ka = k.extract(range(7, 14), range(k.shape[1]))
        ba = b.extract(range(7, 14), range(b.shape[1]))
    with stage("original_rows_acceleration_state"):
        state_rows = e0.hstack(e1).add(e2.matmul(ka))
    with stage("original_rows_state_product"):
        if not state_rows.matmul(m).is_zero_matrix:
            raise ValueError("original retained rows do not vanish on the reduced state embedding")
    with stage("original_rows_control_product"):
        if not e2.matmul(ba).is_zero_matrix:
            raise ValueError("original retained rows depend on remaining controls")


def reconstruct_dynamics(k, b, m, free_rows, *, stage=nullcontext):
    """Extract the dynamics on free rows, then check the two lifts on every row."""
    with stage("terminal_reduced_dynamics"):
        with stage("reduced_KM"):
            km = k.matmul(m)
        with stage("reduced_DM"):
            dm = field_derivative(m)
        with stage("reduced_difference"):
            w = km.sub(dm)
        with stage("reduced_free_rows"):
            reduced_k = w.extract(free_rows, range(w.shape[1]))
            reduced_b = b.extract(free_rows, range(b.shape[1]))

    with stage("terminal_reconstruct_dynamics"):
        with stage("reconstruct_state_product"):
            if not m.matmul(reduced_k).sub(w).is_zero_matrix:
                raise ValueError("free-coordinate K does not reconstruct the embedded dynamics")
        with stage("reconstruct_control_product"):
            if not m.matmul(reduced_b).sub(b).is_zero_matrix:
                raise ValueError("free-coordinate controls do not reconstruct embedded controls")
        with stage("reduced_to_expressions"):
            return expression_matrix(reduced_k), expression_matrix(reduced_b)

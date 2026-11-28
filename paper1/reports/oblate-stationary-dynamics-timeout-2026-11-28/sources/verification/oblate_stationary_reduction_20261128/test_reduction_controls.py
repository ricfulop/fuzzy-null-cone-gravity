"""Manufactured controls for the stationary DAE reducer.

These tests use synthetic 23-row, 7-column gauge-fixed systems.  They exercise
the DAE machinery without importing the stationary extractor or recomputing any
Einstein tensors.
"""

from contextlib import contextmanager

import sympy as s


class MemoryRecorder:
    """Stage/write adapter for reducer tests that must not touch the filesystem."""

    def __init__(self):
        self.stages = []
        self.outputs = {}

    @contextmanager
    def stage(self, name, **details):
        self.stages.append((name, details))
        yield

    def write_json(self, name, value):
        self.outputs[name] = value


def blank_system():
    return [s.zeros(23, 7) for _ in range(3)]


def exact_zero_matrix(red, matrix):
    return all(red.exact_zero(value) for value in matrix)


def assert_retained_rows_vanish(red, gauge_matrices, closure):
    retained_state, retained_control = red.retained_parts(
        gauge_matrices, closure["state"]
    )
    assert exact_zero_matrix(red, retained_state * closure["embedding"])
    assert exact_zero_matrix(red, retained_control)


def run_closure(red, gauge_matrices):
    recorder = MemoryRecorder()
    closure = red.reduce_dae(
        gauge_matrices, recorder, max_steps=red.STATE_COUNT + red.CONTROL_COUNT
    )
    assert recorder.outputs
    return closure, recorder


def check_hidden_constraints_keep_free_response_control(red):
    matrices = blank_system()
    equation_v0_second = red.EXPECTED_SELECTED[0]
    equation_v0_zero = red.EXPECTED_SELECTED[1]

    # v0'' + v1 = 0.
    matrices[2][equation_v0_second, 0] = 1
    matrices[0][equation_v0_second, 1] = 1

    # v0 = 0.  Preservation gives v0'=0, then v1=0, then v1'=0.
    matrices[0][equation_v0_zero, 0] = 1

    # D = v2'' should remain a free response control.
    matrices[2][red.RESPONSE_ROW, 2] = 1

    closure, recorder = run_closure(red, matrices)

    assert closure["state"].constraints.rows == 4
    assert closure["embedding"].cols == 10
    assert len(closure["state"].control_labels) == 5
    assert_retained_rows_vanish(red, matrices, closure)
    assert exact_zero_matrix(red, closure["response_state"])
    assert closure["response_control"] == s.Matrix([[1, 0, 0, 0, 0]])

    return {
        "name": "hidden_constraints_keep_free_response_control",
        "constraint_count": closure["state"].constraints.rows,
        "free_state_count": closure["embedding"].cols,
        "remaining_control_count": len(closure["state"].control_labels),
        "response_state": red.as_strings(closure["response_state"]),
        "response_control": red.as_strings(closure["response_control"]),
        "recorded_stage_count": len(recorder.stages),
        "recorded_output_count": len(recorder.outputs),
    }


def check_variable_constraint_preservation_kills_response(red):
    matrices = blank_system()
    equation_constraint = red.EXPECTED_SELECTED[0]
    equation_v1_second = red.EXPECTED_SELECTED[1]

    # v0 - r v1 = 0.
    matrices[0][equation_constraint, 0] = 1
    matrices[0][equation_constraint, 1] = -red.R

    # v1'' = 0.
    matrices[2][equation_v1_second, 1] = 1

    # The differentiated constraint chain forces v0'' = 2 v1',
    # so this response must reduce to zero exactly.
    matrices[2][red.RESPONSE_ROW, 0] = 1
    matrices[1][red.RESPONSE_ROW, 1] = -2

    closure, recorder = run_closure(red, matrices)

    assert closure["state"].constraints.rows == 2
    assert closure["embedding"].cols == 12
    assert len(closure["state"].control_labels) == 5
    assert_retained_rows_vanish(red, matrices, closure)
    assert exact_zero_matrix(red, closure["response_state"])
    assert exact_zero_matrix(red, closure["response_control"])

    return {
        "name": "variable_constraint_preservation_kills_response",
        "constraint_count": closure["state"].constraints.rows,
        "free_state_count": closure["embedding"].cols,
        "remaining_control_count": len(closure["state"].control_labels),
        "response_state": red.as_strings(closure["response_state"]),
        "response_control": red.as_strings(closure["response_control"]),
        "recorded_stage_count": len(recorder.stages),
        "recorded_output_count": len(recorder.outputs),
    }


def verify_controls(reducer_module):
    cases = [
        check_hidden_constraints_keep_free_response_control(reducer_module),
        check_variable_constraint_preservation_kills_response(reducer_module),
    ]
    return {
        "state": "manufactured DAE controls passed",
        "case_count": len(cases),
        "cases": cases,
    }


if __name__ == "__main__":
    import reduce_stationary as red

    verify_controls(red)

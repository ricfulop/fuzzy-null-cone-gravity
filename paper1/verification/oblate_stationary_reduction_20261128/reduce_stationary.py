"""Exact JSON-only reduction of the exported oblate stationary system.

This script consumes only a completed stationary-operator JSON export.  It does
not import the extractor, recompute curvature, choose a radial boundary model,
or certify a nonlinear obstruction.
"""

from dataclasses import dataclass
from pathlib import Path
import argparse
import json
import os
import platform
import sys

import sympy as s

import reduction_progress as progress


EPSILON0 = s.Rational(1, 100)
R, EPSILON, T = s.symbols("r epsilon T", positive=True)
ATAN_FIXED = s.atan(s.Rational(1, 10) / R)
T_PRIME = -s.Rational(10, 1) / (100 * R**2 + 1)

EXPECTED_COLUMNS = [
    ["lambda", 0],
    ["lambda", 2],
    ["X", 0],
    ["X", 2],
    ["V_polar", 2],
    ["C_polar", 2],
    ["a", 0],
    ["a", 2],
]
EXPECTED_ROWS = [
    ["lambda", 0],
    ["X", 0],
    ["lambda", 1],
    ["X", 1],
    ["V_polar", 1],
    ["V_axial", 1],
    ["lambda", 2],
    ["X", 2],
    ["V_polar", 2],
    ["V_axial", 2],
    ["C_polar", 2],
    ["C_axial", 2],
    ["a", 0],
    ["a", 1],
    ["a", 2],
    ["b", 0],
    ["b", 1],
    ["W_polar", 1],
    ["W_axial", 1],
    ["b", 2],
    ["W_polar", 2],
    ["W_axial", 2],
]
EXPECTED_SELECTED = [0, 6, 1, 7, 8, 10, 12, 14]
EXPECTED_AUTOMATIC = [2, 3, 4, 5, 9, 11, 13, 15, 16, 17, 18, 19, 20, 21]
X0_COLUMN = 2
RESPONSE_ROW = 22
STATE_COUNT = 14
CONTROL_COUNT = 7
EXPECTED_COMPLETION_SHA256 = "640ae2f7c9e1cd3702f8330b074df361ee3fdff701ee2a566ef2b84c4e3a6b4d"


SAFE_LOCALS = {
    "r": R,
    "epsilon": EPSILON,
    "atan": s.atan,
    "sqrt": s.sqrt,
}


@dataclass
class ReductionState:
    K: s.Matrix
    B: s.Matrix
    control_labels: list
    constraints: s.Matrix


def exact_zero(expr):
    return s.cancel(s.together(expr)) == 0


def exact_matrix(matrix):
    return matrix.applyfunc(lambda value: s.cancel(s.together(value)))


def field_diff(expr, order=1):
    result = expr
    for _ in range(order):
        result = s.diff(result, R) + T_PRIME * s.diff(result, T)
        result = s.cancel(s.together(result))
    return result


def field_diff_matrix(matrix):
    return matrix.applyfunc(field_diff)


def as_strings(matrix):
    return [[str(s.cancel(s.together(matrix[i, j]))) for j in range(matrix.cols)]
            for i in range(matrix.rows)]


def matrix_from_rows(rows, columns=STATE_COUNT):
    if not rows:
        return s.zeros(0, columns)
    return s.Matrix(rows)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def load_completed_export(input_dir):
    input_dir = Path(input_dir)
    completion_path = input_dir / "completion.json"
    operators_path = input_dir / "stationary-operators.json"
    reproducibility_path = input_dir / "stationary-reproducibility.json"
    for path in (completion_path, operators_path, reproducibility_path):
        require(path.exists(), f"missing required completed-export file: {path}")
    require(progress.sha256_file(completion_path) == EXPECTED_COMPLETION_SHA256,
            "completion marker differs from the independently reviewed upstream export")

    completion = load_json(completion_path)
    operators = load_json(operators_path)
    reproducibility = load_json(reproducibility_path)

    require(
        completion.get("all_checker_controls_passed") is True,
        "completion marker does not assert passed controls",
    )
    require(completion.get("column_count") == 8, "unexpected column count")
    require(completion.get("retained_row_count") == 22, "unexpected row count")
    require(
        completion.get("source_sha256_before") == completion.get("source_sha256_after"),
        "source hashes changed during upstream export",
    )
    exports = completion.get("export_sha256", {})
    require(
        exports.get("stationary-operators.json")
        == progress.sha256_file(operators_path),
        "stationary-operators.json hash does not match completion marker",
    )
    require(
        exports.get("stationary-reproducibility.json")
        == progress.sha256_file(reproducibility_path),
        "stationary-reproducibility.json hash does not match completion marker",
    )
    require(
        operators.get("radial_DAE_reduced") is False,
        "input is not the raw unreduced operator export",
    )
    require(
        operators.get("response_evaluated_on_admissible_non_gauge_seed") is False,
        "unexpected upstream response evaluation",
    )
    require(operators.get("column_labels") == EXPECTED_COLUMNS, "column labels changed")
    require(
        operators.get("retained_row_labels") == EXPECTED_ROWS,
        "retained row labels changed",
    )
    require(
        operators.get("selected_row_indices") == EXPECTED_SELECTED,
        "selected row ordering changed",
    )
    require(
        operators.get("automatic_zero_row_indices") == EXPECTED_AUTOMATIC,
        "automatic row set changed",
    )
    require(operators.get("last_row") ==
            "discarded response D=6 E_X4-E_Cpolar4; not a retained equation",
            "discarded response convention changed")
    require(operators.get("operator_convention") ==
            "sum_j M_j(r,epsilon) u^(j)(r); no derivative acts on M_j here",
            "radial derivative convention changed")

    matrices = operators.get("matrices_by_derivative_order")
    require(isinstance(matrices, list) and len(matrices) == 3, "expected 3 matrices")
    for order, matrix in enumerate(matrices):
        require(len(matrix) == 23, f"order {order} row count changed")
        for row in matrix:
            require(len(row) == 8, f"order {order} column count changed")

    return {
        "input_dir": str(input_dir),
        "completion": completion,
        "operators": operators,
        "reproducibility": reproducibility,
        "input_sha256": {
            "completion.json": progress.sha256_file(completion_path),
            "stationary-operators.json": progress.sha256_file(operators_path),
            "stationary-reproducibility.json": progress.sha256_file(
                reproducibility_path
            ),
        },
    }


def parse_expression(text):
    return s.sympify(text, locals=SAFE_LOCALS)


def convert_expression(expr):
    fixed = s.cancel(s.together(expr.subs(EPSILON, EPSILON0)))
    converted = fixed.xreplace({ATAN_FIXED: T})
    require(not converted.has(EPSILON), "epsilon survived fixed-epsilon conversion")
    require(not converted.has(s.atan), "atan survived Q(r,T) conversion")
    require(not converted.has(s.sqrt), "sqrt survived Q(r,T) conversion")
    require(converted.free_symbols <= {R, T}, "unexpected symbol in expression")
    require(not converted.atoms(s.Float), "inexact coefficient in formal field")
    numerator, denominator = s.fraction(s.cancel(converted))
    s.Poly(numerator, R, T, domain=s.QQ)
    denominator_poly = s.Poly(denominator, R, T, domain=s.QQ)
    require(not denominator_poly.is_zero, "zero formal denominator")
    recovered = converted.xreplace({T: ATAN_FIXED})
    require(
        exact_zero(recovered - fixed),
        "Q(r,T) conversion failed its reversible substitution check",
    )
    return s.cancel(s.together(converted))


def convert_matrices(matrix_strings):
    converted = []
    for order in range(3):
        matrix = s.zeros(23, 8)
        for row in range(23):
            for column in range(8):
                matrix[row, column] = convert_expression(
                    parse_expression(matrix_strings[order][row][column])
                )
        converted.append(exact_matrix(matrix))
    return converted


def apply_operator(matrices, profiles):
    profile_matrix = s.Matrix(profiles)
    result = s.zeros(matrices[0].rows, 1)
    for order, matrix in enumerate(matrices):
        result += matrix * profile_matrix.applyfunc(
            lambda value, order=order: field_diff(value, order)
        )
    return exact_matrix(result)


def verify_radial_gauge(matrices):
    f = s.Function("f")(R)
    rho = R**2 + EPSILON0
    lambda0 = (R**2 + EPSILON0 / 3) / rho
    lambda2 = 2 * EPSILON0 / (3 * rho)
    profiles = [
        field_diff(lambda0) * f + 2 * lambda0 * field_diff(f),
        field_diff(lambda2) * f + 2 * lambda2 * field_diff(f),
        2 * R * f,
        0,
        0,
        0,
        0,
        0,
    ]
    residual = apply_operator(matrices, profiles)
    require(
        all(exact_zero(value) for value in residual),
        "radial X0 gauge profile is not a right-kernel",
    )
    return residual


def remove_x0_column(matrices):
    kept_columns = [i for i in range(8) if i != X0_COLUMN]
    return [matrix.extract(range(matrix.rows), kept_columns) for matrix in matrices]


def rows(matrix, indices):
    return matrix.extract(indices, range(matrix.cols))


def zero_row(row):
    return all(exact_zero(row[0, column]) for column in range(row.cols))


def rref_with_side(left, right):
    left = exact_matrix(s.Matrix(left))
    right = exact_matrix(s.Matrix(right))
    pivot_row = 0
    pivots = []
    for column in range(left.cols):
        found = None
        for row in range(pivot_row, left.rows):
            if not exact_zero(left[row, column]):
                found = row
                break
        if found is None:
            continue
        if found != pivot_row:
            left.row_swap(found, pivot_row)
            right.row_swap(found, pivot_row)

        pivot = left[pivot_row, column]
        left.row_op(pivot_row, lambda value, _: s.cancel(s.together(value / pivot)))
        right.row_op(pivot_row, lambda value, _: s.cancel(s.together(value / pivot)))

        for row in range(left.rows):
            if row == pivot_row:
                continue
            factor = left[row, column]
            if exact_zero(factor):
                continue
            left.row_op(
                row,
                lambda value, j, factor=factor, pivot_row=pivot_row:
                    s.cancel(s.together(value - factor * left[pivot_row, j])),
            )
            right.row_op(
                row,
                lambda value, j, factor=factor, pivot_row=pivot_row:
                    s.cancel(s.together(value - factor * right[pivot_row, j])),
            )
        pivots.append(column)
        pivot_row += 1
        if pivot_row == left.rows:
            break
    return exact_matrix(left), exact_matrix(right), pivots


def row_basis(matrix):
    matrix = s.Matrix(matrix)
    if matrix.rows == 0:
        return s.zeros(0, matrix.cols), []
    reduced, _, pivots = rref_with_side(matrix, s.zeros(matrix.rows, 0))
    basis = [reduced.row(i) for i in range(reduced.rows) if not zero_row(reduced.row(i))]
    return matrix_from_rows(basis, matrix.cols), pivots


def matrix_rank(matrix):
    return row_basis(matrix)[0].rows


def independent_extension(existing, candidates):
    existing_basis, _ = row_basis(existing)
    accepted = []
    rank = existing_basis.rows
    for row_index in range(candidates.rows):
        row = candidates.row(row_index)
        if zero_row(row):
            continue
        trial = existing_basis.col_join(matrix_from_rows(accepted)).col_join(row)
        trial_rank = matrix_rank(trial)
        if trial_rank > rank + len(accepted):
            accepted.append(row)
    if accepted:
        return matrix_from_rows(accepted, existing.cols)
    return s.zeros(0, existing.cols)


def kernel_embedding(constraints):
    state_count = constraints.cols
    basis, pivots = row_basis(constraints)
    free_rows = [index for index in range(state_count) if index not in pivots]
    embedding = s.zeros(state_count, len(free_rows))
    for free_index, state_index in enumerate(free_rows):
        embedding[state_index, free_index] = 1
    for row, pivot in enumerate(pivots):
        for free_index, state_index in enumerate(free_rows):
            embedding[pivot, free_index] = -basis[row, state_index]
    require(
        all(exact_zero(value) for value in exact_matrix(constraints * embedding)),
        "constraint kernel embedding failed",
    )
    return exact_matrix(embedding), free_rows, basis, pivots


def initial_state():
    K = s.zeros(STATE_COUNT, STATE_COUNT)
    for index in range(CONTROL_COUNT):
        K[index, CONTROL_COUNT + index] = 1
    B = s.zeros(STATE_COUNT, CONTROL_COUNT)
    for index in range(CONTROL_COUNT):
        B[CONTROL_COUNT + index, index] = 1
    controls = [f"d2_{label[0]}_{label[1]}" for label in EXPECTED_COLUMNS
                if label != ["X", 0]]
    return ReductionState(K=K, B=B, control_labels=controls,
                          constraints=s.zeros(0, STATE_COUNT))


def split_selected_system(gauge_matrices):
    selected = [rows(matrix, EXPECTED_SELECTED) for matrix in gauge_matrices]
    acceleration = selected[2]
    state = selected[0].row_join(selected[1])
    return acceleration, state


def solve_controls(state, acceleration, state_matrix):
    controls_in = list(state.control_labels)
    reduced_acc, reduced_state, pivots = rref_with_side(acceleration, state_matrix)
    free = [index for index in range(acceleration.cols) if index not in pivots]

    substitution_state = s.zeros(acceleration.cols, STATE_COUNT)
    substitution_controls = s.zeros(acceleration.cols, len(free))
    for row, pivot in enumerate(pivots):
        substitution_state[pivot, :] = -reduced_state.row(row)
        for free_index, free_column in enumerate(free):
            substitution_controls[pivot, free_index] = -reduced_acc[row, free_column]
    for free_index, free_column in enumerate(free):
        substitution_controls[free_column, free_index] = 1

    residuals = []
    for row in range(len(pivots), reduced_acc.rows):
        require(
            zero_row(reduced_acc.row(row)),
            "control elimination left a nonsolved acceleration row",
        )
        if not zero_row(reduced_state.row(row)):
            residuals.append(reduced_state.row(row))

    return {
        "pivots": pivots,
        "pivot_labels": [controls_in[index] for index in pivots],
        "free_labels": [controls_in[index] for index in free],
        "substitution_state": exact_matrix(substitution_state),
        "substitution_controls": exact_matrix(substitution_controls),
        "residual_constraints": matrix_from_rows(residuals, STATE_COUNT),
    }


def apply_control_solution(state, solution):
    return ReductionState(
        K=exact_matrix(state.K + state.B * solution["substitution_state"]),
        B=exact_matrix(state.B * solution["substitution_controls"]),
        control_labels=list(solution["free_labels"]),
        constraints=state.constraints,
    )


def preservation_equations(constraints, state):
    return (
        exact_matrix(constraints * state.B),
        exact_matrix(field_diff_matrix(constraints) + constraints * state.K),
    )


def response_parts(gauge_matrices, state):
    bottom_rows = list(range(CONTROL_COUNT, STATE_COUNT))
    acceleration_state = state.K.extract(bottom_rows, range(STATE_COUNT))
    acceleration_controls = state.B.extract(bottom_rows, range(state.B.cols))
    response0 = gauge_matrices[0].row(RESPONSE_ROW)
    response1 = gauge_matrices[1].row(RESPONSE_ROW)
    response2 = gauge_matrices[2].row(RESPONSE_ROW)
    return (
        exact_matrix(response0.row_join(response1) + response2 * acceleration_state),
        exact_matrix(response2 * acceleration_controls),
    )


def retained_parts(gauge_matrices, state):
    bottom_rows = list(range(CONTROL_COUNT, STATE_COUNT))
    acceleration_state = state.K.extract(bottom_rows, range(STATE_COUNT))
    acceleration_controls = state.B.extract(bottom_rows, range(state.B.cols))
    retained0 = gauge_matrices[0].extract(range(22), range(CONTROL_COUNT))
    retained1 = gauge_matrices[1].extract(range(22), range(CONTROL_COUNT))
    retained2 = gauge_matrices[2].extract(range(22), range(CONTROL_COUNT))
    return (
        exact_matrix(retained0.row_join(retained1) + retained2 * acceleration_state),
        exact_matrix(retained2 * acceleration_controls),
    )


def reduce_dae(gauge_matrices, recorder, max_steps):
    acceleration, state_matrix = split_selected_system(gauge_matrices)
    state = initial_state()
    pending_acceleration = acceleration
    pending_state = state_matrix
    steps = []

    for step in range(max_steps):
        with recorder.stage("dae_step", step=step, control_count=state.B.cols):
            solution = solve_controls(state, pending_acceleration, pending_state)
            state = apply_control_solution(state, solution)
            new_constraints = independent_extension(
                state.constraints, solution["residual_constraints"]
            )
            if new_constraints.rows:
                state = ReductionState(
                    K=state.K,
                    B=state.B,
                    control_labels=state.control_labels,
                    constraints=state.constraints.col_join(new_constraints),
                )
            steps.append({
                "step": step,
                "pivot_labels": solution["pivot_labels"],
                "free_labels": solution["free_labels"],
                "new_constraint_count": new_constraints.rows,
                "total_constraint_count": state.constraints.rows,
                "remaining_control_count": state.B.cols,
            })
            recorder.write_json(
                f"dae-step-{step:02d}.json",
                {
                    "state": "intermediate DAE elimination stage; not a seed certificate",
                    **steps[-1],
                    "substitution_state": as_strings(solution["substitution_state"]),
                    "substitution_controls": as_strings(
                        solution["substitution_controls"]
                    ),
                    "constraints": as_strings(state.constraints),
                },
            )

            if new_constraints.rows == 0:
                break
            pending_acceleration, pending_state = preservation_equations(
                new_constraints, state
            )
    else:
        raise RuntimeError("DAE prolongation did not close within the step bound")

    embedding, free_rows, constraints_rref, constraint_pivots = kernel_embedding(
        state.constraints
    )
    preservation_acc, preservation_state = preservation_equations(
        state.constraints, state
    )
    require(
        all(exact_zero(value) for value in preservation_acc),
        "remaining controls violate constraint preservation",
    )
    require(
        all(exact_zero(value) for value in exact_matrix(preservation_state * embedding)),
        "constraint subspace is not invariant",
    )

    retained_state, retained_control = retained_parts(gauge_matrices, state)
    require(
        all(exact_zero(value) for value in exact_matrix(retained_state * embedding)),
        "original retained rows do not vanish on the reduced state embedding",
    )
    require(
        all(exact_zero(value) for value in retained_control),
        "original retained rows depend on remaining controls",
    )

    km_minus_mprime = exact_matrix(state.K * embedding - field_diff_matrix(embedding))
    reduced_K = exact_matrix(km_minus_mprime.extract(free_rows, range(embedding.cols)))
    reduced_B = exact_matrix(state.B.extract(free_rows, range(state.B.cols)))
    require(
        all(
            exact_zero(value)
            for value in exact_matrix(
                field_diff_matrix(embedding) + embedding * reduced_K
                - state.K * embedding
            )
        ),
        "free-coordinate K does not reconstruct the embedded dynamics",
    )
    require(
        all(
            exact_zero(value)
            for value in exact_matrix(embedding * reduced_B - state.B)
        ),
        "free-coordinate controls do not reconstruct embedded controls",
    )

    response_state, response_control = response_parts(gauge_matrices, state)
    return {
        "state": state,
        "steps": steps,
        "constraints_rref": constraints_rref,
        "constraint_pivots": constraint_pivots,
        "embedding": embedding,
        "free_rows": free_rows,
        "reduced_K": reduced_K,
        "reduced_B": reduced_B,
        "response_state": exact_matrix(response_state * embedding),
        "response_control": response_control,
    }


def run(input_dir, output_dir, stage, max_steps):
    if not __debug__:
        raise RuntimeError("Scientific assertions must remain enabled")
    progress.report(phase="setup", step=1, of=1, unit="environment")
    export = load_completed_export(input_dir)
    source_paths = sorted(Path(__file__).resolve().parent.glob("*.py"))
    source_hashes = {p.name: progress.sha256_file(p) for p in source_paths}
    identity = {
        "python": platform.python_version(),
        "python_implementation": sys.version,
        "sympy": s.__version__,
        "epsilon0": "1/100",
        "T": "atan(1/(10*r))",
        "T_prime": str(T_PRIME),
        "input_dir": export["input_dir"],
        "input_sha256": export["input_sha256"],
        "upstream_source_commit": export["reproducibility"].get("source_commit"),
        "upstream_identity_sha256": export["completion"].get("identity_sha256"),
        "scope": "JSON-only fixed-epsilon stationary DAE reduction",
        "source_commit": os.environ.get("ORX_COMMIT"),
        "source_sha256": source_hashes,
        "stage": stage,
        "max_steps": max_steps,
    }
    with progress.ReductionRecorder(output_dir, identity) as recorder:
        with recorder.stage("manufactured_controls"):
            from test_reduction_controls import verify_controls
            controls = verify_controls(sys.modules[__name__])
            recorder.write_json("manufactured-controls.json", controls)
        progress.report(m={"manufactured_controls_passed": 1})
        with recorder.stage("fixed_epsilon_field_conversion"):
            matrices = convert_matrices(
                export["operators"]["matrices_by_derivative_order"]
            )
        with recorder.stage("radial_x0_gauge_right_kernel"):
            gauge_residual = verify_radial_gauge(matrices)
        with recorder.stage("principal_block_probe"):
            gauge_matrices = remove_x0_column(matrices)
            acceleration, state_matrix = split_selected_system(gauge_matrices)
            _, _, principal_pivots = rref_with_side(
                acceleration, s.zeros(acceleration.rows, 0)
            )
            recorder.write_json(
                "principal-block.json",
                {
                    "state": "principal rank probe; not a DAE closure",
                    "epsilon0": "1/100",
                    "T": "atan(1/(10*r))",
                    "selected_row_indices": EXPECTED_SELECTED,
                    "removed_column": ["X", 0],
                    "kept_columns": [
                        label for label in EXPECTED_COLUMNS if label != ["X", 0]
                    ],
                    "principal_rank": len(principal_pivots),
                    "principal_pivots": principal_pivots,
                    "principal_block": as_strings(acceleration),
                    "state_block": as_strings(state_matrix),
                    "radial_gauge_residual_zero": all(
                        exact_zero(value) for value in gauge_residual
                    ),
                    "smooth_solution_germ_certified": False,
                    "response_evaluated_on_admissible_non_gauge_seed": False,
                },
            )
        progress.report(m={"principal_rank": len(principal_pivots), "gauge_verified": 1})

        if stage == "principal":
            progress.report(phase="export", step=0, of=1, unit="bundle")
            require(source_hashes == {p.name: progress.sha256_file(p) for p in source_paths},
                    "reducer sources changed during the invocation")
            recorder.write_json(
                "summary.json",
                {
                    "state": "principal probe complete; DAE closure not requested",
                    "smooth_solution_germ_certified": False,
                    "admissible_non_gauge_seed_supplied": False,
                    "nonlinear_consistency_established": False,
                },
            )
            progress.report(phase="export", step=1, of=1, unit="bundle", m={"principal_probe_passed": 1})
            print("PASS: exact principal block and radial gauge; no DAE closure or solution germ.", flush=True)
            return

        with recorder.stage("dae_closure"):
            closure = reduce_dae(gauge_matrices, recorder, max_steps)
        progress.report(phase="export", step=0, of=1, unit="bundle")
        require(source_hashes == {p.name: progress.sha256_file(p) for p in source_paths},
                "reducer sources changed during the invocation")
        with recorder.stage("closure_export"):
            recorder.write_json(
                "dae-closure.json",
                {
                    "state": (
                        "formal DAE closure over Q(r,T) at fixed epsilon; "
                        "no interval pivot certificate or responding seed"
                    ),
                    "constraint_count": closure["state"].constraints.rows,
                    "remaining_control_labels": closure["state"].control_labels,
                    "free_rows": closure["free_rows"],
                    "constraint_pivots": closure["constraint_pivots"],
                    "constraints_rref": as_strings(closure["constraints_rref"]),
                    "embedding": as_strings(closure["embedding"]),
                    "reduced_K": as_strings(closure["reduced_K"]),
                    "reduced_B": as_strings(closure["reduced_B"]),
                    "response_state": as_strings(closure["response_state"]),
                    "response_control": as_strings(closure["response_control"]),
                    "constant_control_augmented_state_supported": True,
                    "smooth_solution_germ_certified": False,
                    "admissible_non_gauge_seed_supplied": False,
                    "nonlinear_consistency_established": False,
                },
            )
            result = {
                "formal_constraint_count": closure["state"].constraints.rows,
                "formal_free_state_count": closure["embedding"].cols,
                "formal_free_control_count": closure["state"].B.cols,
                "formal_response_nonzero": int(any(not exact_zero(v) for v in
                    list(closure["response_state"]) + list(closure["response_control"]))),
                "dae_closure_passed": 1,
            }
            recorder.write_json("summary.json", {
                "state": "formal DAE closure; actual regular interval not certified",
                **result,
                "smooth_solution_germ_certified": False,
                "nonlinear_consistency_established": False,
            })
        progress.report(phase="export", step=1, of=1, unit="bundle", m=result)
        print("PASS: formal DAE closure and all 22 original rows; no actual solution germ yet.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-dir",
        type=Path,
        required=True,
        help="Completed stationary export outputs directory",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Fresh reduction output directory",
    )
    parser.add_argument(
        "--stage",
        choices=["principal", "closure"],
        default="principal",
        help="Stop after the principal probe or attempt full DAE closure",
    )
    parser.add_argument(
        "--max-steps",
        type=int,
        default=STATE_COUNT + CONTROL_COUNT,
        help="Bound on DAE preservation steps",
    )
    args = parser.parse_args()
    run(args.input_dir, args.output_dir, args.stage, args.max_steps)


if __name__ == "__main__":
    main()

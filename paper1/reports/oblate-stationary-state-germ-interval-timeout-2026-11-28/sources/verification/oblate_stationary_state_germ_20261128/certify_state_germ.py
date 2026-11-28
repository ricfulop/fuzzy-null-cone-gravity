"""Direct 14-state germ certificate for the oblate stationary system.

This JSON-only checker consumes:

* the original completed stationary-operator export;
* the terminal state snapshot saved before preservation; and
* a config that pins both input SHA-256 digests.

It does not import the extractor, the reducer, or the previous germ certifier.
It proves a local stationary germ directly in the fourteen state coordinates by
checking the terminal constraint kernel, the invariant constraint identity, the
physical lift v'=w, all 22 retained Euler rows, and a same-radius interval
certificate for regularity, rank six, and any discarded response witness.
"""

import argparse
from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import sys
import tempfile
import time

import sympy as s
from sympy.polys.matrices import DomainMatrix


EPSILON0 = s.Rational(1, 100)
R, EPSILON, T = s.symbols("r epsilon T", positive=True)
FIELD = s.QQ.frac_field(R, T)
FIELD_ZERO = FIELD.from_sympy(s.Integer(0))
FIELD_ONE = FIELD.from_sympy(s.Integer(1))
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
ORIGINAL_COLUMN_COUNT = 8
PROFILE_COUNT = ORIGINAL_COLUMN_COUNT - 1
RETAINED_ROW_COUNT = 22
RESPONSE_ROW = 22
STATE_COUNT = 2 * PROFILE_COUNT
FREE_STATE_COUNT = 8

KEPT_ORIGINAL_COLUMNS = [
    index for index in range(ORIGINAL_COLUMN_COUNT) if index != X0_COLUMN
]
EXPECTED_STATE_LABELS = [
    {"derivative_order": order, "coefficient_label": label}
    for order in (0, 1)
    for label in EXPECTED_COLUMNS
    if label != ["X", 0]
]

SAFE_RAW_LOCALS = {
    "r": R,
    "epsilon": EPSILON,
    "atan": s.atan,
    "sqrt": s.sqrt,
}
SAFE_Q_LOCALS = {"r": R, "T": T}

DEFAULT_R0S = [
    s.Rational(1),
    s.Rational(2),
    s.Rational(1, 2),
    s.Rational(3),
    s.Rational(1, 3),
    s.Rational(3, 2),
    s.Rational(2, 3),
    s.Rational(4, 3),
    s.Rational(3, 4),
    s.Rational(5, 4),
    s.Rational(4, 5),
]


class NullProgress:
    @contextmanager
    def stage(self, _name, **_details):
        yield

    def metric(self, **_metrics):
        pass


class StageProgress:
    """Small flushed progress reporter copied locally from certify_germ.py."""

    def __init__(self):
        self.started = time.monotonic()
        self.sequence = 0

    def report(self, **fields):
        print("::orx:: " + json.dumps(fields, sort_keys=True, allow_nan=False), flush=True)

    def event(self, event, **fields):
        self.sequence += 1
        record = {
            "sequence": self.sequence,
            "event": event,
            "seconds": round(time.monotonic() - self.started, 6),
            **fields,
        }
        print(
            "PROGRESS " + json.dumps(record, sort_keys=True, allow_nan=False),
            file=sys.stderr,
            flush=True,
        )

    @contextmanager
    def stage(self, name, unit="check", **details):
        self.report(phase=name, step=0, of=1, unit=unit)
        started = time.monotonic()
        self.event("begin", stage=name, **details)
        try:
            yield
        except BaseException as error:
            self.event(
                "error",
                stage=name,
                duration_seconds=round(time.monotonic() - started, 6),
                error_type=type(error).__name__,
                error=str(error),
                **details,
            )
            raise
        else:
            self.event(
                "end",
                stage=name,
                duration_seconds=round(time.monotonic() - started, 6),
                **details,
            )
            self.report(phase=name, step=1, of=1, unit=unit)

    def metric(self, **metrics):
        self.event("metric", **metrics)
        numeric = {}
        for key, value in metrics.items():
            if isinstance(value, bool):
                numeric[key] = int(value)
            elif isinstance(value, (int, float)):
                numeric[key] = value
        if numeric:
            self.report(m=numeric)


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def load_json_bytes(path):
    data = Path(path).read_bytes()
    return json.loads(data.decode("utf-8")), data, sha256_bytes(data)


def load_pinned_json(path, expected_sha256, label):
    data = Path(path).read_bytes()
    actual_sha256 = sha256_bytes(data)
    require(actual_sha256 == expected_sha256, f"{label} SHA-256 mismatch")
    return json.loads(data.decode("utf-8")), data, actual_sha256


def canonical_sha256(value):
    return sha256_bytes(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    )


def atomic_json_new(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            dir=path.parent,
            prefix=".partial-",
            delete=False,
        ) as stream:
            temporary = Path(stream.name)
            json.dump(value, stream, indent=2, allow_nan=False)
            stream.write("\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def input_pin(config, key):
    inputs = config.get("inputs")
    require(isinstance(inputs, dict), "config.inputs missing")
    entry = inputs.get(key)
    require(isinstance(entry, dict), f"config.inputs.{key} missing")
    source_path = entry.get("path")
    require(
        isinstance(source_path, str) and source_path.strip(),
        f"config.inputs.{key}.path missing",
    )
    sha256 = entry.get("sha256")
    require(
        isinstance(sha256, str)
        and len(sha256) == 64
        and all(char in "0123456789abcdef" for char in sha256),
        f"config.inputs.{key}.sha256 must be a 64-character lowercase hex string",
    )
    return {"configured_path": source_path, "sha256": sha256}


def exact_zero(expr):
    return s.cancel(s.together(expr)) == 0


def exact_matrix(matrix):
    return matrix.applyfunc(lambda value: s.cancel(s.together(value)))


def assert_zero_matrix(matrix, label):
    for row in range(matrix.rows):
        for col in range(matrix.cols):
            value = s.cancel(s.together(matrix[row, col]))
            require(exact_zero(value), f"{label} nonzero at ({row}, {col}): {value}")


def matrix_shape(matrix):
    return [matrix.rows, matrix.cols]


def domain_shape(matrix):
    return list(matrix.shape)


def as_strings(matrix):
    return [
        [str(s.cancel(s.together(matrix[row, col]))) for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]


def field_diff(expr):
    result = s.diff(expr, R) + T_PRIME * s.diff(expr, T)
    return s.cancel(s.together(result))


def field_diff_matrix(matrix):
    return matrix.applyfunc(field_diff)


def field_value(expr, name):
    require(
        not expr.has(s.Float, s.nan, s.oo, -s.oo, s.zoo),
        f"{name} is inexact or nonfinite",
    )
    return FIELD.from_sympy(expr)


def field_matrix(matrix, name):
    entries = {}
    for (row, col), expr in matrix.todok().items():
        value = field_value(expr, f"{name}[{row},{col}]")
        if value:
            entries.setdefault(row, {})[col] = value
    return DomainMatrix.from_dod(entries, matrix.shape, FIELD)


def field_zero_matrix(rows, cols):
    return DomainMatrix.from_dod({}, (rows, cols), FIELD)


def field_eye(size):
    return DomainMatrix.from_dod(
        {index: {index: FIELD_ONE} for index in range(size)},
        (size, size),
        FIELD,
    )


def field_selector(rows, column_count):
    return DomainMatrix.from_dod(
        {
            out_row: {source_row: FIELD_ONE}
            for out_row, source_row in enumerate(rows)
        },
        (len(rows), column_count),
        FIELD,
    )


def field_to_sympy(value):
    return FIELD.to_sympy(value)


def field_matrix_to_sympy(matrix):
    rows, cols = matrix.shape
    converted = s.zeros(rows, cols)
    for row, values in matrix.to_dod().items():
        for col, value in values.items():
            converted[row, col] = field_to_sympy(value)
    return converted


def assert_zero_field_matrix(matrix, label):
    require(matrix.is_zero_matrix, f"{label} is not zero")


def assert_field_identity(matrix, size, label):
    assert_zero_field_matrix(matrix.sub(field_eye(size)), label)


def field_derivative_matrix(matrix):
    field_r, field_t = FIELD.gens
    field_t_prime = -10 / (100 * field_r**2 + 1)
    entries = {}
    for row, values in matrix.to_dod().items():
        for col, value in values.items():
            derivative = value.diff(field_r) + field_t_prime * value.diff(field_t)
            if derivative:
                entries.setdefault(row, {})[col] = derivative
    return DomainMatrix.from_dod(entries, matrix.shape, FIELD)


def field_terminal_state(terminal):
    return {
        **terminal,
        "C": field_matrix(terminal["C"], "C"),
        "K": field_matrix(terminal["K"], "K"),
        "B": field_matrix(terminal["B"], "B"),
        "M": field_matrix(terminal["M"], "M"),
    }


def field_operator_matrices(matrices):
    return [
        field_matrix(matrix, f"operator_order_{order}")
        for order, matrix in enumerate(matrices)
    ]


def check_qq_expression(expr, name):
    expr = s.cancel(s.together(expr))
    require(not expr.has(s.atan), f"{name} contains atan")
    require(not expr.has(s.sqrt), f"{name} contains sqrt")
    require(not expr.has(EPSILON), f"{name} contains epsilon")
    require(not expr.has(s.Float, s.nan, s.oo, -s.oo, s.zoo), f"{name} is inexact or nonfinite")
    require(expr.free_symbols <= {R, T}, f"{name} contains unexpected symbols")
    numerator, denominator = s.fraction(expr)
    s.Poly(numerator, R, T, domain=s.QQ)
    denominator_poly = s.Poly(denominator, R, T, domain=s.QQ)
    require(not denominator_poly.is_zero, f"{name} has a zero denominator")
    return expr


def parse_raw_expression(text, name):
    require(isinstance(text, str), f"{name} is not a string expression")
    parsed = s.sympify(text, locals=SAFE_RAW_LOCALS)
    require(
        not parsed.has(s.Float, s.nan, s.oo, -s.oo, s.zoo),
        f"{name} is inexact or nonfinite",
    )
    fixed = s.cancel(s.together(parsed.subs(EPSILON, EPSILON0)))
    converted = s.cancel(s.together(fixed.xreplace({ATAN_FIXED: T})))
    recovered = converted.xreplace({T: ATAN_FIXED})
    require(
        exact_zero(recovered - fixed),
        f"{name} failed reversible Q(r,T) conversion",
    )
    return check_qq_expression(converted, name)


def parse_q_expression(text, name):
    require(isinstance(text, str), f"{name} is not a string expression")
    return check_qq_expression(s.sympify(text, locals=SAFE_Q_LOCALS), name)


def require_rectangular(data, rows, cols, name):
    require(isinstance(data, list), f"{name} is not a matrix")
    require(len(data) == rows, f"{name} row count changed")
    for row, values in enumerate(data):
        require(isinstance(values, list), f"{name} row {row} is not a list")
        require(len(values) == cols, f"{name} column count changed in row {row}")


def parse_raw_matrix(data, rows, cols, name):
    require_rectangular(data, rows, cols, name)
    matrix = s.zeros(rows, cols)
    for row in range(rows):
        for col in range(cols):
            matrix[row, col] = parse_raw_expression(
                data[row][col], f"{name}[{row},{col}]"
            )
    return exact_matrix(matrix)


def parse_q_matrix(data, rows, cols, name):
    require_rectangular(data, rows, cols, name)
    matrix = s.zeros(rows, cols)
    for row in range(rows):
        for col in range(cols):
            matrix[row, col] = parse_q_expression(data[row][col], f"{name}[{row},{col}]")
    return exact_matrix(matrix)


def parse_snapshot_matrix(payload, rows, cols, name):
    require(isinstance(payload, dict), f"{name} matrix payload missing")
    require(payload.get("shape") == [rows, cols], f"{name} shape changed")
    return parse_q_matrix(payload.get("entries"), rows, cols, f"{name}.entries")


def parse_operator_export(operators):
    require(
        operators.get("radial_DAE_reduced") is False,
        "operator export was unexpectedly pre-reduced",
    )
    require(
        operators.get("response_evaluated_on_admissible_non_gauge_seed") is False,
        "operator export unexpectedly pre-evaluated the discarded response",
    )
    require(operators.get("column_labels") == EXPECTED_COLUMNS, "operator column labels changed")
    require(operators.get("retained_row_labels") == EXPECTED_ROWS, "retained row labels changed")
    require(
        operators.get("selected_row_indices") == EXPECTED_SELECTED,
        "selected row ordering changed",
    )
    require(
        operators.get("automatic_zero_row_indices") == EXPECTED_AUTOMATIC,
        "automatic zero row set changed",
    )
    require(
        operators.get("operator_convention")
        == "sum_j M_j(r,epsilon) u^(j)(r); no derivative acts on M_j here",
        "operator radial derivative convention changed",
    )
    require(
        operators.get("last_row")
        == "discarded response D=6 E_X4-E_Cpolar4; not a retained equation",
        "discarded response convention changed",
    )
    matrices = operators.get("matrices_by_derivative_order")
    require(isinstance(matrices, list) and len(matrices) == 3, "expected 3 operator matrices")
    return [
        parse_raw_matrix(
            matrix,
            RESPONSE_ROW + 1,
            ORIGINAL_COLUMN_COUNT,
            f"operator_order_{order}",
        )
        for order, matrix in enumerate(matrices)
    ]


def parse_terminal_state(snapshot):
    require(
        snapshot.get("state") == "write-only intermediate terminal state before preservation",
        "terminal snapshot has unexpected state",
    )
    require(
        snapshot.get("state_coordinate_labels") == EXPECTED_STATE_LABELS,
        "terminal state coordinate labels changed",
    )
    require(snapshot.get("remaining_control_labels") == [], "remaining controls are not zero")
    matrices = snapshot.get("matrices")
    require(isinstance(matrices, dict), "terminal snapshot matrices missing")
    free_rows = snapshot.get("free_rows")
    require(isinstance(free_rows, list), "free_rows missing")
    require(len(free_rows) == FREE_STATE_COUNT, "free row count changed")
    require(
        sorted(free_rows) == sorted(set(free_rows))
        and all(type(row) is int and 0 <= row < STATE_COUNT for row in free_rows),
        "free_rows must be unique integer state row indices",
    )
    parsed = {
        "C": parse_snapshot_matrix(matrices.get("C"), 6, STATE_COUNT, "C"),
        "K": parse_snapshot_matrix(matrices.get("K"), STATE_COUNT, STATE_COUNT, "K"),
        "B": parse_snapshot_matrix(matrices.get("B"), STATE_COUNT, 0, "B"),
        "M": parse_snapshot_matrix(matrices.get("M"), STATE_COUNT, FREE_STATE_COUNT, "M"),
        "free_rows": list(free_rows),
        "ignored_publication_flags": {
            "original_constraint_product_zero": snapshot.get("original_constraint_product_zero"),
            "terminal_preservation_checked": snapshot.get("terminal_preservation_checked"),
            "original_22_retained_rows_checked": snapshot.get("original_22_retained_rows_checked"),
            "reduced_dynamics_reconstructed": snapshot.get("reduced_dynamics_reconstructed"),
            "discarded_response_checked": snapshot.get("discarded_response_checked"),
            "dae_closure_passed": snapshot.get("dae_closure_passed"),
            "smooth_solution_germ_certified": snapshot.get("smooth_solution_germ_certified"),
            "nonlinear_consistency_established": snapshot.get("nonlinear_consistency_established"),
        },
    }
    return parsed


def verify_identity_linkage(snapshot, terminal_state_sha256, operators_sha256):
    identity = snapshot.get("identity")
    require(isinstance(identity, dict), "terminal snapshot identity missing")
    require(
        snapshot.get("identity_sha256") == canonical_sha256(identity),
        "terminal snapshot identity SHA-256 mismatch",
    )
    input_hashes = identity.get("input_sha256")
    require(isinstance(input_hashes, dict), "terminal snapshot input hashes missing")
    require(
        input_hashes.get("stationary-operators.json") == operators_sha256,
        "terminal snapshot was not produced from this stationary-operators.json",
    )
    require(identity.get("epsilon0") == "1/100", "terminal snapshot epsilon changed")
    require(identity.get("T") == "atan(1/(10*r))", "terminal snapshot T convention changed")
    require(identity.get("T_prime") == str(T_PRIME), "terminal snapshot derivative changed")
    return {
        "terminal_state_sha256": terminal_state_sha256,
        "identity_sha256": snapshot.get("identity_sha256"),
        "source_commit": identity.get("source_commit"),
        "upstream_source_commit": identity.get("upstream_source_commit"),
        "upstream_identity_sha256": identity.get("upstream_identity_sha256"),
    }


def gauge_blocks(matrices):
    return [
        matrix.extract(range(RESPONSE_ROW + 1), KEPT_ORIGINAL_COLUMNS)
        for matrix in matrices
    ]


def row_selector(rows, column_count):
    selector = s.zeros(len(rows), column_count)
    for out_row, source_row in enumerate(rows):
        selector[out_row, source_row] = 1
    return selector


def complement(indices, size):
    present = set(indices)
    return [index for index in range(size) if index not in present]


def state_projectors():
    pv = s.eye(PROFILE_COUNT).row_join(s.zeros(PROFILE_COUNT, PROFILE_COUNT))
    pw = s.zeros(PROFILE_COUNT, PROFILE_COUNT).row_join(s.eye(PROFILE_COUNT))
    return pv, pw


def build_kernel_identities(terminal):
    c = terminal["C"]
    m = terminal["M"]
    return {
        "J": field_selector(terminal["free_rows"], STATE_COUNT),
        "CM": c.matmul(m),
        "JM": m.extract(terminal["free_rows"], range(FREE_STATE_COUNT)),
    }


def build_preservation_identities(terminal):
    c = terminal["C"]
    k = terminal["K"]
    m = terminal["M"]
    dc = field_derivative_matrix(c)
    s_matrix = dc.add(c.matmul(k))
    return {
        "DC": dc,
        "S": s_matrix,
        "SM": s_matrix.matmul(m),
    }


def build_physical_lift_identities(terminal):
    k = terminal["K"]
    m = terminal["M"]
    pw = field_selector(range(PROFILE_COUNT, STATE_COUNT), STATE_COUNT)
    lift_operator = k.extract(range(PROFILE_COUNT), range(STATE_COUNT)).sub(pw)
    return {
        "Pw": pw,
        "profile_lift": lift_operator.matmul(m),
    }


def state_operator(gauge_operators, k, rows):
    e0, e1, e2 = gauge_operators
    k_lower = k.extract(range(PROFILE_COUNT, STATE_COUNT), range(STATE_COUNT))
    return (
        e0.extract(rows, range(PROFILE_COUNT)).hstack(
            e1.extract(rows, range(PROFILE_COUNT))
        )
        .add(e2.extract(rows, range(PROFILE_COUNT)).matmul(k_lower))
    )


def build_retained_euler_identities(gauge_operators, terminal):
    retained_operator = state_operator(
        gauge_operators,
        terminal["K"],
        range(RETAINED_ROW_COUNT),
    )
    return {
        "retained_operator": retained_operator,
        "retained_on_kernel": retained_operator.matmul(terminal["M"]),
    }


def build_discarded_response_identities(gauge_operators, terminal):
    response_operator = state_operator(
        gauge_operators,
        terminal["K"],
        [RESPONSE_ROW],
    )
    return {
        "response_operator": response_operator,
        "response_on_kernel": response_operator.matmul(terminal["M"]),
    }


def build_state_identities(gauge_operators, terminal):
    identities = {}
    identities.update(build_kernel_identities(terminal))
    identities.update(build_preservation_identities(terminal))
    identities.update(build_physical_lift_identities(terminal))
    identities.update(build_retained_euler_identities(gauge_operators, terminal))
    identities.update(build_discarded_response_identities(gauge_operators, terminal))
    identities["proof_matrix_flags"] = {
        "CM_zero": True,
        "JM_identity": True,
        "SM_zero": True,
        "profile_lift_zero": True,
        "all_22_retained_rows_zero": True,
    }
    return {
        **identities,
    }


def assert_state_identities(identities):
    assert_zero_field_matrix(identities["CM"], "C M")
    assert_field_identity(identities["JM"], FREE_STATE_COUNT, "J M")
    assert_zero_field_matrix(identities["SM"], "(D C + C K) M")
    assert_zero_field_matrix(identities["profile_lift"], "(Pv K - Pw) M")
    assert_zero_field_matrix(
        identities["retained_on_kernel"],
        "retained Euler rows on M",
    )


def denominator_poly(expr):
    _, denominator = s.fraction(s.cancel(s.together(expr)))
    poly = s.Poly(denominator, R, T, domain=s.QQ)
    require(not poly.is_zero, "zero denominator polynomial")
    return poly.monic()


def field_denominator_poly(value):
    _, denominator = s.fraction(field_to_sympy(value))
    poly = s.Poly(denominator, R, T, domain=s.QQ)
    require(not poly.is_zero, "zero denominator polynomial")
    return poly.monic()


def collect_field_denominators(named_matrices):
    seen = {}
    for group, matrix in named_matrices:
        for row, values in matrix.to_dod().items():
            for col, value in values.items():
                poly = field_denominator_poly(value)
                expr = s.cancel(poly.as_expr())
                if expr == 1:
                    continue
                key = str(expr)
                seen.setdefault(
                    key,
                    {
                        "poly_r_T": key,
                        "locations": [],
                    },
                )
                seen[key]["locations"].append([group, row, col])
    derivative_poly = s.Poly(100 * R**2 + 1, R, T, domain=s.QQ).monic()
    seen.setdefault(
        str(derivative_poly.as_expr()),
        {
            "poly_r_T": str(derivative_poly.as_expr()),
            "locations": [["field_derivative", 0, 0]],
        },
    )
    return list(seen.values())


def interval_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def interval_sub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def interval_mul(a, b):
    values = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(values), max(values))


def interval_div(a, b):
    require(interval_sign(b) != 0, "interval denominator contains zero")
    reciprocals = [s.Rational(1, 1) / b[0], s.Rational(1, 1) / b[1]]
    return interval_mul(a, (min(reciprocals), max(reciprocals)))


def interval_horner(poly, interval):
    poly = s.Poly(poly, T, domain=s.QQ)
    result = (s.Rational(0), s.Rational(0))
    for coeff in poly.all_coeffs():
        result = interval_add(interval_mul(result, interval), (coeff, coeff))
    return result


def interval_sign(interval):
    if interval[0] > 0:
        return 1
    if interval[1] < 0:
        return -1
    return 0


def atan_bracket(argument, term_count):
    require(term_count > 0, "atan bracket needs at least one term")
    require(argument > 0 and argument <= 1, "atan bracket expects 0 < x <= 1")
    partial = s.Rational(0)
    for index in range(term_count):
        partial += (-1) ** index * argument ** (2 * index + 1) / (2 * index + 1)
    next_partial = partial + (-1) ** term_count * argument ** (2 * term_count + 1) / (
        2 * term_count + 1
    )
    return (min(partial, next_partial), max(partial, next_partial))


def atan_controls(argument, bracket):
    control_t = interval_horner(s.Poly(T, T, domain=s.QQ), bracket)
    control_t_minus_argument = interval_horner(
        s.Poly(T - argument, T, domain=s.QQ), bracket
    )
    require(control_t[0] > 0, "atan bracket did not prove T > 0")
    require(
        control_t_minus_argument[1] < 0,
        "atan bracket did not prove atan(x) < x",
    )
    return {
        "T_positive": [str(control_t[0]), str(control_t[1])],
        "T_minus_argument_negative": [
            str(control_t_minus_argument[0]),
            str(control_t_minus_argument[1]),
        ],
    }


def poly_at_r0(poly_or_expr, r0):
    return s.Poly(
        s.Poly(poly_or_expr, R, T, domain=s.QQ).as_expr().subs(R, r0),
        T,
        domain=s.QQ,
    )


def eval_poly_interval(poly_or_expr, r0, bracket):
    return interval_horner(poly_at_r0(poly_or_expr, r0), bracket)


def eval_expr_interval(expr, r0, bracket):
    numerator, denominator = s.fraction(s.cancel(s.together(expr)))
    numerator_interval = eval_poly_interval(numerator, r0, bracket)
    denominator_interval = eval_poly_interval(denominator, r0, bracket)
    return interval_div(numerator_interval, denominator_interval)


def certify_denominators_at_r0(denominators, r0, bracket):
    certificates = []
    for item in denominators:
        poly = s.sympify(item["poly_r_T"], locals=SAFE_Q_LOCALS)
        value_interval = eval_poly_interval(poly, r0, bracket)
        sign = interval_sign(value_interval)
        require(sign != 0, f"denominator not separated: {item['poly_r_T']}")
        certificates.append({
            **item,
            "r0": str(r0),
            "poly_at_r0": str(poly_at_r0(poly, r0).as_expr()),
            "interval": [str(value_interval[0]), str(value_interval[1])],
            "sign": sign,
        })
    return certificates


def interval_gaussian_det(matrix, r0, bracket):
    require(matrix.rows == matrix.cols, "minor is not square")
    entries = [
        [eval_expr_interval(matrix[row, col], r0, bracket) for col in range(matrix.cols)]
        for row in range(matrix.rows)
    ]
    determinant = (s.Rational(1), s.Rational(1))
    row_permutation = list(range(matrix.rows))
    pivots = []
    row_swaps = 0
    for col in range(matrix.cols):
        pivot = None
        for row in range(col, matrix.rows):
            if interval_sign(entries[row][col]) != 0:
                pivot = row
                break
        require(pivot is not None, f"minor pivot {col} was interval-singular")
        if pivot != col:
            entries[col], entries[pivot] = entries[pivot], entries[col]
            row_permutation[col], row_permutation[pivot] = (
                row_permutation[pivot],
                row_permutation[col],
            )
            determinant = (-determinant[1], -determinant[0])
            row_swaps += 1
        pivot_interval = entries[col][col]
        pivots.append({
            "column": col,
            "source_row": row_permutation[col],
            "interval": [str(pivot_interval[0]), str(pivot_interval[1])],
            "sign": interval_sign(pivot_interval),
        })
        determinant = interval_mul(determinant, pivot_interval)
        for row in range(col + 1, matrix.rows):
            factor = interval_div(entries[row][col], pivot_interval)
            entries[row][col] = (s.Rational(0), s.Rational(0))
            for target_col in range(col + 1, matrix.cols):
                entries[row][target_col] = interval_sub(
                    entries[row][target_col],
                    interval_mul(factor, entries[col][target_col]),
                )
    sign = interval_sign(determinant)
    require(sign != 0, "minor determinant interval contains zero")
    return {
        "determinant_interval": [str(determinant[0]), str(determinant[1])],
        "determinant_sign": sign,
        "pivot_intervals": pivots,
        "row_permutation": row_permutation,
        "row_swap_count": row_swaps,
    }


def numerator_poly(expr):
    numerator, _ = s.fraction(s.cancel(s.together(expr)))
    return s.Poly(numerator, R, T, domain=s.QQ)


def make_r0_attempt(r0, term_count):
    r0 = s.Rational(r0)
    require(r0 > 0, "r0 must be positive")
    argument = s.Rational(1, 10) / r0
    bracket = atan_bracket(argument, term_count)
    return r0, argument, bracket, atan_controls(argument, bracket)


def certify_at_r0(denominators, minor, r0, term_count):
    r0, argument, bracket, controls = make_r0_attempt(r0, term_count)
    denominator_certificates = certify_denominators_at_r0(denominators, r0, bracket)
    minor_certificate = interval_gaussian_det(minor, r0, bracket)
    return {
        "r0": str(r0),
        "terms": term_count,
        "atan_argument": str(argument),
        "atan_interval": [str(bracket[0]), str(bracket[1])],
        "controls": controls,
        "regularity_certificate": {
            "r0": str(r0),
            "terms": term_count,
            "atan_interval": [str(bracket[0]), str(bracket[1])],
            "denominator_count": len(denominator_certificates),
            "certified_denominators": denominator_certificates,
        },
        "minor_certificate": minor_certificate,
        "bracket": bracket,
    }


def response_candidates(response):
    candidates = []
    response_entries = response.to_dod().get(0, {})
    for col, value in sorted(response_entries.items()):
        candidates.append({
            "column": col,
            "poly": numerator_poly(field_to_sympy(value)),
        })
    return candidates


def certify_regularity_and_response(
    denominators,
    c_matrix,
    free_rows,
    response,
    max_atan_terms,
    r0_candidates=DEFAULT_R0S,
):
    minor_columns = complement(free_rows, STATE_COUNT)
    require(len(minor_columns) == 6, "natural pivot complement is not rank six")
    minor = field_matrix_to_sympy(c_matrix.extract(range(6), minor_columns))
    candidates = response_candidates(response)
    first_regular = None
    failures = []
    for r0 in r0_candidates:
        r0 = s.Rational(r0)
        for term_count in range(1, max_atan_terms + 1):
            try:
                regular = certify_at_r0(denominators, minor, r0, term_count)
            except ValueError as error:
                if term_count == max_atan_terms:
                    failures.append({
                        "r0": str(r0),
                        "stage": "regularity_or_rank_minor",
                        "error": str(error),
                    })
                continue
            if first_regular is None:
                first_regular = regular
            if not candidates:
                return first_regular, {
                    "status": "silent_on_forward_state_germ_subfamily",
                    "formal_response_nonzero": False,
                    "response_nonzero": False,
                    "certified_response_nonzero": False,
                    "stationary_universal_zero_claim": False,
                    "terminal_subfamily_response_zero": True,
                    "converse_certificate_supplied": False,
                }
            for candidate in candidates:
                interval = eval_poly_interval(candidate["poly"], r0, regular["bracket"])
                sign = interval_sign(interval)
                if sign == 0:
                    continue
                initial = ["0"] * FREE_STATE_COUNT
                initial[candidate["column"]] = "1"
                numerator_certificate = {
                    "column": candidate["column"],
                    "r0": regular["r0"],
                    "poly_r_T": str(candidate["poly"].as_expr()),
                    "poly_at_r0": str(poly_at_r0(candidate["poly"], r0).as_expr()),
                    "interval": [str(interval[0]), str(interval[1])],
                    "sign": sign,
                }
                return regular, {
                    "status": "response_nonzero_at_certified_r0",
                    "formal_response_nonzero": True,
                    "response_nonzero": True,
                    "certified_response_nonzero": True,
                    "stationary_universal_zero_claim": False,
                    "responding_initial_vector_supplied": True,
                    "converse_certificate_supplied": False,
                    "r0": regular["r0"],
                    "atan_interval": regular["atan_interval"],
                    "certified_column": candidate["column"],
                    "unit_vector_index": candidate["column"],
                    "unit_vector_kind": "free_state",
                    "certified_numerator": numerator_certificate,
                    "free_state_initial": initial,
                    "regularity_certificate": regular["regularity_certificate"],
                }
        if candidates:
            failures.append({
                "r0": str(r0),
                "stage": "response_numerator",
                "candidate_count": len(candidates),
                "error": "no response numerator separated at this radius",
            })
    require(first_regular is not None, "no regular rank-six radius certified")
    return first_regular, {
        "status": "formal_response_nonzero_unresolved",
        "formal_response_nonzero": True,
        "response_nonzero": None,
        "certified_response_nonzero": False,
        "stationary_universal_zero_claim": False,
        "responding_initial_vector_supplied": False,
        "converse_certificate_supplied": False,
        "candidate_count": len(candidates),
        "sampled_r0": [str(s.Rational(value)) for value in r0_candidates],
        "failures": failures,
    }


def blank_operator_system():
    return [s.zeros(RESPONSE_ROW + 1, ORIGINAL_COLUMN_COUNT) for _ in range(3)]


def manufactured_operator_export(matrices):
    return {
        "radial_DAE_reduced": False,
        "response_evaluated_on_admissible_non_gauge_seed": False,
        "column_labels": EXPECTED_COLUMNS,
        "retained_row_labels": EXPECTED_ROWS,
        "selected_row_indices": EXPECTED_SELECTED,
        "automatic_zero_row_indices": EXPECTED_AUTOMATIC,
        "operator_convention": (
            "sum_j M_j(r,epsilon) u^(j)(r); no derivative acts on M_j here"
        ),
        "last_row": "discarded response D=6 E_X4-E_Cpolar4; not a retained equation",
        "matrices_by_derivative_order": [as_strings(matrix) for matrix in matrices],
    }


def manufactured_terminal_state(c_factor=s.Rational(1)):
    c = s.zeros(6, STATE_COUNT)
    pivot_rows = [0, 1, 2, 7, 8, 9]
    for row, column in enumerate(pivot_rows):
        c[row, column] = c_factor if row == 0 else 1
    k = s.zeros(STATE_COUNT, STATE_COUNT)
    for index in range(PROFILE_COUNT):
        k[index, PROFILE_COUNT + index] = 1
    b = s.zeros(STATE_COUNT, 0)
    free_rows = complement(pivot_rows, STATE_COUNT)
    m = s.zeros(STATE_COUNT, len(free_rows))
    for col, row in enumerate(free_rows):
        m[row, col] = 1
    return {
        "C": c,
        "K": k,
        "B": b,
        "M": m,
        "free_rows": free_rows,
        "ignored_publication_flags": {},
    }


def manufactured_snapshot(state):
    return {
        "state": "write-only intermediate terminal state before preservation",
        "state_coordinate_labels": EXPECTED_STATE_LABELS,
        "remaining_control_labels": [],
        "free_rows": state["free_rows"],
        "matrices": {
            "C": {"shape": [6, STATE_COUNT], "entries": as_strings(state["C"])},
            "K": {"shape": [STATE_COUNT, STATE_COUNT], "entries": as_strings(state["K"])},
            "B": {"shape": [STATE_COUNT, 0], "entries": as_strings(state["B"])},
            "M": {
                "shape": [STATE_COUNT, FREE_STATE_COUNT],
                "entries": as_strings(state["M"]),
            },
        },
    }


def singleton_field_row(col, value):
    return DomainMatrix.from_dod(
        {0: {col: FIELD.from_sympy(value)}},
        (1, FREE_STATE_COUNT),
        FIELD,
    )


def require_rejects(callback, name):
    try:
        callback()
    except ValueError as error:
        return {"name": name, "rejected": True, "error": str(error)}
    raise AssertionError(f"{name} did not reject")


def verify_manufactured_controls():
    state = manufactured_terminal_state()
    field_state = field_terminal_state(state)
    operators = blank_operator_system()
    operators[0][RESPONSE_ROW, KEPT_ORIGINAL_COLUMNS[3]] = 1
    gauge_operators = gauge_blocks(field_operator_matrices(operators))
    identities = build_state_identities(gauge_operators, field_state)
    assert_state_identities(identities)
    denominators = collect_field_denominators([
        ("response", identities["response_on_kernel"])
    ])
    regular, response = certify_regularity_and_response(
        denominators,
        field_state["C"],
        field_state["free_rows"],
        identities["response_on_kernel"],
        max_atan_terms=4,
        r0_candidates=[s.Rational(1)],
    )
    require(response["certified_column"] == 0, "manufactured witness chose wrong column")
    top_level = certify(
        manufactured_operator_export(operators),
        manufactured_snapshot(state),
        max_atan_terms=4,
    )
    require(
        top_level["state"] == "state_germ_with_nonzero_discarded_response",
        "manufactured top-level dispatch lost the response witness",
    )
    require(top_level["B_shape"] == [STATE_COUNT, 0], "manufactured B shape changed")

    cm_bad = {**state, "M": state["M"].copy()}
    cm_bad["M"][0, 0] = 1
    jm_bad = {**state, "M": state["M"].copy()}
    jm_bad["M"][state["free_rows"][0], 0] = 2
    sm_bad = {**state, "K": state["K"].copy()}
    sm_bad["K"][0, state["free_rows"][0]] = 1
    lift_bad = {**state, "K": state["K"].copy()}
    lift_bad["K"][3, 3] = 1
    retained_bad = [matrix.copy() for matrix in operators]
    retained_bad[0][1, KEPT_ORIGINAL_COLUMNS[3]] = 1
    retained_bad_gauge = gauge_blocks(field_operator_matrices(retained_bad))

    minor_factor_state = manufactured_terminal_state(c_factor=R - 1)
    field_minor_factor_state = field_terminal_state(minor_factor_state)
    minor_identities = build_state_identities(
        gauge_operators,
        field_minor_factor_state,
    )
    minor_denominators = collect_field_denominators([
        ("response", minor_identities["response_on_kernel"])
    ])
    minor_regular, _ = certify_regularity_and_response(
        minor_denominators,
        field_minor_factor_state["C"],
        field_minor_factor_state["free_rows"],
        minor_identities["response_on_kernel"],
        max_atan_terms=4,
        r0_candidates=[s.Rational(1), s.Rational(2)],
    )
    require(minor_regular["r0"] == "2", "minor factor did not reject r0=1")

    row_swap_minor = s.Matrix([[0, 1], [R, 0]])
    row_swap = interval_gaussian_det(
        row_swap_minor, s.Rational(1), atan_bracket(s.Rational(1, 10), 1)
    )
    require(row_swap["row_swap_count"] == 1, "interval row swap did not occur")

    invalid_b_snapshot = {
        "state": "write-only intermediate terminal state before preservation",
        "state_coordinate_labels": EXPECTED_STATE_LABELS,
        "remaining_control_labels": [],
        "free_rows": state["free_rows"],
        "matrices": {
            "C": {"shape": [6, STATE_COUNT], "entries": as_strings(state["C"])},
            "K": {"shape": [STATE_COUNT, STATE_COUNT], "entries": as_strings(state["K"])},
            "B": {"shape": [STATE_COUNT, 1], "entries": [["0"] for _ in range(STATE_COUNT)]},
            "M": {
                "shape": [STATE_COUNT, FREE_STATE_COUNT],
                "entries": as_strings(state["M"]),
            },
        },
    }
    invalid_free_row_snapshot = manufactured_snapshot(state)
    invalid_free_row_snapshot["free_rows"] = [
        True,
        *state["free_rows"][1:],
    ]

    f = T / R
    expected_df = -T / R**2 - s.Rational(10) / (R * (100 * R**2 + 1))
    derivative_input = DomainMatrix.from_dod(
        {0: {0: FIELD.from_sympy(f)}},
        (1, 1),
        FIELD,
    )
    derivative = field_derivative_matrix(derivative_input)
    require(
        derivative.to_dod().get(0, {}).get(0, FIELD_ZERO)
        == FIELD.from_sympy(expected_df),
        "field derivative control failed",
    )
    zero_derivative = field_derivative_matrix(field_zero_matrix(2, 3))
    require(
        zero_derivative.shape == (2, 3) and zero_derivative.is_zero_matrix,
        "zero field derivative lost sparse shape",
    )

    zero_response = field_zero_matrix(1, FREE_STATE_COUNT)
    _, zero_response_certificate = certify_regularity_and_response(
        [],
        field_state["C"],
        field_state["free_rows"],
        zero_response,
        max_atan_terms=4,
        r0_candidates=[s.Rational(1)],
    )
    require(
        zero_response_certificate["status"]
        == "silent_on_forward_state_germ_subfamily",
        "zero response did not certify silent forward subfamily",
    )
    require(
        zero_response_certificate["stationary_universal_zero_claim"] is False,
        "zero response claimed universal silence",
    )

    unresolved_response = singleton_field_row(0, R - 1)
    _, unresolved_certificate = certify_regularity_and_response(
        [],
        field_state["C"],
        field_state["free_rows"],
        unresolved_response,
        max_atan_terms=4,
        r0_candidates=[s.Rational(1)],
    )
    require(
        unresolved_certificate["response_nonzero"] is None
        and unresolved_certificate["certified_response_nonzero"] is False,
        "unresolved response was not tri-state null",
    )

    singular_denominators = [{
        "poly_r_T": "r - 2",
        "locations": [["manufactured_denominator", 0, 0]],
    }]
    singular_then_regular, singular_then_response = certify_regularity_and_response(
        singular_denominators,
        field_state["C"],
        field_state["free_rows"],
        singleton_field_row(0, s.Integer(1)),
        max_atan_terms=4,
        r0_candidates=[s.Rational(2), s.Rational(1)],
    )
    require(
        singular_then_regular["r0"] == "1"
        and singular_then_response["r0"] == "1"
        and singular_then_response["certified_numerator"]["r0"] == "1",
        "singular first radius was accepted as a response witness",
    )

    cases = [
        {
            "name": "canonical_positive_response",
            "r0": regular["r0"],
            "certified_column": response["certified_column"],
        },
        {
            "name": "top_level_dispatch_positive_response",
            "state": top_level["state"],
            "unit_vector_index": top_level["response_certificate"][
                "unit_vector_index"
            ],
        },
        require_rejects(
            lambda: assert_state_identities(
                build_state_identities(gauge_operators, field_terminal_state(cm_bad))
            ),
            "CM_defect",
        ),
        require_rejects(
            lambda: assert_state_identities(
                build_state_identities(gauge_operators, field_terminal_state(jm_bad))
            ),
            "JM_defect",
        ),
        require_rejects(
            lambda: assert_state_identities(
                build_state_identities(gauge_operators, field_terminal_state(sm_bad))
            ),
            "SM_defect",
        ),
        require_rejects(
            lambda: assert_state_identities(
                build_state_identities(gauge_operators, field_terminal_state(lift_bad))
            ),
            "physical_lift_defect",
        ),
        require_rejects(
            lambda: assert_state_identities(
                build_state_identities(retained_bad_gauge, field_state)
            ),
            "retained_X0_row_defect",
        ),
        {
            "name": "minor_factor_rejects_r1_and_accepts_r2",
            "certified_r0": minor_regular["r0"],
        },
        {
            "name": "interval_row_swap_control",
            "row_swap_count": row_swap["row_swap_count"],
        },
        {
            "name": "field_derivative_T_over_r",
            "output_shape": list(derivative.shape),
        },
        {
            "name": "zero_field_derivative_sparse_shape",
            "output_shape": list(zero_derivative.shape),
        },
        {
            "name": "zero_response_silent_forward_subfamily",
            "status": zero_response_certificate["status"],
        },
        {
            "name": "formal_response_unresolved",
            "status": unresolved_certificate["status"],
            "response_nonzero": unresolved_certificate["response_nonzero"],
        },
        {
            "name": "response_skips_singular_radius",
            "certified_r0": singular_then_response["r0"],
        },
        require_rejects(
            lambda: interval_gaussian_det(
                s.Matrix([[0, 0], [0, 1]]),
                s.Rational(1),
                atan_bracket(s.Rational(1, 10), 1),
            ),
            "interval_singular_minor",
        ),
        require_rejects(
            lambda: parse_terminal_state(invalid_b_snapshot),
            "invalid_nonzero_control_shape",
        ),
        require_rejects(
            lambda: parse_terminal_state(invalid_free_row_snapshot),
            "invalid_bool_free_row",
        ),
        require_rejects(
            lambda: parse_q_expression("0.5", "manufactured_float"),
            "float_rejected",
        ),
    ]
    return {
        "state": "direct state-germ manufactured controls passed",
        "case_count": len(cases),
        "cases": cases,
    }


def certify(operators_json, snapshot_json, max_atan_terms, progress=None):
    progress = progress or NullProgress()
    with progress.stage(
        "original_operator_conversion",
        unit="convert",
        operator_orders=3,
        row_count=RESPONSE_ROW + 1,
        column_count=ORIGINAL_COLUMN_COUNT,
    ):
        raw_operators = parse_operator_export(operators_json)
        progress.metric(
            operator_order_count=len(raw_operators),
            retained_row_count=RETAINED_ROW_COUNT,
            original_column_count=ORIGINAL_COLUMN_COUNT,
        )
    with progress.stage("terminal_state_parsing", unit="parse"):
        terminal = parse_terminal_state(snapshot_json)
        progress.metric(
            constraint_rows=terminal["C"].rows,
            state_count=terminal["K"].rows,
            free_state_count=terminal["M"].cols,
            remaining_control_count=terminal["B"].cols,
        )
    with progress.stage("field_conversion", unit="convert"):
        raw_field_operators = field_operator_matrices(raw_operators)
        gauge_operators = gauge_blocks(raw_field_operators)
        field_terminal = field_terminal_state(terminal)
        progress.metric(
            field_operator_count=len(raw_field_operators),
            field_state_count=field_terminal["K"].shape[0],
            field_free_state_count=field_terminal["M"].shape[1],
        )
    with progress.stage("kernel_identities", unit="algebra"):
        identities = build_kernel_identities(field_terminal)
        assert_zero_field_matrix(identities["CM"], "C M")
        assert_field_identity(identities["JM"], FREE_STATE_COUNT, "J M")
        progress.metric(
            cm_zero=True,
            jm_identity=True,
            free_state_count=FREE_STATE_COUNT,
        )
    with progress.stage("preservation_identity", unit="algebra"):
        identities.update(build_preservation_identities(field_terminal))
        assert_zero_field_matrix(identities["SM"], "(D C + C K) M")
        progress.metric(sm_zero=True, preservation_rows=field_terminal["C"].shape[0])
    with progress.stage("physical_lift_identity", unit="algebra"):
        identities.update(build_physical_lift_identities(field_terminal))
        assert_zero_field_matrix(identities["profile_lift"], "(Pv K - Pw) M")
        progress.metric(profile_lift_zero=True, profile_count=PROFILE_COUNT)
    with progress.stage("retained_euler_identity", unit="algebra"):
        identities.update(build_retained_euler_identities(gauge_operators, field_terminal))
        assert_zero_field_matrix(identities["retained_on_kernel"], "retained Euler rows")
        progress.metric(all_22_retained_rows_zero=True, retained_row_count=22)
    with progress.stage("discarded_response_lift", unit="algebra"):
        identities.update(
            build_discarded_response_identities(gauge_operators, field_terminal)
        )
        response = identities["response_on_kernel"]
        formal_response_nonzero = not response.is_zero_matrix
        progress.metric(
            formal_discarded_response_nonzero=formal_response_nonzero,
            response_columns=response.shape[1],
        )
    with progress.stage("denominator_collection", unit="interval"):
        denominators = collect_field_denominators([
            ("C", field_terminal["C"]),
            ("K", field_terminal["K"]),
            ("M", field_terminal["M"]),
            ("raw_operators_M0", raw_field_operators[0]),
            ("raw_operators_M1", raw_field_operators[1]),
            ("raw_operators_M2", raw_field_operators[2]),
            ("DC", identities["DC"]),
            ("S", identities["S"]),
            ("profile_lift", identities["profile_lift"]),
            ("retained_operator", identities["retained_operator"]),
            ("retained_operator_on_M", identities["retained_on_kernel"]),
            ("response_operator", identities["response_operator"]),
            ("response_on_M", response),
        ])
        progress.metric(denominator_count=len(denominators))
    with progress.stage(
        "regularity_rank_and_response_search",
        unit="interval",
        denominator_count=len(denominators),
        max_atan_terms=max_atan_terms,
    ):
        regularity, response_certificate = certify_regularity_and_response(
            denominators,
            field_terminal["C"],
            field_terminal["free_rows"],
            response,
            max_atan_terms,
        )
        progress.metric(
            r0=regularity["r0"],
            terms=regularity["terms"],
            rank=6,
            denominator_count=len(denominators),
        )
    with progress.stage(
        "response_outcome_reporting",
        unit="interval",
        status=response_certificate["status"],
    ):
        progress.metric(
            status=response_certificate["status"],
            formal_response_nonzero=response_certificate["formal_response_nonzero"],
            certified_response_nonzero=response_certificate[
                "certified_response_nonzero"
            ],
        )

    rank_certificate = {
        "rank": 6,
        "minor_columns": complement(field_terminal["free_rows"], STATE_COUNT),
        "r0": regularity["r0"],
        "terms": regularity["terms"],
        "atan_interval": regularity["atan_interval"],
        **regularity["minor_certificate"],
    }
    response_nonzero = response_certificate["response_nonzero"]
    if response_certificate["certified_response_nonzero"]:
        state = "state_germ_with_nonzero_discarded_response"
    elif response_certificate["formal_response_nonzero"]:
        state = "state_germ_response_unresolved"
    else:
        state = "state_germ_forward_subfamily_silent"
    return {
        "state": state,
        "epsilon0": "1/100",
        "T": "atan(1/(10*r))",
        "T_prime": str(T_PRIME),
        "state_germ_certificate_passed": True,
        "smooth_solution_germ_certified": True,
        "free_state_count": FREE_STATE_COUNT,
        "remaining_control_count": 0,
        "B_shape": [STATE_COUNT, 0],
        "same_radius_witness_required": True,
        "formal_discarded_response_nonzero": response_certificate[
            "formal_response_nonzero"
        ],
        "discarded_response_nonzero": response_nonzero,
        "certified_discarded_response_nonzero": response_certificate[
            "certified_response_nonzero"
        ],
        "admissible_non_gauge_seed_supplied": response_certificate[
            "certified_response_nonzero"
        ],
        "stationary_universal_zero_claim": False,
        "physical_boundary_domain_certified": False,
        "final_dirac_classification_established": False,
        "nonlinear_consistency_established": False,
        "proof_matrix_flags": {
            "CM_zero": True,
            "JM_identity": True,
            "SM_zero": True,
            "profile_lift_zero": True,
            "all_22_retained_rows_zero": True,
        },
        "constraint_rank_certificate": rank_certificate,
        "regularity_certificate": regularity["regularity_certificate"],
        "retained_euler_certificate": {
            "retained_row_count": RETAINED_ROW_COUNT,
            "x0_equation_included_in_22_retained_rows": True,
            "operator_identity": "(E0 Pv + E1 Pw + E2 Pw K) M",
            "all_22_retained_rows_zero": True,
        },
        "response_certificate": response_certificate,
        "snapshot_publication_flags_ignored": field_terminal[
            "ignored_publication_flags"
        ],
    }


def run(config_path, operators_path, state_path, output_path, max_atan_terms):
    if not __debug__:
        raise RuntimeError("scientific assertions must remain enabled")
    config_path = Path(config_path)
    operators_path = Path(operators_path)
    state_path = Path(state_path)
    progress = StageProgress()
    progress.event(
        "invocation_begin",
        python=sys.version,
        executable=sys.executable,
        sympy=s.__version__,
        max_atan_terms=max_atan_terms,
        config_path=str(config_path),
        operators_path=str(operators_path),
        state_path=str(state_path),
    )
    with progress.stage("input_pin_verification", unit="input"):
        config, config_bytes, config_sha256 = load_json_bytes(config_path)
        operators_pin = input_pin(config, "stationary_operators")
        state_pin = input_pin(config, "terminal_state")
        operators, operator_bytes, operators_sha256 = load_pinned_json(
            operators_path,
            operators_pin["sha256"],
            "stationary operators",
        )
        terminal_state, state_bytes, state_sha256 = load_pinned_json(
            state_path,
            state_pin["sha256"],
            "terminal state",
        )
        progress.metric(
            config_bytes=len(config_bytes),
            operator_bytes=len(operator_bytes),
            terminal_state_bytes=len(state_bytes),
            config_sha256=config_sha256,
            stationary_operators_sha256=operators_sha256,
            terminal_state_sha256=state_sha256,
        )
    with progress.stage("identity_linkage", unit="check"):
        identity_linkage = verify_identity_linkage(
            terminal_state, state_sha256, operators_sha256
        )
        progress.metric(identity_link_verified=True)
    with progress.stage("manufactured_controls", unit="check"):
        manufactured_controls = verify_manufactured_controls()
        progress.metric(manufactured_case_count=manufactured_controls["case_count"])

    result = certify(
        operators,
        terminal_state,
        max_atan_terms=max_atan_terms,
        progress=progress,
    )
    input_sha256 = {
        "config": config_sha256,
        "stationary_operators": operators_sha256,
        "terminal_state": state_sha256,
    }
    with progress.stage("final_export", unit="bundle"):
        atomic_json_new(output_path, {
            **result,
            "runtime_environment": {
                "python": sys.version,
                "executable": sys.executable,
                "sympy": s.__version__,
            },
            "input_sha256": input_sha256,
            "configured_input_paths": {
                "stationary_operators": operators_pin["configured_path"],
                "terminal_state": state_pin["configured_path"],
            },
            "identity_linkage": identity_linkage,
            "manufactured_controls": manufactured_controls,
        })
        progress.metric(output_path=str(output_path))
    progress.report(
        phase="complete",
        step=1,
        of=1,
        unit="certificate",
        m={"state_germ_certificate_passed": 1},
    )
    progress.event("invocation_return")
    print("PASS: direct 14-state stationary germ certified.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--operators", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-atan-terms", type=int, default=256)
    args = parser.parse_args()
    run(
        args.config,
        args.operators,
        args.state,
        args.output,
        args.max_atan_terms,
    )


if __name__ == "__main__":
    main()

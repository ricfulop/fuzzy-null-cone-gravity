"""Independent terminal germ certificate for the oblate stationary reduction.

This script consumes only:

* the original completed stationary-operator export;
* a terminal dae-closure.json from the fixed-epsilon reducer; and
* a small config pinning the two input SHA-256 digests.

It does not import the extractor or the reducer.  The certificate substitutes the
terminal reduced system back into the 23 original rows and then proves local
regularity at r=1 by rational interval arithmetic for T=atan(1/10).
"""

from pathlib import Path
import argparse
import hashlib
import json
import os
import tempfile

import sympy as s


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
ORIGINAL_COLUMN_COUNT = 8
STATE_COUNT = 14

SAFE_RAW_LOCALS = {
    "r": R,
    "epsilon": EPSILON,
    "atan": s.atan,
    "sqrt": s.sqrt,
}
SAFE_Q_LOCALS = {"r": R, "T": T}
KEPT_ORIGINAL_COLUMNS = [
    index for index in range(ORIGINAL_COLUMN_COUNT) if index != X0_COLUMN
]
DEFAULT_RESPONSE_R0S = [
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


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256_bytes(data):
    return hashlib.sha256(data).hexdigest()


def load_json_bytes(path):
    data = Path(path).read_bytes()
    return json.loads(data.decode("utf-8")), data, sha256_bytes(data)


def load_pinned_json(path, expected, label):
    data = Path(path).read_bytes()
    actual = sha256_bytes(data)
    require(actual == expected, f"{label} SHA-256 mismatch")
    return json.loads(data.decode("utf-8")), data, actual


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


def expected_sha256(config, name):
    aliases = {
        "stationary_operators": ["stationary_operators", "stationary-operators.json"],
        "dae_closure": ["dae_closure", "dae-closure.json"],
    }[name]
    candidates = [
        config.get(f"{alias}_sha256")
        for alias in aliases
    ] + [
        config.get(section, {}).get(alias)
        for section in ["sha256", "expected_sha256"]
        for alias in aliases
    ] + [
        config.get("inputs", {}).get(alias, {}).get("sha256")
        for alias in aliases
    ]
    for candidate in candidates:
        if candidate is not None:
            require(
                isinstance(candidate, str) and len(candidate) == 64,
                f"{name} SHA-256 must be a 64-character hex string",
            )
            return candidate
    raise ValueError(f"missing pinned SHA-256 for {name}")


def exact_zero(expr):
    return s.cancel(s.together(expr)) == 0


def exact_matrix(matrix):
    return matrix.applyfunc(lambda value: s.cancel(s.together(value)))


def assert_zero_matrix(matrix, label):
    for i in range(matrix.rows):
        for j in range(matrix.cols):
            value = s.cancel(s.together(matrix[i, j]))
            require(exact_zero(value), f"{label} nonzero at ({i}, {j}): {value}")


def field_diff(expr, order=1):
    result = expr
    for _ in range(order):
        result = s.diff(result, R) + T_PRIME * s.diff(result, T)
        result = s.cancel(s.together(result))
    return result


def field_diff_matrix(matrix):
    return matrix.applyfunc(field_diff)


def as_strings(matrix):
    return [
        [str(s.cancel(s.together(matrix[i, j]))) for j in range(matrix.cols)]
        for i in range(matrix.rows)
    ]


def check_qq_expression(expr, name):
    expr = s.cancel(s.together(expr))
    require(not expr.has(s.atan), f"{name} contains atan")
    require(not expr.has(s.sqrt), f"{name} contains sqrt")
    require(not expr.has(EPSILON), f"{name} contains epsilon")
    require(not expr.atoms(s.Float), f"{name} contains inexact Float atoms")
    require(expr.free_symbols <= {R, T}, f"{name} contains unexpected symbols")
    numerator, denominator = s.fraction(expr)
    s.Poly(numerator, R, T, domain=s.QQ)
    denominator_poly = s.Poly(denominator, R, T, domain=s.QQ)
    require(not denominator_poly.is_zero, f"{name} has a zero denominator")
    return expr


def parse_raw_expression(text, name):
    require(isinstance(text, str), f"{name} is not a string expression")
    parsed = s.sympify(text, locals=SAFE_RAW_LOCALS)
    require(not parsed.atoms(s.Float), f"{name} contains inexact Float atoms")
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


def parse_q_matrix(data, rows, cols, name):
    require_rectangular(data, rows, cols, name)
    matrix = s.zeros(rows, cols)
    for row in range(rows):
        for col in range(cols):
            matrix[row, col] = parse_q_expression(data[row][col], f"{name}[{row},{col}]")
    return exact_matrix(matrix)


def parse_raw_matrix(data, rows, cols, name):
    require_rectangular(data, rows, cols, name)
    matrix = s.zeros(rows, cols)
    for row in range(rows):
        for col in range(cols):
            matrix[row, col] = parse_raw_expression(
                data[row][col], f"{name}[{row},{col}]"
            )
    return exact_matrix(matrix)


def parse_operator_export(operators):
    require(
        operators.get("radial_DAE_reduced") is False,
        "operator export was unexpectedly pre-reduced",
    )
    require(
        operators.get("response_evaluated_on_admissible_non_gauge_seed") is False,
        "operator export unexpectedly pre-evaluated the discarded response",
    )
    require(
        operators.get("column_labels") == EXPECTED_COLUMNS,
        "operator column labels changed",
    )
    require(
        operators.get("retained_row_labels") == EXPECTED_ROWS,
        "operator retained row labels changed",
    )
    require(
        operators.get("selected_row_indices") == EXPECTED_SELECTED,
        "operator selected row ordering changed",
    )
    require(
        operators.get("automatic_zero_row_indices") == EXPECTED_AUTOMATIC,
        "operator automatic row set changed",
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
    require(isinstance(matrices, list) and len(matrices) == 3, "expected 3 matrices")
    return [
        parse_raw_matrix(
            matrix,
            RESPONSE_ROW + 1,
            ORIGINAL_COLUMN_COUNT,
            f"operator_order_{order}",
        )
        for order, matrix in enumerate(matrices)
    ]


def infer_closure_widths(closure):
    embedding = closure.get("embedding")
    require(isinstance(embedding, list), "embedding is not a matrix")
    require(len(embedding) == STATE_COUNT, "embedding row count changed")
    free_count = None
    for row, values in enumerate(embedding):
        require(isinstance(values, list), f"embedding row {row} is not a list")
        if free_count is None:
            free_count = len(values)
        require(len(values) == free_count, "embedding is not rectangular")
    control_labels = closure.get("remaining_control_labels")
    require(isinstance(control_labels, list), "remaining controls changed")
    return free_count or 0, len(control_labels)


def parse_closure(closure):
    require(
        closure.get("constant_control_augmented_state_supported") is True,
        "closure does not preserve the constant-control convention",
    )
    require(
        closure.get("smooth_solution_germ_certified") is False,
        "closure unexpectedly claims an actual germ",
    )
    require(
        closure.get("nonlinear_consistency_established") is False,
        "closure unexpectedly claims nonlinear consistency",
    )
    free_count, control_count = infer_closure_widths(closure)
    embedding = parse_q_matrix(closure["embedding"], STATE_COUNT, free_count, "embedding")
    reduced_k = parse_q_matrix(
        closure.get("reduced_K"), free_count, free_count, "reduced_K"
    )
    reduced_b = parse_q_matrix(
        closure.get("reduced_B"), free_count, control_count, "reduced_B"
    )
    response_state = parse_q_matrix(
        closure.get("response_state"), 1, free_count, "response_state"
    )
    response_control = parse_q_matrix(
        closure.get("response_control"), 1, control_count, "response_control"
    )
    return {
        "free_count": free_count,
        "control_count": control_count,
        "embedding": embedding,
        "reduced_K": reduced_k,
        "reduced_B": reduced_b,
        "response": response_state.row_join(response_control),
    }


def augment_constant_controls(reduced_k, reduced_b):
    free_count, control_count = reduced_b.rows, reduced_b.cols
    return reduced_k.row_join(reduced_b).col_join(
        s.zeros(control_count, free_count + control_count)
    )


def gauge_blocks(matrices):
    return [
        matrix.extract(range(RESPONSE_ROW + 1), KEPT_ORIGINAL_COLUMNS)
        for matrix in matrices
    ]


def split_profile_blocks(embedding):
    columns = range(embedding.cols)
    return (
        embedding.extract(range(0, len(KEPT_ORIGINAL_COLUMNS)), columns),
        embedding.extract(
            range(len(KEPT_ORIGINAL_COLUMNS), 2 * len(KEPT_ORIGINAL_COLUMNS)),
            columns,
        ),
    )


def build_augmented_profiles(state_block, control_block):
    free_count = state_block.cols
    control_count = control_block.cols
    profiles = s.zeros(ORIGINAL_COLUMN_COUNT, free_count + control_count)
    for gauge_row, original_column in enumerate(KEPT_ORIGINAL_COLUMNS):
        for col in range(free_count):
            profiles[original_column, col] = state_block[gauge_row, col]
        for col in range(control_count):
            profiles[original_column, free_count + col] = control_block[gauge_row, col]
    return profiles


def direct_terminal_residual(matrices, closure):
    e0, e1, e2 = gauge_blocks(matrices)
    p_block, s_block = split_profile_blocks(closure["embedding"])
    reduced_k = closure["reduced_K"]
    reduced_b = closure["reduced_B"]
    acceleration_state = exact_matrix(field_diff_matrix(s_block) + s_block * reduced_k)
    acceleration_control = exact_matrix(s_block * reduced_b)
    residual_state = exact_matrix(
        e0 * p_block + e1 * s_block + e2 * acceleration_state
    )
    residual_control = exact_matrix(e2 * acceleration_control)
    return {
        "P": p_block,
        "S": s_block,
        "acceleration_state": acceleration_state,
        "acceleration_control": acceleration_control,
        "profile_derivative_state_defect": exact_matrix(
            field_diff_matrix(p_block) + p_block * reduced_k - s_block
        ),
        "profile_derivative_control_defect": exact_matrix(p_block * reduced_b),
        "residual_state": residual_state,
        "residual_control": residual_control,
        "response": residual_state.row(RESPONSE_ROW).row_join(
            residual_control.row(RESPONSE_ROW)
        ),
    }


def augmented_terminal_residual(matrices, direct, closure):
    zeros = s.zeros(len(KEPT_ORIGINAL_COLUMNS), closure["control_count"])
    u0 = build_augmented_profiles(direct["P"], zeros)
    expected_u1 = build_augmented_profiles(direct["S"], zeros)
    expected_u2 = build_augmented_profiles(
        direct["acceleration_state"], direct["acceleration_control"]
    )
    kaug = augment_constant_controls(closure["reduced_K"], closure["reduced_B"])
    u1 = exact_matrix(field_diff_matrix(u0) + u0 * kaug)
    u2 = exact_matrix(field_diff_matrix(u1) + u1 * kaug)
    residual = exact_matrix(matrices[0] * u0 + matrices[1] * u1 + matrices[2] * u2)
    return {
        "U": u0,
        "U1": u1,
        "U2": u2,
        "Kaug": kaug,
        "residual": residual,
        "expected_residual": direct["residual_state"].row_join(
            direct["residual_control"]
        ),
        "U1_defect": exact_matrix(u1 - expected_u1),
        "U2_defect": exact_matrix(u2 - expected_u2),
    }


def assert_terminal_residuals(matrices, closure):
    direct = direct_terminal_residual(matrices, closure)
    assert_zero_matrix(
        direct["profile_derivative_state_defect"],
        "profile derivative state defect",
    )
    assert_zero_matrix(
        direct["profile_derivative_control_defect"],
        "profile derivative control defect",
    )
    assert_zero_matrix(
        direct["residual_state"].extract(range(RESPONSE_ROW), range(
            direct["residual_state"].cols
        )),
        "retained terminal state residual",
    )
    assert_zero_matrix(
        direct["residual_control"].extract(range(RESPONSE_ROW), range(
            direct["residual_control"].cols
        )),
        "retained terminal control residual",
    )
    response_difference = exact_matrix(
        direct["response"] - closure["response"]
    )
    assert_zero_matrix(response_difference, "discarded response mismatch")
    augmented = augmented_terminal_residual(matrices, direct, closure)
    assert_zero_matrix(augmented["U1_defect"], "augmented U1 defect")
    assert_zero_matrix(augmented["U2_defect"], "augmented U2 defect")
    assert_zero_matrix(
        exact_matrix(augmented["residual"] - augmented["expected_residual"]),
        "augmented residual defect",
    )
    return {**direct, **augmented}


def denominator_poly(expr):
    _, denominator = s.fraction(s.cancel(s.together(expr)))
    poly = s.Poly(denominator, R, T, domain=s.QQ)
    if poly.is_zero:
        raise ValueError("zero denominator polynomial")
    return poly.monic()


def collect_denominators(named_matrices):
    seen = {}
    for group, matrix in named_matrices:
        for row in range(matrix.rows):
            for col in range(matrix.cols):
                poly = denominator_poly(matrix[row, col])
                expr = s.cancel(poly.as_expr())
                if expr == 1:
                    continue
                key = str(expr)
                seen.setdefault(
                    key,
                    {
                        "poly_r_T": key,
                        "poly_at_r1": str(s.Poly(expr.subs(R, 1), T).as_expr()),
                        "locations": [],
                    },
                )
                seen[key]["locations"].append([group, row, col])
    for expr in [100 * R**2 + 1]:
        poly = s.Poly(expr, R, T, domain=s.QQ).monic()
        key = str(poly.as_expr())
        seen.setdefault(
            key,
            {
                "poly_r_T": key,
                "poly_at_r1": str(s.Poly(poly.as_expr().subs(R, 1), T).as_expr()),
                "locations": [["field_derivative", 0, 0]],
            },
        )
    return list(seen.values())


def interval_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def interval_mul(a, b):
    products = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(products), max(products))


def interval_horner(poly, interval):
    poly = s.Poly(poly, T, domain=s.QQ)
    result = (s.Rational(0), s.Rational(0))
    for coeff in poly.all_coeffs():
        result = interval_add(interval_mul(result, interval), (coeff, coeff))
    return result


def atan_bracket(argument, term_count):
    require(term_count > 0, "atan bracket needs at least one term")
    require(argument > 0 and argument <= 1, "atan bracket expects 0 < x <= 1")
    partial = s.Rational(0)
    for n in range(term_count):
        partial += (-1) ** n * argument ** (2 * n + 1) / (2 * n + 1)
    next_partial = partial + (-1) ** term_count * argument ** (2 * term_count + 1) / (
        2 * term_count + 1
    )
    return (min(partial, next_partial), max(partial, next_partial))


def interval_sign(interval):
    if interval[0] > 0:
        return 1
    if interval[1] < 0:
        return -1
    return 0


def certify_polys_at_r0(polys, r0, max_terms):
    r0 = s.Rational(r0)
    argument = s.Rational(1, 10) / r0
    last_certificates = None
    for term_count in range(1, max_terms + 1):
        bracket = atan_bracket(argument, term_count)
        control_t = interval_horner(s.Poly(T, T, domain=s.QQ), bracket)
        control_t_minus_argument = interval_horner(
            s.Poly(T - argument, T, domain=s.QQ), bracket
        )
        if not (control_t[0] > 0 and control_t_minus_argument[1] < 0):
            continue
        certificates = []
        failed = False
        for item in polys:
            poly_r_t = s.Poly(
                s.sympify(item["poly_r_T"], locals=SAFE_Q_LOCALS),
                R,
                T,
                domain=s.QQ,
            )
            poly = s.Poly(poly_r_t.as_expr().subs(R, r0), T, domain=s.QQ)
            value_interval = interval_horner(poly, bracket)
            sign = interval_sign(value_interval)
            if sign == 0:
                failed = True
                break
            certificates.append({
                **item,
                "r0": str(r0),
                "atan_argument": str(argument),
                "poly_at_r0": str(poly.as_expr()),
                "interval": [str(value_interval[0]), str(value_interval[1])],
                "sign": sign,
            })
        last_certificates = certificates
        if not failed:
            return {
                "r0": str(r0),
                "terms": term_count,
                "atan_argument": str(argument),
                "atan_interval": [str(bracket[0]), str(bracket[1])],
                "controls": {
                    "T_positive": [str(control_t[0]), str(control_t[1])],
                    "T_minus_argument_negative": [
                        str(control_t_minus_argument[0]),
                        str(control_t_minus_argument[1]),
                    ],
                },
                "certificates": certificates,
            }
    raise ValueError(
        f"could not separate all polynomials at atan(1/(10*{r0})); "
        f"last certificate count {len(last_certificates or [])}"
    )


def response_polynomial(poly_expr):
    return {
        "poly_r_T": str(s.Poly(poly_expr, R, T, domain=s.QQ).as_expr()),
        "poly_at_r1": str(s.Poly(poly_expr.subs(R, 1), T, domain=s.QQ).as_expr()),
    }


def denominator_certificate_summary(denominators, certificate):
    return {
        "r0": certificate["r0"],
        "terms": certificate["terms"],
        "atan_interval": certificate["atan_interval"],
        "controls": certificate["controls"],
        "coefficient_polynomial_count": len(denominators),
        "coefficients": certificate["certificates"],
    }


def certify_response(
    closure, denominators, max_terms, r0_candidates=DEFAULT_RESPONSE_R0S
):
    response = closure["response"]
    if all(exact_zero(value) for value in response):
        return {
            "status": "silent_on_terminal_subfamily",
            "formal_response_nonzero": False,
            "response_nonzero": False,
            "certified_response_nonzero": False,
            "admissible_non_gauge_seed_supplied": False,
            "converse_certificate_supplied": False,
            "terminal_subfamily_response_zero": True,
            "stationary_universal_zero_claim": False,
            "universal_zero_claim": False,
        }

    candidates = []
    for column in range(response.cols):
        value = s.cancel(s.together(response[0, column]))
        if exact_zero(value):
            continue
        numerator, _ = s.fraction(value)
        candidates.append({
            "column": column,
            "kind": "free_state" if column < closure["free_count"] else "constant_control",
            "poly": response_polynomial(numerator),
        })

    failures = []
    for r0 in r0_candidates:
        try:
            regularity_certificate = certify_polys_at_r0(
                denominators, r0, max_terms
            )
        except ValueError as error:
            failures.append({
                "r0": str(s.Rational(r0)),
                "stage": "terminal_denominators",
                "error": str(error),
            })
            continue

        for candidate in candidates:
            try:
                certificate = certify_polys_at_r0(
                    [candidate["poly"]], r0, max_terms
                )
            except ValueError as error:
                failures.append({
                    "r0": str(s.Rational(r0)),
                    "stage": "response_numerator",
                    "column": candidate["column"],
                    "error": str(error),
                })
                continue
            free_initial = [0] * closure["free_count"]
            control_initial = [0] * closure["control_count"]
            column = candidate["column"]
            if column < closure["free_count"]:
                free_initial[column] = 1
            else:
                control_initial[column - closure["free_count"]] = 1
            return {
                "status": "response_nonzero_at_certified_r0",
                "formal_response_nonzero": True,
                "response_nonzero": True,
                "certified_response_nonzero": True,
                "admissible_non_gauge_seed_supplied": True,
                "responding_initial_vector_supplied": True,
                "converse_certificate_supplied": False,
                "stationary_universal_zero_claim": False,
                "universal_zero_claim": False,
                "certified_column": candidate["column"],
                "certified_column_kind": candidate["kind"],
                "certified_numerator": {
                    **candidate["poly"],
                    **certificate["certificates"][0],
                },
                "regularity_certificate": denominator_certificate_summary(
                    denominators, regularity_certificate
                ),
                "atan_certificate": {
                    key: certificate[key]
                    for key in ["r0", "terms", "atan_interval", "controls"]
                },
                "free_state_initial": [str(value) for value in free_initial],
                "constant_control_initial": [str(value) for value in control_initial],
            }

    return {
        "status": "formal_response_nonzero_but_unresolved_at_sampled_r0",
        "formal_response_nonzero": True,
        "response_nonzero": None,
        "certified_response_nonzero": False,
        "admissible_non_gauge_seed_supplied": False,
        "responding_initial_vector_supplied": False,
        "converse_certificate_supplied": False,
        "stationary_universal_zero_claim": False,
        "universal_zero_claim": False,
        "sampled_r0": [str(value) for value in r0_candidates],
        "candidate_count": len(candidates),
        "failures": failures,
    }


def blank_operator_system():
    return [s.zeros(RESPONSE_ROW + 1, ORIGINAL_COLUMN_COUNT) for _ in range(3)]


def manufactured_closure():
    embedding = s.zeros(STATE_COUNT, 1)
    embedding[0, 0] = 1
    return {
        "free_count": 1,
        "control_count": 0,
        "embedding": embedding,
        "reduced_K": s.zeros(1, 1),
        "reduced_B": s.zeros(1, 0),
        "response": s.zeros(1, 1),
    }


def require_rejects(callback, name):
    try:
        callback()
    except ValueError as error:
        return {"name": name, "rejected": True, "error": str(error)}
    raise AssertionError(f"{name} did not reject")


def verify_response_radius_requires_regular_denominators():
    closure = {
        "free_count": 1,
        "control_count": 0,
        "response": s.Matrix([[1]]),
    }
    singular_at_two = [{
        "poly_r_T": str(R - 2),
        "locations": [["manufactured", 0, 0]],
    }]

    unresolved = certify_response(
        closure, singular_at_two, max_terms=4, r0_candidates=[s.Rational(2)]
    )
    require(
        unresolved["status"] == "formal_response_nonzero_but_unresolved_at_sampled_r0",
        "singular response radius was not rejected",
    )

    certified = certify_response(
        closure,
        singular_at_two,
        max_terms=4,
        r0_candidates=[s.Rational(2), s.Rational(1)],
    )
    require(
        certified["status"] == "response_nonzero_at_certified_r0",
        "regular fallback response radius was not certified",
    )
    require(
        certified["atan_certificate"]["r0"] == "1",
        "response certified at a singular radius",
    )
    require(
        certified["regularity_certificate"]["r0"] == "1",
        "response regularity certificate used a different radius",
    )
    return {
        "name": "response_radius_requires_regular_denominators",
        "singular_r0_status": unresolved["status"],
        "certified_r0": certified["atan_certificate"]["r0"],
    }


def verify_control_only_response():
    closure = {
        "free_count": 2,
        "control_count": 1,
        "embedding": s.zeros(STATE_COUNT, 2),
        "reduced_K": s.Matrix([[0, 1], [0, 0]]),
        "reduced_B": s.Matrix([[0], [1]]),
        "response": s.Matrix([[0, 0, 1]]),
    }
    closure["embedding"][0, 0] = 1
    closure["embedding"][len(KEPT_ORIGINAL_COLUMNS), 1] = 1

    matrices = blank_operator_system()
    matrices[2][RESPONSE_ROW, 0] = 1
    terminal = assert_terminal_residuals(matrices, closure)
    response = certify_response(
        closure, denominators=[], max_terms=4, r0_candidates=[s.Rational(1)]
    )

    require(
        terminal["acceleration_control"][0, 0] == 1,
        "manufactured control did not survive into v''",
    )
    require(
        response["certified_column_kind"] == "constant_control",
        "control-only response did not select a control column",
    )
    require(
        response["constant_control_initial"] == ["1"],
        "control-only response did not choose a unit control datum",
    )
    return {
        "name": "control_only_response",
        "certified_column_kind": response["certified_column_kind"],
        "constant_control_initial": response["constant_control_initial"],
    }


def verify_manufactured_controls():
    closure = manufactured_closure()
    matrices = blank_operator_system()
    assert_terminal_residuals(matrices, closure)

    bad_retained = [matrix.copy() for matrix in matrices]
    bad_retained[0][0, 0] = 1

    bad_response = [matrix.copy() for matrix in matrices]
    bad_response[0][RESPONSE_ROW, 0] = 1

    return {
        "state": "terminal germ manufactured controls passed",
        "cases": [
            {
                "name": "zero_system",
                "retained_rows_zero": True,
                "response_matches": True,
            },
            require_rejects(
                lambda: assert_terminal_residuals(bad_retained, closure),
                "wrong_nonzero_retained_residual",
            ),
            require_rejects(
                lambda: assert_terminal_residuals(bad_response, closure),
                "wrong_nonzero_response_residual",
            ),
            verify_response_radius_requires_regular_denominators(),
            verify_control_only_response(),
        ],
    }


def certify(operators, closure, max_atan_terms):
    matrices = parse_operator_export(operators)
    parsed_closure = parse_closure(closure)
    terminal = assert_terminal_residuals(matrices, parsed_closure)

    denominators = collect_denominators([
        ("operators_M0", matrices[0]),
        ("operators_M1", matrices[1]),
        ("operators_M2", matrices[2]),
        ("P", terminal["P"]),
        ("S", terminal["S"]),
        ("profile_acceleration_state", terminal["acceleration_state"]),
        ("profile_acceleration_control", terminal["acceleration_control"]),
        ("U", terminal["U"]),
        ("U1", terminal["U1"]),
        ("U2", terminal["U2"]),
        ("Kaug", terminal["Kaug"]),
        ("response", parsed_closure["response"]),
    ])
    denominator_certificate = certify_polys_at_r0(
        denominators, s.Rational(1), max_atan_terms
    )
    response_certificate = certify_response(
        parsed_closure, denominators, max_atan_terms
    )

    return {
        "state": "terminal germ certificate passed",
        "epsilon0": "1/100",
        "T": "atan(1/(10*r))",
        "T_prime": str(T_PRIME),
        "free_state_count": parsed_closure["free_count"],
        "constant_control_count": parsed_closure["control_count"],
        "all_22_retained_rows_zero": True,
        "x0_equation_included_in_22_retained_rows": True,
        "discarded_response_matches_reducer_serialization": True,
        "constant_control_augmented_state_supported": True,
        "smooth_solution_germ_certified": True,
        "formal_discarded_response_nonzero": response_certificate[
            "formal_response_nonzero"
        ],
        "discarded_response_nonzero": response_certificate["response_nonzero"],
        "certified_discarded_response_nonzero": response_certificate[
            "certified_response_nonzero"
        ],
        "nonintegrability_conclusion_supported": response_certificate[
            "certified_response_nonzero"
        ],
        "nonintegrability_target_if_response_nonzero": (
            "no C2 amplitude family of exact retained bulk solutions on a common "
            "spacetime neighborhood, with topology controlling the required Euler "
            "jets, can have the specified spacetime tangent e+g_f"
        ),
        "instantaneous_configuration_only_claim": False,
        "nonlinear_consistency_established": False,
        "search_limitation": (
            "terminal denominator regularity at r=1 is required before "
            "alternate response radii are tried"
        ),
        "denominator_certificate": denominator_certificate_summary(
            denominators, denominator_certificate
        ),
        "response_certificate": response_certificate,
    }


def run(config_path, operators_path, closure_path, output_path, max_atan_terms):
    if not __debug__:
        raise RuntimeError("scientific assertions must remain enabled")
    config_path = Path(config_path)
    operators_path = Path(operators_path)
    closure_path = Path(closure_path)

    config, _, config_sha256 = load_json_bytes(config_path)
    operators, _, operators_sha256 = load_pinned_json(
        operators_path,
        expected_sha256(config, "stationary_operators"),
        "stationary operators",
    )
    closure, _, closure_sha256 = load_pinned_json(
        closure_path,
        expected_sha256(config, "dae_closure"),
        "DAE closure",
    )
    actual_hashes = {
        "config": config_sha256,
        "stationary_operators": operators_sha256,
        "dae_closure": closure_sha256,
    }
    manufactured_controls = verify_manufactured_controls()
    result = certify(operators, closure, max_atan_terms)
    atomic_json_new(output_path, {
        **result,
        "input_sha256": actual_hashes,
        "manufactured_controls": manufactured_controls,
    })
    print("PASS: terminal retained rows and local germ interval certified.", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--operators", type=Path, required=True)
    parser.add_argument("--closure", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-atan-terms", type=int, default=256)
    args = parser.parse_args()
    run(
        args.config,
        args.operators,
        args.closure,
        args.output,
        args.max_atan_terms,
    )


if __name__ == "__main__":
    main()

"""Independent coefficient controls, executed only inside the tracked atlas check.

The callback supplies already reconstructed columns; no projection is rerun here.
Expected low harmonics and the nine archived axisymmetric matrix entries are
written explicitly rather than obtained from the production generator.
"""
import sympy as s

from atlas_core import Label, harmonics, x


AXES = ("x", "y", "z")
FIELD_MINIMUM = {
    "a": 0, "b": 0, "lambda": 0, "X": 0,
    "V_polar": 1, "V_axial": 1, "W_polar": 1, "W_axial": 1,
    "C_polar": 2, "C_axial": 2,
}
ARCHIVED_SOURCE = {
    "path": "paper1/reports/oblate-null-map-check-2026-11-27/sources/"
            "verification/oblate_nonlinear_null_20261127/null_compatibility.py",
    "function": "expected_discarded_generator",
    "convention": "unit w; eight output rows, nine nonzero matrix entries",
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def equal(actual, expected, message):
    require(s.cancel(actual - expected) == 0, message)


def expected_columns():
    return [Label(field, ell, m)
            for field, lower in FIELD_MINIMUM.items()
            for ell in range(lower, 3) for m in range(-ell, ell + 1)]


def expected_rows():
    return {Label(field, ell, m) for field in FIELD_MINIMUM
            for ell in (3, 4) for m in range(-ell, ell + 1)}


def conjugate_label(label):
    return Label(label.field, label.ell, -label.m)


def check_low_harmonics():
    """Pin real and complex polarizations in the factored amplitude convention."""
    d = 1 - x**2
    zero = s.zeros(2)
    fixtures = [
        (1, 0, (x, s.Matrix([d, 0]), s.Matrix([0, -1]), zero, zero)),
        (2, 0, ((3*x*x - 1)/2, s.Matrix([3*x*d, 0]),
                s.Matrix([0, -3*x]), s.diag(s.Rational(3, 2), -3*d*d/2),
                s.Matrix([[0, -3*d/2], [-3*d/2, 0]]))),
    ]
    for m in (-1, 1):
        fixtures.append((1, m, (
            -s.S.One, s.Matrix([x, -s.I*m/d]),
            s.Matrix([-s.I*m, -x/d]), zero, zero)))
    for m in (-2, 2):
        fixtures.append((2, m, (
            3*d, s.Matrix([-6*x*d, 3*s.I*m]),
            s.Matrix([3*s.I*m*d, 6*x]),
            s.Matrix([[3*(1+x*x)/d, -3*s.I*m*x],
                      [-3*s.I*m*x, -3*d*(1+x*x)]]),
            s.Matrix([[-3*s.I*m*x/d, -3*(1+x*x)],
                      [-3*(1+x*x), 3*s.I*m*x*d]]))))
    for ell, m, expected in fixtures:
        actual = harmonics(ell, m)
        equal(actual[0], expected[0], f"scalar basis anchor {ell},{m}")
        for kind, a, b in zip(("polar", "axial", "E", "B"), actual[1:], expected[1:]):
            require(a.shape == b.shape, f"basis shape {ell},{m},{kind}")
            for av, bv in zip(a, b):
                equal(av, bv, f"basis anchor {ell},{m},{kind}")
    return {"fixture_count": len(fixtures), "signed_modes": [-2, -1, 0, 1, 2]}


def archived_entries():
    # (input, output, coefficient), transcribed from the archived eight-row map.
    return (
        (Label("C_axial", 2, 0), Label("X", 4, 0), -s.Rational(12, 35)),
        (Label("C_axial", 2, 0), Label("C_polar", 4, 0), s.Rational(2, 35)),
        (Label("X", 1, 0), Label("C_axial", 3, 0), -s.Rational(2, 15)),
        (Label("X", 2, 0), Label("C_axial", 4, 0), -s.Rational(2, 35)),
        (Label("C_polar", 2, 0), Label("C_axial", 4, 0), -s.Rational(2, 35)),
        (Label("V_polar", 1, 0), Label("V_axial", 3, 0), -s.Rational(2, 15)),
        (Label("V_polar", 2, 0), Label("V_axial", 4, 0), -s.Rational(6, 35)),
        (Label("W_polar", 1, 0), Label("W_axial", 3, 0), -s.Rational(2, 15)),
        (Label("W_polar", 2, 0), Label("W_axial", 4, 0), -s.Rational(6, 35)),
    )


def check_axisymmetric(mapped_column):
    columns = [label for label in expected_columns() if label.m == 0]
    expected = {label: {} for label in columns}
    for column, row, value in archived_entries():
        expected[column][row] = value
    require(len(columns) == 22, "axisymmetric input count")
    for column in columns:
        actual = mapped_column(column, "z")
        require(set(actual) == set(expected[column]), f"archived map support {column}")
        for row, value in expected[column].items():
            equal(actual[row], value, f"archived map coefficient {column} -> {row}")
    return {"input_columns": 22, "nonzero_output_rows": 8, "nonzero_entries": 9,
            "all_other_entries_zero": True, "temporal_W_rows_included": True,
            "round_trace_X_from_C_axial_included": True}


def check_catalog_and_selection(mapped_column):
    columns, rows = expected_columns(), expected_rows()
    require(len(columns) == len(set(columns)) == 78, "input label count")
    require(len(rows) == 160, "discarded label count")
    entry_count = 0
    for column in columns:
        for axis in AXES:
            actual = mapped_column(column, axis)
            require(isinstance(actual, dict), "map callback must return a dictionary")
            require(set(actual) <= rows, f"discarded label outside catalog {column},{axis}")
            allowed_modes = {column.m} if axis == "z" else {column.m - 1, column.m + 1}
            for row, value in actual.items():
                require(row.m in allowed_modes, f"Fourier selection {column},{axis},{row}")
                real, imag = s.expand_complex(value).as_real_imag()
                require(isinstance(real, s.Rational) and isinstance(imag, s.Rational),
                        "coefficient outside exact Q(i)")
                require(value != 0 and s.cancel(value) != 0, "explicit zero in sparse map")
                entry_count += 1
    return {"input_columns": 78, "discarded_labels": 160,
            "generator_columns": 234, "nonzero_entries": entry_count}


def check_reality(mapped_column):
    comparisons = 0
    for column in expected_columns():
        for axis in AXES:
            actual = mapped_column(column, axis)
            opposite = mapped_column(conjugate_label(column), axis)
            expected = {conjugate_label(row): s.conjugate(value)
                        for row, value in actual.items()}
            require(set(opposite) == set(expected), f"reality support {column},{axis}")
            for row, value in expected.items():
                equal(opposite[row], value, f"reality coefficient {column},{axis},{row}")
                comparisons += 1
    return {"column_pairs_checked": 234, "coefficient_comparisons": comparisons,
            "negative_m_convention": "conjugate positive mode; no factorial reweighting"}


def check_transverse_phases(mapped_column):
    comparisons = 0
    for column in expected_columns():
        x_map, y_map = mapped_column(column, "x"), mapped_column(column, "y")
        require(set(x_map) == set(y_map), f"transverse support {column}")
        for row, value in x_map.items():
            k = row.m - column.m
            require(k in (-1, 1), f"transverse branch {column},{row}")
            equal(y_map[row], -s.I*k*value, f"Cartesian y phase {column},{row}")
            comparisons += 1
    require(comparisons > 0, "transverse comparison is vacuous")
    return {"input_columns": 78, "coefficient_comparisons": comparisons,
            "raising_factor": "-I", "lowering_factor": "I"}


def check_rotations(mapped_column):
    count = 0
    for column in expected_columns():
        for axis in AXES:
            require(mapped_column(column, axis, rotation=True) == {},
                    f"rotation has discarded rows {column},{axis}")
            count += 1
    return {"rotation_columns_checked": count, "all_discarded_zero": True}


def expect_rejected(check, callback, name):
    try:
        check(callback)
    except AssertionError:
        return {"case": name, "passed": True}
    raise AssertionError(f"Deliberate defect accepted: {name}")


def mutated_column(mapped_column, target_column, target_axis, target_row, transform):
    def lookup(column, axis, rotation=False):
        result = mapped_column(column, axis, rotation=rotation).copy()
        if column == target_column and axis == target_axis and not rotation:
            transform(result, target_row)
        return result
    return lookup


def check_defect_rejection(mapped_column):
    def perturb(result, row):
        result[row] += 1

    def remove(result, row):
        del result[row]

    def flip(result, row):
        result[row] = -result[row]

    defects = [
        expect_rejected(check_axisymmetric, mutated_column(
            mapped_column, Label("C_axial", 2, 0), "z", Label("X", 4, 0), perturb),
            "round_trace_from_axial_tensor_defect_rejected"),
        expect_rejected(check_axisymmetric, mutated_column(
            mapped_column, Label("W_polar", 1, 0), "z", Label("W_axial", 3, 0), remove),
            "missing_temporal_W_row_rejected"),
        expect_rejected(check_transverse_phases, mutated_column(
            mapped_column, Label("a", 2, 2), "y", Label("a", 3, 3), flip),
            "transverse_y_phase_defect_rejected"),
    ]
    return {"case_count": len(defects), "cases": defects}


def verify_independent_controls(mapped_column):
    cases = []
    for name, control in (
        ("explicit_low_harmonic_basis", lambda: check_low_harmonics()),
        ("catalog_and_fourier_selection", lambda: check_catalog_and_selection(mapped_column)),
        ("archived_axisymmetric_eight_rows", lambda: check_axisymmetric(mapped_column)),
        ("all_column_complex_conjugacy", lambda: check_reality(mapped_column)),
        ("all_column_cartesian_transverse_phases", lambda: check_transverse_phases(mapped_column)),
        ("all_column_rotation_zero", lambda: check_rotations(mapped_column)),
        ("deliberate_coefficient_defects", lambda: check_defect_rejection(mapped_column)),
    ):
        evidence = control()
        cases.append({"case": name, "passed": True, "evidence": evidence})
    return {"status": "independent atlas coefficient controls passed",
            "case_count": len(cases), "cases": cases,
            "archived_source": ARCHIVED_SOURCE,
            "scope": "exact coefficient maps and conventions; no Euler pairing or seed"}

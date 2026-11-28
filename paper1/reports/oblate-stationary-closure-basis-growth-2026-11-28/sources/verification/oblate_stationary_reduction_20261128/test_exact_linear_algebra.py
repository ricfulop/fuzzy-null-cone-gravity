"""Prepared manufactured controls for the optional exact row operations."""
import unittest
import io

from sympy import Float, Matrix, cancel, zeros

from exact_linear_algebra import R, T, independent_extension, rref_with_side


class ExactRowControls(unittest.TestCase):
    def assert_matrix_equal(self, actual, expected):
        self.assertEqual(actual.shape, expected.shape)
        self.assertTrue(all(cancel(value) == 0 for value in actual - expected))

    def test_control_pivot_state_pivot_and_free_control(self):
        a, b = (R + 1)/(T + 2), (T + 1)/(R + 2)
        left = Matrix([[a, 2*a], [0, 0], [2*a, 4*a]])
        side = Matrix([[b, 1], [0, 2], [2*b, 4]])
        controls, states, pivots = rref_with_side(left, side)
        self.assertEqual(pivots, [0])
        self.assert_matrix_equal(controls, Matrix([[1, 2], [0, 0], [0, 0]]))
        self.assert_matrix_equal(states, Matrix([[b/a, 0], [0, 1], [0, 0]]))
        # On the remaining constraint z_2=0 the original and reduced systems
        # admit the same two arbitrary data (u_2,z_1).
        lift = Matrix([[-2, -b/a], [1, 0], [0, 1], [0, 0]])
        self.assert_matrix_equal(left.row_join(side)*lift, zeros(3, 2))
        self.assert_matrix_equal(controls.row_join(states)*lift, zeros(3, 2))

    def test_no_controls_keeps_all_state_constraints(self):
        a, b, pivots = rref_with_side(zeros(2, 0), Matrix([[0, T], [0, 2*T]]))
        self.assertEqual((a.shape, pivots), ((2, 0), []))
        self.assert_matrix_equal(b, Matrix([[0, 1], [0, 0]]))

    def test_zero_state_width_and_empty_shapes(self):
        a, b, pivots = rref_with_side(Matrix([[R], [0]]), zeros(2, 0))
        self.assert_matrix_equal(a, Matrix([[1], [0]]))
        self.assertEqual((b.shape, pivots), ((2, 0), [0]))
        for rows, nc, ns in [(0, 2, 3), (4, 0, 0), (2, 1, 2)]:
            a, b, pivots = rref_with_side(zeros(rows, nc), zeros(rows, ns))
            self.assertEqual((a.shape, b.shape, pivots), ((rows, nc), (rows, ns), []))

    def test_original_candidate_rows_in_order(self):
        existing = Matrix([[R, 0, 0], [T, 0, 0]])
        candidates = Matrix([[1, 0, 0], [0, T, 0], [0, 2*T, 0], [1, 2, R]])
        self.assert_matrix_equal(independent_extension(existing, candidates),
                                 candidates.extract([1, 3], range(3)))
        self.assert_matrix_equal(independent_extension(existing, existing), zeros(0, 3))
        self.assert_matrix_equal(independent_extension(zeros(0, 3), candidates),
                                 candidates.extract([0, 1, 3], range(3)))
        self.assert_matrix_equal(independent_extension(zeros(0, 0), zeros(2, 0)), zeros(0, 0))

    def test_inexact_input_rejected(self):
        with self.assertRaises(ValueError):
            rref_with_side(Matrix([[Float("0.5")]]), zeros(1, 0))


def verify_controls():
    stream = io.StringIO()
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(ExactRowControls)
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    if not result.wasSuccessful():
        raise AssertionError("Exact row-operation controls failed:\n" + stream.getvalue())
    return {"state": "exact row-operation controls passed", "test_count": result.testsRun,
            "control_output": stream.getvalue()}


if __name__ == "__main__":
    unittest.main()

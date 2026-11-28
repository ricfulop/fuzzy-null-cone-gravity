"""Stdlib serialization controls; no scientific module or calculation is loaded."""

import json
import sys
import unittest

from exact_decimal import integer_decimal_string, rational_decimal_string


class ExactDecimalTests(unittest.TestCase):
    def test_zero_sign_and_chunk_boundaries(self):
        for value in (0, 1, -1, 999_999_999, 1_000_000_000,
                      1_000_000_001, -1_000_000_001, 10**27 + 7):
            self.assertEqual(integer_decimal_string(value), str(value))

    def test_above_default_limit_preserves_exact_digits(self):
        value = 10**5000 + 7
        expected = "1" + "0" * 4999 + "7"
        self.assertEqual(integer_decimal_string(value), expected)
        self.assertEqual(integer_decimal_string(-value), "-" + expected)
        self.assertEqual(integer_decimal_string(10**5000 - 1), "9" * 5000)

    def test_repeated_decimal_pattern(self):
        value = 1234567890 * ((10**5000 - 1) // (10**10 - 1))
        self.assertEqual(integer_decimal_string(value), "1234567890" * 500)

    def test_large_rational_text_and_json(self):
        value = 10**5000 + 7
        digits = "1" + "0" * 4999 + "7"
        text = rational_decimal_string(-value, value + 2)
        self.assertEqual(text, "-" + digits + "/1" + "0" * 4999 + "9")
        payload = {"determinant_interval": [text, text], "pivot_intervals": [[text, text]]}
        self.assertEqual(json.loads(json.dumps(payload)), payload)
        self.assertEqual(rational_decimal_string(value, 1), digits)

    def test_digit_limit_is_unchanged(self):
        before = sys.get_int_max_str_digits() if hasattr(sys, "get_int_max_str_digits") else None
        integer_decimal_string(10**6000)
        after = sys.get_int_max_str_digits() if hasattr(sys, "get_int_max_str_digits") else None
        self.assertEqual(before, after)

    def test_inexact_boolean_and_invalid_denominator_rejected(self):
        for value in (True, 1.0, "1", None):
            with self.assertRaises(TypeError):
                integer_decimal_string(value)
            with self.assertRaises(TypeError):
                rational_decimal_string(value, 1)
            with self.assertRaises(TypeError):
                rational_decimal_string(1, value)
        for denominator in (0, -1):
            with self.assertRaises(ValueError):
                rational_decimal_string(1, denominator)


if __name__ == "__main__":
    unittest.main()

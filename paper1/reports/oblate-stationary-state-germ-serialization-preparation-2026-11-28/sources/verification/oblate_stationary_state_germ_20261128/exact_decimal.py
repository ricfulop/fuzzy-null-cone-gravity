"""Exact decimal output for already-computed integers, using bounded chunks.

This module does not alter Python's integer input/output digit-limit setting.
It accepts only integers supplied by the certificate's exact arithmetic.
"""


def integer_decimal_string(value):
    if type(value) is not int:
        raise TypeError("Expected a Python integer")
    if value == 0:
        return "0"
    sign = "-" if value < 0 else ""
    value = abs(value)
    chunks = []
    while value:
        value, remainder = divmod(value, 1_000_000_000)
        chunks.append(remainder)
    return sign + str(chunks[-1]) + "".join(
        f"{chunk:09d}" for chunk in reversed(chunks[:-1])
    )


def rational_decimal_string(numerator, denominator):
    if type(numerator) is not int or type(denominator) is not int:
        raise TypeError("Expected integer numerator and denominator")
    if denominator <= 0:
        raise ValueError("Expected a positive denominator")
    result = integer_decimal_string(numerator)
    return result if denominator == 1 else result + "/" + integer_decimal_string(denominator)

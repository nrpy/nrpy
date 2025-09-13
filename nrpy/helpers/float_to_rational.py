"""
Store helper function f2r needed by sebob to convert floating point numbers to sympy rationals.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed modules:
from typing import Union

import sympy as sp


def f2r(
    input_float: Union[float, str],
    zpad: int = 60,
) -> Union[float, sp.Rational]:
    """
    Convert a float-like number to a high-precision sympy.Rational number.

    This helper reduces decimal-to-binary rounding noise before constructing
    a Rational. It is useful when building closed-form targets or attachment
    conditions derived from the BOB model, where tiny float errors can leak
    into symbolic solves and generated C code.

    :param input_float: Float-like input. Accepts float, mpmath.mpf or string.
    :param zpad: Number of zeros to append after the decimal point (default is 60).
    :return: A sympy Rational number with high precision.

    >>> f2r(0.1)
    1/10
    >>> f2r(1.5)
    3/2
    >>> f2r("2.0",zpad=2)
    2
    """
    float_as_string = str(input_float)

    # Ensure the string has a decimal point
    if "." not in float_as_string:
        float_as_string = f"{float_as_string}."

    # Append 60 zeros after the decimal of the floating point number to increase precision
    return sp.Rational(float_as_string + "0" * zpad)

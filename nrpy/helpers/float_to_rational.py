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


def f2r(input_float: float, do_nothing: bool = False) -> Union[float, sp.Rational]:
    """
    Convert a floating-point number to a high-precision rational number.

    This function takes a floating-point number, converts it to a string,
    and appends 60 zeros to increase the precision of the conversion to a rational number.

    :param input_float: The floating-point number to convert.
    :param do_nothing: Boolean flag to return the input float (for debugging, default is False).
    :return: Original float if do_nothing is True, else a sympy Rational number with high precision.

    >>> f2r(0.1)
    1/10
    >>> f2r(1.5)
    3/2
    >>> f2r(2.0,do_nothing=True)
    2.0
    """
    # if do_nothing is True, return the input float
    if do_nothing:
        return input_float
    # Convert the input float to a string
    float_as_string = str(input_float)

    # Ensure the string has a decimal point
    if "." not in float_as_string:
        float_as_string = f"{float_as_string}."

    # Append 60 zeros after the decimal of the floating point number to increase precision
    return sp.Rational(float_as_string + "0" * 60)

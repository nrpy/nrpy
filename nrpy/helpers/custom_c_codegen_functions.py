# nrpy/helpers/custom_c_codegen_functions.py
"""
Custom functions to override SymPy's standard use of POW.

In place of the standard C POW function, we substitute POW in favor of using
computationally more efficient C functions sqrt, cbrt, and standard multiplication
when possible. This mapping also selects the correct C function based on the
floating point type (float, double, long double), including support for
double complex via the C99 complex math functions (e.g., csqrt, cpow, cabs).
Note: there is no standard complex cbrt in C99, so cube roots of complex values
fall back to cpow.
Note: a negative integer power -n is printed as 1.0/(x*...*x) with n factors, so
the product overflows for |x| above DBL_MAX^(1/n) (about 1.3e154 for n=2), and the
quotient then returns 0 where the exact result is a representable number. The
single-precision limit for n=2 is sqrt(FLT_MAX), about 1.8e19.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import sympy as sp


def _fits_c_int(b: sp.Basic) -> bool:
    """
    Return True if b is a SymPy Integer that fits in a 32-bit C int, the argument type of abs().

    :param b: Argument of the absolute value.
    :return: True if the C function abs() can take b without narrowing.
    """
    return isinstance(b, sp.Integer) and abs(b) < 2**31


custom_functions_for_SymPy_ccode = {
    "double": {
        "nrpyAbs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabs({b})"),
        ],
        "fabs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabs({b})"),
        ],
        "Pow": [
            (lambda b, e: e == sp.Rational(1, 2), lambda b, e: f"sqrt({b})"),
            (lambda b, e: e == 0.5, lambda b, e: f"sqrt({b})"),
            (lambda b, e: e == -sp.Rational(1, 2), lambda b, e: f"(1.0/sqrt({b}))"),
            (lambda b, e: e == -0.5, lambda b, e: f"(1.0/sqrt({b}))"),
            (lambda b, e: e == sp.S.One / 3, lambda b, e: f"cbrt({b})"),
            (lambda b, e: e == -sp.S.One / 3, lambda b, e: f"(1.0/cbrt({b}))"),
            (lambda b, e: e == sp.Rational(1, 6), lambda b, e: f"(sqrt(cbrt({b})))"),
            (
                lambda b, e: e == -sp.Rational(1, 6),
                lambda b, e: f"(1.0/(sqrt(cbrt({b}))))",
            ),
            (lambda b, e: e == 2, lambda b, e: f"(({b})*({b}))"),
            (lambda b, e: e == 3, lambda b, e: f"(({b})*({b})*({b}))"),
            (lambda b, e: e == 4, lambda b, e: f"(({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == 5, lambda b, e: f"(({b})*({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == -1, lambda b, e: f"(1.0/({b}))"),
            (lambda b, e: e == -2, lambda b, e: f"(1.0/(({b})*({b})))"),
            (lambda b, e: e == -3, lambda b, e: f"(1.0/(({b})*({b})*({b})))"),
            (lambda b, e: e == -4, lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})))"),
            (
                lambda b, e: e == -5,
                lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})*({b})))",
            ),
            (lambda b, e: e != -5, "pow"),
        ],
    },
    "float": {
        "nrpyAbs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsf({b})"),
        ],
        "fabs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsf({b})"),
        ],
        "Pow": [
            (lambda b, e: e == sp.Rational(1, 2), lambda b, e: f"sqrtf({b})"),
            (lambda b, e: e == 0.5, lambda b, e: f"sqrtf({b})"),
            (lambda b, e: e == -sp.Rational(1, 2), lambda b, e: f"(1.0f/sqrtf({b}))"),
            (lambda b, e: e == -0.5, lambda b, e: f"(1.0f/sqrtf({b}))"),
            (lambda b, e: e == sp.S.One / 3, lambda b, e: f"cbrtf({b})"),
            (lambda b, e: e == -sp.S.One / 3, lambda b, e: f"(1.0f/cbrtf({b}))"),
            (lambda b, e: e == 2, lambda b, e: f"(({b})*({b}))"),
            (lambda b, e: e == 3, lambda b, e: f"(({b})*({b})*({b}))"),
            (lambda b, e: e == 4, lambda b, e: f"(({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == 5, lambda b, e: f"(({b})*({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == -1, lambda b, e: f"(1.0f/({b}))"),
            (lambda b, e: e == -2, lambda b, e: f"(1.0f/(({b})*({b})))"),
            (lambda b, e: e == -3, lambda b, e: f"(1.0f/(({b})*({b})*({b})))"),
            (lambda b, e: e == -4, lambda b, e: f"(1.0f/(({b})*({b})*({b})*({b})))"),
            (
                lambda b, e: e == -5,
                lambda b, e: f"(1.0f/(({b})*({b})*({b})*({b})*({b})))",
            ),
            (lambda b, e: e != -5, "powf"),
        ],
    },
    "long double": {
        "nrpyAbs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsl({b})"),
        ],
        "fabs": [
            (_fits_c_int, lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsl({b})"),
        ],
        "Pow": [
            (lambda b, e: e == sp.Rational(1, 2), lambda b, e: f"sqrtl({b})"),
            (lambda b, e: e == 0.5, lambda b, e: f"sqrtl({b})"),
            (lambda b, e: e == -sp.Rational(1, 2), lambda b, e: f"(1.0l/sqrtl({b}))"),
            (lambda b, e: e == -0.5, lambda b, e: f"(1.0l/sqrtl({b}))"),
            (lambda b, e: e == sp.S.One / 3, lambda b, e: f"cbrtl({b})"),
            (lambda b, e: e == -sp.S.One / 3, lambda b, e: f"(1.0l/cbrtl({b}))"),
            (lambda b, e: e == 2, lambda b, e: f"(({b})*({b}))"),
            (lambda b, e: e == 3, lambda b, e: f"(({b})*({b})*({b}))"),
            (lambda b, e: e == 4, lambda b, e: f"(({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == 5, lambda b, e: f"(({b})*({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == -1, lambda b, e: f"(1.0l/({b}))"),
            (lambda b, e: e == -2, lambda b, e: f"(1.0l/(({b})*({b})))"),
            (lambda b, e: e == -3, lambda b, e: f"(1.0l/(({b})*({b})*({b})))"),
            (lambda b, e: e == -4, lambda b, e: f"(1.0l/(({b})*({b})*({b})*({b})))"),
            (
                lambda b, e: e == -5,
                lambda b, e: f"(1.0l/(({b})*({b})*({b})*({b})*({b})))",
            ),
            (lambda b, e: e != -5, "powl"),
        ],
    },
    "double complex": {
        # magnitude for complex -> returns double
        "nrpyAbs": "cabs",
        "fabs": "cabs",
        "Abs": "cabs",
        # unary transcendentals -> force complex variants
        "sin": "csin",
        "cos": "ccos",
        "tan": "ctan",
        "asin": "casin",
        "acos": "cacos",
        "atan": "catan",
        "exp": "cexp",
        "log": "clog",
        "sinh": "csinh",
        "cosh": "ccosh",
        "tanh": "ctanh",
        "asinh": "casinh",
        "acosh": "cacosh",
        "atanh": "catanh",
        # powers (keep the predicate-list form only for Pow)
        "Pow": [
            (lambda b, e: e == sp.Rational(1, 2), lambda b, e: f"csqrt({b})"),
            (lambda b, e: e == 0.5, lambda b, e: f"csqrt({b})"),
            (lambda b, e: e == -sp.Rational(1, 2), lambda b, e: f"(1.0/csqrt({b}))"),
            (lambda b, e: e == -0.5, lambda b, e: f"(1.0/csqrt({b}))"),
            # no standard complex cbrt; use cpow
            (lambda b, e: e == sp.S.One / 3, lambda b, e: f"cpow({b}, 1.0/3.0)"),
            (lambda b, e: e == -sp.S.One / 3, lambda b, e: f"(1.0/cpow({b}, 1.0/3.0))"),
            (
                lambda b, e: e == sp.Rational(1, 6),
                lambda b, e: f"csqrt(cpow({b}, 1.0/3.0))",
            ),
            (
                lambda b, e: e == -sp.Rational(1, 6),
                lambda b, e: f"(1.0/csqrt(cpow({b}, 1.0/3.0)))",
            ),
            # small integer powers
            (lambda b, e: e == 2, lambda b, e: f"(({b})*({b}))"),
            (lambda b, e: e == 3, lambda b, e: f"(({b})*({b})*({b}))"),
            (lambda b, e: e == 4, lambda b, e: f"(({b})*({b})*({b})*({b}))"),
            (lambda b, e: e == 5, lambda b, e: f"(({b})*({b})*({b})*({b})*({b}))"),
            # negative integer powers
            (lambda b, e: e == -1, lambda b, e: f"(1.0/({b}))"),
            (lambda b, e: e == -2, lambda b, e: f"(1.0/(({b})*({b})))"),
            (lambda b, e: e == -3, lambda b, e: f"(1.0/(({b})*({b})*({b})))"),
            (lambda b, e: e == -4, lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})))"),
            (
                lambda b, e: e == -5,
                lambda b, e: f"(1.0/(({b})*({b})*({b})*({b})*({b})))",
            ),
            # general fallback
            (lambda b, e: True, "cpow"),
        ],
    },
}

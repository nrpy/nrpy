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

Author(s):  Zachariah B. Etienne
            zachetie **at** gmail **dot* com
            Samuel D. Tootle
            sdtootle **at** gmail **dot* com
"""

import sympy as sp

custom_functions_for_SymPy_ccode = {
    "double": {
        "nrpyAbs": [
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabs({b})"),
        ],
        "fabs": [
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
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
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsf({b})"),
        ],
        "fabs": [
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
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
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
            (lambda b: True, lambda b: f"fabsl({b})"),
        ],
        "fabs": [
            (lambda b: isinstance(b, sp.Integer), lambda b: f"abs({b})"),
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
        # unary transcendentals -> force complex variants
        "sin": "csin",
        "cos": "ccos",
        "tan": "ctan",
        "exp": "cexp",
        "log": "clog",
        "sinh": "csinh",
        "cosh": "ccosh",
        "tanh": "ctanh",
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

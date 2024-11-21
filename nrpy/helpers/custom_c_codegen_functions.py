"""
Custom functions to override SymPy's standard use of POW.
In place of the standard C POW function, we substitute POW
in favor of using computationally more efficient C functions sqrt, cbrt, and standard multiplication.
This has been expanded to use the correct C function based on the floating point
type (e.g. float, double, long double)

Author(s):  Zachariah B. Etienne
            zachetie **at** gmail **dot* com
            Samuel D. Tootle
            sdtootle **at** gmail **dot* com
"""

import sympy as sp

custom_functions_for_SymPy_ccode = {
    "double": {
        "nrpyAbs": "fabs",
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
        "nrpyAbs": "fabsf",
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
        "nrpyAbs": "fabsl",
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
}

"""
Construct symbolic expressions of min and max functions, later replacing.

calls of nrpyAbs with fabs calls

"""

from typing import Any, Dict

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.params as par  # NRPy: parameter interface

_nrpyAbs = sp.Function("nrpyAbs")


def register_TINYDOUBLE_if_needed() -> sp.Expr:
    """
    Register TINYDOUBLE if it hasn't been registered, and return the symbol.

    :return: TINYDOUBLE sympy symbol.
    """
    if "TINYDOUBLE" not in par.glb_code_params_dict:
        par.register_CodeParameter(
            "#define", __name__, "TINYDOUBLE", 1e-100, commondata=False
        )
    return sp.Symbol("TINYDOUBLE", real=True)


def min_noif(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    """
    Return the minimum of a and b.

    :param a: sympy expression or number
    :param b: sympy expression or number.

    :return: minimum of the two in symbolic form.
    """
    if a == sp.sympify(0):
        return sp.Rational(1, 2) * (b - _nrpyAbs(b))
    if b == sp.sympify(0):
        return sp.Rational(1, 2) * (a - _nrpyAbs(a))
    return sp.Rational(1, 2) * (a + b - _nrpyAbs(a - b))


def max_noif(a: sp.Expr, b: sp.Expr) -> sp.Expr:
    """
    Return the maximum of a and b.

    :param a: sympy expression or number
    :param b: sympy expression or number.

    :return: maximum of the two in symbolic form
    """
    if a == sp.sympify(0):
        return sp.Rational(1, 2) * (b + _nrpyAbs(b))
    if b == sp.sympify(0):
        return sp.Rational(1, 2) * (a + _nrpyAbs(a))
    return sp.Rational(1, 2) * (a + b + _nrpyAbs(a - b))


def coord_leq_bound(x: sp.Expr, xstar: sp.Expr) -> sp.Expr:
    """
    Return 1.0 if x <= xstar, 0.0 otherwise.

    :param x: sympy variable
    :param xstar: bound for variable.

    :return: symbolic form x <= xstar
    """
    TINYDOUBLE = register_TINYDOUBLE_if_needed()
    return min_noif(x - xstar - TINYDOUBLE, 0.0) / (x - xstar - TINYDOUBLE)


def coord_geq_bound(x: sp.Expr, xstar: sp.Expr) -> sp.Expr:
    """
    Return 1.0 if x >= xstar, 0.0 otherwise.

    :param x: sympy variable
    :param xstar: bound for variable.

    :return: symbolic form x >= xstar
    """
    TINYDOUBLE = register_TINYDOUBLE_if_needed()
    return max_noif(x - xstar + TINYDOUBLE, 0.0) / (x - xstar + TINYDOUBLE)


def coord_less_bound(x: sp.Expr, xstar: sp.Expr) -> sp.Expr:
    """
    Return 1.0 if x < xstar, 0.0 otherwise.

    :param x: sympy variable
    :param xstar: bound for variable.

    :return: symbolic form x > xstar
    """
    TINYDOUBLE = register_TINYDOUBLE_if_needed()
    return min_noif(x - xstar, 0.0) / (x - xstar - TINYDOUBLE)


def coord_greater_bound(x: sp.Expr, xstar: sp.Expr) -> sp.Expr:
    """
    Return 1.0 if x > xstar, 0.0 otherwise.

    :param x: sympy variable
    :param xstar: bound for variable.

    :return: symbolic form x > xstar
    """
    TINYDOUBLE = register_TINYDOUBLE_if_needed()
    return max_noif(x - xstar, 0.0) / (x - xstar + TINYDOUBLE)


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    exprs_dict: Dict[str, Any] = {}

    y = sp.symbols("y", real=True)
    y_star = sp.sympify(0.00125)

    z = sp.symbols("z", real=True)
    z_star = sp.sympify(0.125)

    exprs_dict["min_noif_exp"] = min_noif(y, z).subs(sp.Function("nrpyAbs"), sp.sin)
    exprs_dict["max_noif_exp"] = max_noif(z, y).subs(sp.Function("nrpyAbs"), sp.cos)

    exprs_dict["coord_leq_bound_exp"] = coord_leq_bound(y, y_star).subs(
        sp.Function("nrpyAbs"), sp.exp
    )
    exprs_dict["coord_geq_bound_exp"] = coord_geq_bound(z, z_star).subs(
        sp.Function("nrpyAbs"), sp.sin
    )

    exprs_dict["coord_less_bound_exp"] = coord_less_bound(z, y).subs(
        sp.Function("nrpyAbs"), sp.cos
    )
    exprs_dict["coord_greater_bound_exp"] = coord_greater_bound(y, z_star).subs(
        sp.Function("nrpyAbs"), sp.sin
    )

    results_dict = ve.process_dictionary_of_expressions(
        exprs_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )

"""
Utilities for Parsing or manipulating symbolic expressions.

Author: Samuel Tootle
Email:  sdtootle *at* gmail *dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import List, Tuple, Union

import sympy as sp

import nrpy.params as par  # NRPy+: Parameter interface


def get_unique_expression_symbols_as_strings(
    expr: sp.Expr, exclude: Union[List[str], None] = None
) -> List[str]:
    """
    Get a unique list of expression symbols.

    :param expr: Sympy expression
    :param exclude: List of symbol names to exclude
    :returns: List of unique symbol names from the expression

    DOCTEST:
    >>> from sympy.abc import a, b
    >>> from sympy import cos, sin
    >>> xx0 = sp.Symbol('xx0')
    >>> x = cos(a + b)**2 ++ xx0
    >>> get_unique_expression_symbols_as_strings(x, ["xx0"])
    ['a', 'b']
    """
    exclude = [] if not exclude else exclude

    symbols = {sym.name for sym in expr.free_symbols}
    return sorted(symbols - set(exclude))


def get_params_commondata_symbols_from_expr_list(
    expr_list: List[sp.Expr], exclude: Union[List[str], None] = None
) -> Tuple[List[str], List[str]]:
    """
    Get the param and commondata symbols from a list of expressions.

    :param expr_list: List of sympy expressions
    :param exclude: List of symbol names to exclude
    :returns: List of params symbols, List of commondata symbols

    """
    exclude = [] if not exclude else exclude
    unique_symbols = []
    for tmp_expr in expr_list:
        unique_symbols += get_unique_expression_symbols_as_strings(
            tmp_expr, exclude=exclude
        )
    unique_symbols = list(set(unique_symbols))

    param_symbols = sorted(
        list(
            set(unique_symbols)
            & set(
                k if not v.commondata else ""
                for k, v in par.glb_code_params_dict.items()
            )
        )
    )

    commondata_symbols = sorted(
        list(
            set(unique_symbols)
            & set(
                k if v.commondata else "" for k, v in par.glb_code_params_dict.items()
            )
        )
    )
    return param_symbols, commondata_symbols


def generate_definition_header(
    str_list: List[str], enable_intrinsics: bool = False, var_access: str = "params->"
) -> str:
    """
    Generate the string header for parameter definitions.

    :param str_list: List of Symbol strings to include in the header
    :param enable_intrinsics: Whether to modify str based on hardware intrinsics.
    :param var_access: The variable access string
    :returns: The definition string
    """
    return "\n".join(
        [
            (
                f"const REAL {p} = {var_access}{p};"
                if not enable_intrinsics
                else f"const REAL NOSIMD{p} = {var_access}{p};\n"
                f"MAYBE_UNUSED const REAL_SIMD_ARRAY {p} = ConstSIMD(NOSIMD{p});\n"
            )
            for p in str_list
        ]
    )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

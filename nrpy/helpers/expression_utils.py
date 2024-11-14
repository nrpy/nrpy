"""
Utilities for Parsing or manipulating symbolic expressions.

Author: Samuel Tootle
Email:  sdtootle *at* gmail *dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import List, Union

import sympy as sp


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


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

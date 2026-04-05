# nrpy/helpers/expression_utils.py
"""
Utilities for Parsing or manipulating symbolic expressions.

Author: Samuel Tootle
Email:  sdtootle *at* gmail *dot* com
"""

import sys  # Standard Python module for multiplatform OS-level functions
from typing import Dict, FrozenSet, List, Set, Tuple, Union, cast

import sympy as sp

import nrpy.params as par  # NRPy: Parameter interface


def _get_free_symbol_names_cached(
    expr: sp.Basic, exclude: Set[str], memo: Dict[int, FrozenSet[str]]
) -> FrozenSet[str]:
    """
    Return free-symbol names for `expr` using a memoized DAG traversal.

    This preserves SymPy's `free_symbols` semantics, including removal of
    locally bound symbols for objects like `Integral`, `Lambda`, and `Subs`,
    while avoiding repeated whole-tree traversals for large shared expression
    graphs.

    :param expr: SymPy expression or subexpression to inspect.
    :param exclude: Symbol names to exclude from the result.
    :param memo: Cache from object id to computed free-symbol names.
    :return: Frozen set of free symbol names appearing in ``expr``.
    """
    expr_id = id(expr)
    if expr_id in memo:
        return memo[expr_id]

    stack: List[Tuple[sp.Basic, bool]] = [(expr, False)]
    while stack:
        node, expanded = stack.pop()
        node_id = id(node)
        if node_id in memo:
            continue

        if expanded:
            names: Set[str] = set()
            for arg in node.args:
                names.update(memo[id(arg)])

            bound_symbols = getattr(node, "bound_symbols", ())
            if bound_symbols:
                names.difference_update(
                    cast(str, getattr(sym, "name", str(sym))) for sym in bound_symbols
                )

            memo[node_id] = frozenset(names)
            continue

        if node.is_Atom:
            names = set()
            for sym in node.free_symbols:
                sym_name = cast(str, getattr(sym, "name", str(sym)))
                if sym_name not in exclude:
                    names.add(sym_name)
            memo[node_id] = frozenset(names)
            continue

        stack.append((node, True))
        for arg in node.args:
            if id(arg) not in memo:
                stack.append((arg, False))

    return memo[expr_id]


def get_unique_expression_symbols_as_strings(
    expr: sp.Expr, exclude: Union[List[str], None] = None
) -> List[str]:
    """
    Get a unique list of expression symbols.

    :param expr: Sympy expression
    :param exclude: List of symbol names to exclude
    :return: List of unique symbol names from the expression

    DOCTEST:
    >>> from sympy.abc import a, b
    >>> from sympy import cos, sin
    >>> xx0 = sp.Symbol('xx0')
    >>> x = cos(a + b)**2 ++ xx0
    >>> get_unique_expression_symbols_as_strings(x, ["xx0"])
    ['a', 'b']
    """
    exclude_set = set(exclude or [])
    return sorted(_get_free_symbol_names_cached(expr, exclude_set, {}))


def get_params_commondata_symbols_from_expr_list(
    expr_list: List[sp.Expr], exclude: Union[List[str], None] = None
) -> Tuple[List[str], List[str]]:
    """
    Get the param and commondata symbols from a list of expressions.

    :param expr_list: List of sympy expressions
    :param exclude: List of symbol names to exclude
    :return: List of params symbols, List of commondata symbols

    """
    exclude_set = set(exclude or [])
    memo: Dict[int, FrozenSet[str]] = {}
    unique_symbols: Set[str] = set()
    for tmp_expr in expr_list:
        unique_symbols.update(
            _get_free_symbol_names_cached(tmp_expr, exclude_set, memo)
        )

    param_symbols: List[str] = []
    commondata_symbols: List[str] = []
    for symbol_name in sorted(unique_symbols):
        code_param = par.glb_code_params_dict.get(symbol_name)
        if code_param is None:
            continue
        if code_param.commondata:
            commondata_symbols.append(symbol_name)
        else:
            param_symbols.append(symbol_name)
    return param_symbols, commondata_symbols


def generate_definition_header(
    str_list: List[str], enable_intrinsics: bool = False, var_access: str = "params->"
) -> str:
    """
    Generate the string header for parameter definitions.

    :param str_list: List of Symbol strings to include in the header
    :param enable_intrinsics: Whether to modify str based on hardware intrinsics.
    :param var_access: The variable access string
    :return: The definition string
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

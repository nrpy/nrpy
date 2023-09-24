""" CSE Partial Factorization and Post-Processing

    The following script will perform partial factorization on SymPy expressions,
    which should occur before common subexpression elimination (CSE) to prevent the
    identification of undesirable patterns, and perform post-processing on the
    the resulting replaced/reduced expressions after the CSE procedure was applied.

Author: Ken Sible
Email:  ksible *at* outlook *dot* com
"""
import sys  # Standard Python module for multiplatform OS-level functions
from collections import OrderedDict
from typing import Union, List, Tuple, Dict, cast
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
from nrpy.helpers.expr_tree import ExprTree


def cse_preprocess(
    expr_list: Union[
        sp.Basic,
        List[sp.Basic],
    ],
    # expecting either a single sympy expression or a list of sympy expressions
    prefix: str = "",  # this string prefix will be used for variable names (i.e. rational symbols)
    declare_neg1_as_symbol: bool = False,  # if true, declares a symbol for negative one (i.e. _NegativeOne_)
    factor: bool = True,  # if true, performs partial factorization (excluding negative symbol)
    negative: bool = False,  # if true, includes negative symbol in partial factorization
    debug: bool = False,  # if true, back-substitutes and checks difference for debugging
) -> Tuple[
    List[sp.Basic],
    Dict[sp.Basic, sp.Rational],
]:  # returns a tuple of modified sympy expression(s) and a dictionary mapping symbols to rational numbers
    """Perform CSE Preprocessing

    :arg:    single SymPy expression or list of SymPy expressions
    :arg:    string prefix for variable names (i.e. rational symbols)
    :arg:    declare symbol for negative one (i.e. _NegativeOne_)
    :arg:    perform partial factorization (excluding negative symbol)
    :arg:    include negative symbol in partial factorization
    :arg:    back-substitute and check difference for debugging
    :return: modified SymPy expression(s) where all integers and rationals were replaced
                with temporary placeholder variables that allow for partial factorization

    Doctests:
    >>> from sympy.abc import x, y, z
    >>> expr = -x/12 - y/12 + z
    >>> cse_preprocess(expr, debug=True)
    ([_Rational_1_12*(-x - y) + z], OrderedDict([(_Rational_1_12, 1/12)]))

    >>> cse_preprocess(expr, declare_neg1_as_symbol=True, debug=True)
    ([_Rational_1_12*(_NegativeOne_*x + _NegativeOne_*y) + z], OrderedDict([(_Rational_1_12, 1/12), (_NegativeOne_, -1)]))

    >>> expr = -x/12 - y/12 + z
    >>> cse_preprocess(expr, declare_neg1_as_symbol=True, negative=True, debug=True)
    ([_NegativeOne_*_Rational_1_12*(x + y) + z], OrderedDict([(_Rational_1_12, 1/12), (_NegativeOne_, -1)]))

    >>> cse_preprocess(expr, factor=False, debug=True)
    ([(-_Rational_1_12)*x + (-_Rational_1_12)*y + z], OrderedDict([(_Rational_1_12, 1/12)]))

    >>> cse_preprocess(expr, prefix='FD', debug=True)
    ([FD_Rational_1_12*(-x - y) + z], OrderedDict([(FD_Rational_1_12, 1/12)]))

    >>> from sympy import exp
    >>> expr = exp(3*x + 3*y)
    >>> cse_preprocess(expr, debug=True)
    ([exp(_Integer_3*(x + y))], OrderedDict([(_Integer_3, 3)]))

    >>> from sympy import Mul
    >>> expr = Mul((-1)**3, (3*x + 3*y), evaluate=False)
    >>> cse_preprocess(expr, declare_neg1_as_symbol=True, debug=True)
    ([_Integer_3*_NegativeOne_*(x + y)], OrderedDict([(_NegativeOne_, -1), (_Integer_3, 3)]))
    """
    # Convert input to list if it's not a list
    if not isinstance(expr_list, list):
        expr_list = [expr_list]

    # Make a copy of the list
    expr_list = expr_list[:]

    # Define negative one symbol
    _NegativeOne_ = sp.Symbol(prefix + "_NegativeOne_")

    # Maps to hold symbol to rational and rational to symbol mappings
    symbol_to_Rational_dict: Dict[sp.Basic, sp.Rational] = OrderedDict()
    map_rat_to_sym: Dict[sp.Rational, sp.Basic] = OrderedDict()

    # Loop over each expression in the input list
    for i, expr in enumerate(expr_list):
        # Create an expression tree for each expression
        tree = ExprTree(expr)

        # Preorder traversal of the expression tree, searching for rational numbers
        for subtree in tree.preorder():
            subexpr = subtree.expr

            # If the subexpression is a Rational type and it's not equal to Negative One
            if isinstance(subexpr, sp.Rational) and subexpr != sp.S.NegativeOne:
                # mypy is clueless here; subexpr is sp.Rational, and mypy complains that it could be something else!

                # Continue loop if the subexpression is an exponent in a power function, we don't want to replace it
                if subtree.func == sp.Pow:
                    continue

                # If rational < 0, factor out negative, leaving positive rational
                sign = 1 if subexpr >= 0 else -1
                subexpr *= sign

                # If rational number hasn't been encountered before, replace
                # it with a symbol using get() to avoid try-except;
                # typehinting note: subexpr is guaranteed to be Rational.
                repl = map_rat_to_sym.get(subexpr)  # type: ignore
                if repl is None:
                    p, q = subexpr.p, subexpr.q  # type: ignore

                    # Name the variable based on its value and whether it's an integer or a rational number
                    var_name = (
                        f"{prefix}_Rational_{p}_{q}"
                        if q != 1
                        else f"{prefix}_Integer_{p}"
                    )

                    # Replace the rational number with the symbol in the expression
                    repl = sp.Symbol(var_name)

                    # Add mapping of symbol to rational and rational to symbol
                    symbol_to_Rational_dict[repl], map_rat_to_sym[subexpr] = subexpr, repl  # type: ignore

                # Update subexpression in the subtree
                subtree.expr = repl * sign

                if sign < 0:
                    tree.build(subtree, clear=False)

            # If declare_neg1_as_symbol is True, replace negative one with symbol
            elif declare_neg1_as_symbol and subexpr == sp.S.NegativeOne:
                # using get() to avoid try-except
                subtree.expr = map_rat_to_sym.get(
                    sp.S.NegativeOne, sp.Symbol("didnotfind_subtree_expr")
                )
                if subtree.expr == sp.Symbol("didnotfind_subtree_expr"):
                    symbol_to_Rational_dict[_NegativeOne_] = sp.S.NegativeOne
                    map_rat_to_sym[subexpr] = _NegativeOne_  # type: ignore
                    subtree.expr = _NegativeOne_
        # Update expression from reconstructed tree
        expr = tree.reconstruct()

        # Perform partial factoring if factor is True
        if factor:
            # Get set of symbols to factor, excluding _NegativeOne_
            var_set = [
                var for var in symbol_to_Rational_dict if var != _NegativeOne_
            ]  # using list comprehension

            # Handle factoring of function argument(s)
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    arg = subtree.children[0]
                    for var in var_set:
                        if var in arg.expr.free_symbols:
                            arg.expr = sp.collect(arg.expr, var)
                    tree.build(arg)

            # Update expression from reconstructed tree
            expr = tree.reconstruct()

            # Perform partial factoring on entire expression
            # This collect is very expensive, so make sure var exists in expr.free_symbols!
            free_symbols = expr.free_symbols
            needed_var_set: List[sp.Basic] = []
            for var in var_set:
                if var in free_symbols:
                    needed_var_set += [var]
            expr = sp.collect(expr, needed_var_set)
            tree.root.expr = expr
            tree.build(tree.root)

        # If negative is True, perform partial factoring on _NegativeOne_
        if negative:
            for subtree in tree.preorder():
                if isinstance(subtree.expr, sp.Function):
                    arg = subtree.children[0]
                    arg.expr = sp.collect(arg.expr, _NegativeOne_)
                    tree.build(arg)
            expr = sp.collect(tree.reconstruct(), _NegativeOne_)
            tree.root.expr = expr
            tree.build(tree.root)

        # If declare is True, simplify (-1)^n
        if declare_neg1_as_symbol:
            changed_expr = False
            _One_ = sp.Symbol(prefix + "_Integer_1")
            for subtree in tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == sp.Pow:
                    base, exponent = subexpr.args[0], subexpr.args[1]
                    if base == _NegativeOne_ and isinstance(exponent, int):
                        subtree.expr = _One_ if exponent % 2 == 0 else _NegativeOne_
                        tree.build(subtree)
                        changed_expr = True
            if changed_expr:
                expr = tree.reconstruct()

        # Replace any remaining ones with symbol _One_ after partial factoring
        if factor or negative:
            _One_ = sp.Symbol(prefix + "_Integer_1")
            for subtree in tree.preorder():
                if subtree.expr == sp.S.One:
                    subtree.expr = _One_
            tmp_expr = tree.reconstruct()
            if tmp_expr != expr:
                # using get() to avoid try-except:
                if map_rat_to_sym.get(sp.S.One) is None:
                    symbol_to_Rational_dict[_One_], map_rat_to_sym[sp.S.One] = (
                        sp.S.One,
                        _One_,
                    )
                    subtree.expr = _One_
                expr = tmp_expr

        # If debug is True, back-substitute everything and check difference
        if debug:
            # Helper function to replace symbols with their corresponding rational numbers
            def lookup_rational(arg: sp.Basic) -> sp.Basic:
                if isinstance(arg, sp.Symbol):
                    arg = symbol_to_Rational_dict.get(arg, arg)
                return arg

            # Create new tree for debugging
            debug_tree = ExprTree(expr)

            # Replace symbols with their corresponding rational numbers in the debug tree
            for subtree in debug_tree.preorder():
                subexpr = subtree.expr
                if subexpr.func == sp.Symbol:
                    subtree.expr = lookup_rational(subexpr)
            debug_expr = tree.reconstruct()

            # Calculate the difference between the original and the debug expression
            if sp.simplify(cast(sp.Expr, expr) - debug_expr) != 0:
                # If the difference is not zero, it means something went wrong in the replacement
                raise ValueError(
                    "Debugging error: Difference in expressions is non-zero."
                )

            # Replace the expression in the list with the new processed expression
        expr_list[i] = expr

    # At the end, return the modified expressions and the dictionary mapping symbols to rational numbers
    return expr_list, symbol_to_Rational_dict


def cse_postprocess(
    cse_output: Tuple[List[Tuple[sp.Symbol, sp.Expr]], List[sp.Expr]]
) -> Tuple[List[Tuple[sp.Symbol, sp.Expr]], List[sp.Expr]]:
    """
    Perform CSE Postprocessing

    This function takes the output from SymPy's Common Subexpression Elimination (CSE), and applies post-processing to the replaced and reduced expressions.
    The post-processing includes handling scalar temporary variables, ordering the replaced expressions, handling negative symbols, and back-substituting in some cases.

    :param cse_output: Output from SymPy CSE with tuple format: (list of ordered pairs that contain substituted symbols and their replaced expressions, reduced SymPy expression)
    :return: Output from SymPy CSE where postprocessing, such as back-substitution of addition/product of symbols, has been applied to the replaced/reduced expression(s)

    Doctests:
    >>> from sympy.abc import x, y
    >>> from sympy import cse, cos, sin

    >>> cse_out = cse(3 + x + cos(3 + x))
    >>> cse_postprocess(cse_out)
    ([], [x + cos(x + 3) + 3])

    >>> cse_out = cse(3 + x + y + cos(3 + x + y))
    >>> cse_postprocess(cse_out)
    ([(x0, x + y + 3)], [x0 + cos(x0)])

    >>> cse_out = cse(3*x + cos(3*x))
    >>> cse_postprocess(cse_out)
    ([], [3*x + cos(3*x)])

    >>> cse_out = cse(3*x*y + cos(3*x*y))
    >>> cse_postprocess(cse_out)
    ([(x0, 3*x*y)], [x0 + cos(x0)])

    >>> cse_out = cse(x**2 + cos(x**2))
    >>> cse_postprocess(cse_out)
    ([], [x**2 + cos(x**2)])

    >>> cse_out = cse(x**3 + cos(x**3))
    >>> cse_postprocess(cse_out)
    ([(x0, x**3)], [x0 + cos(x0)])

    >>> cse_out = cse(x*y + cos(x*y) + sin(x*y))
    >>> cse_postprocess(cse_out)
    ([(x0, x*y)], [x0 + sin(x0) + cos(x0)])

    >>> from sympy import exp, log
    >>> expr = -x + exp(-x) + log(-x)
    >>> cse_pre = cse_preprocess(expr, declare_neg1_as_symbol=True)
    >>> cse_out = cse(cse_pre[0])
    >>> cse_postprocess(cse_out)
    ([], [_NegativeOne_*x + exp(_NegativeOne_*x) + log(_NegativeOne_*x)])
    """
    replaced, reduced = cse_output
    replaced, reduced = replaced[:], reduced[:]

    # SCALAR_TMP's are coming in as sp.Eq. They are
    # effectively "replace" expressions put in by hand.
    # Move them to the replaced array.
    reduced2 = []
    for element in reduced:
        if isinstance(element, sp.Equality):
            replaced += [(element.lhs, element.rhs)]
        else:
            reduced2 += [element]
    reduced = reduced2

    # Sort the replaced expressions
    # so that none are evaluated before
    # they are set.
    # Create a lookup dictionary from the replaced expressions for fast access
    lookup = {}
    for repl in replaced:
        lookup[str(repl[0])] = repl
    replaced2 = []
    # Loop until all expressions are evaluated
    while len(lookup) > 0:
        # Create a new lookup dictionary to store unevaluated expressions
        new_lookup = {}
        for repl in lookup.values():
            # Search through the expression symbols to find any unevaluated ones
            found = False
            for sym in repl[1].free_symbols:
                if str(sym) in lookup:
                    found = True
                    break
            # If found any unevaluated symbols, add this expression to the new lookup dictionary
            if found:
                new_lookup[str(repl[0])] = repl
            else:
                replaced2 += [repl]
        # Verify the new lookup dictionary is smaller than the previous one
        assert len(new_lookup) < len(lookup)
        lookup = new_lookup
    # Verify all expressions have been evaluated
    assert len(replaced) == len(replaced2)
    replaced = replaced2

    # Start a while loop to iterate through all replaced expressions
    i = 0
    while i < len(replaced):
        sym, expr = replaced[i]
        args = expr.args
        # Check if the expression is a multiplication of two terms, where one of them is a negative symbol
        # If found, substitute this expression back into all further replaced expressions and reduced expressions
        if (
            expr.func == sp.Mul
            and len(expr.args) == 2
            and any(
                a1.func == sp.Symbol
                and (a2 == sp.S.NegativeOne or "_NegativeOne_" in str(a2))
                for a1, a2 in [args, reversed(args)]
            )
        ):
            for k in range(i + 1, len(replaced)):
                if sym in replaced[k][1].free_symbols:
                    replaced[k] = (replaced[k][0], replaced[k][1].subs(sym, expr))
            for k, element in enumerate(reduced):
                if sym in element.free_symbols:
                    reduced[k] = reduced[k].subs(sym, expr)
            # Remove the current expression from the replaced list as its substitutions are done
            replaced.pop(i)
            if i != 0:
                i -= 1
        # Check if the expression is an addition or product of 2 or less symbols, or a square of a symbol
        # If found, count the number of occurrences of the substituted symbol in all further replaced expressions and reduced expressions
        # If the count is 2 or less, substitute this expression back into all further replaced expressions and reduced expressions
        if (
            (expr.func in {sp.Add, sp.Mul})
            and 0 < len(expr.args) <= 2
            and all(
                (arg.func == sp.Symbol or arg.is_Integer or arg.is_Rational)
                for arg in expr.args
            )
        ) or (
            expr.func == sp.Pow and expr.args[0].func == sp.Symbol and expr.args[1] == 2
        ):
            sym_count = 0  # Count the number of occurrences of the substituted symbol
            for k in range(len(replaced) - i):
                # Check if the substituted symbol appears in the replaced expressions
                if sym in replaced[i + k][1].free_symbols:
                    for arg in sp.preorder_traversal(replaced[i + k][1]):
                        if arg.func == sp.Symbol and str(arg) == str(sym):
                            sym_count += 1
            for element in reduced:
                # Check if the substituted symbol appears in the reduced expression
                if sym_count <= 2 and sym in element.free_symbols:
                    str_sym = str(sym)
                    for arg in sp.preorder_traversal(element):
                        if arg.func == sp.Symbol and str(arg) == str_sym:
                            sym_count += 1
            # If the number of occurrences of the substituted symbol is 2 or less, back-substitute
            if 0 < sym_count <= 2:
                for k in range(i + 1, len(replaced)):
                    if sym in replaced[k][1].free_symbols:
                        replaced[k] = (replaced[k][0], replaced[k][1].subs(sym, expr))
                for k, element in enumerate(reduced):
                    if sym in element.free_symbols:
                        reduced[k] = reduced[k].xreplace({sym: expr})
                # Remove the current expression from the replaced list as its substitutions are done
                replaced.pop(i)
                i -= 1
        i += 1
    # Return the final processed replaced and reduced expressions
    return replaced, reduced


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

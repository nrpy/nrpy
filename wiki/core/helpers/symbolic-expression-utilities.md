# Symbolic Expression Utilities

> Reusable helpers for inspecting, transforming, and composing symbolic expressions. · Status: confirmed
> Up: [Helper APIs](index.md)

## Summary

These helpers inspect and transform SymPy expressions without defining the equations or public code-generation interfaces. `ExprTree` supplies mutable expression nodes. `expression_utils` collects free-symbol names and emits parameter aliases. `f2r` constructs rationals from decimal strings, including scientific notation. `functional.py` supplies composition and iterable utilities.

## Detail

`ExprTree` wraps a SymPy expression in `ExprTree.Node` objects. Each node stores its expression and mutable children. `build()` puts the parent expression function in each child's `.func` field; the root field is `None`. `reconstruct()` instead calls each non-leaf node's current `expr.func` on its child expressions, optionally enabling SymPy evaluation; a function that rejects the `evaluate` keyword with `TypeError`, such as `ExprCondPair` (the pair holding each `Piecewise` branch), `Integral`, `Sum`, or `Indexed`, is called again without it. `Piecewise` itself accepts the keyword. Changing the stored `.func` field therefore does not change reconstruction. `preorder()` and `postorder()` traverse the nodes.

Evidence (descriptive behavior): `expr_tree.py::ExprTree.build`; `expr_tree.py::ExprTree.reconstruct`.

`_get_free_symbol_names_cached()` traverses expression arguments with an object-id memo and unions their symbol names. A node with `bound_symbols`, such as `Integral`, `Lambda`, or `Subs`, and a `Derivative`, `Limit`, or `MatrixSymbol` take their names from SymPy's `free_symbols` for that node, less the excluded names, so the result keeps `x` for `Integral(x, (x, 0, x))`, where `x` is free in the upper limit. `get_unique_expression_symbols_as_strings()` sorts the resulting names after exclusions. Parameter grouping through `get_params_commondata_symbols_from_expr_list()` uses the same names.

Evidence (descriptive behavior): `expression_utils.py::_get_free_symbol_names_cached`.

`get_params_commondata_symbols_from_expr_list()` uses the same memoized traversal across a list of expressions, then looks up every sorted symbol name in `par.glb_code_params_dict`. Registered `CodeParameter` symbols with `commondata=True` go into the commondata result list; other registered code parameters go into the params result list; unregistered free symbols are ignored by this bucketing helper.

`generate_definition_header()` emits one definition line per requested symbol by default, using `const REAL <name> = <var_access><name>;`. When `enable_intrinsics=True`, it emits a scalar `NOSIMD<name>` alias and a `REAL_SIMD_ARRAY` alias through `ConstSIMD(...)`. The helper only builds the definition text; the calling codegen or infrastructure module places that text into generated functions.

`f2r()` passes the string form of its input to `sympy.Rational`, which keeps the decimal value for fixed-point and scientific notation alike: `1e-7` gives `1/10000000` and `1.2e-7` gives `3/25000000`. It does not recover precision already lost in a float.

Evidence (descriptive behavior): `float_to_rational.py::f2r`.

`functional.py` provides composition and iterable helpers. `pipe()` applies functions in sequence; `repeat()` applies one function repeatedly; `chain()` and `flatten()` flatten iterable content; `reduce()` combines values; and `uniquify()` preserves first occurrences.

## Sources

- [nrpy/helpers/expr_tree.py](../../../nrpy/helpers/expr_tree.py) - `ExprTree`, `ExprTree.Node`, `build`, `preorder`, `postorder`, `reconstruct`
- [nrpy/helpers/expression_utils.py](../../../nrpy/helpers/expression_utils.py) - `_get_free_symbol_names_cached`, `get_unique_expression_symbols_as_strings`, `get_params_commondata_symbols_from_expr_list`, `generate_definition_header`
- [nrpy/helpers/float_to_rational.py](../../../nrpy/helpers/float_to_rational.py) - `f2r`
- [nrpy/helpers/functional.py](../../../nrpy/helpers/functional.py) - `pipe`, `repeat`, `chain`, `flatten`, `reduce`, `uniquify`

## See Also

- See also: [CSE And Printer Support](cse-and-printer-support.md)
- See also: [Maintenance And Validation Helpers](maintenance-and-validation-helpers.md)
- See also: [C Codegen](../c-codegen.md)

# Symbolic Expression Utilities

> Reusable helpers for inspecting, transforming, and composing symbolic expressions. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Helper APIs](index.md)

## Summary

These helpers support symbolic plumbing around SymPy expressions without owning equation semantics or public code generation. `ExprTree` gives callers a mutable n-ary view of a SymPy expression tree. `expression_utils` extracts free-symbol names, separates registered code parameters from commondata symbols, and emits C definition headers for parameter aliases. `f2r` converts float-like inputs to high-precision `sympy.Rational` values, and `functional.py` provides small composition and iterable utilities.

## Detail

`ExprTree` wraps a SymPy expression in `ExprTree.Node` objects. Each node stores the current expression, the expression function used to rebuild it, and a mutable list of children. `build()` expands children from `expr.args`; `preorder()` and `postorder()` walk the tree; and `reconstruct()` rebuilds parent expressions from child expressions, optionally asking SymPy to evaluate during reconstruction. This is a helper for local tree manipulation, not a replacement for all SymPy expression APIs.

`_get_free_symbol_names_cached()` traverses a SymPy expression graph with an object-id memo and returns frozen sets of symbol names. It preserves SymPy free-symbol behavior for bound-symbol constructs by removing names listed in a node's `bound_symbols`. `get_unique_expression_symbols_as_strings()` exposes the result as a sorted list after applying an optional name exclusion set.

`get_params_commondata_symbols_from_expr_list()` uses the same memoized traversal across a list of expressions, then looks up every sorted symbol name in `par.glb_code_params_dict`. Registered `CodeParameter` symbols with `commondata=True` go into the commondata result list; other registered code parameters go into the params result list; unregistered free symbols are ignored by this bucketing helper.

`generate_definition_header()` emits one definition line per requested symbol by default, using `const REAL <name> = <var_access><name>;`. When `enable_intrinsics=True`, it emits a scalar `NOSIMD<name>` alias and a `REAL_SIMD_ARRAY` alias through `ConstSIMD(...)`. The helper only builds the definition text; the owning codegen or infrastructure path decides where that text is used.

`f2r()` turns a float-like object into a string, ensures the string has a decimal point, appends `zpad` zeros, and passes the padded string to `sympy.Rational`. The intent is to reduce decimal-to-binary floating-point noise before symbolic solves or generated expressions consume constants.

`functional.py` is a compact pure-Python helper group. `pipe()` applies functions in sequence; `repeat()` applies one function a fixed number of times; `chain()` and `flatten()` yield flattened iterable content; `reduce()` applies a binary function cumulatively; `uniquify()` preserves first occurrence order; and `product()` yields Cartesian products, including a `repeat` keyword for repeating one input iterable.

## Sources

- [nrpy/helpers/expr_tree.py](../../../nrpy/helpers/expr_tree.py) - `ExprTree`, `ExprTree.Node`, `build`, `preorder`, `postorder`, `reconstruct`
- [nrpy/helpers/expression_utils.py](../../../nrpy/helpers/expression_utils.py) - `_get_free_symbol_names_cached`, `get_unique_expression_symbols_as_strings`, `get_params_commondata_symbols_from_expr_list`, `generate_definition_header`
- [nrpy/helpers/float_to_rational.py](../../../nrpy/helpers/float_to_rational.py) - `f2r`
- [nrpy/helpers/functional.py](../../../nrpy/helpers/functional.py) - `pipe`, `repeat`, `chain`, `flatten`, `reduce`, `uniquify`, `product`

## See Also

- [Helper APIs](index.md)
- [CSE And Printer Support](cse-and-printer-support.md)
- [Maintenance And Validation Helpers](maintenance-and-validation-helpers.md)
- [C Codegen](../c-codegen.md)
- [Gridfunctions And Parameters](../gridfunctions-and-parameters.md)

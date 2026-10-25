# Symbolic Expression Utilities

> Reusable helpers for inspecting, transforming, and composing symbolic expressions. · Status: confirmed
> Up: [Helper APIs](index.md)

## Summary

These helpers inspect and transform SymPy expressions without defining the equations or public code-generation interfaces. `ExprTree` supplies mutable expression nodes. `expression_utils` collects free-symbol names and emits parameter aliases. `f2r` constructs rationals from decimal strings, including scientific notation. `functional.py` supplies composition and iterable utilities, with `product` equal to `itertools.product`.

## Detail

`ExprTree` wraps a SymPy expression in `ExprTree.Node` objects. Each node stores its expression and mutable children. `build()` puts the parent expression function in each child's `.func` field; the root field is `None`. `reconstruct()` instead calls each non-leaf node's current `expr.func` on its child expressions, optionally enabling SymPy evaluation; a function that rejects the `evaluate` keyword with `TypeError`, such as `Piecewise`, `Integral`, `Sum`, or `Indexed`, is called again without it. Changing the stored `.func` field therefore does not change reconstruction. `preorder()` and `postorder()` traverse the nodes.

Claim evidence:
- Claim: Reconstruction uses `node.expr.func`, not the `.func` field assigned by `build()`, and retries a function without the `evaluate` keyword when that keyword raises `TypeError`.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/expr_tree.py](../../../nrpy/helpers/expr_tree.py), `ExprTree.build, ExprTree.reconstruct`
- Corroboration: `none available`; the named implementation decides this behavior.

`_get_free_symbol_names_cached()` traverses expression arguments with an object-id memo and unions their symbol names. A node with `bound_symbols`, such as `Integral`, `Lambda`, or `Subs`, and a `Derivative` take their names from SymPy's `free_symbols` for that node, less the excluded names, so the result keeps `x` for `Integral(x, (x, 0, x))`, where `x` is free in the upper limit. `get_unique_expression_symbols_as_strings()` sorts the resulting names after exclusions. Parameter grouping through `get_params_commondata_symbols_from_expr_list()` uses the same names.

Claim evidence:
- Claim: `_get_free_symbol_names_cached()` takes the names of a node with `bound_symbols`, or of a `Derivative`, from SymPy's `free_symbols` for that node and every other node's names from the union over its arguments, less the excluded names, so `Integral(x, (x, 0, x))` keeps `x`.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/expression_utils.py](../../../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached`
- Corroboration: `none available`; the named implementation decides this behavior.

`get_params_commondata_symbols_from_expr_list()` uses the same memoized traversal across a list of expressions, then looks up every sorted symbol name in `par.glb_code_params_dict`. Registered `CodeParameter` symbols with `commondata=True` go into the commondata result list; other registered code parameters go into the params result list; unregistered free symbols are ignored by this bucketing helper.

`generate_definition_header()` emits one definition line per requested symbol by default, using `const REAL <name> = <var_access><name>;`. When `enable_intrinsics=True`, it emits a scalar `NOSIMD<name>` alias and a `REAL_SIMD_ARRAY` alias through `ConstSIMD(...)`. The helper only builds the definition text; the owning codegen or infrastructure path decides where that text is used.

`f2r()` converts its input to a lower-case string, splits off any exponent at `e`, adds a decimal point to the mantissa if none is present, appends `zpad` zeros to the mantissa only, and passes the mantissa, zeros, and exponent to `sympy.Rational`. The decimal value is preserved for fixed-point and scientific notation alike: with `zpad=1`, `1e-7` gives `1/10000000` and `1.2e-7` gives `3/25000000`. Padding does not recover precision already lost in a float.

Claim evidence:
- Claim: `f2r()` appends the `zpad` zeros to the mantissa before any exponent, so the rational equals the decimal value of the input; padding adds no information.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/float_to_rational.py](../../../nrpy/helpers/float_to_rational.py), `f2r`
- Corroboration: `none available`; the named implementation decides this behavior.

`functional.py` provides composition and iterable helpers. `pipe()` applies functions in sequence; `repeat()` applies one function repeatedly; `chain()` and `flatten()` flatten iterable content; `reduce()` combines values; and `uniquify()` preserves first occurrences. `product(*iterable, repeat=1)` yields the tuples of `itertools.product(*iterable, repeat=repeat)`: iterable-valued elements stay whole, one-shot iterators are consumed once, and `repeat=0` yields a single empty tuple.

Claim evidence:
- Claim: `product()` delegates to `itertools.product(*iterable, repeat=repeat)` and yields its tuples unchanged.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/functional.py](../../../nrpy/helpers/functional.py), `product`
- Corroboration: `none available`; the named implementation decides this behavior.

## Sources

- [nrpy/helpers/expr_tree.py](../../../nrpy/helpers/expr_tree.py) - `ExprTree`, `ExprTree.Node`, `build`, `preorder`, `postorder`, `reconstruct`
- [nrpy/helpers/expression_utils.py](../../../nrpy/helpers/expression_utils.py) - `_get_free_symbol_names_cached`, `get_unique_expression_symbols_as_strings`, `get_params_commondata_symbols_from_expr_list`, `generate_definition_header`
- [nrpy/helpers/float_to_rational.py](../../../nrpy/helpers/float_to_rational.py) - `f2r`
- [nrpy/helpers/functional.py](../../../nrpy/helpers/functional.py) - `pipe`, `repeat`, `chain`, `flatten`, `reduce`, `uniquify`, `product`

## See Also

- Parent: [Helper APIs](index.md)
- See also: [CSE And Printer Support](cse-and-printer-support.md)
- See also: [Maintenance And Validation Helpers](maintenance-and-validation-helpers.md)
- See also: [C Codegen](../c-codegen.md)
- See also: [Gridfunctions And Parameters](../gridfunctions-and-parameters.md)

# Symbolic Expression Utilities

> Reusable helpers for inspecting, transforming, and composing symbolic expressions. · Status: contested
> Up: [Helper APIs](index.md)

## Summary

These helpers inspect and transform SymPy expressions without defining the equations or public code-generation interfaces. `ExprTree` supplies mutable expression nodes. `expression_utils` collects symbol names and emits parameter aliases, with a bound-symbol limitation described below. `f2r` constructs rationals from fixed-point decimal strings; scientific notation is unsafe. `functional.py` supplies composition and iterable utilities with restricted product behavior.

## Detail

`ExprTree` wraps a SymPy expression in `ExprTree.Node` objects. Each node stores its expression and mutable children. `build()` puts the parent expression function in each child's `.func` field; the root field is `None`. `reconstruct()` instead calls each non-leaf node's current `expr.func` on its child expressions, optionally enabling SymPy evaluation. Changing the stored `.func` field therefore does not change reconstruction. `preorder()` and `postorder()` traverse the nodes.

Claim evidence:
- Claim: Reconstruction uses `node.expr.func`, not the `.func` field assigned by `build()`.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/expr_tree.py](../../../nrpy/helpers/expr_tree.py), `ExprTree.build, ExprTree.reconstruct`
- Corroboration: `none available`; the named implementation decides this behavior.

`_get_free_symbol_names_cached()` traverses expression arguments with an object-id memo, unions their symbol names, then removes names in each node's `bound_symbols`. This does not fully preserve SymPy's free-symbol semantics: `Integral(x, (x, 0, x))` has free symbol `x` in its upper limit, but this helper removes it. `get_unique_expression_symbols_as_strings()` sorts the resulting names after exclusions. Parameter grouping through `get_params_commondata_symbols_from_expr_list()` inherits this omission.

Claim status: contested; contradiction: CONTR-0016. See [CONTR-0016](../../contradictions.md#contr-0016).

Claim evidence:
- Claim: Bound-symbol subtraction can remove free occurrences in integration limits or substitution points; full SymPy free-symbol equivalence is not guaranteed.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/expression_utils.py](../../../nrpy/helpers/expression_utils.py), `_get_free_symbol_names_cached`
- Corroboration: `none available`; the named implementation decides this behavior.

`get_params_commondata_symbols_from_expr_list()` uses the same memoized traversal across a list of expressions, then looks up every sorted symbol name in `par.glb_code_params_dict`. Registered `CodeParameter` symbols with `commondata=True` go into the commondata result list; other registered code parameters go into the params result list; unregistered free symbols are ignored by this bucketing helper.

`generate_definition_header()` emits one definition line per requested symbol by default, using `const REAL <name> = <var_access><name>;`. When `enable_intrinsics=True`, it emits a scalar `NOSIMD<name>` alias and a `REAL_SIMD_ARRAY` alias through `ConstSIMD(...)`. The helper only builds the definition text; the owning codegen or infrastructure path decides where that text is used.

`f2r()` converts its input to a string, adds a decimal point if none is present, appends `zpad` zeros, and passes the result to `sympy.Rational`. Fixed-point decimal strings preserve the decimal value; padding does not recover precision already lost in a float. Scientific notation is unsafe: with `zpad=1`, `1e-7` produces an invalid string, while `1.2e-7` becomes the value of `1.2e-70`. Supply a fixed-point string for the intended rational.

Claim status: contested; contradiction: CONTR-0017. See [CONTR-0017](../../contradictions.md#contr-0017).

Claim evidence:
- Claim: `f2r()` appends zeros to the whole input string; scientific notation may fail or change the exponent, and padding adds no information.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/float_to_rational.py](../../../nrpy/helpers/float_to_rational.py), `f2r`
- Corroboration: `none available`; the named implementation decides this behavior.

`functional.py` provides composition and iterable helpers. `pipe()` applies functions in sequence; `repeat()` applies one function repeatedly; `chain()` and `flatten()` flatten iterable content; `reduce()` combines values; and `uniquify()` preserves first occurrences. `product()` combines inputs through repeated reduction and flattening. It is not a general replacement for `itertools.product`: a single integer-valued input raises `TypeError`, a one-shot second iterator loses later combinations, and iterable-valued elements are flattened. Its `repeat` handling repeats one input only when `repeat > 1`; `repeat=0` does not produce the conventional empty product. Use reusable scalar-valued inputs in the supported multi-input form, or verify the required edge case separately.

Claim status: contested; contradiction: CONTR-0018. See [CONTR-0018](../../contradictions.md#contr-0018).

Claim evidence:
- Claim: `product()` reduces and flattens its inputs and supports only a restricted subset of Cartesian-product behavior.
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

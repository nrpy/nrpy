# CSE And Printer Support

> Helper route for CSE preprocessing, CSE postprocessing, register-pressure statement ordering, custom C printer mappings, and JAX printer behavior. · Status: contested
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.cse_preprocess_postprocess` wraps SymPy common-subexpression elimination with NRPy-specific preprocessing, postprocessing, and deterministic sorting helpers. `nrpy.helpers.register_pressure_ordering` reorders the statements of an already-emitted finite-difference kernel body to shorten live ranges. `nrpy.helpers.custom_c_codegen_functions` supplies the user-function table passed to SymPy's C printer from C codegen. `nrpy.helpers.jax_printer` supplies the `NRPyJaxPrinter` used by JAX-oriented Python codegen. These helpers own expression and printer mechanics; `c_codegen()` and `py_codegen()` own the public code-generation interfaces that call them.

## Detail

`cse_preprocess()` accepts either one SymPy expression or a list of expressions and returns a list plus a symbol-to-rational dictionary. It walks each expression through `ExprTree`, replaces non-power integer and rational constants with generated symbols, can declare `_NegativeOne_`, can factor on the generated rational symbols, and can optionally factor negative-one symbols. Before any replacement it raises `ValueError` when an input symbol is spelled like a constant it generates: `<prefix>_Integer_N`, `<prefix>_Rational_N_M`, or `<prefix>_NegativeOne_`. The factoring passes over function arguments skip a function that has no arguments or whose first argument is not a SymPy expression. With `declare_neg1_as_symbol`, a power of `<prefix>_NegativeOne_` folds to `<prefix>_NegativeOne_` or `<prefix>_Integer_1`, by parity, only when its exponent is an integer; an even power registers `<prefix>_Integer_1` among the generated constants. With `debug=True`, it builds a separate `debug_tree` from the processed expression, replaces the generated symbols by their rational values, reconstructs that tree, and raises `ValueError` when `sp.simplify()` of its difference from the input expression is nonzero.

`cse_postprocess()` consumes the `(replacements, reduced_exprs)` shape produced by SymPy CSE. It moves `sympy.Eq` scalar temporaries from the reduced-expression list into the replacement list, substitutes constant replacements into later expressions, orders replacements so dependencies are available before use, back-substitutes small or cheap temporaries when their use count stays low, handles products involving `_NegativeOne_`-style symbols, and removes replacement entries that are no longer used by the final reduced expressions.

`sort_cse_output_deterministically()` sorts dependency-ready CSE replacements and renames them to `symbol_prefix + "tmpN"`, numbering them in order and skipping every name already in use: a name in its `reserved_names` argument, or the name of a free symbol of the replacement or reduced expressions that is not itself a replaced temporary. `c_codegen()` passes the identifiers of its output strings as `reserved_names`, so a temporary never takes the name of an output variable or an input symbol.

Claim evidence:
- Claim: `cse_preprocess()` raises `ValueError` for an input symbol named `<prefix>_Integer_N`, `<prefix>_Rational_N_M`, or `<prefix>_NegativeOne_`; `sort_cse_output_deterministically()` gives no temporary a name in `reserved_names` or the name of a free symbol that is not a replaced temporary; `c_codegen()` passes the output identifiers as `reserved_names`.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/cse_preprocess_postprocess.py](../../../nrpy/helpers/cse_preprocess_postprocess.py), `cse_preprocess, sort_cse_output_deterministically`; [nrpy/c_codegen.py](../../../nrpy/c_codegen.py), `c_codegen`
- Corroboration: `none available`; the named implementations decide this behavior.

The structural key of `sort_cse_output_deterministically()` sorts arguments whenever the resulting expression is commutative, including powers: `x**y` and `y**x` can have equal keys. Ties retain input order, so replacement-list permutations need not produce identical emitted text.

Claim status: contested; contradiction: CONTR-0014. See [CONTR-0014](../../contradictions.md#contr-0014).

Claim evidence:
- Claim: CSE ordering is not invariant under tied-key input permutations, because the structural key sorts the arguments of commutative-valued expressions such as powers and ties retain input order.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/cse_preprocess_postprocess.py](../../../nrpy/helpers/cse_preprocess_postprocess.py), `sort_cse_output_deterministically`
- Corroboration: `none available`; the named implementation decides this behavior.

`order_statements_for_register_pressure()` takes the C text of a finite-difference kernel body as `c_codegen()` assembles it (every gridfunction read, then every derivative, then every common subexpression, then every output), splits it into statements, and re-emits the same statements in a greedy free-most-first list-schedule order: compile-time constants, and `ConstSIMD`/`ConstCUDA` broadcasts of them, first in their original order, then ready loads and definitions with the statement whose operands die chosen next, then stores last in their original order for alias safety. Statement text is unchanged, so source arithmetic operations are unchanged; reordering can change the compiler's choice of which multiply/add pairs to contract and alter final rounding, including for nonzero results. The helper is conservative: a body containing braces, preprocessor lines, repeated assignment to one name, a `const` definition that follows a store, or any statement that is neither a single-assignment `const` (or `constexpr`) definition nor a store is returned unchanged. `c_codegen()` calls it after `gridfunction_management_and_FD_codegen()` on both the finite-difference and automatic-memory-read routes; there is no separate scheduling option to turn it off.

Claim evidence:
- Claim: `order_statements_for_register_pressure()` schedules constants first, ready loads and definitions next, and original-order stores last for alias safety; source statements remain unchanged, but changed compiler FMA contraction can alter final rounding for nonzero results; `c_codegen()` applies the helper after the FD/memory-read helper on both routes without a separate scheduling opt-out.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/register_pressure_ordering.py](../../../nrpy/helpers/register_pressure_ordering.py), `order_statements_for_register_pressure`; [nrpy/c_codegen.py](../../../nrpy/c_codegen.py), `c_codegen`
- Corroboration: none available; the owning implementations and targeted local validation decide the claim.

`custom_functions_for_SymPy_ccode` maps selected SymPy printer function names by floating-point mode. For `double`, `float`, and `long double`, it maps `nrpyAbs` and `fabs` to the matching absolute-value functions, using the integer `abs()` only for an integer literal of magnitude below `2**31` and `fabs`, `fabsf`, or `fabsl` otherwise, and rewrites common powers to `sqrt`, `cbrt`, repeated multiplication, or reciprocal forms before falling back to `pow`, `powf`, or `powl`. A reciprocal form prints `1.0/(x*...*x)`, so for a base of magnitude above the `1/n` power of the largest finite value (about `1.3e154` for `n = 2` in `double`, and about `1.8e19` in `float`) the product overflows and the quotient returns 0, although the exact result is representable. For `double complex`, it maps absolute value (`nrpyAbs`, `fabs`, and `Abs` to `cabs`) and the transcendental functions `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `exp`, `log`, `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, and `atanh` to C99 complex forms and maps powers to `csqrt`, `cpow`, repeated multiplication, or reciprocal forms. `c_codegen()` passes the table to `sp.ccode()` for ordinary assignments and CSE temporaries.

`NRPyJaxPrinter` subclasses SymPy's JAX printer when available and falls back to the NumPy printer on older SymPy installs. It rewrites known function and constant mappings to use the `jnp` namespace, sets `_module = "jnp"`, and mirrors NRPy's C power simplifications for JAX output: square roots and cube roots become `jnp.sqrt()` and `jnp.cbrt()`, small positive integer powers become repeated multiplication, and small negative integer powers become reciprocals of repeated multiplication. Unsupported powers fall back to the parent printer.

The real C mappings and JAX printer rewrite a one-third power to a real cube root without checking the base's domain. For a negative real base, this differs from SymPy's principal complex power: at `-8`, the real cube root is `-2`, while the principal one-third power is `1 + sqrt(3)*I`. Use these rewrites only when the intended equation calls for the real root. The C `double complex` mapping uses complex power functions instead.

Claim evidence:
- Claim: Real cube-root printing does not preserve principal-complex-power semantics on negative real bases.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/custom_c_codegen_functions.py](../../../nrpy/helpers/custom_c_codegen_functions.py), `custom_functions_for_SymPy_ccode`
- Corroboration: `none available`; the named implementation decides this behavior.
- Additional deciding authority: [nrpy/helpers/jax_printer.py](../../../nrpy/helpers/jax_printer.py), `NRPyJaxPrinter._print_Pow`.

Two further rules keep the printed expressions evaluable by JAX. An integer that does not fit in a signed 32-bit integer is printed as a float, and a rational whose numerator or denominator does not fit in a signed 32-bit integer is printed as the float value of the rational, with `ValueError` raised when that value exceeds the range of a double-precision float; smaller integers and rationals keep their exact form, such as `3` or `1/3`. In operations with arrays, JAX rejects Python integers that do not fit in a signed 32-bit integer (signed 64-bit with `jax_enable_x64`). `Max` and `Min` are printed as nested two-argument `jnp.maximum` and `jnp.minimum` calls, so the printed code needs no import of `functools`, which SymPy 1.14 and later use for these functions.

Claim evidence:
- Claim: `NRPyJaxPrinter` prints integers of magnitude at least `2**31`, and rationals whose numerator magnitude or denominator is at least `2**31`, as Python float literals, and prints `Max`/`Min` as nested `jnp.maximum`/`jnp.minimum` calls. Float printing rounds such values to double precision, and a rational beyond the double range raises `ValueError`.
- Role: descriptive behavior
- Deciding authority: [nrpy/helpers/jax_printer.py](../../../nrpy/helpers/jax_printer.py), `NRPyJaxPrinter._print_Integer`, `_print_Rational`, `_print_Max`, `_print_Min`
- Corroboration: `none available`; the only targeted check is the `NRPyJaxPrinter` class doctest in the deciding module, which validates this behavior locally but is not a separate source.

The JAX printer also implements `_print_ArrayElementwiseApplyFunc()`. Unary lambda elementwise application is lowered by printing the array operand once, printing the scalar lambda body with a sentinel symbol, and replacing each whole-identifier occurrence of that sentinel with the parenthesized array expression. The sentinel name gains trailing underscores until it differs from every free symbol of the body, and a longer name that contains the sentinel text is left unchanged. Multi-argument lambdas fall back to `jnp.vectorize(...)`, and non-lambda callables are printed as callable applications to the array string.

`py_codegen()` is the local integration point for `NRPyJaxPrinter`: it constructs a module-level printer, requires the `Infrastructure` parameter to be `JAX`, optionally runs SymPy CSE through `cse_postprocess()`, prints each unexpanded right-hand side through `printer.doprint()`, and writes the assignment itself; see [Python Codegen](../python-codegen.md) for why the output name is not passed to the printer. Unlike the C path, this Python path does not call `cse_preprocess()` or the deterministic `order="none"` post-sort helper in the current source.

## Sources

- [nrpy/helpers/cse_preprocess_postprocess.py](../../../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_preprocess`, `cse_postprocess`, `sort_cse_output_deterministically`
- [nrpy/helpers/register_pressure_ordering.py](../../../nrpy/helpers/register_pressure_ordering.py) - `order_statements_for_register_pressure`
- [nrpy/helpers/custom_c_codegen_functions.py](../../../nrpy/helpers/custom_c_codegen_functions.py) - `custom_functions_for_SymPy_ccode`
- [nrpy/helpers/jax_printer.py](../../../nrpy/helpers/jax_printer.py) - `NRPyJaxPrinter`, `_print_Pow`, `_print_Integer`, `_print_Rational`, `_print_Max`, `_print_Min`, `_print_ArrayElementwiseApplyFunc`
- [nrpy/c_codegen.py](../../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`
- [nrpy/py_codegen.py](../../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`

## See Also

- Parent: [Helper APIs](index.md)
- See also: [C Codegen](../c-codegen.md)
- See also: [Python Codegen](../python-codegen.md)
- See also: [Finite Difference](../finite-difference.md)
- See also: [SIMD And Intrinsic Support](simd-and-intrinsic-support.md)
- See also: [Symbolic Expression Utilities](symbolic-expression-utilities.md)

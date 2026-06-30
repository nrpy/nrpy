# CSE And Printer Support

> Helper route for CSE preprocessing, CSE postprocessing, custom C printer mappings, and JAX printer behavior. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [Helper APIs](index.md)

## Summary

`nrpy.helpers.cse_preprocess_postprocess` wraps SymPy common-subexpression elimination with NRPy-specific preprocessing, postprocessing, and deterministic sorting helpers. `nrpy.helpers.custom_c_codegen_functions` supplies the user-function table passed to SymPy's C printer from C codegen. `nrpy.helpers.jax_printer` supplies the `NRPyJaxPrinter` used by JAX-oriented Python codegen. These helpers own expression and printer mechanics; `c_codegen()` and `py_codegen()` own the public code-generation interfaces that call them.

## Detail

`cse_preprocess()` accepts either one SymPy expression or a list of expressions and returns a list plus a symbol-to-rational dictionary. It walks each expression through `ExprTree`, replaces non-power integer and rational constants with generated symbols, can declare `_NegativeOne_`, can factor on the generated rational symbols, and can optionally factor negative-one symbols. The current debug branch is intended to back-substitute generated symbols with rational values and compare the result, but the living source mutates a separate `debug_tree` and then reconstructs `tree`; this page therefore does not treat debug mode as effective back-substitution validation until the source behavior is reconciled.

`cse_postprocess()` consumes the `(replacements, reduced_exprs)` shape produced by SymPy CSE. It moves `sympy.Eq` scalar temporaries from the reduced-expression list into the replacement list, substitutes constant replacements into later expressions, orders replacements so dependencies are available before use, back-substitutes small or cheap temporaries when their use count stays low, handles products involving `_NegativeOne_`-style symbols, and removes replacement entries that are no longer used by the final reduced expressions.

`sort_cse_output_deterministically()` is a post-sort and rename pass for fast `sp.cse(..., order="none")` call sites. It builds replacement dependencies, repeatedly selects dependency-ready temporaries, sorts each ready set using deterministic expression keys, and renames the chosen symbols to a stable `symbol_prefix + "tmpN"` sequence. This is NRPy's deterministic helper behavior around an already-computed CSE result; it is not a guarantee about upstream SymPy's unsorted CSE order.

`custom_functions_for_SymPy_ccode` maps selected SymPy printer function names by floating-point mode. For `double`, `float`, and `long double`, it maps `nrpyAbs` and `fabs` to the matching absolute-value functions and rewrites common powers to `sqrt`, `cbrt`, repeated multiplication, or reciprocal forms before falling back to `pow`, `powf`, or `powl`. For `double complex`, it maps absolute value and common transcendental functions to C99 complex forms and maps powers to `csqrt`, `cpow`, repeated multiplication, or reciprocal forms. `c_codegen()` passes the table to `sp.ccode()` for ordinary assignments and CSE temporaries.

`NRPyJaxPrinter` subclasses SymPy's JAX printer when available and falls back to the NumPy printer on older SymPy installs. It rewrites known function and constant mappings to use the `jnp` namespace, sets `_module = "jnp"`, and mirrors NRPy's C power simplifications for JAX output: square roots and cube roots become `jnp.sqrt()` and `jnp.cbrt()`, small positive integer powers become repeated multiplication, and small negative integer powers become reciprocals of repeated multiplication. Unsupported powers fall back to the parent printer.

The JAX printer also implements `_print_ArrayElementwiseApplyFunc()`. Unary lambda elementwise application is lowered by printing the array operand once, printing the scalar lambda body with a sentinel symbol, and replacing that sentinel with the parenthesized array expression. Multi-argument lambdas fall back to `jnp.vectorize(...)`, and non-lambda callables are printed as callable applications to the array string.

`py_codegen()` is the local integration point for `NRPyJaxPrinter`: it constructs a module-level printer, requires the `Infrastructure` parameter to be `JAX`, optionally runs SymPy CSE through `cse_postprocess()`, and emits assignments through `printer.doprint()`. Unlike the C path, this Python path does not call `cse_preprocess()` or the deterministic `order="none"` post-sort helper in the current source.

## Sources

- [nrpy/helpers/cse_preprocess_postprocess.py](../../../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_preprocess`, `cse_postprocess`, `sort_cse_output_deterministically`
- [nrpy/helpers/custom_c_codegen_functions.py](../../../nrpy/helpers/custom_c_codegen_functions.py) - `custom_functions_for_SymPy_ccode`
- [nrpy/helpers/jax_printer.py](../../../nrpy/helpers/jax_printer.py) - `NRPyJaxPrinter`, `_print_Pow`, `_print_ArrayElementwiseApplyFunc`
- [nrpy/c_codegen.py](../../../nrpy/c_codegen.py) - `CCodeGen`, `c_codegen`
- [nrpy/py_codegen.py](../../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`

## See Also

- [Helper APIs](index.md)
- [C Codegen](../c-codegen.md)
- [Python Codegen](../python-codegen.md)
- [Finite Difference](../finite-difference.md)
- [SIMD And Intrinsic Support](simd-and-intrinsic-support.md)
- [Symbolic Expression Utilities](symbolic-expression-utilities.md)

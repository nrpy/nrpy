# Python Codegen

> Core route for turning SymPy expressions into JAX-compatible Python assignment text. · Status: contested
> Up: [Core APIs](index.md)

## Summary

`py_codegen()` emits JAX-compatible Python assignment strings from SymPy expressions. `PyCodeGen` is the option object behind that public function: a stripped-down Python counterpart to `CCodeGen` for JAX-oriented code, without the finite-difference memory-read path used by C kernels. CSE combined with postprocessing substitutions can raise an exception; the limitation below qualifies this interface. An output name that a later printed assignment reads is rejected with `ValueError`.

## Detail

`PyCodeGen` accepts `prestring`, `poststring`, `verbose`, `enable_cse`, `cse_sorting`, `cse_varprefix`, and `postproc_substitution_dict`. The final returned string is the optional verbose comment block followed by `prestring`, generated assignment text, and `poststring`.

The constructor reads `par.parval_from_str("Infrastructure")` and requires the value to be exactly `"JAX"`. If the value differs, it raises `ValueError("Infrastructure must be 'jax' for py_codegen")`. This guard matches the implementation boundary: the module constructs a module-level `NRPyJaxPrinter`, and this Python codegen path targets JAX-compatible assignment output rather than general infrastructure-specific generated packages.

`py_codegen()` rejects tuple inputs for both the SymPy expression argument and the output-name argument. Non-list scalar inputs are normalized into one-element lists, list inputs are copied so later processing does not mutate caller-owned lists, and the expression list length must match the output-name list length.

When `verbose` is enabled, `py_codegen()` emits a Python comment block before generated code. The block records the original SymPy expression or expressions paired with their output assignment names; plural output uses bracketed `"[name = expression]"` lines.

With `enable_cse=False`, each expression is emitted independently. If `postproc_substitution_dict` is nonempty, `apply_substitution_dict()` first rewrites matching free-symbol names by appending the configured suffix. The right-hand side is then printed with `printer.doprint(expr)`, and the line `output_name = <printed expression>` is appended to the output string.

With `enable_cse=True`, `py_codegen()` collects the expressions and output names, then calls SymPy CSE with numbered temporaries from `cse_varprefix + "tmp"` and the configured `cse_sorting` order. For SymPy versions before 1.3, the implementation prints a warning and uses the raw `sp.cse()` result. For SymPy 1.3 and newer, it passes the `sp.cse()` result through `cse_postprocess()`. Without postprocessing substitutions, CSE temporaries and final reduced expressions are emitted without `sp.expand()`: `printer.doprint()` prints each right-hand side and `py_codegen()` writes the assignment. In the postprocessed CSE path, replacement entries are tuples. A nonempty `postproc_substitution_dict` makes `py_codegen()` attempt to assign into a tuple whenever a temporary survives, raising `TypeError`, even if no substitution key matches. Substitution works on the no-CSE path and on final reduced expressions when no temporary survives.

Claim status: contested; contradiction: CONTR-0012.
See [CONTR-0012](../contradictions.md#contr-0012) for the source discrepancy and
required resolution check. Combining CSE and postprocessing substitutions is
not guaranteed to return code.

Both paths emit output assignments in list order into one Python namespace, so
an output assignment replaces any earlier value bound to that name. Before
returning code, `py_codegen()` parses each printed output right-hand side with
Python's `ast` module and raises `ValueError` if an earlier output name is a
Python name that a later printed right-hand side reads. Because the check uses
the printed code, after `postproc_substitution_dict` renaming and after the
module-level `NRPyJaxPrinter` renames reserved words, it covers an input symbol
such as `x`; a CSE temporary such as `tmp0` for expressions
`[x**3 + sin(x**3), x**3]` with output names `["tmp0", "b"]`; a symbol named
`lambda`, which prints as `lambda_`; and names that the printer emits itself,
such as `jnp` or `abs`. Without this check, the later output would read the
overwritten value and the generated code would return wrong numbers or fail at
run time. An output name that no later right-hand side reads, such as a final
output that reuses an input symbol's name, is accepted and printed unchanged.
This check resolved [CONTR-0013](../contradictions.md#contr-0013).

Claim evidence:
- Claim: `py_codegen()` raises `ValueError` when an earlier output name is a Python name read by a later printed output right-hand side, with CSE enabled or disabled, and otherwise emits the same assignments as without the check. This claims nothing about two outputs sharing one name, or about `prestring` and `poststring`.
- Role: descriptive behavior
- Deciding authority: [nrpy/py_codegen.py](../../nrpy/py_codegen.py), `_check_no_read_after_overwrite` and its two calls in `py_codegen`, and the module-level `printer`
- Corroboration: `none available`; the `py_codegen` doctests in the deciding module show the rejected `tmp0`, input-symbol, and printer-renamed `lambda_` cases and the accepted final-output case, but are not a separate source.

`py_codegen()` never passes the output name to `doprint()` as the `assign_to` argument. In SymPy 1.11 to 1.14, whose printer is `JaxPrinter`, and in some SymPy development commits, that path converts the right-hand side to an array expression and silently drops negative powers, printing, for example, `m2/m1` as `m2` and `1/x**2` as `jnp.einsum("", )`.

Claim evidence:
- Claim: `py_codegen()` prints each temporary and output as `name = <printed right-hand side>`, calling `printer.doprint()` on the right-hand side only; quotients and negative integer powers are therefore printed as divisions or reciprocals on every supported SymPy version. This makes no claim about output formatting beyond the printer's own conventions.
- Role: descriptive behavior
- Deciding authority: [nrpy/py_codegen.py](../../nrpy/py_codegen.py), `py_codegen`
- Corroboration: `none available`; the only targeted check is the `py_codegen` doctest in the deciding module, which prints `m2/m1`, `1/x**2`, and `y/(x + z)` with CSE disabled and enabled but is not a separate source.

Neither path calls `sp.expand()`, matching `c_codegen()`. Expanding the SEOBNRv5
final mass `M_f` places its rational coefficients over common denominators and
prints coefficients larger than the largest finite float32 value. JAX evaluates
such a Python float literal in the float32 precision of a float32 input array,
so the coefficient overflows before the terms cancel, and `M_f` and the QNM
values computed from it become `nan`. The unexpanded expression prints no such
coefficient.

Claim evidence:
- Claim: `py_codegen()` prints CSE temporaries and outputs without calling `sp.expand()`. For the SEOBNRv5 final mass that `register_PyFunction_SEOBNRv5_aligned_spin_coefficients()` passes to `py_codegen()`, the expanded form prints coefficients larger than the largest finite float32 value and the unexpanded form does not. This claims no float32 accuracy bound and no result for other expressions or inputs.
- Role: descriptive behavior
- Deciding authority: [nrpy/py_codegen.py](../../nrpy/py_codegen.py), `py_codegen`; [SEOBNRv5_aligned_spin_constants.py](../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py), `SEOBNR_aligned_spin_constants.final_mass_non_precessing_UIB2016`; [SEOBNRv5_aligned_spin_coefficients.py](../../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py), `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`
- Corroboration: `none available`; no repository test or CI job evaluates the generated coefficient function with a float32 input

Unlike `c_codegen()`, this path does not run `sort_cse_output_deterministically()` when `cse_sorting="none"`. That option requests SymPy's unsorted CSE path; this page makes no deterministic-output guarantee for it.

This page owns the public `py_codegen()` and `PyCodeGen` contract. Helper pages own CSE helper and printer internals; JAX infrastructure pages own generated package assembly.

Import these APIs from `nrpy.py_codegen`. The empty `nrpy/__init__.py` does not re-export them.

## Sources

- [cse_preprocess_postprocess.py](../../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_postprocess` replacement tuples

- [nrpy/py_codegen.py](../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`, `printer`
- [nrpy/c_codegen.py](../../nrpy/c_codegen.py) - `apply_substitution_dict`
- [nrpy/__init__.py](../../nrpy/__init__.py) - empty package initializer; no Python-codegen re-exports

## See Also

- Parent: [Core APIs](index.md)
- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Contrasts with: [C Codegen](c-codegen.md)
- See also: [JAX Project Generation Lifecycle](../infrastructures/jax/project-generation-lifecycle.md)

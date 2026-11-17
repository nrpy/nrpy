# Python Codegen

> SymPy-to-JAX Python assignment code generation via `py_codegen()` and `PyCodeGen`, CSE/substitution limits, and `float32`/overwrite safeguards. · Status: contested
> Up: [Core APIs](index.md)

## Summary

`nrpy.py_codegen.py_codegen()` (configured by `PyCodeGen`) emits JAX-compatible Python assignment blocks from SymPy expressions using `NRPyJaxPrinter`, omitting `c_codegen`'s finite-difference memory-read lowering. It rejects read-after-overwrite output names via `ast` inspection, avoids `assign_to` and `sp.expand()` to preserve negative powers and `float32` coefficient range, and raises `TypeError` when `enable_cse=True` survives with a nonempty `postproc_substitution_dict`.

## Detail

### Options, `Infrastructure == "JAX"` Guard, And CSE Substitution Bug

`PyCodeGen` configures `prestring`, `poststring`, `verbose`, `enable_cse`, `cse_sorting`, `cse_varprefix`, and `postproc_substitution_dict`, returning optional `# [name = expr]` verbose comments + `prestring` + assignments + `poststring`. Its constructor requires `par.parval_from_str("Infrastructure") == "JAX"` (raising `ValueError("Infrastructure must be 'jax' for py_codegen")` otherwise). `py_codegen()` rejects tuple inputs, wraps scalar expressions/names into 1-element lists, copies list arguments without mutating caller lists, and requires equal list lengths (import directly from `nrpy.py_codegen`; `nrpy/__init__.py` is empty).

With `enable_cse=False`, optional `apply_substitution_dict()` appends suffixes to matching free-symbol names before `printer.doprint(expr)` emits `output_name = <rhs>`. With `enable_cse=True`, `sp.cse()` runs with prefix `cse_varprefix + "tmp"` and `cse_sorting` (passing through `cse_postprocess()` on SymPy `>= 1.3`, without `c_codegen`'s `sort_cse_output_deterministically()` when `cse_sorting="none"`). Because `cse_postprocess()` returns tuple entries, a nonempty `postproc_substitution_dict` attempts tuple item assignment whenever any CSE temporary survives, raising `TypeError` even if no key matches ([CONTR-0012](../contradictions.md#contr-0012)).

Claim status: contested; contradiction: CONTR-0012.

### Read-After-Overwrite Check, `doprint` Without `assign_to`, And Unexpanded `float32` Literals

Before returning, `_check_no_read_after_overwrite` parses every printed output RHS with Python's `ast` module (after `postproc_substitution_dict` and `NRPyJaxPrinter` reserved-word renaming such as `lambda -> lambda_`) and raises `ValueError` if an earlier output name is read by a later output RHS—catching overwritten input symbols (`x`), CSE temporaries (`tmp0`), renamed keywords (`lambda_`), or printer identifiers (`jnp`, `abs`), while permitting final outputs whose names are not read afterward.

`py_codegen()` always formats `name = printer.doprint(rhs)` without passing `assign_to` to `doprint()`: in SymPy `1.11`–`1.14`, `JaxPrinter.doprint(..., assign_to=...)` converts expressions to array contractions and silently drops negative powers (printing `m2 / m1` as `m2` and `1 / x**2` as `jnp.einsum("", )`).

Like `c_codegen()`, `py_codegen()` never calls `sp.expand()`: expanding `SEOBNR_aligned_spin_constants.final_mass_non_precessing_UIB2016` (`M_f` in `SEOBNRv5_aligned_spin_coefficients`) combines rational fractions over common denominators exceeding `float32` `max`, overflowing to `nan` on `float32` JAX arrays, whereas the unexpanded expression stays within `float32` range.

Evidence (descriptive behavior): `nrpy/py_codegen.py::_check_no_read_after_overwrite`; `nrpy/py_codegen.py::py_codegen`; `cse_preprocess_postprocess.py::cse_postprocess`; `SEOBNRv5_aligned_spin_constants.py`.

## Sources

- [py_codegen.py](../../nrpy/py_codegen.py) - `PyCodeGen`, `py_codegen`, `_check_no_read_after_overwrite`, `printer`
- [c_codegen.py](../../nrpy/c_codegen.py) - `apply_substitution_dict`
- [cse_preprocess_postprocess.py](../../nrpy/helpers/cse_preprocess_postprocess.py) - `cse_postprocess`
- [SEOBNRv5_aligned_spin_constants.py](../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_constants.py) - `final_mass_non_precessing_UIB2016`
- [SEOBNRv5_aligned_spin_coefficients.py](../../nrpy/infrastructures/JAX/sebob/SEOBNRv5_aligned_spin_coefficients.py) - `register_PyFunction_SEOBNRv5_aligned_spin_coefficients`
- [__init__.py](../../nrpy/__init__.py) - empty top-level package initializer

## See Also

- Depends on: [CSE And Printer Support](helpers/cse-and-printer-support.md)
- Contrasts with: [C Codegen](c-codegen.md)
- See also: [JAX Project Layout And Module Emission](../infrastructures/jax/project-generation-lifecycle.md)

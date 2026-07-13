# Trusted Expression Pipeline

> Explain how symbolic equation outputs become trusted numerical validation files. · Status: confirmed · Last reconciled: 07-12-2026
> Up: [Equations](index.md)

## Summary

Equation modules fingerprint symbolic expressions by evaluating them at a
deterministic high-precision sample point and comparing those values with a
sibling trusted file under `tests/`. This can detect expression drift under the
sampled conditions; it is not a formal symbolic-equality proof. The same
mechanics cover BSSN, GR conversions and diagnostics, GRHD, wave, elliptic,
TOV, SEOBNR/BOB, and most geometry-support helpers. Quaternion tensor rotation
is a confirmed doctest-only exception.

## Detail

The common equation-module flow is:

1. Build a dictionary from object state or explicit expression names.
2. Call `process_dictionary_of_expressions(...)`, usually with
   `fixed_mpfs_for_free_symbols=True`.
3. Call `compare_or_generate_trusted_results(...)` with the owning module path,
   working directory, trusted-file basename, and processed results.

Equation entry points currently pass
`fixed_mpfs_for_free_symbols=True`. For each ordinary free symbol,
`convert_free_symbols_set_to_mpf_dict` seeds Python's `random` module from the
MD5 text digest of that symbol and assigns one `mpf` value in `[0, 1)`. This is
deterministic sampling, not domain-aware test-data selection. `PI` and `M_PI`
instead receive `mp.pi`; `SQRT1_2` and `M_SQRT1_2` receive `1/sqrt(2)`.

`process_dictionary_of_expressions` sorts dictionary items, ignores keys
containing `funcform`, recursively flattens lists, and names flattened entries
`<key>_<flat-index>`. Each SymPy expression is reduced with
`sp.cse(..., order="none")` at the module's default `precision = 30` decimal
digits. Conversion tries `mpf` first and falls back to `mpc` for non-real
results. A small nonzero value below `10**(-4*mp.dps/5)` is retried at twice the
precision and is collapsed to exact zero only when it remains below the tighter
threshold.

`compare_or_generate_trusted_results` derives the owning module's
`tests/<trusted_file_basename>.py` path. If the file exists, it imports
`trusted_dict`, checks dictionary length, then looks up and compares every
trusted key. If the file does not exist, it writes a new trusted dictionary;
that generation branch does not compare against an independent expected value.
With matching dictionary lengths, a renamed or absent result key surfaces as a
`KeyError` during lookup rather than the missing-key diagnostic used for a
length mismatch.

`compare_against_trusted` uses relative tolerance
`10**(-4*mp.dps/5)` against the magnitude of the trusted value after converting
that value to `mpc`. Consequently a trusted zero has a zero relative-error
bound and must match exactly after conversion's near-zero handling. There is no
separate absolute tolerance.

`assert_equal` is a related sampled-numerical check, despite its name. It
processes both inputs with fixed substitutions and compares zipped result values
using relative tolerance against their average. It does not establish a SymPy
identity, and for dictionary inputs it does not independently assert matching
key names or equal lengths. The GR `nrpylatex/test_parse_BSSN.py` cross-check
uses this helper, so that test is also a deterministic sampled comparison of the
parsed and handwritten expression sets.

`output_trusted` writes only the needed `mpmath` imports plus `trusted_dict`, and
formats the file with Black. A count or value mismatch tells maintainers to
delete the stale trusted file and rerun the owning module only when the new
result is independently understood and accepted.

Trusted-value files under `*/tests/*.py` are treated specially. They should
contain only generated `mpf` or `mpc` data, no module docstrings, no functions,
and no classes. The preserved agent rules say not to hand-edit trusted values;
regenerate them from the owning module and explain the reason in the commit
message.

Family pages own the implementation-specific validation inventory. In compact
form, current coverage includes BSSN quantities/RHSs/constraints, ADM/BSSN and
four-metric conversions, initial data, Psi4/tetrads, geodesics, horizon
diagnostics, Fishbone-Moncrief data, GRHD equations/speeds/HLL fluxes, scalar
wave RHSs and initial data, conformally flat elliptic RHS/source terms, TOV ODE
RHSs, SEOBNRv5/BOB dynamics and waveform quantities, basis transforms,
GeneralRFM fisheye maps, SO(3) rotations, and spin-weighted spherical
harmonics. The owning leaves link to representative trusted files and stable
symbols; this page owns only the common processing and comparison mechanics.

These evidence levels are distinct:

- `assert_equal` or `check_zero` is a sampled numerical expression comparison,
  not formal symbolic equality.
- Trusted-value matching says current sampled values match stored sampled
  values for exported keys and selected options.
- Successful Python construction or trusted matching does not prove generated
  C/C++ builds.
- A successful generated-project build does not prove the executable runs or
  its output is scientifically correct.
- Runtime and numerical-accuracy claims require separate infrastructure,
  example, or validation evidence with exact inputs and assertions.

## Sources

- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `precision`, `convert_free_symbols_set_to_mpf_dict`, `convert_one_expression_to_mpfmpc`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `assert_equal`, `check_zero`, `output_trusted`, `compare_against_trusted`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `Required Checks`, `Expression Validation`
- [test_parse_BSSN.py](../../nrpy/equations/general_relativity/nrpylatex/test_parse_BSSN.py) - `test_example_BSSN`
- [BSSN_RHSs.py](../../nrpy/equations/general_relativity/BSSN_RHSs.py) - `BSSNRHSs`, `BSSN_RHSs`
- [BSSN_RHSs_Cartesian.py](../../nrpy/equations/general_relativity/tests/BSSN_RHSs_Cartesian.py) - `trusted_dict`
- [GRHD_equations.py](../../nrpy/equations/grhd/GRHD_equations.py) - `GRHD_Equations`, `construct_all_equations`
- [WaveEquation_RHSs.py](../../nrpy/equations/wave_equation/WaveEquation_RHSs.py) - `WaveEquation_RHSs`
- [SEOBNRv5_aligned_spin_Hamiltonian.py](../../nrpy/equations/seobnr/SEOBNRv5_aligned_spin_Hamiltonian.py) - `SEOBNRv5_aligned_spin_Hamiltonian_quantities`
- [jacobians.py](../../nrpy/equations/basis_transforms/jacobians.py) - `BasisTransforms`
- [tensor_rotation.py](../../nrpy/equations/quaternion_rotations/tensor_rotation.py) - `rotate`

## See Also

- Parent: [Equations](index.md)
- Validated by: [Expression Validation Helpers](../validation/expression-validation-helpers.md)
- Example: [BSSN Family](general-relativity/bssn-family.md)
- Example: [GRHD](grhd.md)
- Example: [Wave Equation](wave-equation.md)
- Example: [SEOBNR And BOB](seobnr/index.md)
- Example: [Geometry And Special-Function Support](geometry-and-special-function-support.md)

# Trusted Expression Pipeline

> Explain how symbolic equation outputs become trusted numerical validation files. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [Equations](index.md)

## Summary

Equation modules validate symbolic expressions by building a results dictionary,
converting it to deterministic high-precision numerical values, and comparing
or generating a sibling trusted file under `tests/`. The same mechanics cover
BSSN, GR conversions and diagnostics, GRHD, wave, elliptic, TOV, SEOBNR/BOB,
and most geometry-support helpers. Quaternion tensor rotation is a confirmed
doctest-only exception.

## Detail

The common equation-module flow is:

1. Build a dictionary from object state or explicit expression names.
2. Call `process_dictionary_of_expressions(...)`, usually with
   `fixed_mpfs_for_free_symbols=True`.
3. Call `compare_or_generate_trusted_results(...)` with the owning module path,
   working directory, trusted-file basename, and processed results.

`process_dictionary_of_expressions` sorts dictionary items, ignores keys
containing `funcform`, flattens tensor lists, and converts every SymPy expression
to an `mpf` or `mpc` result. `compare_or_generate_trusted_results` derives the
owning module's `tests/<trusted_file_basename>.py` path. If the file exists, it
loads `trusted_dict` and compares expression count and values; if it does not
exist, it writes a new trusted dictionary.

`output_trusted` writes only the needed `mpmath` imports plus `trusted_dict`, and
formats the file with Black. `compare_against_trusted` raises on missing keys or
value mismatch and tells maintainers to delete the stale trusted file and rerun
the owning module only when the new result is trusted.

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

## Sources

- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `output_trusted`, `compare_against_trusted`
- [original-agents.md](../../raw/source-docs/original-agents.md) - `Required Checks`, `Expression Validation`
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

# Trusted Expression Pipeline

> Explain how symbolic equation outputs become trusted numerical validation files. · Status: confirmed · Last reconciled: 07-20-2026
> Up: [Equations](index.md)

## Summary

When exact or semantic invariants are impractical, equation modules can sample
symbolic expressions at a deterministic high-precision point and compare those
values with a sibling trusted file under `tests/`. This can detect expression
drift under sampled conditions; it is not a formal symbolic-equality proof. The
same flow covers BSSN, GR conversions and diagnostics, GRHD, wave, elliptic,
TOV, SEOBNR/BOB, and most geometry-support helpers. Quaternion tensor rotation
is a confirmed doctest-only exception.

## Detail

Prefer an exact analytic, symbolic, or semantic invariant. When sampled
regression is appropriate, the common equation-module flow is:

1. Build a dictionary with stable descriptive keys from object state or
   explicit expression names.
2. Call `process_dictionary_of_expressions(...)`, usually with
   `fixed_mpfs_for_free_symbols=True`.
3. Call `compare_or_generate_trusted_results(...)` with the owning module path,
   working directory, trusted-file basename, and processed results.

This flow is sampled numerical regression, not a symbolic identity proof.
[Expression Validation
Helpers](../validation/expression-validation-helpers.md) owns substitution,
conversion, structure and collision handling, tolerances, non-finite behavior,
and trusted-file mechanics. [Test Oracles And Safe
Updates](../validation/test-oracles-and-safe-updates.md) owns store admission
and safe regeneration; [Code Test Policy](../validation/code-test-policy.md)
owns test placement and meaningfulness.

Family pages own the implementation-specific validation inventory. In compact
form, current coverage includes BSSN quantities/RHSs/constraints, ADM/BSSN and
four-metric conversions, initial data, Psi4/tetrads, geodesics, horizon
diagnostics, Fishbone-Moncrief data, GRHD equations/speeds/HLL fluxes, scalar
wave RHSs and initial data, conformally flat elliptic RHS/source terms, TOV ODE
RHSs, SEOBNRv5/BOB dynamics and waveform quantities, basis transforms,
GeneralRFM fisheye maps, SO(3) rotations, and spin-weighted spherical
harmonics. The owning leaves link to representative trusted files and stable
symbols; this page owns the common caller flow and family inventory.

## Sources

- [validate_expressions.py](../../nrpy/validate_expressions/validate_expressions.py) - `process_dictionary_of_expressions`, `compare_or_generate_trusted_results`
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
- Depends on: [Test Oracles And Safe Updates](../validation/test-oracles-and-safe-updates.md)
- Depends on: [Code Test Policy](../validation/code-test-policy.md)
- Example: [BSSN Family](general-relativity/bssn-family.md)
- Example: [GRHD](grhd.md)
- Example: [Wave Equation](wave-equation.md)
- Example: [SEOBNR And BOB](seobnr/index.md)
- Example: [Geometry And Special-Function Support](geometry-and-special-function-support.md)

# GR Poststep Repairs

> ETLegacy lapse flooring, initial determinant/trace projection, and determinant-only MoL repair. · Status: confirmed · Last reconciled: 08-26-2026
> Up: [ETLegacy](index.md)

## Summary

This leaf owns ETLegacy GR repair kernels. Registration places the lapse floor
in `MoL_PostStep`, the combined determinant/trace projection in
`CCTK_INITIAL`, and a determinant-only repair in `MoL_PostStep` and
`MoL_PseudoEvolution`. The split preserves ETLegacy's established evolved
solution while making freshly converted initial data trace-free. Generated
Baikal and BaikalVacuum thorns were built and all five ET tests passed; restart
validation was not done.

These functions are repair kernels around the evolved state. They do not own
the symbolic BSSN equations; see
[GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
and [BSSN Family](../../equations/general-relativity/bssn-family.md) for the
evolution and constraint expressions.

## Detail

`register_CFunction_floor_the_lapse()` emits `<thorn>_floor_the_lapse`. Its
generated body defines `NRPYMAX` if needed, then loops over all points with
SIMD disabled and assigns `alphaGF = NRPYMAX(alphaGF, lapse_floor)`. The
registered schedule entry is in `MoL_PostStep` before
`<thorn>_enforce_detgbar_equals_detghat`, reads and writes `alphaGF`
everywhere, and records `lapse_floor` as a code parameter used by the thorn.

`register_CFunction_enforce_detgbar_equals_detghat_trAzero()` emits
`<thorn>_enforce_detgbar_equals_detghat_trAzero`. It selects the BSSN and
reference-metric objects for either the plain coordinate system or the
`_rfm_precompute` variant. In one all-points loop it loads all six independent
`hDD` and six independent `aDD` components, reconstructs `gammabarDD` and
`AbarDD`, and rescales `gammabarDD` by
`(|detgammahat|/detgammabar)^(1/3)`. It then inverts that corrected metric and
removes `gammabar^{ij} Abar_ij / 3` from `AbarDD`. Only after both corrected
tensors are formed does it write all twelve evolved components. Thus the
trace projection uses the post-enforced conformal metric and cannot observe
partially updated input.

The combined projection runs only in `CCTK_INITIAL`, after every registered
ADM-to-BSSN producer. The same registration function also emits
`<thorn>_enforce_detgbar_equals_detghat`, which rescales only the six `hDD`
components. This determinant-only repair runs after evolved-variable boundary
conditions in `MoL_PostStep` and before BSSN constraints in
`MoL_PseudoEvolution`. BSSN-to-ADM consumers run after that repair. SIMD
remains disabled because both kernels use a cube root and write all points.

Claim evidence:
- Claim: ETLegacy registers the combined determinant/trace projection in `CCTK_INITIAL` and a determinant-only repair in `MoL_PostStep` and `MoL_PseudoEvolution`; exact ET 2024-06 beta generation, build, and five-test execution pass, while recovery and restart behavior remain untested.
- Role: generated evidence
- Deciding authority: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: [floor_the_lapse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py), `register_CFunction_floor_the_lapse`
- Validation: `inspected=pass; generated=pass; built=pass; run=pass; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.12.1, ET 2024-06 beta; backend=ETLegacy generated Baikal, BaikalVacuum, and WaveToyNRPy thorns; precision=exact determinant/trace identities, generated-source structure, and ET trusted-output tolerances; GPU=not-run; restart=not-run; distributed=2 MPI processes; error_path=not-run; options=Cartesian ET tests plus Cartesian, SinhSpherical precompute, and GeneralRFM generator validation; date=08-26-2026`

## Sources

- [floor_the_lapse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py) - `register_CFunction_floor_the_lapse`
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

## See Also

- [ETLegacy](index.md)
- [GR ADM/BSSN, Slicing, And Matter Coupling](gr-adm-bssn-slicing-and-matter-coupling.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Contrasts with: [CarpetX GR Poststep Repairs](../carpetx/gr-poststep-repairs.md)

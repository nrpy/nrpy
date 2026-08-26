# GR Poststep Repairs

> ETLegacy lapse flooring and combined conformal determinant/trace projection across initial and MoL lifecycle bins. · Status: confirmed · Last reconciled: 08-26-2026
> Up: [ETLegacy](index.md)

## Summary

This leaf owns ETLegacy GR repair kernels. Registration places the lapse floor
in `MoL_PostStep` and emits combined-projection entries for `CCTK_INITIAL`,
`MoL_PostStep`, and `MoL_PseudoEvolution`. Those entries restore
`det(gammabar)=det(gammahat)` and `tr(Abar)=0` when executed. No generated
thorn build, initial/recovery run, or restart validation was done.

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
`<thorn>_enforce_detgbar_equals_detghat_trAzero`, reads and writes `alphaGF`
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

The combined projection has entries in `CCTK_INITIAL`, `MoL_PostStep`, and
`MoL_PseudoEvolution`. Poststep it runs after the lapse floor; pseudo-evolution
it runs after auxiliary BCs and before BSSN constraints.
Schedule metadata reads and writes all twelve gridfunctions everywhere. SIMD
remains disabled because the kernel uses a cube root and writes all points.

Claim evidence:
- Claim: ETLegacy registers the combined projection in `CCTK_INITIAL`, `MoL_PostStep`, and `MoL_PseudoEvolution`, with emitted ordering after the named BC aliases and before pseudo-evolution constraints; source registration does not prove fresh-data, recovery, or restart execution.
- Role: generated evidence
- Deciding authority: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: [floor_the_lapse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py), `register_CFunction_floor_the_lapse`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=ETLegacy generated projection kernel and schedule metadata; precision=exact determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=Cartesian and SinhSpherical reference-metric precompute kernel variants; date=08-26-2026`

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

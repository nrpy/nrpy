# CarpetX GR Poststep Repairs

> CarpetX lapse flooring and combined conformal determinant/trace projection across initial and ODESolvers poststep bins. · Status: confirmed · Last reconciled: 08-26-2026
> Up: [CarpetX](index.md)

## Summary

This leaf owns the CarpetX lapse floor and combined conformal
determinant/trace projection. Registration emits entries in `CCTK_INITIAL` and
`ODESolvers_PostStep`. At poststep, lapse flooring comes first; emitted ADM
BSSN-constraint schedule text names projection as a predecessor. ADM export
ordering is owned by [GR ADM/BSSN And Matter Coupling](gr-adm-bssn-and-matter-coupling.md).
No generated thorn build, initial/recovery run, or restart validation was done.

## Detail

`register_CFunction_floor_the_lapse()` emits
`<thorn>_floor_the_lapse` in `ODESolvers_PostStep` before
`<thorn>_enforce_detgbar_equals_detghat_trAzero`. The generated loop covers all points,
reads and writes `alphaGF(everywhere)`, defines `NRPYMAX` when needed, and
assigns `alphaGF = NRPYMAX(alphaGF, lapse_floor)`. The registration also marks
`lapse_floor` in `ET_current_thorn_CodeParams_used`, so CarpetX parameter CCL
generation can see that runtime code parameter.

`register_CFunction_enforce_detgbar_equals_detghat_trAzero()` emits
`<thorn>_enforce_detgbar_equals_detghat_trAzero` in `CCTK_INITIAL` and
`ODESolvers_PostStep`. It
selects coordinate or `_rfm_precompute` objects and uses one all-points loop.
The loop first loads all six independent `hDD` and six independent `aDD`
components. It reconstructs both tensors, rescales `gammabarDD` by
`(abs(detgammahat)/detgammabar)^(1/3)`, inverts the corrected metric, and then
subtracts `gammabar_ij tr(Abar)/3` from `AbarDD`. It writes all twelve corrected
components only after forming both outputs, so the trace uses the post-enforced
metric and no store can corrupt a later input.

Schedule metadata reads and writes all twelve gridfunctions everywhere. SIMD
is disabled for the cube-root, all-points kernel. CarpetX schedules BSSN
constraint diagnostics after poststep projection.

Claim evidence:
- Claim: CarpetX registers the combined projection in `CCTK_INITIAL` and `ODESolvers_PostStep`, with poststep schedule text after lapse flooring and before BSSN constraints; source registration does not prove fresh-data, recovery, or restart execution.
- Role: generated evidence
- Deciding authority: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/CarpetX/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: [floor_the_lapse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/floor_the_lapse.py), `register_CFunction_floor_the_lapse`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=CarpetX generated projection kernel and schedule metadata; precision=exact determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=Cartesian and SinhSpherical reference-metric precompute kernel variants; date=08-26-2026`

## Sources

- [floor_the_lapse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/floor_the_lapse.py) - `register_CFunction_floor_the_lapse`
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/CarpetX/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

## See Also

- [CarpetX](index.md)
- [GR ADM/BSSN And Matter Coupling](gr-adm-bssn-and-matter-coupling.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- Contrasts with: [ETLegacy GR Poststep Repairs](../etlegacy/gr-poststep-repairs.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)

# CarpetX GR Poststep Repairs

> CarpetX lapse flooring and conformal-metric determinant repair in `ODESolvers_PostStep`. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [CarpetX](index.md)

## Summary

This leaf owns two CarpetX GR repair kernels that run after ODESolvers steps:
the lapse floor and the conformal-metric determinant enforcement. The lapse
floor is scheduled before determinant enforcement, so `alphaGF` is repaired
before the `hDD*` conformal-metric perturbations are rescaled.

## Detail

`register_CFunction_floor_the_lapse()` emits
`<thorn>_floor_the_lapse` in `ODESolvers_PostStep before
<thorn>_enforce_detgammahat_constraint`. The generated loop covers all points,
reads and writes `alphaGF(everywhere)`, defines `NRPYMAX` when needed, and
assigns `alphaGF = NRPYMAX(alphaGF, lapse_floor)`. The registration also marks
`lapse_floor` in `ET_current_thorn_CodeParams_used`, so CarpetX parameter CCL
generation can see that runtime code parameter.

`register_CFunction_enforce_detgammahat_constraint()` emits
`<thorn>_enforce_detgammahat_constraint` in `ODESolvers_PostStep`. It selects
the BSSN and reference-metric objects for either the coordinate system or the
`_rfm_precompute` variant, computes `detgammabar` from `gammabarDD`, and writes
the six independent `hDD*` gridfunctions with corrected values. The correction
rescales `delta + hDD` by `(abs(detgammahat) / detgammabar)^(1/3)` and then
subtracts `delta`, which restores `det(gammabar) = det(gammahat)` without
duplicating the symbolic BSSN derivation here.

The determinant repair loops over all points, reads
`hDD00GF`, `hDD01GF`, `hDD02GF`, `hDD11GF`, `hDD12GF`, and `hDD22GF`, and
writes the same six gridfunctions everywhere. SIMD is disabled both in the
`c_codegen()` call and the CarpetX `simple_loop()` call. The source comments
explain that choice with cube-root support and all-points memory-safety
concerns.

Within `ODESolvers_PostStep`, this repair pair is the local source-backed
precondition for the later BSSN-to-ADM export path: the lapse is floored first,
then determinant enforcement repairs `hDD*`, and the ADM/BSSN handoff page owns
the downstream export details.

## Sources

- [floor_the_lapse.py](../../../nrpy/infrastructures/CarpetX/general_relativity/floor_the_lapse.py) - `register_CFunction_floor_the_lapse`
- [enforce_detgammahat_constraint.py](../../../nrpy/infrastructures/CarpetX/general_relativity/enforce_detgammahat_constraint.py) - `register_CFunction_enforce_detgammahat_constraint`

## See Also

- [CarpetX](index.md)
- [GR ADM/BSSN And Matter Coupling](gr-adm-bssn-and-matter-coupling.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- [ETLegacy GR Poststep Repairs](../etlegacy/gr-poststep-repairs.md)
- [BSSN Family](../../equations/general-relativity/bssn-family.md)

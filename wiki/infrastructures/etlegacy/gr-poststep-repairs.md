# GR Poststep Repairs

> ETLegacy lapse flooring and conformal-metric determinant correction after MoL steps. · Status: confirmed · Last reconciled: 2026-06-30
> Up: [ETLegacy](index.md)

## Summary

This leaf owns the corrective ETLegacy GR kernels that run in `MoL_PostStep`.
The lapse-floor kernel applies a runtime floor to `alpha` before the
determinant-enforcement kernel rescales conformal-metric perturbations so the
BSSN determinant condition is restored.

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
`<thorn>_enforce_detgammahat_constraint`, reads and writes `alphaGF`
everywhere, and records `lapse_floor` as a code parameter used by the thorn.

`register_CFunction_enforce_detgammahat_constraint()` emits
`<thorn>_enforce_detgammahat_constraint`. It selects the BSSN and
reference-metric objects for either the plain coordinate system or the
`_rfm_precompute` variant, computes `detgammabar` from `gammabarDD`, and
updates each independent `hDD` component with
`(|detgammahat| / detgammabar)^(1/3) * (delta + hDD) - delta`. That rescaling
is the source-level repair applied to the conformal-metric perturbation so the
reconstructed `gammabarDD` satisfies `det(gammabar) = det(gammahat)`.

The determinant repair is also an all-points `MoL_PostStep` kernel with SIMD
disabled. It reads the six independent `hDD` gridfunctions and writes the same
six gridfunctions everywhere. The implementation comments explain the SIMD
choice in terms of cube-root support and all-points memory safety; this page
therefore documents the actual generated expression and schedule, not unrelated
or misleading docstring text in neighboring functions.

## Sources

- [floor_the_lapse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py) - `register_CFunction_floor_the_lapse`
- [enforce_detgammahat_constraint.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgammahat_constraint.py) - `register_CFunction_enforce_detgammahat_constraint`
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md) - gridfunction and code-parameter ownership

## See Also

- [ETLegacy](index.md)
- [GR ADM/BSSN, Slicing, And Matter Coupling](gr-adm-bssn-slicing-and-matter-coupling.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- [Gridfunctions And Parameters](../../core/gridfunctions-and-parameters.md)

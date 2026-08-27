# GR Poststep Repairs

> ETLegacy lapse flooring and combined conformal determinant/trace projection across initial and MoL lifecycle bins. · Status: confirmed · Last reconciled: 08-26-2026
> Up: [ETLegacy](index.md)

## Summary

This leaf owns ETLegacy GR repair kernels. Registration places the lapse floor
in `MoL_PostStep` and emits combined-projection entries for `CCTK_INITIAL`,
`MoL_PostStep`, and `MoL_PseudoEvolution`. For finite, nonsingular states where
the projection expressions are defined, those entries restore
`det(gammabar)=det(gammahat)` and `tr(Abar)=0` to roundoff. Registration emits
no determinant-only ETLegacy repair entry. No generated thorn build,
initial/recovery run, restart validation, or singular/nonfinite error-path test
was done.

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
`(|detgammahat|/detgammabar)^(1/3)`. For GeneralRFM it explicitly inverts that
corrected metric before removing the trace from `AbarDD`. Standard orthogonal
reference metrics instead invert the reduced `H=I+h` matrix and use algebraic
cancellation of the rescaling factor, yielding the equivalent trace removal
relative to the corrected metric. The generated block loads all six `hDD` and
six `aDD` components before its first store, computes six corrected metric and
six corrected curvature outputs, then stores all twelve, so the projection
cannot observe partially updated input.

The combined projection has entries in `CCTK_INITIAL`, `MoL_PostStep`, and
`MoL_PseudoEvolution`. Poststep it runs after evolved-variable boundary
conditions and after the lapse floor; pseudo-evolution it runs before BSSN
constraints without an auxiliary-boundary dependency.
Schedule metadata lists all twelve gridfunctions as reads and marks all twelve
writes `everywhere`. SIMD remains disabled; source comments cite cube-root
support plus all-points memory and thread-safety concerns. All three repair
schedules invoke this combined projector; none invokes a determinant-only
repair.

Claim evidence:
- Claim: The inspected ETLegacy projector registration assigns the combined determinant/trace projector to `CCTK_INITIAL`, after evolved-variable boundary conditions in `MoL_PostStep`, and before constraints without an auxiliary-boundary dependency in `MoL_PseudoEvolution`; that registration emits no determinant-only repair entry, and source inspection does not prove fresh-data, recovery, restart, or singular/nonfinite error-path behavior.
- Role: descriptive behavior
- Deciding authority: [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/enforce_detgbar_equals_detghat_trAzero.py), `register_CFunction_enforce_detgbar_equals_detghat_trAzero`
- Corroboration: [floor_the_lapse.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/floor_the_lapse.py), `register_CFunction_floor_the_lapse`
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, SymPy 1.14.0; backend=ETLegacy generated projection kernel and schedule metadata; precision=exact determinant/trace identities and generated-source structure; GPU=not-run; restart=not-run; distributed=not-run; error_path=not-run; options=Cartesian, SinhSpherical reference-metric precompute, and GeneralRFM kernel variants; date=08-26-2026`

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

# GR ADM/BSSN, Slicing, And Matter Coupling

> ETLegacy initialization, ADMBase export, slicing registration, and stress-energy raising. · Status: confirmed · Last reconciled: 07-06-2026
> Up: [ETLegacy](index.md)

## Summary

This leaf owns ETLegacy runtime glue around the symbolic GR equation modules.
It covers conversion between Einstein Toolkit ADMBase variables and NRPy BSSN
evolved variables, the startup slicing registration hook, and the TmunuBase
stress-energy raise used before matter-coupled RHS and constraint kernels.

Symbolic ADM/BSSN and matter formulas remain owned by
[Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
and [BSSN Family](../../equations/general-relativity/bssn-family.md). This
page documents where ETLegacy schedules those formulas and which Cactus
gridfunction groups they read or write.

## Detail

`register_CFunction_ADM_to_BSSN()` emits a thorn-local C function named
`<thorn>_ADM_to_BSSN_order_<fd_order>`. Its first all-points loop reads
ADMBase metric, shift, extrinsic-curvature, time-derivative-of-shift, and lapse
data, copies the ADM lapse into the BSSN `alpha` gridfunction, and uses the
symbolic `ADM_to_BSSN` helper to write `cf`, `trK`, `vetU*`, `betU*`, `hDD*`,
and `aDD*`. The generated schedule runs at `CCTK_INITIAL` after
`ADMBase_PostInitial`, guarded by the requested finite-difference order, and
writes/syncs ETLegacy `evol_variables`.

After that ADM-to-BSSN conversion, the same registration function computes
`LambdabarU` from `gammabarUU`, `GammabarUDD`, and the reference-metric
`GammahatUDD`, rescales it into `lambdaU*`, evaluates those components on the
interior, and calls `ExtrapolateGammas` for `lambdaU0GF`, `lambdaU1GF`, and
`lambdaU2GF`. This makes the lambda gridfunctions a post-conversion derived
state, not independent ADMBase input.

`register_CFunction_BSSN_to_ADM()` emits `<thorn>_BSSN_to_ADM`. It copies BSSN
`alpha`, `vetU*`, and `betU*` back to ADMBase lapse, shift, and `dtshift`, then
uses the symbolic `BSSN_to_ADM` helper to reconstruct ADMBase metric and
extrinsic curvature from `cf`, `trK`, `hDD*`, and `aDD*`. Its normal schedule
entry is in `MoL_PostStep`, after `<thorn>_evol_ApplyBCs` and before
`ADMBase_SetADMVars`, where it writes ADMBase metric, shift, curvature,
`dtshift`, and lapse everywhere. When `thorn_name == "Baikal"`, the function
also registers a `MoL_PseudoEvolution` entry after `<thorn>_aux_ApplyBCs`; the
source schedule text marks that path as needed for HydroBase integration.

`register_CFunction_T4DD_to_T4UU()` emits `<thorn>_T4DD_to_T4UU`. It builds the
contravariant four-metric through `BSSN_to_g4UU`, reads TmunuBase
stress-energy components from the `eT*` gridfunctions, and raises both tensor
indices into ETLegacy `T4UU*` gridfunctions. The function is scheduled twice:
in `MoL_CalcRHS` before `<thorn>_RHS`, so matter-coupled RHS code sees raised
stress-energy data, and in `MoL_PseudoEvolution` before
`<thorn>_BSSN_constraints`, so constraint code sees the same raised data. The
implementation disables SIMD in both `c_codegen` and the all-points loop; its
description warns that total local grid size is not guaranteed to be a multiple
of SIMD width.

`register_CFunction_RegisterSlicing()` emits `<thorn>_RegisterSlicing`, an
`int` function scheduled at `STARTUP` with Cactus `meta` options. Its body calls
`Einstein_RegisterSlicing("<thorn>")` and returns zero, registering the
generated thorn as a 3+1 slicing provider before later schedule bins run.

## Sources

- [ADM_to_BSSN.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/ADM_to_BSSN.py) - `register_CFunction_ADM_to_BSSN`
- [BSSN_to_ADM.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/BSSN_to_ADM.py) - `register_CFunction_BSSN_to_ADM`
- [T4DD_to_T4UU.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/T4DD_to_T4UU.py) - `register_CFunction_T4DD_to_T4UU`
- [RegisterSlicing.py](../../../nrpy/infrastructures/ETLegacy/general_relativity/RegisterSlicing.py) - `register_CFunction_RegisterSlicing`

## See Also

- [ETLegacy](index.md)
- [GR BSSN RHS, Ricci, Constraints, And Validation](gr-bssn-rhs-ricci-constraints-and-validation.md)
- [GR Poststep Repairs](gr-poststep-repairs.md)
- [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- [BSSN Family](../../equations/general-relativity/bssn-family.md)

# CarpetX GR ADM/BSSN And Matter Coupling

> CarpetX ADMBaseX import/export and TmunuBaseX stress-energy raising around BSSN evolution. · Status: confirmed
> Up: [CarpetX](index.md)

## Summary

This leaf covers CarpetX runtime transfers between external Einstein Toolkit
interfaces and NRPy BSSN variables. `ADM_to_BSSN` imports ADMBaseX initial data
into CarpetX `evol_variables`, `BSSN_to_ADM` exports poststep BSSN state back
to ADMBaseX, and `T4DD_to_T4UU` raises TmunuBaseX covariant stress-energy data
for matter-coupled RHS and constraint kernels.

Symbolic ADM/BSSN and stress-energy formulas are defined in
[Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
and [BSSN Family](../../equations/general-relativity/bssn-family.md); the
`CarpetX/general_relativity/` registrars add ODESolvers schedule bins,
read/write group declarations, and tile/subgrid loop transfers.

## Detail

`register_CFunction_ADM_to_BSSN()` emits
`<thorn>_ADM_to_BSSN_order_<fd_order>` and temporarily sets the global
`fd_order` parameter while building finite-difference code for the requested
order. The schedule entry is guarded by `if(fd_order == <fd_order>)` and runs
at `CCTK_INITIAL after ADMBaseX_PostInitial`. It reads ADMBaseX metric, shift,
curvature (`ADMBaseX::curv`), `dtshift`, and lapse groups, writes
`evol_variables`, and syncs `evol_variables`.

Inside the first all-points loop, `ADM_to_BSSN` copies ADMBaseX lapse `alp`
into BSSN `alpha`, reads ADMBaseX shift, `dtshift`, metric, and curvature
components, calls the symbolic ADM-to-BSSN helper in Cartesian coordinates, and
writes `cf`, `trK`, `vetU*`, `betU*`, `hDD*`, and `aDD*`. It then computes
`LambdabarU` from `gammabarUU`, `GammabarUDD`, and reference-metric
`GammahatUDD`, rescales by `ReU`, writes `lambdaU0..2` on the interior, and
leaves the CarpetX `ExtrapolateGammas` calls commented in the generated body.

`register_CFunction_BSSN_to_ADM()` emits `<thorn>_BSSN_to_ADM` for the reverse
handoff. The function copies BSSN `alpha`, `vetU*`, and `betU*` into ADMBaseX
lapse, shift, and `dtshift`, reconstructs ADM metric and extrinsic curvature
from `cf`, `trK`, `hDD*`, and `aDD*`, and writes the ADMBaseX groups
everywhere. Its schedule runs in `ODESolvers_PostStep` after
`<thorn>_enforce_detgbar_equals_detghat_trAzero` and before
`ADMBaseX_SetADMVars`, so exported ADM data consumes both the determinant- and
trace-projected conformal tensors.

Evidence (generated evidence): `CarpetX/general_relativity/BSSN_to_ADM.py::register_CFunction_BSSN_to_ADM`; corroboration: `CarpetX/general_relativity/enforce_detgbar_equals_detghat_trAzero.py::register_CFunction_enforce_detgbar_equals_detghat_trAzero`.

`register_CFunction_T4DD_to_T4UU()` emits `<thorn>_T4DD_to_T4UU`. It builds
`g4UU` from BSSN variables through the symbolic BSSN-to-four-metric conversion,
reads `TmunuBaseX::stress_energy_scalar`,
`TmunuBaseX::stress_energy_vector`, and
`TmunuBaseX::stress_energy_tensor` as `eTtt`, `eTtx`, and related components,
raises both indices, and writes `T4UU00`, `T4UU01`, `T4UU02`, `T4UU03`,
`T4UU11`, `T4UU12`, `T4UU13`, `T4UU22`, `T4UU23`, and `T4UU33` everywhere.

The matter conversion schedules once before `<thorn>_RHS` so matter-coupled
RHS code can consume raised stress-energy data, and once before
`<thorn>_BSSN_constraints` so constraint code sees the same raised data. The
source disables SIMD in both `c_codegen` and the all-points loop because
`cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]` is not guaranteed to be a multiple of
`SIMD_width`.

In `T4DD_to_T4UU.py`, the Python tuple keys are `ODESolvers_CalcRHS` and
`ODESolvers_PseudoEvolution`, whereas the emitted `schedule.ccl` strings place
the routines in `ODESolvers_RHS` and `ODESolvers_PostStep`.

## Sources

- [ADM_to_BSSN.py](../../../nrpy/infrastructures/CarpetX/general_relativity/ADM_to_BSSN.py) - `register_CFunction_ADM_to_BSSN`
- [BSSN_to_ADM.py](../../../nrpy/infrastructures/CarpetX/general_relativity/BSSN_to_ADM.py) - `register_CFunction_BSSN_to_ADM`
- [T4DD_to_T4UU.py](../../../nrpy/infrastructures/CarpetX/general_relativity/T4DD_to_T4UU.py) - `register_CFunction_T4DD_to_T4UU`
- [enforce_detgbar_equals_detghat_trAzero.py](../../../nrpy/infrastructures/CarpetX/general_relativity/enforce_detgbar_equals_detghat_trAzero.py) - `register_CFunction_enforce_detgbar_equals_detghat_trAzero`

## See Also

- Depends on: [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [Reference Metrics](../../core/reference-metrics.md)

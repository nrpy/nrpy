# BHaHAHA Horizon Runtime

> Orient the BHaH BHaHAHA runtime and route detailed apparent-horizon work to the owning leaves. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

BHaHAHA is the BHaH apparent-horizon finder integration. A generated BHaH
application registers per-horizon data and diagnostics storage, prepares metric
input for each active horizon, and calls the BHaHAHA solver through its public C
boundary.

The runtime has two layers. `bhahaha_find_horizons` is the BHaH-facing
orchestrator for multiple horizons, generated controls, metric adaptation, and
persistent horizon history. `bah_find_horizon` is the public single-horizon
entry point backed by `find_horizon`, which runs BHaHAHA's multiresolution
relaxation and returns a solver status.

## Detail

Use this page to choose the next owner:

- [BHaHAHA Public API And Input Contract](bhahaha-public-api-and-input-contract.md)
  owns `BHaHAHA_header.h`, caller-owned inputs, generated BHaH controls, public
  helper prototypes, poisoning checks, and solver status vocabulary.
- [BHaHAHA Grid, Interpolation, And Boundaries](bhahaha-grid-interpolation-and-boundaries.md)
  owns the BHaH BSSN-to-ADM adapter, spherical input-grid setup, interpolation
  stages, BHaHAHA-local grid setup, and boundary handling.
- [BHaHAHA Relaxation, Diagnostics, And State](bhahaha-relaxation-diagnostics-and-state.md)
  owns the single-horizon relaxation loop, multiresolution progression, initial
  surface history, BBH/common-horizon state, diagnostics, output handoff, and
  temporary-memory cleanup.

The overview boundary is intentionally narrow. Facts about exact struct fields,
input-grid layout, interpolation enums, boundary algorithms, residual tests,
diagnostic quantities, or checkpoint serialization belong in the task leaves or
their BHaH neighbors, not here.

## Sources

- [BHaH_implementation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `bhahaha_find_horizons`, `initialize_bhahaha_solver_params_and_shapes`, `BHaHAHA_interpolate_metric_data_nrpy`, `BHaHAHA_BSSN_to_ADM_Cartesian`, `free_bhahaha_horizon_shape_data_all_horizons`
- [find_horizon.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/find_horizon.py) - `register_CFunction_find_horizon`, `find_horizon`, `free_all_but_external_input_gfs`
- [BHaHAHA_header.h](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`, `NUM_EXT_INPUT_CARTESIAN_GFS`, `bah_find_horizon`

## See Also

- [BHaH](index.md)
- [BHaHAHA Public API And Input Contract](bhahaha-public-api-and-input-contract.md)
- [BHaHAHA Grid, Interpolation, And Boundaries](bhahaha-grid-interpolation-and-boundaries.md)
- [BHaHAHA Relaxation, Diagnostics, And State](bhahaha-relaxation-diagnostics-and-state.md)

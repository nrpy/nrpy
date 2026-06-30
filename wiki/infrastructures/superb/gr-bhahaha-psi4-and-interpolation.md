# GR, BHaHAHA, Psi4, And Interpolation

> Explain superB GR service wiring for BHaHAHA horizons, Psi4 shell extraction, and distributed interpolation. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [superB](index.md)

## Summary

superB wires GR observables through service hooks around the generated
timestepping chares. The BHaHAHA path registers horizon-finding parameters and a
two-phase `bhahaha_find_horizons` wrapper, then uses `Horizon_finder` and
`Interpolator3d` chares to stage target points, distribute interpolation work,
return packed BSSN samples, transform those samples into the ADM metric data
expected by BHaHAHA, and run the horizon solve.

Psi4 has a related shell workflow. The spectroscopy example enables Psi4 only
with BHaHAHA, registers the equation-side `psi4` diagnostic and spin-weight
-2 spherical-harmonic helper, registers the superB shell decomposition
function, and asks `Interpolator3d` to handle Psi4 shell interpolation when
`enable_psi4` is true.

## Detail

The examples gate these services explicitly. `superB_two_blackholes_collide.py`
sets `enable_BHaHAHA`, builds or imports the BHaHAHA library path, registers
`bhahaha_find_horizons`, registers the uniform-grid interpolator, and emits the
`Interpolator3d` and `Horizon_finder` service files. The spectroscopy example
adds `enable_psi4`, rejects Psi4 without BHaHAHA, registers the BHaH Psi4
diagnostic and spin-weighted spherical-harmonic C functions, registers the
superB `psi4_spinweightm2_decomposition` function, and passes `enable_psi4` to
`output_interpolator3d_h_cpp_ci`.

`register_bhahaha_commondata_and_params` adds BHaHAHA state to `commondata`:
per-horizon `bhahaha_params_and_data`, diagnostics, search radii, convergence
tolerances, multigrid angular resolutions, initial center guesses, BBH-mode
indices, and active-horizon bookkeeping. `register_CFunction_bhahaha_find_horizons`
also registers `free_bhahaha_horizon_shape_data_all_horizons`, builds a
prefunction block, and publishes the generated C entry point with the
`BHAHAHA_FIND_HORIZONS_SETUP` and `BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE`
phases defined in `superB.h`.

The BHaHAHA prefunction block is the superB-side adapter, not the solver math.
It defines the BSSN gridfunction list used for interpolation, validates
multigrid resolution inputs, initializes per-horizon solver data and shape
history arrays, prepares spherical interpolation targets around each current
horizon guess, allocates temporary BSSN sample buffers, and converts
interpolated BSSN data into Cartesian ADM `gamma_ij` and `K_ij` components.
`generate_bssn_to_adm_codegen` supplies the symbolic code for that final
BSSN-to-ADM packing step using the active coordinate system and conformal-factor
parameter.

The horizon service flow is split across `Horizon_finder` and `Interpolator3d`.
Each `Horizon_finder` instance calls `bhahaha_find_horizons` in setup mode for
its horizon index, which computes guesses, radial interpolation points, and the
destination coordinate array. It then waits for `ready_for_interpolation` and
asks the `Interpolator3d` array to interpolate `BHAHAHA_NUM_INTERP_GFS` values
with request type `INTERP_REQUEST_BHAHAHA`. After an `InterpBufMsg` returns,
`process_interpolation_results` unpacks the concatenated point-indexed payload
into `dst_data_ptrs`, and the second `bhahaha_find_horizons` phase transforms
the returned BSSN samples to ADM, calls `bah_find_horizon`, writes BHaHAHA
diagnostics on success, records failure state when needed, and frees temporary
metric data.

`Interpolator3d` is the shared distributed interpolation service. Its
`perform_interpolation` method filters requested destination points to the
local chare's grid bounds, allocates per-gridfunction local result arrays, calls
`interpolation_3d_general__uniform_src_grid`, and packs each contribution as
the original destination index plus the requested gridfunction values.
`contribute_interpolation_results` sends those packed buffers to interpolator
chare `(0,0,0)`, `recv_interp_msg` checks that all incoming buffers match the
same request type and id, and `send_interp_concat` concatenates the full result.
For BHaHAHA requests the aggregate message is routed to the matching
`Horizon_finder`; for Psi4 requests it is unpacked locally into shell arrays and
decomposed.

`register_CFunction_psi4_spinweightm2_decomposition` owns the non-chare C
function for shell extraction and mode output. It registers angular shell sizes,
the number of extraction radii, and the extraction-radius array; initializes a
`psi4_shell_angular_grid_t`; fills shell points through
`psi4_spinweightm2_shell_fill_points`; interpolates the real and imaginary
diagnostic gridfunctions when the shell lies inside the grid interior; and calls
`psi4_spinweightm2_decompose_shell` to append `R_ext * psi4_{l,m}` time series.
When `Interpolator3d` is generated with `enable_psi4`, its `.ci` flow instead
uses request type `INTERP_REQUEST_PSI4`, loops over
`commondata.num_psi4_extraction_radii`, builds each shell on interpolator chare
zero, distributes interpolation across the array, unpacks two gridfunctions, and
calls the same shell decomposition helper.

This page intentionally treats BHaHAHA and Psi4 equations as neighboring
context only. The equation-side horizon page owns apparent-horizon geometry and
spin diagnostics; the Psi4 page owns the Weyl scalar and tetrad construction;
this superB page owns how generated services schedule, transport, unpack, and
adapt those quantities at runtime.

## Sources

- [BHaH_implementation.py](../../../nrpy/infrastructures/superB/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `register_bhahaha_commondata_and_params`, `build_bhahaha_prefunc`, `generate_bssn_to_adm_codegen`, `register_CFunction_free_bhahaha_horizon_shape_data_all_horizons`
- [psi4_spinweightm2_decomposition.py](../../../nrpy/infrastructures/superB/general_relativity/psi4_spinweightm2_decomposition.py) - `register_CFunction_psi4_spinweightm2_decomposition`, `psi4_spinweightm2_shell_init`, `psi4_spinweightm2_shell_fill_points`, `psi4_spinweightm2_decompose_shell`
- [interpolator3d_chare.py](../../../nrpy/infrastructures/superB/interpolator3d_chare.py) - `InterpBufMsg`, `Interpolator3d::perform_interpolation`, `Interpolator3d::contribute_interpolation_results`, `output_interpolator3d_h_cpp_ci`
- [horizon_finder_chare.py](../../../nrpy/infrastructures/superB/horizon_finder_chare.py) - `Horizon_finder::process_interpolation_results`, `output_horizon_finder_ci`
- [superB_two_blackholes_collide.py](../../../nrpy/examples/superB_two_blackholes_collide.py) - `enable_BHaHAHA`
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - `enable_psi4`, `enable_BHaHAHA`
- [superB.h](../../../nrpy/infrastructures/superB/superB/superB.h) - `BHAHAHA_FIND_HORIZONS_SETUP`, `BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE`, `psi4_shell_angular_grid_t`, `unpack_interpolation_buffer`

## See Also

- [superB](index.md)
- [BHaH Lifecycle](../bhah-lifecycle.md)
- [Chare Entrypoints And Runtime](chare-entrypoints-and-runtime.md)
- [Diagnostics And Observables](diagnostics-and-observables.md)
- [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- [Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)
- [Horizon Diagnostics](../../equations/general-relativity/horizon-diagnostics.md)
- [Geometry And Special-Function Support](../../equations/geometry-and-special-function-support.md)

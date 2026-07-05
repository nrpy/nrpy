# GR Application Wiring

> Map how BHaH registers generated CFunctions that connect GR equations, initial data, diagnostics, and basis transforms. Status: confirmed. Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

BHaH GR wiring is a code-generation layer over the symbolic GR modules. It does
not rederive the equations. It chooses coordinate systems and feature flags,
pulls symbolic expressions from the BSSN, ADM, Psi4, and initial-data modules,
wraps them in BHaH loop/kernel infrastructure, and registers concrete
CFunctions such as `rhs_eval`, `Ricci_eval`, `constraints_eval`,
`initial_data`, `diagnostic_gfs_set`, and `psi4`.

The black-hole and GRHD examples show the dataflow shape: register initial-data
import/conversion, grid and diagnostic helpers, reference-metric precompute,
Ricci/RHS/determinant enforcement/constraints, Method of Lines step glue,
coordinate and basis transforms, wrapper dispatchers, headers, parser, `main`,
and cleanup.

## Detail

`register_CFunction_rhs_eval` generates the BSSN RHS CFunction. It pulls
non-gauge RHS expressions from `BSSN_RHSs[...]`, gauge RHS expressions from
`BSSN_gauge_RHSs`, and optional constraint damping terms from
`BSSN_constraints[...]`. It builds a sorted local name-to-expression dictionary,
maps each RHS name to the matching `rhs_gfs` gridfunction with
`BHaHGridFunction.access_gf`, and emits an interior `simple_loop` with finite
difference codegen, optional SIMD/CUDA intrinsics, optional finite-difference
helper functions, and upwinding controlled by beta. The generated function
signature changes with `enable_rfm_precompute`: it receives either
`rfmstruct` or coordinate arrays. It can include `RbarDD` gridfunctions,
`T4munu`, Kreiss-Oliger dissipation, curvature-aware KO, Hamiltonian-constraint
damping, and slow-start lapse through code-generation flags and commondata
parameters.

`register_CFunction_Ricci_eval` emits `Ricci_eval` from
`BSSN_quantities[CoordSystem + "_rfm_precompute"].Ricci_exprs`. It always uses
the rfm-precompute expression family, stores either ordinary `RBARDD*GF`
auxiliary gridfunctions or `DIAG_RBARDD*GF` channels for the host-only
diagnostics version, and wraps the interior loop with BHaH kernel/launch code.
CUDA generation is rejected for `GeneralRFM`; `host_only_version=True`
temporarily forces OpenMP generation so CUDA applications can still compute
host-side diagnostic Ricci data as `Ricci_eval_host`.

`register_CFunction_constraints_eval` emits the diagnostics-side Hamiltonian
and momentum-constraint evaluator. It temporarily forces OpenMP, reads
`BSSN_constraints[CoordSystem + "_rfm_precompute_RbarDD_gridfunctions" +
optional "_T4munu"]`, writes `DIAG_HAMILTONIANGF` and `DIAG_MSQUAREDGF`, and
places the function in the `diagnostics/` subdirectory. When the original
parallelization is CUDA, generated references to `RBARDD` and optional `T4UU`
auxiliary gridfunctions are rewritten to diagnostic channels so host-side
constraint evaluation consumes the diagnostic buffer filled for output.

`register_CFunction_enforce_detgammabar_equals_detgammahat` emits the
determinant-enforcement CFunction used after RHS/MoL updates. It reconstructs
`detgammabar` from `BSSN_quantities[...]`, builds corrected `hDD` components so
`det(gammabar)` matches `det(gammahat)`, writes those components back into
`in_gfs`, and runs over all points rather than only interiors. The generated
kernel accepts either `rfmstruct` or coordinate arrays according to
`enable_rfm_precompute`, plus read-only `auxevol_gfs`.

`register_CFunction_initial_data` is the application-level initial-data
assembler. For built-in exact data it instantiates `InitialData_Cartesian` or
`InitialData_Spherical`, registers the exact ADM provider by
`register_CFunction_exact_ADM_ID_function`, then registers per-coordinate
ADM-to-BSSN readers through
`register_CFunctions_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`.
The generated `initial_data()` optionally attempts `read_checkpoint()` first;
on restart it applies inner boundary conditions and optional interpatch
interpolation before returning. Without restart it creates an `ID_persist`
struct, optionally populates it, loops over grids, optionally runs
`generalrfm_precompute`, calls
the selected initial-data reader conversion function, applies
outer-extrapolation plus inner boundary conditions, and optionally frees
persistent initial-data storage.

The ADM reader converter is a generated prefunc chain. `register_BHaH_defines_h`
adds `initial_data_struct` and `ID_persist_struct` to `BHaH_defines.h`.
`Cfunction_ADM_SphorCart_to_Cart` transforms ADM variables from the input
spherical, Cartesian, or GeneralRFM basis to Cartesian. `Cfunction_ADM_Cart_to_BSSN_Cart`
converts Cartesian ADM data to Cartesian BSSN fields. `Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm`
transforms those BSSN tensors/vectors to the destination reference-metric basis
and applies BSSN rescalings. `build_initial_data_conversion_loop` writes
`alpha`, `cf`, `trK`, `hDD`, `aDD`, `vetU`, `betU`, and optional `T4UU` into
MoL gridfunction arrays. `build_lambdaU_zeroing_block` initializes `lambdaU`,
`build_apply_inner_bcs_block` applies parity-sensitive inner boundary
conditions, and `Cfunction_initial_data_lambdaU_grid_interior` computes
`lambdaU` by finite differencing the initialized conformal metric.

`register_CFunction_diagnostic_gfs_set` bridges evolved GR state to diagnostics.
It registers diagnostic gridfunctions, builds parity metadata for them, and
generates `diagnostic_gfs_set(commondata, griddata, diagnostic_gfs)`. Runtime
flow is per grid: compute Ricci into `auxevol_gfs` or diagnostics, evaluate
constraints into diagnostics, optionally compute Psi4 and apply inner boundary
conditions to `DIAG_PSI4_RE/IM`, optionally apply inner boundary conditions to
constraint diagnostics before interpolation, then copy lapse, conformal factor,
and grid index into diagnostic channels. `register_CFunction_diagnostics_nearest`
and `register_CFunction_diagnostics_volume_integration` consume these
diagnostic buffers without owning their memory.

The nearest and volume diagnostic wiring is deliberately generic after
`diagnostic_gfs_set`. `diagnostics_nearest()` exposes user-editable `which_gfs`
arrays and dispatches to 0D, 1D, and 2D nearest samplers. Volume diagnostics
build recipes from diagnostic enum tokens and call
`diags_integration_execute_recipes`. The volume-element helper is
coordinate-specialized by `register_CFunction_sqrt_detgammahat_d3xx_volume_element`
when diagnostics are registered.

Basis transforms are registered through
`basis_transforms.register_all.register_CFunctions`. The two production modules
emit private per-coordinate single-point kernels with coordinate-system suffixes.
The public unsuffixed runtime dispatchers are emitted later by
`rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs`. The
rfm-to-Cartesian path reads rescaled BSSN storage, reconstructs barred
quantities, transforms `gammabarDD`, `AbarDD`, beta, B, and Lambdabar to
Cartesian, and writes Cartesian-basis storage. The Cartesian-to-rfm path reads
Cartesian storage, transforms tensors/vectors into the destination basis,
subtracts the destination reference metric where needed, rescales by `ReDD` and
`ReU`, and writes native BSSN storage.

TwoPunctures is wired as an external compact-object initial-data family.
`TwoPunctures_lib.register_C_functions_explicit` registers
`initialize_ID_persist_struct`, `TP_CoordTransf`, `TP_Equations`,
`TP_FuncAndJacobian`, `TP_Newton`, `TP_Interp`, `TP_solve`, and
`TP_utilities`. `ID_persist_str` registers commondata controls for binary
description, spectral resolution, and optional bare masses, then contributes
the persistent spectral-solve fields. `initialize_ID_persist_struct` sets
defaults, copies or derives binary masses, separation, momenta, spin, center
offset, spectral grid sizes, and orientation. The explicit orientation choices
are `native_cartesian_xy_plane` and `legacy_swap_xz`; the older
`register_C_functions(enable_xy_plane=...)` shim maps the boolean interface to
those names. If TwoPunctures/TOVola setup detail grows beyond routing and
dataflow, split it into a future `compact-object-initial-data.md` leaf.

TOVola is the corresponding single-star initial-data path. `TOVola.ID_persist_str`
registers central density, polytropic EOS constants, ODE controls, and
interpolation stencil sizes, then contributes radial-table pointers and count
fields to `ID_persist_struct`. `register_CFunction_TOVola_solve` registers the
GSL ODE integration driver, stores the solved radial data in the persistent
arrays, and frees temporary solve storage. `register_CFunction_TOVola_interp`
registers the pointwise interpolation provider that maps Cartesian points to
isotropic radius, interpolates the radial table, and fills ADM-like
`initial_data_struct` fields including lapse, spherical spatial metric,
zero extrinsic curvature, and stress-energy components. The GRHD TOV example
uses this by registering `TOVola_interp`, `TOVola_solve`, and a spherical
ADM-to-BSSN reader whose `ID_persist_struct_str` comes from TOVola.

Psi4 wiring has two layers. `psi4.register_CFunction_psi4` temporarily forces
OpenMP generation, requires the grid origin to be zero, calls
`generate_CFunction_psi4_tetrad` and
`generate_CFunction_psi4_metric_deriv_quantities` to create local helper
kernels, then loops over interior points to write `DIAG_PSI4_REGF` and
`DIAG_PSI4_IMGF`. The tetrad helper uses `Psi4Tetrads` to compute
`mre4U`, `mim4U`, and `n4U` from metric gridfunction values; the derivative
helper computes `gammaDDdDD`, `GammaUDD`, and `KDDdD` arrays used by the
symbolic `Psi4` expression. `psi4_spinweightm2_decomposition` then interpolates
diagnostic Psi4 data from grid 0 onto spherical extraction shells, calls
`spin_weight_minus2_sph_harmonics` for each `(l,m)`, numerically integrates
over the shell, and appends `R_ext * psi4_{l,m}` time series to radius- and
mode-tagged text files.

`spin_weight_minus2_sph_harmonics` is registered in the BHaH special-functions
branch. It registers `swm2sh_maximum_l_mode_to_compute` in commondata and emits
a switch over generated `l,m` pairs. Out-of-range requests print an error and
exit, so the decomposition generator must register enough modes for the
requested extraction.

## Sources

- [rhs_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/rhs_eval.py) - `register_CFunction_rhs_eval`
- [Ricci_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/Ricci_eval.py) - `register_CFunction_Ricci_eval`
- [constraints_eval.py](../../../nrpy/infrastructures/BHaH/general_relativity/constraints_eval.py) - `register_CFunction_constraints_eval`
- [enforce_detgammabar_equals_detgammahat.py](../../../nrpy/infrastructures/BHaH/general_relativity/enforce_detgammabar_equals_detgammahat.py) - `register_CFunction_enforce_detgammabar_equals_detgammahat`
- [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - `register_CFunction_exact_ADM_ID_function`, `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`, `Cfunction_ADM_SphorCart_to_Cart`, `Cfunction_ADM_Cart_to_BSSN_Cart`, `Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm`, `Cfunction_initial_data_lambdaU_grid_interior`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics_nearest.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [diagnostics_volume_integration.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_volume_integration.py) - `register_CFunction_diagnostics_volume_integration`
- [register_all.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/register_all.py) - `register_CFunctions`
- [basis_transform_BSSN_rfm_to_Cartesian_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_rfm_to_Cartesian_single_point.py) - `register_CFunction_basis_transform_BSSN_rfm_to_Cartesian_single_point`
- [basis_transform_BSSN_Cartesian_to_rfm_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_Cartesian_to_rfm_single_point.py) - `register_CFunction_basis_transform_BSSN_Cartesian_to_rfm_single_point`
- [TwoPunctures_lib.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TwoPunctures_lib.py) - `register_C_functions_explicit`, `register_C_functions`
- [ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/ID_persist_struct.py) - `ID_persist_str`, `register_CFunction_initialize_ID_persist_struct`
- [TOVola/ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/ID_persist_struct.py) - `ID_persist_str`
- [TOVola_solve.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_solve.py) - `register_CFunction_TOVola_solve`
- [TOVola_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_interp.py) - `register_CFunction_TOVola_interp`
- [psi4.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/psi4.py) - `register_CFunction_psi4`
- [compute_psi4_metric_deriv.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_metric_deriv.py) - `generate_CFunction_psi4_metric_deriv_quantities`
- [compute_psi4_tetrad.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_tetrad.py) - `generate_CFunction_psi4_tetrad`
- [psi4_spinweightm2_decomposition.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4_spinweightm2_decomposition.py) - `register_CFunction_psi4_spinweightm2_decomposition`, `lowlevel_decompose_psi4_into_swm2_modes`
- [spin_weight_minus2_spherical_harmonics.py](../../../nrpy/infrastructures/BHaH/special_functions/spin_weight_minus2_spherical_harmonics.py) - `register_CFunction_spin_weight_minus2_sph_harmonics`
- [two_blackholes_collide.py](../../../nrpy/examples/two_blackholes_collide.py) - `BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval`, `BHaH.general_relativity.basis_transforms.register_all.register_CFunctions`
- [groovy_TOV_BSSN.py](../../../nrpy/examples/groovy_TOV_BSSN.py) - `BHaH.general_relativity.TOVola.TOVola_interp.register_CFunction_TOVola_interp`, `BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`
- [SOURCES.md](../../../raw/SOURCES.md) - `infrastructure-modules-and-embedded-headers`

## See Also

- Parent: [BHaH](index.md)
- Depends on: [BSSN Family](../../equations/general-relativity/bssn-family.md)
- Depends on: [Initial Data](../../equations/general-relativity/initial-data.md)
- Depends on: [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- Depends on: [Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)
- See also: [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)

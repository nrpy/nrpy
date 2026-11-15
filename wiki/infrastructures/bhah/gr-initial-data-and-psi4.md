# GR Initial Data, Diagnostics, And Psi4

> BHaH CFunction registration for ADM-to-BSSN/fCCZ4 initial-data conversion, TwoPunctures, TOVola, basis transforms, diagnostic buffers, and Psi4 mode extraction. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH's initial-data and diagnostic modules register `initial_data`, ADM-to-BSSN/fCCZ4 coordinate and rescaling converters, compact-object solvers (`TwoPunctures` and `TOVola`), curvilinear-to-Cartesian basis transforms, `diagnostic_gfs_set` with nearest/volume samplers, and `psi4` spin-weight `-2` spherical-harmonic extraction.

## Detail

### Initial-Data Assembly, ADM Conversion, And Basis Transforms

`initial_data.py::register_CFunction_initial_data` assembles exact (`InitialData_Cartesian` / `InitialData_Spherical` via `register_CFunction_exact_ADM_ID_function`) or external initial data and per-coordinate ADM converters (`ADM_Initial_Data_Reader__BSSN_Converter.py`). Generated `initial_data()` optionally calls `read_checkpoint()` first, applying inner boundary conditions, optional interpatch interpolation, and optional `enable_conformal_projection` before returning—preserving checkpointed `Theta_fCCZ4` when `enable_fCCZ4=True`. On fresh runs it allocates and populates `ID_persist_struct`, optionally runs `generalrfm_precompute`, converts ADM data on each grid, applies `apply_bcs_outerextrap_and_inner` (and optional conformal projection), and frees `ID_persist_struct`.

`ADM_Initial_Data_Reader__BSSN_Converter.py` emits a `prefunc` chain alongside `initial_data_struct` and `ID_persist_struct` in `BHaH_defines.h`: `Cfunction_ADM_SphorCart_to_Cart` -> `Cfunction_ADM_Cart_to_BSSN_Cart` -> `Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm`, writing `alpha`, `cf`, `trK`, `hDD`, `aDD`, `vetU`, `betU`, optional `T4UU`, and fresh-data `Theta_fCCZ4 = 0` (`enable_fCCZ4=True`; YBS-MOM adds no state). It then zeroes `lambdaU`, applies parity-aware inner BCs, and finite-differences the conformal metric in `Cfunction_initial_data_lambdaU_grid_interior`.

`basis_transforms/register_all.py::register_CFunctions` emits coordinate-suffixed single-point kernels (`rfm_to_Cartesian` reconstructing barred fields and transforming `gammabarDD`, `AbarDD`, `beta`, `B`, and `Lambdabar` to Cartesian; `Cartesian_to_rfm` transforming to the destination basis, subtracting `gammahatDD`, and dividing by `ReDD`/`ReU`), dispatched at runtime by `rfm_wrapper_functions.py`.

Evidence (descriptive behavior): `initial_data.py::register_CFunction_initial_data`; `ADM_Initial_Data_Reader__BSSN_Converter.py`; `basis_transforms/register_all.py`; `blackhole_spectroscopy.py`.

### TwoPunctures And TOVola Solvers

`TwoPunctures/TwoPunctures_lib.py::register_C_functions_explicit` registers `initialize_ID_persist_struct`, `TP_CoordTransf`, `TP_Equations`, `TP_FuncAndJacobian`, `TP_Newton`, `TP_Interp`, `TP_solve`, and `TP_utilities` (`register_C_functions(enable_xy_plane=...)` maps booleans to `native_cartesian_xy_plane` and `legacy_swap_xz`). `ID_persist_struct.py` registers binary mass, separation, momentum, spin, offset, and spectral-grid `commondata` parameters. While the shared initial-lapse default is `psi^n`, `TP_Interp` also supports `"W"` (`alpha = 1 / (psi1 / static_psi)^2`, including spectral correction `U` in `psi1`), selected explicitly by `blackhole_spectroscopy.py` and `superB_blackhole_spectroscopy.py` for slow-start lapse while `superb_paper2` scripts retain `psi^n`.

`TOVola` (`ID_persist_struct.py`, `TOVola_solve.py`, `TOVola_interp.py`) registers polytropic EOS and GSL ODE controls, integrates the radial TOV system into `ID_persist_struct` tables, and interpolates isotropic-radius spherical ADM metric, lapse, zero `K_ij`, and `T^(mu nu)` into `initial_data_struct` (used by `groovy_TOV_BSSN.py`).

Evidence (descriptive behavior): `TwoPunctures/TwoPunctures_lib.py`; `TwoPunctures/ID_persist_struct.py`; `TwoPunctures/TP_interp.py`; `TOVola/TOVola_solve.py`; `TOVola/TOVola_interp.py`.

### Diagnostic Buffers And Psi4 Mode Extraction

`diagnostic_gfs_set.py::register_CFunction_diagnostic_gfs_set` populates `diagnostic_gfs` per grid: it computes Ricci and `constraints_eval`, optionally evaluates `psi4` and applies `apply_bcs_inner_only_specific_gfs` to `DIAG_PSI4_REGF`/`DIAG_PSI4_IMGF` (and to constraint channels when interpolation diagnostics are enabled), and copies `alpha`, `cf`, and `grid_idx`. `diagnostics_nearest.py` (0D/1D/2D `which_gfs` samplers) and `diagnostics_volume_integration.py` (`diags_integration_execute_recipes` with `sqrt_detgammahat_d3xx_volume_element`) read `diagnostic_gfs` without managing its allocation.

`psi4/psi4.py::register_CFunction_psi4` forces OpenMP, generates `compute_psi4_tetrad` (`Psi4Tetrads` `mre4U`, `mim4U`, `n4U`) and `compute_psi4_metric_deriv_quantities` (`gammaDDdDD`, `GammaUDD`, `KDDdD`), and writes interior `DIAG_PSI4_REGF` and `DIAG_PSI4_IMGF`. `psi4_spinweightm2_decomposition.py` interpolates grid `0` onto origin-centered extraction spheres, projects against `spin_weight_minus2_sph_harmonics` (`swm2sh_maximum_l_mode_to_compute`, inclusive `m in [-l, +l]` in C), and appends `R_ext * psi4_{l,m}`.
Because the tetrad uses grid-local coordinates while shells center on `(0, 0, 0)`, `Cart_origin{x,y,z}` must be zero (unchecked in C; satisfied in BHaH/superB spectroscopy where `Cart_origin*` default to `0` and are excluded from parfile parsing). Shell interpolation runs only when every shell point's unclamped nearest index lies in the grid interior (writing `0` otherwise) and reads inner-BC-filled ghost zones without outer-BC fills; see [Inner Ghost Points In Interpolation Stencils](curvilinear-boundaries-and-stencils.md#inner-ghost-points-in-interpolation-stencils).

Evidence (descriptive behavior): `diagnostic_gfs_set.py`; `psi4/psi4.py`; `BHaH/general_relativity/psi4_spinweightm2_decomposition.py`; `spin_weight_minus2_spherical_harmonics.py`; `grid.py::Cart_origin`.

## Sources

- [initial_data.py](../../../nrpy/infrastructures/BHaH/general_relativity/initial_data.py) - `register_CFunction_initial_data`
- [ADM_Initial_Data_Reader__BSSN_Converter.py](../../../nrpy/infrastructures/BHaH/general_relativity/ADM_Initial_Data_Reader__BSSN_Converter.py) - `register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN`, `Cfunction_initial_data_lambdaU_grid_interior`
- [diagnostic_gfs_set.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostic_gfs_set.py) - `register_CFunction_diagnostic_gfs_set`
- [diagnostics_nearest.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_nearest.py) - `register_CFunction_diagnostics_nearest`
- [diagnostics_volume_integration.py](../../../nrpy/infrastructures/BHaH/general_relativity/diagnostics_volume_integration.py) - `register_CFunction_diagnostics_volume_integration`
- [register_all.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/register_all.py) - `register_CFunctions`
- [basis_transform_BSSN_rfm_to_Cartesian_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_rfm_to_Cartesian_single_point.py) - `register_CFunction_basis_transform_BSSN_rfm_to_Cartesian_single_point`
- [basis_transform_BSSN_Cartesian_to_rfm_single_point.py](../../../nrpy/infrastructures/BHaH/general_relativity/basis_transforms/basis_transform_BSSN_Cartesian_to_rfm_single_point.py) - `register_CFunction_basis_transform_BSSN_Cartesian_to_rfm_single_point`
- [TwoPunctures_lib.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TwoPunctures_lib.py) - `register_C_functions_explicit`, `register_C_functions`
- [ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/ID_persist_struct.py) - `ID_persist_str`, `register_CFunction_initialize_ID_persist_struct`
- [TP_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TwoPunctures/TP_interp.py) - `register_CFunction_TP_Interp`
- [TOVola/ID_persist_struct.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/ID_persist_struct.py) - `ID_persist_str`
- [TOVola_solve.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_solve.py) - `register_CFunction_TOVola_solve`
- [TOVola_interp.py](../../../nrpy/infrastructures/BHaH/general_relativity/TOVola/TOVola_interp.py) - `register_CFunction_TOVola_interp`
- [psi4.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/psi4.py) - `register_CFunction_psi4`
- [compute_psi4_metric_deriv.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_metric_deriv.py) - `generate_CFunction_psi4_metric_deriv_quantities`
- [compute_psi4_tetrad.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4/compute_psi4_tetrad.py) - `generate_CFunction_psi4_tetrad`
- [psi4_spinweightm2_decomposition.py](../../../nrpy/infrastructures/BHaH/general_relativity/psi4_spinweightm2_decomposition.py) - `register_CFunction_psi4_spinweightm2_decomposition`
- [superB psi4_spinweightm2_decomposition.py](../../../nrpy/infrastructures/superB/general_relativity/psi4_spinweightm2_decomposition.py) - `register_CFunction_psi4_spinweightm2_decomposition`
- [spin_weight_minus2_spherical_harmonics.py](../../../nrpy/infrastructures/BHaH/special_functions/spin_weight_minus2_spherical_harmonics.py) - `register_CFunction_spin_weight_minus2_sph_harmonics`
- [grid.py](../../../nrpy/grid.py) - `Cart_origin`
- [blackhole_spectroscopy.py](../../../nrpy/examples/blackhole_spectroscopy.py) - `--fccz4` initial-data and TwoPunctures `"W"` lapse wiring
- [superB_blackhole_spectroscopy.py](../../../nrpy/examples/superB_blackhole_spectroscopy.py) - TwoPunctures `"W"` lapse wiring
- [superB_blackhole_spectroscopy_8Mseparation.py](../../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_8Mseparation.py) - default `psi^n` TwoPunctures lapse
- [superB_blackhole_spectroscopy_last_orbit.py](../../../nrpy/examples/superb_paper2/superB_blackhole_spectroscopy_last_orbit.py) - default `psi^n` TwoPunctures lapse
- [groovy_TOV_BSSN.py](../../../nrpy/examples/groovy_TOV_BSSN.py) - TOVola + spherical ADM-to-BSSN wiring

## See Also

- Depends on: [Initial Data](../../equations/general-relativity/initial-data.md)
- Depends on: [Metric Conversions And Matter](../../equations/general-relativity/metric-conversions-and-matter.md)
- Depends on: [Psi4 And Tetrads](../../equations/general-relativity/psi4-and-tetrads.md)

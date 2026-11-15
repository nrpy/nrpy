# Grids, Coordinates, And Reference Metrics

> BHaH CFunction registration for numerical grid parameters, CFL timesteps, `rfm_precompute`, `GeneralRFM` fisheye maps, coordinate converters, and `CoordSystem_hash` wrappers. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH grid initialization registers per-grid coordinate parameters, allocates cell-centered `xx[3]` arrays for independent or multipatch grids, computes CFL-limited timesteps, populates ordinary `rfm_struct` or `GeneralRFM` `AUXEVOL` reference-metric data, and registers `xx_to_Cart` / `Cart_to_xx` converters with runtime `CoordSystem_hash` dispatchers.

## Detail

### Grid Parameters, Multi-Grid Topologies, And CFL Timestep

`numerical_grids_and_timestep.py::register_CFunctions` registers per-grid parameters `Nxx*`, `Nxx_plus_2NGHOSTS*`, `xxmin*`, `xxmax*`, `dxx*`, `invdxx*`, `convergence_factor`, `CoordSystem_hash`, `grid_idx`, and `gridname`. `numerical_grid_params_Nxx_dxx_xx()` initializes default `Nxx` from `Nxx_dict`, applies optional runtime `Nx[]` overrides, multiplies counts by `commondata->convergence_factor` except along symmetry axes with `Nxx == 2`, evaluates `xxmin*`/`xxmax*` from reference-metric parameters, computes `dxx*` and `invdxx*`, and allocates cell-centered `xx[3]` arrays (allocating device `xx` and copying `d_params` on CUDA).

`numerical_grids_and_timestep()` handles `"independent grid(s)"` directly (setting `gridname`, `CoordSystem_hash`, `grid_physical_size`, per-coordinate initializers, `NUMGRIDS`, `grid_idx`, and initial `nn`, `nn_0`, `t_0`, `time`) and delegates `"multipatch"` to `multipatch_grids_set_up()`.
When CFL calculation is enabled, `cfl_limited_timestep()` sets `commondata->dt = min(commondata->dt, ds_min * commondata->CFL_FACTOR)` across active grids by evaluating `ds_min_single_pt()` from reference-metric scale factors on standard coordinate systems or from the fisheye metric diagonal on `GeneralRFM_fisheyeN*` (non-fisheye `GeneralRFM` reports `ds_min` unsupported). Multipatch setups (`gridding_approach == "multipatch"`) also register `ds_min_single_pt()` alongside `ds_min_radial_like_dirns_single_pt()` (`register_CFunction_ds_min_radial_like_dirns_single_pt`), which evaluates minimum spacing across `rfm.radial_like_dirns` (non-angular coordinate axes: `[0]` spherical/wedge, `[0, 1]` SymTP/prolate-spheroidal, `[0, 2]` cylindrical, `[0, 1, 2]` Cartesian and `GeneralRFM_fisheyeN*`) for regrid and masking support.

Evidence (descriptive behavior): `numerical_grids_and_timestep.py::register_CFunctions`; `numerical_grids_and_timestep.py::register_CFunction_cfl_limited_timestep`; `numerical_grids_and_timestep.py::register_CFunction_ds_min_radial_like_dirns_single_pt`.

### Reference-Metric Precompute And Fisheye Parameters

Standard non-`GeneralRFM` precompute (`rfm_precompute.py`) identifies coordinate-dependent reference-metric subexpressions, defines `rfm_struct`, and registers host/CUDA `rfm_precompute_malloc`, `rfm_precompute_defines`, and `rfm_precompute_free`.

For `GeneralRFM_fisheyeN*`, `generalrfm_precompute.py` stores `ghatDD`, `ghatUU`, `detgammahat`, `ghatDDdD`, and `ghatDDdDD` in `AUXEVOL` gridfunctions, rejecting CUDA and unsupported providers. `generalrfm_cart_to_xx.py` registers a coordinate-specialized `generalrfm_Cart_to_xx` radial Newton solver returning nonzero on non-convergence. `fisheye/phys_params_to_fisheye.py` registers physical fisheye `CodeParameter` entries in `commondata_struct`, computes internal `fisheye_a*`, `fisheye_R*`, `fisheye_s*`, and `fisheye_c`, and emits a post-`params_struct_set_to_default()` conversion hook.

Evidence (descriptive behavior): `rfm_precompute.py::ReferenceMetricPrecompute`; `generalrfm_precompute.py`; `generalrfm_cart_to_xx.py`; `fisheye/phys_params_to_fisheye.py`.

### Coordinate Converters And Hash Dispatch Wrappers

`xx_tofrom_Cart.py` registers `register_CFunction_xx_to_Cart` and `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid` (`GeneralRFM` CUDA is unsupported). On independent grids, the inverse subtracts `Cart_origin*` and computes `xx`; if `Cart_to_i0i1i2 == NULL` it returns before float-to-int conversion, otherwise it computes nearest indices from `xxmin`, `dxx`, and `NGHOSTS` assuming an in-domain point without coordinate bounds checks (which remain prohibited unless explicitly requested; see [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)). On multipatch grids with `params->grid_rotates`, inverse conversion applies cumulative transpose `R^T` before origin/coordinate inversion and forward `xx_to_Cart` applies `R` afterward; each multipatch registrar registers only its duplicate-safe SO(3) matrix/`R` or `R^T` helper closure, while independent grids register none.

Full-source doctest baselines in `xx_tofrom_Cart.py` cover default independent converters and Cartesian multipatch OpenMP/CUDA variants, whereas OpenMP `GeneralRFM` fisheye inverse (`_run_generalrfm_fisheye_inverse_roundtrip_check`) compiles `REAL=float` and `REAL=double` C checks for finite round-trips and non-finite input failure. `rfm_wrapper_functions.py` emits unsuffixed dispatchers switching on `params->CoordSystem_hash` and writes uppercase coordinate hash macros to `BHaH_defines.h`.

Evidence (generated evidence): `xx_tofrom_Cart.py::register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`; `xx_tofrom_Cart.py::_run_generalrfm_fisheye_inverse_roundtrip_check`; `rfm_wrapper_functions.py`.

## Sources

- [coding_style.md](../../../coding_style.md) - `## Coordinate Bounds-Check Prohibition`
- [reference_metric.py](../../../nrpy/reference_metric.py) - `ReferenceMetric.radial_like_dirns`
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `register_CFunctions`, `register_CFunction_numerical_grid_params_Nxx_dxx_xx`, `register_CFunction_numerical_grids_and_timestep`, `register_CFunction_cfl_limited_timestep`, `register_CFunction_ds_min_single_pt`, `register_CFunction_ds_min_radial_like_dirns_single_pt`
- [rfm_precompute.py](../../../nrpy/infrastructures/BHaH/rfm_precompute.py) - `ReferenceMetricPrecompute`, `register_CFunctions_rfm_precompute`
- [generalrfm_precompute.py](../../../nrpy/infrastructures/BHaH/generalrfm_precompute.py) - `register_CFunction_generalrfm_precompute`, `register_CFunctions_generalrfm_support`
- [generalrfm_cart_to_xx.py](../../../nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py) - `register_CFunction_generalrfm_Cart_to_xx`
- [phys_params_to_fisheye.py](../../../nrpy/infrastructures/BHaH/fisheye/phys_params_to_fisheye.py) - `register_CFunction_fisheye_params_from_physical_N`, `build_post_params_struct_set_to_default_hook`
- [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py) - `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, `register_CFunction_xx_to_Cart`, `_run_generalrfm_fisheye_inverse_roundtrip_check`
- [rfm_wrapper_functions.py](../../../nrpy/infrastructures/BHaH/rfm_wrapper_functions.py) - `get_CoordSystem_hash`, `register_CFunctions_CoordSystem_wrapper_funcs`

## See Also

- Depends on: [Reference Metrics](../../core/reference-metrics.md)
- Depends on: [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)
- See also: [Curvilinear Boundaries And Stencils](curvilinear-boundaries-and-stencils.md)

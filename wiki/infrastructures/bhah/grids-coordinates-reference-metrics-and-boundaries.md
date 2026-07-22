# Grids, Coordinates, Reference Metrics, And Boundaries

> BHaH route for numerical grid setup, coordinate wrappers, reference-metric precompute, GeneralRFM, fisheye hooks, and curvilinear boundaries. Status: confirmed. Last reconciled: 07-20-2026
> Up: [BHaH](index.md)

## Summary

BHaH grid setup starts by registering per-grid numerical parameters, allocating
cell-centered coordinate arrays, and selecting either independent grids or a
multipatch setup. It then optionally populates reference-metric precompute data,
sets curvilinear boundary metadata, chooses a CFL-limited timestep, and registers
coordinate-system wrappers so generated code can dispatch by
`params->CoordSystem_hash`.

## Detail

`register_CFunctions()` in `numerical_grids_and_timestep.py` registers the core
grid parameters: `Nxx*`, `Nxx_plus_2NGHOSTS*`, `xxmin*`, `xxmax*`, `dxx*`,
`invdxx*`, `convergence_factor`, `CoordSystem_hash`, `grid_idx`, and
`gridname`. `numerical_grid_params_Nxx_dxx_xx()` sets default `Nxx` from the
generator's `Nxx_dict`, allows a runtime `Nx[]` override, applies
`commondata->convergence_factor` except along symmetry axes with `Nxx == 2`,
sets `Nxx_plus_2NGHOSTS`, evaluates coordinate extents from reference-metric
parameters, computes `dxx` and `invdxx`, and allocates cell-centered `xx[3]`
arrays. CUDA setup allocates these arrays on device and copies per-grid params
to `d_params`.

`numerical_grids_and_timestep()` supports `"independent grid(s)"` directly and
delegates `"multipatch"` to `multipatch_grids_set_up()`. Independent-grid setup
assigns `gridname`, `CoordSystem_hash`, `grid_physical_size`, calls the
coordinate-grid initializer for each coordinate system, records `NUMGRIDS`, and
sets `grid_idx`. First-time setup initializes `nn`, `nn_0`, `t_0`, and `time`.
When CFL setup is enabled, `cfl_limited_timestep()` computes
`commondata->dt = min(commondata->dt, ds_min * commondata->CFL_FACTOR)` across
active grids. For ordinary coordinate systems, `ds_min_single_pt()` uses
reference-metric scale factors; for `GeneralRFM_fisheyeN*`, it uses the fisheye
metric diagonal. Non-fisheye `GeneralRFM` reports unsupported `ds_min`.

Reference-metric precompute has two paths. `rfm_precompute.py` handles ordinary
non-`GeneralRFM` precompute by discovering coordinate-dependent expressions,
creating `rfm_struct`, and registering `rfm_precompute_malloc`,
`rfm_precompute_defines`, and `rfm_precompute_free`; generated host and CUDA
paths allocate, populate, read, and free lookup arrays. `generalrfm_precompute.py`
handles `GeneralRFM_fisheyeN*` by storing `ghatDD`, `ghatUU`, `detgammahat`,
`ghatDDdD`, and `ghatDDdDD` into `AUXEVOL` gridfunctions. That GeneralRFM
precompute path rejects CUDA and unsupported providers. `generalrfm_cart_to_xx.py`
registers a coordinate-specialized `generalrfm_Cart_to_xx` Newton solve that
inverts a supported GeneralRFM fisheye map and returns nonzero on failure.
`phys_params_to_fisheye.py` registers physical fisheye CodeParameters in
`commondata_struct`, computes internal `fisheye_a*`, `fisheye_R*`,
`fisheye_s*`, and `fisheye_c`, and provides a post-params hook that runs the
conversion after `params_struct_set_to_default()`.

Coordinate wrappers live in two layers. `xx_tofrom_Cart.py` registers the Python
registrars `register_CFunction_xx_to_Cart` and
`register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, which emit the
coordinate-specific C functions `xx_to_Cart` and
`Cart_to_xx_and_nearest_i0i1i2_assume_valid`. For an independent grid, the
inverse unshifts Cartesian points by `Cart_origin*` for local-grid calculations
and maps them to `xx`. When the caller supplies `Cart_to_i0i1i2`, the emitted C
function then converts `xx` to nearest indices with `xxmin`, `dxx`, and
`NGHOSTS`; this conversion assumes a valid, in-domain, cell-centered point and
does not perform a bounds check. When `Cart_to_i0i1i2 == NULL`, it returns after
computing `xx`, before any floating-to-integer index conversion. For a
multipatch grid with `params->grid_rotates`, inverse conversion builds the
cumulative rotation matrix and applies `R^T` to the Cartesian vector before
local-origin handling and native-coordinate inversion; forward `xx_to_Cart`
applies `R` after native-to-Cartesian mapping and origin handling. GeneralRFM
fisheye wrappers use a radial Newton solve for inverse mapping and a closed-form
fisheye radius map for forward mapping. Each multipatch registrar also registers
its exact SO(3) dependency closure: inverse registers the matrix builder and
`R^T` vector helper, while forward registers the matrix builder and `R` vector
helper. Registration remains duplicate-safe in either converter order and
independent-grid registration adds no SO(3) helper.

Maintenance rule: coordinate-admission checks are request-gated by
[Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md);
the nullable coordinates-only path does not authorize adding them.

The 30 current `Cart_to_xx_and_nearest_i0i1i2_assume_valid` and 30 current
`xx_to_Cart` default independent-grid baselines and four explicitly named
Cartesian multipatch OpenMP/CUDA baselines passed isolated candidate
generation, review, byte comparison, and a second fresh-process source comparison. No C/CUDA
compilation, runtime inverse check, or numerical-result guarantee was
established. Then
`rfm_wrapper_functions.py` creates non-coordinate-specific wrapper functions that
switch on `params->CoordSystem_hash`, calls the matching coordinate-specific
function, and registers uppercase coordinate hash macros in `BHaH_defines.h`.

Claim evidence:
- Claim: `Cart_to_xx_and_nearest_i0i1i2_assume_valid` returns logical coordinates without index conversion when `Cart_to_i0i1i2 == NULL`; each multipatch converter self-registers its exact SO(3) dependency closure; and rotating-multipatch inverse/forward conversion applies `R^T` and `R`, respectively.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/BHaH/xx_tofrom_Cart.py` - `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, `register_CFunction_xx_to_Cart`
- Corroboration: none available; owner-derived emitted-source comparisons are not independent evidence.
- Validation: `inspected=pass; generated=pass; built=not-run; run=not-run; result_checked=pass`
- Dimensions: `platform=Linux; tool_version=Python 3.12.3, clang-format 22.1.8; backend=OpenMP C and CUDA source; precision=not-applicable; GPU=not-run; restart=not-applicable; distributed=not-applicable; error_path=not-run; options=15 default independent coordinate systems plus Cartesian rotating multipatch; date=07-20-2026`

Curvilinear boundary registration starts in
`CurviBoundaryConditions/register_all.py`. It registers `outer_bc_type` with
default `"radiation"`, the coordinate-specific `bcstruct_set_up()`, radiation
outer plus inner BCs, inner-only BCs, specific-gridfunction inner-only BCs,
extrapolation outer plus inner BCs, `bc_struct`'s `griddata_struct`
contribution, and BHaH header definitions. CUDA builds also register
`cpyHosttoDevice_bc_struct()`.

The boundary header defines `innerpt_bc_struct`, `outerpt_bc_struct`,
`bc_info_struct`, and `bc_struct`. `bc_struct` owns `inner_bc_array`,
`pure_outer_bc_array[NGHOSTS*3]`, and `bc_info`. Parity arrays are inferred from
registered gridfunction names for evolved gridfunctions, optionally auxiliary
and auxevol gridfunctions; high-rank GeneralRFM auxevol metric fields receive a
placeholder parity only to preserve table indexing.

`bcstruct_set_up()` classifies boundary points by mapping each ghost-zone point
through the coordinate map to an inbounds point. Points that map away from
themselves become inner boundary points with destination index, source index,
and ten parity values. Points that map to themselves become pure outer boundary
points, stored by ghost-zone layer and direction with face signs. Optional masks
skip unset points and identify outer-boundary points; CUDA masking is rejected
by numerical-grid setup.

Boundary application has three main variants. `apply_bcs_inner_only()` uses
`bcstruct->inner_bc_array` and parity tables to copy from source to destination
inner boundary points. `apply_bcs_outerradiation_and_inner()` fills pure outer
points first with radiation BCs using wavespeed and `f_infinity`, then applies
inner BCs to RHS gridfunctions. `apply_bcs_outerextrap_and_inner()` fills pure
outer points by second-order extrapolation, then applies inner BCs; a
specific-gridfunction variant performs the same sequence for selected grid
functions.

## Sources

- [coding_style.md](../../../coding_style.md) - `## Coordinate Bounds-Check Prohibition`
- [numerical_grids_and_timestep.py](../../../nrpy/infrastructures/BHaH/numerical_grids_and_timestep.py) - `register_CFunctions`, `register_CFunction_numerical_grid_params_Nxx_dxx_xx`, `register_CFunction_numerical_grids_and_timestep`, `register_CFunction_cfl_limited_timestep`, `register_CFunction_ds_min_single_pt`, `register_CFunction_ds_min_radial_like_dirns_single_pt`
- [rfm_precompute.py](../../../nrpy/infrastructures/BHaH/rfm_precompute.py) - `ReferenceMetricPrecompute`, `register_CFunctions_rfm_precompute`
- [rfm_wrapper_functions.py](../../../nrpy/infrastructures/BHaH/rfm_wrapper_functions.py) - `get_CoordSystem_hash`, `register_CFunctions_CoordSystem_wrapper_funcs`
- [xx_tofrom_Cart.py](../../../nrpy/infrastructures/BHaH/xx_tofrom_Cart.py) - `register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid`, `register_CFunction_xx_to_Cart`
- [generalrfm_precompute.py](../../../nrpy/infrastructures/BHaH/generalrfm_precompute.py) - `register_CFunction_generalrfm_precompute`, `register_CFunctions_generalrfm_support`
- [generalrfm_cart_to_xx.py](../../../nrpy/infrastructures/BHaH/generalrfm_cart_to_xx.py) - `register_CFunction_generalrfm_Cart_to_xx`
- [phys_params_to_fisheye.py](../../../nrpy/infrastructures/BHaH/fisheye/phys_params_to_fisheye.py) - `register_CFunction_fisheye_params_from_physical_N`, `build_post_params_struct_set_to_default_hook`
- [register_all.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py) - `register_C_functions`
- [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`
- [BHaH_defines.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/BHaH_defines.py) - `register_BHaH_defines_h`, `BHaH_defines_set_gridfunction_defines_with_parity_types`, `register_griddata_commondata`
- [apply_bcs_outerradiation_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py) - `register_CFunction_apply_bcs_outerradiation_and_inner`, `setup_Cfunction_radiation_bcs`
- [apply_bcs_outerextrap_and_inner.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerextrap_and_inner.py) - `register_CFunction_apply_bcs_outerextrap_and_inner`, `register_CFunction_apply_bcs_outerextrap_and_inner_specific_gfs`
- [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`, `generate_apply_bcs_inner_only__kernel_body`

## See Also

- [BHaH](index.md)
- Depends on: [Contribution Style And Static Analysis](../../architecture/contribution-style-and-static-analysis.md)
- [Runtime Data, Parameters, Headers, And CLI](runtime-data-parameters-headers-and-cli.md)
- [Reference Metrics](../../core/reference-metrics.md)

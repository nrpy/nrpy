# BHaHAHA Grid, Interpolation, And Boundaries

> Explain how BHaHAHA constructs spherical grids, interpolates metric data, and applies local eigencoordinate boundary parity. · Status: confirmed
> Up: [BHaHAHA](index.md)

## Summary

In generated BHaH evolutions, `BHaHAHA_interpolate_metric_data_nrpy` samples BSSN gridfunctions from the evolution grid, transforms them to Cartesian basis, and writes the `12` Cartesian ADM components (`gamma_{ij}`, `K_{ij}`) expected by `BHaHAHA_header.h`. Inside `find_horizon`, BHaHAHA maps that Cartesian input across three spherical grids: the external input grid, the interpolation-source grid (carrying metric derivatives), and the 2D evolved horizon grid (`h(theta, phi)` at one active radial index).

## Detail

### BHaH Evolution-Grid Sampling And Radial Grid Setup

`BHaHAHA_interpolate_metric_data_nrpy` samples `aDD`, `cf`, `hDD`, and `trK` using private `INTERP_BSSN_GF_INDICES` and writes `FINAL_ADM_METRIC_INDICES` matching `INTERP_GAMMADDXXGF`..`INTERP_KDDZZGF`. Destination points on a cell-centered `(r, theta, phi)` sphere around the horizon center pass through `Cart_to_xx_and_nearest_i0i1i2_assume_valid` and `interpolation_3d_general__uniform_src_grid`. The source grid must obey [angular inner-boundary requirements](../curvilinear-boundaries-and-stencils.md#hard-requirements-for-angular-inner-boundaries), and stencil centering follows [inner ghost points in interpolation stencils](../curvilinear-boundaries-and-stencils.md#inner-ghost-points-in-interpolation-stencils); `BHaHAHA_interpolate_metric_data_nrpy` does not check the interpolator's return status, leaving out-of-bounds destinations unset. Sampled points convert via `basis_transform_BSSN_rfm_to_Cartesian_single_point` and `BHaHAHA_BSSN_to_ADM_Cartesian`.

`bah_radial_grid_cell_centered_set_up` clamps `input_r_min >= 0` and `input_r_max <= max_search_radius`, enforces at least `4` interior radial points plus ghost zones, outputs `(Nr_external_input, r_min_external_input, dr_external_input, radii)`, and falls back to `r_min_external_input = 0` whenever ghost zones would cross `r = 0`.

### External, Interpolation-Source, And Evolved Spherical Grids

1. **External input grid (`bah_numgrid__external_input_set_up`)**: Allocates `external_input_gfs` at the finest multigrid `(Ntheta, Nphi)` resolution, copies caller Cartesian ADM data into ghosted storage (shifting active radial indices by `NGHOSTS` when `r_min_external_input == 0`), and converts to rescaled spherical fields `external_spherical_{WW,trK,hDD*,aDD*}`.
2. **Interpolation-source grid (`bah_numgrid__interp_src_set_up`)**: Shares the external radial grid at the current multigrid angular resolution, runs 2D Lagrange angular interpolation (`bah_interpolation_2d_external_input_to_interp_src_grid`) into `src_*` gridfunctions, evaluates interior derivatives (`bah_hDD_dD_and_W_dD_in_interp_src_grid_interior` for `partial_D_hDD` and `partial_D_WW`), and fills radial ghost derivatives via `bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding` (outer boundaries always; inner boundaries only when `r_min_external_input != 0`).
3. **Evolved horizon surface grid (`bah_numgrid__evol_set_up`)**: Sets `NUMGRIDS = 1`, `Nxx0` from evolved shape (`1` active radial zone at `NGHOSTS`), `theta in [0, pi]`, `phi in [-pi, pi]`, evolved `(hh, vv)`, auxiliary metric/derivative fields, and local `bcstruct`.

### Radial-Spoke Interpolation And Local Parity Boundaries

`bah_interpolation_1d_radial_spokes_on_3d_src_grid` interpolates metric and derivative fields along radial spokes onto `r = hh(theta, phi)` at radial index `NGHOSTS`, returning error codes for null pointers, excessive order, `hh > r_max`, or `hh < r_min_external_input`. `bah_interpolation_2d_general__uniform_src_grid` handles arbitrary `(theta, phi)` transfers such as coarse-to-fine prolongation; all 1D/2D/3D interpolators validate pointers, order vs grid size, and stencil bounds.

`bah_bcstruct_set_up` (`bcstruct_set_up.py`, which rejects `GeneralRFM` coordinates at codegen time) emits BHaHAHA-local `innerpt_bc_struct`, `outerpt_bc_struct`, `bc_info_struct`, and `bc_struct`, mapping inner ghost points via spherical eigencoordinates with `10` parity signs plus a signed-permutation derivative Jacobian (returning a boundary error if unrepresentable). `bah_apply_bcs_inner_only` applies inner parity conditions on the evolved grid at radial index `NGHOSTS`, while external and interpolation-source setups apply `bcstruct` inner mappings across all radial shells.

## Sources

- [BHaHAHA_header.h](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `NUM_EXT_INPUT_CARTESIAN_GFS`, `INTERP_GAMMADDXXGF`, `INTERP_KDDZZGF`
- [BHaH_implementation.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `FINAL_ADM_METRIC_INDICES`, `INTERP_BSSN_GF_INDICES`, `BHaHAHA_BSSN_to_ADM_Cartesian`, `BHaHAHA_interpolate_metric_data_nrpy`
- [radial_grid_cell_centered_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/radial_grid_cell_centered_set_up.py) - `register_CFunction_radial_grid_cell_centered_set_up`
- [numgrid__external_input_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__external_input_set_up.py) - `register_CFunction_numgrid__external_input_set_up`
- [numgrid__interp_src_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__interp_src_set_up.py) - `register_CFunction_numgrid__interp_src_set_up`
- [numgrid__evol_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__evol_set_up.py) - `register_CFunction_numgrid__evol_set_up`
- [interpolation_3d_general__uniform_src_grid.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py) - `register_CFunction_interpolation_3d_general__uniform_src_grid`
- [interpolation_2d_external_input_to_interp_src_grid.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_external_input_to_interp_src_grid.py) - `register_CFunction_interpolation_2d_external_input_to_interp_src_grid`
- [interpolation_2d_general__uniform_src_grid.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_general__uniform_src_grid.py) - `register_CFunction_interpolation_2d_general__uniform_src_grid`
- [interpolation_1d_radial_spokes_on_3d_src_grid.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_1d_radial_spokes_on_3d_src_grid.py) - `register_CFunction_interpolation_1d_radial_spokes_on_3d_src_grid`
- [hDD_dD_and_W_dD_in_interp_src_grid_interior.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/hDD_dD_and_W_dD_in_interp_src_grid_interior.py) - `register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior`
- [bcstruct_set_up.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`, `BHaH_defines_set_gridfunction_defines_with_parity_types`
- [apply_bcs_inner_only.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`
- [apply_bcs_r_maxmin_partial_r_hDD_upwinding.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_r_maxmin_partial_r_hDD_upwinding.py) - `register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding`

## See Also

- Depends on: [BHaHAHA Public C Header And Input Fields](bhahaha-public-api-and-input-contract.md)
- Depends on: [Curvilinear Boundaries And Stencils](../curvilinear-boundaries-and-stencils.md)
- Depends on: [Reference Metrics](../../../core/reference-metrics.md)

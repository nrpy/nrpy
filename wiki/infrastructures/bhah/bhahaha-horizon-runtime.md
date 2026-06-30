# BHaHAHA Horizon Runtime

> Explain how the BHaH BHaHAHA runtime registers, feeds, solves, diagnoses, and cleans up apparent-horizon searches. · Status: confirmed · Last reconciled: 2026-06-29
> Up: [BHaH](index.md)

## Summary

BHaHAHA is the BHaH runtime leaf for apparent-horizon finding. The generated
application registers per-horizon commondata arrays, builds a spherical metric
input grid around each candidate horizon, converts interpolated BSSN data into
the ADM Cartesian metric components required by the public BHaHAHA C interface,
then calls `bah_find_horizon` for each active horizon.

The runtime has two layers. `bhahaha_find_horizons` is the NRPy-side
orchestrator that runs inside a BHaH evolution and manages multiple horizons,
BBH/common-horizon state, BSSN-to-ADM interpolation, diagnostics output, and
persistent previous-shape storage. `find_horizon` is the single-horizon
hyperbolic-relaxation driver that owns BHaHAHA's multiresolution grid cycle,
boundary setup, RHS/KO evaluation, over-relaxation, convergence checks, final
diagnostics, and temporary memory cleanup.

## Detail

`register_CFunction_bhahaha_find_horizons` registers two commondata arrays:
`bhahaha_params_and_data[max_horizons]` and
`bhahaha_diagnostics[max_horizons]`. It also registers the public control
parameters used by the generated code: per-horizon search radii, mass scales,
CFL factors, damping strengths, convergence tolerances, KO strength, initial
grid centers, interpolation radial capacity, verbosity, maximum iterations,
multigrid angular-resolution arrays, and BBH/common-horizon switches.

At runtime, `bhahaha_find_horizons` exits when `bah_max_num_horizons <= 0`,
validates `bah_Ntheta_array_multigrid` and `bah_Nphi_array_multigrid`, and
initializes persistent state on the first call. For each horizon it sets
initial center history, radius history, failed-find timestamps, and
`use_fixed_radius_guess_on_full_sphere`, then
`initialize_bhahaha_solver_params_and_shapes` copies commondata controls into
the horizon's `bhahaha_params_and_data_struct`. On the first call it allocates
`prev_horizon_m1`, `prev_horizon_m2`, and `prev_horizon_m3` at the finest
configured angular resolution; `free_bhahaha_horizon_shape_data_all_horizons`
later frees those arrays.

BBH mode requires `bah_max_num_horizons == 3`. The two indices in
`bah_BBH_mode_inspiral_BH_idxs` start active and
`bah_BBH_mode_common_horizon_idx` starts inactive. Once both individual
horizons have recent successful finds, the orchestrator tests whether center
separation plus both stored maximum radii fits inside twice the common
horizon's configured maximum search radius. If so, it activates the common
horizon, seeds its center with a mass-weighted average of the two individual
centers, resets its history to force a full-sphere fixed-radius guess, and
uses the common horizon's `bah_max_search_radius`.

Before each solve, `bah_xyz_center_r_minmax` extrapolates the next center and
radial search interval from the stored `m1/m2/m3` centers, radii, and times.
Unreliable history reduces the horizon CFL factor, forces a full-sphere guess,
and, for an initial BBH common-horizon attempt, promotes the first multigrid
level to the second configured angular resolution and doubles the maximum
iteration count. `bah_radial_grid_cell_centered_set_up` then clamps the search
interval to the configured maximum radius, enforces a minimum number of
interior radial points, includes radial ghost zones, and fills the radial
sampling array passed to interpolation.

`BHaHAHA_interpolate_metric_data_nrpy` is the BHaH-to-BHaHAHA metric bridge. It
builds destination spherical points around the current center guess, maps those
Cartesian locations to the source grid with `Cart_to_xx_and_nearest_i0i1i2`,
interpolates the BSSN gridfunctions listed in `bhahaha_gf_interp_indices` with
`interpolation_3d_general__uniform_src_grid`, transforms each point from the
source reference-metric basis to Cartesian with
`basis_transform_BSSN_rfm_to_Cartesian_single_point`, and writes the final ADM
metric components using `BHaHAHA_BSSN_to_ADM_Cartesian`. The output order is
the private `FINAL_ADM_METRIC_INDICES` enum; the stable public input order is
the `NUM_EXT_INPUT_CARTESIAN_GFS`/`INTERP_*` layout in `BHaHAHA_header.h`.

`find_horizon` is the single-horizon solver driver behind public
`bah_find_horizon`. It creates local `commondata`, copies in the caller-owned
`bhahaha_params_and_data_struct` and `bhahaha_diagnostics_struct`, sets damping,
CFL, KO strength, and diagnostic counters, then calls
`bah_numgrid__external_input_set_up`. That setup allocates
`external_input_gfs`, creates uniform cell-centered `(r, theta, phi)` arrays,
converts external Cartesian ADM data into rescaled spherical conformal
gridfunctions, builds a BHaHAHA boundary structure, and applies inner
boundary/parity data.

For each multigrid angular resolution, `find_horizon` calls
`bah_numgrid__interp_src_set_up` to allocate the radial-spoke interpolation
source grid, perform angular interpolation from the external input grid with
`bah_interpolation_2d_external_input_to_interp_src_grid`, compute radial and
angular derivatives with `bah_hDD_dD_and_W_dD_in_interp_src_grid_interior`, and
apply radial/min-max/upwind boundary handling through
`bah_apply_bcs_r_maxmin_partial_r_hDD_upwinding`. It then allocates the evolved
grid through `bah_numgrid__evol_set_up`, which registers/evolves `hh` and
`vv`, stores metric and derivative auxevol fields, creates the spherical
evolution grid with `Nxx0 = 1`, allocates reference-metric precompute data, and
sets up the evolved-grid `bcstruct`.

`bah_initial_data` seeds `hh(theta, phi)` either by prolongating a coarser
horizon with `bah_interpolation_2d_general__uniform_src_grid`, by using a
fixed fraction of the full radial search sphere, or by quadratic extrapolation
from `prev_horizon_m1/m2/m3`. It sets `vv = eta_damping * hh` so the initial
surface has zero relaxation velocity, then calls `bah_apply_bcs_inner_only`.

Each relaxation step starts by interpolating metric data onto the current
surface with `bah_interpolation_1d_radial_spokes_on_3d_src_grid`, updates the
CFL-limited timestep with `bah_cfl_limited_timestep_based_on_h_equals_r`, and
may call `bah_over_relaxation`. Over-relaxation stores a previous surface in
`h_p`, scans extrapolated overstep factors, evaluates
`bah_diagnostics_area_centroid_and_Theta_norms` after trial metric
interpolations, accepts an improved overstep, resets `vv`, and refreshes the
surface interpolation and timestep.

The evolution RHS is registered by `register_CFunction_rhs_eval`: `hh_rhs` is
`vv - eta_damping * hh`, and `vv_rhs` is the negative BHaHAHA expansion
function `Theta` from the equation-side horizon diagnostics. `register_CFunction_KO_apply`
adds optional angular Kreiss-Oliger dissipation when `KO_diss_strength` is
nonzero. Both use BHaH simple loops over the surface interior with spherical
reference-metric precompute support.

Diagnostics run during the relaxation and once again at final resolution.
`bah_diagnostics` refreshes surface metric data when needed, calls
`bah_diagnostics_area_centroid_and_Theta_norms`, prints radius/Theta progress,
computes centroid-relative min/max/mean radii, computes coordinate-plane and
general spin-axis proper circumferences, cycles successful `m1/m2/m3` center,
radius, and shape history, and `bah_diagnostics_file_output` writes
AHFinderDirect-style output. On convergence, `find_horizon` stores the final
surface in `prev_horizon_m1`, shifts older shapes, and frees all temporary
per-resolution data through `free_all_but_external_input_gfs`; it also frees
external input grids before returning the final BHaHAHA error code.

`BHaHAHA_header.h` is the stable external C boundary. It defines
`bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`,
`NUM_EXT_INPUT_CARTESIAN_GFS`, `IDX2`, `MAX_RESOLUTIONS`, public poisoning
helpers, `bah_radial_grid_cell_centered_set_up`, `bah_xyz_center_r_minmax`,
`bah_find_horizon`, and `bah_diagnostics_file_output`. Generated NRPy helper
enums such as `INTERP_BSSN_GF_INDICES` and `FINAL_ADM_METRIC_INDICES` are
internal to the BHaH orchestrator rather than part of that public header API.

## Sources

- [BHaH_implementation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `bhahaha_find_horizons`, `initialize_bhahaha_solver_params_and_shapes`, `BHaHAHA_interpolate_metric_data_nrpy`, `BHaHAHA_BSSN_to_ADM_Cartesian`, `free_bhahaha_horizon_shape_data_all_horizons`
- [find_horizon.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/find_horizon.py) - `register_CFunction_find_horizon`, `find_horizon`, `free_all_but_external_input_gfs`
- [initial_data.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/initial_data.py) - `register_CFunction_initial_data`, `initial_data`
- [rhs_eval_KO_apply.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/rhs_eval_KO_apply.py) - `register_CFunction_rhs_eval`, `register_CFunction_KO_apply`
- [over_relaxation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/over_relaxation.py) - `register_CFunction_over_relaxation`, `over_relaxation`
- [numgrid__external_input_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__external_input_set_up.py) - `register_CFunction_numgrid__external_input_set_up`, `numgrid__external_input_set_up`
- [numgrid__interp_src_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__interp_src_set_up.py) - `register_CFunction_numgrid__interp_src_set_up`, `numgrid__interp_src_set_up`
- [numgrid__evol_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/numgrid__evol_set_up.py) - `register_CFunction_numgrid__evol_set_up`, `numgrid__evol_set_up`
- [radial_grid_cell_centered_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/radial_grid_cell_centered_set_up.py) - `register_CFunction_radial_grid_cell_centered_set_up`, `radial_grid_cell_centered_set_up`
- [xyz_center_r_minmax.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/xyz_center_r_minmax.py) - `register_CFunction_xyz_center_r_minmax`, `xyz_center_r_minmax`
- [bcstruct_set_up.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/bcstruct_set_up.py) - `register_CFunction_bcstruct_set_up`, `BHaH_defines_set_gridfunction_defines_with_parity_types`
- [apply_bcs_inner_only.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_inner_only.py) - `register_CFunction_apply_bcs_inner_only`, `apply_bcs_inner_only`
- [apply_bcs_r_maxmin_partial_r_hDD_upwinding.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/apply_bcs_r_maxmin_partial_r_hDD_upwinding.py) - `register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding`, `apply_bcs_r_maxmin_partial_r_hDD_upwinding`
- [interpolation_1d_radial_spokes_on_3d_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_1d_radial_spokes_on_3d_src_grid.py) - `register_CFunction_interpolation_1d_radial_spokes_on_3d_src_grid`, `interpolation_1d_radial_spokes_on_3d_src_grid`
- [interpolation_2d_external_input_to_interp_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_external_input_to_interp_src_grid.py) - `register_CFunction_interpolation_2d_external_input_to_interp_src_grid`, `interpolation_2d_external_input_to_interp_src_grid`
- [interpolation_2d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_2d_general__uniform_src_grid.py) - `register_CFunction_interpolation_2d_general__uniform_src_grid`, `interpolation_2d_general__uniform_src_grid`
- [interpolation_3d_general__uniform_src_grid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/interpolation_3d_general__uniform_src_grid.py) - `register_CFunction_interpolation_3d_general__uniform_src_grid`, `interpolation_3d_general__uniform_src_grid`
- [hDD_dD_and_W_dD_in_interp_src_grid_interior.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/hDD_dD_and_W_dD_in_interp_src_grid_interior.py) - `register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior`, `hDD_dD_and_W_dD_in_interp_src_grid_interior`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics.py) - `register_CFunction_diagnostics`, `diagnostics`
- [diagnostics_area_centroid_and_Theta_norms.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_area_centroid_and_Theta_norms.py) - `register_CFunction_diagnostics_area_centroid_and_Theta_norms`, `diagnostics_area_centroid_and_Theta_norms`
- [diagnostics_min_max_mean_radii_wrt_centroid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_min_max_mean_radii_wrt_centroid.py) - `register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid`, `diagnostics_min_max_mean_radii_wrt_centroid`
- [diagnostics_proper_circumferences.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences.py) - `register_CFunction_diagnostics_proper_circumferences`, `diagnostics_proper_circumferences`
- [diagnostics_proper_circumferences_general.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences_general.py) - `register_CFunction_diagnostics_proper_circumferences_general`, `diagnostics_proper_circumferences_general`
- [diagnostics_file_output.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - `register_CFunction_diagnostics_file_output`, `diagnostics_file_output`
- [BHaHAHA_header.h](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_params_and_data_struct`, `bhahaha_diagnostics_struct`, `NUM_EXT_INPUT_CARTESIAN_GFS`, `bah_find_horizon`

## See Also

- [BHaH](index.md)
- [Lifecycle And Project Assembly](lifecycle-and-project-assembly.md)
- [Horizon Diagnostics](../../equations/general-relativity/horizon-diagnostics.md)
- [Reference Metrics](../../core/reference-metrics.md)

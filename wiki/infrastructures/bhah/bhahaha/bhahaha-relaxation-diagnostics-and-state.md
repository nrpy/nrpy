# BHaHAHA Relaxation, Diagnostics, And State

> Explain multi-horizon BHaH state, BBH common-horizon activation, multiresolution hyperbolic relaxation, diagnostics, and cleanup. · Status: confirmed
> Up: [BHaHAHA](index.md)

## Summary

`bhahaha_find_horizons` manages persistent per-horizon center/radius/time/shape history and binary-black-hole (BBH) common-horizon activation in generated BHaH evolutions, while `find_horizon` solves one horizon at a time across configured multigrid `(Ntheta, Nphi)` levels using `hh_rhs = vv - eta_damping * hh` and `vv_rhs = -Theta` (with constant wavespeed `1` and optional angular Kreiss-Oliger dissipation).

## Detail

### Multi-Horizon History, Extrapolation, And BBH Mode

On iteration `0`, `bhahaha_find_horizons` seeds centers from `bah_initial_grid_{x,y,z}_center`, sets `t_m1/m2/m3 = -1.0`, sets radius bounds from `bah_max_search_radius`, forces a full-sphere initial guess, and allocates `prev_horizon_m1/m2/m3` at `max(Ntheta) x max(Nphi)`; later steps update non-persistent controls in place. In BBH mode (`3` horizons: `2` inspiral initially active, `1` common initially inactive), once both inspiral horizons succeed and their center separation plus maximum radii fit inside the common search diameter, the common search activates at the mass-weighted centroid with reset history and a full-sphere guess; once found, the two inspiral horizons deactivate.

Before each solve, `bah_xyz_center_r_minmax` calls `bah_quadratic_extrapolation` (`3` times -> quadratic, `2` -> linear, `<2` -> zeroth order); missing history lowers CFL, forces full-sphere mode, and on initial common-horizon attempts starts at the second multigrid resolution with `2x` `max_iterations`. `bhahaha_find_horizons` then interpolates Cartesian ADM metric data into `input_metric_data` and calls `bah_find_horizon`.

### Multiresolution Initialization, Relaxation Loop, And RHS

`bah_initial_data` initializes `hh(theta, phi)` by prolongating a coarser solved level when available, setting `0.8 * r_max_interior` in full-sphere mode, or extrapolating `prev_horizon_m1/m2/m3` in time; it sets `vv = eta_damping * hh` (`hh_rhs = 0`) and applies inner parity boundary conditions.

Each relaxation iteration runs radial-spoke metric interpolation, computes `dt = CFL_FACTOR * min(hh * dtheta, hh * dphi * sin(theta))` (`bah_cfl_limited_timestep_based_on_h_equals_r`), optionally tests linear overstep factors (`bah_over_relaxation` using `h_p`, accepting only sufficient `L_inf(Theta)` reductions and resetting `vv = eta_damping * hh`), runs scheduled diagnostics, checks `max_iterations` and `Theta` tolerances, and steps Method of Lines. `register_CFunction_rhs_eval` emits `hh_rhs = vv - eta_damping * hh` and `vv_rhs = -Theta` (plus `register_CFunction_KO_apply` when `KO_diss_strength != 0`; `variable_wavespeed.py` is unregistered).

### Diagnostics, Output Files, And Memory Cleanup

`bah_diagnostics` refreshes radial-spoke data and calls `bah_diagnostics_area_centroid_and_Theta_norms` (horizon area, centroid, `L2`/`L_inf` `Theta` norms, grid-relative radii, `Theta_eval_points_counter`, and `FIND_HORIZON_HORIZON_TOO_SMALL` if `r_min < 3 * dr`). On convergence at the finest level, `find_horizon` shifts `prev_horizon_m1/m2/m3`, updates `m1/m2/m3` center/radius/time history, and runs `bah_diagnostics_min_max_mean_radii_wrt_centroid` plus `bah_diagnostics_proper_circumferences` (`xy`/`xz`/`yz` circumferences and spin ratios; `bah_diagnostics_proper_circumferences_general` handles arbitrary spin axes).

If `bah_enable_spectre_spin_diagnostic` is enabled, one three-mode AKV solve
supplies the SpECTRE-style and Gram-matrix vectors and the lowest-shear
AKV-pole vector. A shared solve or integration failure leaves all three
vectors unavailable. A pole-fit failure leaves the AKV-pole vector
unavailable and preserves the other two vectors. The generalized
proper-circumference calculation uses the SpECTRE direction to obtain a
separate spin magnitude and vector.

`bah_diagnostics_file_output` appends to `BHaHAHA_diagnostics.ah%d.gp`. A new
file labels all 34 columns. Columns 22–24 contain the SpECTRE-style components,
25–27 the Gram-matrix components, 28–30 the lowest-shear AKV-pole components,
and 31–34 the circumference-ratio magnitude followed by its three components.
The writer emits `NaN` for unavailable spin values, retains any existing
header, and writes the Cartesian horizon coordinates to `h.t%07d.ah%d.gp`.

`free_all_but_external_input_gfs` frees per-level evolved/source grids, `rfmstruct`, `bcstruct`, MoL buffers, and `h_p`; `find_horizon` then frees `external_input_gfs`. `bhahaha_find_horizons` frees `input_metric_data` per solve and `prev_horizon_m1/m2/m3` via `free_bhahaha_horizon_shape_data_all_horizons` at shutdown.

## Sources

- [BHaH_implementation.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `bhahaha_find_horizons`, `initialize_bhahaha_solver_params_and_shapes`, `free_bhahaha_horizon_shape_data_all_horizons`
- [find_horizon.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/find_horizon.py) - `register_CFunction_find_horizon`, `find_horizon`, `free_all_but_external_input_gfs`
- [initial_data.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/initial_data.py) - `register_CFunction_initial_data`
- [xyz_center_r_minmax.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/xyz_center_r_minmax.py) - `register_CFunction_xyz_center_r_minmax`
- [quadratic_extrapolation.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/quadratic_extrapolation.py) - `register_CFunction_quadratic_extrapolation`
- [cfl_limited_timestep_based_on_h_equals_r.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/cfl_limited_timestep_based_on_h_equals_r.py) - `register_CFunction_cfl_limited_timestep_based_on_h_equals_r`
- [variable_wavespeed.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/variable_wavespeed.py) - `register_CFunction_variable_wavespeed`
- [rhs_eval_KO_apply.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/rhs_eval_KO_apply.py) - `register_CFunction_rhs_eval`, `register_CFunction_KO_apply`
- [over_relaxation.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/over_relaxation.py) - `register_CFunction_over_relaxation`
- [diagnostics.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics.py) - `register_CFunction_diagnostics`
- [diagnostics_area_centroid_and_Theta_norms.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_area_centroid_and_Theta_norms.py) - `register_CFunction_diagnostics_area_centroid_and_Theta_norms`
- [diagnostics_min_max_mean_radii_wrt_centroid.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_min_max_mean_radii_wrt_centroid.py) - `register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid`
- [diagnostics_proper_circumferences.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences.py) - `register_CFunction_diagnostics_proper_circumferences`
- [diagnostics_proper_circumferences_general.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences_general.py) - `register_CFunction_diagnostics_proper_circumferences_general`
- [diagnostics_integration_weights.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_integration_weights.py) - `register_CFunction_diagnostics_integration_weights`
- [diagnostics_file_output.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - `register_CFunction_diagnostics_file_output`

- [spectre_spin_integrator.py](../../../../nrpy/infrastructures/BHaH/BHaHAHA/spectre_spin_integrator.py) - `bah_compute_spectre_spin_potentials`, `diagnostics_spectre_spin`
- [BHaHAHA_header.h](../../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h) - `bhahaha_diagnostics_struct`, `bah_initialize_diagnostics_struct`

## See Also

- Depends on: [BHaHAHA Public C Header And Input Fields](bhahaha-public-api-and-input-contract.md)
- Depends on: [BHaHAHA Grid, Interpolation, And Boundaries](bhahaha-grid-interpolation-and-boundaries.md)
- Depends on: [Horizon Diagnostics](../../../equations/general-relativity/horizon-diagnostics.md)

# BHaHAHA Relaxation, Diagnostics, And State

> Explain the single-horizon relaxation loop, persistent horizon state, BBH common-horizon state, diagnostics, file output, and cleanup paths. · Status: confirmed · Last reconciled: 06-30-2026
> Up: [BHaH](index.md)

## Summary

The BHaHAHA relaxation/state leaf starts at the BHaH multi-horizon
orchestrator and narrows to the single-horizon `find_horizon` solver. The
orchestrator owns persistent per-horizon center, radius, time, previous-shape,
and BBH active/inactive state. The solver owns one candidate horizon at a
time: it accepts prepared Cartesian ADM metric input, sets up BHaHAHA-local
external/interpolation/evolution grids, advances the hyperbolic relaxation
through the configured multiresolution angular sequence, records diagnostics,
and returns a BHaHAHA error code.

The active RHS uses `hh_rhs = vv - eta_damping * hh` and
`vv_rhs = -Theta`, with optional angular Kreiss-Oliger dissipation. The
`variable_wavespeed.py` file still defines a registration helper for a variable
wavespeed routine, but the current BHaHAHA generator does not register it and
the active RHS path hard-codes wavespeed `1`, so it is not active solver
behavior.

## Detail

`bhahaha_find_horizons` is the persistent state owner for generated BHaH
evolutions. On the first evolution iteration it initializes each horizon's
center history from the configured initial centers, marks `t_m1/t_m2/t_m3` as
missing with `-1.0`, initializes radius history from the maximum search radius,
forces a fixed full-sphere first guess, and allocates the
`prev_horizon_m1/m2/m3` shape arrays at the finest configured angular
resolution. Later calls refresh non-persistent solver controls without
reallocating the shape history.

BBH mode is stateful in the same orchestrator. It requires three horizons:
two individual inspiral horizons start active and the common horizon starts
inactive. Once both individual horizons have recent successful finds, the
orchestrator tests whether their center separation plus both stored maximum
radii fits within the common horizon's configured search diameter. If it
does, it activates the common search, seeds the common center with a
mass-weighted average of the two individual centers, resets the common
history to missing, and forces a fixed full-sphere search. After the common
horizon has been found and remains active, the individual horizons are
deactivated.

Before each active horizon solve, `bhahaha_find_horizons` extrapolates the
next center and radial search interval through `bah_xyz_center_r_minmax`.
That helper uses `bah_quadratic_extrapolation`, which performs quadratic
extrapolation when three history times exist, linear extrapolation when two
exist, and a zeroth-order estimate when history is missing. Unreliable history
reduces the horizon CFL factor, forces a full-sphere guess, and, for the
initial BBH common-horizon attempts, promotes the first multigrid resolution
to the second configured level and doubles the maximum iteration budget.

The BHaH orchestrator then builds the per-horizon radial interpolation range,
allocates `input_metric_data`, fills it through the BHaH metric interpolation
bridge, and calls `bah_find_horizon`. The solver treats that metric buffer as
external input. Its first handoff is `bah_numgrid__external_input_set_up`,
which constructs BHaHAHA's external input grid from the supplied Cartesian ADM
data. Each multiresolution level then calls `bah_numgrid__interp_src_set_up`
for the interpolation-source grid and `bah_numgrid__evol_set_up` for the
surface evolution grid. Detailed grid, interpolation, and boundary mechanics
belong to the BHaHAHA grid/interpolation owner; this page owns the lifecycle
handoff and solver use.

`bah_initial_data` seeds `hh(theta, phi)` for each resolution from one of
three sources. If a coarser resolved horizon exists, it prolongates that
surface onto the finer angular grid. If the horizon is in fixed-radius
full-sphere mode, it initializes to `0.8` times the current interior outer
search radius. Otherwise it extrapolates each angular point from
`prev_horizon_m1/m2/m3` at the current external input time. It sets
`vv = eta_damping * hh`, making the initial `hh_rhs` zero, then applies inner
boundary conditions and marks the coarse-horizon path available for the next
resolution.

Each relaxation step interpolates the metric data onto the current surface,
sets `dt` with `bah_cfl_limited_timestep_based_on_h_equals_r`, optionally
runs `bah_over_relaxation`, runs diagnostics on the configured cadence, tests
the maximum-iteration and Theta-norm convergence conditions, and advances the
Method-of-Lines state when no stop condition is met. The CFL helper computes
the minimum angular surface spacing from `hh * dtheta` and
`hh * dphi * sin(theta)` and multiplies it by `CFL_FACTOR`.

`bah_over_relaxation` stores a previous surface in `h_p`, periodically scans
linear extrapolation overstep factors, evaluates area/centroid/Theta
diagnostics on each trial surface after refreshing radial-spoke interpolation,
accepts an improved overstep only when it reduces the L-infinity Theta norm
enough, resets `vv = eta_damping * hh`, refreshes the surface interpolation,
updates the CFL timestep, and clears the stored overstep surface.

`register_CFunction_rhs_eval` is the active relaxation RHS registration. It
builds `hh_rhs` from `vv - eta_damping * hh` and `vv_rhs` from the negative
equation-side expansion function `Theta`. `register_CFunction_KO_apply` adds
optional angular Kreiss-Oliger terms to both RHS components only when
`KO_diss_strength` is nonzero. `register_CFunction_variable_wavespeed`
generates a routine that can populate an auxiliary variable wavespeed field,
but the active RHS explicitly comments out the symbolic variable and uses
constant wavespeed `1`.

`bah_diagnostics` refreshes radial-spoke interpolation before its cadence
work. `bah_diagnostics_area_centroid_and_Theta_norms` selects
cell-centered integration weights, accumulates horizon area, area-weighted
coordinate centroid, area-weighted L2 norm of `Theta`, L-infinity norm of
`Theta`, and internal min/max radii relative to the current grid center. It
also increments the Theta evaluation point counter. If the grid-center
minimum radius falls below three radial spacings, diagnostics set the
`FIND_HORIZON_HORIZON_TOO_SMALL` error.

Final diagnostics add centroid-relative radius and circumference work.
`bah_diagnostics_min_max_mean_radii_wrt_centroid` computes minimum, maximum,
and area-weighted mean coordinate radius relative to the diagnostic centroid.
`bah_diagnostics_proper_circumferences` computes coordinate-plane proper
circumferences in the `xy`, `xz`, and `yz` planes, estimates spin magnitudes
from those circumference ratios when the ratio inversion is valid, and returns
diagnostic allocation or interpolation errors to the caller. The separately
registered `bah_diagnostics_proper_circumferences_general` computes polar and
equatorial proper circumferences for a supplied spin axis and stores the
general-axis ratio-based spin estimate when that routine is invoked.

On a successful final-resolution solve, `find_horizon` shifts the previous
surface arrays so the new surface becomes `prev_horizon_m1`, marks final
diagnostics, and calls `bah_diagnostics` once with cadence forced to every
iteration. Final diagnostics also cycle the center, radius, and time history:
older `m1/m2` entries shift to `m2/m3`, new `m1` center is the extrapolated
search center plus the measured centroid offset, and new `m1` radii come from
the centroid-relative min/max diagnostics.

`bah_diagnostics_file_output` is the local AHFinderDirect-style output writer
used after a successful orchestrator solve. It appends one row to
`BHaHAHA_diagnostics.ah%d.gp`, writing iteration, time, centroid,
centroid-relative radii, coordinate-plane circumferences, area, irreducible
mass, Theta norms, and ratio-based spin estimates. It also writes a
gnuplot-compatible surface file `h.t%07d.ah%d.gp` containing Cartesian horizon
surface points from `prev_horizon_m1`.

Cleanup is split by layer. Inside `find_horizon`,
`free_all_but_external_input_gfs` frees each per-resolution evolution grid,
reference-metric precompute data, boundary arrays, MoL storage, coordinate
arrays, interpolation-source storage, and over-relaxation history. The solver
then frees external input coordinate arrays and gridfunctions after all
resolutions finish or after an error path exits the loop. The BHaH
orchestrator frees the temporary `input_metric_data` buffer after each
per-horizon solve and uses
`free_bhahaha_horizon_shape_data_all_horizons` to release persistent
`prev_horizon_m1/m2/m3` arrays when the generated application cleans up.

## Sources

- [BHaH_implementation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/BHaH_implementation.py) - `register_CFunction_bhahaha_find_horizons`, `bhahaha_find_horizons`, `initialize_bhahaha_solver_params_and_shapes`, `free_bhahaha_horizon_shape_data_all_horizons`
- [find_horizon.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/find_horizon.py) - `register_CFunction_find_horizon`, `find_horizon`, `free_all_but_external_input_gfs`
- [initial_data.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/initial_data.py) - `register_CFunction_initial_data`, `initial_data`
- [xyz_center_r_minmax.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/xyz_center_r_minmax.py) - `register_CFunction_xyz_center_r_minmax`, `xyz_center_r_minmax`
- [quadratic_extrapolation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/quadratic_extrapolation.py) - `register_CFunction_quadratic_extrapolation`, `quadratic_extrapolation`
- [cfl_limited_timestep_based_on_h_equals_r.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/cfl_limited_timestep_based_on_h_equals_r.py) - `register_CFunction_cfl_limited_timestep_based_on_h_equals_r`, `cfl_limited_timestep_based_on_h_equals_r`
- [variable_wavespeed.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/variable_wavespeed.py) - `register_CFunction_variable_wavespeed`, `variable_wavespeed`
- [rhs_eval_KO_apply.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/rhs_eval_KO_apply.py) - `register_CFunction_rhs_eval`, `register_CFunction_KO_apply`
- [over_relaxation.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/over_relaxation.py) - `register_CFunction_over_relaxation`, `over_relaxation`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics.py) - `register_CFunction_diagnostics`, `diagnostics`
- [diagnostics_area_centroid_and_Theta_norms.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_area_centroid_and_Theta_norms.py) - `register_CFunction_diagnostics_area_centroid_and_Theta_norms`, `diagnostics_area_centroid_and_Theta_norms`
- [diagnostics_min_max_mean_radii_wrt_centroid.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_min_max_mean_radii_wrt_centroid.py) - `register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid`, `diagnostics_min_max_mean_radii_wrt_centroid`
- [diagnostics_proper_circumferences.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences.py) - `register_CFunction_diagnostics_proper_circumferences`, `diagnostics_proper_circumferences`
- [diagnostics_proper_circumferences_general.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_proper_circumferences_general.py) - `register_CFunction_diagnostics_proper_circumferences_general`, `diagnostics_proper_circumferences_general`
- [diagnostics_integration_weights.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_integration_weights.py) - `register_CFunction_diagnostics_integration_weights`, `diagnostics_integration_weights`
- [diagnostics_file_output.py](../../../nrpy/infrastructures/BHaH/BHaHAHA/diagnostics_file_output.py) - `register_CFunction_diagnostics_file_output`, `diagnostics_file_output`

## See Also

- [BHaH](index.md)
- [BHaHAHA Horizon Runtime](bhahaha-horizon-runtime.md)
- [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
- [GR Application Wiring](gr-application-wiring.md)
- [Horizon Diagnostics](../../equations/general-relativity/horizon-diagnostics.md)

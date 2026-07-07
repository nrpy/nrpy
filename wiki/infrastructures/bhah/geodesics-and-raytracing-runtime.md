# Geodesics And Raytracing Runtime

> BHaH runtime pieces for standalone geodesics and evolution-time raytracing export. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

BHaH has two distinct raytracing-facing paths. The standalone photon geodesic
runtime builds a generated `main`, initializes camera/window/source parameters,
tiles the camera window, integrates batches of photon geodesics with RKF45, and
serializes per-tile light-blueprint results. The evolution-time export path is
part of diagnostics: it writes Cartesian `g4DD` and `Gamma4UDD` time-slice data
from a live BSSN evolution, then a combiner validates and stacks those slices
for later numerical-spacetime interpolation.

## Detail

The standalone photon entrypoint is `main` in the photon geodesics package. It
registers window tiling parameters, initializes `commondata`, parses command
line and parfile input, computes a camera basis with fallback logic for
near-degenerate up vectors, loops over `(tx, ty)` tiles, shifts the active
window center per tile, calls `batch_integrator_numerical`, shifts local tile
hits back to the global window coordinate system, and writes one result set per
tile. This path is a standalone geodesic program; it is not the same as
diagnostics emitted by an evolving BHaH spacetime.

`batch_integrator_numerical` is the host orchestrator for photon batches. It
registers integration limits and RKF45 controls in `commondata`, allocates the
Structure-of-Arrays photon state, history, step-size, status, event-lock, and
result buffers, sets up CUDA streams or CPU equivalents, and uses the
TimeSlotManager to process active rays by coordinate-time slots. The split
pipeline stages exchange flat bundles for state, metric, connection, derivative
stages, affine parameter, retries, and termination state.

The RKF45 kernels are deliberately split. `interpolation_kernel` evaluates
spacetime-specialized `g4DD_metric` and `connections` helpers for each ray and
writes 10 metric and 40 connection components into bundles. `calculate_ode_rhs_kernel`
unpacks coordinates, momenta, metric entries, and Christoffels, evaluates the
nine photon RHS expressions, and writes the selected RK stage into
`d_k_bundle`. `rkf45_stage_update` reads the base state, stage derivatives, and
per-ray step size, applies the RKF45 Butcher coefficients for stages 1-5, skips
stage 6, and writes the temporary state for the next stage. The finalization and
control kernel, event manager, and time-slot helpers then decide which rays
remain active and which step sizes advance.

`event_detection_manager_kernel` owns geometric termination and result capture.
It rejects rays whose temporal momentum exceeds `p_t_max`, marks rays outside
`r_escape`, checks crossings of the immutable global window plane and source
plane using current and two historical states, calls
`find_event_time_and_state` to reconstruct the crossing, delegates physical
coordinate extraction to `handle_window_plane_intersection` or
`handle_source_plane_intersection`, and shifts the state history only for rays
that remain active.

Metric and connection generation is shared by photon and massive geodesic
runtimes. `g4DD_metric` writes the upper-triangular 10-component covariant
metric into a thread-local array, while `connections` writes the 40 unique
Christoffel components. Both unpack only coordinates used by the generated
SymPy expressions and validate the particle-state width by `PARTICLE`.
Massive geodesics use the GSL path instead of the batched photon RKF45 path:
the spacetime-specialized massive-geodesic GSL wrapper casts the GSL parameter
pointer to `commondata_struct`, evaluates metric and connection locally, calls
`calculate_ode_rhs_massive`, and returns `GSL_SUCCESS`.

Numerical-spacetime interpolation is separate from analytic metric evaluation.
`register_CFunction_numerical_interpolation` emits a CPU wrapper that assumes a
mapped `NumericalTimeWindowManager`; for every ray in a chunk it asks the time
window for the temporal stencil, performs spherical azimuthal-symmetry spatial
Lagrange interpolation on every mapped slice, then calls temporal Lagrange
interpolation at the photon's coordinate time. The time-window manager owns the
read-only mmap window over a combined numerical-spacetime container and widens
each slot by temporal interpolation halo plus backward RKF45 lookahead.

Evolution-time raytracing export lives under BHaH diagnostics. When enabled by
`register_all_diagnostics`, the diagnostics function calls
`output_raytracing_data` on scheduled output steps. The exporter requires a
host/OpenMP build, `enable_rfm_precompute=True`,
`enable_RbarDD_gridfunctions=True`, one active grid, and source coordinate
system `Cartesian` or `Spherical`. It refreshes same-slice Ricci/RHS scratch
data, evaluates Cartesian-basis `g4DD` and `Gamma4UDD` symbolic recipes on
interior points, writes fixed-width metadata and binary64 point records, and
records that ghost-zone points are excluded. `combine_raytracing_time_slices.py`
then parses those stage-1 files, validates headers, sorts by simulation time,
and writes a read-only stacked container for downstream interpolation.

## Sources

- [main.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main.py) - `main`
- [batch_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py) - `batch_integrator_numerical`
- [time_slot_manager_helpers.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/time_slot_manager_helpers.py) - `time_slot_manager_helpers`, `TimeSlotManager`
- [interpolation_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/interpolation_kernel.py) - `interpolation_kernel`
- [calculate_ode_rhs_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_ode_rhs_kernel.py) - `calculate_ode_rhs_kernel`
- [rkf45_stage_update.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_stage_update.py) - `rkf45_stage_update`
- [rkf45_finalize_and_control_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py) - `rkf45_finalize_and_control_kernel`, `rkf45_finalize_and_control`
- [event_detection_manager_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py) - `event_detection_manager_kernel`
- [find_event_time_and_state.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py) - `find_event_time_and_state`
- [handle_window_plane_intersection.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_window_plane_intersection.py) - `handle_window_plane_intersection`
- [handle_source_plane_intersection.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_source_plane_intersection.py) - `handle_source_plane_intersection`
- [g4DD_metric.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/g4DD_metric.py) - `g4DD_metric`
- [connections.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/connections.py) - `connections`
- [ode_gsl_wrapper_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/ode_gsl_wrapper_massive.py) - `ode_gsl_wrapper_massive`
- [calculate_ode_rhs_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/calculate_ode_rhs_massive.py) - `calculate_ode_rhs_massive`
- [numerical_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py) - `register_CFunction_numerical_interpolation`
- [time_window_manager_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py) - `time_window_manager_numerical`, `NumericalTimeWindowManager`
- [azimuthal_symmetry_spatial_lagrange_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py) - `register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation`
- [output_raytracing_data.py](../../../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py) - `register_CFunction_output_raytracing_data`, `raytracing_data_point_index_from_logical_indices`
- [combine_raytracing_time_slices.py](../../../nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py) - `Stage1Info`, `Layout`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`, `enable_raytracing_data_output`

## See Also

- [BHaH](index.md)
- [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
- [Geodesics](../../equations/general-relativity/geodesics.md)
- [Black Hole Evolution](../../examples/black-hole-evolution.md)

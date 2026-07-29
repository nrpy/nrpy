# Geodesics And Raytracing Runtime

> BHaH runtime pieces for standalone geodesics and evolution-time raytracing export. · Status: confirmed · Last reconciled: 06-29-2026
> Up: [BHaH](index.md)

## Summary

BHaH has two distinct raytracing-facing paths. The standalone photon geodesic
runtime builds a generated `main`, initializes observer and independent plane
parameters, tiles the angular pixel grid, integrates batches of photon geodesics with RKF45,
and serializes per-tile light-blueprint results. The evolution-time export path is
part of diagnostics: it writes Cartesian `g4DD` and `Gamma4UDD` time-slice data
from a live BSSN evolution, then a combiner validates and stacks those slices
for later numerical-spacetime interpolation.

## Detail

The standalone photon entrypoint is `main` in the photon geodesics package. It
registers angular tile-grid parameters, initializes `commondata`, parses
command-line and parfile input, computes an observer tetrad with fallback logic
for near-degenerate up vectors, loops over `(tx, ty)` tile indices, and calls
the selected batch integrator. Each record receives normalized image-sample
coordinates before serialization; image placement does not depend on integer
pixel identity fields.
This path is a standalone geodesic program; it is not the same as
diagnostics emitted by an evolving BHaH spacetime.

Claim evidence:
- Claim: The standalone photon entrypoint registers angular tile-grid controls, initializes an observer tetrad, assigns normalized image-sample coordinates, and is distinct from evolution-time diagnostics.
- Role: generated-output boundary
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py` — `main`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_and_fill_blueprint_data_universal.py` — normalized record fields
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

The shared `set_initial_conditions_kernel` receives one metric evaluated at the
observer event and constructs one validated metric-orthonormal tetrad per tile
batch call. It uses that tetrad for every ray in the call, writes complete
contravariant `p^mu`, and only then performs any requested normalized-variable
conversion. The event-coordinate-independent plane bases serve event-coordinate
diagnostics; they are not the metric tetrad used to initialize momentum.

Claim evidence:
- Claim: `set_initial_conditions_kernel` constructs one validated metric-orthonormal observer tetrad per batch call, writes complete contravariant momentum, and performs normalized-variable conversion afterward when requested.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py` — observer initialization
- Corroboration: `nrpy/examples/photon_batch_geodesic_integrator_numerical.py` — shared observer initialization arguments
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`batch_integrator_numerical` is the host orchestrator for photon batches. It
registers integration limits and RKF45 controls in `commondata`, allocates the
Structure-of-Arrays photon state, history, step-size, status, event-lock, and
result buffers, sets up CUDA streams or CPU equivalents, and uses the
TimeSlotManager to process active rays by coordinate-time slots. The split
pipeline stages exchange flat bundles for state, metric, connection, derivative
stages, affine parameter, retries, and termination state.

Claim evidence:
- Claim: `batch_integrator_numerical` orchestrates host-side photon batches, manages flat state/history/geometry/result bundles, and processes active rays through coordinate-time slots.
- Role: generated-output boundary
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py` — `batch_integrator_numerical`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/time_slot_manager_helpers.py` — `TimeSlotManager`
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

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

Claim evidence:
- Claim: The split RKF45 pipeline writes ten metric and forty Christoffel components, evaluates the selected photon RHS at each stage, skips the stage-6 intermediate update, and delegates acceptance/control and termination decisions to the remaining kernels.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/interpolation_kernel.py`, `calculate_ode_rhs_kernel.py`, and `rkf45_stage_update.py` — generated kernel interfaces
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py` — acceptance/control
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

`event_detection_manager_kernel` owns geometric termination and result capture.
It rejects rays whose temporal momentum exceeds `p_t_max`, marks rays outside
`r_escape`, checks crossings of the independent nonterminal and terminal planes
using current and two historical states, calls `find_event_time_and_state` to
reconstruct the crossing, delegates physical coordinate extraction to
`handle_non_terminal_plane_intersection` or `handle_terminal_plane_intersection`,
and shifts the state history only for rays that remain active.

Claim evidence:
- Claim: `event_detection_manager_kernel` owns photon termination and result capture for momentum limits, escape radius, independent plane crossings, event-state reconstruction, and active-ray history shifts.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py` — `event_detection_manager_kernel`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py` — event reconstruction
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Blueprint headers carry tile identity/counts, `alpha_w`, `alpha_h`, and schema
version 6. Records carry plane diagnostics, final angles, termination times,
and normalized image-sample fractions.
The renderer places rays from the normalized fractions and preserves the
vertical raster flip, while plane diagnostics remain available to
`blueprint_analysis.py`.

Claim evidence:
- Claim: Blueprint schema version 6 records tile/header geometry, plane diagnostics, termination times, final angles, and normalized image-sample fractions; the renderer preserves the documented vertical raster flip.
- Role: generated-output boundary
- Deciding authority: `nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py` — schema constants and dtype
- Corroboration: `nrpy/examples/geodesic_visualizations/visualize_lensed_image.py` — image placement
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Metric and connection generation is shared by photon and massive geodesic
runtimes. `g4DD_metric` writes the upper-triangular 10-component covariant
metric into a thread-local array, while `connections` writes the 40 unique
Christoffel components. Both unpack only coordinates used by the generated
SymPy expressions and validate the particle-state width by `PARTICLE`.
Massive geodesics use the GSL path instead of the batched photon RKF45 path:
the spacetime-specialized massive-geodesic GSL wrapper casts the GSL parameter
pointer to `commondata_struct`, evaluates metric and connection locally, calls
`calculate_ode_rhs_massive`, and returns `GSL_SUCCESS`.

Claim evidence:
- Claim: Analytic metric/connection helpers write ten metric and forty Christoffel components, while the massive GSL wrapper evaluates the specialized RHS and returns `GSL_SUCCESS`.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/g4DD_metric.py`, `connections.py`, and `massive/ode_gsl_wrapper_massive.py` — generated helper registrations
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/calculate_ode_rhs_massive.py` — massive RHS
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Numerical-spacetime interpolation is separate from analytic metric evaluation.
`register_CFunction_numerical_interpolation` emits a CPU wrapper that assumes a
mapped `NumericalTimeWindowManager`; for every ray in a chunk it asks the time
window for the temporal stencil, performs `SinhCylindricalv2n2` azimuthal-
symmetry spatial Lagrange interpolation on every mapped slice, then calls
temporal Lagrange interpolation at the photon's coordinate time. The time-window
manager owns the read-only mmap window over a combined numerical-spacetime
container and widens each slot by temporal interpolation halo plus backward
RKF45 lookahead. The exporter and combiner preserve the full logical grid,
including ghost-zone points; `InputSliceInfo` describes source slices in the
combiner. Stored slice-table times are authoritative: the first stored time is
loaded into runtime `t_numerical_initial` and printed at startup, while the
user supplies only `t_numerical_end`. The nominal `dt_numerical_spacetime_data`
describes approximate evolution spacing and is used only when synthetic
stencil edge times are required; stored output times may be nearby, nonuniform
values.

Claim evidence:
- Claim: Numerical interpolation uses authoritative combined-container slice times, preserves the full logical grid including ghost zones, and uses nominal spacing only for approximate synthetic temporal-stencil edge times.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py` — slice-table loading and stencil construction
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py` — combined layout and metadata
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Numerical endpoint dispatch is piecewise constant. At or below the first stored
time, the first slice is spatially interpolated; at or above the selected final
slice, that final slice is reused. `g4DD` obtains temporal metric derivatives
from temporal interpolation, `g4DD_d0` uses stored metric derivatives and
zeroes their temporal slots only for direct endpoint dispatch, and `GammaUDD`
reuses stored Christoffels at endpoints. Mixed `g4DD_d0` temporal padding
remains a separate interpolation behavior; it is not controlled by
`--raytracing-static-christoffels`. That flag changes the exported final
Christoffel payload only.

Claim evidence:
- Claim: Numerical endpoint dispatch reuses the first/final stored slices as piecewise-constant continuations; `g4DD`, `g4DD_d0`, and `GammaUDD` obtain their geometry according to their selected payload contracts, while `--raytracing-static-christoffels` affects only exported GammaUDD data.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py` — endpoint dispatch; `output_raytracing_data.py` — payload mode
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/temporal_lagrange_interpolation.py` — temporal endpoint handling
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Startup reads `.bin` `NGHOSTS` metadata and terminates if the configured spatial
half-width requires more ghost zones than available. It also terminates when
observer or positive Cartesian `r_escape` probes cannot support the requested
native stencil; the relevant `.par` controls are
`observer_x/y/z`, `r_escape`, and
`numerical_spacetime_spatial_interp_half_width`.

Claim evidence:
- Claim: Numerical startup rejects insufficient ghost-zone capacity and observer or positive Cartesian `r_escape` probes that cannot support the requested native stencil.
- Role: user-facing commands and interfaces
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py` — startup domain validation
- Corroboration: `nrpy/infrastructures/BHaH/xx_tofrom_Cart.py` — Cartesian-to-native inversion and admission checks
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

Evolution-time raytracing export lives under BHaH diagnostics. When enabled by
`register_all_diagnostics`, the diagnostics function calls
`output_raytracing_data` on scheduled output steps. The exporter requires a
host/OpenMP build, `enable_rfm_precompute=True`,
`enable_RbarDD_gridfunctions=True`, one active grid, and the raytracing
example's `SinhCylindricalv2n2` source coordinate system. It refreshes
same-slice Ricci/RHS scratch data, evaluates Cartesian-basis `g4DD` and
`Gamma4UDD` symbolic recipes on interior points, fills the full logical-grid
payload including ghost zones, writes fixed-width metadata and binary64 point
records, and preserves the serialized logical-grid bounds. Optional
`--raytracing-static-christoffels` changes only the selected GammaUDD values in
the qualifying final output slice; interpolation consumes whichever values were
stored and applies its own endpoint policy. `combine_raytracing_time_slices.py`
then parses those stage-1 files, validates headers, sorts by simulation time,
and writes a read-only stacked container for downstream interpolation.

Claim evidence:
- Claim: Evolution-time export is host/OpenMP, single-grid, SinhCylindricalv2n2-only in the documented path; it writes mode-selected Cartesian metric/geometry payloads over the full logical grid, and optional static Christoffels affect only GammaUDD values on the last scheduled output while the default remains nonstatic.
- Role: generated-output boundary
- Deciding authority: `nrpy/infrastructures/BHaH/diagnostics/diagnostics.py` and `output_raytracing_data.py` — diagnostics registration and exporter
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py` — stage-1 validation and combined container
- Validation: `inspected=pass; generated=not-run; built=not-run; run=not-run; result_checked=not-run`
- Dimensions: `platform=not-applicable; tool_version=not-applicable; backend=not-applicable; precision=not-applicable; GPU=not-applicable; restart=not-applicable; distributed=not-applicable; error_path=not-applicable; options=not-applicable; date=07-28-2026`

## Sources

- [main_batch.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py) - `main`
- [batch_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py) - `batch_integrator_numerical`
- [time_slot_manager_helpers.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/time_slot_manager_helpers.py) - `time_slot_manager_helpers`, `TimeSlotManager`
- [interpolation_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/interpolation_kernel.py) - `interpolation_kernel`
- [calculate_ode_rhs_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_ode_rhs_kernel.py) - `calculate_ode_rhs_kernel`
- [rkf45_stage_update.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_stage_update.py) - `rkf45_stage_update`
- [rkf45_finalize_and_control_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py) - `rkf45_finalize_and_control_kernel`, `rkf45_finalize_and_control`
- [event_detection_manager_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py) - `event_detection_manager_kernel`
- [find_event_time_and_state.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py) - `find_event_time_and_state`
- [handle_non_terminal_plane_intersection.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_non_terminal_plane_intersection.py) - `handle_non_terminal_plane_intersection`
- [handle_terminal_plane_intersection.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/handle_terminal_plane_intersection.py) - `handle_terminal_plane_intersection`
- [g4DD_metric.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/g4DD_metric.py) - `g4DD_metric`
- [connections.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/connections.py) - `connections`
- [ode_gsl_wrapper_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/ode_gsl_wrapper_massive.py) - `ode_gsl_wrapper_massive`
- [calculate_ode_rhs_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/calculate_ode_rhs_massive.py) - `calculate_ode_rhs_massive`
- [numerical_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py) - `register_CFunction_numerical_interpolation`
- [time_window_manager_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py) - `time_window_manager_numerical`, `NumericalTimeWindowManager`
- [azimuthal_symmetry_spatial_lagrange_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py) - `register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation`
- [output_raytracing_data.py](../../../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py) - `register_CFunction_output_raytracing_data`, `raytracing_data_point_index_from_logical_indices`
- [combine_raytracing_time_slices.py](../../../nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py) - `InputSliceInfo`, `Layout`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`, `enable_raytracing_data_output`

## See Also

- [BHaH](index.md)
- [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
- [Geodesics](../../equations/general-relativity/geodesics.md)
- [Black Hole Evolution](../../examples/black-hole-evolution.md)

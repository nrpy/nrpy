# Geodesics And Raytracing Runtime

> BHaH runtime pieces for standalone geodesics and evolution-time raytracing export. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH has standalone single- and batch-photon geodesic programs plus an
evolution-time raytracing-data export. Both analytical and numerical
single-photon programs support optional terminal and nonterminal planes. The
batch program initializes observer and independent plane parameters, tiles the
angular pixel grid, integrates photon groups with RKF45, and writes per-tile
light-blueprint binary files. The evolution-time export writes mode-selected
Cartesian `g4DD`, `g4DD_d0`, and `Gamma4UDD` time-slice data from a live BSSN
evolution; a combiner validates and stacks those slices for later
numerical-spacetime interpolation.

## Detail

The standalone batch-photon entrypoint is `main` in the photon geodesics package. It
registers angular tile-grid parameters, initializes `commondata`, parses
command-line and parfile input, loops over `(tx, ty)` tile indices, and calls
the selected batch integrator. The shared initializer invoked by that integrator
constructs the observer tetrad, including fallback logic for near-degenerate up
vectors. Each record receives normalized image-sample
coordinates before serialization; image placement does not depend on integer
pixel identity fields.
This path is a standalone geodesic program; it is not the same as
diagnostics emitted by an evolving BHaH spacetime. Single-photon generators use
the particle-independent forwarding `main` from `geodesics/main_single.py`.
Both single-photon generators register optional event-plane parameters and the
crossing functions used by their integrators.

Claim evidence:
- Claim: The standalone batch-photon entrypoint registers angular tile-grid controls, delegates observer-tetrad construction to the shared initializer, assigns normalized image-sample coordinates, and is distinct from both single-photon programs and evolution-time diagnostics.
- Role: generated evidence
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py` — `main`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/calculate_and_fill_blueprint_data_universal.py` — normalized record fields

The shared `set_initial_conditions_kernel` receives one metric evaluated at the
observer event and constructs one validated metric-orthonormal tetrad per call.
It uses that tetrad for every requested ray, writes complete contravariant
`p^mu`, and only then performs any requested normalized-variable conversion.
Batch generators and both single-photon generators request initialization of
event-plane-side history when registering their shared initializer.
Event-plane bases used by crossing handlers are not the metric tetrad used to
initialize momentum.
The generated ray construction sets the initial camera-tetrad photon energy
magnitude to `E_camera = 1` before any normalized-variable conversion. This
affine normalization is distinct from the hypersurface-normal measure
`|alpha p^0|` used by later termination diagnostics.

Claim evidence:
- Claim: `set_initial_conditions_kernel` constructs one validated metric-orthonormal camera tetrad per call, initializes unit camera-tetrad energy magnitude `E_camera=1`, writes complete contravariant momentum, optionally initializes batch or single-photon event-plane-side history, and performs normalized-variable conversion afterward when requested; the later hypersurface-normal energy magnitude is `|alpha p^0|`.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py` — observer initialization
- Corroboration: `nrpy/examples/photon_batch_geodesic_integrator_numerical.py`, `nrpy/examples/photon_single_geodesic_integrator_analytical.py`, `nrpy/examples/photon_single_geodesic_integrator_numerical.py`, `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py`, and `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py` — shared observer and event-history initialization arguments

`batch_integrator_numerical` is the host orchestrator for photon batches. It
registers integration limits and RKF45 controls in `commondata`, allocates the
Structure-of-Arrays photon state, history, step-size, status, event-lock, and
result buffers, sets up CUDA streams or CPU equivalents, and uses the
TimeSlotManager to process active rays by coordinate-time slots. The split
pipeline stages exchange flat bundles for state, metric, connection, derivative
stages, affine parameter, retries, and termination state. After an accepted RKF45
step, direct mode refreshes the metric at the accepted state using the same
trial-locked spatial stencil centers used by the RK stages. The integrator then
fills a common per-ray log-energy bundle: normalized mode copies `u`, while
direct mode computes `ln|alpha p^0|` with `normal_observer_log_energy`. The event
manager applies the shared upper-only `evolution_measure_max` cutoff. Its
normalization sidecar stores `|C_direct|` for direct EOM and
`|exp(2u)(C_normalized - 1)|` for normalized EOM.

Claim evidence:
- Claim: `batch_integrator_numerical` orchestrates host-side photon batches, manages flat state/history/geometry/result bundles, reuses trial-locked spatial centers for accepted-state metric refresh, computes the common log-energy bundle, and processes active rays through coordinate-time slots.
- Role: generated evidence
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py` — `batch_integrator_numerical`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/time_slot_manager_helpers.py` — `TimeSlotManager`

The RKF45 kernels are deliberately split. `interpolation_kernel` evaluates
spacetime-specialized `g4DD_metric` and `connections` helpers for each ray and
writes 10 metric and 40 connection components into bundles. `calculate_ode_rhs_kernel`
unpacks coordinates, momenta, metric entries, and Christoffels, evaluates the
nine photon RHS expressions, and writes the selected RK stage into
`d_k_bundle`. `rkf45_stage_update` reads the base state, stage derivatives, and
per-ray step size, applies the RKF45 Butcher coefficients for stages 1-5, skips
stage 6, and writes the temporary state for the next stage. The finalization and
control kernel commits accepted states; the accepted-state metric refresh,
log-energy calculation, event manager, and time-slot helpers then decide which
rays remain active and which step sizes advance.

Claim evidence:
- Claim: The split RKF45 pipeline writes ten metric and forty Christoffel components, evaluates the selected photon RHS at each stage, skips the stage-6 intermediate update, commits accepted states, and delegates downstream metric, energy, control, and termination decisions to the remaining kernels.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/interpolation_kernel.py`, `calculate_ode_rhs_kernel.py`, and `rkf45_stage_update.py` — generated kernel interfaces
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py` — acceptance/control

`event_detection_manager_kernel` applies log-energy and coordinate-radius
termination rules and checks enabled planes. Both single-photon integrators call
it after each accepted step when a plane is enabled. The manager marks a sign
change of the plane function as pending, counts subsequent accepted states, and
calls `find_event_time_and_state_centered` when the requested centered stencil
is available or a stop requires an earlier fit. Rejected RKF45 steps do not
advance this count. The interpolator tries the requested polynomial degree,
then lower even degrees through 2. It checks distinct integration parameters
in each contiguous stencil and returns the nine-component crossing state and
the degree used. The integration parameter is affine parameter for direct EOM
and coordinate time for normalized EOM; state component zero stores the other
quantity. The plane handlers compute local coordinates and apply the terminal
plane's radial bounds. An accepted terminal crossing stops the photon; a
nonterminal crossing does not.

For numerical batches, the manager applies the coordinate-time slot limit on
the accepted state before resolving pending crossings. This preserves a
crossing when that state reaches the last allowed time slot. The host time-slot
manager then removes stopped photons from the active list.

Batch integrators shift accepted-state history and write separate sparse
`light_blueprint_non_terminal_crossings_XX_YY.bin` and
`light_blueprint_terminal_crossings_XX_YY.bin` files. Each native record holds
the tile-local photon index, polynomial degree, crossing integration parameter,
two local plane coordinates, and nine interpolated state components. The image
blueprint also stores the axial angular momentum, coordinate time, and signed
plane distance at the first accepted state after a nonterminal crossing.
Single-photon integrators write
`plane_crossings.txt` with plane type, coordinate time, affine parameter, local
coordinates, nine state components, and `interpolation_degree`.

Claim evidence:
- Claim: The event manager delays detected plane crossings until a centered stencil or physical stop is available, selects the highest usable generated degree, and stores full crossing states and actual degrees in separate batch files or `plane_crossings.txt` for single photons.
- Role: descriptive behavior
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py` — `event_detection_manager_kernel`
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_event_time_and_state.py` — stencil selection and crossing state; `single_integrator_analytical.py` and `single_integrator_numerical.py` — text output; `batch_integrator_analytical.py` and `batch_integrator_numerical.py` — sparse crossing files; `handle_non_terminal_plane_intersection.py` and `handle_terminal_plane_intersection.py` — local coordinates and radius filtering

Blueprint headers carry tile identity/counts, `alpha_w`, `alpha_h`, and binary-layout
version 7. Each native same-build record is 124 bytes. Records carry nonterminal
crossing coordinates, affine parameter, and coordinate time; terminal-plane
texture coordinates; final angles, termination affine parameter and coordinate
time; normalized image-sample fractions; and the first accepted state's axial
angular momentum, coordinate time, and signed distance after a nonterminal
crossing. Final records admit spatial- and
temporal-interpolation failure statuses in addition to the existing physical
stops and numerical failures; internal `ACTIVE` and `REJECTED` states remain
invalid serialized outcomes.
The renderer places rays from the normalized fractions and preserves the
vertical raster flip, while plane diagnostics remain available to
`blueprint_analysis.py`. That analysis also consumes optional matching
normalization sidecars and plots their magnitudes by termination status.

Claim evidence:
- Claim: Blueprint binary-layout version 7 uses 124-byte records for terminal texture coordinates, nonterminal crossing and post-step diagnostics, final affine parameter and coordinate time, final angles, and normalized image-sample fractions; the renderer preserves the documented vertical raster flip.
- Role: generated evidence
- Deciding authority: `nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py` — binary-layout constants and dtype
- Corroboration: `nrpy/examples/geodesic_visualizations/visualize_lensed_image.py` — image placement

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
Callers may supply trial-locked native spatial-stencil centers; when provided,
the wrapper reuses those center indices for every temporal-stencil payload and
endpoint refresh. Normalized-EOM calls also provide the integration parameter,
trial step, and RK stage so interpolated geometry uses the RK stage coordinate
time; direct-EOM calls continue to read coordinate time from the state.

The numerical batch integrator can enable `perform_synthetic_slice_check` to
record synthetic temporal-stencil nodes whose times lie strictly below `0` or
above `t_numerical_end`. Direct first/final-slice endpoint dispatch does not
count, and normalization or `L_z` diagnostic interpolations do not count. The
wrapper records the first successful RKF45 interpolation request time for each
photon and boundary. It writes one separate 26-byte tile sidecar record per
photon: a `uint64_t` photon index, two `uint8_t` lower/upper-use flags, then
lower and upper `f64` request times. A false flag has a `NaN` time. Records
have no file header and use native byte order.
The numerical single-photon integrator uses the same wrapper check for RKF45
stage interpolations and prints its lower/upper flags and first request times
once the run finishes; it does not write a sidecar file.

The spatial helper performs tensor-product Lagrange interpolation only in the
two non-azimuthal native coordinates. It obtains spatial metric derivatives
from analytic derivatives of the same uniform-grid basis, rotates tensors from
the stored reference-azimuth plane to the target azimuth, and transforms native
derivatives to Cartesian coordinates with the reference-metric Jacobian. The
temporal helper builds a barycentric Lagrange basis from the actual, potentially
nonuniform slice times. For `g4DD`, analytic derivatives of that temporal basis
provide metric time derivatives; `g4DD_d0` and `GammaUDD` instead interpolate
their stored secondary payloads.

Ray-local interpolation errors have distinct terminal statuses. Coordinate
inversion, spatial-helper, and nonfinite spatial-output failures produce
`FAILURE_SPATIAL_INTERPOLATION`; temporal-window, temporal-helper, and
nonfinite final temporal-output failures produce
`FAILURE_TEMPORAL_INTERPOLATION`. The wrapper fills that ray's interpolation
scratch outputs with `NAN`. RKF45 finalization and event detection preserve the
status and last accepted persistent state, after which the existing host
routing serializes the ray as completed and continues the batch. Observer
interpolation remains fatal because no rays can be initialized without the
observer metric. Optional terminal and nonterminal normalization diagnostics
use temporary statuses, skip failed diagnostic samples, and leave their
sidecar entries as `NAN` without replacing a physical termination status.

Claim evidence:
- Claim: Numerical interpolation uses authoritative combined-container slice times, preserves the full logical grid including ghost zones, performs two-dimensional native spatial Lagrange interpolation with analytic basis derivatives and azimuthal tensor rotation, uses a nonuniform temporal barycentric basis, reuses caller-supplied trial-locked spatial centers, supplies RK stage coordinate time for normalized EOM, classifies ray-local spatial and temporal interpolation failures separately, and uses nominal spacing only for approximate synthetic temporal-stencil edge times.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py`, `temporal_lagrange_interpolation.py`, and `time_window_manager_numerical.py` — spatial interpolation, temporal interpolation, and time-window behavior
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py` — combined layout and metadata; `nrpy/infrastructures/BHaH/interpolation/differentiate_interpolation_lagrange_uniform.h` — uniform-basis derivative helper

Claim evidence:
- Claim: Optional batch tracking records the first successful RKF45 request that uses synthetic stencil nodes outside `[0, t_numerical_end]`; direct endpoint and diagnostic interpolations are excluded, and each tile writes fixed-order per-photon records with `NaN` times for false flags.
- Role: runtime diagnostic output
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py` and `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py` — stencil classification, per-photon storage, and binary output
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py` — optional parameter, per-tile filename, and runtime telemetry

Claim evidence:
- Claim: The numerical single-photon integrator passes tracking storage only to RKF45 stage interpolations and prints each boundary flag and first request time during cleanup when `perform_synthetic_slice_check` is enabled; diagnostic interpolations pass null tracking pointers.
- Role: runtime diagnostic output
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py` — single-photon tracker storage, evolution calls, and terminal reporting
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py` — synthetic-node detection and first-request recording

Numerical endpoint dispatch is piecewise constant. At or below the first stored
time, the first slice is spatially interpolated; at or above the selected final
slice, that final slice is reused. `g4DD` obtains temporal metric derivatives
from temporal interpolation. `g4DD_d0` uses stored metric derivatives on real
slices, but zeroes their temporal slots for direct endpoint dispatch and for
synthetic lower or upper nodes used to extend a frozen endpoint; the metric and
spatial-derivative slots remain copied from that endpoint. `GammaUDD` reuses
stored Christoffels at endpoints. The `g4DD_d0` padding policy is not controlled
by `--raytracing-static-christoffels`; that flag changes the exported final
Christoffel payload only.

Claim evidence:
- Claim: Numerical endpoint dispatch reuses the first/final stored slices as piecewise-constant continuations; frozen `g4DD_d0` nodes retain endpoint metric and spatial-derivative data but zero metric time derivatives, while `--raytracing-static-christoffels` affects only exported GammaUDD data.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py` — endpoint dispatch; `output_raytracing_data.py` — payload mode
- Corroboration: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/temporal_lagrange_interpolation.py` — temporal endpoint handling

Startup reads `.bin` `NGHOSTS` metadata and terminates if the configured spatial
half-width requires more ghost zones than available. It also terminates when
observer or positive Cartesian `r_escape` probes cannot support the requested
native stencil; the relevant `.par` controls are
`observer_x/y/z`, `r_escape`, and
`numerical_spacetime_spatial_interp_half_width`.

Claim evidence:
- Claim: Numerical startup rejects insufficient ghost-zone capacity and observer or positive Cartesian `r_escape` probes that cannot support the requested native stencil.
- Role: public/scientific contract
- Deciding authority: `nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py` — startup domain validation
- Corroboration: `nrpy/infrastructures/BHaH/xx_tofrom_Cart.py` — Cartesian-to-native inversion and admission checks

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
- Role: generated evidence
- Deciding authority: `nrpy/infrastructures/BHaH/diagnostics/diagnostics.py` and `output_raytracing_data.py` — diagnostics registration and exporter
- Corroboration: `nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py` — stage-1 validation and combined container

## Sources

- [main_batch.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py) - `main`
- [batch_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py) - `batch_integrator_numerical`
- [batch_integrator_analytical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_analytical.py) - `batch_integrator_analytical`
- [single_integrator_analytical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py) - `single_integrator_analytical`
- [single_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py) - `single_integrator_numerical`
- [normal_observer_log_energy.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/normal_observer_log_energy.py) - `normal_observer_log_energy`
- [set_initial_conditions_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py) - `set_initial_conditions_kernel`
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
- [temporal_lagrange_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/temporal_lagrange_interpolation.py) - `register_CFunction_temporal_lagrange_interpolation`
- [differentiate_interpolation_lagrange_uniform.h](../../../nrpy/infrastructures/BHaH/interpolation/differentiate_interpolation_lagrange_uniform.h) - `compute_lagrange_basis_derivative_coeffs_xi`
- [output_raytracing_data.py](../../../nrpy/infrastructures/BHaH/diagnostics/output_raytracing_data.py) - `register_CFunction_output_raytracing_data`, `raytracing_data_point_index_from_logical_indices`
- [combine_raytracing_time_slices.py](../../../nrpy/infrastructures/BHaH/diagnostics/combine_raytracing_time_slices.py) - `InputSliceInfo`, `Layout`
- [diagnostics.py](../../../nrpy/infrastructures/BHaH/diagnostics/diagnostics.py) - `register_all_diagnostics`, `enable_raytracing_data_output`

## See Also

- [BHaH](index.md)
- [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)
- [Geodesics](../../equations/general-relativity/geodesics.md)
- [Black Hole Evolution](../../examples/black-hole-evolution.md)

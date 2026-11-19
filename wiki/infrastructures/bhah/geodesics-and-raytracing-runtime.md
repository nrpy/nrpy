# Geodesics And Raytracing Runtime

> Explain BHaH geodesic integration, numerical spacetime interpolation, and photon diagnostics. · Status: confirmed
> Up: [BHaH](index.md)

## Summary

BHaH generates standalone massive and photon geodesic integrators. Massive particles use GSL. The photon integrators use split RKF45 kernels with analytic or numerical metrics. Numerical evolution reads combined spacetime slices and reuses a fixed interpolation stencil and spatial mask through each RKF45 attempt.

## Detail

### Observer And Photon Evolution

The shared photon initializer constructs a metric-orthonormal observer tetrad and sets the initial camera-frame energy magnitude to one. It converts the resulting direct four-momentum to normalized variables when requested. The [geodesic equations](../../equations/general-relativity/geodesics.md) define both state vectors and null constraints.

Numerical batches group active photons into coordinate-time slots and process chunks with separate arrays for each state component. The RKF45 controller accepts or rejects each photon individually. On acceptance, the direct EOM refreshes the metric with the attempt's interpolation selection and evaluates `ln|alpha p^0|`. The normalized EOM reads the same quantity from `u`. The evolution-measure stop uses an upper bound on this quantity. A separate `rkf45_max_consecutive_min_steps` limit counts attempts at the configured minimum step, including accepted attempts. The retry limit counts consecutive rejections.

The optional batch norm files express the residual on the direct momentum scale. They store `|C_direct|` for direct evolution or `|exp(2u)(C_normalized - 1)|` for normalized evolution. Interpolation failures during delayed diagnostics skip the affected measurement. The terminal progress display counts retryable rejections by RKF45 error, spatial interpolation failure, and temporal interpolation failure for the latest chunk.

Evidence (descriptive behavior): `set_initial_conditions_kernel.py::set_initial_conditions_kernel`, `batch_integrator_numerical.py::batch_integrator_numerical`, and `rkf45_finalize_and_control_kernel.py::rkf45_finalize_and_control_kernel`.

The generated massive-particle GSL wrapper evaluates the metric and Christoffels at the current state before calling the eight-component RHS.

Evidence (descriptive behavior): `ode_gsl_wrapper_massive.py`; corroboration: `calculate_ode_rhs_massive.py`.

### Numerical-Spacetime Interpolation

The time-window manager maps combined spacetime slices for the current scheduling slot, including the interpolation halo and backward lookahead. Stored times specify the actual temporal nodes. The nominal spacing contributes to the spacing of synthetic endpoint nodes. See [Raytracing Export](diagnostics-output-and-checkpointing.md#raytracing-export-and-progress-indicator) for the three stored data formats.

The reader uses the numerical file's ghost-zone count to index the stored grid. Startup rejects a spatial half-width that exceeds that count. The batch startup checks that the observer and escape-radius probes support the requested stencil.

On each selected slice, spatial Lagrange interpolation uses the nonazimuthal native coordinates of a `SinhCylindricalv2n2` grid with two stored azimuthal planes separated by pi. The stencil reads the full stored grid, including ghost records. It accumulates the Cartesian components stored at those nodes, then rotates the reconstructed tensors from the selected reference plane to the target azimuth. Metric interpolation differentiates the same spatial polynomial and azimuthal rotation and converts the native derivatives to Cartesian derivatives using the inverse reference-metric Jacobian.

Evidence (descriptive behavior): `azimuthal_symmetry_spatial_lagrange_interpolation.py::register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation`; corroboration: `time_window_manager_numerical.py::time_window_manager_numerical`.

The wrapper then interpolates in time. For `g4DD`, it differentiates the temporal metric polynomial. For `g4DD_d0`, it interpolates the stored coordinate-time derivatives. For `GammaUDD`, it interpolates the stored Christoffels. Outside the retained time range it samples the nearest endpoint. The `g4DD_d0` path sets the coordinate-time derivative slots to zero for frozen endpoint evaluations and synthetic endpoint nodes. The `GammaUDD` path reuses the stored endpoint Christoffels.

Evidence (descriptive behavior): `numerical_interpolation.py::register_CFunction_numerical_interpolation`; corroboration: `temporal_lagrange_interpolation.py::register_CFunction_temporal_lagrange_interpolation`.

With `perform_synthetic_slice_check` enabled, the numerical integrators record the first RKF45 geometry request using a synthetic node on either time boundary. The batch output is `light_blueprint_synthetic_slice_usage_XX_YY.bin`.

Evidence (descriptive behavior): `numerical_interpolation.py::register_CFunction_numerical_interpolation` and `main_batch.py::main`.

### Fixed Stencils And The Temporal Step Limit

The numerical integrators choose the native spatial center at the start of an attempt. They choose the retained temporal node nearest the predicted midpoint, `t + h/2` in normalized evolution or `t + h*p^0/2` in direct evolution. The stages evaluate geometry at their own coordinates using those centers. A retry chooses the centers again after changing the step size.

For temporal half-width `n > 0`, the runtime coordinate-time cap `H_t` is just below `(2*n-1)*delta_safe`, where `delta_safe` is the minimum spacing used by the retained and synthetic nodes. Normalized evolution caps `|h|` by `H_t`. Direct evolution caps it by `0.9*H_t/|p^0|` using the accepted start momentum. This predicts the time span when `p^0` varies. The initial and subsequent steps also obey `rkf45_h_max` and the mapped-window limit. A required cap below `rkf45_h_min` is reported as a failure. Temporal half-width zero uses one node and disables the temporal travel cap.

Evidence (descriptive behavior): `time_window_manager_numerical.py::time_window_manager_numerical_temporal_step_limit` and `rkf45_finalize_and_control_kernel.py::rkf45_cap_numerical_step`.

### Spatial Node Exclusion For Numerical Photons

On the stored azimuthal plane selected at the start of an RKF45 attempt, the mask compares each spatial metric's Frobenius norm with the median norm of its neighboring stencil nodes multiplied by `numerical_spacetime_spatial_mask_threshold`. A norm strictly above this product on any selected time slice excludes that spatial node across the temporal stencil. Synthetic time nodes use the corresponding endpoint slice. Spatial half-width zero leaves the mask empty.

Evidence (descriptive behavior): `time_window_manager_numerical.py::time_window_manager_numerical_build_spatial_mask`.

With any exclusions, the spatial interpolator forms Lagrange weights along rho within each retained row, then along z across rows. Each row needs at least two retained nodes, and the z interpolation needs at least two usable rows. The derivatives come from that same reconstruction. An empty mask selects ordinary spatial interpolation. The temporal interpolation method stays the same for either reconstruction.

Evidence (descriptive behavior): `azimuthal_symmetry_spatial_lagrange_interpolation_masked.py::spatial_lagrange_masked_weights`.

For positive spatial half-width, the numerical batch integrator prepares outlier flags for selected time slices as it maps them and reuses flags for slices shared by consecutive mappings. Each node on the two stored azimuthal planes has a flag for each possible stencil position at that node, such as a corner, an edge, or the interior. During an attempt, the mask reads the flag for each node's stencil position and combines the exclusions across the selected slices. The single-photon integrator computes the norms and medians directly on its selected stencil.

Evidence (descriptive behavior): `time_window_manager_numerical.py::time_window_manager_numerical_prepare_slice_outlier_flags`, `time_window_manager_numerical.py::time_window_manager_numerical_mmap_for_slot`, and `time_window_manager_numerical.py::time_window_manager_numerical_build_spatial_mask`.

Single and batch integrations retain the mask through the stages and accepted-step metric evaluations. The batch integrator rebuilds the mask for delayed diagnostics from the saved accepted stencil selection or the saved plane-crossing candidate. Numerical Brent refinement uses one selection for the full crossing bracket. A change of stencil or mask between attempts can change the reconstructed metric.

Evidence (descriptive behavior): `single_integrator_numerical.py::single_integrator_numerical`, `batch_integrator_numerical.py::batch_integrator_numerical`, and `find_plane_crossings_brent.py::plane_crossing_brent_refine_candidate`.

### Fixed-Time Geometry Sampling

With fast light, the numerical integrators map geometry at `t_start` while coordinate time and the integration parameter continue to evolve. Interior metric derivatives or Christoffels are evaluated at that sampling time. The same sampling time is used for refinement and diagnostics. The temporal travel cap is disabled.

Evidence (descriptive behavior): `numerical_interpolation.py::register_CFunction_numerical_interpolation`, `batch_integrator_numerical.py::batch_integrator_numerical`, and `single_integrator_numerical.py::single_integrator_numerical`.

### Metric-Gradient Tolerances

The `g4DD` and `g4DD_d0` modes can scale the absolute and relative RKF45 tolerances. Let `G` be the maximum absolute Cartesian metric derivative over successful stages of the attempted step. The multiplier is `max(tol_scale_min, 1/(1+L_sense*G))`. This measures coordinate-dependent variation. The normalized log-energy tolerance and the Brent tolerances are separate. With `tol_scale_min = 1` or `L_sense = 0`, the integrators skip the derivative scan.

Evidence (descriptive behavior): `rkf45_finalize_and_control_kernel.py::rkf45_finalize_and_control_kernel`, `batch_integrator_numerical.py::batch_integrator_numerical`, and `single_integrator_numerical.py::single_integrator_numerical`.

### Planes And Apparent Horizons

Plane side changes on accepted steps save the starting state and integration-parameter interval. After evolution, Brent refinement evaluates fifth-order RKF45 trials from that saved state. The first refined terminal crossing inside the specified radial bounds supplies the reported terminal event. A nonterminal crossing is kept if it precedes the selected terminal event.

Evidence (descriptive behavior): `event_detection_manager_kernel.py::event_detection_manager_kernel` and `find_plane_crossings_brent.py::find_plane_crossings_brent`.

If an apparent-horizon `.bin` is supplied, the numerical program validates its slice table against the spacetime data. It checks the observer state and accepted photon states using the photon coordinate time or `t_start` in fast light. A nonfinite lookup time returns `HORIZON_NOT_CHECKED`. A finite negative lookup time is set to zero before comparison with the retained slice times. A lookup before the first retained slice or after the last returns `HORIZON_NOT_CHECKED`.

Evidence (descriptive behavior): `puncture_diagnostics.py::apparent_horizons_open`, `puncture_diagnostics.py::apparent_horizons_classify`, and `event_detection_manager_kernel.py::event_detection_manager_kernel`.

The lookup checks the individual horizons before common-horizon detection and the common horizon afterward. For each active horizon, it subtracts the angularly interpolated surface radius from the photon radius measured from the saved surface origin. Between retained slices, it interpolates the signed radial gap and containment tolerance linearly in coordinate time. An exact slice match uses that slice alone. A horizon check with missing bracketing data is unavailable. A containing surface sets `STOP_CONDITION_INSIDE_APPARENT_HORIZON`. A refined earlier terminal plane can still supply the reported terminal event.

Evidence (descriptive behavior): `puncture_diagnostics.py::apparent_horizons_classify` and `event_detection_manager_kernel.py::event_detection_manager_kernel`.

The text horizon diagnostics classify the final accepted state in single-photon runs and evolution-measure stops in batch runs when binary horizon stopping is disabled.

Evidence (descriptive behavior): `puncture_diagnostics.py::puncture_diagnostics` and `single_integrator_numerical.py::single_integrator_numerical`.

When either horizon input is enabled, the numerical batch integrator writes a puncture ID byte, a horizon status byte, and a double-precision radial gap for each photon. With a combined apparent-horizon input file, it writes the saved horizon check and prints the record count. For directory diagnostics alone, it prints the fraction inside a horizon among evolution-measure stops with a containment result.

Evidence (descriptive behavior): `puncture_diagnostics.py::write_puncture_diagnostics` and `batch_integrator_numerical.py::batch_integrator_numerical`.

## Sources

- [g4DD_metric.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/g4DD_metric.py) - g4DD_metric
- [ode_gsl_wrapper_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/ode_gsl_wrapper_massive.py) - ode_gsl_wrapper_massive
- [calculate_ode_rhs_massive.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/calculate_ode_rhs_massive.py) - calculate_ode_rhs_massive
- [time_window_manager_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/time_window_manager_numerical.py) - time_window_manager_numerical and its temporal cap and mask helpers
- [azimuthal_symmetry_spatial_lagrange_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation.py) - spatial Lagrange interpolation
- [azimuthal_symmetry_spatial_lagrange_interpolation_masked.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/azimuthal_symmetry_spatial_lagrange_interpolation_masked.py) - spatial_lagrange_masked_weights
- [numerical_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/numerical_interpolation.py) - register_CFunction_numerical_interpolation
- [temporal_lagrange_interpolation.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/interpolation/temporal_lagrange_interpolation.py) - register_CFunction_temporal_lagrange_interpolation
- [set_initial_conditions_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py) - set_initial_conditions_kernel
- [main_batch.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py) - `main`, optional batch output files
- [batch_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/batch_integrator_numerical.py) - batch_integrator_numerical
- [single_integrator_numerical.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_numerical.py) - single_integrator_numerical
- [rkf45_finalize_and_control_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/rkf45_finalize_and_control_kernel.py) - rkf45_finalize_and_control_kernel
- [event_detection_manager_kernel.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/event_detection_manager_kernel.py) - event_detection_manager_kernel
- [find_plane_crossings_brent.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/find_plane_crossings_brent.py) - find_plane_crossings_brent
- [puncture_diagnostics.py](../../../nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/puncture_diagnostics.py) - horizon lookup, containment checks, and `write_puncture_diagnostics`

## See Also

- Depends on: [Geodesic Equations](../../equations/general-relativity/geodesics.md)
- Example: [Geodesic Raytracing](../../examples/geodesic-raytracing.md)
- See also: [Diagnostics Output And Checkpointing](diagnostics-output-and-checkpointing.md)

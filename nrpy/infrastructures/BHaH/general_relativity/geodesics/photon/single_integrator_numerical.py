"""
Generate a standalone CPU project for one photon in a numerical spacetime.

The generated C executable follows the numerical batch integrator's numerical
interpolation, time-window management, initial-condition geometry, and RKF45
control logic, but removes batching, double buffering, device memory, streams,
and blueprint output. It writes one trajectory row after every accepted RKF45
step, including a signed normalization diagnostic. Optional RKF45 debugging
writes trial-level and stage-level records from this single-photon host loop.
Optional terminal and nonterminal planes use the shared event-time interpolation
and write crossing times, local plane coordinates, and all nine interpolated
state components to ``plane_crossings.txt``.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

# pylint: disable=too-many-lines

import os
import shutil

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    main_single,
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation import (
    numerical_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_non_terminal_plane_intersection,
    handle_terminal_plane_intersection,
    normal_observer_log_energy,
    photon_momentum_to_normalized_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.set_initial_conditions_kernel import (
    register_photon_batch_structs,
    set_initial_conditions_kernel,
)


def single_integrator_numerical(  # pylint: disable=invalid-name,too-many-locals
    spacetime_name: str,
    dataset_coord_system: str,
    interpolation_method: str = "g4DD",
    maximum_degree: int = 4,
    normalized_eom: bool = False,
    enable_rkf45_trial_debug: bool = False,
    axisymmetric_about_z: bool = False,
    track_synthetic_slice_usage: bool = False,
) -> None:
    """
    Register the standalone numerical-spacetime single-photon C integrator.

    :param spacetime_name: Spacetime identifier used to select photon equations.
    :param dataset_coord_system: Coordinate system used by the numerical dataset.
    :param interpolation_method: Numerical geometry payload method used by the generated project.
    :param maximum_degree: Largest generated plane interpolation polynomial degree.
    :param normalized_eom: Whether to evolve normalized coordinate-time equations.
    :param enable_rkf45_trial_debug: Whether to write one diagnostic row for every
        RKF45 trial to ``rkf45_trials.txt`` and one row for each of its six
        stages to ``rkf45_stages.txt``. When axial symmetry is enabled, this
        also adds accepted-state $L_z$ values to ``trajectory.txt``.
    :param axisymmetric_about_z: Whether the numerical spacetime has rotational
        symmetry about the Cartesian ``z`` axis, allowing $L_z$ diagnostics.
    :param track_synthetic_slice_usage: Whether to enable reporting the first
        RKF45 request that uses synthetic temporal-stencil nodes outside the
        numerical-spacetime time range. The numerical interpolation function
        must be registered with the same option.
    :raises ValueError: If the interpolation method or dataset coordinate system
        is unsupported.

    Doctests:
    >>> import os
    >>> import tempfile
    >>> from unittest.mock import patch
    >>> import nrpy.c_function as cfc
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     single_integrator_numerical("Numerical", "SinhCylindricalv2n2")
    ...     generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "# lambda t x y z p^t p^x p^y p^z L_normal norm" in generated
    True
    >>> "const double trajectory_norm = normalization.C;" in generated
    True
    >>> "fabs(normalization.C)" not in generated
    True
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     single_integrator_numerical("Numerical", "SinhCylindricalv2n2", normalized_eom=True)
    ...     generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> ("# lambda t x y z u Pi_1 Pi_2 Pi_3 L_normal norm" in generated and
    ...  "const double trajectory_norm = normalization.C - 1.0;" in generated)
    True
    >>> "fabs(normalization.C - 1.0)" not in generated
    True
    >>> "trial_debug_file" not in generated
    True
    >>> "stage_debug_file" not in generated
    True
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     single_integrator_numerical(
    ...         "Numerical", "SinhCylindricalv2n2", enable_rkf45_trial_debug=True
    ...     )
    ...     generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "rkf45_trials.txt" in generated
    True
    >>> "trial_debug_file" in generated
    True
    >>> "rkf45_stages.txt" in generated
    True
    >>> "stage_debug_file" in generated
    True
    >>> "commondata.tiles_width = 1;" not in generated
    True
    >>> "commondata.tile_index_width = 0;" not in generated
    True
    >>> par.glb_code_params_dict["tiles_width"].add_to_parfile
    True
    >>> par.glb_code_params_dict["tile_index_width"].add_to_parfile
    True
    >>> "Cart_to_xx_and_nearest_i0i1i2_assume_valid" in generated
    True
    >>> "k_bundle[(stage - 1) * 9 + 5]" in generated
    True
    >>> with tempfile.TemporaryDirectory() as cache_dir, patch.dict(os.environ, {"XDG_CACHE_HOME": cache_dir}):
    ...     cfc.CFunction_dict.clear()
    ...     single_integrator_numerical(
    ...         "Numerical", "SinhCylindricalv2n2", track_synthetic_slice_usage=True
    ...     )
    ...     generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "Synthetic temporal-stencil nodes outside [0, t_numerical_end]:" in generated
    True
    >>> "used_lower_endpoint: %s, photon_request_time: %.15e" in generated
    True
    >>> par.glb_code_params_dict["perform_synthetic_slice_check"].add_to_parfile
    True
    """
    if interpolation_method not in ("g4DD", "g4DD_d0", "GammaUDD"):
        raise ValueError(
            "interpolation_method must be one of ('g4DD', 'g4DD_d0', 'GammaUDD'); "
            f"found '{interpolation_method}'."
        )
    if maximum_degree < 3:
        raise ValueError("Plane interpolation degree must be at least three")
    if dataset_coord_system != "SinhCylindricalv2n2":
        raise ValueError(
            "single_integrator_numerical supports only "
            "dataset_coord_system='SinhCylindricalv2n2'; "
            f"found '{dataset_coord_system}'."
        )
    if not isinstance(enable_rkf45_trial_debug, bool):
        raise ValueError(
            "enable_rkf45_trial_debug must be a bool, got "
            f"{type(enable_rkf45_trial_debug).__name__}."
        )
    if not isinstance(axisymmetric_about_z, bool):
        raise ValueError(
            f"axisymmetric_about_z must be a bool, got {type(axisymmetric_about_z).__name__}."
        )
    if not isinstance(track_synthetic_slice_usage, bool):
        raise ValueError(
            "track_synthetic_slice_usage must be a bool, got "
            f"{type(track_synthetic_slice_usage).__name__}."
        )
    enable_accepted_Lz_diagnostic = enable_rkf45_trial_debug and axisymmetric_about_z
    phi_dim = 1

    # Register the shared batch state definitions and metric-tetrad initializer.
    register_photon_batch_structs()
    BHaH_defines_h.register_BHaH_defines(
        "single_photon_macros",
        """
    #undef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    """,
    )
    set_initial_conditions_kernel(
        normalized_eom=normalized_eom,
        initialize_event_history=True,
        tile_indices_add_to_parfile=True,
    )

    # The shared initializer uses the batch tile-sampling rules. Expose the tile
    # counts and active tile indices so one single-ray process can reproduce
    # any batch-camera sample exactly. The defaults remain the center ray of
    # one tile. The shared initializer registers the tile indices.
    par.register_CodeParameters(
        "int",
        __name__,
        [
            "tiles_width",
            "tiles_height",
        ],
        [1, 1],
        commondata=True,
        add_to_parfile=True,
    )

    # Step 1: Register single-photon escape and numerical-dataset parameters.
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["r_escape"],
        [150.0],
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameter(
        "char[4096]",
        __name__,
        "numerical_spacetime_bin_path",
        "",
        commondata=True,
        add_to_parfile=True,
        description=(
            "Path to the validated combined numerical raytracing .bin file used by "
            "the numerical single-photon integrator."
        ),
    )
    par.register_CodeParameter(
        "bool",
        __name__,
        "perform_normalization_check",
        False,
        commondata=True,
        add_to_parfile=True,
    )
    if track_synthetic_slice_usage:
        par.register_CodeParameter(
            "bool",
            __name__,
            "perform_synthetic_slice_check",
            False,
            commondata=True,
            add_to_parfile=True,
            description=(
                "Record first RKF45 photon request times that use synthetic "
                "temporal-stencil nodes below t=0 or above t_numerical_end, "
                "then print the two flags and request times after integration."
            ),
        )
    trial_spatial_center_setup = r"""
    NumericalSpatialStencilCenter trial_spatial_center;
    const REAL trial_cartesian[3] = {
        (REAL)f_start[1], (REAL)f_start[2], (REAL)f_start[3]};
    REAL trial_native[3];
    int trial_automatic_center_idx[3];
    int trial_selected_center_idx[3];
    if (time_window_manager_numerical_resolve_spatial_target_and_stencil(
            &numerical_params, trial_cartesian, NULL, trial_native,
            trial_automatic_center_idx, trial_selected_center_idx) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
      // Preserve the last accepted state and route this ray through the same
      // per-photon spatial-interpolation failure path as the wrapper.
      *status = FAILURE_SPATIAL_INTERPOLATION;
      trial_spatial_center.i0 = 0;
      trial_spatial_center.i2 = 0;
    } else {
      trial_spatial_center.i0 = trial_selected_center_idx[0];
      trial_spatial_center.i2 = trial_selected_center_idx[2];
    } // END ELSE: RKF45 trial spatial center resolved
"""

    # Step 2: Select emitted C expressions for the two state conventions.
    if normalized_eom:
        # Normalized state stores affine parameter in f[0].  Its origin is
        # always lambda=0; coordinate time remains the integration parameter.
        initial_state_time = "0.0"
        initial_integration_value = "commondata.t_start"
        coordinate_time_expression = "*integration_param"
        trajectory_lambda_expression = "f[0]"
        trajectory_time_expression = "*integration_param"
        trajectory_header = "# lambda t x y z u Pi_1 Pi_2 Pi_3 L_normal norm\\n"
        interpolation_stage_arguments = (
            "&trial_spatial_center.i0, &trial_spatial_center.i2, "
            "integration_param, h, stage,"
        )
        interpolation_initial_arguments = "NULL, NULL, integration_param, h, 1,"
        rhs_integration_arguments = "integration_param, h,"
        momentum_conversion_call = (
            "photon_momentum_to_normalized_kernel(f, metric, chunk_size);"
        )
        normalization_kernel_name = "normalization_constraint_photon_normalized"
        normalization_diagnostic_expression = "normalization.C - 1.0"
    else:
        # Direct state stores coordinate time in f[0].  Affine integration
        # parameter starts at lambda=0.
        initial_state_time = "commondata.t_start"
        initial_integration_value = "0.0"
        coordinate_time_expression = "f[0]"
        trajectory_lambda_expression = "*integration_param"
        trajectory_time_expression = "f[0]"
        trajectory_header = "# lambda t x y z p^t p^x p^y p^z L_normal norm\\n"
        interpolation_stage_arguments = (
            "&trial_spatial_center.i0, &trial_spatial_center.i2,"
        )
        interpolation_initial_arguments = "NULL, NULL,"
        rhs_integration_arguments = ""
        momentum_conversion_call = ""
        normalization_kernel_name = "normalization_constraint_photon"
        normalization_diagnostic_expression = "normalization.C"

    event_state_columns = (
        "interpolated_lambda interpolated_x interpolated_y interpolated_z "
        "interpolated_u interpolated_Pi_1 interpolated_Pi_2 "
        "interpolated_Pi_3 interpolated_L_normal"
        if normalized_eom
        else "interpolated_t interpolated_x interpolated_y interpolated_z "
        "interpolated_p^t interpolated_p^x interpolated_p^y "
        "interpolated_p^z interpolated_L_normal"
    )
    event_state_format = " ".join(["%.17e"] * 9)
    non_terminal_event_state_arguments = ", ".join(
        f"non_terminal_plane_event_state[{component}]" for component in range(9)
    )
    terminal_event_state_arguments = ", ".join(
        f"terminal_plane_event_state[{component}]" for component in range(9)
    )

    if enable_accepted_Lz_diagnostic:
        trajectory_header = trajectory_header[:-2] + " L_z\\n"
        axial_angular_momentum_metric_argument = "" if normalized_eom else "metric, "
        initial_Lz_evaluation = f"""
  axial_angular_momentum_t initial_angular_momentum;
  axial_angular_momentum_z(
      f,
      {axial_angular_momentum_metric_argument}&initial_angular_momentum,
      chunk_size,
      stream_idx);
  const double initial_Lz = initial_angular_momentum.Lz;
  if (!isfinite(initial_Lz)) {{
    fprintf(stderr, "ERROR: initial axial angular momentum was not finite.\\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial L_z diagnostic was invalid
  printf("Initial axial angular momentum L_z = %.17e\\n", initial_Lz);
"""
        accepted_Lz_evaluation = f"""
      axial_angular_momentum_t accepted_angular_momentum;
      axial_angular_momentum_z(
          f,
          {axial_angular_momentum_metric_argument}&accepted_angular_momentum,
          chunk_size,
          stream_idx);
      const double accepted_Lz = accepted_angular_momentum.Lz;
      if (!isfinite(accepted_Lz)) {{
        fprintf(
            stderr,
            "ERROR: accepted-state axial angular momentum was not finite.\\n");
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: accepted-state L_z diagnostic was invalid
"""
        accepted_Lz_before_interpolation = (
            accepted_Lz_evaluation if normalized_eom else ""
        )
        accepted_Lz_after_interpolation = (
            "" if normalized_eom else accepted_Lz_evaluation
        )
        accepted_Lz_format = " %.17e"
        accepted_Lz_argument = ", accepted_Lz"
        failed_interpolation_Lz_format = " %.17e"
        failed_interpolation_Lz_argument = (
            ", accepted_Lz" if normalized_eom else ", NAN"
        )
    else:
        initial_Lz_evaluation = ""
        accepted_Lz_evaluation = ""
        accepted_Lz_before_interpolation = ""
        accepted_Lz_after_interpolation = ""
        accepted_Lz_format = ""
        accepted_Lz_argument = ""
        failed_interpolation_Lz_format = ""
        failed_interpolation_Lz_argument = ""

    accepted_metric_interpolation_arguments = (
        interpolation_initial_arguments
        if normalized_eom
        else "&trial_spatial_center.i0, &trial_spatial_center.i2,"
    )
    synthetic_slice_usage_declarations = (
        r"""
  const long int single_photon_indices[1] = {0};
  synthetic_slice_usage_t synthetic_slice_usage_by_photon[1] = {
      {false, NAN, false, NAN}};
"""
        if track_synthetic_slice_usage
        else ""
    )
    synthetic_slice_usage_stage_arguments = (
        "single_photon_indices,\n          synthetic_slice_usage_by_photon,"
        if track_synthetic_slice_usage
        else ""
    )
    synthetic_slice_usage_null_arguments = (
        "NULL,\n          NULL," if track_synthetic_slice_usage else ""
    )
    synthetic_slice_usage_report = (
        r"""
  if (commondata.perform_synthetic_slice_check) {
    const synthetic_slice_usage_t *usage =
        &synthetic_slice_usage_by_photon[0];
    printf("Synthetic temporal-stencil nodes outside [0, t_numerical_end]:\n");
    printf(
        "  used_lower_endpoint: %s, photon_request_time: %.15e\n",
        usage->used_lower_endpoint ? "true" : "false",
        usage->lower_request_time);
    printf(
        "  used_upper_endpoint: %s, photon_request_time: %.15e\n",
        usage->used_upper_endpoint ? "true" : "false",
        usage->upper_request_time);
  } // END IF: synthetic temporal-stencil reporting enabled
"""
        if track_synthetic_slice_usage
        else ""
    )
    log_energy_evaluation = (
        "const double log_energy_measure = f[4];"
        if normalized_eom
        else "double log_energy_bundle[1];\n"
        "      normal_observer_log_energy(\n"
        "          f, metric, log_energy_bundle, chunk_size, stream_idx);\n"
        "      const double log_energy_measure = log_energy_bundle[0];"
    )

    stage_normalization_diagnostic_expression = (
        normalization_diagnostic_expression.replace(
            "normalization.", "stage_normalization."
        )
    )
    if normalized_eom:
        stage_time_expression = """*integration_param +
          rkf45_stage_time_fractions[stage - 1] * *h"""
    else:
        stage_time_expression = "f_temp[0]"

    if enable_rkf45_trial_debug:
        trial_component_names = (
            r"""
    "lambda",
    "x",
    "y",
    "z",
    "u",
    "Pi_1",
    "Pi_2",
    "Pi_3",
    "L_normal"
"""
            if normalized_eom
            else r"""
    "t",
    "x",
    "y",
    "z",
    "p^0",
    "p^1",
    "p^2",
    "p^3",
    "L_normal"
"""
        )
        trial_debug_declarations = r"""
  FILE *trial_debug_file = NULL;
  rkf45_trial_diagnostic_t trial_debug;
  const char *trial_component_names[] = {
{trial_component_names}
  };
"""
        trial_debug_declarations = trial_debug_declarations.replace(
            "{trial_component_names}", trial_component_names
        )
        if normalized_eom:
            stage_debug_header = """# accepted_step trial_number retry_number stage h_trial t_start stage_time lambda x y z r_stage xx0 xx1 xx2 i0 i1 i2 stencil_i0_low stencil_i0_high stencil_i2_low stencil_i2_high stage_norm_error u Pi_1 Pi_1_derivative
"""
        else:
            stage_debug_header = """# accepted_step trial_number retry_number stage h_trial t_start stage_time t x y z r_stage xx0 xx1 xx2 i0 i1 i2 stencil_i0_low stencil_i0_high stencil_i2_low stencil_i2_high stage_norm_error p^0 p^1 p^1_derivative
"""
        stage_debug_header_c = (
            stage_debug_header.rstrip("\n").replace("\\", "\\\\").replace('"', '\\"')
            + "\\n"
        )
        stage_debug_declarations = (
            r"""
  // Host-side stage logging preserves the sequential order of every trial.
  FILE *stage_debug_file = NULL;
  const double rkf45_stage_time_fractions[] = {
      0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};
"""
            if normalized_eom
            else r"""
  // Host-side stage logging preserves the sequential order of every trial.
  FILE *stage_debug_file = NULL;
"""
        )
        stage_debug_open = rf"""
  stage_debug_file = fopen("rkf45_stages.txt", "w");
  if (stage_debug_file == NULL) {{
    fprintf(stderr, "ERROR: could not open rkf45_stages.txt for writing.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: RKF45 stage diagnostics unavailable
  fprintf(stage_debug_file, "{stage_debug_header_c}");
"""
        stage_debug_record = rf"""
      // A failed interpolation intentionally leaves NaN RK scratch data. Do
      // not treat that expected sentinel as a separate diagnostics failure.
      if (*status != FAILURE_SPATIAL_INTERPOLATION &&
          *status != FAILURE_TEMPORAL_INTERPOLATION) {{
      // Record the interpolated stage before the next RKF45 stage update.
      normalization_constraint_t stage_normalization;
      {normalization_kernel_name}(
          f_temp,
          metric,
          &stage_normalization,
          chunk_size,
          stream_idx);
      const double stage_norm_error =
          {stage_normalization_diagnostic_expression};
      if (!isfinite(stage_norm_error)) {{
        fprintf(
            stderr,
            "ERROR: stage %d normalization error was not finite.\n",
            stage);
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: stage normalization error invalid

      const REAL stage_cartesian[3] = {{
          (REAL)f_temp[1], (REAL)f_temp[2], (REAL)f_temp[3]}};
      REAL stage_xx[3];
      int stage_indices[3];
      Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__SinhCylindricalv2n2(
          &numerical_params,
          stage_cartesian,
          stage_xx,
          stage_indices);
      const double stage_radius = sqrt(
          f_temp[1] * f_temp[1] +
          f_temp[2] * f_temp[2] +
          f_temp[3] * f_temp[3]);
      const double stage_time = {stage_time_expression};
      const double stage_pi1_derivative =
          k_bundle[(stage - 1) * 9 + 5];
      fprintf(
          stage_debug_file,
          "%ld %ld %d %d "
          "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
          "%d %d %d %d %d %d %d "
          "%.17g %.17g %.17g %.17g\n",
          accepted_step_before_trial,
          trial_number,
          retry_number_before,
          stage,
          h_trial,
          t_start,
          stage_time,
          f_temp[0],
          f_temp[1],
          f_temp[2],
          f_temp[3],
          stage_radius,
          (double)stage_xx[0],
          (double)stage_xx[1],
          (double)stage_xx[2],
          stage_indices[0],
          stage_indices[1],
          stage_indices[2],
          stage_indices[0] - commondata.numerical_spacetime_spatial_interp_half_width,
          stage_indices[0] + commondata.numerical_spacetime_spatial_interp_half_width,
          stage_indices[2] - commondata.numerical_spacetime_spatial_interp_half_width,
          stage_indices[2] + commondata.numerical_spacetime_spatial_interp_half_width,
          stage_norm_error,
          f_temp[4],
          f_temp[5],
          stage_pi1_derivative);
      fflush(stage_debug_file);
      }} // END IF: stage interpolation valid
"""
        trial_debug_open = r"""
  trial_debug_file = fopen("rkf45_trials.txt", "w");
  if (trial_debug_file == NULL) {
    fprintf(stderr, "ERROR: could not open rkf45_trials.txt for writing.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  } // END IF: RKF45 trial diagnostics unavailable
  fprintf(
      trial_debug_file,
      "# accepted_step trial_number retry_number_before retry_number_after "
      "status_name status_value t_start h_trial h_error_controller h_proposed "
      "err_norm limiting_component limiting_component_name "
      "limiting_delta_5_minus_4 limiting_error_absolute limiting_scale "
      "limiting_error_normalized trial_result x_start y_start z_start r_start\n");
"""
        trial_debug_trial_metadata = r"""
    const long int accepted_step_before_trial = accepted_steps;
    const long int trial_number = rkf45_attempts + 1;
    const int retry_number_before = *rejection_retries;
    const double t_start = coordinate_time;
    const double h_trial = *h;
    const double x_start = f_start[1];
    const double y_start = f_start[2];
    const double z_start = f_start[3];
    const double r_start = sqrt(
        x_start * x_start + y_start * y_start + z_start * z_start);
"""
        trial_debug_call_argument = "&trial_debug,\n        "
        trial_debug_record = r"""
    const int trial_status_value = (int)*status;
    const char *trial_status_name =
        (trial_status_value >= 0 && trial_status_value < 11)
            ? status_names[trial_status_value]
            : "UNKNOWN_STATUS";
    const int trial_component_value = trial_debug.limiting_component;
    const char *trial_component_name =
        (trial_component_value >= 0 && trial_component_value < 9)
            ? trial_component_names[trial_component_value]
            : "UNKNOWN_COMPONENT";
    const int retry_number_after = *rejection_retries;
    const char *trial_result = "FAILED_OTHER";
    if (*status == ACTIVE) {
      trial_result = "ACCEPTED";
    } else if (*status == REJECTED) {
      trial_result = "REJECTED";
    } else if (*status == FAILURE_RKF45_REJECTION_LIMIT) {
      trial_result = "FAILED_REJECTION_LIMIT";
    } else if (*status == FAILURE_SPATIAL_INTERPOLATION) {
      trial_result = "FAILED_SPATIAL_INTERPOLATION";
    } else if (*status == FAILURE_TEMPORAL_INTERPOLATION) {
      trial_result = "FAILED_TEMPORAL_INTERPOLATION";
    } // END ELSE IF: classify RKF45 trial result
    fprintf(
        trial_debug_file,
        "%ld %ld %d %d %s %d %.17g %.17g %.17g %.17g %.17g %d %s "
        "%.17g %.17g %.17g %.17g %s %.17g %.17g %.17g %.17g\n",
        accepted_step_before_trial,
        trial_number,
        retry_number_before,
        retry_number_after,
        trial_status_name,
        trial_status_value,
        t_start,
        h_trial,
        trial_debug.h_error_controller,
        *h,
        trial_debug.err_norm,
        trial_component_value,
        trial_component_name,
        trial_debug.limiting_delta_5_minus_4,
        trial_debug.limiting_error_absolute,
        trial_debug.limiting_scale,
        trial_debug.limiting_error_normalized,
        trial_result,
        x_start,
        y_start,
        z_start,
        r_start);
    fflush(trial_debug_file);
"""
        trial_debug_cleanup = r"""
  if (trial_debug_file != NULL)
    fclose(trial_debug_file);
"""
        stage_debug_cleanup = r"""
  if (stage_debug_file != NULL)
    fclose(stage_debug_file);
"""
    else:
        trial_debug_declarations = ""
        stage_debug_header = ""
        stage_debug_declarations = ""
        stage_debug_open = ""
        stage_debug_record = ""
        stage_debug_cleanup = ""
        trial_debug_open = ""
        trial_debug_trial_metadata = ""
        trial_debug_call_argument = ""
        trial_debug_record = ""
        trial_debug_cleanup = ""

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdbool.h>",
        "<stdio.h>",
        "<stdlib.h>",
        "<string.h>",
    ]

    desc = rf"""Integrate one photon through a numerical spacetime on the CPU.

The executable initializes one photon from direct position and momentum parameters,
solves the initial null constraint, maps numerical time windows by coordinate-time
slot, and advances the state with the shared six-stage RKF45 pipeline. A trajectory
row, including signed normalization deviation, is written only after an accepted
step. ``initial_state.txt`` stores the initial event and direct contravariant
momentum at full precision before any normalized-state conversion. If the
numerical spacetime is axisymmetric about the Cartesian ``z`` axis
and RKF45 trial debugging is enabled, the initial $L_z$ is printed and each
accepted-state $L_z$ is appended to that row. When RKF45 trial debugging is
enabled, ``rkf45_trials.txt`` records every trial and ``rkf45_stages.txt``
records all six interpolation/RHS stages of every trial in execution order.
When a terminal or nonterminal plane is enabled, accepted steps pass through the
shared event detector. It interpolates each crossing state and writes one row per
detected plane to ``plane_crossings.txt`` with plane type, coordinate time,
affine parameter, local plane coordinates, and all nine interpolated state
components. A terminal-plane hit stops the integration only when its local radius
is within the configured bounds.

For normalized equations, the state layout is
``(lambda, x, y, z, u, Pi_1, Pi_2, Pi_3, L_normal)``: ``f[0]`` is lambda and
the RKF45 integration parameter is coordinate time. For non-normalized
equations, the state layout is ``(t, x, y, z, p^t, p^x, p^y, p^z, L_normal)``:
the RKF45 integration parameter is lambda and ``f[0]`` is coordinate time.

@param argc Number of command-line arguments.
@param[in] argv Command-line argument array.
@return EXIT_SUCCESS on success, or EXIT_FAILURE after a setup or integration error.

@note The numerical dataset coordinate system is {dataset_coord_system}.
@note The selected analytic photon equation family is {spacetime_name}.
"""

    cfunc_type = "int"
    name = "single_integrator_numerical"
    params = "int argc, const char *argv[]"

    body = rf"""
  //==========================================
  // 1. COMMONDATA AND PARAMETER SETUP
  //==========================================
  commondata_struct commondata;
  commondata_struct_set_to_default(&commondata);
  cmdline_input_and_parfile_parser(&commondata, argc, argv);

  const long int num_rays = 1;
  const long int chunk_size = 1;
  const int stream_idx = 0;
  const long int max_accepted_steps = 200000;
  const long int max_rkf45_attempts = 2000000;

  const char *status_names[] = {{
    "STOP_CONDITION_COORD_RADIUS_EXCEEDED",
    "STOP_CONDITION_TERMINAL_PLANE",
    "STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED",
    "FAILURE_RKF45_REJECTION_LIMIT",
    "STOP_CONDITION_T_MAX_EXCEEDED",
    "FAILURE_SLOT_MANAGER_ERROR",
    "FAILURE_GENERIC",
    "ACTIVE",
    "REJECTED",
    "FAILURE_SPATIAL_INTERPOLATION",
    "FAILURE_TEMPORAL_INTERPOLATION",
    "FAILURE_PLANE_INTERPOLATION_HISTORY"
  }};

  int exit_status = EXIT_SUCCESS;
  FILE *trajectory_file = NULL;
  FILE *plane_crossings_file = NULL;
{trial_debug_declarations}
{stage_debug_declarations}
{synthetic_slice_usage_declarations}
  bool slot_manager_initialized = false;
  bool numerical_window_initialized = false;

  double *f = NULL;
  double *f_start = NULL;
  double *f_temp = NULL;
  double *metric = NULL;
  double observer_metric[10];
  double observer_tetrad[4][4];
  double *rhs_geometry = NULL;
  double *k_bundle = NULL;
  double *integration_param = NULL;
  double *h = NULL;
  int *rejection_retries = NULL;
  termination_type_t *status = NULL;
  PhotonStateSoA initial_photon = {{0}};
  double f_event_history[{maximum_degree} * 9] = {{0.0}};
  double integration_param_event_history[{maximum_degree}] = {{0.0}};
  bool on_positive_side_of_non_terminal_plane_prev = false;
  bool on_positive_side_of_terminal_plane_prev = false;
  bool non_terminal_plane_event_found = false;
  bool terminal_plane_event_found = false;
  bool non_terminal_plane_crossing_pending = false;
  bool terminal_plane_crossing_pending = false;
  int non_terminal_plane_steps_past = 0;
  int terminal_plane_steps_past = 0;
  int non_terminal_plane_event_degree = 0;
  int terminal_plane_event_degree = 0;
  double non_terminal_plane_event_state[9] = {{0.0}};
  double terminal_plane_event_state[9] = {{0.0}};
  blueprint_data_t plane_crossing_results = {{0}};
  const long int event_chunk_index[1] = {{0}};
  const bool event_planes_enabled =
      commondata.non_terminal_plane_enabled || commondata.terminal_plane_enabled;

  //==========================================
  // 2. SINGLE-RAY CPU MEMORY
  //==========================================
  BHAH_MALLOC(f, sizeof(double) * 9);
  BHAH_MALLOC(f_start, sizeof(double) * 9);
  BHAH_MALLOC(f_temp, sizeof(double) * 9);
  BHAH_MALLOC(metric, sizeof(double) * 10);
  BHAH_MALLOC(rhs_geometry, sizeof(double) * 40);
  BHAH_MALLOC(k_bundle, sizeof(double) * 6 * 9);
  BHAH_MALLOC(integration_param, sizeof(double));
  BHAH_MALLOC(h, sizeof(double));
  BHAH_MALLOC(rejection_retries, sizeof(int));
  BHAH_MALLOC(status, sizeof(termination_type_t));

  if (f == NULL || f_start == NULL || f_temp == NULL || metric == NULL ||
      rhs_geometry == NULL || k_bundle == NULL || integration_param == NULL ||
      h == NULL || rejection_retries == NULL || status == NULL) {{
    fprintf(stderr, "ERROR: failed to allocate single-photon CPU buffers.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: single-photon CPU allocation failed

  //==========================================
  // 3. TIME-SLOT AND NUMERICAL-WINDOW SETUP
  //==========================================
  TimeSlotManager tsm;
  const double slot_manager_t_max = commondata.t_start + 1.0e-5;
  slot_manager_init(
      &tsm,
      commondata.slot_manager_t_min,
      slot_manager_t_max,
      commondata.slot_manager_delta_t,
      num_rays);
  slot_manager_initialized = true;

  if (tsm.photon_next_ptrs == NULL || tsm.slot_heads == NULL ||
      tsm.slot_counts == NULL) {{
    fprintf(stderr, "ERROR: failed to allocate TimeSlotManager storage.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: TimeSlotManager allocation failed

  NumericalTimeWindowManager numerical_window;
  time_window_manager_numerical_set_inert(&numerical_window);

  commondata_struct commondata_for_params_defaults = commondata;
  griddata_struct dummy_griddata[MAXNUMGRIDS];
  params_struct_set_to_default(&commondata_for_params_defaults, dummy_griddata);
  params_struct numerical_params = dummy_griddata[0].params;

  if (commondata.numerical_spacetime_bin_path[0] == '\0') {{
    fprintf(
        stderr,
        "ERROR: numerical_spacetime_bin_path is empty. Set it in the generated "
        "parameter file.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: numerical_spacetime_bin_path was empty

  if (time_window_manager_numerical_init(
          &numerical_window,
          commondata.numerical_spacetime_bin_path,
          &commondata,
          commondata.numerical_spacetime_temporal_interp_half_width,
          &numerical_params) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
    fprintf(
        stderr,
        "ERROR: failed to initialize the numerical time-window manager from '%s'.\n",
        commondata.numerical_spacetime_bin_path);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: numerical time-window manager initialization failed
  numerical_window_initialized = true;
  printf("Numerical spacetime first stored slice time: %.15e\n",
         (double)commondata.t_numerical_initial);

  const REAL observer_position[3] = {{
      (REAL)commondata.observer_x,
      (REAL)commondata.observer_y,
      (REAL)commondata.observer_z}};
  if (time_window_manager_numerical_validate_startup_domain(
          &numerical_window,
          &commondata,
          &numerical_params,
          observer_position,
          "observer_x, observer_y, and observer_z") !=
      TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
    fprintf(stderr,
            "ERROR: numerical observer/r_escape startup stencil validation failed.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: startup domain validation failed

  if (numerical_params.Nxx{phi_dim} != 2) {{
    fprintf(
        stderr,
        "ERROR: numerical spatial interpolation expects exactly two stored phi "
        "planes in native dimension {phi_dim}; got Nxx{phi_dim}=%d.\n",
        numerical_params.Nxx{phi_dim});
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: dataset lacked two phi planes

  azimuthal_symmetry_spatial_lagrange_context_struct spatial_context;
  spatial_context.stored_phi_samples[0] =
      numerical_params.xxmin{phi_dim} + 0.5 * numerical_params.dxx{phi_dim};
  spatial_context.stored_phi_samples[1] =
      numerical_params.xxmin{phi_dim} + 1.5 * numerical_params.dxx{phi_dim};

  //==========================================
  // 4. OBSERVER-TETRAD INITIAL CONDITIONS
  //==========================================
  // Reuse the batch initializer's observer parameters. The single-ray
  // command-line momentum is the observer look-forward direction seed, while q
  // sets the initial observer-frame energy to one.
  // Single-ray execution is the one-sample-per-tile specialization of the
  // shared angular sample mapping. Runtime tile counts and indices permit
  // exact reproduction of any batch-camera pixel; their defaults select the
  // center ray of a one-tile camera.
  commondata.scan_density = 1;
  // The temporary state supplies the observer event to one metric
  // interpolation. The shared initializer overwrites f and h with the
  // validated complete direct tetrad ray.
  f[0] = {initial_state_time};
  f[1] = commondata.observer_x;
  f[2] = commondata.observer_y;
  f[3] = commondata.observer_z;
  for (int component = 4; component < 9; ++component) {{
    f[component] = 0.0;
  }} // END LOOP: for component over temporary momentum

  *integration_param = {initial_integration_value};
  *h = commondata.initial_h;

  for (int component = 0; component < 9; ++component) {{
    if (!isfinite(f[component])) {{
      fprintf(stderr, "ERROR: initial photon state was not finite.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: one initial photon state component
  }} // END LOOP: for component over initial photon
  // The temporary observer state intentionally has zero momentum. The shared
  // tetrad initializer supplies and validates the complete momentum after the
  // one observer-metric interpolation below.
  if (!isfinite(*h) || *h == 0.0) {{
    fprintf(stderr, "ERROR: initial_h must be finite and nonzero.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial integration step was invalid

  int mapped_slot_index = -1;
  const int initial_slot_index = slot_get_index(&tsm, commondata.t_start);
  if (initial_slot_index < 0) {{
    fprintf(
        stderr,
        "ERROR: t_start=%e is outside the configured TimeSlotManager range.\n",
        (double)commondata.t_start);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial coordinate time outside bounds

  if (time_window_manager_numerical_mmap_for_slot(
          &numerical_window, &tsm, initial_slot_index) !=
      TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
    fprintf(
        stderr,
        "ERROR: failed to map the initial numerical time window for slot %d.\n",
        initial_slot_index);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial numerical time-window mapping failed
  mapped_slot_index = initial_slot_index;

  // Observer interpolation uses a temporary status because failure here means
  // no photon can be initialized; it is therefore a project-level error.
  termination_type_t observer_interpolation_status = ACTIVE;
  numerical_interpolation(
      &commondata,
      &numerical_params,
      &spatial_context,
      &numerical_window,
      f,
      &observer_interpolation_status,
      {interpolation_initial_arguments}
      {synthetic_slice_usage_null_arguments}
      metric,
      NULL,
      chunk_size,
      stream_idx);

  if (observer_interpolation_status == FAILURE_SPATIAL_INTERPOLATION ||
      observer_interpolation_status == FAILURE_TEMPORAL_INTERPOLATION) {{
    fprintf(
        stderr,
        "ERROR: observer metric interpolation failed with status %d.\n",
        (int)observer_interpolation_status);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: observer metric interpolation failed

  // Preserve the one interpolated observer metric while the shared initializer
  // constructs and validates the observer tetrad.
  for (int component = 0; component < 10; ++component) {{
    observer_metric[component] = metric[component];
  }} // END LOOP: for component over observer metric

  initial_photon.f = f;
  initial_photon.h = h;
  initial_photon.integration_param = integration_param;
  initial_photon.f_history = f_event_history;
  initial_photon.integration_param_history = integration_param_event_history;
  initial_photon.on_positive_side_of_non_terminal_plane_prev =
      &on_positive_side_of_non_terminal_plane_prev;
  initial_photon.on_positive_side_of_terminal_plane_prev =
      &on_positive_side_of_terminal_plane_prev;
  initial_photon.non_terminal_plane_crossing_pending =
      &non_terminal_plane_crossing_pending;
  initial_photon.terminal_plane_crossing_pending =
      &terminal_plane_crossing_pending;
  initial_photon.non_terminal_plane_steps_past = &non_terminal_plane_steps_past;
  initial_photon.terminal_plane_steps_past = &terminal_plane_steps_past;
  initial_photon.non_terminal_plane_event_degree =
      &non_terminal_plane_event_degree;
  initial_photon.terminal_plane_event_degree = &terminal_plane_event_degree;
  set_initial_conditions_kernel(
      &commondata, num_rays, &initial_photon, observer_metric, observer_tetrad);
  // Retain the direct tangent before normalized-EOM conversion. It defines
  // the inertial straight-line reference for an analytic spacetime.
  const double initial_direct_momentum[4] = {{f[4], f[5], f[6], f[7]}};
  {momentum_conversion_call}

  // Seed accepted-state history after normalized-momentum conversion.
  for (int history_step = 0; history_step < {maximum_degree}; ++history_step) {{
    for (int component = 0; component < 9; ++component)
      f_event_history[history_step * 9 + component] = f[component];
    integration_param_event_history[history_step] = *integration_param;
  }} // END LOOP: initialize accepted-state history

  normalization_constraint_t initial_normalization;
  {normalization_kernel_name}(
      f, metric, &initial_normalization, chunk_size, stream_idx);
  const double initial_constraint_error =
      fabs(initial_normalization.C - {"1.0" if normalized_eom else "0.0"});
  if (!isfinite(initial_constraint_error) ||
      initial_constraint_error > 1.0e-9) {{
    fprintf(
        stderr,
        "ERROR: initial tetrad-ray constraint residual=%e exceeds tolerance.\n",
        initial_constraint_error);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial tetrad-ray constraint invalid

  for (int component = 0; component < 9; ++component) {{
    if (!isfinite(f[component])) {{
      fprintf(
          stderr,
          "ERROR: initial constrained state component %d was not finite.\n",
          component);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: one initial constrained state component
  }} // END LOOP: for component over initial state
{initial_Lz_evaluation}

  *rejection_retries = 0;
  *status = ACTIVE;

  printf("Initial State:\n");
  printf("  Pos (%.4f, %.4f, %.4f)\n", f[1], f[2], f[3]);
  printf(
      "  Energy/momentum state (f[4], f[5], f[6], f[7]) = "
      "(%.4f, %.4f, %.4f, %.4f)\n",
      f[4],
      f[5],
      f[6],
      f[7]);
  printf("  Integration parameter = %.15e\n", {trajectory_lambda_expression});
  printf("  Coordinate time = %.15e\n", {trajectory_time_expression});

  trajectory_file = fopen("trajectory.txt", "w");
  if (trajectory_file == NULL) {{
    fprintf(stderr, "ERROR: could not open trajectory.txt for writing.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: trajectory output unavailable
  fprintf(
      trajectory_file,
      "{trajectory_header}");
  FILE *initial_state_file = fopen("initial_state.txt", "w");
  if (initial_state_file == NULL) {{
    fprintf(stderr, "ERROR: could not open initial_state.txt for writing.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial-state output unavailable
  fprintf(initial_state_file, "# t x y z p^t p^x p^y p^z\n");
  fprintf(
      initial_state_file,
      "%.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
      (double)commondata.t_start, f[1], f[2], f[3],
      initial_direct_momentum[0], initial_direct_momentum[1],
      initial_direct_momentum[2], initial_direct_momentum[3]);
  if (fclose(initial_state_file) != 0) {{
    fprintf(stderr, "ERROR: failed to close initial_state.txt.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial-state output failed
{trial_debug_open}
{stage_debug_open}

  printf("Starting CPU numerical single-photon integration.\n");
  printf("  spacetime equations: {spacetime_name}\n");
  printf("  dataset coordinates: {dataset_coord_system}\n");
  printf("  normalized_eom: %s\n", {str(normalized_eom).lower()} ? "true" : "false");
  printf("  numerical data: %s\n", commondata.numerical_spacetime_bin_path);

  //==========================================
  // 5. SINGLE-PHOTON RKF45 LOOP
  //==========================================
  long int accepted_steps = 0;
  long int rkf45_attempts = 0;

  while (accepted_steps < max_accepted_steps &&
         rkf45_attempts < max_rkf45_attempts) {{
    const double coordinate_time = {coordinate_time_expression};
    const int slot_index = slot_get_index(&tsm, coordinate_time);
    if (slot_index < 0) {{
      *status = STOP_CONDITION_T_MAX_EXCEEDED;
      printf(
          "Coordinate time %.15e left the configured numerical slot range.\n",
          coordinate_time);
      break;
    }} // END IF: current time left data window

    if (slot_index != mapped_slot_index) {{
      if (time_window_manager_numerical_mmap_for_slot(
              &numerical_window, &tsm, slot_index) !=
          TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
        fprintf(
            stderr,
            "ERROR: failed to map numerical time window for slot %d at t=%e.\n",
            slot_index,
            coordinate_time);
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: numerical time-window mapping failed
      mapped_slot_index = slot_index;
    }} // END IF: photon moved to different slot

    memcpy(f_start, f, sizeof(double) * 9);
    memcpy(f_temp, f, sizeof(double) * 9);
{trial_debug_trial_metadata}
{trial_spatial_center_setup}

    for (int stage = 1; stage <= 6; ++stage) {{
      numerical_interpolation(
          &commondata,
          &numerical_params,
          &spatial_context,
          &numerical_window,
          f_temp,
          status,
          {interpolation_stage_arguments}
          {synthetic_slice_usage_stage_arguments}
          metric,
          rhs_geometry,
          chunk_size,
          stream_idx);

      calculate_ode_rhs_kernel(
          f_temp,
          metric,
          rhs_geometry,
          {rhs_integration_arguments}
          k_bundle,
          stage,
          chunk_size,
          stream_idx);
{stage_debug_record}
      if (stage < 6) {{
        rkf45_stage_update(
            f_start,
            k_bundle,
            h,
            stage,
            chunk_size,
            f_temp,
            stream_idx);
      }} // END IF: skip stage-6 intermediate state update
    }} // END LOOP: for stage over RKF45 stages

    rkf45_finalize_and_control(
        &commondata,
        f,
        f_start,
        k_bundle,
        h,
        status,
        integration_param,
        rejection_retries,
        {trial_debug_call_argument}chunk_size,
        stream_idx);
    rkf45_attempts++;
{trial_debug_record}

    if (*status == ACTIVE) {{
      for (int component = 0; component < 9; ++component) {{
        if (!isfinite(f[component])) {{
          fprintf(
              stderr,
              "ERROR: accepted state component %d was not finite after attempt %ld.\n",
              component,
              rkf45_attempts);
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: accepted state component invalid
      }} // END LOOP: for component over accepted state

{accepted_Lz_before_interpolation}
      numerical_interpolation(
          &commondata,
          &numerical_params,
          &spatial_context,
          &numerical_window,
          f,
          status,
          {accepted_metric_interpolation_arguments}
          {synthetic_slice_usage_null_arguments}
          metric,
          NULL,
          chunk_size,
          stream_idx);

      if (*status == FAILURE_SPATIAL_INTERPOLATION ||
          *status == FAILURE_TEMPORAL_INTERPOLATION) {{
        // RKF45 already accepted and committed this state. Preserve it in the
        // trajectory even though its normalization diagnostic is unavailable.
        fprintf(
            trajectory_file,
            "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e{failed_interpolation_Lz_format}\n",
            {trajectory_lambda_expression},
            {trajectory_time_expression},
            f[1],
            f[2],
            f[3],
            f[4],
            f[5],
            f[6],
            f[7],
            f[8],
            NAN{failed_interpolation_Lz_argument});
        fflush(trajectory_file);
        accepted_steps++;
        printf(
            "Accepted-state interpolation failed with status %d; "
            "preserving the accepted state.\n",
            (int)*status);
        break;
      }} // END IF: accepted-state interpolation failed

{accepted_Lz_after_interpolation}
      {log_energy_evaluation}
      if (!isfinite(log_energy_measure)) {{
        fprintf(stderr, "ERROR: accepted-state log-energy measure was not finite.\n");
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: accepted-state log-energy measure invalid

      if (event_planes_enabled) {{
        event_detection_manager_kernel(
            &commondata,
            f,
            &log_energy_measure,
            f_event_history,
            integration_param,
            integration_param_event_history,
            &plane_crossing_results,
            status,
            &on_positive_side_of_non_terminal_plane_prev,
            &on_positive_side_of_terminal_plane_prev,
            &non_terminal_plane_event_found,
            &terminal_plane_event_found,
            &non_terminal_plane_crossing_pending,
            &terminal_plane_crossing_pending,
            &non_terminal_plane_steps_past,
            &terminal_plane_steps_past,
            &non_terminal_plane_event_degree,
            &terminal_plane_event_degree,
            non_terminal_plane_event_state,
            terminal_plane_event_state,
            accepted_steps + 1 >= max_accepted_steps ||
                rkf45_attempts >= max_rkf45_attempts ||
                slot_get_index(&tsm, {coordinate_time_expression}) < 0,
            event_chunk_index,
            chunk_size,
            stream_idx);
      }} // END IF: a physical plane is enabled

      normalization_constraint_t normalization;
      {normalization_kernel_name}(
          f, metric, &normalization, chunk_size, stream_idx);
      const double trajectory_norm = {normalization_diagnostic_expression};
      if (!isfinite(trajectory_norm)) {{
        fprintf(stderr, "ERROR: accepted-state norm was not finite.\n");
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: accepted-state norm invalid

      fprintf(
          trajectory_file,
          "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e{accepted_Lz_format}\n",
          {trajectory_lambda_expression},
          {trajectory_time_expression},
          f[1],
          f[2],
          f[3],
          f[4],
          f[5],
          f[6],
          f[7],
          f[8],
          trajectory_norm{accepted_Lz_argument});
      fflush(trajectory_file);
      accepted_steps++;

      if (event_planes_enabled && *status != ACTIVE) {{
        if (*status == STOP_CONDITION_TERMINAL_PLANE) {{
          printf("Photon crossed the configured terminal plane.\n");
        }} else if (*status == STOP_CONDITION_COORD_RADIUS_EXCEEDED) {{
          printf("Photon escaped to r > %.15e.\n", (double)commondata.r_escape);
        }} else if (*status == STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED) {{
          printf(
              "Evolution measure exceeded %.15e.\n",
              commondata.evolution_measure_max);
        }} else if (*status == FAILURE_PLANE_INTERPOLATION_HISTORY ||
                   *status == FAILURE_GENERIC) {{
          fprintf(stderr, "Plane crossing failed with status %d.\n", (int)*status);
          exit_status = EXIT_FAILURE;
        }} // END ELSE IF: event manager established a stop status
        break;
      }} // END IF: event manager stopped the photon

      const double radius_squared = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
      if (radius_squared > commondata.r_escape * commondata.r_escape) {{
        *status = STOP_CONDITION_COORD_RADIUS_EXCEEDED;
        printf("Photon escaped to r > %.15e.\n", (double)commondata.r_escape);
        break;
      }} // END IF: photon state crossed boundary

      if (log_energy_measure > commondata.evolution_measure_max) {{
        *status = STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED;
        printf("Evolution measure exceeded %.15e.\n", commondata.evolution_measure_max);
        break;
      }} // END IF: evolution measure exceeded limit
    }} else if (*status == REJECTED) {{
      *status = ACTIVE;
      continue;
    }} // END ELSE IF: retry rejected RKF45 step
    else if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
      printf("RKF45 reached its consecutive-rejection limit.\n");
      break;
    }} else if (*status == FAILURE_SPATIAL_INTERPOLATION ||
               *status == FAILURE_TEMPORAL_INTERPOLATION) {{
      printf("Interpolation failed with status %d.\n", (int)*status);
      break;
    }} else {{
      fprintf(stderr, "ERROR: unexpected integration status %d.\n", (int)*status);
      exit_status = EXIT_FAILURE;
      goto cleanup;
  }} // END ELSE: unexpected RKF45 finalization status
  }} // END WHILE: evolve photon through accepted steps

  if (event_planes_enabled) {{
    plane_crossings_file = fopen("plane_crossings.txt", "w");
    if (plane_crossings_file == NULL) {{
      fprintf(stderr, "ERROR: could not open plane_crossings.txt for writing.\n");
      exit_status = EXIT_FAILURE;
    }} else {{
      fprintf(
          plane_crossings_file,
          "# plane_type coordinate_time affine_parameter local_y local_z {event_state_columns} interpolation_degree\n");
      if (non_terminal_plane_event_found) {{
        fprintf(
            plane_crossings_file,
            "nonterminal %.17e %.17e %.17e %.17e {event_state_format} %d\n",
            plane_crossing_results.non_terminal_plane_t,
            plane_crossing_results.non_terminal_plane_lambda,
            plane_crossing_results.y_nt,
            plane_crossing_results.z_nt,
            {non_terminal_event_state_arguments},
            non_terminal_plane_event_degree);
        printf(
            "Nonterminal-plane crossing: t=%.15e, lambda=%.15e, "
            "local=(%.15e, %.15e)\n",
            plane_crossing_results.non_terminal_plane_t,
            plane_crossing_results.non_terminal_plane_lambda,
            plane_crossing_results.y_nt,
            plane_crossing_results.z_nt);
      }} // END IF: nonterminal plane was crossed
      if (terminal_plane_event_found) {{
        fprintf(
            plane_crossings_file,
            "terminal %.17e %.17e %.17e %.17e {event_state_format} %d\n",
            plane_crossing_results.t_f,
            plane_crossing_results.L_f,
            plane_crossing_results.y_t,
            plane_crossing_results.z_t,
            {terminal_event_state_arguments},
            terminal_plane_event_degree);
        printf(
            "Terminal-plane crossing: t=%.15e, lambda=%.15e, "
            "local=(%.15e, %.15e)\n",
            plane_crossing_results.t_f,
            plane_crossing_results.L_f,
            plane_crossing_results.y_t,
            plane_crossing_results.z_t);
      }} // END IF: terminal plane was crossed
      if (fclose(plane_crossings_file) != 0) {{
        fprintf(stderr, "ERROR: failed to close plane_crossings.txt.\n");
        exit_status = EXIT_FAILURE;
      }}
      plane_crossings_file = NULL;
    }} // END ELSE: crossing output file opened
  }} // END IF: at least one physical plane is enabled

  if ((*status == ACTIVE || *status == REJECTED) &&
      accepted_steps >= max_accepted_steps) {{
    *status = FAILURE_GENERIC;
    printf("Integration stopped at the accepted-step safety limit.\n");
  }} else if ((*status == ACTIVE || *status == REJECTED) &&
             rkf45_attempts >= max_rkf45_attempts) {{
    *status = FAILURE_GENERIC;
    printf("Integration stopped at the RKF45-attempt safety limit.\n");
  }} // END ELSE IF: integration reached the RKF45-attempt safety

  const int final_status_index = (int)*status;
  const char *final_status_name =
      (final_status_index >= 0 && final_status_index < 12)
          ? status_names[final_status_index]
          : "UNKNOWN_STATUS";
  printf(
      "Integration finished after %ld accepted steps and %ld RKF45 attempts.\n",
      accepted_steps,
      rkf45_attempts);
  printf(
      "Final status: %s (%d), lambda=%.15e, t=%.15e\n",
      final_status_name,
      final_status_index,
      {trajectory_lambda_expression},
      {trajectory_time_expression});

  //==========================================
  // 6. OPTIONAL TERMINAL NORMALIZATION CHECK
  //==========================================
  if (commondata.perform_normalization_check) {{
    const double terminal_coordinate_time = {coordinate_time_expression};
    const int terminal_slot_index = slot_get_index(&tsm, terminal_coordinate_time);
    if (terminal_slot_index < 0) {{
      printf(
          "Terminal normalization skipped: t=%.15e is outside the numerical "
          "slot range.\n",
          terminal_coordinate_time);
    }} else {{
      bool terminal_normalization_ready = true;
      if (terminal_slot_index != mapped_slot_index) {{
        if (time_window_manager_numerical_mmap_for_slot(
                &numerical_window, &tsm, terminal_slot_index) !=
            TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
          printf(
              "Terminal normalization skipped: failed to map its numerical "
              "time window.\n");
          terminal_normalization_ready = false;
        }} else {{
          mapped_slot_index = terminal_slot_index;
        }} // END ELSE: terminal normalization time window mapped
      }} // END IF: terminal state changed slot

      if (terminal_normalization_ready) {{
        // A separate status prevents this optional diagnostic from replacing
        // the physical termination status established by the evolution.
        termination_type_t terminal_interpolation_status = ACTIVE;
        numerical_interpolation(
            &commondata,
            &numerical_params,
            &spatial_context,
            &numerical_window,
            f,
            &terminal_interpolation_status,
            {interpolation_initial_arguments}
            {synthetic_slice_usage_null_arguments}
            metric,
            NULL,
            chunk_size,
            stream_idx);

        if (terminal_interpolation_status == FAILURE_SPATIAL_INTERPOLATION ||
            terminal_interpolation_status == FAILURE_TEMPORAL_INTERPOLATION) {{
          printf(
              "Terminal normalization skipped: interpolation failed with "
              "status %d.\n",
              (int)terminal_interpolation_status);
        }} else {{
          normalization_constraint_t normalization;
          {normalization_kernel_name}(
              f, metric, &normalization, chunk_size, stream_idx);
          const double normalization_deviation =
              {normalization_diagnostic_expression};
          if (!isfinite(normalization_deviation)) {{
            printf(
                "Terminal normalization skipped: deviation was not finite.\n");
          }} else {{
            printf(
                "Final signed normalization deviation: %.15e\n",
                normalization_deviation);
          }} // END ELSE: terminal normalization deviation is finite
        }} // END ELSE: terminal normalization interpolation succeeded
      }} // END IF: terminal normalization window available
    }} // END ELSE: terminal state inside window
  }} // END IF: terminal normalization diagnostics were requested

  cleanup:
{synthetic_slice_usage_report}
{trial_debug_cleanup}
{stage_debug_cleanup}
  if (trajectory_file != NULL)
    fclose(trajectory_file);
  if (plane_crossings_file != NULL)
    fclose(plane_crossings_file);
  if (numerical_window_initialized)
    time_window_manager_numerical_free(&numerical_window);
  if (slot_manager_initialized)
    slot_manager_free(&tsm);

  BHAH_FREE(f);
  BHAH_FREE(f_start);
  BHAH_FREE(f_temp);
  BHAH_FREE(metric);
  BHAH_FREE(rhs_geometry);
  BHAH_FREE(k_bundle);
  BHAH_FREE(integration_param);
  BHAH_FREE(h);
  BHAH_FREE(rejection_retries);
  BHAH_FREE(status);

  return exit_status;
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    # Step 1: Select the CPU/OpenMP BHaH backend.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    # Step 2: Select the project, numerical dataset, and equation convention.
    PROJECT_NAME = "single_integrator_numerical"
    project_dir = os.path.abspath(os.path.join("project", PROJECT_NAME))
    parfile_path = os.path.join(project_dir, f"{PROJECT_NAME}.par")

    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "photon"
    DATASET_COORD_SYSTEM = "SinhCylindricalv2n2"
    INTERPOLATION_METHOD = "g4DD"
    NORMALIZED_EOM = False
    rhs_uses_metric_derivatives = INTERPOLATION_METHOD != "GammaUDD"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"

    # Step 3: Recreate the generated project directory.
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Step 4: Acquire the photon equations used by the shared RKF45 RHS kernel.
    print(f"Acquiring symbolic data for {GEO_KEY}...")
    geodesic_data = Geodesic_Equations[GEO_KEY]

    # Step 5: Register the numerical interpolation and shared RKF45 pipeline.
    print("Registering numerical single-photon kernels...")
    if rhs_uses_metric_derivatives:
        geodesic_rhs = (
            geodesic_data.geodesic_eom_rhs_photon_normalized()
            if NORMALIZED_EOM
            else geodesic_data.geodesic_eom_rhs_photon()
        )
    else:
        geodesic_rhs = (
            geodesic_data.geodesic_eom_rhs_photon_normalized_christoffel()
            if NORMALIZED_EOM
            else geodesic_data.geodesic_eom_rhs_photon_christoffel()
        )

    u_expr, PiD_exprs = geodesic_data.photon_momentum_to_normalized_quantities()
    if NORMALIZED_EOM:
        photon_momentum_to_normalized_kernel.photon_momentum_to_normalized_kernel(
            u_expr, PiD_exprs
        )
    else:
        normal_observer_log_energy.normal_observer_log_energy(u_expr)

    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )
    numerical_interpolation.register_CFunction_numerical_interpolation(
        DATASET_COORD_SYSTEM,
        interpolation_method=INTERPOLATION_METHOD,
        enable_simd=False,
        project_dir=project_dir,
        normalized_eom=NORMALIZED_EOM,
        skip_non_active_status=True,
    )
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_rhs,
        geodesic_data.xx,
        rhs_uses_metric_derivatives=rhs_uses_metric_derivatives,
        normalized_eom=NORMALIZED_EOM,
    )
    rkf45_stage_update.rkf45_stage_update()
    find_event_time_and_state.find_event_time_and_state(4)
    handle_terminal_plane_intersection.handle_terminal_plane_intersection()
    handle_non_terminal_plane_intersection.handle_non_terminal_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel(
        maximum_degree=4,
        normalized_eom=NORMALIZED_EOM,
        capture_event_state=True,
        single_photon_step_limit=True,
    )
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True,
        retry_interpolation_failures=True,
    )
    single_integrator_numerical(
        SPACETIME,
        DATASET_COORD_SYSTEM,
        interpolation_method=INTERPOLATION_METHOD,
        normalized_eom=NORMALIZED_EOM,
    )
    main_single.main_single("single_integrator_numerical")

    for internal_func in [
        "find_event_time_and_state_centered",
        "handle_terminal_plane_intersection",
        "handle_non_terminal_plane_intersection",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 7: Generate parameter headers, the default parfile, and CPU definitions.
    print("Generating headers, parameter handling, and Makefile...")
    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=PROJECT_NAME
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=PROJECT_NAME,
        cmdline_inputs=[
            "tiles_width",
            "tiles_height",
            "tile_index_width",
            "tile_index_height",
        ],
    )

    # Shared RKF45 kernels retain architecture-neutral intrinsic names. These
    # scalar definitions execute entirely on the CPU and require no GPU runtime.
    cpu_scalar_intrinsics = {
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))\n",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))\n",
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))\n",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))\n",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))\n",
        "FusedMulAddCUDA(a, b, c)": (
            "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))\n"
        ),
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)\n",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)\n",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)\n",
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC\n",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline\n",
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": (
            "#define BHAH_WARP_ATOMIC_ADD(ptr, val) "
            '_Pragma("omp atomic") *(ptr) += (val)\n'
        ),
        "GLOBAL_COMMONDATA_EXTERN": (
            "// CPU execution passes commondata explicitly.\n"
        ),
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while (0)\n",
    }
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        supplemental_defines_dict=cpu_scalar_intrinsics,
    )

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=PROJECT_NAME,
        exec_or_library_name=PROJECT_NAME,
        compiler_opt_option="fast",
        addl_CFLAGS=[
            "-fopenmp",
            "-O3",
            "-DDEBUG",
            "-Wno-stringop-truncation",
        ],
        addl_libraries=["-lm"],
        CC="gcc",
        src_code_file_ext="c",
    )

    # Step 8: Copy the packaged trajectory visualizer.
    copy_files(
        package="nrpy.examples.geodesic_visualizations",
        filenames_list=["visualize_trajectory.py"],
        project_dir=project_dir,
        subdirectory="",
    )

    print(f"Finished generating {project_dir}.")
    print(f"Set numerical_spacetime_bin_path in {parfile_path} before running.")
    print(f"Build with: cd {project_dir} && make")
    print(f"Run with: ./{PROJECT_NAME} {PROJECT_NAME}.par")
    print("Accepted trajectory samples will be written to trajectory.txt.")

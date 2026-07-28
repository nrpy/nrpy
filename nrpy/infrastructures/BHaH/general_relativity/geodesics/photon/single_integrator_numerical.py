"""
Generate a standalone CPU project for one photon in a numerical spacetime.

The generated C executable follows the numerical batch integrator's numerical
interpolation, time-window management, initial-condition geometry, and RKF45
control logic, but removes batching, double buffering, device memory, streams,
and blueprint output. It writes one trajectory row after every accepted RKF45
step, including a signed normalization diagnostic. Optional RKF45 debugging
writes trial-level and stage-level records from this single-photon host loop.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

# pylint: disable=too-many-lines

import os
import shutil

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.equations.general_relativity.geodesics.geodesics import Geodesic_Equations
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation import (
    numerical_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    calculate_ode_rhs_kernel,
    main_single,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
)


def single_integrator_numerical(  # pylint: disable=invalid-name,too-many-locals
    spacetime_name: str,
    dataset_coord_system: str,
    interpolation_method: str = "g4DD",
    normalized_eom: bool = False,
    enable_rkf45_trial_debug: bool = False,
) -> None:
    """
    Register the standalone numerical-spacetime single-photon C integrator.

    :param spacetime_name: Spacetime identifier used to select photon equations.
    :param dataset_coord_system: Coordinate system used by the numerical dataset.
    :param interpolation_method: Numerical geometry payload method used by the generated project.
    :param normalized_eom: Whether to evolve normalized coordinate-time equations.
    :param enable_rkf45_trial_debug: Whether to write one diagnostic row for every
        RKF45 trial to ``rkf45_trials.txt`` and one row for each of its six
        stages to ``rkf45_stages.txt``.
    :raises ValueError: If the interpolation method or dataset coordinate system
        is unsupported.

    Doctests:
    >>> import os
    >>> import nrpy.c_function as cfc
    >>> os.environ["XDG_CACHE_HOME"] = "/tmp"
    >>> cfc.CFunction_dict.clear()
    >>> single_integrator_numerical("Numerical", "SinhCylindricalv2n2")
    >>> generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "# lambda t x y z energy_measure p_x p_y p_z aux norm" in generated
    True
    >>> "const double trajectory_norm = normalization.C;" in generated
    True
    >>> "fabs(normalization.C)" not in generated
    True
    >>> cfc.CFunction_dict.clear()
    >>> single_integrator_numerical("Numerical", "SinhCylindricalv2n2", normalized_eom=True)
    >>> generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "const double trajectory_norm = normalization.C - 1.0;" in generated
    True
    >>> "fabs(normalization.C - 1.0)" not in generated
    True
    >>> "trial_debug_file" not in generated
    True
    >>> "stage_debug_file" not in generated
    True
    >>> cfc.CFunction_dict.clear()
    >>> single_integrator_numerical(
    ...     "Numerical", "SinhCylindricalv2n2", enable_rkf45_trial_debug=True
    ... )
    >>> generated = cfc.CFunction_dict["single_integrator_numerical"].full_function
    >>> "rkf45_trials.txt" in generated
    True
    >>> "trial_debug_file" in generated
    True
    >>> "rkf45_stages.txt" in generated
    True
    >>> "stage_debug_file" in generated
    True
    >>> "Cart_to_xx_and_nearest_i0i1i2_assume_valid" in generated
    True
    >>> "k_bundle[(stage - 1) * 9 + 5]" in generated
    True
    """
    if interpolation_method not in ("g4DD", "g4DD_d0", "GammaUDD"):
        raise ValueError(
            "interpolation_method must be one of ('g4DD', 'g4DD_d0', 'GammaUDD'); "
            f"found '{interpolation_method}'."
        )
    if dataset_coord_system != "SinhCylindricalv2n2":
        raise ValueError(
            "single_integrator_numerical supports only "
            "dataset_coord_system='SinhCylindricalv2n2'; "
            f"found '{dataset_coord_system}'."
        )
    phi_dim = 1

    single_photon_definitions = r"""
    #ifndef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    #endif

    typedef enum {
      TERMINATION_TYPE_COORD_RADIUS_EXCEEDED,
      TERMINATION_TYPE_SOURCE_PLANE,
      FAILURE_EVOLUTION_MEASURE_EXCEEDED,
      FAILURE_RKF45_REJECTION_LIMIT,
      FAILURE_T_MAX_EXCEEDED,
      FAILURE_SLOT_MANAGER_ERROR,
      TERMINATION_TYPE_FAILURE,
      ACTIVE,
      REJECTED
    } termination_type_t; // END ENUM: termination_type_t
    """
    BHaH_defines_h.register_BHaH_defines(
        "single_photon_numerical_definitions", single_photon_definitions
    )

    # Step 1: Register direct single-photon initial-state parameters.
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "initial_t",
            "initial_x",
            "initial_y",
            "initial_z",
            "initial_p_x",
            "initial_p_y",
            "initial_p_z",
            "initial_integration_param",
            "initial_eulerian_distance",
            "initial_h",
        ],
        [28.0, 10.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.05],
        commondata=True,
        add_to_parfile=True,
    )
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

    # Step 2: Select emitted C expressions for the two state conventions.
    if normalized_eom:
        initial_state_time = "commondata.initial_integration_param"
        initial_integration_param = "commondata.initial_t"
        coordinate_time_expression = "*integration_param"
        trajectory_lambda_expression = "f[0]"
        trajectory_time_expression = "*integration_param"
        interpolation_stage_arguments = (
            "&trial_spatial_center.i0, &trial_spatial_center.i2, "
            "integration_param, h, stage,"
        )
        interpolation_initial_arguments = "NULL, NULL, integration_param, h, 1,"
        rhs_integration_arguments = "integration_param, h,"
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
      fprintf(stderr, "ERROR: could not resolve the RKF45 trial spatial center.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    } // END IF: RKF45 trial spatial center was invalid
    trial_spatial_center.i0 = trial_selected_center_idx[0];
    trial_spatial_center.i2 = trial_selected_center_idx[2];
"""
        momentum_conversion_call = (
            "photon_momentum_to_normalized_kernel("
            "f, metric, chunk_size, stream_idx);"
        )
        normalization_kernel_name = "normalization_constraint_photon_normalized"
        normalization_diagnostic_expression = "normalization.C - 1.0"
    else:
        initial_state_time = "commondata.initial_t"
        initial_integration_param = "commondata.initial_integration_param"
        coordinate_time_expression = "f[0]"
        trajectory_lambda_expression = "*integration_param"
        trajectory_time_expression = "f[0]"
        interpolation_stage_arguments = ""
        interpolation_initial_arguments = ""
        rhs_integration_arguments = ""
        trial_spatial_center_setup = ""
        momentum_conversion_call = ""
        normalization_kernel_name = "normalization_constraint_photon"
        normalization_diagnostic_expression = "normalization.C"

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
    "L_Euler"
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
    "L_Euler"
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
        (trial_status_value >= 0 && trial_status_value < 9)
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
step. When RKF45 trial debugging is enabled, ``rkf45_trials.txt`` records every
trial and ``rkf45_stages.txt`` records all six interpolation/RHS stages of every
trial in execution order.

For normalized equations, the state layout is
``(lambda, x, y, z, u, Pi_1, Pi_2, Pi_3, L_Euler)``: ``f[0]`` is lambda and
the RKF45 integration parameter is coordinate time. For non-normalized
equations, the state layout is ``(t, x, y, z, p_0, p_1, p_2, p_3, L_Euler)``:
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
    "TERMINATION_TYPE_COORD_RADIUS_EXCEEDED",
    "TERMINATION_TYPE_SOURCE_PLANE",
    "FAILURE_EVOLUTION_MEASURE_EXCEEDED",
    "FAILURE_RKF45_REJECTION_LIMIT",
    "FAILURE_T_MAX_EXCEEDED",
    "FAILURE_SLOT_MANAGER_ERROR",
    "TERMINATION_TYPE_FAILURE",
    "ACTIVE",
    "REJECTED"
  }};

  int exit_status = EXIT_SUCCESS;
  FILE *trajectory_file = NULL;
{trial_debug_declarations}
{stage_debug_declarations}
  bool slot_manager_initialized = false;
  bool numerical_window_initialized = false;

  double *f = NULL;
  double *f_start = NULL;
  double *f_temp = NULL;
  double *metric = NULL;
  double *rhs_geometry = NULL;
  double *k_bundle = NULL;
  double *integration_param = NULL;
  double *h = NULL;
  int *rejection_retries = NULL;
  termination_type_t *status = NULL;

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
  const double slot_manager_t_max = commondata.initial_t + 1.0e-5;
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

  const REAL camera_position[3] = {{
      (REAL)commondata.initial_x,
      (REAL)commondata.initial_y,
      (REAL)commondata.initial_z}};
  if (time_window_manager_numerical_validate_startup_domain(
          &numerical_window,
          &commondata,
          &numerical_params,
          camera_position,
          "initial_x, initial_y, and initial_z") !=
      TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
    fprintf(stderr,
            "ERROR: numerical camera/r_escape startup stencil validation failed.\n");
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
  }} // END IF: numerical dataset did not contain

  azimuthal_symmetry_spatial_lagrange_context_struct spatial_context;
  spatial_context.stored_phi_samples[0] =
      numerical_params.xxmin{phi_dim} + 0.5 * numerical_params.dxx{phi_dim};
  spatial_context.stored_phi_samples[1] =
      numerical_params.xxmin{phi_dim} + 1.5 * numerical_params.dxx{phi_dim};

  //==========================================
  // 4. DIRECT INITIAL CONDITIONS
  //==========================================
  f[0] = {initial_state_time};
  f[1] = commondata.initial_x;
  f[2] = commondata.initial_y;
  f[3] = commondata.initial_z;
  f[4] = 0.0;
  f[5] = commondata.initial_p_x;
  f[6] = commondata.initial_p_y;
  f[7] = commondata.initial_p_z;
  f[8] = commondata.initial_eulerian_distance;

  *integration_param = {initial_integration_param};
  *h = commondata.initial_h;
  *rejection_retries = 0;
  *status = ACTIVE;

  for (int component = 0; component < 9; ++component) {{
    if (!isfinite(f[component])) {{
      fprintf(stderr, "ERROR: initial photon state was not finite.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: one initial photon state component
  }} // END LOOP: for component over initial photon
  const double spatial_momentum_squared =
      f[5] * f[5] + f[6] * f[6] + f[7] * f[7];
  if (!isfinite(spatial_momentum_squared) || spatial_momentum_squared <= 0.0) {{
    fprintf(stderr, "ERROR: initial spatial photon momentum must be finite and nonzero.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial spatial momentum invalid
  if (!isfinite(*h) || *h == 0.0) {{
    fprintf(stderr, "ERROR: initial_h must be finite and nonzero.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial integration step was invalid

  int mapped_slot_index = -1;
  const int initial_slot_index = slot_get_index(&tsm, commondata.initial_t);
  if (initial_slot_index < 0) {{
    fprintf(
        stderr,
        "ERROR: initial_t=%e is outside the configured TimeSlotManager range.\n",
        (double)commondata.initial_t);
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

  numerical_interpolation(
      &commondata,
      &numerical_params,
      &spatial_context,
      &numerical_window,
      f,
      {interpolation_initial_arguments}
      metric,
      NULL,
      chunk_size,
      stream_idx);

  for (int component = 0; component < 10; ++component) {{
    if (!isfinite(metric[component])) {{
      fprintf(
          stderr,
          "ERROR: initial numerical metric component %d was not finite.\n",
          component);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: initial metric component invalid
  }} // END LOOP: for component over initial metric

  p0_reverse_kernel(f, metric, chunk_size, stream_idx);
  {momentum_conversion_call}

  for (int component = 0; component < 9; ++component) {{
    if (!isfinite(f[component])) {{
      fprintf(
          stderr,
          "ERROR: initial constrained state component %d was not finite.\n",
          component);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: one initial constrained state component
  }} // END LOOP: for component over the initial

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
      "# lambda t x y z energy_measure p_x p_y p_z aux norm\n");
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
      *status = FAILURE_T_MAX_EXCEEDED;
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
    }} // END IF: photon moved to a different

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
          {interpolation_stage_arguments}
          metric,
          rhs_geometry,
          chunk_size,
          stream_idx);

      for (int component = 0; component < 10; ++component) {{
        if (!isfinite(metric[component])) {{
          fprintf(
              stderr,
              "ERROR: stage %d metric component %d was not finite at t=%e.\n",
              stage,
              component,
              coordinate_time);
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: stage metric component invalid
      }} // END LOOP: for component over stage metric

      for (int component = 0; component < 40; ++component) {{
        if (!isfinite(rhs_geometry[component])) {{
          fprintf(
              stderr,
              "ERROR: stage %d geometry component %d was not finite "
              "at t=%e.\n",
              stage,
              component,
              coordinate_time);
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: stage geometry component invalid
      }} // END LOOP: for component over stage geometry

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
      }} // END LOOP: for component over the accepted

      numerical_interpolation(
          &commondata,
          &numerical_params,
          &spatial_context,
          &numerical_window,
          f,
          {interpolation_initial_arguments}
          metric,
          NULL,
          chunk_size,
          stream_idx);

      for (int component = 0; component < 10; ++component) {{
        if (!isfinite(metric[component])) {{
          fprintf(
              stderr,
              "ERROR: accepted-state metric component %d was not finite.\n",
              component);
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: accepted-state metric component invalid
      }} // END LOOP: for component over accepted-state metric

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
          "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
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
          trajectory_norm);
      fflush(trajectory_file);
      accepted_steps++;

      const double radius_squared = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
      if (radius_squared > commondata.r_escape * commondata.r_escape) {{
        *status = TERMINATION_TYPE_COORD_RADIUS_EXCEEDED;
        printf("Photon escaped to r > %.15e.\n", (double)commondata.r_escape);
        break;
      }} // END IF: photon state crossed boundary

      // f[4] is p^0 for direct evolution and the normalized log-energy measure otherwise.
      if (fabs(f[4]) > commondata.evolution_measure_max) {{
        *status = FAILURE_EVOLUTION_MEASURE_EXCEEDED;
        printf("Evolution measure exceeded %.15e.\n", commondata.evolution_measure_max);
        break;
      }} // END IF: evolution measure exceeded limit
    }} else if (*status == REJECTED)
      continue;
    else if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
      printf("RKF45 reached its consecutive-rejection limit.\n");
      break;
    }} else {{
      fprintf(stderr, "ERROR: unexpected integration status %d.\n", (int)*status);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END ELSE: RKF45 finalization returned an unexpected
  }} // END WHILE: evolve one photon through accepted

  if ((*status == ACTIVE || *status == REJECTED) &&
      accepted_steps >= max_accepted_steps) {{
    *status = TERMINATION_TYPE_FAILURE;
    printf("Integration stopped at the accepted-step safety limit.\n");
  }} else if ((*status == ACTIVE || *status == REJECTED) &&
             rkf45_attempts >= max_rkf45_attempts) {{
    *status = TERMINATION_TYPE_FAILURE;
    printf("Integration stopped at the RKF45-attempt safety limit.\n");
  }} // END ELSE IF: integration reached the RKF45-attempt safety

  const int final_status_index = (int)*status;
  const char *final_status_name =
      (final_status_index >= 0 && final_status_index < 9)
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
      if (terminal_slot_index != mapped_slot_index) {{
        if (time_window_manager_numerical_mmap_for_slot(
                &numerical_window, &tsm, terminal_slot_index) !=
            TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
          fprintf(
              stderr,
              "ERROR: failed to map the terminal normalization time window.\n");
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: terminal normalization time-window mapping failed
        mapped_slot_index = terminal_slot_index;
      }} // END IF: terminal state occupied a different

      numerical_interpolation(
          &commondata,
          &numerical_params,
          &spatial_context,
          &numerical_window,
          f,
          {interpolation_initial_arguments}
          metric,
          NULL,
          chunk_size,
          stream_idx);

      for (int component = 0; component < 10; ++component) {{
        if (!isfinite(metric[component])) {{
          fprintf(
              stderr,
              "ERROR: terminal metric component %d was not finite.\n",
              component);
          exit_status = EXIT_FAILURE;
          goto cleanup;
        }} // END IF: terminal metric component invalid
      }} // END LOOP: for component over terminal metric

      normalization_constraint_t normalization;
      {normalization_kernel_name}(
          f, metric, &normalization, chunk_size, stream_idx);
      const double normalization_deviation = {normalization_diagnostic_expression};
      if (!isfinite(normalization_deviation)) {{
        fprintf(stderr, "ERROR: terminal normalization deviation was not finite.\n");
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: terminal normalization deviation invalid
      printf(
          "Final signed normalization deviation: %.15e\n",
          normalization_deviation);
    }} // END ELSE: terminal state inside window
  }} // END IF: terminal normalization diagnostics were requested

  cleanup:
{trial_debug_cleanup}
{stage_debug_cleanup}
  if (trajectory_file != NULL)
    fclose(trajectory_file);
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
    if geodesic_data.p0_photon is None:
        raise ValueError(f"p0_photon is None for {GEO_KEY}")
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

    p0_reverse_kernel.p0_reverse_kernel(geodesic_data.p0_photon)
    normalization_constraint.normalization_constraint(
        geodesic_data.norm_constraint_expr, PARTICLE
    )
    numerical_interpolation.register_CFunction_numerical_interpolation(
        DATASET_COORD_SYSTEM,
        interpolation_method=INTERPOLATION_METHOD,
        enable_simd=False,
        project_dir=project_dir,
        normalized_eom=NORMALIZED_EOM,
    )
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_rhs,
        geodesic_data.xx,
        rhs_uses_metric_derivatives=rhs_uses_metric_derivatives,
        normalized_eom=NORMALIZED_EOM,
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True,
        register_numerical_initial_h=False,
    )
    single_integrator_numerical(
        SPACETIME,
        DATASET_COORD_SYSTEM,
        interpolation_method=INTERPOLATION_METHOD,
        normalized_eom=NORMALIZED_EOM,
    )
    main_single.main_single("single_integrator_numerical")

    # Step 7: Generate parameter headers, the default parfile, and CPU definitions.
    print("Generating headers, parameter handling, and Makefile...")
    CPs.write_CodeParameters_h_files(set_commondata_only=True, project_dir=project_dir)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=PROJECT_NAME
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=PROJECT_NAME
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

    # Step 8: Copy the trajectory visualizer when this script resides in nrpy/examples.
    script_dir = os.path.dirname(os.path.abspath(__file__))
    visualization_source = os.path.join(
        script_dir, "geodesic_visualizations", "visualize_trajectory.py"
    )
    if os.path.exists(visualization_source):
        shutil.copy(visualization_source, project_dir)
    else:
        print(
            "Warning: visualize_trajectory.py was not found beside this example; "
            "the integrator project was still generated."
        )

    print(f"Finished generating {project_dir}.")
    print(f"Set numerical_spacetime_bin_path in {parfile_path} before running.")
    print(f"Build with: cd {project_dir} && make")
    print(f"Run with: ./{PROJECT_NAME} {PROJECT_NAME}.par")
    print("Accepted trajectory samples will be written to trajectory.txt.")

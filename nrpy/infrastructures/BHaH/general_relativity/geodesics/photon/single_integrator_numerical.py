"""
Generate a standalone CPU project for one photon in a numerical spacetime.

The generated C executable follows the numerical batch integrator's numerical
interpolation, time-window management, initial-condition geometry, and RKF45
control logic, but removes batching, double buffering, device memory, streams,
and blueprint output. It writes one trajectory row after every accepted RKF45
step.

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
) -> None:
    """
    Register the standalone numerical-spacetime single-photon C integrator.

    :param spacetime_name: Spacetime identifier used to select photon equations.
    :param dataset_coord_system: Coordinate system used by the numerical dataset.
    :param interpolation_method: Numerical geometry payload method used by the generated project.
    :param normalized_eom: Whether to evolve normalized coordinate-time equations.
    :raises ValueError: If the interpolation method or dataset coordinate system
        is unsupported.
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

    # Step 1: Register the batch-authoritative single-ray input parameters.
    par.register_CodeParameter(
        "int",
        __name__,
        "scan_density",
        500,
        commondata=True,
        add_to_parfile=True,
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        ["window_center_x", "window_center_y", "window_center_z"],
        [50.0, 0.0, 0.0],
        commondata=True,
        add_to_parfile=False,
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "original_window_center_x",
            "original_window_center_y",
            "original_window_center_z",
            "camera_pos_x",
            "camera_pos_y",
            "camera_pos_z",
            "window_up_vec_x",
            "window_up_vec_y",
            "window_up_vec_z",
            "window_width",
            "window_height",
            "t_start",
        ],
        [50.0, 0.0, 0.0, 51.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 100.0],
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
        initial_state_time = "0.0"
        initial_integration_param = "commondata.t_start"
        coordinate_time_expression = "*integration_param"
        trajectory_lambda_expression = "f[0]"
        trajectory_time_expression = "*integration_param"
        interpolation_stage_arguments = "integration_param, h, stage,"
        interpolation_initial_arguments = "integration_param, h, 1,"
        rhs_integration_arguments = "integration_param, h,"
        momentum_conversion_call = (
            "photon_momentum_to_normalized_kernel("
            "f, metric, chunk_size, stream_idx);"
        )
        normalization_kernel_name = "normalization_constraint_photon_normalized"
        normalization_error_expression = "fabs(normalization.C - 1.0)"
    else:
        initial_state_time = "commondata.t_start"
        initial_integration_param = "0.0"
        coordinate_time_expression = "f[0]"
        trajectory_lambda_expression = "*integration_param"
        trajectory_time_expression = "f[0]"
        interpolation_stage_arguments = ""
        interpolation_initial_arguments = ""
        rhs_integration_arguments = ""
        momentum_conversion_call = ""
        normalization_kernel_name = "normalization_constraint_photon"
        normalization_error_expression = "fabs(normalization.C)"

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

The executable initializes ray zero with the same camera-window geometry used by
the numerical batch integrator, solves the initial null constraint, maps numerical
time windows by coordinate-time slot, and advances the state with the shared
six-stage RKF45 pipeline. A trajectory row is written only after an accepted step.

For normalized equations, the output lambda column is f[0] and the output time
column is the RKF45 integration parameter. For non-normalized equations, lambda is
the integration parameter and time is f[0].

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

  // The batch driver seeds each active tile from the configured master window.
  commondata.window_center_x = commondata.original_window_center_x;
  commondata.window_center_y = commondata.original_window_center_y;
  commondata.window_center_z = commondata.original_window_center_z;

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
  }} // END IF: one or more single-photon CPU allocations failed

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
          commondata.numerical_spacetime_temporal_interp_order,
          &numerical_params) != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {{
    fprintf(
        stderr,
        "ERROR: failed to initialize the numerical time-window manager from '%s'.\n",
        commondata.numerical_spacetime_bin_path);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: numerical time-window manager initialization failed
  numerical_window_initialized = true;

  if (numerical_params.Nxx{phi_dim} != 2) {{
    fprintf(
        stderr,
        "ERROR: numerical spatial interpolation expects exactly two stored phi "
        "planes in native dimension {phi_dim}; got Nxx{phi_dim}=%d.\n",
        numerical_params.Nxx{phi_dim});
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: numerical dataset did not contain two azimuthal planes

  azimuthal_symmetry_spatial_lagrange_context_struct spatial_context;
  spatial_context.stored_phi_samples[0] =
      numerical_params.xxmin{phi_dim} + 0.5 * numerical_params.dxx{phi_dim};
  spatial_context.stored_phi_samples[1] =
      numerical_params.xxmin{phi_dim} + 1.5 * numerical_params.dxx{phi_dim};

  //==========================================
  // 4. BATCH-AUTHORITATIVE INITIAL CONDITIONS
  //==========================================
  const double camera_x = commondata.camera_pos_x;
  const double camera_y = commondata.camera_pos_y;
  const double camera_z = commondata.camera_pos_z;

  double n_z[3] = {{
    commondata.original_window_center_x - camera_x,
    commondata.original_window_center_y - camera_y,
    commondata.original_window_center_z - camera_z
  }};
  const double magnitude_n_z =
      sqrt(n_z[0] * n_z[0] + n_z[1] * n_z[1] + n_z[2] * n_z[2]);
  if (!isfinite(magnitude_n_z) || magnitude_n_z <= 0.0) {{
    fprintf(stderr, "ERROR: camera and original window center must differ.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: camera line-of-sight vector was invalid
  for (int component = 0; component < 3; ++component)
    n_z[component] /= magnitude_n_z;

  const double guide_up[3] = {{
    commondata.window_up_vec_x,
    commondata.window_up_vec_y,
    commondata.window_up_vec_z
  }};
  double n_x[3] = {{
    guide_up[1] * n_z[2] - guide_up[2] * n_z[1],
    guide_up[2] * n_z[0] - guide_up[0] * n_z[2],
    guide_up[0] * n_z[1] - guide_up[1] * n_z[0]
  }};
  double magnitude_n_x =
      sqrt(n_x[0] * n_x[0] + n_x[1] * n_x[1] + n_x[2] * n_x[2]);

  if (magnitude_n_x < 1.0e-9) {{
    double alternative_up[3] = {{0.0, 1.0, 0.0}};
    if (fabs(n_z[1]) > 0.999) {{
      alternative_up[1] = 0.0;
      alternative_up[2] = 1.0;
    }} // END IF: line of sight was parallel to the first fallback axis
    n_x[0] = alternative_up[1] * n_z[2] - alternative_up[2] * n_z[1];
    n_x[1] = alternative_up[2] * n_z[0] - alternative_up[0] * n_z[2];
    n_x[2] = alternative_up[0] * n_z[1] - alternative_up[1] * n_z[0];
    magnitude_n_x =
        sqrt(n_x[0] * n_x[0] + n_x[1] * n_x[1] + n_x[2] * n_x[2]);
  }} // END IF: camera up vector required the batch fallback basis

  if (!isfinite(magnitude_n_x) || magnitude_n_x <= 0.0) {{
    fprintf(stderr, "ERROR: failed to construct the camera horizontal basis.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: camera horizontal basis remained singular
  for (int component = 0; component < 3; ++component)
    n_x[component] /= magnitude_n_x;

  const double n_y[3] = {{
    n_z[1] * n_x[2] - n_z[2] * n_x[1],
    n_z[2] * n_x[0] - n_z[0] * n_x[2],
    n_z[0] * n_x[1] - n_z[1] * n_x[0]
  }};

  if (commondata.scan_density <= 0) {{
    fprintf(stderr, "ERROR: scan_density must be positive.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: scan_density was not positive

  const long int ray_index = 0;
  const int row = (int)(ray_index / commondata.scan_density);
  const int column = (int)(ray_index % commondata.scan_density);
  const double pixel_x =
      -0.5 * commondata.window_width +
      (column + 0.5) * (commondata.window_width / commondata.scan_density);
  const double pixel_y =
      -0.5 * commondata.window_height +
      (row + 0.5) * (commondata.window_height / commondata.scan_density);

  const double target_position[3] = {{
    commondata.window_center_x + pixel_x * n_x[0] + pixel_y * n_y[0],
    commondata.window_center_y + pixel_x * n_x[1] + pixel_y * n_y[1],
    commondata.window_center_z + pixel_x * n_x[2] + pixel_y * n_y[2]
  }};
  const double direction_x = target_position[0] - camera_x;
  const double direction_y = target_position[1] - camera_y;
  const double direction_z = target_position[2] - camera_z;
  const double direction_magnitude = sqrt(
      direction_x * direction_x + direction_y * direction_y +
      direction_z * direction_z);
  if (!isfinite(direction_magnitude) || direction_magnitude <= 0.0) {{
    fprintf(stderr, "ERROR: initial photon direction was invalid.\n");
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial photon direction was invalid

  f[0] = {initial_state_time};
  f[1] = camera_x;
  f[2] = camera_y;
  f[3] = camera_z;
  f[4] = 0.0;
  f[5] = direction_x / direction_magnitude;
  f[6] = direction_y / direction_magnitude;
  f[7] = direction_z / direction_magnitude;
  f[8] = 0.0;

  *integration_param = {initial_integration_param};
  *h = commondata.numerical_initial_h;
  *rejection_retries = 0;
  *status = ACTIVE;

  int mapped_slot_index = -1;
  const int initial_slot_index = slot_get_index(&tsm, commondata.t_start);
  if (initial_slot_index < 0) {{
    fprintf(
        stderr,
        "ERROR: t_start=%e is outside the configured TimeSlotManager range.\n",
        (double)commondata.t_start);
    exit_status = EXIT_FAILURE;
    goto cleanup;
  }} // END IF: initial coordinate time was outside the slot lattice

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
    }} // END IF: one initial metric component was not finite
  }} // END LOOP: for component over initial metric components

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
    }} // END IF: one initial constrained state component was not finite
  }} // END LOOP: for component over the initial constrained state

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
  }} // END IF: trajectory output file could not be opened
  fprintf(trajectory_file, "# lambda t x y z energy_measure p_x p_y p_z aux\n");

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
    }} // END IF: current coordinate time left the slot lattice

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
      }} // END IF: numerical time-window mapping failed during integration
      mapped_slot_index = slot_index;
    }} // END IF: photon moved to a different coordinate-time slot

    memcpy(f_start, f, sizeof(double) * 9);
    memcpy(f_temp, f, sizeof(double) * 9);

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
        }} // END IF: one stage metric component was not finite
      }} // END LOOP: for component over stage metric components

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
        }} // END IF: one stage geometry component was not finite
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
      rkf45_stage_update(
          f_start,
          k_bundle,
          h,
          stage,
          chunk_size,
          f_temp,
          stream_idx);
    }} // END LOOP: for stage over all six Cash-Karp stages

    rkf45_finalize_and_control(
        &commondata,
        f,
        f_start,
        k_bundle,
        h,
        status,
        integration_param,
        rejection_retries,
        chunk_size,
        stream_idx);
    rkf45_attempts++;

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
        }} // END IF: one accepted state component was not finite
      }} // END LOOP: for component over the accepted photon state

      fprintf(
          trajectory_file,
          "%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          {trajectory_lambda_expression},
          {trajectory_time_expression},
          f[1],
          f[2],
          f[3],
          f[4],
          f[5],
          f[6],
          f[7],
          f[8]);
      fflush(trajectory_file);
      accepted_steps++;

      const double radius_squared = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
      if (radius_squared > commondata.r_escape * commondata.r_escape) {{
        *status = TERMINATION_TYPE_COORD_RADIUS_EXCEEDED;
        printf("Photon escaped to r > %.15e.\n", (double)commondata.r_escape);
        break;
      }} // END IF: accepted photon state crossed the escape sphere

      // f[4] is p^0 for direct evolution and the normalized log-energy measure otherwise.
      if (fabs(f[4]) > commondata.evolution_measure_max) {{
        *status = FAILURE_EVOLUTION_MEASURE_EXCEEDED;
        printf("Evolution measure exceeded %.15e.\n", commondata.evolution_measure_max);
        break;
      }} // END IF: accepted evolution measure exceeded the configured limit
    }} else if (*status == REJECTED)
      continue;
    else if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
      printf("RKF45 reached its consecutive-rejection limit.\n");
      break;
    }} else {{
      fprintf(stderr, "ERROR: unexpected integration status %d.\n", (int)*status);
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END ELSE: RKF45 finalization returned an unexpected terminal status
  }} // END WHILE: evolve one photon through accepted and rejected RKF45 attempts

  if ((*status == ACTIVE || *status == REJECTED) &&
      accepted_steps >= max_accepted_steps) {{
    *status = TERMINATION_TYPE_FAILURE;
    printf("Integration stopped at the accepted-step safety limit.\n");
  }} else if ((*status == ACTIVE || *status == REJECTED) &&
             rkf45_attempts >= max_rkf45_attempts) {{
    *status = TERMINATION_TYPE_FAILURE;
    printf("Integration stopped at the RKF45-attempt safety limit.\n");
  }} // END ELSE IF: integration reached the RKF45-attempt safety limit

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
      }} // END IF: terminal state occupied a different coordinate-time slot

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
        }} // END IF: one terminal metric component was not finite
      }} // END LOOP: for component over terminal metric components

      normalization_constraint_t normalization;
      {normalization_kernel_name}(
          f, metric, &normalization, chunk_size, stream_idx);
      const double normalization_error = {normalization_error_expression};
      if (!isfinite(normalization_error)) {{
        fprintf(stderr, "ERROR: terminal normalization error was not finite.\n");
        exit_status = EXIT_FAILURE;
        goto cleanup;
      }} // END IF: terminal normalization error was not finite
      printf("Final normalization absolute error: %.15e\n", normalization_error);
    }} // END ELSE: terminal state was inside the numerical time domain
  }} // END IF: terminal normalization diagnostics were requested

  cleanup:
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
        geodesic_data.geodesic_rhs,
        geodesic_data.xx,
        normalized_eom=NORMALIZED_EOM,
        interpolation_method=INTERPOLATION_METHOD,
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True
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

# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/single_integrator_analytical.py
"""
Define the reusable single-photon analytical geodesic integrator orchestrator.

This module registers the C function that evolves one massless test particle from a
symbolic spacetime recipe using the split Runge-Kutta-Fehlberg 4(5) photon pipeline. Its
initial state is constructed from the metric-driven observer tetrad shared with
the numerical integrator. The registered function writes trajectory samples and
reports normalization and conserved-quantity diagnostics while preserving the
Structure of Arrays layout expected by the shared geodesic kernels. Optional
terminal and nonterminal planes record crossing time, affine parameter, local
plane coordinates, and all nine interpolated state components in
``plane_crossings.txt``. Optional RKF45 debugging writes trial-controller
records and analytic metric, connection, and right-hand-side records from the
single-photon integration loop.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon.set_initial_conditions_kernel import (
    register_photon_batch_structs,
)


def single_integrator_analytical(
    spacetime: str,
    particle: str,
    maximum_degree: int,
    normalized_eom: bool = False,
    enable_rkf45_trial_debug: bool = False,
) -> None:
    """
    Register the single-photon analytical geodesic integrator C function.

    The generated C function initializes one photon state, evolves it with the
    split RKF45 pipeline, writes trajectory samples, and reports final
    normalization and conserved-quantity diagnostics. Optional plane crossings
    are reconstructed after accepted steps. ``plane_crossings.txt`` records
    crossing time, affine parameter, local plane coordinates, and all nine
    interpolated state components.

    :param spacetime: The background spacetime descriptor.
    :param particle: The test-particle type.
    :param maximum_degree: Largest generated plane interpolation polynomial degree.
    :param normalized_eom: Whether to use normalized photon evolution.
    :param enable_rkf45_trial_debug: Whether to write one diagnostic row for every
        RKF45 trial to ``rkf45_trials.txt`` and one analytic metric, connection,
        and right-hand-side row for each of its six stages to
        ``rkf45_stages.txt``.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> import os
    >>> import tempfile
    >>> from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import rkf45_finalize_and_control_kernel
    >>> cache_dir = tempfile.TemporaryDirectory()
    >>> os.environ["XDG_CACHE_HOME"] = cache_dir.name
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
    ...     normalized_eom=True, enable_rkf45_trial_debug=True
    ... )
    >>> single_integrator_analytical(
    ...     "BrillLindquist_InitialData_Static_Cartesian",
    ...     "photon",
    ...     4,
    ...     normalized_eom=True,
    ...     enable_rkf45_trial_debug=True,
    ... )
    >>> generated = cfc.CFunction_dict["single_integrator_analytical"].full_function
    >>> ("rkf45_trials.txt" in generated and "rkf45_stages.txt" in generated and
    ...  "# lambda t x y z u Pi_1 Pi_2 Pi_3 L_normal" in generated)
    True
    >>> cache_dir.cleanup()
    """
    if maximum_degree < 3:
        raise ValueError("Plane interpolation degree must be at least three")
    register_photon_batch_structs()

    macro_defs = r"""
    #undef BUNDLE_CAPACITY
    #define BUNDLE_CAPACITY 1
    """
    BHaH_defines_h.register_BHaH_defines("single_photon_macros", macro_defs)

    # The shared initializer expects the batch tile-sampling parameters. This
    # standalone path is fixed to one tile containing one ray; the active
    # indices and scan density are registered by set_initial_conditions_kernel.
    par.register_CodeParameters(
        "int",
        __name__,
        ["tiles_width", "tiles_height"],
        [1, 1],
        commondata=True,
        add_to_parfile=False,
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "r_escape",
        ],
        [150.0],
        commondata=True,
        add_to_parfile=True,
    )

    includes = [
        "math.h",
        "stdio.h",
        "stdlib.h",
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]

    desc = """Integrate one photon geodesic with the split RKF45 pipeline.

Initializes one photon state, evolves it with the split RKF45 pipeline,
writes trajectory samples, and reports final normalization and
conserved-quantity diagnostics. Optional terminal and nonterminal planes record
crossing time, affine parameter, local plane coordinates, and all nine
interpolated state components in ``plane_crossings.txt``.

When RKF45 trial debugging is enabled, ``rkf45_trials.txt`` records every
adaptive-step trial and ``rkf45_stages.txt`` records all six analytic metric,
connection, and right-hand-side stages of each trial in execution order.

@return EXIT_SUCCESS on successful completion; EXIT_FAILURE if setup fails.
"""

    cfunc_type = "int"
    name = "single_integrator_analytical"
    params = "int argc, const char *argv[]"

    if spacetime == "KerrSchild_Cartesian":
        conserved_quantity_q_report = r"""
    printf("  Delta Q  = %.4e\n", fabs(cq_final.Q - cq_init.Q));
"""
    else:
        conserved_quantity_q_report = r"""
    printf("  Delta Q  = not defined for this spacetime.\n");
"""

    if normalized_eom:
        conserved_quantity_q_report = ""
        conservation_error_report = (
            '    printf("Conservation Absolute Errors are not reported for '
            'normalized EOM.\\n");'
        )
        conserved_quantity_initialization = ""
        conserved_quantity_finalization = ""
        initial_integration_parameter = "commondata.t_start"
        trajectory_lambda_expression = "f[0]"
        trajectory_time_expression = "*integration_param"
        trajectory_header = "# lambda t x y z u Pi_1 Pi_2 Pi_3 L_normal\\n"
        rhs_integration_arguments = "integration_param, h,"
        log_energy_evaluation = "const double log_energy_measure = f[4];"
        normalization_kernel_name = "normalization_constraint_photon_normalized"
        normalization_error_expression = "fabs(norm_final.C - 1.0)"
        initial_normalization_check = r"""
    normalization_constraint_t norm_initial;
    normalization_constraint_photon_normalized(
      f, metric, &norm_initial, chunk_size, stream_idx
    );
    if (!isfinite(norm_initial.C) || fabs(norm_initial.C - 1.0) > 1.0e-9) {
      fprintf(stderr, "ERROR: initial normalized photon constraint is invalid.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }
"""
    else:
        conservation_error_report = (
            '    printf("Conservation Absolute Errors:\\n");\n'
            '    printf("  Delta E  = %.4e\\n", fabs(cq_final.E - cq_init.E));\n'
            '    printf("  Delta Lz = %.4e\\n", fabs(cq_final.Lz - cq_init.Lz));'
        )
        conserved_quantity_initialization = rf"""
    conserved_quantities_t cq_init;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(
      &commondata, &all_photons, num_rays, &cq_init
    );
"""
        conserved_quantity_finalization = rf"""
    conserved_quantities_t cq_final;
    calculate_conserved_quantities_universal_{spacetime}_{particle}(
      &commondata, &all_photons, num_rays, &cq_final
    );
"""
        initial_integration_parameter = "0.0"
        trajectory_lambda_expression = "*integration_param"
        trajectory_time_expression = "f[0]"
        trajectory_header = "# lambda t x y z p^t p^x p^y p^z L_normal\\n"
        rhs_integration_arguments = ""
        log_energy_evaluation = r"""
        normal_observer_log_energy(
          f, metric, log_energy_bundle, chunk_size, stream_idx
        );
        const double log_energy_measure = log_energy_bundle[0];
"""
        normalization_kernel_name = "normalization_constraint_photon"
        normalization_error_expression = "fabs(norm_final.C)"
        initial_normalization_check = ""

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

    stage_normalization_diagnostic_expression = (
        "stage_normalization.C - 1.0" if normalized_eom else "stage_normalization.C"
    )
    trial_start_time_expression = (
        "*integration_param" if normalized_eom else "f_base[0]"
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
""".replace("{trial_component_names}", trial_component_names)
        stage_debug_declarations = (
            r"""
    FILE *stage_debug_file = NULL;
    const double rkf45_stage_time_fractions[] = {
      0.0, 1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0};
"""
            if normalized_eom
            else r"""
    FILE *stage_debug_file = NULL;
"""
        )
        stage_debug_header = (
            """# accepted_step trial_number retry_number stage h_trial t_start stage_time lambda x y z r_stage stage_norm_error u Pi_1 Pi_1_derivative
"""
            if normalized_eom
            else """# accepted_step trial_number retry_number stage h_trial t_start stage_time t x y z r_stage stage_norm_error p^0 p^1 p^1_derivative
"""
        )
        stage_debug_header_c = (
            stage_debug_header.rstrip("\n").replace("\\", "\\\\").replace('"', '\\"')
            + "\\n"
        )
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
        stage_debug_open = rf"""
    stage_debug_file = fopen("rkf45_stages.txt", "w");
    if (stage_debug_file == NULL) {{
      fprintf(stderr, "ERROR: could not open rkf45_stages.txt for writing.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: RKF45 stage diagnostics unavailable
    fprintf(stage_debug_file, "{stage_debug_header_c}");
"""
        trial_debug_trial_metadata = rf"""
      const long int accepted_step_before_trial = steps;
      const long int trial_number = rkf45_attempts + 1;
      const int retry_number_before = *rejection_retries;
      const double t_start = {trial_start_time_expression};
      const double h_trial = *h;
      const double x_start = f_base[1];
      const double y_start = f_base[2];
      const double z_start = f_base[3];
      const double r_start = sqrt(
        x_start * x_start + y_start * y_start + z_start * z_start);
"""
        stage_debug_record = rf"""
        normalization_constraint_t stage_normalization;
        {normalization_kernel_name}(
          f_temp, metric, &stage_normalization, chunk_size, stream_idx
        );
        const double stage_norm_error =
          {stage_normalization_diagnostic_expression};
        // Preserve nonfinite diagnostic values in the stage record. The
        // RKF45 finalizer handles a nonfinite candidate through its ordinary
        // rejection path, so enabling diagnostics does not change evolution.

        const double stage_radius = sqrt(
          f_temp[1] * f_temp[1] +
          f_temp[2] * f_temp[2] +
          f_temp[3] * f_temp[3]);
        const double stage_time = {stage_time_expression};
        const double stage_p1_derivative = k_bundle[(stage - 1) * 9 + 5];
        fprintf(
          stage_debug_file,
          "%ld %ld %d %d "
          "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g "
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
          stage_norm_error,
          f_temp[4],
          f_temp[5],
          stage_p1_derivative);
        fflush(stage_debug_file);
"""
        trial_debug_call_argument = "&trial_debug,\n        "
        trial_debug_attempt_declaration = "    long int rkf45_attempts = 0;\n"
        trial_debug_attempt_increment = "      rkf45_attempts++;\n"
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
        stage_debug_declarations = ""
        trial_debug_open = ""
        stage_debug_open = ""
        trial_debug_trial_metadata = ""
        stage_debug_record = ""
        trial_debug_call_argument = ""
        trial_debug_attempt_declaration = ""
        trial_debug_attempt_increment = ""
        trial_debug_record = ""
        trial_debug_cleanup = ""
        stage_debug_cleanup = ""

    normalized_momentum_conversion = (
        """    photon_momentum_to_normalized_kernel(
      f, metric, chunk_size
    );"""
        if normalized_eom
        else ""
    )
    initial_state_report = (
        '    printf("  Normalized state (u=%.4f, Pi_1=%.4f, Pi_2=%.4f, Pi_3=%.4f)\\n", f[4], f[5], f[6], f[7]);'
        if normalized_eom
        else '    printf("  Mom (%.4f, %.4f, %.4f, %.4f)\\n", f[4], f[5], f[6], f[7]);'
    )
    final_parameter_report = (
        """    printf(
      "Integration finished after %d steps. Final lambda = %.4f, final coordinate time = %.4f\\n",
      steps,
      f[0],
      *integration_param
    );"""
        if normalized_eom
        else """    printf(
      "Integration finished after %d steps. Final lambda = %.4f\\n",
      steps,
      *integration_param
    );"""
    )

    body = rf"""
    // ==========================================
    // STRUCTURAL SETUP & PARAMETERS
    // ==========================================
    commondata_struct commondata;
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

    const long int num_rays = 1;
    const long int chunk_size = 1;
    const int stream_idx = 0;
    // These values match the shared termination_type_t through REJECTED.
    // Numerical-only interpolation failures are omitted because this path
    // evaluates the analytic metric and connection directly. ACTIVE and
    // REJECTED are intentionally 7 and 8, not the first two enum values.
    const char *status_names[] = {{
      "STOP_CONDITION_COORD_RADIUS_EXCEEDED",
      "STOP_CONDITION_TERMINAL_PLANE",
      "STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED",
      "FAILURE_RKF45_REJECTION_LIMIT",
      "STOP_CONDITION_T_MAX_EXCEEDED",
      "FAILURE_SLOT_MANAGER_ERROR",
      "FAILURE_GENERIC",
      "ACTIVE",
      "REJECTED"
    }};

    int exit_status = EXIT_SUCCESS;
    FILE *fp = NULL;
    FILE *plane_crossings_file = NULL;
{trial_debug_declarations}
{stage_debug_declarations}
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
    double f_event_history[{maximum_degree} * 9] = {{0.0}};
    double integration_param_event_history[{maximum_degree}] = {{0.0}};
    blueprint_data_t plane_crossing_results = {{0}};
    const long int event_chunk_index[1] = {{0}};
    const bool event_planes_enabled =
        commondata.non_terminal_plane_enabled || commondata.terminal_plane_enabled;

    printf("Starting Split-Pipeline Geodesic Integrator...\n");
    printf("spacetime: {spacetime}\n");

    // ==========================================
    // GLOBAL MEMORY ALLOCATION
    // ==========================================
    double *f = NULL;
    double *f_base = NULL;
    double *f_temp = NULL;
    double *metric = NULL;
    double *connection = NULL;
    double *k_bundle = NULL;
    double *integration_param = NULL;
    double *h = NULL;
    int *rejection_retries = NULL;
    termination_type_t *status = NULL;
    double *log_energy_bundle = NULL;

    BHAH_MALLOC(f, 9 * sizeof(double));
    BHAH_MALLOC(f_base, 9 * sizeof(double));
    BHAH_MALLOC(f_temp, 9 * sizeof(double));
    BHAH_MALLOC(metric, 10 * sizeof(double));
    BHAH_MALLOC(connection, 40 * sizeof(double));
    BHAH_MALLOC(k_bundle, 6 * 9 * sizeof(double));
    BHAH_MALLOC(integration_param, sizeof(double));
    BHAH_MALLOC(h, sizeof(double));
    BHAH_MALLOC(rejection_retries, sizeof(int));
    BHAH_MALLOC(status, sizeof(termination_type_t));
    BHAH_MALLOC(log_energy_bundle, sizeof(double));

    if (f == NULL || f_base == NULL || f_temp == NULL || metric == NULL ||
        connection == NULL || k_bundle == NULL || integration_param == NULL ||
        h == NULL || rejection_retries == NULL || status == NULL ||
        log_energy_bundle == NULL) {{
      fprintf(stderr, "Error: failed to allocate photon state buffers.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: photon buffer allocation failed

    PhotonStateSoA all_photons = {{0}};
    all_photons.f = f;
    all_photons.h = h;
    all_photons.integration_param = integration_param;
    all_photons.f_history = f_event_history;
    all_photons.integration_param_history = integration_param_event_history;
    all_photons.on_positive_side_of_non_terminal_plane_prev =
        &on_positive_side_of_non_terminal_plane_prev;
    all_photons.on_positive_side_of_terminal_plane_prev =
        &on_positive_side_of_terminal_plane_prev;
    all_photons.non_terminal_plane_crossing_pending =
        &non_terminal_plane_crossing_pending;
    all_photons.terminal_plane_crossing_pending =
        &terminal_plane_crossing_pending;
    all_photons.non_terminal_plane_steps_past = &non_terminal_plane_steps_past;
    all_photons.terminal_plane_steps_past = &terminal_plane_steps_past;
    all_photons.non_terminal_plane_event_degree =
        &non_terminal_plane_event_degree;
    all_photons.terminal_plane_event_degree = &terminal_plane_event_degree;

    // ==========================================
    // INITIAL CONDITIONS
    // ==========================================
    f[0] = commondata.t_start;
    f[1] = commondata.observer_x;
    f[2] = commondata.observer_y;
    f[3] = commondata.observer_z;
    f[4] = 0.0;
    f[5] = 0.0;
    f[6] = 0.0;
    f[7] = 0.0;
    f[8] = 0.0;

    // Single-ray execution is the one-tile, one-sample specialization of the
    // shared angular sample mapping. Tile origins and pixel dimensions are
    // deliberately not part of commondata.
    commondata.tiles_width = 1;
    commondata.tiles_height = 1;
    commondata.tile_index_width = 0;
    commondata.tile_index_height = 0;
    commondata.scan_density = 1;
    *integration_param = {initial_integration_parameter};
    *h = commondata.initial_h;
    *rejection_retries = 0;
    *status = ACTIVE;

    // ==========================================
    // INITIAL ANALYTIC METRIC EVALUATION
    // ==========================================
    interpolation_kernel_{spacetime}(
      &commondata, f, metric, NULL, chunk_size, stream_idx
    );
    double observer_tetrad[4][4];
    set_initial_conditions_kernel(
      &commondata, num_rays, &all_photons, metric, observer_tetrad
    );
{normalized_momentum_conversion}
    for (int history_step = 0; history_step < {maximum_degree}; ++history_step) {{
      for (int component = 0; component < 9; ++component)
        f_event_history[history_step * 9 + component] = f[component];
      integration_param_event_history[history_step] = *integration_param;
    }} // END LOOP: initialize accepted-state history
{initial_normalization_check}

    printf("Initial State:\n");
    printf("  Pos (%.4f, %.4f, %.4f)\n", f[1], f[2], f[3]);
{initial_state_report}

    // ==========================================
    // PRE-INTEGRATION DIAGNOSTICS
    // ==========================================
{conserved_quantity_initialization}

    fp = fopen("trajectory.txt", "w");
    if (fp == NULL) {{
      fprintf(stderr, "Error: could not open trajectory.txt for writing.\n");
      exit_status = EXIT_FAILURE;
      goto cleanup;
    }} // END IF: trajectory output unavailable
    fprintf(fp, "{trajectory_header}");
{trial_debug_open}
{stage_debug_open}

    // ==========================================
    // MODULAR SPLIT-PIPELINE INTEGRATION LOOP
    // ==========================================
    int steps = 0;
{trial_debug_attempt_declaration}

    const long int max_accepted_steps = 200000;
    while (steps < max_accepted_steps) {{
      for (int i = 0; i < 9; i++) {{
        f_base[i] = f[i];
        f_temp[i] = f[i];
      }} // END LOOP: copy current state
{trial_debug_trial_metadata}

      for (int stage = 1; stage <= 6; stage++) {{
        // Evaluate the analytic metric and connection at the stage state.
        interpolation_kernel_{spacetime}(
          &commondata, f_temp, metric, connection, chunk_size, stream_idx
        );
        calculate_ode_rhs_kernel(
          f_temp, metric, connection, {rhs_integration_arguments}
          k_bundle, stage, chunk_size, stream_idx
        );
{stage_debug_record}

        if (stage < 6)
          rkf45_stage_update(
            f, k_bundle, h, stage, chunk_size, f_temp, stream_idx
          );
      }} // END LOOP: for stage over RKF45 stages

      rkf45_finalize_and_control(
        &commondata,
        f,
        f_base,
        k_bundle,
        h,
        status,
        integration_param,
        rejection_retries,
        {trial_debug_call_argument}chunk_size,
        stream_idx
      );
{trial_debug_attempt_increment}
{trial_debug_record}

      if (*rejection_retries == 0) {{
        interpolation_kernel_{spacetime}(
          &commondata, f, metric, NULL, chunk_size, stream_idx
        );
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
            steps + 1 >= max_accepted_steps,
            event_chunk_index,
            chunk_size,
            stream_idx);
        }} // END IF: a physical plane is enabled

        fprintf(
          fp,
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
          f[8]
        );
        steps++;

        if (event_planes_enabled && *status != ACTIVE) {{
          if (*status == STOP_CONDITION_TERMINAL_PLANE) {{
            printf("Photon crossed the configured terminal plane.\n");
          }} else if (*status == STOP_CONDITION_COORD_RADIUS_EXCEEDED) {{
            printf("Particle escaped to r > %.2f.\n", commondata.r_escape);
          }} else if (*status == STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED) {{
            printf("Evolution measure exceeded numerical limit.\n");
          }} else if (*status == FAILURE_PLANE_INTERPOLATION_HISTORY ||
                     *status == FAILURE_GENERIC) {{
            fprintf(stderr, "Plane crossing failed with status %d.\n", *status);
            exit_status = EXIT_FAILURE;
          }} // END ELSE IF: event manager established a stop status
          break;
        }} // END IF: event manager stopped the photon

        if (log_energy_measure > commondata.evolution_measure_max) {{
          *status = STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED;
          printf("Evolution measure exceeded numerical limit.\n");
          break;
        }} // END IF: evolution measure exceeded limit
      }} // END IF: accepted step written

      const double r2 = f[1] * f[1] + f[2] * f[2] + f[3] * f[3];
      if (r2 > commondata.r_escape * commondata.r_escape) {{
        printf("Particle escaped to r > %.2f.\n", commondata.r_escape);
        break;
      }} // END IF: particle crossed the escape radius

      if (*status == FAILURE_RKF45_REJECTION_LIMIT) {{
        printf(
          "Integration terminated with status: %s (%d)\n",
          status_names[*status],
          *status
        );
        break;
      }} // END IF: RKF45 rejection limit was reached

      if (steps >= max_accepted_steps)
        printf(
          "Termination: reached maximum accepted step count (%ld).\n",
          max_accepted_steps
        );
    }} // END WHILE: integrate the photon geodesic

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

{final_parameter_report}

    // ==========================================
    // POST-INTEGRATION DIAGNOSTICS
    // ==========================================
{conserved_quantity_finalization}

    normalization_constraint_t norm_final;
    interpolation_kernel_{spacetime}(
      &commondata, f, metric, NULL, chunk_size, stream_idx
    );
    {normalization_kernel_name}(
      f, metric, &norm_final, chunk_size, stream_idx
    );

    printf("\nFinal normalization constraint error = %.4e\n",
           {normalization_error_expression});
{conservation_error_report}
{conserved_quantity_q_report}

    cleanup:
{trial_debug_cleanup}
{stage_debug_cleanup}
    if (fp != NULL)
      fclose(fp);
    if (plane_crossings_file != NULL)
      fclose(plane_crossings_file);
    BHAH_FREE(f);
    BHAH_FREE(f_base);
    BHAH_FREE(f_temp);
    BHAH_FREE(metric);
    BHAH_FREE(connection);
    BHAH_FREE(k_bundle);
    BHAH_FREE(integration_param);
    BHAH_FREE(h);
    BHAH_FREE(rejection_retries);
    BHAH_FREE(status);
    BHAH_FREE(log_energy_bundle);

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

    par.adjust_CodeParam_default("rkf45_absolute_error_tolerance", 1e-17)
    par.adjust_CodeParam_default("rkf45_error_tolerance", 1e-17)
    par.adjust_CodeParam_default("rkf45_h_max", 10.0)
    par.adjust_CodeParam_default("rkf45_h_min", 1e-20)
    par.adjust_CodeParam_default("rkf45_max_retries", 15)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

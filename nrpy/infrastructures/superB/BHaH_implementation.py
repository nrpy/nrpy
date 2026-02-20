"""
Register C++ function for finding horizons with superB BHaHAHA.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.jacobians as jac
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    griddata_commondata,
)


def string_for_step7e_to_g_main_loop_for_each_horizon(
    single_horizon: bool = False,
    enable_timing: bool = True,
) -> str:
    r"""
    Generate the C string for STEP 7 e to g: main loop for each horizon.

    :param single_horizon: If True, emit code for a single horizon (`which_horizon`). Otherwise, emit the standard for-loop.
    :param enable_timing: If False, skip emitting `gettimeofday` lines.
    :return: Raw C string.
    """
    outstring = r"""
    if (commondata->bah_verbosity_level >= 2) {"""
    if enable_timing:
        outstring += r"""
        struct timeval temp_time;
        gettimeofday(&temp_time, NULL);
        printf("%.6f s : FINISH interpolate H%d (%.3fs)\n", timeval_to_seconds(start_time_total, temp_time), h,
               timeval_to_seconds(iter_time_tracker, temp_time));
        gettimeofday(&iter_time_tracker, NULL); // Reset for solver timing"""
    else:
        outstring += r"""
        printf("FINISH interpolate\n");"""
    outstring += r"""
    } // END IF: verbosity for finish interpolate

    // STEP 7.e: Call `bah_find_horizon` solver if metric data available.
    if (current_horizon_params->input_metric_data && current_horizon_params->Nr_external_input > 0) {
      bhahaha_diagnostics_struct diags_local; // Local diagnostics struct per thread/horizon

      const int rc = bah_find_horizon(current_horizon_params, &diags_local);

      if (rc == BHAHAHA_SUCCESS) {
        // STEP 7.e.i: On success, update state and output diagnostics.
        current_horizon_params->use_fixed_radius_guess_on_full_sphere = 0;

        bah_diagnostics_file_output(&diags_local, current_horizon_params, commondata->bah_max_num_horizons, x_guess[h], y_guess[h], z_guess[h], ".");
      } else { // Horizon find failed
        // STEP 7.e.ii: On failure, log warning and set flags for next time.
        fprintf(stderr,
                "WARNING: it=%d t=%.2e Horizon %d (idx %d) find FAILED (error code=%d): %s. Forcing fixed-radius, full-sphere guess next time.\n",
                commondata->nn, commondata->time, h + 1, h, rc, bah_error_message((bhahaha_error_codes)rc));
        current_horizon_params->use_fixed_radius_guess_on_full_sphere = 1;
        current_horizon_params->t_m1 = -1.0; // Mark as not found for extrapolation logic.
      } // END ELSE: solver failed
    } else { // Skipped solve (due to no input_metric_data or Nr_external_input == 0)
      // STEP 7.f: If solve was skipped, set flags for next time.
      if (commondata->bah_verbosity_level > 0 && current_horizon_params->Nr_external_input == 0) {
        printf("INFO: Horizon %d (idx %d) skipped solve due to Nr_external_input=0.\n", h + 1, h);
      } // END IF: verbosity for skipped solve
      current_horizon_params->use_fixed_radius_guess_on_full_sphere = 1;
      current_horizon_params->t_m1 = -1.0;
    } // END ELSE: skipped solve

    // STEP 7.g: Free the allocated `input_metric_data` buffer.
    BHAH_FREE(current_horizon_params->input_metric_data);

    if (commondata->bah_verbosity_level >= 1) {"""
    if enable_timing:
        outstring += r"""
      struct timeval temp_time;
      gettimeofday(&temp_time, NULL);
      printf("%.6f s : FINISH find H%d (%.3fs for solve section)\n", timeval_to_seconds(start_time_total, temp_time), h,
             timeval_to_seconds(iter_time_tracker, temp_time));"""
    else:
        outstring += r"""
      printf("FINISH find\n");"""
    outstring += r"""
    } // END IF: verbosity for finish find horizon"""
    if not single_horizon:
        outstring += r"""
  } // END LOOP: for h (main find loop)
"""
    return outstring


def string_for_step7a_to_c_main_loop_for_each_horizon(
    single_horizon: bool = False,
    allocate_radii_interp: bool = False,
    enable_timing: bool = True,
    enable_BBH_mode: bool = True,
) -> str:
    """
    Generate the C string for STEP 7 a to c: main loop for each horizon.

    :param single_horizon: If True, emit code for a single horizon (`which_horizon`). Otherwise, emit the standard for-loop.
    :param allocate_radii_interp: If True, allocate `radii_interp` on the heap
                                   with malloc; otherwise, use stack array.
    :param enable_timing: If False, skip emitting `gettimeofday` lines.
    :param enable_BBH_mode: **Only used when `single_horizon` is True.**
                            If True, emit the BBH activity guard
                            (`if (!commondata->bah_BBH_mode_horizon_active[h]) return;`)
                            for the selected horizon. If False, **omit** that guard.
                            (No effect when `single_horizon` is False.)
    :return: Raw C string.
    """
    if single_horizon:
        outstring = r"""
      // STEP 7: Main computation for one horizon.
        int h = which_horizon;
        """
        if enable_BBH_mode:
            outstring += r"""
      // STEP 7.a: Skip if horizon is not active.
        if (!commondata->bah_BBH_mode_horizon_active[h])
          return; // END IF: not active, return
"""
    else:
        outstring = r"""
  // STEP 7: Main loop for each horizon.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
  // STEP 7.a: Skip if horizon is not active.
    if (!commondata->bah_BBH_mode_horizon_active[h])
      continue; // END IF: not active, continue loop
"""
    if allocate_radii_interp:
        radii_interp_decl = "*radii_interp = (REAL*)malloc(sizeof(REAL) * commondata->bah_Nr_interp_max);"
        radii_interp_arg = "*radii_interp"
    else:
        radii_interp_decl = "REAL radii_interp[commondata->bah_Nr_interp_max];"
        radii_interp_arg = "radii_interp"

    if enable_timing:
        outstring += r"""
    gettimeofday(&iter_time_tracker, NULL); // Per-horizon timer start
"""

    outstring += rf"""
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];

    // STEP 7.b: Set up the radial grid for interpolation.
    {radii_interp_decl}                           // Max capacity from commondata (scalar)
    bah_radial_grid_cell_centered_set_up(commondata->bah_Nr_interp_max,         // Max capacity of radii_interp array
                                         commondata->bah_max_search_radius[h],  // Overall max search radius for this horizon
                                         r_min_guess[h], r_max_guess[h], // Current guess for search range from extrapolation
                                         &current_horizon_params->Nr_external_input, // Output: Actual Nr used for this horizon
                                         &current_horizon_params->r_min_external_input, // Output: Actual r_min used
                                         &current_horizon_params->dr_external_input, // Output: Actual dr used
                                         {radii_interp_arg});                         // Output: populated radii for this horizon

    // STEP 7.c: Allocate buffer for interpolated ADM metric data.
    const int Ntheta_interp_eff = current_horizon_params->Ntheta_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
    const int Nphi_interp_eff = current_horizon_params->Nphi_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
    const size_t npts_metric_adm =
        (size_t)BHAHAHA_NUM_METRIC_COMPONENTS * current_horizon_params->Nr_external_input * Ntheta_interp_eff * Nphi_interp_eff;

    if (current_horizon_params->Nr_external_input > 0) {{
      current_horizon_params->input_metric_data = (REAL *)malloc(npts_metric_adm * sizeof(REAL));
      if (!current_horizon_params->input_metric_data) {{
        fprintf(stderr, "ERROR: malloc() failed for input_metric_data for H%d (%zu REALs).\n", h, npts_metric_adm);
        exit(EXIT_FAILURE);
      }} // END IF: malloc failed for input_metric_data
    }} else {{ // Nr_external_input <= 0
      current_horizon_params->input_metric_data = NULL;
    }} // END ELSE: Nr_external_input <= 0

    if (commondata->bah_verbosity_level >= 2) {{
      struct timeval temp_time;"""
    if enable_timing:
        outstring += r"""
      gettimeofday(&temp_time, NULL);
      printf("%.6f s : START interpolate H%d\n", timeval_to_seconds(start_time_total, temp_time), h);"""
    else:
        outstring += r"""
      printf("START interpolate\n");"""
    outstring += r"""
    } // END IF: verbosity for start interpolate
"""
    return outstring


def string_for_step6_apply_robustness_improv_and_extrap_horizon_guesses(
    enable_BBH_mode: bool = True,
) -> str:
    r"""
    Generate the C string for STEP 6: apply robustness improvements and extrapolate horizon guesses.


    :param enable_BBH_mode: If True, include the extra robustness block
                            for the BBH common horizon (STEP 6.a.iii).
                            If False, skip it.
    :return: Raw C string.
    """
    outstring = r"""
  // STEP 6: Apply robustness improvements and extrapolate horizon guesses.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];

    // STEP 6.a: For each horizon, if not found reliably, adjust parameters.
    if (current_horizon_params->t_m1 == -1.0 || current_horizon_params->t_m2 == -1.0 || current_horizon_params->t_m3 == -1.0) {
      // STEP 6.a.i: Reduce CFL.
      current_horizon_params->cfl_factor *= 0.9;
      // STEP 6.a.ii: Force fixed radius guess.
      current_horizon_params->use_fixed_radius_guess_on_full_sphere = 1; // Force fixed radius if not found reliably.
"""
    if enable_BBH_mode:
        outstring += r"""
      // STEP 6.a.iii: If BBH common horizon, adjust resolution/iterations.
      if (commondata->bah_enable_BBH_mode && h == commondata->bah_BBH_mode_common_horizon_idx &&
          current_horizon_params->num_resolutions_multigrid >= 2) {
        current_horizon_params->Ntheta_array_multigrid[0] = commondata->bah_Ntheta_array_multigrid[1]; // Use 2nd level res for 1st.
        current_horizon_params->Nphi_array_multigrid[0] = commondata->bah_Nphi_array_multigrid[1];
        current_horizon_params->max_iterations *= 2;
      } // END IF: BBH mode common horizon robustness for resolution/iterations"""
    outstring += r"""
    } // END IF: robustness improvements for not-yet-found horizons

    // STEP 6.b: Extrapolate current guesses.
    bah_xyz_center_r_minmax(current_horizon_params, &x_guess[h], &y_guess[h], &z_guess[h], &r_min_guess[h], &r_max_guess[h]);
    if (current_horizon_params->use_fixed_radius_guess_on_full_sphere) { // This flag is key.
      r_min_guess[h] = 0.0;
      r_max_guess[h] = commondata->bah_max_search_radius[h]; // Use per-horizon parameter from commondata.
    } // END IF: use_fixed_radius_guess_on_full_sphere for extrapolation override
  } // END LOOP: for h (Robustness & Extrapolation)
"""
    return outstring


def string_for_step4_populate_iteration_metadata() -> str:
    r"""
    Generate the C string for STEP 4: populate current iteration metadata.

    :return: Raw C string.
    """
    outstring = r"""
  // STEP 4: Populate current iteration metadata.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];
    current_horizon_params->which_horizon = h + 1; // 1-based for BHaHAHA library
    current_horizon_params->num_horizons = commondata->bah_max_num_horizons;
    current_horizon_params->iteration_external_input = commondata->nn;
    current_horizon_params->time_external_input = commondata->time;
  } // END LOOP: for h (populating current iteration metadata)
"""
    return outstring


def string_for_step3_initialize_bhahaha_data_structs_and_solver_params(
    enable_BBH_mode: bool = True,
) -> str:
    r"""
    Generate the C string for STEP 3: initialize BHaHAHA data structures and solver parameters.

    :param enable_BBH_mode: If True, include BBH mode initialization logic.
                            If False, skip it and always activate all horizons.
    :return: Raw C string.
    """
    outstring = r"""
  // STEP 3: Initialize BHaHAHA data structures and solver parameters.
  if (commondata->nn == 0) {
    // STEP 3.a: If it's the first call (nn=0):
    for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
      // STEP 3.a.i: Initialize persistent fields.
      bhahaha_params_and_data_struct *current_horizon_data = &commondata->bhahaha_params_and_data[h];
      current_horizon_data->x_center_m1 = commondata->bah_initial_grid_x_center[h];
      current_horizon_data->y_center_m1 = commondata->bah_initial_grid_y_center[h];
      current_horizon_data->z_center_m1 = commondata->bah_initial_grid_z_center[h];
      current_horizon_data->x_center_m2 = current_horizon_data->x_center_m3 = current_horizon_data->x_center_m1;
      current_horizon_data->y_center_m2 = current_horizon_data->y_center_m3 = current_horizon_data->y_center_m1;
      current_horizon_data->z_center_m2 = current_horizon_data->z_center_m3 = current_horizon_data->z_center_m1;
      current_horizon_data->t_m1 = current_horizon_data->t_m2 = current_horizon_data->t_m3 = -1.0;
      current_horizon_data->r_min_m1 = current_horizon_data->r_min_m2 = current_horizon_data->r_min_m3 = 0.0;
      current_horizon_data->r_max_m1 = current_horizon_data->r_max_m2 = current_horizon_data->r_max_m3 = commondata->bah_max_search_radius[h];
      current_horizon_data->use_fixed_radius_guess_on_full_sphere = 1; // Default to fixed radius guess

      // STEP 3.a.ii: Call `initialize_bhahaha_solver_params_and_shapes` for allocation and non-persistent params.
      initialize_bhahaha_solver_params_and_shapes(commondata, h, true); // true for shape alloc

    } // END LOOP: for h (iteration 0 persistent data initialization)
"""
    if enable_BBH_mode:
        outstring += r"""
    // STEP 3.a.iii: Initialize `commondata->bah_BBH_mode_horizon_active`.
    if (commondata->bah_enable_BBH_mode) {
      if (commondata->bah_max_num_horizons != 3) {
        fprintf(stderr,
                "ERROR: bah_enable_BBH_mode requires bah_max_num_horizons==3, to account for common horizon, plus two individual horizons.\n");
        exit(EXIT_FAILURE);
      } // END IF: incorrect num_horizons for BBH mode
      const int bh1 = commondata->bah_BBH_mode_inspiral_BH_idxs[0];
      const int bh2 = commondata->bah_BBH_mode_inspiral_BH_idxs[1];
      const int com = commondata->bah_BBH_mode_common_horizon_idx;
      if (bh1 < 0 || bh1 >= commondata->bah_max_num_horizons || bh2 < 0 || bh2 >= commondata->bah_max_num_horizons || com < 0 ||
          com >= commondata->bah_max_num_horizons || bh1 == bh2 || bh1 == com || bh2 == com) {
        fprintf(stderr, "ERROR: Invalid or overlapping BBH horizon indices: bh1=%d, bh2=%d, com=%d. Max idx is %d.\n", bh1, bh2, com,
                commondata->bah_max_num_horizons - 1);
        exit(EXIT_FAILURE);
      } // END IF: invalid BBH indices
      commondata->bah_BBH_mode_horizon_active[bh1] = 1;
      commondata->bah_BBH_mode_horizon_active[bh2] = 1;
      commondata->bah_BBH_mode_horizon_active[com] = 0;
    } else { // Not BBH mode
      for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
        commondata->bah_BBH_mode_horizon_active[h] = 1; // Activate all configured horizons
      } // END LOOP: for h (activating all horizons if not BBH mode)
    } // END ELSE: not enable_BBH_mode (iteration 0 horizon activity)"""
    else:
        outstring += r"""
    // STEP 3.a.iii: Initialize `commondata->bah_BBH_mode_horizon_active` (non-BBH mode only).
    for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
      commondata->bah_BBH_mode_horizon_active[h] = 1; // Activate all configured horizons
    } // END LOOP: for h (activating all horizons if BBH mode disabled in codegen)
"""
    outstring += r"""
  } else { // Not the first call (nn != 0)
    // STEP 3.b: If not the first call:
    for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
      // STEP 3.b.i: Call `initialize_bhahaha_solver_params_and_shapes` for non-persistent params.
      initialize_bhahaha_solver_params_and_shapes(commondata, h, false); // false: don't re-alloc shapes
    } // END LOOP: for h (non-first call initialization of solver params)
  } // END ELSE: not commondata->nn == 0
"""
    return outstring


def string_for_step2_validate_multigrid_inputs() -> str:
    r"""
    Generate the C string for STEP 2: validate multigrid resolution inputs.

    :return: Raw C string.
    """
    outstring = r"""
// STEP 2: Validate multigrid resolution inputs.
  check_multigrid_resolution_inputs(commondata); // Validates Ntheta, Nphi arrays in commondata
"""
    return outstring


def string_for_static_func_initialize_bhahaha_solver_params_and_shapes() -> str:
    r"""
    Generate the C string for initialize_bhahaha_solver_params_and_shapes static function.

    :return: Raw C string.
    """
    outstring = r"""
/**
 * Initializes non-persistent BHaHAHA solver parameters for a specific horizon and
 * allocates memory for horizon shape history arrays if it's the first call.
 *
 * The function performs the following steps:
 * 1. Gets a pointer to the current horizon's `bhahaha_params_and_data_struct`.
 * 2. If `is_first_call_for_shapes` is true:
 *    a. Calculates max Ntheta and Nphi for this horizon.
 *    b. Allocates memory for `prev_horizon_m1`, `prev_horizon_m2`, and `prev_horizon_m3`.
 *    c. Exits if memory allocation fails.
 * 3. Copies non-persistent solver parameters (verbosity, max_iterations, CFL, tolerances,
 *    M_scale, eta_damping, KO_strength, multigrid resolutions) from `commondata`
 *    to the `current_horizon_params`.
 * 4. Sets `current_horizon_params->input_metric_data` to NULL.
 *
 * @param commondata - Pointer to `commondata_struct` containing global settings and the
 *                     `bhahaha_params_and_data` array.
 * @param h - Index of the horizon to initialize.
 * @param is_first_call_for_shapes - True if shape history arrays need allocation, false otherwise.
 * @return - None (`void`).
 */
static void initialize_bhahaha_solver_params_and_shapes(commondata_struct *restrict commondata, int h, bool is_first_call_for_shapes) {
  // STEP 1: Get a pointer to the current horizon's `bhahaha_params_and_data_struct`.
  bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];

  // STEP 2: If `is_first_call_for_shapes` is true:
  if (is_first_call_for_shapes) {
    // STEP 2.a: Calculates max Ntheta and Nphi for this horizon.
    const int Ntheta_max_this_h = commondata->bah_Ntheta_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
    const int Nphi_max_this_h = commondata->bah_Nphi_array_multigrid[commondata->bah_num_resolutions_multigrid - 1];
    // STEP 2.b: Allocates memory for `prev_horizon_m1`, `prev_horizon_m2`, and `prev_horizon_m3`.
    current_horizon_params->prev_horizon_m1 = (REAL *)malloc(Ntheta_max_this_h * Nphi_max_this_h * sizeof(REAL));
    current_horizon_params->prev_horizon_m2 = (REAL *)malloc(Ntheta_max_this_h * Nphi_max_this_h * sizeof(REAL));
    current_horizon_params->prev_horizon_m3 = (REAL *)malloc(Ntheta_max_this_h * Nphi_max_this_h * sizeof(REAL));
    // STEP 2.c: Exits if memory allocation fails.
    if (!current_horizon_params->prev_horizon_m1 || !current_horizon_params->prev_horizon_m2 || !current_horizon_params->prev_horizon_m3) {
      fprintf(stderr, "ERROR: Memory allocation failed for prev_horizon shape data for horizon %d.\n", h);
      exit(EXIT_FAILURE);
    } // END IF: memory allocation failed for shape data
  } // END IF: is_first_call_for_shapes

  // STEP 3: Copies non-persistent solver parameters from `commondata` to `current_horizon_params`.
  current_horizon_params->enable_eta_varying_alg_for_precision_common_horizon = commondata->bah_enable_eta_varying_alg_for_precision_common_horizon;
  current_horizon_params->verbosity_level = commondata->bah_verbosity_level;
  current_horizon_params->max_iterations = commondata->bah_max_iterations;
  current_horizon_params->cfl_factor = commondata->bah_cfl_factor[h];
  current_horizon_params->Theta_L2_times_M_tolerance = commondata->bah_Theta_L2_times_M_tolerance[h];
  current_horizon_params->Theta_Linf_times_M_tolerance = commondata->bah_Theta_Linf_times_M_tolerance[h];
  current_horizon_params->M_scale = commondata->bah_M_scale[h];
  current_horizon_params->eta_damping_times_M = commondata->bah_eta_damping_times_M[h];
  current_horizon_params->KO_strength = commondata->bah_KO_strength[h];
  current_horizon_params->num_resolutions_multigrid = commondata->bah_num_resolutions_multigrid;
  memcpy(current_horizon_params->Ntheta_array_multigrid, commondata->bah_Ntheta_array_multigrid,
         sizeof(int) * commondata->bah_num_resolutions_multigrid);
  memcpy(current_horizon_params->Nphi_array_multigrid, commondata->bah_Nphi_array_multigrid, sizeof(int) * commondata->bah_num_resolutions_multigrid);

  // STEP 4: Sets `current_horizon_params->input_metric_data` to NULL.
  current_horizon_params->input_metric_data = NULL;
} // END FUNCTION: initialize_bhahaha_solver_params_and_shapes
"""
    return outstring


def string_for_static_func_check_multigrid_resolution_inputs() -> str:
    r"""
    Generate the C string for static function: check_multigrid_resolution_inputs.

    :return: Raw C string.
    """
    outstring = r"""
/**
 * Validates BHaHAHA multigrid resolution inputs.
 * Checks if `bah_num_resolutions_multigrid` is positive and if all
 * `bah_Ntheta_array_multigrid` and `bah_Nphi_array_multigrid` entries are positive.
 * If validation fails, prints an error and exits.
 *
 * @param commondata - Pointer to `commondata_struct` containing multigrid settings.
 * @return - None (`void`).
 */
static void check_multigrid_resolution_inputs(const commondata_struct *restrict commondata) {
  int trigger_error = 0;
  if (commondata->bah_num_resolutions_multigrid <= 0) {
    trigger_error = 1;
  } else { // num_resolutions_multigrid > 0
    for (int res = 0; res < commondata->bah_num_resolutions_multigrid; res++) {
      if (commondata->bah_Ntheta_array_multigrid[res] <= 0 || commondata->bah_Nphi_array_multigrid[res] <= 0) {
        trigger_error = 1;
        break;
      } // END IF: invalid resolution
    } // END LOOP: for res
  } // END ELSE: num_resolutions_multigrid > 0
  if (trigger_error) {
    fprintf(stderr, "ERROR: BHaHAHA multigrid resolutions are unset or invalid. Please specify "
                    "non-zero values, e.g.:\n"
                    "   bah_num_resolutions_multigrid = 3\n"
                    "   bah_Ntheta_array_multigrid   = [8, 16, 32]\n"
                    "   bah_Nphi_array_multigrid     = [16, 32, 64]\n");
    exit(EXIT_FAILURE);
  } // END IF: trigger_error
} // END FUNCTION: check_multigrid_resolution_inputs
"""
    return outstring


def string_for_prefunc_enums_and_interp_indices(
    add_bhahaha_gf_interp_ind_to_bhah_defines: bool = False,
) -> str:
    r"""
    Generate the C string for prefunc enums and interpolator index arrays.

    :param add_bhahaha_gf_interp_ind_to_bhah_defines: add to BHaH_defines.h instead of defining here.
    :return: Raw C string.
    """
    outstring = r"""
// Enum for indexing the final ADM metric components BHaHAHA expects.
enum FINAL_ADM_METRIC_INDICES {
  FINAL_INTERP_GAMMADDXXGF,
  FINAL_INTERP_GAMMADDXYGF,
  FINAL_INTERP_GAMMADDXZGF,
  FINAL_INTERP_GAMMADDYYGF,
  FINAL_INTERP_GAMMADDYZGF,
  FINAL_INTERP_GAMMADDZZGF,
  FINAL_INTERP_KDDXXGF,
  FINAL_INTERP_KDDXYGF,
  FINAL_INTERP_KDDXZGF,
  FINAL_INTERP_KDDYYGF,
  FINAL_INTERP_KDDYZGF,
  FINAL_INTERP_KDDZZGF,
  BHAHAHA_NUM_METRIC_COMPONENTS
}; // END ENUM: FINAL_ADM_METRIC_INDICES
"""
    enum_INTERP_BSSN_GF_INDICES_parts = """INTERP_ADD00GF_IDX,
  INTERP_ADD01GF_IDX,
  INTERP_ADD02GF_IDX,
  INTERP_ADD11GF_IDX,
  INTERP_ADD12GF_IDX,
  INTERP_ADD22GF_IDX,
  INTERP_CFGF_IDX,
  INTERP_HDD00GF_IDX,
  INTERP_HDD01GF_IDX,
  INTERP_HDD02GF_IDX,
  INTERP_HDD11GF_IDX,
  INTERP_HDD12GF_IDX,
  INTERP_HDD22GF_IDX,
  INTERP_TRKGF_IDX,"""

    if not add_bhahaha_gf_interp_ind_to_bhah_defines:
        outstring += f"""
// Enum for indexing interpolated BSSN gridfunctions. (NRPy-specific)
enum INTERP_BSSN_GF_INDICES {{
  {enum_INTERP_BSSN_GF_INDICES_parts}
  BHAHAHA_NUM_INTERP_GFS
}}; // END ENUM: INTERP_BSSN_GF_INDICES
"""
    else:
        outstring += f"""
// Enum for indexing interpolated BSSN gridfunctions. (NRPy-specific)
enum INTERP_BSSN_GF_INDICES {{
  {enum_INTERP_BSSN_GF_INDICES_parts}
}}; // END ENUM: INTERP_BSSN_GF_INDICES
"""

    bhahaha_gf_interp_indices_string = r"""// BSSN gridfunctions input into the interpolator. Must be in the same order as the enum list below.
const int bhahaha_gf_interp_indices[BHAHAHA_NUM_INTERP_GFS] = {
    ADD00GF, ADD01GF, ADD02GF, ADD11GF, ADD12GF, ADD22GF, // Traceless, rescaled extrinsic curvature components.
    CFGF,                                                 // Conformal factor.
    HDD00GF, HDD01GF, HDD02GF, HDD11GF, HDD12GF, HDD22GF, // Rescaled conformal 3-metric components.
    TRKGF                                                 // Trace of extrinsic curvature.
};"""

    if not add_bhahaha_gf_interp_ind_to_bhah_defines:
        outstring += "\n" + bhahaha_gf_interp_indices_string + "\n"
    else:
        BHaH_defines_h.register_BHaH_defines(
            __name__,
            "\n#define INTERP_REQUEST_BHAHAHA 0\n#define INTERP_REQUEST_PSI4 1\n\n#define BHAHAHA_NUM_INTERP_GFS 14\n"
            + bhahaha_gf_interp_indices_string,
        )

    return outstring


def string_for_func_free_bhahaha_horizon_shape_data_all_horizons() -> str:
    r"""
    Generate the C string for function: free_bhahaha_horizon_shape_data_all_horizons.

    :return: Raw C string.
    """
    outstring = r"""
/**
 * Frees memory allocated for horizon shape history arrays (`prev_horizon_m1/m2/m3`)
 * for all horizons.
 *
 * @param commondata - Pointer to `commondata_struct` containing the BHaHAHA data.
 * @return - None (`void`).
 */
void free_bhahaha_horizon_shape_data_all_horizons(commondata_struct *restrict commondata) {
  for (int h = 0; h < commondata->bah_max_num_horizons; ++h) {
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];
    if (current_horizon_params != NULL) {
      BHAH_FREE(current_horizon_params->prev_horizon_m1);
      BHAH_FREE(current_horizon_params->prev_horizon_m2);
      BHAH_FREE(current_horizon_params->prev_horizon_m3);
    } // END IF: current_horizon_params != NULL
  } // END LOOP: for h
} // END FUNCTION: free_bhahaha_horizon_shape_data_all_horizons
"""
    return outstring


def generate_bssn_to_adm_codegen(CoordSystem: str) -> str:
    r"""
    Generate the C code string for BSSN→ADM Cartesian transformation.

    Sets up rfm tensors, computes $\bar{\gamma}_{ij}$ and $\bar{A}_{ij}$
    in Cartesian basis, applies conformal factor choice (phi/chi/W),
    and constructs ADM $\gamma_{ij}$ and $K_{ij}$. Returns C code
    generated by NRPy’s c_codegen.

    :param CoordSystem: Coordinate system key in NRPy reference_metric.
    :return: Raw C code string implementing the transformation.
    :raises ValueError: If ``EvolvedConformalFactor_cf`` is unsupported.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    rfm_aDD = ixp.declarerank2("rfm_aDD", symmetry="sym01")
    rfm_hDD = ixp.declarerank2("rfm_hDD", symmetry="sym01")
    rfm_gammabarDD = ixp.zerorank2()
    rfm_AbarDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            rfm_gammabarDD[i][j] = rfm_hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            rfm_AbarDD[i][j] = rfm_aDD[i][j] * rfm.ReDD[i][j]
    Cart_gammabarDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        CoordSystem, rfm_gammabarDD
    )
    Cart_AbarDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        CoordSystem, rfm_AbarDD
    )
    exp4phi = sp.sympify(0)
    cf = sp.Symbol("cf", real=True)
    EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    if EvolvedConformalFactor_cf == "phi":
        exp4phi = sp.exp(4 * cf)
    elif EvolvedConformalFactor_cf == "chi":
        exp4phi = 1 / cf
    elif EvolvedConformalFactor_cf == "W":
        exp4phi = 1 / cf**2
    else:
        raise ValueError(
            f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
        )
    expr_list: List[sp.Expr] = []
    name_list: List[str] = []
    Cart_gammaDD = ixp.zerorank2()
    labels = ["X", "Y", "Z"]
    for i in range(3):
        for j in range(i, 3):
            Cart_gammaDD[i][j] = exp4phi * Cart_gammabarDD[i][j]
            expr_list += [Cart_gammaDD[i][j]]
            name_list += [
                f"input_metric_data_target_array[FINAL_INTERP_GAMMADD{labels[i]}{labels[j]}GF * offset + idx3]"
            ]
    Cart_KDD = ixp.zerorank2()
    trK_sym = sp.Symbol("trK", real=True)
    for i in range(3):
        for j in range(i, 3):
            Cart_KDD[i][j] = (
                exp4phi * Cart_AbarDD[i][j]
                + sp.Rational(1, 3) * Cart_gammaDD[i][j] * trK_sym
            )
            expr_list += [Cart_KDD[i][j]]
            name_list += [
                f"input_metric_data_target_array[FINAL_INTERP_KDD{labels[i]}{labels[j]}GF * offset + idx3]"
            ]
    codegen_str = ccg.c_codegen(
        expr_list, name_list, include_braces=True, verbose=False
    )
    return codegen_str


def register_bhahaha_commondata_and_params(max_horizons: int) -> None:
    """
    Register BHaHAHA commondata arrays and all related CodeParameters.

    :param max_horizons: Maximum number of horizons to support.
    """
    griddata_commondata.register_griddata_commondata(
        __name__,
        f"bhahaha_params_and_data_struct bhahaha_params_and_data[{max_horizons}]",
        "BHaHAHA parameters and data, including previous horizon data",
        is_commondata=True,
    )
    griddata_commondata.register_griddata_commondata(
        __name__,
        f"bhahaha_diagnostics_struct bhahaha_diagnostics[{max_horizons}]",
        "BHaHAHA diagnostics",
        is_commondata=True,
    )

    # Register commondata parameters
    # Dictionary for REAL type parameters
    real_params = {
        "bah_max_search_radius": ("Maximum radius to search for horizons", 1.5),
        "bah_cfl_factor": ("CFL factor for timestep. 0.98 recommended.", 0.98),
        "bah_Theta_Linf_times_M_tolerance": (
            "Convergence criterion: Linf norm of residual times M_scale, 1e-2 recommended",
            1e-2,
        ),
        "bah_Theta_L2_times_M_tolerance": (
            "Convergence criterion: L2 norm of residual times M_scale, 2e-5 recommended",
            2e-5,
        ),
        "bah_eta_damping_times_M": (
            "Exponential damping parameter, 3.5 recommended",
            7.0,
        ),
        "bah_M_scale": ("Mass scale for horizon", 1.0),
        "bah_KO_strength": (
            "Kreiss-Oliger dissipation parameter, 0.0 strongly recommended",
            0.0,
        ),
        "bah_initial_grid_x_center": ("Initial x-center coordinates", 0.0),
        "bah_initial_grid_y_center": ("Initial y-center coordinates", 0.0),
        "bah_initial_grid_z_center": ("Initial z-center coordinates", 0.0),
    }
    for key, value in real_params.items():
        _ = par.register_CodeParameter(
            cparam_type=f"REAL[{max_horizons}]",
            module="BHaHAHA",
            name=key,
            defaultvalue=value[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=value[0],
        )

    # Dictionary for INT type parameters
    int_params = {
        "bah_Nr_interp_max": (
            "Maximum number of radial grid points used in 3D interpolation",
            48,
        ),
        "bah_verbosity_level": ("Verbosity level", 2),
        "bah_max_iterations": (
            "Maximum iterations before giving up. 10000 recommended",
            5000,
        ),
        "bah_max_num_horizons": (
            "Maximum number of horizons to find",
            max_horizons,
        ),
        "bah_num_resolutions_multigrid": (
            "Number of resolutions for low-to-high multigrid pass",
            3,
        ),
        "bah_enable_BBH_mode": (
            "BBH mode? Enable=1 ; Disable=0",
            0,
        ),
        "bah_BBH_mode_common_horizon_idx": (
            "BBH mode: Horizon index (zero-offset) for common (merged BH) horizon",
            2,
        ),
        "bah_enable_eta_varying_alg_for_precision_common_horizon": (
            "Enable varying-eta prescription for precision common-horizon finding",
            0,
        ),
    }
    for key, value in int_params.items():
        _ = par.register_CodeParameter(
            cparam_type="int",
            module="BHaHAHA",
            name=key,
            defaultvalue=value[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=value[0],
        )

    # Dictionary for INT array type parameters
    int_array_params: Dict[str, Tuple[str, List[int]]] = {
        "bah_Ntheta_array_multigrid": (
            "Array of theta resolutions for multigrid cycle",
            [8, 16, 32],
        ),
        "bah_Nphi_array_multigrid": (
            "Array of phi resolutions for multigrid cycle",
            [16, 32, 64],
        ),
        "bah_BBH_mode_inspiral_BH_idxs": (
            "BBH mode: Horizon indices (zero-offset) for inspiralling BHs",
            [0, 1],
        ),
    }
    for key, intvalue in int_array_params.items():
        _ = par.register_CodeParameter(
            cparam_type=f"int[{len(intvalue[1])}]",
            module="BHaHAHA",
            name=key,
            defaultvalue=intvalue[1],
            commondata=True,
            add_to_parfile=True,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=True,
            description=intvalue[0],
        )

    _ = par.register_CodeParameter(
        cparam_type=f"int[{max_horizons}]",
        module="BHaHAHA",
        name="bah_BBH_mode_horizon_active",
        defaultvalue=[1] * max_horizons,
        commondata=True,
        add_to_parfile=False,
        add_to_set_CodeParameters_h=False,
        add_to_glb_code_params_dict=True,
        description="BBH mode: Record of which horizons are active.",
    )


def string_for_spherical_interp_setup_steps_1_to_4() -> str:
    r"""
    Generate the C string for spherical interpolation setup (steps 1–4).

    :return: Raw C string.
    """
    outstring = r"""
  // STEP 1: Determine spherical grid parameters and total interpolation points.
  const int Ntheta_interp = current_horizon_params->Ntheta_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
  const int Nphi_interp = current_horizon_params->Nphi_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
  const REAL dtheta_interp = M_PI / ((REAL)Ntheta_interp);
  const REAL dphi_interp = 2.0 * M_PI / ((REAL)Nphi_interp);

  const int actual_Nr_interp = current_horizon_params->Nr_external_input;
  *total_interp_points = actual_Nr_interp * Ntheta_interp * Nphi_interp;

  // STEP 3: Allocate memory for destination reference-metric coordinates.
  *dst_x0x1x2_interp = (REAL(*)[3])malloc(*total_interp_points * 3 * sizeof(REAL));
  if (dst_x0x1x2_interp == NULL) {
    fprintf(stderr, "ERROR: Failed to allocate memory for dst_x0x1x2_interp.\n");
    exit(EXIT_FAILURE);
  } // END IF: dst_x0x1x2_interp == NULL

#define IDX3_SPH_INTERP_LOCAL(ir, itheta, iphi) ((ir) + actual_Nr_interp * ((itheta) + Ntheta_interp * (iphi)))

  // STEP 4: Populate `dst_x0x1x2_interp`.
#pragma omp parallel for
  for (int iphi = 0; iphi < Nphi_interp; iphi++) {
    const REAL phi = -M_PI + ((REAL)iphi + 0.5) * dphi_interp;
    const REAL sinphi = sin(phi);
    const REAL cosphi = cos(phi);
    for (int itheta = 0; itheta < Ntheta_interp; itheta++) {
      const REAL theta = ((REAL)itheta + 0.5) * dtheta_interp;
      const REAL sintheta = sin(theta);
      const REAL costheta = cos(theta);
      for (int ir = 0; ir < actual_Nr_interp; ir++) {
        const REAL r = radii[ir];
        const int idx3 = IDX3_SPH_INTERP_LOCAL(ir, itheta, iphi);
        const REAL xCart[3] = {x_center + r * sintheta * cosphi, y_center + r * sintheta * sinphi, z_center + r * costheta};
        int Cart_to_i0i1i2_not_stored_to_save_memory[3];
        Cart_to_xx_and_nearest_i0i1i2(params, xCart, (*dst_x0x1x2_interp)[idx3], Cart_to_i0i1i2_not_stored_to_save_memory);
      } // END LOOP: for ir (spherical grid setup)
    } // END LOOP: for itheta (spherical grid setup)
  } // END LOOP: for iphi (#pragma omp parallel for, spherical grid setup)
"""
    return outstring


def string_for_spherical_interp_setup_step_6_allocate_tmp_bssn() -> str:
    r"""
    Generate the C string for spherical interpolation setup (step 6).

    :return: Raw C string.
    """
    outstring = r"""
// STEP 6: Allocate temporary memory for interpolated BSSN variables.
  *dst_data_ptrs_bssn = (REAL **)malloc(BHAHAHA_NUM_INTERP_GFS * sizeof(REAL *));
  for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
    (*dst_data_ptrs_bssn)[i] = (REAL *)malloc(*total_interp_points * sizeof(REAL));
    if ((*dst_data_ptrs_bssn)[i] == NULL) {
      fprintf(stderr, "ERROR: Failed to allocate memory for dst_data_ptrs_bssn[%d].\n", i);
      for (int k = 0; k < i; ++k)
        BHAH_FREE((*dst_data_ptrs_bssn)[k]);
      BHAH_FREE(*dst_x0x1x2_interp);
      exit(EXIT_FAILURE);
    } // END IF: dst_data_ptrs_bssn[i] == NULL
  } // END LOOP: for i (allocating dst_data_ptrs_bssn)
"""
    return outstring


def string_for_bssn_to_adm_transformation_block(CoordSystem: str) -> str:
    r"""
    Generate the C string for STEP 8: BSSN→ADM Cartesian transformation block.

    :param CoordSystem: CoordSystem of project, where horizon finding will take place.
    :return: Raw C string to be injected as the transformation block.
    """
    outstring = r"""
{                             // Start of BSSN to ADM transformation block
#include "set_CodeParameters.h" // NRPy-specific include for coordinate transformations and symbolic expressions

    // STEP 8: Transform interpolated BSSN data to ADM Cartesian components.
#pragma omp parallel for
    for (int iphi = 0; iphi < Nphi_interp; iphi++) {
      for (int itheta = 0; itheta < Ntheta_interp; itheta++) {
        for (int ir = 0; ir < actual_Nr_interp; ir++) {
          const int offset = total_interp_points;
          const int idx3 = IDX3_SPH_INTERP_LOCAL(ir, itheta, iphi);
          const REAL xx0 = (*dst_x0x1x2_interp)[idx3][0];
          const REAL xx1 = (*dst_x0x1x2_interp)[idx3][1];
          const REAL xx2 = (*dst_x0x1x2_interp)[idx3][2];

          const REAL cf = (*dst_data_ptrs_bssn)[INTERP_CFGF_IDX][idx3];
          const REAL trK = (*dst_data_ptrs_bssn)[INTERP_TRKGF_IDX][idx3];
"""
    defines_list: List[str] = []
    for i in range(3):
        for j in range(i, 3):
            defines_list += [
                f"          const REAL rfm_hDD{i}{j} = (*dst_data_ptrs_bssn)[INTERP_HDD{i}{j}GF_IDX][idx3];\n"
            ]
            defines_list += [
                f"          const REAL rfm_aDD{i}{j} = (*dst_data_ptrs_bssn)[INTERP_ADD{i}{j}GF_IDX][idx3];\n"
            ]
    outstring += "".join(sorted(defines_list, key=str.casefold))

    outstring += generate_bssn_to_adm_codegen(CoordSystem)

    outstring += r"""
        } // END LOOP: for ir (BSSN to ADM transformation)
      } // END LOOP: for itheta (BSSN to ADM transformation)
    } // END LOOP: for iphi (#pragma omp parallel for, BSSN to ADM transformation)
  } // End of BSSN to ADM transformation block
"""
    return outstring


def build_bhahaha_prefunc(
    CoordSystem: str,
    add_bhahaha_gf_interp_ind_to_bhah_defines: bool = False,
) -> str:
    """
    Construct the C prelude used by the generated horizon-finder function.
    Includes:
      * FINAL_ADM_METRIC_INDICES and INTERP_BSSN_GF_INDICES enums
      * bhahaha_gf_interp_indices[] mapping
      * input validation (check_multigrid_resolution_inputs)
      * per-horizon init/free helpers
      * interpolation + BSSN→ADM transformation routine

    :param CoordSystem: CoordSystem of project, where horizon finding will take place.
    :param add_bhahaha_gf_interp_ind_to_bhah_defines: add to BHaH_defines.h instead of defining here.
    :return: Raw C string to be injected as the function preamble.
    """
    prefunc = string_for_prefunc_enums_and_interp_indices(
        add_bhahaha_gf_interp_ind_to_bhah_defines=add_bhahaha_gf_interp_ind_to_bhah_defines
    )

    prefunc += string_for_static_func_check_multigrid_resolution_inputs()

    prefunc += string_for_static_func_initialize_bhahaha_solver_params_and_shapes()

    prefunc += string_for_func_free_bhahaha_horizon_shape_data_all_horizons()

    prefunc += r"""

/**
 * Prepares spherical interpolation points and allocates storage
 * for interpolated BSSN gridfunctions.
 *
 * Steps:
 * 1. Determine spherical grid parameters (Ntheta, Nphi, dtheta, dphi from
 *    `current_horizon_params`) and compute total number of interpolation points.
 * 2. Return if `*total_interp_points` is zero.
 * 3. Allocate memory for destination reference-metric coordinates
 *    (`dst_x0x1x2_interp`). Exits on failure.
 * 4. Populate `dst_x0x1x2_interp`: For each spherical grid point
 *    (r, theta, phi around x_center, y_center, z_center), convert
 *    Cartesian coordinates to reference-metric coordinates using
 *    `Cart_to_xx_and_nearest_i0i1i2`.
 * 5. Allocate arrays for interpolated BSSN variables
 *    (`dst_data_ptrs_bssn`). Exits on failure.
 *
 * NOTE: This function only sets up spherical target points and allocates
 *       buffers. Interpolation and ADM transformation are done elsewhere.
 */
static void BHaHAHA_compute_dst_pts_metric_data_nrpy(const params_struct *restrict params,
                                                 bhahaha_params_and_data_struct *restrict current_horizon_params, const REAL x_center,
                                                 const REAL y_center, const REAL z_center, const REAL radii[],
                                                 int *total_interp_points, REAL (**dst_x0x1x2_interp)[3], REAL ***dst_data_ptrs_bssn) {



"""

    prefunc += string_for_spherical_interp_setup_steps_1_to_4()

    prefunc += string_for_spherical_interp_setup_step_6_allocate_tmp_bssn()

    prefunc += r"""
} // END FUNCTION: BHaHAHA_compute_dst_pts_metric_data_nrpy
"""

    prefunc += r"""
/**
 * Transforms pre-interpolated BSSN data (on spherical target points) into
 * ADM Cartesian components and writes them into the flat array for BHaHAHA.
 *
 * Steps:
 * 8. Using the precomputed reference-metric coords (`*dst_x0x1x2_interp`) and
 *    pre-interpolated BSSN variables (`*dst_data_ptrs_bssn`), compute ADM
 *    Cartesian g_ij and K_ij at each target point (parallelized), and store into
 *    `input_metric_data_target_array` in the expected flat layout.
 * 10. Free allocated temporary buffers (`*dst_x0x1x2_interp`, `*dst_data_ptrs_bssn`).
 *
 * @param commondata - Pointer to `commondata_struct` (used by set_CodeParameters.h).
 * @param params - Pointer to `params_struct` (used by set_CodeParameters.h).
 * @param current_horizon_params - Horizon grid sizes (Nr, Ntheta, Nphi).
 * @param dst_x0x1x2_interp - Ref-metric coords at target points (size: Nr*Ntheta*Nphi × 3).
 * @param dst_data_ptrs_bssn - Pre-interpolated BSSN values at target points.
 * @param input_metric_data_target_array - Output flat array for ADM g_ij and K_ij.
 * @return - None (`void`).
 */
static void BHaHAHA_transform_BSSN_to_ADM(
    const commondata_struct *restrict commondata,
    const params_struct *restrict params,
    bhahaha_params_and_data_struct *restrict current_horizon_params,
    REAL (**dst_x0x1x2_interp)[3],
    REAL ***dst_data_ptrs_bssn,
    REAL *restrict input_metric_data_target_array)
{

  //Determine spherical grid parameters and total interpolation points.
  const int Ntheta_interp = current_horizon_params->Ntheta_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
  const int Nphi_interp = current_horizon_params->Nphi_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
  const REAL dtheta_interp = M_PI / ((REAL)Ntheta_interp);
  const REAL dphi_interp = 2.0 * M_PI / ((REAL)Nphi_interp);

  const int actual_Nr_interp = current_horizon_params->Nr_external_input;
  const int total_interp_points = actual_Nr_interp * Ntheta_interp * Nphi_interp;

#define IDX3_SPH_INTERP_LOCAL(ir, itheta, iphi) ((ir) + actual_Nr_interp * ((itheta) + Ntheta_interp * (iphi)))
"""

    prefunc += string_for_bssn_to_adm_transformation_block(CoordSystem)

    prefunc += r"""
  // STEP 10: Free allocated temporary memory.
  BHAH_FREE(*dst_x0x1x2_interp);
  for (int i = 0; i < BHAHAHA_NUM_INTERP_GFS; i++) {
    BHAH_FREE((*dst_data_ptrs_bssn)[i]);
  } // END LOOP: for i (freeing (*dst_data_ptrs_bssn))
  BHAH_FREE(*dst_data_ptrs_bssn);
} // END FUNCTION: BHaHAHA_interpolate_metric_data_nrpy
"""
    return prefunc


def register_CFunction_bhahaha_find_horizons(
    CoordSystem: str,
    max_horizons: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 3D Lagrange interpolation.

    :param CoordSystem: CoordSystem of project, where horizon finding will take place.
    :param max_horizons: Maximum number of horizons to search for.
    :return: None if in registration phase, else the updated NRPy environment.

    >>> env = register_CFunction_bhahaha_find_horizons()
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    register_bhahaha_commondata_and_params(max_horizons)

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "sys/time.h",
    ]

    prefunc = build_bhahaha_prefunc(
        CoordSystem, add_bhahaha_gf_interp_ind_to_bhah_defines=True
    )

    desc = r"""Main driver function for finding apparent horizons using the BHaHAHA library.
It orchestrates initialization, extrapolation, spherical-target setup/staging, solving,
and result updates for multiple horizons.

- The function performs the following steps:
- 1. Checks if a horizon find is scheduled for the current iteration based on
-    `diagnostics_output_every` and `dt`. Returns if not scheduled.
- 2. Validates multigrid resolution inputs using `check_multigrid_resolution_inputs`. Exits on failure.
- 3. Initializes BHaHAHA data structures and solver parameters:
-    a. If it's the first call (nn=0):
-       i.  Initializes persistent fields (centers, times, radii) in `commondata->bhahaha_params_and_data`
-           for each horizon using initial guess parameters.
-       ii. Calls `initialize_bhahaha_solver_params_and_shapes` with `is_first_call_for_shapes = true`
-           to allocate shape history arrays and set non-persistent solver params.
-       iii. Marks all configured horizons active.
-    b. If not the first call (nn!=0):
-       i.  Calls `initialize_bhahaha_solver_params_and_shapes` with `is_first_call_for_shapes = false`
-           to refresh non-persistent solver params (shape arrays are not re-allocated).
- 4. Populates current iteration metadata (iteration number `nn`, simulation `time`) into each
-    horizon's `bhahaha_params_and_data_struct`.
- 5. Applies robustness improvements and extrapolates horizon guesses for the current find:
-    a. For each horizon, if it was not found reliably in the previous three attempts (t_m1, t_m2, or t_m3 is -1.0):
-       i.  Reduces its `cfl_factor` by 10%.
-       ii. Sets `use_fixed_radius_guess_on_full_sphere` to true.
-    b. Extrapolates current guesses for center (x_guess, y_guess, z_guess) and radial search range
-       (r_min_guess, r_max_guess) using `bah_xyz_center_r_minmax` based on historical data.
-       If `use_fixed_radius_guess_on_full_sphere` is true for a horizon, its r_min_guess is set to 0.0 and
-       r_max_guess is set to `commondata->bah_max_search_radius[h]`.
- 6. Main computation (single-horizon entry; `h = which_horizon`):
-    a. Sets up the radial grid for interpolation using `bah_radial_grid_cell_centered_set_up`,
-       passing the current `r_min_guess[h]`, `r_max_guess[h]`, and storing the resulting
-       actual `Nr_external_input`, `r_min_external_input`, `dr_external_input`, and `radii_interp`
-       array in/via `current_horizon_params`.
-    b. Allocates a buffer for interpolated ADM metric data (`current_horizon_params->input_metric_data`)
-       if `Nr_external_input` > 0. Exits on allocation failure.
-    c. If the metric data buffer was allocated and `Nr_external_input` > 0, calls
-       `BHaHAHA_compute_dst_pts_metric_data_nrpy` to prepare spherical target points and per-point BSSN storage
-       using the current `x_guess[h]`, `y_guess[h]`, `z_guess[h]`, and the populated `radii_interp`.
-       The BSSN→ADM packing occurs later via `BHaHAHA_transform_BSSN_to_ADM`.
- 7. Solve path:
-    a. If metric data is available (`input_metric_data` is not NULL and `Nr_external_input` > 0):
-       i.  Calls `BHaHAHA_transform_BSSN_to_ADM` and then the solver `bah_find_horizon`, passing
-           `current_horizon_params` and a local `bhahaha_diagnostics_struct`.
-       ii. If the solver returns `BHAHAHA_SUCCESS`:
-           - Updates `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 0 (false).
-           - Outputs diagnostics using `bah_diagnostics_file_output`.
-       iii. If the solver fails:
-            - Prints a warning.
-            - Sets `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 1 (true).
-            - Sets `current_horizon_params->t_m1` to -1.0 to indicate failure for future extrapolation.
-    b. If the solve was skipped (e.g., `Nr_external_input == 0` or no metric data buffer):
-       - Sets `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 1.
-       - Sets `current_horizon_params->t_m1` to -1.0.
-    c. Frees the allocated `input_metric_data` buffer for the current horizon.
- 8. Emits start/finish messages when verbosity is enabled (timing may be reported elsewhere).

@param commondata - Pointer to `commondata_struct` holding global BHaHAHA settings,
                    persistent horizon data, and simulation state. This struct is
                    extensively read from and modified.
@param griddata - Pointer to `griddata_struct` containing simulation grid information
                  and gridfunctions for metric data interpolation.
@return - None (`void`).
"""

    cfunc_type = "void"
    name = "bhahaha_find_horizons"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata, REAL x_guess[], REAL y_guess[], REAL z_guess[],
                           REAL **radii_interp, int *total_elements, REAL (**dst_x0x1x2)[3], REAL ***dst_data_ptrs, const int which_horizon,
                           const int which_bhahaha_part"""
    body = """
    switch (which_bhahaha_part) {
    case BHAHAHA_FIND_HORIZONS_SETUP: {
    """
    body += """
    const int max_num_horizons = commondata->bah_max_num_horizons;"""

    body += string_for_step2_validate_multigrid_inputs()

    body += string_for_step3_initialize_bhahaha_data_structs_and_solver_params(
        enable_BBH_mode=False
    )

    body += r"""
  // Local arrays for per-horizon guesses for the current find.
  REAL r_min_guess[max_num_horizons], r_max_guess[max_num_horizons];
  """

    body += string_for_step4_populate_iteration_metadata()

    body += string_for_step6_apply_robustness_improv_and_extrap_horizon_guesses(
        enable_BBH_mode=False
    )

    body += r"""
  if (commondata->bah_verbosity_level > 1) {
    printf("--- NRPy BHaHAHA Pre-Interpolation State (Iter %d, Time %.3f) ---\n", commondata->nn, commondata->time);
    for (int h_diag = 0; h_diag < commondata->bah_max_num_horizons; h_diag++) {
      const bhahaha_params_and_data_struct *p = &commondata->bhahaha_params_and_data[h_diag];
      printf("  H%d: act=%d, fixed_sph=%d, ctr_guess=(%.3e,%.3e,%.3e), r_guess=[%.3e,%.3e]\n"
             "     state_vars: x_c_m1=%.3e, y_c_m1=%.3e, z_c_m1=%.3e, r_max_m1=%.3e, r_min_m1=%.3e, t_m1=%.3f\n",
             h_diag, commondata->bah_BBH_mode_horizon_active[h_diag], // Activity from commondata array
             p->use_fixed_radius_guess_on_full_sphere,                // Flag from BHaHAHA struct
             x_guess[h_diag], y_guess[h_diag], z_guess[h_diag], r_min_guess[h_diag], r_max_guess[h_diag], p->x_center_m1, p->y_center_m1,
             p->z_center_m1, p->r_max_m1, p->r_min_m1, p->t_m1);
    } // END LOOP: for h_diag (pre-interpolation print)
    printf("-----------------------------------------------------------\n");
  } // END IF: verbosity_level > 1 (pre-interpolation print)
  """

    body += string_for_step7a_to_c_main_loop_for_each_horizon(
        single_horizon=True,
        allocate_radii_interp=True,
        enable_timing=False,
        enable_BBH_mode=False,
    )

    body += """
    // STEP 7.d: Interpolate metric data.
    if (current_horizon_params->input_metric_data && current_horizon_params->Nr_external_input > 0) {
      BHaHAHA_compute_dst_pts_metric_data_nrpy(&griddata[0].params, current_horizon_params, x_guess[h], y_guess[h],
                                           z_guess[h],   // Use per-horizon guess
                                           *radii_interp, // Pass the populated radii_interp array for this horizon
                                           total_elements, dst_x0x1x2, dst_data_ptrs);
    } // END IF: call BHaHAHA_compute_dst_pts_metric_data_nrpy
    """

    body += r"""
    break;
    }
    case BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE: {
      int h = which_horizon;
      bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];
      if (current_horizon_params->input_metric_data && current_horizon_params->Nr_external_input > 0) {
        BHaHAHA_transform_BSSN_to_ADM(commondata, &griddata[0].params, current_horizon_params, dst_x0x1x2, dst_data_ptrs, current_horizon_params->input_metric_data);
      } // END IF: call BHaHAHA_transform_BSSN_to_ADM
  """

    body += string_for_step7e_to_g_main_loop_for_each_horizon(
        single_horizon=True, enable_timing=False
    )

    body += r"""
    break;
    }
    }
  """

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

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
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.infrastructures.BHaH.BHaHAHA.BHaH_implementation import (
    build_bhahaha_prefunc,
    register_bhahaha_commondata_and_params,
)


def register_CFunction_bhahaha_find_horizons(
    CoordSystem: str,
    max_horizons: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for general-purpose 3D Lagrange interpolation.

    :param CoordSystem: CoordSystem of project, where horizon finding will take place.
    :param max_horizons: Maximum number of horizons to search for.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If EvolvedConformalFactor_cf set to unsupported value.

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
It orchestrates initialization, BBH logic, extrapolation, interpolation, solving,
and result updates for multiple horizons.

- The function performs the following steps:
- 1. Checks if a horizon find is scheduled for the current iteration based on
-    `diagnostics_output_every` and `dt`. Returns if not scheduled.
- 2. Validates multigrid resolution inputs using `check_multigrid_resolution_inputs`. Exits on failure.
- 3. Initializes BHaHAHA data structures and solver parameters:
-    a. If it's the first call (nn=0):
-       i. Initializes persistent fields (centers, times, radii) in `commondata->bhahaha_params_and_data`
-          for each horizon using initial guess parameters.
-       ii. Calls `initialize_bhahaha_solver_params_and_shapes` with `is_first_call_for_shapes = true`
-           to allocate shape history arrays and set non-persistent solver params.
-       iii. Initializes `commondata->bah_BBH_mode_horizon_active` based on `bah_BBH_mode_enable`.
-           Exits if BBH mode is enabled with an invalid `bah_max_num_horizons` or invalid/overlapping indices.
-    b. If not the first call:
-       i. Calls `initialize_bhahaha_solver_params_and_shapes` with `is_first_call_for_shapes = false`
-          to re-initialize non-persistent solver params (shape arrays are not re-allocated).
- 4. Populates current iteration metadata (iteration number `nn`, simulation `time`) into each
-    horizon's `bhahaha_params_and_data_struct`.
- 5. Applies BBH mode logic if `commondata->bah_BBH_mode_enable` is true:
-    a. Deactivates individual inspiral BH searches if the common horizon was found previously and is active.
-    b. Activates the common horizon search if both individual BHs are active, the common horizon is not yet active,
-       and the individual BHs (based on their previous find's center and max radius) meet a proximity criterion
-       (distance + radii sum <= threshold diameter).
-       If activated, initializes the common horizon's persistent state (center from mass-weighted average,
-       time set to -1.0 for fixed radius guess, max radius from parameters) and sets its current find guess.
- 6. Applies robustness improvements and extrapolates horizon guesses for the current find:
-    a. For each horizon, if it was not found reliably in the previous three attempts (t_m1, t_m2, or t_m3 is -1.0):
-       i. Reduces its `cfl_factor` by 10%.
-       ii. Sets `use_fixed_radius_guess_on_full_sphere` to true.
-       iii. If in BBH mode and it's the common horizon making its initial attempts, increases the resolution
-           of its first multigrid level and doubles `max_iterations`.
-    b. Extrapolates current guesses for center (x_guess, y_guess, z_guess) and radial search range
-       (r_min_guess, r_max_guess) using `bah_xyz_center_r_minmax` based on historical data.
-       If `use_fixed_radius_guess_on_full_sphere` is true for a horizon, its r_min_guess is set to 0.0 and
-       r_max_guess is set to `commondata->bah_max_search_radius[h]`.
- 7. Main loop for each horizon `h` from 0 to `max_num_horizons - 1`:
-    a. Skips the horizon if it's not marked active in `commondata->bah_BBH_mode_horizon_active[h]`.
-    b. Sets up the radial grid for interpolation using `bah_radial_grid_cell_centered_set_up`,
-       passing the current `r_min_guess[h]`, `r_max_guess[h]`, and storing the resulting
-       actual `Nr_external_input`, `r_min_external_input`, `dr_external_input`, and `radii_interp`
-       array in/via `current_horizon_params`.
-    c. Allocates a buffer for interpolated ADM metric data (`current_horizon_params->input_metric_data`)
-       if `Nr_external_input` > 0. Exits on allocation failure.
-    d. If the metric data buffer was allocated and `Nr_external_input` > 0, calls
-       `BHaHAHA_interpolate_metric_data_nrpy` to fill the buffer using the current `x_guess[h]`,
-       `y_guess[h]`, `z_guess[h]`, and the populated `radii_interp`.
-    e. If metric data is available (`input_metric_data` is not NULL and `Nr_external_input` > 0):
-       i. Calls the BHaHAHA solver `bah_find_horizon`, passing `current_horizon_params` (which
-          contains all input data and state) and a local `bhahaha_diagnostics_struct`.
-       ii. If the solver returns `BHAHAHA_SUCCESS`:
-           - Updates `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 0 (false).
-           - Outputs diagnostics using `bah_diagnostics_file_output`.
-       iii. If the solver fails:
-            - Prints a warning.
-            - Sets `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 1 (true).
-            - Sets `current_horizon_params->t_m1` to -1.0 to indicate failure for future extrapolation.
-    f. If the solve was skipped (e.g., `Nr_external_input == 0` or no metric data buffer):
-       - Sets `current_horizon_params->use_fixed_radius_guess_on_full_sphere` to 1.
-       - Sets `current_horizon_params->t_m1` to -1.0.
-    g. Frees the allocated `input_metric_data` buffer for the current horizon.
- 8. Prints total elapsed time for the entire horizon finding process if verbosity is enabled.

@param commondata - Pointer to `commondata_struct` holding global BHaHAHA settings,
                    persistent horizon data, simulation state, and BBH mode flags.
                    This struct is extensively read from and modified.
@param griddata - Pointer to `griddata_struct` containing simulation grid information
                  and gridfunctions for metric data interpolation.
@return - None (`void`).
"""

    cfunc_type = "void"
    name = "bhahaha_find_horizons"
    params = (
        """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""
    )

    body = r"""
  // STEP 1: Check if horizon find is scheduled for the current iteration.
  if (commondata->diagnostics_output_every <= 0 || (commondata->nn % (int)(commondata->diagnostics_output_every / commondata->dt + 0.5)) != 0) {
    int bah_find_every = 1;                                      // Placeholder: find every iteration. This should be a commondata param.
    if (commondata->diagnostics_output_every > commondata->dt) { // A basic way to get find_every from time interval
      bah_find_every = (int)(commondata->diagnostics_output_every / commondata->dt + 0.5);
      if (bah_find_every == 0)
        bah_find_every = 1;
    }
    if (bah_find_every <= 0 || (commondata->nn % bah_find_every) != 0) {
      return; // END IF: not scheduled this iteration
    }
  } // END IF: diagnostics_output_every > 0

  const int max_num_horizons = commondata->bah_max_num_horizons;
  if (max_num_horizons <= 0) {
    return; // END IF: no horizons
  }

  struct timeval start_time_total, iter_time_tracker;
  gettimeofday(&start_time_total, NULL);

  // STEP 2: Validate multigrid resolution inputs.
  check_multigrid_resolution_inputs(commondata); // Validates Ntheta, Nphi arrays in commondata

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

    // STEP 3.a.iii: Initialize `commondata->bah_BBH_mode_horizon_active`.
    if (commondata->bah_BBH_mode_enable) {
      if (commondata->bah_max_num_horizons != 3) {
        fprintf(stderr,
                "ERROR: bah_BBH_mode_enable requires bah_max_num_horizons==3, to account for common horizon, plus two individual horizons.\n");
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
    } // END ELSE: not BBH_mode_enable (iteration 0 horizon activity)
  } else { // Not the first call (nn != 0)
    // STEP 3.b: If not the first call:
    for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
      // STEP 3.b.i: Call `initialize_bhahaha_solver_params_and_shapes` for non-persistent params.
      initialize_bhahaha_solver_params_and_shapes(commondata, h, false); // false: don't re-alloc shapes
    } // END LOOP: for h (non-first call initialization of solver params)
  } // END ELSE: not commondata->nn == 0

  // Local arrays for per-horizon guesses for the current find.
  REAL x_guess[max_num_horizons], y_guess[max_num_horizons], z_guess[max_num_horizons];
  REAL r_min_guess[max_num_horizons], r_max_guess[max_num_horizons];

  // STEP 4: Populate current iteration metadata.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];
    current_horizon_params->which_horizon = h + 1; // 1-based for BHaHAHA library
    current_horizon_params->num_horizons = commondata->bah_max_num_horizons;
    current_horizon_params->iteration_external_input = commondata->nn;
    current_horizon_params->time_external_input = commondata->time;
  } // END LOOP: for h (populating current iteration metadata)

  // STEP 5: Apply BBH mode logic.
  if (commondata->bah_BBH_mode_enable) {
    const int bh1_idx = commondata->bah_BBH_mode_inspiral_BH_idxs[0];
    const int bh2_idx = commondata->bah_BBH_mode_inspiral_BH_idxs[1];
    const int com_idx = commondata->bah_BBH_mode_common_horizon_idx;

    const bool com_found_prev = (commondata->bhahaha_params_and_data[com_idx].use_fixed_radius_guess_on_full_sphere == 0);

    // STEP 5.a: Deactivate individual BHs if common horizon was found previously and is active.
    if (commondata->bah_BBH_mode_horizon_active[com_idx] && com_found_prev && commondata->bah_BBH_mode_horizon_active[bh1_idx] &&
        commondata->bah_BBH_mode_horizon_active[bh2_idx]) {
      commondata->bah_BBH_mode_horizon_active[bh1_idx] = 0;
      commondata->bah_BBH_mode_horizon_active[bh2_idx] = 0;
      if (commondata->bah_verbosity_level > 0) {
        printf("NRPy BHaHAHA BBH (Iter %d): Common active and found. Deactivating individual BHs (H%d, H%d).\n", commondata->nn, bh1_idx, bh2_idx);
      } // END IF: verbosity for deactivation
    } // END IF: common horizon found, deactivate individual

    // STEP 5.b: Activate common horizon search if criteria met.
    if (commondata->bah_BBH_mode_horizon_active[bh1_idx] && commondata->bah_BBH_mode_horizon_active[bh2_idx] &&
        !commondata->bah_BBH_mode_horizon_active[com_idx]) {

      const bhahaha_params_and_data_struct *params_bh1 = &commondata->bhahaha_params_and_data[bh1_idx];
      const bhahaha_params_and_data_struct *params_bh2 = &commondata->bhahaha_params_and_data[bh2_idx];

      const REAL x1 = params_bh1->x_center_m1, y1 = params_bh1->y_center_m1, z1 = params_bh1->z_center_m1;
      const REAL rmax1 = params_bh1->r_max_m1, t1 = params_bh1->t_m1;
      const REAL x2 = params_bh2->x_center_m1, y2 = params_bh2->y_center_m1, z2 = params_bh2->z_center_m1;
      const REAL rmax2 = params_bh2->r_max_m1, t2 = params_bh2->t_m1;

      if (t1 >= 0.0 && t2 >= 0.0) { // Both BH1 and BH2 were found previously.
        const REAL dist_centers = sqrt(SQR(x1 - x2) + SQR(y1 - y2) + SQR(z1 - z2));
        const REAL threshold_diam = 2.0 * commondata->bah_max_search_radius[com_idx]; // Diameter for common horizon search area.

        if (commondata->bah_verbosity_level > 0) {
          printf("NRPy BBH Trigger (Iter %d, CommonH_idx %d): sep=%.3f, rmax1(H%d)=%.3f, rmax2(H%d)=%.3f. Sum_dist+radii=%.3f <= Thr_diam=%.3f?\n",
                 commondata->nn, com_idx, dist_centers, bh1_idx, rmax1, bh2_idx, rmax2, (dist_centers + rmax1 + rmax2), threshold_diam);
        } // END IF: verbosity for trigger check

        if (dist_centers + rmax1 + rmax2 <= threshold_diam) {   // Trigger condition met.
          commondata->bah_BBH_mode_horizon_active[com_idx] = 1; // Activate common horizon search.

          const REAL M1 = params_bh1->M_scale, M2 = params_bh2->M_scale;
          const REAL M_sum = (M1 + M2); // This will never sum to zero, or we have a much bigger problem on our hands.
          const REAL com_x = (M1 * x1 + M2 * x2) / M_sum;
          const REAL com_y = (M1 * y1 + M2 * y2) / M_sum;
          const REAL com_z = (M1 * z1 + M2 * z2) / M_sum;

          bhahaha_params_and_data_struct *params_com = &commondata->bhahaha_params_and_data[com_idx];
          params_com->x_center_m1 = com_x;
          params_com->y_center_m1 = com_y;
          params_com->z_center_m1 = com_z;
          params_com->t_m1 = params_com->t_m2 = params_com->t_m3 = -1.0; // Ensures fixed radius guess via extrapolation logic.
          params_com->r_min_m1 = 0.0;
          params_com->r_max_m1 = commondata->bah_max_search_radius[com_idx]; // Set r_max_m1 for extrapolation.
          params_com->use_fixed_radius_guess_on_full_sphere = 1;             // Explicitly force full sphere for its first find.

          x_guess[com_idx] = com_x;
          y_guess[com_idx] = com_y;
          z_guess[com_idx] = com_z;
          r_min_guess[com_idx] = 0.0;
          r_max_guess[com_idx] = commondata->bah_max_search_radius[com_idx];

          if (commondata->bah_verbosity_level > 0) {
            printf("NRPy BBH Activate Common (Iter %d, H%d): Center_guess=(%.3f,%.3f,%.3f), r_max_guess=%.3f\n", commondata->nn, com_idx, com_x,
                   com_y, com_z, r_max_guess[com_idx]);
          } // END IF: verbosity for common activation
        } // END IF: trigger condition met
      } // END IF: both inspiral BHs found previously
    } // END IF: check for activating common horizon
  } // END IF: bah_BBH_mode_enable (BBH Mode Logic)

  // STEP 6: Apply robustness improvements and extrapolate horizon guesses.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];

    // STEP 6.a: For each horizon, if not found reliably, adjust parameters.
    if (current_horizon_params->t_m1 == -1.0 || current_horizon_params->t_m2 == -1.0 || current_horizon_params->t_m3 == -1.0) {
      // STEP 6.a.i: Reduce CFL.
      current_horizon_params->cfl_factor *= 0.9;
      // STEP 6.a.ii: Force fixed radius guess.
      current_horizon_params->use_fixed_radius_guess_on_full_sphere = 1; // Force fixed radius if not found reliably.

      // STEP 6.a.iii: If BBH common horizon, adjust resolution/iterations.
      if (commondata->bah_BBH_mode_enable && h == commondata->bah_BBH_mode_common_horizon_idx &&
          current_horizon_params->num_resolutions_multigrid >= 2) {
        current_horizon_params->Ntheta_array_multigrid[0] = commondata->bah_Ntheta_array_multigrid[1]; // Use 2nd level res for 1st.
        current_horizon_params->Nphi_array_multigrid[0] = commondata->bah_Nphi_array_multigrid[1];
        current_horizon_params->max_iterations *= 2;
      } // END IF: BBH mode common horizon robustness for resolution/iterations
    } // END IF: robustness improvements for not-yet-found horizons

    // STEP 6.b: Extrapolate current guesses.
    bah_xyz_center_r_minmax(current_horizon_params, &x_guess[h], &y_guess[h], &z_guess[h], &r_min_guess[h], &r_max_guess[h]);
    if (current_horizon_params->use_fixed_radius_guess_on_full_sphere) { // This flag is key.
      r_min_guess[h] = 0.0;
      r_max_guess[h] = commondata->bah_max_search_radius[h]; // Use per-horizon parameter from commondata.
    } // END IF: use_fixed_radius_guess_on_full_sphere for extrapolation override
  } // END LOOP: for h (Robustness & Extrapolation)

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

  // STEP 7: Main loop for each horizon.
  for (int h = 0; h < commondata->bah_max_num_horizons; h++) {
    // STEP 7.a: Skip if horizon is not active.
    if (!commondata->bah_BBH_mode_horizon_active[h])
      continue; // END IF: not active, continue loop

    gettimeofday(&iter_time_tracker, NULL); // Per-horizon timer start

    bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];

    // STEP 7.b: Set up the radial grid for interpolation.
    REAL radii_interp[commondata->bah_Nr_interp_max];                           // Max capacity from commondata (scalar)
    bah_radial_grid_cell_centered_set_up(commondata->bah_Nr_interp_max,         // Max capacity of radii_interp array
                                         commondata->bah_max_search_radius[h],  // Overall max search radius for this horizon
                                         r_min_guess[h], r_max_guess[h], // Current guess for search range from extrapolation
                                         &current_horizon_params->Nr_external_input, // Output: Actual Nr used for this horizon
                                         &current_horizon_params->r_min_external_input, // Output: Actual r_min used
                                         &current_horizon_params->dr_external_input, // Output: Actual dr used
                                         radii_interp);                         // Output: populated radii for this horizon

    // STEP 7.c: Allocate buffer for interpolated ADM metric data.
    const int Ntheta_interp_eff = current_horizon_params->Ntheta_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
    const int Nphi_interp_eff = current_horizon_params->Nphi_array_multigrid[current_horizon_params->num_resolutions_multigrid - 1];
    const size_t npts_metric_adm =
        (size_t)BHAHAHA_NUM_METRIC_COMPONENTS * current_horizon_params->Nr_external_input * Ntheta_interp_eff * Nphi_interp_eff;

    if (current_horizon_params->Nr_external_input > 0) {
      current_horizon_params->input_metric_data = (REAL *)malloc(npts_metric_adm * sizeof(REAL));
      if (!current_horizon_params->input_metric_data) {
        fprintf(stderr, "ERROR: malloc() failed for input_metric_data for H%d (%zu REALs).\n", h, npts_metric_adm);
        exit(EXIT_FAILURE);
      } // END IF: malloc failed for input_metric_data
    } else { // Nr_external_input <= 0
      current_horizon_params->input_metric_data = NULL;
    } // END ELSE: Nr_external_input <= 0

    if (commondata->bah_verbosity_level >= 2) {
      struct timeval temp_time;
      gettimeofday(&temp_time, NULL);
      printf("%.6f s : START interpolate H%d\n", timeval_to_seconds(start_time_total, temp_time), h);
    } // END IF: verbosity for start interpolate

    // STEP 7.d: Interpolate metric data.
    if (current_horizon_params->input_metric_data && current_horizon_params->Nr_external_input > 0) {
      BHaHAHA_interpolate_metric_data_nrpy(commondata, &griddata[0].params, griddata[0].xx, griddata[0].gridfuncs.y_n_gfs, current_horizon_params,
                                           x_guess[h], y_guess[h],
                                           z_guess[h],   // Use per-horizon guess
                                           radii_interp, // Pass the populated radii_interp array for this horizon
                                           current_horizon_params->input_metric_data);
    } // END IF: call BHaHAHA_interpolate_metric_data_nrpy

    if (commondata->bah_verbosity_level >= 2) {
      struct timeval temp_time;
      gettimeofday(&temp_time, NULL);
      printf("%.6f s : FINISH interpolate H%d (%.3fs)\n", timeval_to_seconds(start_time_total, temp_time), h,
             timeval_to_seconds(iter_time_tracker, temp_time));
      gettimeofday(&iter_time_tracker, NULL); // Reset for solver timing
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

    if (commondata->bah_verbosity_level >= 1) {
      struct timeval temp_time;
      gettimeofday(&temp_time, NULL);
      printf("%.6f s : FINISH find H%d (%.3fs for solve section)\n", timeval_to_seconds(start_time_total, temp_time), h,
             timeval_to_seconds(iter_time_tracker, temp_time));
    } // END IF: verbosity for finish find horizon
  } // END LOOP: for h (main find loop)

  // STEP 8: Print total elapsed time.
  if (commondata->bah_verbosity_level >= 0) {
    struct timeval end_time_total;
    gettimeofday(&end_time_total, NULL);
    printf("NRPy_BHaHAHA total elapsed time (Iter %d): %.6f s\n", commondata->nn, timeval_to_seconds(start_time_total, end_time_total));
  } // END IF: verbosity for total time print
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

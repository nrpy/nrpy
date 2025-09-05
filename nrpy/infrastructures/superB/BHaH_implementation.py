"""
Register C++ function for finding horizons with superB BHaHAHA.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures.BHaH.BHaHAHA.BHaH_implementation import (
    generate_bssn_to_adm_codegen,
    register_bhahaha_commondata_and_params,
    string_for_func_free_bhahaha_horizon_shape_data_all_horizons,
    string_for_prefunc_enums_and_interp_indices,
    string_for_static_func_check_multigrid_resolution_inputs,
    string_for_static_func_initialize_bhahaha_solver_params_and_shapes,
    string_for_static_func_timeval_to_seconds,
    string_for_step1_horizon_schedule,
    string_for_step2_validate_multigrid_inputs,
    string_for_step3_initialize_bhahaha_data_structs_and_solver_params,
    string_for_step4_populate_iteration_metadata,
    string_for_step5_apply_bbh_mode_logic,
    string_for_step6_apply_robustness_improv_and_extrap_horizon_guesses,
    string_for_step7a_to_c_main_loop_for_each_horizon,
    string_for_step7e_to_g_main_loop_for_each_horizon,
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

  // STEP 2: Return if no points to interpolate.
  if (total_interp_points == 0)
    return; // END IF: total_interp_points == 0, no points to interpolate

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
      * timing helper (timeval_to_seconds)
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

    prefunc += string_for_static_func_timeval_to_seconds()

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
It orchestrates initialization, BBH logic, extrapolation, spherical-target setup/staging, solving,
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
-       `BHaHAHA_compute_dst_pts_metric_data_nrpy` to prepare spherical target points and per-point BSSN storage
-       using the current `x_guess[h]`, `y_guess[h]`, `z_guess[h]`, and the populated `radii_interp`.
-       The BSSN→ADM packing occurs later via `BHaHAHA_transform_BSSN_to_ADM`.
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
- 8. Emits start/finish messages when verbosity is enabled (timing may be reported elsewhere).

@param commondata - Pointer to `commondata_struct` holding global BHaHAHA settings,
                    persistent horizon data, simulation state, and BBH mode flags.
                    This struct is extensively read from and modified.
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
    body += string_for_step1_horizon_schedule()

    body += string_for_step2_validate_multigrid_inputs()

    body += string_for_step3_initialize_bhahaha_data_structs_and_solver_params()

    body += r"""
  // Local arrays for per-horizon guesses for the current find.
  REAL r_min_guess[max_num_horizons], r_max_guess[max_num_horizons];
  """

    body += string_for_step4_populate_iteration_metadata()

    body += string_for_step5_apply_bbh_mode_logic()

    body += string_for_step6_apply_robustness_improv_and_extrap_horizon_guesses()

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

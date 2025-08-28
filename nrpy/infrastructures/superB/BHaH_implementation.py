"""
Register C++ function for finding horizons with superB BHaHAHA.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures.BHaH.BHaHAHA.BHaH_implementation import (
    build_bhahaha_prefunc,
    register_bhahaha_commondata_and_params,
    string_for_step1_horizon_schedule,
    string_for_step2_validate_multigrid_inputs,
    string_for_step3_initialize_bhahaha_data_structs_and_solver_params,
    string_for_step4_populate_iteration_metadata,
    string_for_step5_apply_bbh_mode_logic,
    string_for_step6_apply_robustness_improv_and_extrap_horizon_guesses,
    string_for_step7a_to_d_main_loop_for_each_horizon,
    string_for_step7e_to_g_main_loop_for_each_horizon,
    string_for_step8_print_total_elapsed_time,
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

    prefunc += """
    static struct timeval start_time_total, iter_time_tracker;
    """

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
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata, REAL x_guess[], REAL y_guess[], REAL z_guess[],REAL **radii, int *total_elements, REAL (**dst_x0x1x2)[3], REAL ***dst_data_ptrs,const int which_horizon, const int which_bhahaha_part"""
    body = """
    switch (which_bhahaha_part) {
    case BHAHAHA_FIND_HORIZONS_SETUP: {
    """
    body += string_for_step1_horizon_schedule()

    body += r"""
  gettimeofday(&start_time_total, NULL);
  """

    body += string_for_step2_validate_multigrid_inputs()

    body += string_for_step3_initialize_bhahaha_data_structs_and_solver_params()

    body += r"""
  // Local arrays for per-horizon guesses for the current find.
  REAL x_guess[max_num_horizons], y_guess[max_num_horizons], z_guess[max_num_horizons];
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

    body += string_for_step7a_to_d_main_loop_for_each_horizon(single_horizon=True)

    body += r"""
    break;
    }
    case BHAHAHA_FIND_HORIZONS_FIND_AND_WRITE_TO_FILE: {
      int h = which_horizon;
      bhahaha_params_and_data_struct *current_horizon_params = &commondata->bhahaha_params_and_data[h];
  """

    body += string_for_step7e_to_g_main_loop_for_each_horizon(single_horizon=True)

    body += string_for_step8_print_total_elapsed_time()

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

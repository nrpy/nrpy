"""
Register and configure C functions for BHaHAHA diagnostics.
This includes interpolating 3D grid data to a 2D surface, computing various norms,
coordinate radii, and the proper area of the horizon. Diagnostic data is written
to output files for analysis.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures.BHaH import griddata_commondata


def register_CFunction_diagnostics() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for diagnostics in the simulation.

    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    DocTests:
    >>> env = register_CFunction_diagnostics()
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Add BHaHAHA.h's diagnostics struct to commondata, so we can easily read/write those variables here.
    griddata_commondata.register_griddata_commondata(
        __name__,
        "bhahaha_diagnostics_struct *restrict bhahaha_diagnostics",
        "diagnostics quantities; struct defined in BHaHAHA.h",
        is_commondata=True,
    )
    griddata_commondata.register_griddata_commondata(
        __name__,
        "int is_final_iteration",
        "diagnostics quantities; struct defined in BHaHAHA.h",
        is_commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
/**
 * Displays spin values based on provided circumference ratio comparisons.
 *
 * This helper function prints the spin component labels along with their corresponding
 * ratio-based values. If spin_from_ratios are not applicable (indicated by -10.0), it displays "N/A".
 *
 * @param spin_label The label identifying the spin component (e.g., "spin_x").
 * @param spin_from_ratio1 The first spin_from_ratio value used for calculating the spin component.
 * @param spin_from_ratio2 The second spin_from_ratio value used for calculating the spin component.
 * @param spin_from_ratio1_label The label describing the first spin_from_ratio.
 * @param spin_from_ratio2_label The label describing the second spin_from_ratio.
 */
static void display_spin(const char *spin_label, const double spin_from_ratio1, const double spin_from_ratio2, const char *spin_from_ratio1_label,
                         const char *spin_from_ratio2_label) {
  // Spin values of -10.0 indicate that a spin was not derivable from the circumference ratio.
  if (spin_from_ratio1 == -10.0 && spin_from_ratio2 == -10.0) {
    printf("#%s = (    N/A    ,     N/A    ) based on (%s, %s) proper circumference spin_from_ratios.\n", spin_label, spin_from_ratio1_label,
           spin_from_ratio2_label);
  } else {
    printf("#%s = (%.4g, %.4g) based on (%s, %s) proper circumference spin_from_ratios.\n", spin_label,
           (spin_from_ratio1 != -10.0) ? spin_from_ratio1 : 0.0, (spin_from_ratio2 != -10.0) ? spin_from_ratio2 : 0.0, spin_from_ratio1_label,
           spin_from_ratio2_label);
  } // END IF: Spin derivable from ratios
} // END FUNCTION: display_spin
"""
    desc = """Performs apparent horizon diagnostics for BHaHAHA.

This function conducts various diagnostic checks and computations related to
the apparent horizon in the BHaHAHA simulation. Diagnostics are performed
at specified iteration intervals and may include interpolation, centroid
calculations, norm evaluations, and detailed final iteration analyses.

@param commondata Pointer to the common data structure containing simulation parameters and state.
@param griddata Pointer to the grid data structure containing grid-related parameters and functions.

@note This function updates the error_flag within commondata based on diagnostic outcomes."""
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    # Removed the unused variable 'Th'
    _ = gri.register_gridfunctions(
        "Theta", group="AUX", gf_array_name="diagnostic_output_gfs"
    )
    _ = par.register_CodeParameter(
        "int", __name__, "output_diagnostics_every_nn", 500, commondata=True
    )
    body = r"""
  // Check if diagnostics should be performed at the current iteration.
  if (commondata->nn % commondata->output_diagnostics_every_nn == 0) {
    {
      // Retrieve grid sizes including ghost zones for each dimension.
      const int Nxx_plus_2NGHOSTS0 = griddata[0].params.Nxx_plus_2NGHOSTS0;
      const int Nxx_plus_2NGHOSTS1 = griddata[0].params.Nxx_plus_2NGHOSTS1;
      const int Nxx_plus_2NGHOSTS2 = griddata[0].params.Nxx_plus_2NGHOSTS2;

      // Perform interpolation on the source grid using radial spokes.
      commondata->error_flag = bah_interpolation_1d_radial_spokes_on_3d_src_grid(
          &griddata[0].params, commondata, &griddata[0].gridfuncs.y_n_gfs[IDX4(HHGF, 0, 0, 0)], griddata[0].gridfuncs.auxevol_gfs);

      // Exit diagnostics if interpolation fails.
      if (commondata->error_flag != BHAHAHA_SUCCESS)
        return;
    } // END BLOCK: Interpolation and grid size setup

    // Calculate area centroid and theta norms for the apparent horizon.
    bah_diagnostics_area_centroid_and_Theta_norms(commondata, griddata);

    // Check if verbosity level is set to display detailed diagnostics.
    if (commondata->bhahaha_params_and_data->verbosity_level == 2) {

      // Print diagnostic headers during the first diagnostic iteration.
      if (commondata->nn == 0) {
        printf("#*** Horizon %d / %d : (Nr x Ntheta x Nphi) = (%d x %d x %d)\n#*** Tolerances: (Linf_Theta*M, L2_Theta*M) = (%e, %e)\n",
               commondata->bhahaha_params_and_data->which_horizon, commondata->bhahaha_params_and_data->num_horizons, commondata->external_input_Nxx0,
               griddata[0].params.Nxx1, griddata[0].params.Nxx2, commondata->bhahaha_params_and_data->Theta_Linf_times_M_tolerance,
               commondata->bhahaha_params_and_data->Theta_L2_times_M_tolerance);
        printf("#Iter |min_r |max_r |maxsrch_r|Linf_Theta|L2_Theta|      Area      |   M_irr   |Nth|Nph|N_Theta_eval|\n");
        printf("#-----|------|------|---------|----------|--------|----------------|-----------|---|---|------------|\n");
      } // END IF nn == 0

      // Access the diagnostics structure for current iteration values.
      bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

      // r_max_interior = r_min_external_input + ((Nr_external_input-BHAHAHA_NGHOSTS) + 0.5) * dr
      const REAL r_max_interior =
          commondata->bhahaha_params_and_data->r_min_external_input +
          ((commondata->bhahaha_params_and_data->Nr_external_input - BHAHAHA_NGHOSTS) + 0.5) * commondata->bhahaha_params_and_data->dr_external_input;

      // Display current diagnostic metrics.
      printf("%6d %6.4f %6.4f   %6.4f  %8.4e %6.2e %14.10e %7.5e %3d %3d %12ld\n", commondata->nn, commondata->min_radius_wrt_grid_center,
             commondata->max_radius_wrt_grid_center, r_max_interior, bhahaha_diags->Theta_Linf_times_M, bhahaha_diags->Theta_L2_times_M,
             bhahaha_diags->area, sqrt(bhahaha_diags->area / (16 * M_PI)), griddata[0].params.Nxx1, griddata[0].params.Nxx2,
             bhahaha_diags->Theta_eval_points_counter);
    } // END IF verbosity level == 2

    // Verify that the minimum coordinate radius meets the required threshold.
    if (commondata->min_radius_wrt_grid_center < 3.0 * commondata->external_input_dxx0) {
      commondata->error_flag = FIND_HORIZON_HORIZON_TOO_SMALL;
      return;
    }

    // Perform additional diagnostics during the final iteration.
    if (commondata->is_final_iteration) {
      // Compute minimum, maximum, and mean radii with respect to the centroid.
      bah_diagnostics_min_max_mean_radii_wrt_centroid(commondata, griddata);

      // Calculate proper circumferences and update the error flag based on the outcome.
      commondata->error_flag = bah_diagnostics_proper_circumferences(commondata, griddata);
      if (commondata->error_flag != BHAHAHA_SUCCESS)
        return;

      // Display detailed final iteration diagnostics if verbosity is enabled.
      {
        bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

        if (commondata->bhahaha_params_and_data->verbosity_level > 0) {
          printf("#-={ BHaHAHA Horizon %d / %d Diagnostics, t=%.3f }=-\n", commondata->bhahaha_params_and_data->which_horizon,
                 commondata->bhahaha_params_and_data->num_horizons, commondata->bhahaha_params_and_data->time_external_input);
          printf("#Final iteration: %d, at (Nr x Ntheta x Nphi) = (%d x %d x %d)\n", commondata->nn, commondata->external_input_Nxx0,
                 griddata[0].params.Nxx1, griddata[0].params.Nxx2);
          printf("#(%5.5e, %5.5e) = (L_inf Theta, L2 Theta)\n", bhahaha_diags->Theta_Linf_times_M, bhahaha_diags->Theta_L2_times_M);
          printf("#(%5.5e, %5.5e) = (Area, M_irr)\n", bhahaha_diags->area, sqrt(bhahaha_diags->area / (16 * M_PI)));
          printf("#(%+4.4e, %+4.4e, %+4.4e) = (x, y, z) centroid, wrt input grid origin\n", bhahaha_diags->x_centroid_wrt_coord_origin,
                 bhahaha_diags->y_centroid_wrt_coord_origin, bhahaha_diags->z_centroid_wrt_coord_origin);
        } // END IF verbosity_level > 0, then print the diagnostics

        REAL x_center, y_center, z_center, r_min, r_max;
        bah_xyz_center_r_minmax(commondata->bhahaha_params_and_data, &x_center, &y_center, &z_center, &r_min, &r_max);

        // Cycle data from previous horizon finds. _m1 data depend on centroids being computed in bah_diagnostics()
        commondata->bhahaha_params_and_data->t_m3 = commondata->bhahaha_params_and_data->t_m2;
        commondata->bhahaha_params_and_data->t_m2 = commondata->bhahaha_params_and_data->t_m1;
        commondata->bhahaha_params_and_data->t_m1 = commondata->bhahaha_params_and_data->time_external_input;
        commondata->bhahaha_params_and_data->x_center_m3 = commondata->bhahaha_params_and_data->x_center_m2;
        commondata->bhahaha_params_and_data->y_center_m3 = commondata->bhahaha_params_and_data->y_center_m2;
        commondata->bhahaha_params_and_data->z_center_m3 = commondata->bhahaha_params_and_data->z_center_m2;
        commondata->bhahaha_params_and_data->x_center_m2 = commondata->bhahaha_params_and_data->x_center_m1;
        commondata->bhahaha_params_and_data->y_center_m2 = commondata->bhahaha_params_and_data->y_center_m1;
        commondata->bhahaha_params_and_data->z_center_m2 = commondata->bhahaha_params_and_data->z_center_m1;
        commondata->bhahaha_params_and_data->x_center_m1 = x_center + bhahaha_diags->x_centroid_wrt_coord_origin;
        commondata->bhahaha_params_and_data->y_center_m1 = y_center + bhahaha_diags->y_centroid_wrt_coord_origin;
        commondata->bhahaha_params_and_data->z_center_m1 = z_center + bhahaha_diags->z_centroid_wrt_coord_origin;
        commondata->bhahaha_params_and_data->r_min_m3 = commondata->bhahaha_params_and_data->r_min_m2;
        commondata->bhahaha_params_and_data->r_max_m3 = commondata->bhahaha_params_and_data->r_max_m2;
        commondata->bhahaha_params_and_data->r_min_m2 = commondata->bhahaha_params_and_data->r_min_m1;
        commondata->bhahaha_params_and_data->r_max_m2 = commondata->bhahaha_params_and_data->r_max_m1;
        commondata->bhahaha_params_and_data->r_min_m1 = commondata->bhahaha_diagnostics->min_coord_radius_wrt_centroid;
        commondata->bhahaha_params_and_data->r_max_m1 = commondata->bhahaha_diagnostics->max_coord_radius_wrt_centroid;

        if (commondata->bhahaha_params_and_data->verbosity_level > 0)
          printf("#(%+4.4e, %+4.4e, %+4.4e) = (x, y, z) centroid, wrt global origin\n", commondata->bhahaha_params_and_data->x_center_m1,
                 commondata->bhahaha_params_and_data->y_center_m1, commondata->bhahaha_params_and_data->z_center_m1);

        if (commondata->bhahaha_params_and_data->verbosity_level > 0) {
          printf("#(%5.5e, %5.5e, %5.5e) = (min, max, mean) coord radii, relative to centroid\n", bhahaha_diags->min_coord_radius_wrt_centroid,
                 bhahaha_diags->max_coord_radius_wrt_centroid, bhahaha_diags->mean_coord_radius_wrt_centroid);

          printf("#(%5.5e, %5.5e, %5.5e) = (xy, xz, yz) proper circumferences\n", bhahaha_diags->xy_plane_circumference,
                 bhahaha_diags->xz_plane_circumference, bhahaha_diags->yz_plane_circumference);

          // Display spin_x based on (xy/yz, xz/yz) ratios
          display_spin("spin_x", bhahaha_diags->spin_a_x_from_xy_over_yz_prop_circumfs, bhahaha_diags->spin_a_x_from_xz_over_yz_prop_circumfs, //
                       "xy/yz", "xz/yz");

          // Display spin_y based on (xy/xz, yz/xz) ratios
          display_spin("spin_y", bhahaha_diags->spin_a_y_from_xy_over_xz_prop_circumfs, bhahaha_diags->spin_a_y_from_yz_over_xz_prop_circumfs, //
                       "xy/xz", "yz/xz");

          // Display spin_z based on (xz/xy, yz/xy) ratios
          display_spin("spin_z", bhahaha_diags->spin_a_z_from_xz_over_xy_prop_circumfs, bhahaha_diags->spin_a_z_from_yz_over_xy_prop_circumfs, //
                       "xz/xy", "yz/xy");
        } // END IF verbosity level > 0
      } // END compute, store, and (optionally) print final diagnostics
    } // END IF final iteration
  } // END IF output diagnostics
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


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        raise RuntimeError(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)"
        )
    print(f"Doctest passed: All {results.attempted} test(s) passed")

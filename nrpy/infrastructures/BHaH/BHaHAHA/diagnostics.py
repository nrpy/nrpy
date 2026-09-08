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
from nrpy.infrastructures import BHaH


def register_CFunction_diagnostics() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for diagnostics in the simulation.

    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Add BHaHAHA.h's diagnostics struct to commondata, so we can easily read/write those variables here.
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "bhahaha_diagnostics_struct *restrict bhahaha_diagnostics",
        "diagnostics quantities; struct defined in BHaHAHA.h",
        is_commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
"""
    desc = """Performs apparent horizon diagnostics for BHaHAHA.

This function conducts various diagnostic checks and computations related to
the apparent horizon in the BHaHAHA simulation. Diagnostics are performed
at specified iteration intervals and may include interpolation, centroid
calculations, norm evaluations.

@param[in,out] commondata Pointer to the common data structure containing simulation parameters and state.
@param[in,out] griddata Pointer to the grid data structure containing grid-related parameters and functions.

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
    cfunc_decorators = r"""
#ifdef __CUDACC__
__device__
#endif
"""
    body = r"""
  #ifdef __CUDACC__
  // Set up cooperative gropu 
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();

  // Set up shared memory
  extern __shared__ REAL s[]; 
  #endif
  // Check if diagnostics should be performed at the current iteration.
  if (commondata->nn % commondata->output_diagnostics_every_nn == 0) {
    {
      // Retrieve grid sizes including ghost zones for each dimension.
      const int Nxx_plus_2NGHOSTS0 = griddata[0].params.Nxx_plus_2NGHOSTS0;
      const int Nxx_plus_2NGHOSTS1 = griddata[0].params.Nxx_plus_2NGHOSTS1;
      const int Nxx_plus_2NGHOSTS2 = griddata[0].params.Nxx_plus_2NGHOSTS2;

      // Perform interpolation on the source grid using radial spokes.
      bah_interpolation_1d_radial_spokes_on_3d_src_grid(&griddata[0].params, commondata, &griddata[0].gridfuncs.y_n_gfs[IDX4pt(HHGF, 0)], griddata[0].gridfuncs.auxevol_gfs);
#ifdef __CUDACC__
      gpu_grid.sync();
#endif
      // Exit diagnostics if interpolation fails.
      if (commondata->error_flag != BHAHAHA_SUCCESS)
        return;
    } // END BLOCK: Interpolation and grid size setup

    // Calculate area centroid and theta norms for the apparent horizon.
    bah_diagnostics_area_centroid_and_Theta_norms(commondata, griddata);
#ifdef __CUDACC__
    gpu_grid.sync();
#endif

#ifdef __CUDACC__
    CUDA_ONE_THREAD(gpu_grid) {
#endif
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
#ifdef __CUDACC__
    } END_CUDA_ONE_THREAD //End verbose==2 printing
    gpu_grid.sync();
#endif
    // Verify that the minimum coordinate radius meets the required threshold.
    if (commondata->min_radius_wrt_grid_center < 3.0 * commondata->external_input_dxx0) {
      commondata->error_flag = FIND_HORIZON_HORIZON_TOO_SMALL;
      return;
    }
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
        cfunc_decorators=cfunc_decorators,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

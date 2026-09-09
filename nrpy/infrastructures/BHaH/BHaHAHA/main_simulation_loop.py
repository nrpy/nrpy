"""

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""


import nrpy.c_function as cfc
from nrpy.infrastructures import BHaH


def register_CFunction_main_simulation_loop() -> None:
    """
    Register the C function for the main simulation loop in find_horizon.

    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = r"""
"""
    desc = """Runs the main relaxation while loop for each resolution. This function:
1. Interpolates metric to current best sruface
2. Update the timestep based on the CFL condition.
3. Attempt over-relaxation every once in a while. If an over-relaxation is 
performed, the CFL timestep and metric will be updated to the new surface.
4. Time-varying eta prescription -- reduce eta with residual
5. Output diagnostic information.  
6. Determine if stop conditions are met to exit the simulation loop.
7. Advance the simulation using the Method of Lines with Runge-Kutta-like 
integration.
                                                                                                        
@params commondata - Pointer to the common data structure containing simulation parameters and state.
@params griddata - Pointer to the grid data structure containing grid-related parameters and functions.
@params resolution - int labeling the resolution from coarse(0) to fine(n:default=3);
                                                                                                        
@note this function updates the error_flag within commondata based on simulation outcomes. 
@returns void"""
    cfunc_type = "void"
    name = "main_simulation_loop"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata, int resolution"
    )
    cfunc_decorators = r"""
#ifdef __CUDACC__
__global__
#endif
"""
    body = r"""
#ifdef __CUDACC__
  //Set up cooperative group
  namespace cg = cooperative_groups;
  cg::grid_group gpu_grid = cg::this_grid();

  // Set up shared memory array
  extern __shared__ REAL s[];
#endif

  bhahaha_params_and_data_struct *bhahaha_params_and_data = commondata->bhahaha_params_and_data;
  bhahaha_diagnostics_struct *bhahaha_diags = commondata->bhahaha_diagnostics;
 
  const params_struct *restrict params = &griddata[0].params;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  int stop_condition = 0;
  while (commondata->time < commondata->t_final) { // Main loop to advance the simulation.
    // Step 1.a: Interpolate metric to current best surface. This is done at the start of each
    //           timestep instead of at each RK substep for efficiency reasons.
    bah_interpolation_1d_radial_spokes_on_3d_src_grid(
            &griddata[0].params, commondata, &griddata[0].gridfuncs.y_n_gfs[IDX4pt(HHGF, 0)], griddata[0].gridfuncs.auxevol_gfs);
#ifdef __CUDACC__
    gpu_grid.sync();
#endif
    if (commondata->error_flag != BHAHAHA_SUCCESS) {
      return;
    }

    // Step 1.b: Update the timestep based on the CFL condition.
    bah_cfl_limited_timestep_based_on_h_equals_r(commondata, griddata);

    // Reset Horizon
    if (commondata->nn == 0) {
      PARALLEL_LOOP(i0, NGHOSTS, NGHOSTS+1, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
        commondata->h_p[IDX3(i0,i1,i2)] = 0.0;
      } END_PARALLEL_LOOP 
      #ifdef __CUDACC__
      gpu_grid.sync();
      #endif
    } //End Reset horizon

    // Step 1.c: Attempt over-relaxation every once in a while. If an
    //           over-relaxation is performed, the CFL timestep
    //           and metric will be updated to the new surface.
    bah_over_relaxation(commondata, griddata);
    if (commondata->error_flag != BHAHAHA_SUCCESS) { // could fail due to metric interpolation.
      return;
    }

    // Step 1.d: Time-varying eta prescription -- reduce eta with residual
    if (bhahaha_params_and_data->enable_eta_varying_alg_for_precision_common_horizon && commondata->nn % 10000 == 0) {
        bah_interpolation_1d_radial_spokes_on_3d_src_grid(&griddata[0].params, commondata, &griddata[0].gridfuncs.y_n_gfs[IDX4pt(HHGF, 0)], griddata[0].gridfuncs.auxevol_gfs);
#ifdef __CUDACC__
      gpu_grid.sync();
#endif
      bah_diagnostics_area_centroid_and_Theta_norms(commondata, griddata);
#ifdef __CUDACC__
      CUDA_ONE_THREAD(gpu_grid) {
#endif
        REAL eta_min_times_M = 0.15;
        REAL eta_max_times_M = 3.0;
        if (resolution == 1)
          eta_max_times_M = 3.0;
        if (resolution >= 2)
          eta_max_times_M = 30.0;
  
        const REAL eta_damping_times_M = NRPYMAX(eta_min_times_M, eta_max_times_M * sqrt(bhahaha_diags->Theta_Linf_times_M));
        commondata->eta_damping = eta_damping_times_M / bhahaha_params_and_data->M_scale;
#ifdef __CUDACC__
      } END_CUDA_ONE_THREAD //End global variable modification
      gpu_grid.sync();
#endif
      PARALLEL_LOOP(i0, NGHOSTS, NGHOSTS+1, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
        griddata[0].gridfuncs.y_n_gfs[IDX4(VVGF, i0, i1, i2)] =
            commondata->eta_damping * griddata[0].gridfuncs.y_n_gfs[IDX4(HHGF, i0, i1, i2)];
      } END_PARALLEL_LOOP // END LOOP over all gridpoints on horizon surface.
      #ifdef __CUDACC__
      gpu_grid.sync();
      #endif
    } // END time-varying eta prescription.

    // Step 1.e: Output diagnostic information.  
    bah_diagnostics(commondata, griddata);
#ifdef __CUDACC__
    gpu_grid.sync();
#endif
    if (commondata->error_flag != BHAHAHA_SUCCESS) {
      return;
    } // END IF: Check for diagnostic errors

    // Step 1.f: Determine if stop conditions are met to exit the simulation loop.
    if (commondata->nn > bhahaha_params_and_data->max_iterations) {
      commondata->error_flag = FIND_HORIZON_MAX_ITERATIONS_EXCEEDED;
      stop_condition = 1;
      return;
    } else if (commondata->nn > 2 && // Ensure a minimum number of iterations.
               bhahaha_diags->Theta_Linf_times_M <= bhahaha_params_and_data->Theta_Linf_times_M_tolerance &&
               bhahaha_diags->Theta_L2_times_M <= bhahaha_params_and_data->Theta_L2_times_M_tolerance) {
      stop_condition = 1;
      return;
    } // END IF: Check multiple stop conditions

    // Step 1.g: Advance the simulation using the Method of Lines with Runge-Kutta-like integration.
    if (!stop_condition) {
      bah_MoL_step_forward_in_time(commondata, griddata);
#ifdef __CUDACC__
      gpu_grid.sync();
#endif
    }
    if (commondata->error_flag != BHAHAHA_SUCCESS) {
      return;
    } // END IF: Check for time-stepping errors
#ifdef __CUDACC__
    gpu_grid.sync();
#endif
  } // END LOOP: Main simulation loop
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

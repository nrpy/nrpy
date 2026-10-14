"""
Driver function for finding apparent horizons.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc
from nrpy.infrastructures import BHaH


def register_CFunction_find_horizon() -> None:
    """Register main driver (find horizon) function for BHaHAHA, based on BlackHoles@Home's main_c.py."""
    # Add BHaHAHA.h include near the top of BHaH_defines.h ("general" module):
    BHaH.BHaH_defines_h.register_BHaH_defines(
        "after_general", """#include "BHaHAHA.h"\n"""
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <sys/time.h> // Include sys/time.h for timing functions like gettimeofday()

/**
 * Converts two timeval structures to milliseconds and returns the elapsed time.
 *
 * @param start - The starting time.
 * @param end - The ending time.
 * @return - The elapsed time in milliseconds.
 */
static REAL timeval_to_milliseconds(struct timeval start, struct timeval end) {
  double start_ms = start.tv_sec * 1000.0 + start.tv_usec / 1000.0;
  double end_ms = end.tv_sec * 1000.0 + end.tv_usec / 1000.0;
  return end_ms - start_ms;
}

/**
 * Frees all dynamically allocated memory associated with griddata,
 * except for external input grid functions.
 *
 * @param commondata - Pointer to the common data structure containing shared parameters.
 * @param griddata - Pointer to the grid data structure to be freed.
 */
static void free_all_but_external_input_gfs(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
  const int grid = 0;

  // Free precomputed reference metric arrays.
  bah_rfm_precompute_free(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
  FREE(griddata[grid].rfmstruct);

  // Free inner boundary condition array.
  FREE(griddata[grid].bcstruct.inner_bc_array);

  // Free pure outer boundary condition arrays.
  for (int ng = 0; ng < NGHOSTS * 3; ng++) {
    FREE(griddata[grid].bcstruct.pure_outer_bc_array[ng]);
  } // END LOOP: freeing pure outer boundary condition arrays

  // Free y_n_gfs, intermediate-stage gfs, and auxevol_gfs, needed by MoL.
  #ifdef __CUDACC__
  cudaFree(griddata[grid].gridfuncs.y_n_gfs);
  #else
  BHAH_FREE(griddata[grid].gridfuncs.y_n_gfs);
  #endif
  bah_MoL_free_intermediate_stage_gfs(&griddata[grid].gridfuncs);
  #ifdef __CUDACC__
  cudaFree(griddata[grid].gridfuncs.auxevol_gfs);
  #else
  BHAH_FREE(griddata[grid].gridfuncs.auxevol_gfs);
  #endif

  // Free coordinate arrays for each dimension.
  for (int i = 0; i < 3; i++) {
    FREE(griddata[grid].xx[i]);
  } // END LOOP: freeing coordinate arrays

  // Free the griddata structure itself.
  free(griddata);

  // Free interpolation source grid functions.
  FREE(commondata->interp_src_gfs);

  // Free interpolation source coordinate arrays for each dimension.
  for (int i = 0; i < 3; i++) {
    FREE(commondata->interp_src_r_theta_phi[i]);
  } // END LOOP: freeing interpolation source coordinate arrays

  // Free previous horizon guess array, used for overstep.
  if (commondata->h_p != NULL)
    FREE(commondata->h_p);
} // END FUNCTION: free_all_but_external_input_gfs
"""
    desc = """
Finds the apparent horizon using BHaHAHA.

This driver function initializes necessary data structures, sets up grids, and runs the main simulation loop
to identify the apparent horizon with progressively refined grid resolutions.

@param[in,out] bhahaha_params_and_data Input parameters and data for the algorithm.
@param[out] bhahaha_diags Diagnostics data structure to be updated during execution.
@return Returns BHaHAHA (0) on success or a nonzero error code on failure.
 """
    cfunc_type = "int"
    name = "find_horizon"
    params = "bhahaha_params_and_data_struct *restrict bhahaha_params_and_data, bhahaha_diagnostics_struct *restrict bhahaha_diags"
    body = r"""
  // Step 1.a: Start global timer.
  struct timeval start_time;
  {
    // Verify that gettimeofday() is functional before proceeding.
    if (gettimeofday(&start_time, NULL) != 0) {
      return FIND_HORIZON_GETTIMEOFDAY_BROKEN;
    }
  } // END BLOCK: gettimeofday() sanity check.

  commondata_struct commondata; // Structure containing parameters common to all grids.


  // Step 1.b: Initialize commondata parameters to their default values.
  bah_commondata_struct_set_to_default(&commondata);
  // Assign input diagnostics and parameters after defaults are set so zero-initialization
  // does not clobber these caller-owned pointers. 
  commondata.bhahaha_diagnostics = bhahaha_diags;
  commondata.bhahaha_params_and_data = bhahaha_params_and_data;
  commondata.eta_damping = bhahaha_params_and_data->eta_damping_times_M / bhahaha_params_and_data->M_scale;
  commondata.CFL_FACTOR = bhahaha_params_and_data->cfl_factor;
  commondata.KO_diss_strength = commondata.bhahaha_params_and_data->KO_strength;
  // Initialize counter for total number of points at which Theta is evaluated.
  commondata.bhahaha_diagnostics->Theta_eval_points_counter = 0;

  // Set internal flags for horizon refinement and final iteration diagnostics.
  commondata.use_coarse_horizon = 0;
  commondata.is_final_iteration = 0;


  // Step 1.c: Define angular resolution patterns to optimize the solving process.
  const int n_resolutions = bhahaha_params_and_data->num_resolutions_multigrid;
  int Ntheta[MAX_RESOLUTIONS], Nphi[MAX_RESOLUTIONS];
  memcpy(Ntheta, bhahaha_params_and_data->Ntheta_array_multigrid, sizeof(int) * MAX_RESOLUTIONS);
  memcpy(Nphi, bhahaha_params_and_data->Nphi_array_multigrid, sizeof(int) * MAX_RESOLUTIONS);

  // Step 1.d: Set up external input grids by adding inner ghost zones and applying boundary conditions.
  // A deep copy of bhahaha_params_and_data and bhahaha_diagnostics occurs within the external set up.

#ifdef __CUDACC__
  // Variables for deep copies used after GPU kernels
  bhahaha_params_and_data_struct *h_bhahaha_params_and_data = commondata.bhahaha_params_and_data;
  bhahaha_diagnostics_struct *h_bhahaha_diagnostics = commondata.bhahaha_diagnostics;
#endif

  commondata.external_input_gfs_Cart_basis_no_gzs = bhahaha_params_and_data->input_metric_data;
  bah_numgrid__external_input_set_up(&commondata, n_resolutions, Ntheta, Nphi);
  if (commondata.error_flag != BHAHAHA_SUCCESS) {
    return commondata.error_flag;
  }
  
#ifdef __CUDACC__
  //Variables for Deep copy bhahaha_params_and_data & bhahaha_diagnostics
  bhahaha_params_and_data_struct *d_bhahaha_params_and_data = commondata.bhahaha_params_and_data;
  bhahaha_diagnostics_struct *d_bhahaha_diagnostics = commondata.bhahaha_diagnostics;

  commondata_struct *d_commondata = NULL;
  cudaMalloc((void**)&d_commondata, sizeof(commondata_struct));
#endif
  
  // Step 2: Iterate over different grid resolutions to refine the apparent horizon.
  commondata.error_flag = BHAHAHA_SUCCESS; // Assume success initially.

  for (int resolution = 0; resolution < n_resolutions; resolution++) {
    // Start timing for the current resolution.
    struct timeval res_start_time;
    gettimeofday(&res_start_time, NULL);

    // Adjust diagnostics output frequency based on grid resolution.
    commondata.output_diagnostics_every_nn = 4 * (int)pow(Ntheta[n_resolutions - 1] / Ntheta[resolution], 2.0);

    int Nx_evol_grid[3];
    griddata_struct *restrict griddata; // Structure containing data specific to the current grid.
    const int grid = 0;

    // Define the evolution grid size for the current resolution.
    Nx_evol_grid[0] = 1;
    Nx_evol_grid[1] = Ntheta[resolution];
    Nx_evol_grid[2] = Nphi[resolution];

    // Step 2.b: Set up interpolation source grid and allocate grid functions.
    bah_numgrid__interp_src_set_up(&commondata, Nx_evol_grid);
    if (commondata.error_flag != BHAHAHA_SUCCESS) {
      for (int i = 0; i < 3; i++) {
        FREE(commondata.external_input_r_theta_phi[i]);
      }
      FREE(commondata.external_input_gfs);
      return commondata.error_flag;
    }

    // Step 2.c: Allocate memory for MAXNUMGRIDS griddata structures.
    griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * MAXNUMGRIDS);
#ifdef __CUDACC__
    // Create a device side copy of griddata
    griddata_struct *d_griddata = NULL;
    cudaMalloc((void**)&d_griddata, sizeof(griddata_struct) * MAXNUMGRIDS);
#endif

    // Step 2.d: Initialize griddata parameters to their default values.
    bah_params_struct_set_to_default(&commondata, griddata);

    // Step 2.e: Configure the 2D numerical grid for the Apparent Horizon finder.
    bah_numgrid__evol_set_up(&commondata, griddata, Nx_evol_grid);
    if (commondata.error_flag != BHAHAHA_SUCCESS) {
      for (int i = 0; i < 3; i++) {
        FREE(commondata.external_input_r_theta_phi[i]);
      }
      FREE(commondata.external_input_gfs);
      return commondata.error_flag;
    }
    const params_struct *restrict params = &griddata[grid].params;
    const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
    const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
    const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

    {
      const int grid = 0;

      // Step 3.a: Allocate storage for initial grid functions (y_n_gfs).
#ifdef __CUDACC__
      cudaMalloc((void**)&griddata[grid].gridfuncs.y_n_gfs, NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));
#else
      BHAH_MALLOC(griddata[grid].gridfuncs.y_n_gfs, NUM_EVOL_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));
#endif

      // Step 3.b: Allocate storage for additional grid functions required for time-stepping.
      bah_MoL_malloc_intermediate_stage_gfs(&commondata, params, &griddata[grid].gridfuncs);
#ifdef __CUDACC__
      gpuErrchk( cudaMalloc((void**)&griddata[grid].gridfuncs.auxevol_gfs, NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL)) );
#else
      BHAH_MALLOC(griddata[grid].gridfuncs.auxevol_gfs, NUM_AUXEVOL_GFS * Nxx_plus_2NGHOSTS_tot * sizeof(REAL));
#endif

      // Step 3.c: Initialize commondata.h_p = NULL, so that if Step 5.a (interp 1D) fails,
      //   it doesn't trigger a double free() of h_p.
      commondata.h_p = NULL;
    } // END BLOCK: Allocation of grid functions

    // Step 4: Initialize initial data for the simulation.
    bah_initial_data(&commondata, griddata);
    if (commondata.error_flag != BHAHAHA_SUCCESS) {
      free_all_but_external_input_gfs(&commondata, griddata);
      for (int i = 0; i < 3; i++) {
        FREE(commondata.external_input_r_theta_phi[i]);
      }
      FREE(commondata.external_input_gfs);
      //return INITIAL_DATA_MALLOC_ERROR;
      return commondata.error_flag;
    }

    // Step 5: Execute the main simulation loop to evolve the horizon over time.
#ifdef __CUDACC__
    const int INTERP_ORDER = (2 * NinterpGHOSTS + 1); // Interpolation order corresponds to the number of points in the stencil per dimension.
#endif
    
#ifdef __CUDACC__
    //Allocate h_p
    cudaMalloc((void**)&commondata.h_p,sizeof(REAL) * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2);
#else
    commondata.h_p = (double *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2);
#endif

#ifdef __CUDACC__
    // Allocate struct for intermediate gpu norm calculations
    gpuErrchk( cudaMalloc((void**)&commondata.norms, sizeof(diag_norms_struct)) );
  
    // Update device side copies of commondata and griddata
    cudaMemcpy(d_commondata, &commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice);
    cudaMemcpy(d_griddata, griddata, sizeof(griddata_struct), cudaMemcpyHostToDevice);
#endif

    //Launcing cooperative kernel for main simulation loop
    {
#ifdef __CUDACC__
      void *Args[] ={ &d_commondata, &d_griddata, &resolution};
      int sharesize = 8*THREADSPERBLOCK*sizeof(REAL);
      COOPERATIVE_KERNEL(bah_main_simulation_loop, Args, sharesize);
#else
      bah_main_simulation_loop(&commondata, griddata, resolution);
#endif
    }

#ifdef __CUDACC__
    // Retrieve commondata & griddata from device
    gpuErrchk( cudaMemcpy(&commondata, d_commondata, sizeof(commondata_struct), cudaMemcpyDeviceToHost) );

    gpuErrchk( cudaMemcpy(d_griddata, griddata, sizeof(griddata_struct), cudaMemcpyHostToDevice) );
#endif

    {
      // End timing for the current resolution and display elapsed time.
      struct timeval end_time;
      gettimeofday(&end_time, NULL);

      if (bhahaha_params_and_data->verbosity_level == 2) {
        printf("#Nth x Nph = %d x %d elapsed time = %.1f ms / %.1f ms so far...\n", params->Nxx1, params->Nxx2,
               timeval_to_milliseconds(res_start_time, end_time), timeval_to_milliseconds(start_time, end_time));
      }
    } // END BLOCK: Timing and logging

    int compute_proper_circumferences= 1; 
    // Step 6: Compute final diagnostics if Horizon is found.
    if (commondata.error_flag == BHAHAHA_SUCCESS) {
  
      // Adjust setting for the final iteration, to trigger diagnostics and compute additional diagnostics.
      if (resolution >= n_resolutions - 1)
        commondata.is_final_iteration = 1;

      const int orig_output_diagnostics_every_nn = commondata.output_diagnostics_every_nn;
      commondata.output_diagnostics_every_nn = 1;
      // Compute diagnostics, cycling _m{3,2} and storing _m1 data while we're at it for {x,y,z}_center
      //   as they depend on centroids being computed, and r_{min,max} for good measure.

      // Allocate arrays needed for proper circumference diagnostics only on the final iteration
      if (commondata.is_final_iteration) {
        int NUM_DIAG_GFS = 5;
        int N_angle = griddata[grid].params.Nxx2;
#ifdef __CUDACC__
        cudaMalloc((void**)&commondata.diagnostics_arrays.metric_data_gfs, griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(REAL));
        cudaMalloc((void**)&commondata.diagnostics_arrays.dst_pts, sizeof(REAL) * N_angle*2);
        cudaMalloc((void**)&commondata.diagnostics_arrays.integrand, N_angle  * sizeof(REAL));
        cudaMalloc((void**)&commondata.diagnostics_arrays.theta, N_angle  * sizeof(REAL));
        cudaMalloc((void**)&commondata.diagnostics_arrays.phi, N_angle  * sizeof(REAL));
        // Update device side with diagnostics arrays allocations
        gpuErrchk( cudaMemcpy(d_commondata, &commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice) );
#else
        BHAH_MALLOC(commondata.diagnostics_arrays.metric_data_gfs, griddata[grid].params.Nxx_plus_2NGHOSTS0 * griddata[grid].params.Nxx_plus_2NGHOSTS1 * griddata[grid].params.Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(REAL));
        commondata.diagnostics_arrays.dst_pts = (REAL (*)[2])malloc(sizeof(REAL) * N_angle*2);
        commondata.diagnostics_arrays.integrand = (REAL*)malloc(N_angle  * sizeof(REAL));
        commondata.diagnostics_arrays.theta = (REAL*)malloc(N_angle  * sizeof(REAL));
        commondata.diagnostics_arrays.phi = (REAL*)malloc(N_angle  * sizeof(REAL));
#endif
        if (commondata.diagnostics_arrays.metric_data_gfs == NULL || commondata.diagnostics_arrays.dst_pts == NULL || commondata.diagnostics_arrays.integrand== NULL || commondata.diagnostics_arrays.theta== NULL || commondata.diagnostics_arrays.phi == NULL) {
          if (commondata.diagnostics_arrays.metric_data_gfs != NULL)
            FREE(commondata.diagnostics_arrays.metric_data_gfs);
          if (commondata.diagnostics_arrays.integrand!= NULL)
            FREE(commondata.diagnostics_arrays.integrand);
          if (commondata.diagnostics_arrays.dst_pts != NULL)
            FREE(commondata.diagnostics_arrays.dst_pts);
          if (commondata.diagnostics_arrays.theta!= NULL)
            FREE(commondata.diagnostics_arrays.theta);
          if (commondata.diagnostics_arrays.phi!= NULL)
            FREE(commondata.diagnostics_arrays.phi);
          compute_proper_circumferences= 0;
        }
      }

      {
#ifdef __CUDACC__
        void *Args[] = {&d_commondata, &d_griddata, &compute_proper_circumferences};
        int sharesize = sizeof(REAL)*(INTERP_ORDER + (THREADSPERBLOCK*4)*INTERP_ORDER + THREADSPERBLOCK*INTERP_ORDER*INTERP_ORDER) + sizeof(REAL)*INTERP_ORDER;
        COOPERATIVE_KERNEL(bah_diagnostics_full, Args, sharesize);
#else
        bah_diagnostics_full(&commondata, griddata, compute_proper_circumferences);
#endif
      }

#ifdef __CUDACC__
      // Free diagnostic arrays
      if (commondata.is_final_iteration && compute_proper_circumferences) {
        FREE(commondata.diagnostics_arrays.metric_data_gfs);
        FREE(commondata.diagnostics_arrays.dst_pts);
        FREE(commondata.diagnostics_arrays.integrand);
        FREE(commondata.diagnostics_arrays.theta);
        FREE(commondata.diagnostics_arrays.phi);
      }
#endif

      commondata.output_diagnostics_every_nn = orig_output_diagnostics_every_nn;
    } // END BLOCK: Compute final diagnostics 

    // Step 7: Save the coarse horizon for subsequent resolutions or output final diagnostics.
    if (commondata.error_flag == BHAHAHA_SUCCESS  && !commondata.is_final_iteration) {
      // Allocate memory for storing coarse horizon data.
      const int total_points = Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
#ifdef __CUDACC__
      cudaMalloc((void**)&commondata.coarse_horizon, sizeof(REAL) * total_points);
#else
      commondata.coarse_horizon = (double *)malloc(sizeof(REAL) * total_points);
#endif

      // Save grid parameters for the coarse horizon to maintain consistency.
      commondata.coarse_horizon_dxx1 = params->dxx1;
      commondata.coarse_horizon_dxx2 = params->dxx2;
      commondata.coarse_horizon_Nxx_plus_2NGHOSTS1 = Nxx_plus_2NGHOSTS1;
      commondata.coarse_horizon_Nxx_plus_2NGHOSTS2 = Nxx_plus_2NGHOSTS2;

      // Allocate coordinate arrays for the coarse horizon.
#ifdef __CUDACC__
      cudaMalloc((void**)&commondata.coarse_horizon_r_theta_phi[0], sizeof(REAL) * Nxx_plus_2NGHOSTS0);
      cudaMalloc((void**)&commondata.coarse_horizon_r_theta_phi[1], sizeof(REAL) * Nxx_plus_2NGHOSTS1);
      cudaMalloc((void**)&commondata.coarse_horizon_r_theta_phi[2], sizeof(REAL) * Nxx_plus_2NGHOSTS2);
#else
      commondata.coarse_horizon_r_theta_phi[0] = (double *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
      commondata.coarse_horizon_r_theta_phi[1] = (double *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
      commondata.coarse_horizon_r_theta_phi[2] = (double *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS2);
#endif

#ifdef __CUDACC__
      // Update device side copy of commondata and griddata
      gpuErrchk( cudaMemcpy(d_commondata, &commondata, sizeof(commondata_struct), cudaMemcpyHostToDevice) ); 
      gpuErrchk( cudaMemcpy(d_griddata, griddata, sizeof(griddata_struct), cudaMemcpyHostToDevice) ); 
#endif

#ifdef __CUDACC__
      bah_store_horizon<<<12,64>>>(d_commondata, d_griddata);
      cudaDeviceSynchronize();
#else
      bah_store_horizon(&commondata, griddata);
#endif
    } else { // IF: Horizon found at final resolution or horizon too large
#ifdef __CUDACC__
      //Deep retrieveal of commondata from device if final iteration or 1D interpolation was too large
      gpuErrchk( cudaMemcpy(&commondata, d_commondata, sizeof(commondata_struct), cudaMemcpyDeviceToHost) ); 
      gpuErrchk( cudaMemcpy(griddata, d_griddata, sizeof(griddata_struct), cudaMemcpyDeviceToHost) ); 

      gpuErrchk( cudaMemcpy(h_bhahaha_diagnostics, d_bhahaha_diagnostics, sizeof(bhahaha_diagnostics_struct), cudaMemcpyDeviceToHost) );
      commondata.bhahaha_diagnostics = h_bhahaha_diagnostics;
  
      gpuErrchk( cudaMemcpy(h_bhahaha_params_and_data, d_bhahaha_params_and_data, sizeof(bhahaha_params_and_data_struct), cudaMemcpyDeviceToHost) );
      commondata.bhahaha_params_and_data = h_bhahaha_params_and_data;

      REAL *d_y_n_gfs = griddata[0].gridfuncs.y_n_gfs;
      REAL *h_y_n_gfs = (REAL*)malloc(sizeof(REAL)*Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS);
      gpuErrchk( cudaMemcpy(h_y_n_gfs, d_y_n_gfs, sizeof(REAL)*Nxx_plus_2NGHOSTS_tot*NUM_EVOL_GFS,cudaMemcpyDeviceToHost) );
      griddata[0].gridfuncs.y_n_gfs = h_y_n_gfs;
#endif
      if (commondata.error_flag == BHAHAHA_SUCCESS  && commondata.is_final_iteration) {
        // Store the final horizon data and perform a last diagnostic output.
        const int NUM_THETA = params->Nxx1; // Required for IDX2() macro.
        memcpy(commondata.bhahaha_params_and_data->prev_horizon_m3, commondata.bhahaha_params_and_data->prev_horizon_m2,
               sizeof(REAL) * NUM_THETA * params->Nxx2);
        memcpy(commondata.bhahaha_params_and_data->prev_horizon_m2, commondata.bhahaha_params_and_data->prev_horizon_m1,
               sizeof(REAL) * NUM_THETA * params->Nxx2);
#pragma omp parallel for
        for (int i2 = 0; i2 < params->Nxx2; i2++) {
          for (int i1 = 0; i1 < params->Nxx1; i1++) {
            commondata.bhahaha_params_and_data->prev_horizon_m1[IDX2(i1, i2)] =
                griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, NGHOSTS, i1 + NGHOSTS, i2 + NGHOSTS)];
          } // END LOOP: theta indices
        } // END LOOP: phi indices
      } else if (commondata.error_flag == INTERP1D_HORIZON_TOO_LARGE) {
        // Handle specific error when the horizon exceeds interpolation limits.
        REAL max_radius = -1e10;
#pragma omp parallel for reduction(max : max_radius)
        for (int i2 = 0; i2 < params->Nxx2; i2++) {
          for (int i1 = 0; i1 < params->Nxx1; i1++) {
            REAL current_radius = griddata[grid].gridfuncs.y_n_gfs[IDX4(HHGF, NGHOSTS, i1 + NGHOSTS, i2 + NGHOSTS)];
            if (current_radius > max_radius) {
              max_radius = current_radius;
            }
          } // END LOOP: theta indices
        } // END LOOP: phi indices

        if (commondata.bhahaha_params_and_data->verbosity_level > 0) {
          // r_max_interior = r_min_external_input + ((Nr_external_input-BHAHAHA_NGHOSTS) + 0.5)*dr
          const REAL r_max_interior =
              commondata.bhahaha_params_and_data->r_min_external_input +
              ((commondata.bhahaha_params_and_data->Nr_external_input - BHAHAHA_NGHOSTS) + 0.5) * commondata.bhahaha_params_and_data->dr_external_input;
          printf("ERROR: h_max = %#.4g too close to r_max_search = %#.4g. "
                 "Try either increasing search radius or decreasing cfl_factor.\n",
                 max_radius, r_max_interior);
        }
      } // END IF: Handling specific error conditions
#ifdef __CUDACC__
      //Free host copy of y_n_gfs and reset device side
      free(h_y_n_gfs);
      griddata[grid].gridfuncs.y_n_gfs = d_y_n_gfs;
#endif
    } // END IF: Handling specific error conditions


    // Step 8: Release all allocated memory for the current grid resolution.
    free_all_but_external_input_gfs(&commondata, griddata);
#ifdef __CUDACC__ 
    //Free per iteration device side copies
    gpuErrchk( cudaFree(d_griddata) );
    gpuErrchk( cudaFree(commondata.norms) );
    //Free device side copies 
    if (commondata.is_final_iteration) {
      cudaFree(d_commondata);
      cudaFree(d_bhahaha_diagnostics);
      cudaFree(d_bhahaha_params_and_data);
    }
#endif

    // Report error stopping computation of proper circumferences
    if (compute_proper_circumferences) {
    } else {
      commondata.error_flag = DIAG_PROPER_CIRCUM_MALLOC_ERROR; 
    }

    if (commondata.error_flag != BHAHAHA_SUCCESS) {
      break;
    } // END IF: Check for errors after freeing memory
  } // END LOOP: Iterating over grid resolutions

  // Step 9: After processing all resolutions, release external input memory.
  for (int i = 0; i < 3; i++)
    FREE(commondata.external_input_r_theta_phi[i]);
  FREE(commondata.external_input_gfs);

  // Display final timing information if verbosity is enabled.
  if (commondata.bhahaha_params_and_data->verbosity_level > 0) {
    struct timeval end_time;
    gettimeofday(&end_time, NULL);
    printf("#-={ BHaHAHA finished horizon %d / %d in %#.4g seconds }=-\n", commondata.bhahaha_params_and_data->which_horizon,
           commondata.bhahaha_params_and_data->num_horizons, timeval_to_milliseconds(start_time, end_time) / 1000.0);
  }

  return commondata.error_flag;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )

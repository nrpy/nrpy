#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * -={ main() function }=-
 * Step 1.a: Initialize each CodeParameter in the commondata struct to its default value.
 * Step 1.b: Overwrite the default values with those from the parameter file.
 * Then overwrite the parameter file values with those provided via command line arguments.
 * Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
 * where each structure contains data specific to an individual grid.
 * Step 1.d: Initialize each CodeParameter in griddata.params to its default value.
 * Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
 * bcstruct, rfm_precompute, timestep, and others.
 * Step 2: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
 * Step 3: Set  initial data.
 * Step 4: Allocate storage for non-y_n gridfunctions needed for the Runge-Kutta-like time-stepping.
 * Step 5: MAIN SIMULATION LOOP
 * - Step 5.a: Functions to run prior to diagnostics. E.g., regridding.
 * - Step 5.b: Output diagnostics.
 * - Step 5.c: Prepare to step forward in time.
 * - Step 5.d: Step forward in time using Method of Lines with RK4 algorithm, applying  boundary conditions.
 * - Step 5.e: Finish up step in time.
 * Step 6: Free all allocated memory.
 */
int main(int argc, const char *argv[]) {
  commondata_struct commondata; // commondata contains parameters common to all grids.
  griddata_struct *griddata;    // griddata contains data specific to an individual grid.

  // Step 1.a: Initialize each CodeParameter in the commondata struct to its default value.
  commondata_struct_set_to_default(&commondata);

  // Step 1.b: Overwrite the default values with those from the parameter file.
  //           Then overwrite the parameter file values with those provided via command line arguments.
  cmdline_input_and_parfile_parser(&commondata, argc, argv);

  // Step 1.c: Allocate memory for MAXNUMGRIDS griddata structs,
  //           where each structure contains data specific to an individual grid.
  griddata = (griddata_struct *)malloc(sizeof(griddata_struct) * MAXNUMGRIDS);

  // Step 1.d: Initialize each CodeParameter in griddata.params to its default value.
  params_struct_set_to_default(&commondata, griddata);

  // Step 1.e: Set up numerical grids, including parameters such as NUMGRIDS, xx[3], masks, Nxx, dxx, invdxx,
  //           bcstruct, rfm_precompute, timestep, and others.
  {
    // If this function is being called for the first time, initialize commondata.time, nn, t_0, and nn_0 to 0.
    const bool calling_for_first_time = true;
    numerical_grids_and_timestep(&commondata, griddata, calling_for_first_time);
  }

  for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
    // Step 2: Allocate storage for the initial data (y_n_gfs gridfunctions) on each grid.
    MoL_malloc_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);
  }

  // Step 3: Set up initial data.
  initial_data(&commondata, griddata);

  // Step 4: Allocate storage for non-y_n gridfunctions, needed for the Runge-Kutta-like timestepping.
  for (int grid = 0; grid < commondata.NUMGRIDS; grid++)
    MoL_malloc_non_y_n_gfs(&commondata, &griddata[grid].params, &griddata[grid].gridfuncs);

  // Step 5: MAIN SIMULATION LOOP
  while (commondata.time < commondata.t_final) { // Main loop to progress forward in time.
    // Step 5.a: Main loop, part 1 (pre_diagnostics): Functions to run prior to diagnostics. E.g., regridding.
    // (nothing here; specify by setting pre_diagnostics string in register_CFunction_main_c().)

    // Step 5.b: Main loop, part 2: Output diagnostics
    diagnostics(&commondata, griddata);

    // Step 5.c: Main loop, part 3 (pre_MoL_step_forward_in_time): Prepare to step forward in time
    // (nothing here; specify by setting pre_MoL_step_forward_in_time string in register_CFunction_main_c().)

    // Step 5.d: Main loop, part 4: Step forward in time using Method of Lines with RK4 algorithm,
    //           applying  boundary conditions.
    MoL_step_forward_in_time(&commondata, griddata);

    // Step 5.e: Main loop, part 5 (post_MoL_step_forward_in_time): Finish up step in time
    // (nothing here; specify by setting post_MoL_step_forward_in_time string in register_CFunction_main_c().)

  } // End main loop to progress forward in time.
  BHAH_DEVICE_SYNC();
  // Step 6: Free all allocated memory
  {
    const bool free_non_y_n_gfs_and_core_griddata_pointers = true;
    griddata_free(&commondata, griddata, free_non_y_n_gfs_and_core_griddata_pointers);
  }
  return 0;
} // END FUNCTION main

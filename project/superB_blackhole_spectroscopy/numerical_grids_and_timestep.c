#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Set up a cell-centered grids of size grid_physical_size.
 */
void numerical_grids_and_timestep(commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time) {
  // Step 1.a: Set CoordSystem_hash
  CoordSystem_hash_setup(commondata, griddata);

  // Step 1.b: Set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx);
  }

  // Step 1.c: Allocate memory for and define reference-metric precomputation lookup tables

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    rfm_precompute_malloc(commondata, &griddata[grid].params, &griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, &griddata[grid].rfmstruct, griddata[grid].xx);
  }

  // Step 1.d: Set up curvilinear boundary condition struct (bcstruct)

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
  }

  // Step 1.e: Set timestep based on minimum spacing between neighboring gridpoints.
  commondata->dt = 1e30;
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
  }

  // Step 1.f: Initialize timestepping parameters to zero if this is the first time this function is called.
  if (calling_for_first_time) {
    commondata->nn = 0;
    commondata->nn_0 = 0;
    commondata->t_0 = 0.0;
    commondata->time = 0.0;
  }
}

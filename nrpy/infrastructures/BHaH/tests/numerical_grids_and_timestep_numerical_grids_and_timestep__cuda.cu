#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Set up numerical grids and timestep.
 */
void numerical_grids_and_timestep(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_host,
                                  bool calling_for_first_time) {

  // Step 1.a: Set each CodeParameter in griddata.params to default, for MAXNUMGRIDS grids.
  if (calling_for_first_time)
    params_struct_set_to_default(commondata, griddata);
  // Step 1.b: Set commondata->NUMGRIDS to number of CoordSystems we have
  commondata->NUMGRIDS = 1;

  {
    // Independent grids
    int Nx[3] = {-1, -1, -1};

    // Step 1.c: For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool set_xxmin_xxmax_to_defaults = true;
    int grid = 0;
    // In multipatch, gridname is a helpful alias indicating position of the patch. E.g., "lower Spherical patch"
    snprintf(griddata[grid].params.gridname, 100, "grid_Spherical");
    griddata[grid].params.CoordSystem_hash = SPHERICAL;
    griddata[grid].params.grid_physical_size = 10.0;
    numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, set_xxmin_xxmax_to_defaults);
    memcpy(&griddata_host[grid].params, &griddata[grid].params, sizeof(params_struct));
    grid++;
  }

  // Step 1.d: Allocate memory for and define reference-metric precomputation lookup tables
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    BHAH_MALLOC_DEVICE(griddata[grid].rfmstruct, sizeof(rfm_struct))
    griddata_host->rfmstruct = nullptr;

    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    rfm_precompute_malloc(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, griddata[grid].rfmstruct, griddata[grid].xx);
  }

  cpyDevicetoHost__grid(commondata, griddata_host, griddata);
  BHAH_DEVICE_SYNC();

  // Step 1.e: Set up curvilinear boundary condition struct (bcstruct)
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    bcstruct_set_up(commondata, &griddata[grid].params, griddata_host[grid].xx, &griddata_host[grid].bcstruct, &griddata[grid].bcstruct);
  }

  // Step 1.f: Set timestep based on minimum spacing between neighboring gridpoints.
  commondata->dt = 1e30;
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
  }
  // Step 1.g: Initialize timestepping parameters to zero if this is the first time this function is called.
  if (calling_for_first_time) {
    commondata->nn = 0;
    commondata->nn_0 = 0;
    commondata->t_0 = 0.0;
    commondata->time = 0.0;
  }

  // Step 1.h: Set grid_idx for each grid.
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    griddata[grid].params.grid_idx = grid;
  }
} // END FUNCTION numerical_grids_and_timestep

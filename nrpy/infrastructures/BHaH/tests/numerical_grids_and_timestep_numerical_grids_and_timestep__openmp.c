#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Set up numerical grids and timestep.
 */
void numerical_grids_and_timestep(commondata_struct *restrict commondata, griddata_struct *restrict griddata, bool calling_for_first_time) {
  // Step 1.a: Set up independent grids: first set NUMGRIDS == number of unique CoordSystems we have.
  commondata->NUMGRIDS = 1;
  {
    // Independent grids
    int Nx[3] = {-1, -1, -1};

    // For each grid, set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
    const bool apply_convergence_factor_and_set_xxminmax_defaults = true;
    int grid = 0;

    // In multipatch, gridname is a helpful alias indicating position of the patch. E.g., "lower Spherical patch"
    snprintf(griddata[grid].params.gridname, 100, "grid_Spherical");
    griddata[grid].params.CoordSystem_hash = SPHERICAL;
    griddata[grid].params.grid_physical_size = 10.0;
    numerical_grid_params_Nxx_dxx_xx(commondata, &griddata[grid].params, griddata[grid].xx, Nx, apply_convergence_factor_and_set_xxminmax_defaults);

#ifdef __CUDACC__
    // Now that numerical grid parameters have been set, copy them over to host.
    memcpy(&griddata_host[grid].params, &griddata[grid].params, sizeof(params_struct));
    // After copying params struct, set is_host=true:
    griddata_host[grid].params.is_host = true;
#endif // __CUDACC__
    grid++;
  } // END independent grid setup

  // This populates griddata_host->xx, which (if rfm_precompute is enabled) is needed for host-side rfm_precompute_defines.
  IFCUDARUN(cpyDevicetoHost__grid(commondata, griddata_host, griddata); BHAH_DEVICE_SYNC(););

  // Step 1.b: Allocate memory for and define reference-metric precomputation lookup tables.
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
#ifdef __CUDACC__
    // Set up device rfmstruct
    BHAH_MALLOC_DEVICE(griddata[grid].rfmstruct, sizeof(rfm_struct));
    cpyHosttoDevice_params__constant(&griddata[grid].params, griddata[grid].params.grid_idx % NUM_STREAMS);
    rfm_precompute_malloc(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, griddata[grid].rfmstruct, griddata[grid].xx);
    // Set up host rfmstruct
    BHAH_MALLOC(griddata_host[grid].rfmstruct, sizeof(rfm_struct));
    rfm_precompute_malloc(commondata, &griddata_host[grid].params, griddata_host[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata_host[grid].params, griddata_host[grid].rfmstruct, griddata_host[grid].xx);
#else  // if NOT defined __CUDACC__
    BHAH_MALLOC(griddata[grid].rfmstruct, sizeof(rfm_struct));
    rfm_precompute_malloc(commondata, &griddata[grid].params, griddata[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata[grid].params, griddata[grid].rfmstruct, griddata[grid].xx);
#endif // __CUDACC__
  } // END LOOP over grids

  // Step 1.c: Set up curvilinear boundary condition struct (bcstruct)
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    bcstruct_set_up(commondata, &griddata[grid].params, griddata[grid].xx, &griddata[grid].bcstruct);
  }

  // Step 1.d: Set timestep based on minimum spacing between neighboring gridpoints.
  commondata->dt = 1e30;
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {

    cfl_limited_timestep(commondata, &griddata[grid].params, griddata[grid].xx);
  }
  // Step 1.e: Initialize timestepping parameters to zero if this is the first time this function is called.
  if (calling_for_first_time) {
    commondata->nn = 0;
    commondata->nn_0 = 0;
    commondata->t_0 = 0.0;
    commondata->time = 0.0;
  }

  // Step 1.f: Set grid_idx for each grid.
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    griddata[grid].params.grid_idx = grid;
  }
} // END FUNCTION numerical_grids_and_timestep

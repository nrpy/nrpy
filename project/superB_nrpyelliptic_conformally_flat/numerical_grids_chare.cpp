#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Set up a cell-centered grids of size grid_physical_size.
 */
void numerical_grids_chare(commondata_struct *restrict commondata, griddata_struct *restrict griddata, griddata_struct *restrict griddata_chare,
                           const int chare_index[3]) {
  int grid = 0;
  snprintf(griddata_chare[grid].params.gridname, sizeof griddata_chare[grid].params.gridname, "%s", griddata[grid].params.gridname);
  griddata_chare[grid].params.CoordSystem_hash = griddata[grid].params.CoordSystem_hash;
  griddata_chare[grid].params.grid_physical_size = griddata[grid].params.grid_physical_size;

  // Step 1.b: Set Nxx & Nxx_plus_2NGHOSTS, as well as dxx, invdxx, & xx based on grid_physical_size
  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    numerical_grid_params_Nxx_dxx_xx_chare(commondata, &griddata[grid].params, &griddata_chare[grid].params, griddata_chare[grid].xx, chare_index);
  }

  // Step 1.c: Allocate memory for and define reference-metric precomputation lookup tables

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    griddata_chare[grid].rfmstruct = (rfm_struct *)malloc(sizeof(rfm_struct));
    rfm_precompute_malloc(commondata, &griddata_chare[grid].params, griddata_chare[grid].rfmstruct);
    rfm_precompute_defines(commondata, &griddata_chare[grid].params, griddata_chare[grid].rfmstruct, griddata_chare[grid].xx);
  }

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    charecommstruct_set_up(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct, chare_index);

    // Step 1.d: Set up curvilinear boundary condition struct (bcstruct)

    bcstruct_chare_set_up(commondata, &griddata[grid].params, &griddata_chare[grid].params, &griddata_chare[grid].charecommstruct,
                          griddata_chare[grid].xx, &griddata[grid].bcstruct, &griddata_chare[grid].bcstruct,
                          &griddata_chare[grid].nonlocalinnerbcstruct, chare_index);

    //NEW
    // Initialize the diagnostics struct with zero
    griddata_chare[grid].diagnosticstruct = (diagnostic_struct){0};

  }

  // 1D diagnostics set up
  diagnostics(commondata, griddata, griddata_chare, NULL, chare_index, 0, Ck::IO::Session{}, DIAGNOSTICS_SETUP_1D);

  // 2D diagnostics set up
  diagnostics(commondata, griddata, griddata_chare, NULL, chare_index, 0, Ck::IO::Session{}, DIAGNOSTICS_SETUP_2D);

} // END FUNCTION numerical_grids_chare



#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/*
 * Set up a cell-centered SinhSpherical grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.
 */
void numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical(commondata_struct *restrict commondata, params_struct *restrict params,
                                                          REAL *restrict xx[3]) {
  params->Nxx0 = 800;
  params->Nxx1 = 16;
  params->Nxx2 = 2;

  const REAL grid_physical_size = params->grid_physical_size;
  snprintf(params->CoordSystemName, 50, "SinhSpherical");

  // convergence_factor does not increase resolution across an axis of symmetry:
  if (params->Nxx0 != 2)
    params->Nxx0 *= commondata->convergence_factor;
  if (params->Nxx1 != 2)
    params->Nxx1 *= commondata->convergence_factor;
  if (params->Nxx2 != 2)
    params->Nxx2 *= commondata->convergence_factor;

  params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2 * NGHOSTS;
  params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2 * NGHOSTS;

  // Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
  params->AMPL = grid_physical_size;

  // Set xxmin, xxmax
  params->xxmin0 = 0;
  params->xxmin1 = 0;
  params->xxmin2 = -M_PI;
  params->xxmax0 = 1;
  params->xxmax1 = M_PI;
  params->xxmax2 = M_PI;

  params->dxx0 = (params->xxmax0 - params->xxmin0) / ((REAL)params->Nxx0);
  params->dxx1 = (params->xxmax1 - params->xxmin1) / ((REAL)params->Nxx1);
  params->dxx2 = (params->xxmax2 - params->xxmin2) / ((REAL)params->Nxx2);

  params->invdxx0 = ((REAL)params->Nxx0) / (params->xxmax0 - params->xxmin0);
  params->invdxx1 = ((REAL)params->Nxx1) / (params->xxmax1 - params->xxmin1);
  params->invdxx2 = ((REAL)params->Nxx2) / (params->xxmax2 - params->xxmin2);

  // Set up cell-centered Cartesian coordinate grid, centered at the origin.
  xx[0] = (REAL *restrict)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  xx[1] = (REAL *restrict)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  xx[2] = (REAL *restrict)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
  for (int j = 0; j < params->Nxx_plus_2NGHOSTS0; j++)
    xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx0;
  for (int j = 0; j < params->Nxx_plus_2NGHOSTS1; j++)
    xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx1;
  for (int j = 0; j < params->Nxx_plus_2NGHOSTS2; j++)
    xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS) + (1.0 / 2.0)) * params->dxx2;
}

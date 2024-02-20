#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Set up a cell-centered SinhSpherical grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.
 */
void numerical_grid_params_Nxx_dxx_xx(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    numerical_grid_params_Nxx_dxx_xx__rfm__SinhSpherical(commondata, params, xx);
    break;
  default:
    fprintf(stderr, "ERROR in numerical_grid_params_Nxx_dxx_xx(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

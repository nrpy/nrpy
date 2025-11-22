#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Set up a cell-centered SinhSpherical grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.
 */
void numerical_grid_params_Nxx_dxx_xx_chare(commondata_struct *restrict commondata, const params_struct *restrict params,
                                            params_struct *restrict params_chare, REAL *restrict xx[3], const int chare_index[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    numerical_grid_params_Nxx_dxx_xx_chare__rfm__SinhSpherical(commondata, params, params_chare, xx, chare_index);
    break;
  default:
    fprintf(stderr, "ERROR in numerical_grid_params_Nxx_dxx_xx_chare(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION numerical_grid_params_Nxx_dxx_xx_chare

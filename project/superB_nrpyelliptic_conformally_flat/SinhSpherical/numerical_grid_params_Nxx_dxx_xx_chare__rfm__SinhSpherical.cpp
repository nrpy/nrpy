#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Set up a cell-centered SinhSpherical grid of size grid_physical_size. Set params: Nxx, Nxx_plus_2NGHOSTS, dxx, invdxx, and xx.
 */
void numerical_grid_params_Nxx_dxx_xx_chare__rfm__SinhSpherical(commondata_struct *restrict commondata, const params_struct *restrict params,
                                                                params_struct *restrict params_chare, REAL *restrict xx[3],
                                                                const int chare_index[3]) {
  const int Nchare0 = commondata->Nchare0;
  const int Nchare1 = commondata->Nchare1;
  const int Nchare2 = commondata->Nchare2;
  if (params->Nxx0 % Nchare0 != 0) {
    fprintf(stderr, "Error: Division does not result in an integer value: Nxx0 %% Nchare0 = %d %% %d = %d\n", params->Nxx0, Nchare0,
            params->Nxx0 % Nchare0);
    exit(1);
  }
  if (params->Nxx1 % Nchare1 != 0) {
    fprintf(stderr, "Error: Division does not result in an integer value: Nxx1 %% Nchare1 = %d %% %d = %d\n", params->Nxx1, Nchare1,
            params->Nxx1 % Nchare1);
    exit(1);
  }
  if (params->Nxx2 % Nchare2 != 0) {
    fprintf(stderr, "Error: Division does not result in an integer value: Nxx2 %% Nchare2 = %d %% %d = %d\n", params->Nxx2, Nchare2,
            params->Nxx2 % Nchare2);
    exit(1);
  }
  if (Nchare0 > 1 && params->Nxx0 / Nchare0 < NGHOSTS) {
    fprintf(stderr, "Error: params->Nxx0 / Nchare0 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx0, Nchare0, params->Nxx0 / Nchare0,
            NGHOSTS);
    exit(1);
  }
  if (Nchare1 > 1 && params->Nxx1 / Nchare1 < NGHOSTS) {
    fprintf(stderr, "Error: params->Nxx1 / Nchare1 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx1, Nchare1, params->Nxx1 / Nchare1,
            NGHOSTS);
    exit(1);
  }
  if (Nchare2 > 1 && params->Nxx2 / Nchare2 < NGHOSTS) {
    fprintf(stderr, "Error: params->Nxx2 / Nchare2 is less than NGHOSTS: %d / %d = %d < %d\n", params->Nxx2, Nchare2, params->Nxx2 / Nchare2,
            NGHOSTS);
    exit(1);
  }
  params_chare->Nxx0 = params->Nxx0 / Nchare0;
  params_chare->Nxx1 = params->Nxx1 / Nchare1;
  params_chare->Nxx2 = params->Nxx2 / Nchare2;

  const REAL grid_physical_size = params_chare->grid_physical_size;
  snprintf(params_chare->CoordSystemName, 100, "SinhSpherical");

  params_chare->Nxx_plus_2NGHOSTS0 = params_chare->Nxx0 + 2 * NGHOSTS;
  params_chare->Nxx_plus_2NGHOSTS1 = params_chare->Nxx1 + 2 * NGHOSTS;
  params_chare->Nxx_plus_2NGHOSTS2 = params_chare->Nxx2 + 2 * NGHOSTS;

  // Set grid size to grid_physical_size (set above, based on params->grid_physical_size):
  params_chare->AMPL = grid_physical_size;

  // Set xxmin, xxmax
  params_chare->xxmin0 = params->xxmin0 + (params->dxx0 * (REAL)(params_chare->Nxx0 * chare_index[0]));
  params_chare->xxmin1 = params->xxmin1 + (params->dxx1 * (REAL)(params_chare->Nxx1 * chare_index[1]));
  params_chare->xxmin2 = params->xxmin2 + (params->dxx2 * (REAL)(params_chare->Nxx2 * chare_index[2]));
  params_chare->xxmax0 = params->xxmax0 - (params->dxx0 * (REAL)(params_chare->Nxx0 * (Nchare0 - 1 - chare_index[0])));
  params_chare->xxmax1 = params->xxmax1 - (params->dxx1 * (REAL)(params_chare->Nxx1 * (Nchare1 - 1 - chare_index[1])));
  params_chare->xxmax2 = params->xxmax2 - (params->dxx2 * (REAL)(params_chare->Nxx2 * (Nchare2 - 1 - chare_index[2])));

  params_chare->dxx0 = params->dxx0;
  params_chare->dxx1 = params->dxx1;
  params_chare->dxx2 = params->dxx2;

  params_chare->invdxx0 = params->invdxx0;
  params_chare->invdxx1 = params->invdxx1;
  params_chare->invdxx2 = params->invdxx2;

  // Set up cell-centered Cartesian coordinate grid, centered at the origin.
  xx[0] = (REAL *restrict)malloc(sizeof(REAL) * params_chare->Nxx_plus_2NGHOSTS0);
  xx[1] = (REAL *restrict)malloc(sizeof(REAL) * params_chare->Nxx_plus_2NGHOSTS1);
  xx[2] = (REAL *restrict)malloc(sizeof(REAL) * params_chare->Nxx_plus_2NGHOSTS2);
  for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS0; j++)
    xx[0][j] = params->xxmin0 + ((REAL)(j - NGHOSTS + (params_chare->Nxx0 * chare_index[0])) + (1.0 / 2.0)) * params_chare->dxx0;
  for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS1; j++)
    xx[1][j] = params->xxmin1 + ((REAL)(j - NGHOSTS + (params_chare->Nxx1 * chare_index[1])) + (1.0 / 2.0)) * params_chare->dxx1;
  for (int j = 0; j < params_chare->Nxx_plus_2NGHOSTS2; j++)
    xx[2][j] = params->xxmin2 + ((REAL)(j - NGHOSTS + (params_chare->Nxx2 * chare_index[2])) + (1.0 / 2.0)) * params_chare->dxx2;
} // END FUNCTION numerical_grid_params_Nxx_dxx_xx_chare__rfm__SinhSpherical

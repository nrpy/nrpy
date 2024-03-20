#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Output minimum gridspacing ds_min on a SinhSpherical numerical grid.
 */
void cfl_limited_timestep(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                          bc_struct *restrict bcstruct) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    cfl_limited_timestep__rfm__SinhSpherical(commondata, params, xx, bcstruct);
    break;
  default:
    fprintf(stderr, "ERROR in cfl_limited_timestep(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

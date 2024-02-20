#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Output diagnostic quantities at grid's *physical* center.
 * For example:
 * In Cartesian this will be at i0_mid,i1_mid,i2_mid.
 * In Spherical, this will be at i0_min,i1_mid,i2_mid (i1 and i2 don't matter).
 * In Cylindrical, this will be at i0_min,i1_mid,i2_mid (i1 == phi doesn't matter).
 * In SinhSymTP, this will be at i0_min,i1_mid,i2_mid (i2 == phi doesn't matter).
 */
void diagnostics_nearest_grid_center(commondata_struct *restrict commondata, const params_struct *restrict params,
                                     MoL_gridfunctions_struct *restrict gridfuncs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    diagnostics_nearest_grid_center__rfm__SinhSpherical(commondata, params, gridfuncs);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_grid_center(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

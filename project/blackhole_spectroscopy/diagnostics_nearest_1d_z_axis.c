#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Output diagnostic quantities at gridpoints closest to z axis.
 */
void diagnostics_nearest_1d_z_axis(commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                   MoL_gridfunctions_struct *restrict gridfuncs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    diagnostics_nearest_1d_z_axis__rfm__SinhSpherical(commondata, params, xx, gridfuncs);
    break;
  default:
    fprintf(stderr, "ERROR in diagnostics_nearest_1d_z_axis(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

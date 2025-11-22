#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Call functions that set up all AUXEVOL gridfunctions.
 */
void auxevol_gfs_set_to_constant(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                 MoL_gridfunctions_struct *restrict gridfuncs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    auxevol_gfs_set_to_constant__rfm__SinhSpherical(commondata, params, xx, gridfuncs);
    break;
  default:
    fprintf(stderr, "ERROR in auxevol_gfs_set_to_constant(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION auxevol_gfs_set_to_constant

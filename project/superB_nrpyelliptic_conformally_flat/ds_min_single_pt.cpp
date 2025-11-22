#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min.
 */
void ds_min_single_pt(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    ds_min_single_pt__rfm__SinhSpherical(params, xx0, xx1, xx2, ds_min);
    break;
  default:
    fprintf(stderr, "ERROR in ds_min_single_pt(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION ds_min_single_pt

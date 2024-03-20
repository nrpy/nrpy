#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2),   accounting for the origin of this grid being possibly off-center.
 */
void xx_to_Cart(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const int i0, const int i1,
                const int i2, REAL xCart[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    xx_to_Cart__rfm__SinhSpherical(commondata, params, xx, i0, i1, i2, xCart);
    break;
  default:
    fprintf(stderr, "ERROR in xx_to_Cart(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

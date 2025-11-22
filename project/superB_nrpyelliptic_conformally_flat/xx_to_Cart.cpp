#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {


  // DEBUG: minimal params dump
  //~ fprintf(stderr,
          //~ "[xx_to_Cart] params=%p hash=%d Np2G=(%d,%d,%d) xx=(%.15e, %.15e, %.15e)\n",
          //~ (const void *)params,
          //~ params->CoordSystem_hash,
          //~ params->Nxx_plus_2NGHOSTS0,
          //~ params->Nxx_plus_2NGHOSTS1,
          //~ params->Nxx_plus_2NGHOSTS2,
          //~ xx[0], xx[1], xx[2]);



  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    xx_to_Cart__rfm__SinhSpherical(params, xx, xCart);
    break;
  default:
    fprintf(stderr, "ERROR in xx_to_Cart(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION xx_to_Cart

#include "BHaH_defines.h"

/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} from the
 * local coordinate vector {xx[0], xx[1], xx[2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart__rfm__SymTP(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->Cart_originx + xx0*sin(xx1)*cos(xx2)]"
   *  "[xCart[1] = params->Cart_originy + xx0*sin(xx1)*sin(xx2)]"
   *  "[xCart[2] = params->Cart_originz + sqrt(params->bScale**2 + xx0**2)*cos(xx1)]"
   */
  {
    const REAL tmp0 = xx0 * sin(xx1);
    xCart[0] = params->Cart_originx + tmp0 * cos(xx2);
    xCart[1] = params->Cart_originy + tmp0 * sin(xx2);
    xCart[2] = params->Cart_originz + sqrt(((params->bScale) * (params->bScale)) + ((xx0) * (xx0))) * cos(xx1);
  }
} // END FUNCTION xx_to_Cart__rfm__SymTP

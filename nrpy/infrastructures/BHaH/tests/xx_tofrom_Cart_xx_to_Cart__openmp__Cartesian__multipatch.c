#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center,
 *  and the grid is rotated.
 */
void xx_to_Cart__rfm__Cartesian(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->Cart_originx + xx0]"
   *  "[xCart[1] = params->Cart_originy + xx1]"
   *  "[xCart[2] = params->Cart_originz + xx2]"
   */
  {
    xCart[0] = params->Cart_originx + xx0;
    xCart[1] = params->Cart_originy + xx1;
    xCart[2] = params->Cart_originz + xx2;
  }

  if (params->grid_rotates) {
    REAL R[3][3];
    so3_build_R_from_hats(params->cumulatively_rotated_xhatU, params->cumulatively_rotated_yhatU, params->cumulatively_rotated_zhatU, R);
    so3_apply_R_to_vector(R, xCart);
  } // END IF: grid rotates
} // END FUNCTION: xx_to_Cart__rfm__Cartesian

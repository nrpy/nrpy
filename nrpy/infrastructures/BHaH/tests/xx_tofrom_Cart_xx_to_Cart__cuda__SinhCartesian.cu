#include "BHaH_defines.h"

/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} from the
 * local coordinate vector {xx[0], xx[1], xx[2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
__host__ __device__ void xx_to_Cart__rfm__SinhCartesian(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->AMPLXYZ*(exp(xx0/params->SINHWXYZ) - exp(-xx0/params->SINHWXYZ))/(exp(1/params->SINHWXYZ) - exp(-1/params->SINHWXYZ)) +
   * params->Cart_originx]"
   *  "[xCart[1] = params->AMPLXYZ*(exp(xx1/params->SINHWXYZ) - exp(-xx1/params->SINHWXYZ))/(exp(1/params->SINHWXYZ) - exp(-1/params->SINHWXYZ)) +
   * params->Cart_originy]"
   *  "[xCart[2] = params->AMPLXYZ*(exp(xx2/params->SINHWXYZ) - exp(-xx2/params->SINHWXYZ))/(exp(1/params->SINHWXYZ) - exp(-1/params->SINHWXYZ)) +
   * params->Cart_originz]"
   */
  {
    const REAL tmp0 = (1.0 / (params->SINHWXYZ));
    const REAL tmp2 = params->AMPLXYZ / (exp(tmp0) - exp(-tmp0));
    xCart[0] = params->Cart_originx + tmp2 * (exp(tmp0 * xx0) - exp(-tmp0 * xx0));
    xCart[1] = params->Cart_originy + tmp2 * (exp(tmp0 * xx1) - exp(-tmp0 * xx1));
    xCart[2] = params->Cart_originz + tmp2 * (exp(tmp0 * xx2) - exp(-tmp0 * xx2));
  }
} // END FUNCTION xx_to_Cart__rfm__SinhCartesian

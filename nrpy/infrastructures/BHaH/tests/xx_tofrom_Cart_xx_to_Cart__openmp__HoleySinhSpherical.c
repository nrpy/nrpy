#include "../BHaH_defines.h"
/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart__rfm__HoleySinhSpherical(const params_struct *restrict params, REAL xx[3], REAL xCart[3]) {

  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->AMPL*(exp(xx0/params->SINHW) - exp(-xx0/params->SINHW))*sin(xx1)*cos(xx2)/(exp(1/params->SINHW) - exp(-1/params->SINHW)) +
   * params->Cart_originx]"
   *  "[xCart[1] = params->AMPL*(exp(xx0/params->SINHW) - exp(-xx0/params->SINHW))*sin(xx1)*sin(xx2)/(exp(1/params->SINHW) - exp(-1/params->SINHW)) +
   * params->Cart_originy]"
   *  "[xCart[2] = params->AMPL*(exp(xx0/params->SINHW) - exp(-xx0/params->SINHW))*cos(xx1)/(exp(1/params->SINHW) - exp(-1/params->SINHW)) +
   * params->Cart_originz]"
   */
  {
    const REAL tmp0 = (1.0 / (params->SINHW));
    const REAL tmp2 = params->AMPL * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    const REAL tmp3 = tmp2 * sin(xx1);
    xCart[0] = params->Cart_originx + tmp3 * cos(xx2);
    xCart[1] = params->Cart_originy + tmp3 * sin(xx2);
    xCart[2] = params->Cart_originz + tmp2 * cos(xx1);
  }
} // END FUNCTION xx_to_Cart__rfm__HoleySinhSpherical

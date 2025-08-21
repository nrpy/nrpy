#include "../BHaH_defines.h"
/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart__rfm__SinhSymTP(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {

  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->AMAX*(exp(xx0/params->SINHWAA) - exp(-xx0/params->SINHWAA))*sin(xx1)*cos(xx2)/(exp(1/params->SINHWAA) -
   * exp(-1/params->SINHWAA)) + params->Cart_originx]"
   *  "[xCart[1] = params->AMAX*(exp(xx0/params->SINHWAA) - exp(-xx0/params->SINHWAA))*sin(xx1)*sin(xx2)/(exp(1/params->SINHWAA) -
   * exp(-1/params->SINHWAA)) + params->Cart_originy]"
   *  "[xCart[2] = params->Cart_originz + sqrt(params->AMAX**2*(exp(xx0/params->SINHWAA) - exp(-xx0/params->SINHWAA))**2/(exp(1/params->SINHWAA) -
   * exp(-1/params->SINHWAA))**2 + params->bScale**2)*cos(xx1)]"
   */
  {
    const REAL tmp0 = (1.0 / (params->SINHWAA));
    const REAL tmp1 = exp(tmp0) - exp(-tmp0);
    const REAL tmp3 = exp(tmp0 * xx0) - exp(-tmp0 * xx0);
    const REAL tmp4 = params->AMAX * tmp3 * sin(xx1) / tmp1;
    xCart[0] = params->Cart_originx + tmp4 * cos(xx2);
    xCart[1] = params->Cart_originy + tmp4 * sin(xx2);
    xCart[2] = params->Cart_originz +
               sqrt(((params->AMAX) * (params->AMAX)) * ((tmp3) * (tmp3)) / ((tmp1) * (tmp1)) + ((params->bScale) * (params->bScale))) * cos(xx1);
  }
} // END FUNCTION xx_to_Cart__rfm__SinhSymTP

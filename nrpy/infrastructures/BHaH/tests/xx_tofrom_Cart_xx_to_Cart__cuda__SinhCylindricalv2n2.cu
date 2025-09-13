#include "BHaH_defines.h"
/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} in terms of
 * local grid coordinates {xx[0][i0], xx[1][i1], xx[2][i2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
__host__ __device__ void xx_to_Cart__rfm__SinhCylindricalv2n2(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = params->Cart_originx + (params->rho_slope*xx0 + xx0**2*(params->AMPLRHO - params->rho_slope)*(exp(xx0/params->SINHWRHO) -
   * exp(-xx0/params->SINHWRHO))/(exp(1/params->SINHWRHO) - exp(-1/params->SINHWRHO)))*cos(xx1)]"
   *  "[xCart[1] = params->Cart_originy + (params->rho_slope*xx0 + xx0**2*(params->AMPLRHO - params->rho_slope)*(exp(xx0/params->SINHWRHO) -
   * exp(-xx0/params->SINHWRHO))/(exp(1/params->SINHWRHO) - exp(-1/params->SINHWRHO)))*sin(xx1)]"
   *  "[xCart[2] = params->Cart_originz + params->z_slope*xx2 + xx2**2*(params->AMPLZ - params->z_slope)*(exp(xx2/params->SINHWZ) -
   * exp(-xx2/params->SINHWZ))/(exp(1/params->SINHWZ) - exp(-1/params->SINHWZ))]"
   */
  {
    const REAL tmp0 = (1.0 / (params->SINHWRHO));
    const REAL tmp3 = (1.0 / (params->SINHWZ));
    const REAL tmp2 = params->rho_slope * xx0 +
                      ((xx0) * (xx0)) * (params->AMPLRHO - params->rho_slope) * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    xCart[0] = params->Cart_originx + tmp2 * cos(xx1);
    xCart[1] = params->Cart_originy + tmp2 * sin(xx1);
    xCart[2] = params->Cart_originz + params->z_slope * xx2 +
               ((xx2) * (xx2)) * (params->AMPLZ - params->z_slope) * (exp(tmp3 * xx2) - exp(-tmp3 * xx2)) / (exp(tmp3) - exp(-tmp3));
  }
} // END FUNCTION xx_to_Cart__rfm__SinhCylindricalv2n2

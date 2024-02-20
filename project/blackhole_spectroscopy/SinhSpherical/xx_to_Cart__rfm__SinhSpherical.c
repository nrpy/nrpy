#include "../BHaH_defines.h"
/*
 * Compute Cartesian coordinates given local grid coordinate (xx0,xx1,xx2),   accounting for the origin of this grid being possibly off-center.
 */
void xx_to_Cart__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
                                    const int i0, const int i1, const int i2, REAL xCart[3]) {
#include "../set_CodeParameters.h"

  REAL xx0 = xx[0][i0];
  REAL xx1 = xx[1][i1];
  REAL xx2 = xx[2][i2];
  /*
   *  Original SymPy expressions:
   *  "[xCart[0] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*cos(xx2)/(exp(1/SINHW) - exp(-1/SINHW)) + Cart_originx]"
   *  "[xCart[1] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)*sin(xx2)/(exp(1/SINHW) - exp(-1/SINHW)) + Cart_originy]"
   *  "[xCart[2] = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))*cos(xx1)/(exp(1/SINHW) - exp(-1/SINHW)) + Cart_originz]"
   */
  {
    const REAL tmp0 = (1.0 / (SINHW));
    const REAL tmp2 = AMPL * (exp(tmp0 * xx0) - exp(-tmp0 * xx0)) / (exp(tmp0) - exp(-tmp0));
    const REAL tmp3 = tmp2 * sin(xx1);
    xCart[0] = Cart_originx + tmp3 * cos(xx2);
    xCart[1] = Cart_originy + tmp3 * sin(xx2);
    xCart[2] = Cart_originz + tmp2 * cos(xx1);
  }
}

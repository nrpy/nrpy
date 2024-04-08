#include "../BHaH_defines.h"
/*
 * Given Cartesian point (x,y,z), this function outputs the corresponding
 * (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                       const REAL xCart[3], REAL xx[3], int Cart_to_i0i1i2[3]) {
#include "../set_CodeParameters.h"

  // See comments for description on how coordinates are computed relative to the local grid center.
  const REAL Cartx = xCart[0] - Cart_originx;
  const REAL Carty = xCart[1] - Cart_originy;
  const REAL Cartz = xCart[2] - Cart_originz;
  /*
   *  Original SymPy expressions:
   *  "[xx[0] = SINHW*asinh(sqrt(Cartx**2 + Carty**2 + Cartz**2)*sinh(1/SINHW)/AMPL)]"
   *  "[xx[1] = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2))]"
   *  "[xx[2] = atan2(Carty, Cartx)]"
   */
  const REAL tmp0 = sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
  xx[0] = SINHW * asinh(tmp0 * sinh((1.0 / (SINHW))) / AMPL);
  xx[1] = acos(Cartz / tmp0);
  xx[2] = atan2(Carty, Cartx);

  // Then find the nearest index (i0,i1,i2) on underlying grid to (x,y,z)
  Cart_to_i0i1i2[0] = (int)((xx[0] - (0)) / params->dxx0 + 0.5 + NGHOSTS - 0.5);     // Account for (int) typecast rounding down
  Cart_to_i0i1i2[1] = (int)((xx[1] - (0)) / params->dxx1 + 0.5 + NGHOSTS - 0.5);     // Account for (int) typecast rounding down
  Cart_to_i0i1i2[2] = (int)((xx[2] - (-M_PI)) / params->dxx2 + 0.5 + NGHOSTS - 0.5); // Account for (int) typecast rounding down
}

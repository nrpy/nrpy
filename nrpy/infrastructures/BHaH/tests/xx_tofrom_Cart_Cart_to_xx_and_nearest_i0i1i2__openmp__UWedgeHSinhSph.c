#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2__rfm__UWedgeHSinhSph(const params_struct *restrict params, const REAL xCart[3], REAL xx[3],
                                                        int Cart_to_i0i1i2[3]) {
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];

  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= params->Cart_originx;
  Carty -= params->Cart_originy;
  Cartz -= params->Cart_originz;
  {
    /*
     *  Original SymPy expressions:
     *  "[xx[0] = params->SINHW*asinh(sqrt(Cartx**2 + Carty**2 + Cartz**2)*sinh(1/params->SINHW)/params->AMPL)]"
     *  "[xx[1] = acos(-Cartx/sqrt(Cartx**2 + Carty**2 + Cartz**2))]"
     *  "[xx[2] = atan2(Carty, Cartz)]"
     */
    const REAL tmp0 = sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
    xx[0] = params->SINHW * asinh(tmp0 * sinh((1.0 / (params->SINHW))) / params->AMPL);
    xx[1] = acos(-Cartx / tmp0);
    xx[2] = atan2(Carty, Cartz);

    // Find the nearest grid indices (i0, i1, i2) for the given Cartesian coordinates (x, y, z).
    // Assuming a cell-centered grid, which follows the pattern:
    //   xx0[i0] = params->xxmin0 + ((REAL)(i0 - NGHOSTS) + 0.5) * params->dxx0
    // The index i0 can be derived as:
    //   i0 = (xx0[i0] - params->xxmin0) / params->dxx0 - 0.5 + NGHOSTS
    // Now, including typecasts:
    //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS)
    // C float-to-int conversion truncates toward zero; for nonnegative inputs this matches floor().
    // Assuming (xx - xxmin)/dxx + NGHOSTS is nonnegative (typical for valid interior points), this is safe.
    //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS + 0.5)
    // The 0.5 values cancel out:
    //   i0 =           (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS )
    Cart_to_i0i1i2[0] = (int)((xx[0] - params->xxmin0) / params->dxx0 + (REAL)NGHOSTS);
    Cart_to_i0i1i2[1] = (int)((xx[1] - params->xxmin1) / params->dxx1 + (REAL)NGHOSTS);
    Cart_to_i0i1i2[2] = (int)((xx[2] - params->xxmin2) / params->dxx2 + (REAL)NGHOSTS);
  }
} // END FUNCTION Cart_to_xx_and_nearest_i0i1i2__rfm__UWedgeHSinhSph

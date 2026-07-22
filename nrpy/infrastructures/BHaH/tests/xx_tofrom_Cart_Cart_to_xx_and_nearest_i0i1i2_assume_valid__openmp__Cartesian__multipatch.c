#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Given Cartesian point (x,y,z), this function first unrotates to the co-rotating frame, and then unshifts the grid back to the origin to output the
 * corresponding (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid. The index output may be NULL when only logical coordinates are needed.
 */
void Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__Cartesian(const params_struct *restrict params, const REAL xCart[3], REAL xx[3],
                                                                int Cart_to_i0i1i2[3]) {
  // Set (Cartx, Carty, Cartz) relative to the global (as opposed to local) grid.
  //   This local grid may be offset from the origin by adjusting
  //   (Cart_originx, Cart_originy, Cart_originz) to nonzero values.
  REAL Cartx = xCart[0];
  REAL Carty = xCart[1];
  REAL Cartz = xCart[2];

  if (params->grid_rotates) {
    REAL R[3][3];
    REAL vU[3] = {Cartx, Carty, Cartz};
    so3_build_R_from_hats(params->cumulatively_rotated_xhatU, params->cumulatively_rotated_yhatU, params->cumulatively_rotated_zhatU, R);
    so3_apply_RT_to_vector(R, vU);
    Cartx = vU[0];
    Carty = vU[1];
    Cartz = vU[2];
  } // END IF: grid rotates

  // Set the origin, (Cartx, Carty, Cartz) = (0, 0, 0), to the center of the local grid patch.
  Cartx -= params->Cart_originx;
  Carty -= params->Cart_originy;
  Cartz -= params->Cart_originz;
  {
    /*
     *  Original SymPy expressions:
     *  "[xx[0] = Cartx]"
     *  "[xx[1] = Carty]"
     *  "[xx[2] = Cartz]"
     */
    xx[0] = Cartx;
    xx[1] = Carty;
    xx[2] = Cartz;

    // A NULL index output requests logical coordinates only. In particular,
    // recoverable callers use this path to validate logical coordinates before
    // performing any floating-to-integer conversion.
    if (Cart_to_i0i1i2 == NULL)
      return;

    // Find the nearest grid indices (i0, i1, i2) for the given Cartesian coordinates (x, y, z).
    // Assuming a cell-centered grid, which follows the pattern:
    //   xx0[i0] = params->xxmin0 + ((REAL)(i0 - NGHOSTS) + 0.5) * params->dxx0
    // The index i0 can be derived as:
    //   i0 = (xx0[i0] - params->xxmin0) / params->dxx0 - 0.5 + NGHOSTS
    // Now, including typecasts:
    //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS)
    // C integer conversion truncates toward zero. For valid in-domain cell-centered
    // coordinates, adding 0.5 converts the center-based quantity to the nearest index:
    //   i0 = (int)((xx[0] - params->xxmin0) / params->dxx0 - 0.5 + (REAL)NGHOSTS + 0.5)
    // The 0.5 values cancel out:
    //   i0 =           (int)( ( xx[0] - params->xxmin0 ) / params->dxx0 + (REAL)NGHOSTS )
    Cart_to_i0i1i2[0] = (int)((xx[0] - params->xxmin0) / params->dxx0 + (REAL)NGHOSTS);
    Cart_to_i0i1i2[1] = (int)((xx[1] - params->xxmin1) / params->dxx1 + (REAL)NGHOSTS);
    Cart_to_i0i1i2[2] = (int)((xx[2] - params->xxmin2) / params->dxx2 + (REAL)NGHOSTS);
  } // END BLOCK: Cartesian-to-grid conversion
} // END FUNCTION: Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__Cartesian

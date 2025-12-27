#include "BHaH_defines.h"

/**
 * Compute Cartesian coordinates {x, y, z} = {xCart[0], xCart[1], xCart[2]} from the
 * local coordinate vector {xx[0], xx[1], xx[2]} = {xx0, xx1, xx2},
 * taking into account the possibility that the origin of this grid is off-center.
 */
void xx_to_Cart__rfm__GeneralRFM_FisheyeN2(const params_struct *restrict params, const REAL xx[3], REAL xCart[3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  const REAL r2 = xx0 * xx0 + xx1 * xx1 + xx2 * xx2;

  if (r2 <= (REAL)0.0) {
    // Fisheye map sends the origin to the origin (plus any patch offset).
    xCart[0] = params->Cart_originx;
    xCart[1] = params->Cart_originy;
    xCart[2] = params->Cart_originz;
  } else {
    const REAL r = sqrt(r2);
    REAL rbar;
    /*
     *  Original SymPy expression:
     *  "rbar = params->fisheye_c*(params->fisheye_a2*r + params->fisheye_s1*(params->fisheye_a0 - params->fisheye_a1)*log(cosh((params->fisheye_R1 +
     * r)/params->fisheye_s1)/cosh((-params->fisheye_R1 + r)/params->fisheye_s1))/(2*tanh(params->fisheye_R1/params->fisheye_s1)) +
     * params->fisheye_s2*(params->fisheye_a1 - params->fisheye_a2)*log(cosh((params->fisheye_R2 + r)/params->fisheye_s2)/cosh((-params->fisheye_R2 +
     * r)/params->fisheye_s2))/(2*tanh(params->fisheye_R2/params->fisheye_s2)))"
     */
    const REAL tmp0 = (1.0 / (params->fisheye_s1));
    const REAL tmp2 = (1.0 / (params->fisheye_s2));
    rbar = params->fisheye_c *
           (params->fisheye_a2 * r +
            (1.0 / 2.0) * params->fisheye_s1 * (params->fisheye_a0 - params->fisheye_a1) *
                log(cosh(tmp0 * (params->fisheye_R1 + r)) / cosh(tmp0 * (-params->fisheye_R1 + r))) / tanh(params->fisheye_R1 * tmp0) +
            (1.0 / 2.0) * params->fisheye_s2 * (params->fisheye_a1 - params->fisheye_a2) *
                log(cosh(tmp2 * (params->fisheye_R2 + r)) / cosh(tmp2 * (-params->fisheye_R2 + r))) / tanh(params->fisheye_R2 * tmp2));

    const REAL lam = rbar / r;

    xCart[0] = lam * xx0 + params->Cart_originx;
    xCart[1] = lam * xx1 + params->Cart_originy;
    xCart[2] = lam * xx2 + params->Cart_originz;
  }
} // END FUNCTION xx_to_Cart__rfm__GeneralRFM_FisheyeN2

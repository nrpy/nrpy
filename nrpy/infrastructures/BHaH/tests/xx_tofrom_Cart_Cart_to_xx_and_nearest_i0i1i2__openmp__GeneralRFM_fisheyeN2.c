#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2__rfm__GeneralRFM_fisheyeN2(const params_struct *restrict params, const REAL xCart[3], REAL xx[3],
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

    const REAL rCart = sqrt(Cartx * Cartx + Carty * Carty + Cartz * Cartz);
    if (rCart <= (REAL)0.0) {
      xx[0] = (REAL)0.0;
      xx[1] = (REAL)0.0;
      xx[2] = (REAL)0.0;
    } else {
      const REAL XX_TOLERANCE = (REAL)1e-12;
      const REAL F_TOLERANCE = (REAL)1e-12;
      const int ITER_MAX = 100;

      // Use a robust scale for convergence tests:
      const REAL dxx_scale = (params->dxx0 + params->dxx1 + params->dxx2) / (REAL)3.0;
      const REAL rscale = (rCart > dxx_scale) ? rCart : dxx_scale;

      int iter = 0;
      int tolerance_has_been_met = 0;

      // Two heuristic initial guesses:
      //   Near origin: rbar ~ c*a0*r  => r ~ rCart/(c*a0)
      //   Far field  : rbar ~ c*aN*r  => r ~ rCart/(c*aN)
      REAL r_guess0 = rCart / (params->fisheye_c * params->fisheye_a0);
      REAL r_guessN = rCart / (params->fisheye_c * params->fisheye_a2);
      if (!(r_guess0 > (REAL)0.0))
        r_guess0 = rCart;
      if (!(r_guessN > (REAL)0.0))
        r_guessN = rCart;

      // Pick the initial guess that yields the smaller |f(r)|.
      REAL r = r_guessN;
      REAL f_of_r, fprime_of_r;
      {
        const REAL tmp0 = (1.0 / (params->fisheye_s1));
        const REAL tmp5 = (1.0 / (params->fisheye_s2));
        const REAL tmp1 = tmp0 * (params->fisheye_R1 + r);
        const REAL tmp3 = tmp0 * (-params->fisheye_R1 + r);
        const REAL tmp4 = (1.0 / 2.0) * (params->fisheye_a0 - params->fisheye_a1) / tanh(params->fisheye_R1 * tmp0);
        const REAL tmp6 = tmp5 * (params->fisheye_R2 + r);
        const REAL tmp7 = tmp5 * (-params->fisheye_R2 + r);
        const REAL tmp8 = (1.0 / 2.0) * (params->fisheye_a1 - params->fisheye_a2) / tanh(params->fisheye_R2 * tmp5);
        f_of_r = params->fisheye_c * (params->fisheye_a2 * r + params->fisheye_s1 * tmp4 * log(cosh(tmp1) / cosh(tmp3)) +
                                      params->fisheye_s2 * tmp8 * log(cosh(tmp6) / cosh(tmp7))) -
                 rCart;
        fprime_of_r = params->fisheye_c * (params->fisheye_a2 + tmp4 * (tanh(tmp1) - tanh(tmp3)) + tmp8 * (tanh(tmp6) - tanh(tmp7)));
      }

      REAL fN = fabs(f_of_r);

      r = r_guess0;
      {
        const REAL tmp0 = (1.0 / (params->fisheye_s1));
        const REAL tmp5 = (1.0 / (params->fisheye_s2));
        const REAL tmp1 = tmp0 * (params->fisheye_R1 + r);
        const REAL tmp3 = tmp0 * (-params->fisheye_R1 + r);
        const REAL tmp4 = (1.0 / 2.0) * (params->fisheye_a0 - params->fisheye_a1) / tanh(params->fisheye_R1 * tmp0);
        const REAL tmp6 = tmp5 * (params->fisheye_R2 + r);
        const REAL tmp7 = tmp5 * (-params->fisheye_R2 + r);
        const REAL tmp8 = (1.0 / 2.0) * (params->fisheye_a1 - params->fisheye_a2) / tanh(params->fisheye_R2 * tmp5);
        f_of_r = params->fisheye_c * (params->fisheye_a2 * r + params->fisheye_s1 * tmp4 * log(cosh(tmp1) / cosh(tmp3)) +
                                      params->fisheye_s2 * tmp8 * log(cosh(tmp6) / cosh(tmp7))) -
                 rCart;
        fprime_of_r = params->fisheye_c * (params->fisheye_a2 + tmp4 * (tanh(tmp1) - tanh(tmp3)) + tmp8 * (tanh(tmp6) - tanh(tmp7)));
      }

      REAL f0 = fabs(f_of_r);

      r = (f0 < fN) ? r_guess0 : r_guessN;

      while (iter < ITER_MAX && !tolerance_has_been_met) {

        {
          const REAL tmp0 = (1.0 / (params->fisheye_s1));
          const REAL tmp5 = (1.0 / (params->fisheye_s2));
          const REAL tmp1 = tmp0 * (params->fisheye_R1 + r);
          const REAL tmp3 = tmp0 * (-params->fisheye_R1 + r);
          const REAL tmp4 = (1.0 / 2.0) * (params->fisheye_a0 - params->fisheye_a1) / tanh(params->fisheye_R1 * tmp0);
          const REAL tmp6 = tmp5 * (params->fisheye_R2 + r);
          const REAL tmp7 = tmp5 * (-params->fisheye_R2 + r);
          const REAL tmp8 = (1.0 / 2.0) * (params->fisheye_a1 - params->fisheye_a2) / tanh(params->fisheye_R2 * tmp5);
          f_of_r = params->fisheye_c * (params->fisheye_a2 * r + params->fisheye_s1 * tmp4 * log(cosh(tmp1) / cosh(tmp3)) +
                                        params->fisheye_s2 * tmp8 * log(cosh(tmp6) / cosh(tmp7))) -
                   rCart;
          fprime_of_r = params->fisheye_c * (params->fisheye_a2 + tmp4 * (tanh(tmp1) - tanh(tmp3)) + tmp8 * (tanh(tmp6) - tanh(tmp7)));
        }

        // Unnecessary guard against division by zero in Newton step;
        //   valid coordinate systems must have f'(r) > 0
        // if(fprime_of_r == (REAL)0.0) {
        //  break;
        // }

        const REAL r_np1_unclamped = r - f_of_r / fprime_of_r;

        // Keep r nonnegative (fisheye assumes r >= 0 and rbar(r) is odd/monotone).
        REAL r_np1 = r_np1_unclamped;
        if (r_np1 <= (REAL)0.0)
          r_np1 = (REAL)0.5 * r;

        if (fabs(r - r_np1) <= XX_TOLERANCE * rscale && fabs(f_of_r) <= F_TOLERANCE * rscale) {
          tolerance_has_been_met = 1;
        }
        r = r_np1;
        iter++;
      }

      if (iter >= ITER_MAX || !tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: Newton-Raphson failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n", (double)rCart,
               (double)Cartx, (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: Newton-Raphson failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n", (double)rCart,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
#endif
      }

      const REAL scale = r / rCart;
      xx[0] = scale * Cartx;
      xx[1] = scale * Carty;
      xx[2] = scale * Cartz;
    }

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
} // END FUNCTION Cart_to_xx_and_nearest_i0i1i2__rfm__GeneralRFM_fisheyeN2

#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
__host__ __device__ void Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindricalv2n2(const params_struct *restrict params, const REAL xCart[3],
                                                                                 REAL xx[3], int Cart_to_i0i1i2[3]) {
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
    // First compute analytical coordinate inversions:
    /*
     *  Original SymPy expression:
     *  "xx[1] = atan2(Carty, Cartx)"
     */
    xx[1] = atan2(Carty, Cartx);

    // Next perform Newton-Raphson iterations as needed:
    const REAL XX_TOLERANCE = 1e-12;      // that's 1 part in 1e12 dxxi.
    const REAL F_OF_XX_TOLERANCE = 1e-12; // tolerance of function for which we're finding the root.
    const int ITER_MAX = 100;
    int iter;

    {
      int tolerance_has_been_met = 0;
      iter = 0;
      REAL xx0 = (REAL)0.5 * (params->xxmin0 + params->xxmax0);
      while (iter < ITER_MAX && !tolerance_has_been_met) {
        REAL f_of_xx0, fprime_of_xx0;

        {
          const REAL tmp0 = (1.0 / (params->SINHWRHO));
          const REAL tmp6 = params->AMPLRHO - params->rho_slope;
          const REAL tmp5 = (1.0 / (exp(tmp0) - exp(-tmp0)));
          const REAL tmp2 = exp(tmp0 * xx0);
          const REAL tmp3 = exp(-tmp0 * xx0);
          const REAL tmp7 = tmp5 * tmp6 * ((xx0) * (xx0));
          const REAL tmp4 = tmp2 - tmp3;
          f_of_xx0 = params->rho_slope * xx0 + tmp4 * tmp7 - sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)));
          fprime_of_xx0 = params->rho_slope + 2 * tmp4 * tmp5 * tmp6 * xx0 + tmp7 * (tmp0 * tmp2 + tmp0 * tmp3);
        }

        if (fprime_of_xx0 == (REAL)0.0) {
          break;
        }
        const REAL xx0_np1 = xx0 - f_of_xx0 / fprime_of_xx0;

        if (fabs(xx0 - xx0_np1) <= XX_TOLERANCE * params->dxx0 && fabs(f_of_xx0) <= F_OF_XX_TOLERANCE) {
          tolerance_has_been_met = 1;
        }
        xx0 = xx0_np1;
        iter++;
      } // END Newton-Raphson iterations to compute xx0
      if (iter >= ITER_MAX || !tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: Newton-Raphson failed for SinhCylindricalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0, (double)Cartx,
               (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: Newton-Raphson failed for SinhCylindricalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0, (double)Cartx,
                (double)Carty, (double)Cartz);
        exit(1);
#endif
      }
      xx[0] = xx0;
    }

    {
      int tolerance_has_been_met = 0;
      iter = 0;
      REAL xx2 = (REAL)0.5 * (params->xxmin2 + params->xxmax2);
      while (iter < ITER_MAX && !tolerance_has_been_met) {
        REAL f_of_xx2, fprime_of_xx2;

        {
          const REAL tmp0 = (1.0 / (params->SINHWZ));
          const REAL tmp6 = params->AMPLZ - params->z_slope;
          const REAL tmp5 = (1.0 / (exp(tmp0) - exp(-tmp0)));
          const REAL tmp2 = exp(tmp0 * xx2);
          const REAL tmp3 = exp(-tmp0 * xx2);
          const REAL tmp7 = tmp5 * tmp6 * ((xx2) * (xx2));
          const REAL tmp4 = tmp2 - tmp3;
          f_of_xx2 = -Cartz + params->z_slope * xx2 + tmp4 * tmp7;
          fprime_of_xx2 = params->z_slope + 2 * tmp4 * tmp5 * tmp6 * xx2 + tmp7 * (tmp0 * tmp2 + tmp0 * tmp3);
        }

        if (fprime_of_xx2 == (REAL)0.0) {
          break;
        }
        const REAL xx2_np1 = xx2 - f_of_xx2 / fprime_of_xx2;

        if (fabs(xx2 - xx2_np1) <= XX_TOLERANCE * params->dxx2 && fabs(f_of_xx2) <= F_OF_XX_TOLERANCE) {
          tolerance_has_been_met = 1;
        }
        xx2 = xx2_np1;
        iter++;
      } // END Newton-Raphson iterations to compute xx2
      if (iter >= ITER_MAX || !tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: Newton-Raphson failed for SinhCylindricalv2n2: xx2=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx2, (double)Cartx,
               (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: Newton-Raphson failed for SinhCylindricalv2n2: xx2=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx2, (double)Cartx,
                (double)Carty, (double)Cartz);
        exit(1);
#endif
      }
      xx[2] = xx2;
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
} // END FUNCTION Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindricalv2n2

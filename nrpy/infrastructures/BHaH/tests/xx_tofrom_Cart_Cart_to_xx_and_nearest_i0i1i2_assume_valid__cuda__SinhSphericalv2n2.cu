#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid. The index output may be NULL when only logical coordinates are needed.
 */
__host__ __device__ void Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__SinhSphericalv2n2(const params_struct *restrict params, const REAL xCart[3],
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
     *  Original SymPy expressions:
     *  "[xx[1] = acos(Cartz/sqrt(Cartx**2 + Carty**2 + Cartz**2))]"
     *  "[xx[2] = atan2(Carty, Cartx)]"
     */
    xx[1] = acos(Cartz / sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz))));
    xx[2] = atan2(Carty, Cartx);

    // Next perform safeguarded Newton-Raphson iterations as needed:
    const REAL F_OF_XX_TOLERANCE = 1e-12; // Tolerance of function for which we're finding the root.
    const int ITER_MAX = 100;
    int iter;

    {
      int tolerance_has_been_met = 0;
      iter = 0;
      REAL xx0_lo = params->xxmin0;
      REAL xx0_hi = params->xxmax0;
      REAL xx0 = xx0_lo;
      REAL f_of_xx0, fprime_of_xx0;

      {
        const REAL tmp0 = (1.0 / (params->SINHW));
        const REAL tmp6 = params->AMPL - params->r_slope;
        const REAL tmp5 = (1.0 / (exp(tmp0) - exp(-tmp0)));
        const REAL tmp2 = exp(tmp0 * xx0);
        const REAL tmp3 = exp(-tmp0 * xx0);
        const REAL tmp7 = tmp5 * tmp6 * ((xx0) * (xx0));
        const REAL tmp4 = tmp2 - tmp3;
        f_of_xx0 = params->r_slope * xx0 + tmp4 * tmp7 - sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
        fprime_of_xx0 = params->r_slope + 2 * tmp4 * tmp5 * tmp6 * xx0 + tmp7 * (tmp0 * tmp2 + tmp0 * tmp3);
      }

      REAL f_of_xx0_lo = f_of_xx0;

      xx0 = xx0_hi;

      {
        const REAL tmp0 = (1.0 / (params->SINHW));
        const REAL tmp6 = params->AMPL - params->r_slope;
        const REAL tmp5 = (1.0 / (exp(tmp0) - exp(-tmp0)));
        const REAL tmp2 = exp(tmp0 * xx0);
        const REAL tmp3 = exp(-tmp0 * xx0);
        const REAL tmp7 = tmp5 * tmp6 * ((xx0) * (xx0));
        const REAL tmp4 = tmp2 - tmp3;
        f_of_xx0 = params->r_slope * xx0 + tmp4 * tmp7 - sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
        fprime_of_xx0 = params->r_slope + 2 * tmp4 * tmp5 * tmp6 * xx0 + tmp7 * (tmp0 * tmp2 + tmp0 * tmp3);
      }

      REAL f_of_xx0_hi = f_of_xx0;

      if (fabs(f_of_xx0_lo) <= F_OF_XX_TOLERANCE) {
        xx0 = xx0_lo;
        tolerance_has_been_met = 1;
      } // END IF: lower coordinate is inverse root
      else if (fabs(f_of_xx0_hi) <= F_OF_XX_TOLERANCE) {
        xx0 = xx0_hi;
        tolerance_has_been_met = 1;
      } // END ELSE IF: upper coordinate is inverse root
      else {
        xx0 = (REAL)0.5 * (xx0_lo + xx0_hi);
        while (iter < ITER_MAX && !tolerance_has_been_met) {

          {
            const REAL tmp0 = (1.0 / (params->SINHW));
            const REAL tmp6 = params->AMPL - params->r_slope;
            const REAL tmp5 = (1.0 / (exp(tmp0) - exp(-tmp0)));
            const REAL tmp2 = exp(tmp0 * xx0);
            const REAL tmp3 = exp(-tmp0 * xx0);
            const REAL tmp7 = tmp5 * tmp6 * ((xx0) * (xx0));
            const REAL tmp4 = tmp2 - tmp3;
            f_of_xx0 = params->r_slope * xx0 + tmp4 * tmp7 - sqrt(((Cartx) * (Cartx)) + ((Carty) * (Carty)) + ((Cartz) * (Cartz)));
            fprime_of_xx0 = params->r_slope + 2 * tmp4 * tmp5 * tmp6 * xx0 + tmp7 * (tmp0 * tmp2 + tmp0 * tmp3);
          }

          if (fabs(f_of_xx0) <= F_OF_XX_TOLERANCE) {
            tolerance_has_been_met = 1;
          } // END IF: Newton iterate is inverse root
          else {
            if ((f_of_xx0_lo < (REAL)0.0 && f_of_xx0 < (REAL)0.0) || (f_of_xx0_lo > (REAL)0.0 && f_of_xx0 > (REAL)0.0)) {
              xx0_lo = xx0;
              f_of_xx0_lo = f_of_xx0;
            } // END IF: root lies above current iterate
            else {
              xx0_hi = xx0;
              f_of_xx0_hi = f_of_xx0;
            } // END ELSE: root lies below current iterate

            REAL xx0_trial = (REAL)0.5 * (xx0_lo + xx0_hi);
            if (fprime_of_xx0 != (REAL)0.0) {
              const REAL xx0_newton = xx0 - f_of_xx0 / fprime_of_xx0;
              if (xx0_newton > xx0_lo && xx0_newton < xx0_hi)
                xx0_trial = xx0_newton;
            } // END IF: Newton proposal has derivative

            xx0 = xx0_trial;
          } // END ELSE: update bracket and select trial
          iter++;
        } // END WHILE: safeguarded Newton iterations for xx0
      } // END ELSE: safeguarded Newton-Raphson with bisection fallback
      if (!tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: safeguarded Newton-Raphson failed for SinhSphericalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0, (double)Cartx,
               (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: safeguarded Newton-Raphson failed for SinhSphericalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
#endif
      } // END ELSE IF: safeguarded Newton-Raphson did not converge
      xx[0] = xx0;
    } // END BLOCK: inverse coordinate xx0

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
} // END FUNCTION: Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__SinhSphericalv2n2

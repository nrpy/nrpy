#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid
 */
void Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindricalv2n2(const params_struct *restrict params, const REAL xCart[3], REAL xx[3],
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
    // First compute analytical coordinate inversions:
    /*
     *  Original SymPy expression:
     *  "xx[1] = atan2(Carty, Cartx)"
     */
    xx[1] = atan2(Carty, Cartx);

    // Next perform safeguarded Newton-Raphson iterations as needed:
    const REAL F_OF_XX_TOLERANCE = 1e-12; // Tolerance of function for which we're finding the root.
    const int ITER_MAX = 100;
    int iter;

    {
      int tolerance_has_been_met = 0;
      int bracket_contains_root = 1;
      iter = 0;
      REAL xx0_lo = params->xxmin0;
      REAL xx0_hi = params->xxmax0;
      REAL xx0 = xx0_lo;
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

      REAL f_of_xx0_lo = f_of_xx0;

      xx0 = xx0_hi;

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

      REAL f_of_xx0_hi = f_of_xx0;

      if (fabs(f_of_xx0_lo) <= F_OF_XX_TOLERANCE) {
        xx0 = xx0_lo;
        tolerance_has_been_met = 1;
      } // END IF: lower physical-grid coordinate is the inverse root
      else if (fabs(f_of_xx0_hi) <= F_OF_XX_TOLERANCE) {
        xx0 = xx0_hi;
        tolerance_has_been_met = 1;
      } // END ELSE IF: upper physical-grid coordinate is the inverse root
      else if ((f_of_xx0_lo < (REAL)0.0 && f_of_xx0_hi < (REAL)0.0) || (f_of_xx0_lo > (REAL)0.0 && f_of_xx0_hi > (REAL)0.0)) {
        bracket_contains_root = 0;
      } // END ELSE IF: physical-grid coordinate range does not bracket inverse root
      else {
        xx0 = (REAL)0.5 * (xx0_lo + xx0_hi);
        while (iter < ITER_MAX && !tolerance_has_been_met) {

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

          if (fabs(f_of_xx0) <= F_OF_XX_TOLERANCE) {
            tolerance_has_been_met = 1;
          } // END IF: Newton iterate is the inverse root
          else {
            if ((f_of_xx0_lo < (REAL)0.0 && f_of_xx0 < (REAL)0.0) || (f_of_xx0_lo > (REAL)0.0 && f_of_xx0 > (REAL)0.0)) {
              xx0_lo = xx0;
              f_of_xx0_lo = f_of_xx0;
            } // END IF: inverse root lies above current iterate
            else {
              xx0_hi = xx0;
              f_of_xx0_hi = f_of_xx0;
            } // END ELSE: inverse root lies below current iterate

            REAL xx0_trial = (REAL)0.5 * (xx0_lo + xx0_hi);
            if (fprime_of_xx0 != (REAL)0.0) {
              const REAL xx0_newton = xx0 - f_of_xx0 / fprime_of_xx0;
              if (xx0_newton > xx0_lo && xx0_newton < xx0_hi)
                xx0_trial = xx0_newton;
            } // END IF: Newton proposal has a nonzero derivative

            xx0 = xx0_trial;
          } // END ELSE: update bracket and select safeguarded Newton trial
          iter++;
        } // END WHILE: safeguarded Newton-Raphson iterations to compute xx0
      } // END ELSE: physical-grid coordinate range brackets inverse root
      if (!bracket_contains_root) {
#ifdef __CUDA_ARCH__
        printf("ERROR: inverse root for SinhCylindricalv2n2 is outside xx0=[%.15e,%.15e]: f(min),f(max), x,y,z = %.15e %.15e %.15e %.15e %.15e\n",
               (double)params->xxmin0, (double)params->xxmax0, (double)f_of_xx0_lo, (double)f_of_xx0_hi, (double)Cartx, (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr,
                "ERROR: inverse root for SinhCylindricalv2n2 is outside xx0=[%.15e,%.15e]: f(min),f(max), x,y,z = %.15e %.15e %.15e %.15e %.15e\n",
                (double)params->xxmin0, (double)params->xxmax0, (double)f_of_xx0_lo, (double)f_of_xx0_hi, (double)Cartx, (double)Carty,
                (double)Cartz);
        exit(1);
#endif
      } // END IF: physical-grid coordinate range does not bracket inverse root
      else if (!tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: safeguarded Newton-Raphson failed for SinhCylindricalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0, (double)Cartx,
               (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: safeguarded Newton-Raphson failed for SinhCylindricalv2n2: xx0=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx0,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
#endif
      } // END ELSE IF: safeguarded Newton-Raphson did not converge
      xx[0] = xx0;
    } // END BLOCK: inverse coordinate xx0

    {
      int tolerance_has_been_met = 0;
      int bracket_contains_root = 1;
      iter = 0;
      REAL xx2_lo = params->xxmin2;
      REAL xx2_hi = params->xxmax2;
      REAL xx2 = xx2_lo;
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

      REAL f_of_xx2_lo = f_of_xx2;

      xx2 = xx2_hi;

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

      REAL f_of_xx2_hi = f_of_xx2;

      if (fabs(f_of_xx2_lo) <= F_OF_XX_TOLERANCE) {
        xx2 = xx2_lo;
        tolerance_has_been_met = 1;
      } // END IF: lower physical-grid coordinate is the inverse root
      else if (fabs(f_of_xx2_hi) <= F_OF_XX_TOLERANCE) {
        xx2 = xx2_hi;
        tolerance_has_been_met = 1;
      } // END ELSE IF: upper physical-grid coordinate is the inverse root
      else if ((f_of_xx2_lo < (REAL)0.0 && f_of_xx2_hi < (REAL)0.0) || (f_of_xx2_lo > (REAL)0.0 && f_of_xx2_hi > (REAL)0.0)) {
        bracket_contains_root = 0;
      } // END ELSE IF: physical-grid coordinate range does not bracket inverse root
      else {
        xx2 = (REAL)0.5 * (xx2_lo + xx2_hi);
        while (iter < ITER_MAX && !tolerance_has_been_met) {

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

          if (fabs(f_of_xx2) <= F_OF_XX_TOLERANCE) {
            tolerance_has_been_met = 1;
          } // END IF: Newton iterate is the inverse root
          else {
            if ((f_of_xx2_lo < (REAL)0.0 && f_of_xx2 < (REAL)0.0) || (f_of_xx2_lo > (REAL)0.0 && f_of_xx2 > (REAL)0.0)) {
              xx2_lo = xx2;
              f_of_xx2_lo = f_of_xx2;
            } // END IF: inverse root lies above current iterate
            else {
              xx2_hi = xx2;
              f_of_xx2_hi = f_of_xx2;
            } // END ELSE: inverse root lies below current iterate

            REAL xx2_trial = (REAL)0.5 * (xx2_lo + xx2_hi);
            if (fprime_of_xx2 != (REAL)0.0) {
              const REAL xx2_newton = xx2 - f_of_xx2 / fprime_of_xx2;
              if (xx2_newton > xx2_lo && xx2_newton < xx2_hi)
                xx2_trial = xx2_newton;
            } // END IF: Newton proposal has a nonzero derivative

            xx2 = xx2_trial;
          } // END ELSE: update bracket and select safeguarded Newton trial
          iter++;
        } // END WHILE: safeguarded Newton-Raphson iterations to compute xx2
      } // END ELSE: physical-grid coordinate range brackets inverse root
      if (!bracket_contains_root) {
#ifdef __CUDA_ARCH__
        printf("ERROR: inverse root for SinhCylindricalv2n2 is outside xx2=[%.15e,%.15e]: f(min),f(max), x,y,z = %.15e %.15e %.15e %.15e %.15e\n",
               (double)params->xxmin2, (double)params->xxmax2, (double)f_of_xx2_lo, (double)f_of_xx2_hi, (double)Cartx, (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr,
                "ERROR: inverse root for SinhCylindricalv2n2 is outside xx2=[%.15e,%.15e]: f(min),f(max), x,y,z = %.15e %.15e %.15e %.15e %.15e\n",
                (double)params->xxmin2, (double)params->xxmax2, (double)f_of_xx2_lo, (double)f_of_xx2_hi, (double)Cartx, (double)Carty,
                (double)Cartz);
        exit(1);
#endif
      } // END IF: physical-grid coordinate range does not bracket inverse root
      else if (!tolerance_has_been_met) {
#ifdef __CUDA_ARCH__
        printf("ERROR: safeguarded Newton-Raphson failed for SinhCylindricalv2n2: xx2=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx2, (double)Cartx,
               (double)Carty, (double)Cartz);
        asm("trap;");
#else
        fprintf(stderr, "ERROR: safeguarded Newton-Raphson failed for SinhCylindricalv2n2: xx2=%.15e, x,y,z = %.15e %.15e %.15e\n", (double)xx2,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
#endif
      } // END ELSE IF: safeguarded Newton-Raphson did not converge
      xx[2] = xx2;
    } // END BLOCK: inverse coordinate xx2

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
} // END FUNCTION: Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindricalv2n2

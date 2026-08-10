#include "BHaH_defines.h"

/**
 * Given Cartesian point (x,y,z), this function unshifts the grid back to the origin to output the corresponding
 *             (xx0,xx1,xx2) and the "closest" (i0,i1,i2) for the given grid. The index output may be NULL when only logical coordinates are needed.
 */
__host__ __device__ void Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_fisheyeN2(const params_struct *restrict params,
                                                                                               const REAL xCart[3], REAL xx[3],
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

    const REAL rCart = sqrt((Cartx) * (Cartx) + (Carty) * (Carty) + (Cartz) * (Cartz));
    if (!(isfinite(rCart))) {
      fprintf(stderr, "ERROR: bracketed inverse failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n", (double)rCart,
              (double)Cartx, (double)Carty, (double)Cartz);
      exit(1);
    } else if (rCart <= (REAL)0.0) {
      xx[0] = (REAL)0.0;
      xx[1] = (REAL)0.0;
      xx[2] = (REAL)0.0;
    } else {
      const REAL radial_scale = rCart;
      const REAL residual_tolerance = (REAL)1.0e-12 * radial_scale;
      REAL asymptotic_scale;
      asymptotic_scale = params->fisheye_a2 * params->fisheye_c;

      const REAL inv_asymptotic_scale = (fabs(asymptotic_scale) > (REAL)1.0e-15) ? (REAL)1.0 / asymptotic_scale : (REAL)1.0;
      REAL low = (REAL)0.0;
      REAL high = NRPYMAX(rCart * inv_asymptotic_scale, radial_scale);
      const REAL bracket_tolerance = (REAL)1.0e-12 * NRPYMAX(high, radial_scale);
      REAL radial_seed = (REAL)0.5 * high;
      int bracket_found = 0;
      int converged = 0;
      for (int expand = 0; expand < 80; expand++) {
        REAL high_map;
        const REAL high_tmp0 = (1.0 / (params->fisheye_s1));
        const REAL high_tmp1 = fmin(high, params->fisheye_R1);
        const REAL high_tmp4 = (1.0 / (params->fisheye_s2));
        const REAL high_tmp5 = fmin(high, params->fisheye_R2);
        const REAL high_tmp3 = exp(2 * high_tmp0 * high_tmp1 - 2 * high_tmp0 * fmax(high, params->fisheye_R1));
        const REAL high_tmp7 = exp(2 * high_tmp4 * high_tmp5 - 2 * high_tmp4 * fmax(high, params->fisheye_R2));
        high_map = params->fisheye_c * (high * params->fisheye_a2 +
                                        (1.0 / 2.0) * params->fisheye_s1 * (params->fisheye_a0 - params->fisheye_a1) *
                                            (2 * high_tmp0 * high_tmp1 + log1p(high_tmp3 * expm1(-4 * high_tmp0 * high_tmp1) / (high_tmp3 + 1))) /
                                            tanh(high_tmp0 * params->fisheye_R1) +
                                        (1.0 / 2.0) * params->fisheye_s2 * (params->fisheye_a1 - params->fisheye_a2) *
                                            (2 * high_tmp4 * high_tmp5 + log1p(high_tmp7 * expm1(-4 * high_tmp4 * high_tmp5) / (high_tmp7 + 1))) /
                                            tanh(high_tmp4 * params->fisheye_R2));

        const REAL high_residual = high_map - rCart;
        if (isfinite(high_residual) && high_residual >= (REAL)0.0) {
          bracket_found = 1;
          break;
        }
        high = NRPYMAX(high * (REAL)2.0, (REAL)1.0);
      } // END LOOP: for expand over bracket expansions
      if (!bracket_found) {
        fprintf(stderr, "ERROR: bracketed inverse failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n", (double)rCart,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
      }
      radial_seed = (REAL)0.5 * (low + high);
      for (int iter = 0; iter < 80; iter++) {
        REAL radial_map;
        REAL radial_map_prime;
        const REAL radial_tmp0 = fmin(params->fisheye_R1, radial_seed);
        const REAL radial_tmp1 = (1.0 / (params->fisheye_s1));
        const REAL radial_tmp5 = fmin(params->fisheye_R2, radial_seed);
        const REAL radial_tmp6 = (1.0 / (params->fisheye_s2));
        const REAL radial_tmp4 = (1.0 / 2.0) * (params->fisheye_a0 - params->fisheye_a1) / tanh(params->fisheye_R1 * radial_tmp1);
        const REAL radial_tmp9 = (1.0 / 2.0) * (params->fisheye_a1 - params->fisheye_a2) / tanh(params->fisheye_R2 * radial_tmp6);
        const REAL radial_tmp3 = exp(2 * radial_tmp0 * radial_tmp1 - 2 * radial_tmp1 * fmax(params->fisheye_R1, radial_seed));
        const REAL radial_tmp8 = exp(2 * radial_tmp5 * radial_tmp6 - 2 * radial_tmp6 * fmax(params->fisheye_R2, radial_seed));
        radial_map = params->fisheye_c *
                     (params->fisheye_a2 * radial_seed +
                      params->fisheye_s1 * radial_tmp4 *
                          (2 * radial_tmp0 * radial_tmp1 + log1p(radial_tmp3 * expm1(-4 * radial_tmp0 * radial_tmp1) / (radial_tmp3 + 1))) +
                      params->fisheye_s2 * radial_tmp9 *
                          (2 * radial_tmp5 * radial_tmp6 + log1p(radial_tmp8 * expm1(-4 * radial_tmp5 * radial_tmp6) / (radial_tmp8 + 1))));
        radial_map_prime =
            params->fisheye_c *
            (params->fisheye_a2 +
             radial_tmp4 * (-tanh(radial_tmp1 * (-params->fisheye_R1 + radial_seed)) + tanh(radial_tmp1 * (params->fisheye_R1 + radial_seed))) +
             radial_tmp9 * (-tanh(radial_tmp6 * (-params->fisheye_R2 + radial_seed)) + tanh(radial_tmp6 * (params->fisheye_R2 + radial_seed))));

        const REAL radial_residual = radial_map - rCart;
        REAL trial_seed = (REAL)0.5 * (low + high);
        if (isfinite(radial_map_prime) && fabs(radial_map_prime) > (REAL)1.0e-14) {
          const REAL newton_seed = radial_seed - radial_residual / radial_map_prime;
          if (isfinite(newton_seed) && newton_seed >= low && newton_seed <= high)
            trial_seed = newton_seed;
        }
        REAL trial_map;
        const REAL trial_tmp0 = (1.0 / (params->fisheye_s1));
        const REAL trial_tmp1 = fmin(params->fisheye_R1, trial_seed);
        const REAL trial_tmp4 = (1.0 / (params->fisheye_s2));
        const REAL trial_tmp5 = fmin(params->fisheye_R2, trial_seed);
        const REAL trial_tmp3 = exp(2 * trial_tmp0 * trial_tmp1 - 2 * trial_tmp0 * fmax(params->fisheye_R1, trial_seed));
        const REAL trial_tmp7 = exp(2 * trial_tmp4 * trial_tmp5 - 2 * trial_tmp4 * fmax(params->fisheye_R2, trial_seed));
        trial_map =
            params->fisheye_c * (params->fisheye_a2 * trial_seed +
                                 (1.0 / 2.0) * params->fisheye_s1 * (params->fisheye_a0 - params->fisheye_a1) *
                                     (2 * trial_tmp0 * trial_tmp1 + log1p(trial_tmp3 * expm1(-4 * trial_tmp0 * trial_tmp1) / (trial_tmp3 + 1))) /
                                     tanh(params->fisheye_R1 * trial_tmp0) +
                                 (1.0 / 2.0) * params->fisheye_s2 * (params->fisheye_a1 - params->fisheye_a2) *
                                     (2 * trial_tmp4 * trial_tmp5 + log1p(trial_tmp7 * expm1(-4 * trial_tmp4 * trial_tmp5) / (trial_tmp7 + 1))) /
                                     tanh(params->fisheye_R2 * trial_tmp4));

        REAL trial_residual = trial_map - rCart;
        if (!isfinite(trial_residual)) {
          trial_seed = (REAL)0.5 * (low + high);
          REAL trial_map_fallback;
          const REAL fallback_tmp0 = (1.0 / (params->fisheye_s1));
          const REAL fallback_tmp1 = fmin(params->fisheye_R1, trial_seed);
          const REAL fallback_tmp4 = (1.0 / (params->fisheye_s2));
          const REAL fallback_tmp5 = fmin(params->fisheye_R2, trial_seed);
          const REAL fallback_tmp3 = exp(2 * fallback_tmp0 * fallback_tmp1 - 2 * fallback_tmp0 * fmax(params->fisheye_R1, trial_seed));
          const REAL fallback_tmp7 = exp(2 * fallback_tmp4 * fallback_tmp5 - 2 * fallback_tmp4 * fmax(params->fisheye_R2, trial_seed));
          trial_map_fallback =
              params->fisheye_c *
              (params->fisheye_a2 * trial_seed +
               (1.0 / 2.0) * params->fisheye_s1 * (params->fisheye_a0 - params->fisheye_a1) *
                   (2 * fallback_tmp0 * fallback_tmp1 + log1p(fallback_tmp3 * expm1(-4 * fallback_tmp0 * fallback_tmp1) / (fallback_tmp3 + 1))) /
                   tanh(fallback_tmp0 * params->fisheye_R1) +
               (1.0 / 2.0) * params->fisheye_s2 * (params->fisheye_a1 - params->fisheye_a2) *
                   (2 * fallback_tmp4 * fallback_tmp5 + log1p(fallback_tmp7 * expm1(-4 * fallback_tmp4 * fallback_tmp5) / (fallback_tmp7 + 1))) /
                   tanh(fallback_tmp4 * params->fisheye_R2));

          trial_residual = trial_map_fallback - rCart;
          if (!isfinite(trial_residual)) {
            fprintf(stderr, "ERROR: bracketed inverse failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n",
                    (double)rCart, (double)Cartx, (double)Carty, (double)Cartz);
            exit(1);
          } // END IF: fallback residual is non-finite
        } // END IF: primary trial residual is non-finite
        if (trial_residual >= (REAL)0.0)
          high = trial_seed;
        else
          low = trial_seed;
        if (fabs(trial_residual) < residual_tolerance &&
            (fabs(high - low) < bracket_tolerance || fabs(trial_seed - radial_seed) < bracket_tolerance)) {
          radial_seed = trial_seed;
          converged = 1;
          break;
        }
        radial_seed = trial_seed;
      } // END LOOP: for iter over bracketed Newton iterations
      if (!converged || !isfinite(radial_seed) || radial_seed < (REAL)0.0) {
        fprintf(stderr, "ERROR: bracketed inverse failed for GeneralRFM_fisheyeN2 (fisheye): rCart, x,y,z = %.15e %.15e %.15e %.15e\n", (double)rCart,
                (double)Cartx, (double)Carty, (double)Cartz);
        exit(1);
      }
      const REAL scale = radial_seed / rCart;
      xx[0] = scale * Cartx;
      xx[1] = scale * Carty;
      xx[2] = scale * Cartz;
    } // END ELSE: invert fisheye radius away from the origin

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
} // END FUNCTION: Cart_to_xx_and_nearest_i0i1i2_assume_valid__rfm__GeneralRFM_fisheyeN2

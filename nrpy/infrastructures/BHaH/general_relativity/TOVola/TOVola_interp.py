"""
Register TOVola code TOVola_interp.c.

TOVola creates Tolman-Oppenheimer-Volkoff spherically symmetric initial data,
 typically for single neutron stars.

Authors: David Boyer
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""

import nrpy.c_function as cfc


def register_CFunction_TOVola_interp() -> None:
    """
    Register C function TOVola_interp().

    Provides spectral interpolator to provide data at arbitrary point x,y,z in Cartesian basis.
    """
    includes = ["BHaH_defines.h"]
    desc = "Provide high-order interpolation from TOVola grids onto an arbitrary point xCart[3] = {x,y,z} in the Spherical basis."
    prefunc = r"""
#ifndef TOVOLA_LAGRANGE_MAX_STENCIL_SIZE
#define TOVOLA_LAGRANGE_MAX_STENCIL_SIZE 64
#endif

#ifndef TOVOLA_LAGRANGE_COEFF_SUM_MAX
#define TOVOLA_LAGRANGE_COEFF_SUM_MAX ((REAL)1.0e6)
#endif

/**
 * Find the closest TOVola radial-table index to the requested isotropic radius.
 *
 * @param rr_iso Requested isotropic radius.
 * @param numpoints_arr Number of radial-table points.
 * @param[in] r_iso_arr Monotonic isotropic-radius table.
 * @return Closest radial-table index.
 */
static int TOVola_bisection_idx_finder(const REAL rr_iso, const int numpoints_arr, const REAL *restrict r_iso_arr) {
  int x1 = 0;
  int x2 = numpoints_arr - 1;
  REAL y1 = rr_iso - r_iso_arr[x1];
  REAL y2 = rr_iso - r_iso_arr[x2];
  if (y1 * y2 > 0) {
    fprintf(stderr, "INTERPOLATION BRACKETING ERROR: r_iso_min = %e ?<= r_iso = %.15e ?<= %e = r_iso_max\n", r_iso_arr[0], rr_iso,
            r_iso_arr[numpoints_arr - 1]);
    exit(EXIT_FAILURE);
  }
  for (int i = 0; i < numpoints_arr; i++) {
    int x_midpoint = (x1 + x2) / 2;
    REAL y_midpoint = rr_iso - r_iso_arr[x_midpoint];
    if (y_midpoint * y1 <= 0) {
      x2 = x_midpoint;
      y2 = y_midpoint;
    } else {
      x1 = x_midpoint;
      y1 = y_midpoint;
    }
    if (abs(x2 - x1) == 1) {
      // If r_iso_arr[x1] is closer to rr_iso than r_iso_arr[x2] then return x1:
      if (fabs(rr_iso - r_iso_arr[x1]) < fabs(rr_iso - r_iso_arr[x2])) {
        return x1;
      }
      // Otherwise return x2:
      return x2;
    } // END IF: bracket interval has collapsed to neighboring indices
  } // END LOOP: for i over all table points
  fprintf(stderr, "INTERPOLATION BRACKETING ERROR: r_iso_min = %e ?<= r_iso = %.15e ?<= %e = r_iso_max\n", r_iso_arr[0], rr_iso,
          r_iso_arr[numpoints_arr - 1]);
  exit(EXIT_FAILURE);
} // END FUNCTION: TOVola_bisection_idx_finder

/**
 * Compute one local Lagrange interpolation stencil and its amplification.
 *
 * @param rr_iso Requested isotropic radius.
 * @param numpoints_arr Number of radial-table points.
 * @param[in] r_iso_arr Monotonic isotropic-radius table.
 * @param idxmin Candidate stencil start.
 * @param stencil_size Candidate stencil size.
 * @param[out] coeffs_out Lagrange coefficients for the candidate stencil.
 * @param[out] coeff_sum_out Sum of absolute coefficient values.
 * @param[out] coeff_max_out Maximum absolute coefficient value.
 * @return Nonzero if the candidate is legal and finite; zero otherwise.
 */
static int TOVola_compute_lagrange_coefficients(const REAL rr_iso, const int numpoints_arr, const REAL *restrict r_iso_arr,
                                                const int idxmin, const int stencil_size, REAL *restrict coeffs_out,
                                                REAL *restrict coeff_sum_out, REAL *restrict coeff_max_out) {
  const int idxmax = idxmin + stencil_size - 1;
  if (idxmin < 0 || idxmax >= numpoints_arr || stencil_size < 3 || stencil_size > TOVOLA_LAGRANGE_MAX_STENCIL_SIZE) {
    return 0;
  } // END IF: candidate indices or stencil size are invalid

  if (!isfinite((double)rr_iso) || !isfinite((double)r_iso_arr[idxmin]) || !isfinite((double)r_iso_arr[idxmax]) ||
      r_iso_arr[idxmin] > rr_iso || r_iso_arr[idxmax] < rr_iso) {
    return 0;
  } // END IF: candidate does not bracket the interpolation radius

  for (int i = idxmin + 1; i <= idxmax; i++) {
    if (!isfinite((double)r_iso_arr[i]) || r_iso_arr[i] <= r_iso_arr[i - 1]) {
      return 0;
    } // END IF: candidate nodes are non-finite or not strictly increasing
  } // END LOOP: for i over candidate nodes

  *coeff_sum_out = 0.0;
  *coeff_max_out = 0.0;
  for (int i = 0; i < stencil_size; i++) {
    const REAL x_i = r_iso_arr[idxmin + i];
    REAL numer = 1.0;
    REAL denom = 1.0;
    for (int j = 0; j < stencil_size; j++) {
      if (j != i) {
        const REAL x_j = r_iso_arr[idxmin + j];
        const REAL denom_factor = x_i - x_j;
        if (denom_factor == 0.0) {
          return 0;
        } // END IF: candidate contains duplicate interpolation nodes
        numer *= rr_iso - x_j;
        denom *= denom_factor;
      } // END IF: skip the self Lagrange factor
    } // END LOOP: for j over candidate Lagrange denominator nodes

    coeffs_out[i] = numer / denom;
    if (!isfinite((double)coeffs_out[i])) {
      return 0;
    } // END IF: coefficient construction overflowed or produced NaN
    const REAL abs_coeff = fabs(coeffs_out[i]);
    *coeff_sum_out += abs_coeff;
    *coeff_max_out = NRPYMAX(*coeff_max_out, abs_coeff);
  } // END LOOP: for i over candidate Lagrange coefficients

  return isfinite((double)(*coeff_sum_out)) && isfinite((double)(*coeff_max_out));
} // END FUNCTION: TOVola_compute_lagrange_coefficients

/**
 * Choose a bounded-amplification local Lagrange interpolation stencil.
 *
 * @param rr_iso Requested isotropic radius.
 * @param[in] commondata Runtime parameter structure.
 * @param numpoints_arr Number of radial-table points.
 * @param[in] r_iso_arr Monotonic isotropic-radius table.
 * @param idx_mid Closest radial-table index returned by bisection.
 * @param[out] idxmin_out Selected stencil start.
 * @param[out] stencil_size_out Selected stencil size.
 * @param[out] coeffs_out Selected Lagrange coefficients.
 * @return Zero on success; nonzero on invalid configuration or unsafe interpolation.
 */
static int TOVola_choose_lagrange_stencil(const REAL rr_iso, const commondata_struct *restrict commondata, const int numpoints_arr,
                                          const REAL *restrict r_iso_arr, const int idx_mid, int *restrict idxmin_out,
                                          int *restrict stencil_size_out, REAL *restrict coeffs_out) {
  const int requested_stencil_size = commondata->interpolation_stencil_size;
  const int max_stencil_size = commondata->max_interpolation_stencil_size;
  if (numpoints_arr < 3 || requested_stencil_size < 3 || max_stencil_size < 3) {
    fprintf(stderr,
            "TOVola interpolation stencil error: invalid size metadata "
            "(rr_iso=%.17e idx_mid=%d requested=%d max=%d numpoints_arr=%d).\n",
            (double)rr_iso, idx_mid, requested_stencil_size, max_stencil_size, numpoints_arr);
    return 1;
  } // END IF: size metadata cannot support interpolation

  if (requested_stencil_size > max_stencil_size || max_stencil_size > TOVOLA_LAGRANGE_MAX_STENCIL_SIZE) {
    fprintf(stderr,
            "TOVola interpolation stencil error: requested size exceeds capacity "
            "(rr_iso=%.17e idx_mid=%d requested=%d max=%d compiled_max=%d numpoints_arr=%d).\n",
            (double)rr_iso, idx_mid, requested_stencil_size, max_stencil_size, TOVOLA_LAGRANGE_MAX_STENCIL_SIZE, numpoints_arr);
    return 1;
  } // END IF: requested or maximum stencil size is invalid

  if (requested_stencil_size % 2 == 0) {
    fprintf(stderr,
            "TOVola interpolation stencil error: interpolation_stencil_size must be odd "
            "(rr_iso=%.17e idx_mid=%d requested=%d max=%d numpoints_arr=%d).\n",
            (double)rr_iso, idx_mid, requested_stencil_size, max_stencil_size, numpoints_arr);
    return 1;
  } // END IF: requested stencil size is even

  int requested_size = requested_stencil_size;
  if (requested_size > numpoints_arr) {
    requested_size = numpoints_arr;
    if (requested_size % 2 == 0) {
      requested_size--;
    } // END IF: tiny table requires an odd clamped stencil size
  } // END IF: requested stencil size exceeds a tiny valid table

  int bracket_low = idx_mid;
  int bracket_high = idx_mid;
  if (r_iso_arr[idx_mid] < rr_iso && idx_mid + 1 < numpoints_arr) {
    bracket_high = idx_mid + 1;
  } else if (r_iso_arr[idx_mid] > rr_iso && idx_mid > 0) {
    bracket_low = idx_mid - 1;
  } // END ELSE IF: closest index lies above the interpolation radius

  REAL best_coeff_sum = -1.0;
  REAL trial_coeffs[TOVOLA_LAGRANGE_MAX_STENCIL_SIZE];
  for (int stencil_size = requested_size; stencil_size >= 3; stencil_size -= 2) {
    const int max_idxmin = numpoints_arr - stencil_size;
    int idxmin = NRPYMAX(0, idx_mid - stencil_size / 2 - 1);
    idxmin = NRPYMIN(idxmin, max_idxmin);
    if (idxmin > bracket_low) {
      idxmin = bracket_low;
    } // END IF: current-rule window misses the lower bracketing point
    if (idxmin + stencil_size - 1 < bracket_high) {
      idxmin = bracket_high - stencil_size + 1;
    } // END IF: current-rule window misses the upper bracketing point
    idxmin = NRPYMAX(0, NRPYMIN(idxmin, max_idxmin));

    REAL coeff_sum = 0.0;
    REAL coeff_max = 0.0;
    if (!TOVola_compute_lagrange_coefficients(rr_iso, numpoints_arr, r_iso_arr, idxmin, stencil_size, trial_coeffs, &coeff_sum,
                                              &coeff_max)) {
      continue;
    } // END IF: local candidate is not legal and finite

    if (best_coeff_sum < 0.0 || coeff_sum < best_coeff_sum) {
      best_coeff_sum = coeff_sum;
    } // END IF: local candidate is the least-amplifying stencil so far

    if (coeff_sum <= TOVOLA_LAGRANGE_COEFF_SUM_MAX) {
      *idxmin_out = idxmin;
      *stencil_size_out = stencil_size;
      for (int i = 0; i < stencil_size; i++) {
        coeffs_out[i] = trial_coeffs[i];
      } // END LOOP: for i over selected Lagrange coefficients
      return 0;
    } // END IF: local window is well-conditioned
  } // END LOOP: for stencil_size over odd fallback sizes

  fprintf(stderr,
          "TOVola interpolation stencil error: no safe local Lagrange stencil found "
          "(rr_iso=%.17e idx_mid=%d requested=%d max=%d numpoints_arr=%d best_coeff_sum=%.17e coeff_sum_max=%.17e).\n",
          (double)rr_iso, idx_mid, requested_stencil_size, max_stencil_size, numpoints_arr,
          (double)best_coeff_sum, (double)TOVOLA_LAGRANGE_COEFF_SUM_MAX);
  return 1;
} // END FUNCTION: TOVola_choose_lagrange_stencil

/**
 * Interpolate TOVola one-dimensional TOV data using a stable Lagrange stencil.
 *
 * @param rr_iso Requested isotropic radius.
 * @param[in] commondata Runtime parameter structure.
 * @param numpoints_arr Number of radial-table points.
 * @param[in] r_Schw_arr Schwarzschild-radius table.
 * @param[in] rho_energy_arr Energy-density table.
 * @param[in] rho_baryon_arr Baryon-density table.
 * @param[in] P_arr Pressure table.
 * @param[in] M_arr Enclosed-mass table.
 * @param[in] expnu_arr Lapse-related metric table.
 * @param[in] exp4phi_arr Conformal-factor table.
 * @param[in] r_iso_arr Isotropic-radius table.
 * @param[out] rho_energy Interpolated energy density.
 * @param[out] rho_baryon Interpolated baryon density.
 * @param[out] P Interpolated pressure.
 * @param[out] M Interpolated enclosed mass.
 * @param[out] expnu Interpolated lapse-related metric value.
 * @param[out] exp4phi Interpolated conformal factor.
 */
static void TOVola_TOV_interpolate_1D(REAL rr_iso, const commondata_struct *restrict commondata, const int numpoints_arr,
                                      const REAL *restrict r_Schw_arr, const REAL *restrict rho_energy_arr,
                                      const REAL *restrict rho_baryon_arr, const REAL *restrict P_arr, const REAL *restrict M_arr,
                                      const REAL *restrict expnu_arr, const REAL *restrict exp4phi_arr, const REAL *restrict r_iso_arr,
                                      REAL *restrict rho_energy, REAL *restrict rho_baryon, REAL *restrict P, REAL *restrict M,
                                      REAL *restrict expnu, REAL *restrict exp4phi) {
  if (numpoints_arr < 3) {
    fprintf(stderr, "TOVola interpolation error: numpoints_arr must be at least 3, got %d.\n", numpoints_arr);
    exit(EXIT_FAILURE);
  } // END IF: radial table is too small for interpolation

  const int R_idx = numpoints_arr - 1;
  const REAL M_star = M_arr[R_idx];
  const REAL r_iso_max_inside_star = r_iso_arr[R_idx];
  REAL r_Schw = 0.0;
  if (rr_iso < r_iso_max_inside_star) { // If we are INSIDE the star, we need to interpollate the data to the grid.
    // For this case, we know that for all our scalars, f(r) = f(-r)
    if (rr_iso < 0) {
      rr_iso = -rr_iso;
    } // END IF: reflect negative radius into the positive-radius table

    // First find the central interpolation stencil index:
    const int idx_mid = TOVola_bisection_idx_finder(rr_iso, numpoints_arr, r_iso_arr);

    int idxmin = 0;
    int stencil_size = 0;
    REAL l_i_of_r[TOVOLA_LAGRANGE_MAX_STENCIL_SIZE];
    if (TOVola_choose_lagrange_stencil(rr_iso, commondata, numpoints_arr, r_iso_arr, idx_mid, &idxmin, &stencil_size, l_i_of_r) != 0) {
      exit(EXIT_FAILURE);
    } // END IF: no stable interpolation stencil could be selected

    // Then perform the interpolation:
    *rho_energy = 0.0;
    *rho_baryon = 0.0;
    *P = 0.0;
    *M = 0.0;
    *expnu = 0.0;
    *exp4phi = 0.0;

    for (int i = idxmin; i < idxmin + stencil_size; i++) {
      r_Schw += l_i_of_r[i - idxmin] * r_Schw_arr[i];
      *rho_energy += l_i_of_r[i - idxmin] * rho_energy_arr[i];
      *rho_baryon += l_i_of_r[i - idxmin] * rho_baryon_arr[i];
      *P += l_i_of_r[i - idxmin] * P_arr[i];
      *M += l_i_of_r[i - idxmin] * M_arr[i];
      *expnu += l_i_of_r[i - idxmin] * expnu_arr[i];
      *exp4phi += l_i_of_r[i - idxmin] * exp4phi_arr[i];
    } // END LOOP: for i over selected TOVola interpolation stencil
  } else {
    // If we are OUTSIDE the star, the solution is just Schwarzschild.
    r_Schw = (rr_iso + M_star) + M_star * M_star / (4.0 * rr_iso); // Need to know what r_Schw is at our current grid location.
    *rho_energy = 0;
    *rho_baryon = 0;
    *P = 0;
    *M = M_star;
    *expnu = 1. - 2.0 * (M_star) / r_Schw;
    *exp4phi = (r_Schw * r_Schw) / (rr_iso * rr_iso);
  } // END ELSE: point lies outside the stellar surface
  //printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e hhhh\n", rr_iso, r_Schw, *rho_energy, *rho_baryon, *P, *M, *expnu, *exp4phi);
} // END FUNCTION: TOVola_TOV_interpolate_1D
"""
    name = "TOVola_interp"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
     const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"""
    body = r"""
  const REAL x = xCart[0];
  const REAL y = xCart[1];
  const REAL z = xCart[2];
  const REAL r_iso = sqrt(x * x + y * y + z * z);
  // self.xxSph[1] = sp.acos(self.xx[2] / self.xxSph[0])
  const REAL theta = acos(z / r_iso);

  // Perform pointwise interpolation to radius r using ID_persist data
  REAL rho_energy_val, rho_baryon_val, P_val, M_val, expnu_val, exp4phi_val;
  TOVola_TOV_interpolate_1D(r_iso, commondata, ID_persist->numpoints_arr, ID_persist->r_Schw_arr, ID_persist->rho_energy_arr,
                            ID_persist->rho_baryon_arr, ID_persist->P_arr, ID_persist->M_arr, ID_persist->expnu_arr, ID_persist->exp4phi_arr,
                            ID_persist->r_iso_arr, &rho_energy_val, &rho_baryon_val, &P_val, &M_val, &expnu_val, &exp4phi_val);

  // Assign interpolated values to initial_data_struct
  initial_data->alpha = sqrt(expnu_val);

  // Assuming beta and B fields are zero in this context
  initial_data->betaSphorCartU0 = 0.0;
  initial_data->betaSphorCartU1 = 0.0;
  initial_data->betaSphorCartU2 = 0.0;
  initial_data->BSphorCartU0 = 0.0;
  initial_data->BSphorCartU1 = 0.0;
  initial_data->BSphorCartU2 = 0.0;

  // Metric components (assuming diagonal for simplicity)
  initial_data->gammaSphorCartDD00 = exp4phi_val;
  initial_data->gammaSphorCartDD01 = 0.0;
  initial_data->gammaSphorCartDD02 = 0.0;
  initial_data->gammaSphorCartDD11 = exp4phi_val * r_iso * r_iso;
  initial_data->gammaSphorCartDD12 = 0.0;
  initial_data->gammaSphorCartDD22 = exp4phi_val * r_iso * r_iso * sin(theta) * sin(theta);

  // Extrinsic curvature components set to zero
  initial_data->KSphorCartDD00 = 0.0;
  initial_data->KSphorCartDD01 = 0.0;
  initial_data->KSphorCartDD02 = 0.0;
  initial_data->KSphorCartDD11 = 0.0;
  initial_data->KSphorCartDD12 = 0.0;
  initial_data->KSphorCartDD22 = 0.0;

  initial_data->T4SphorCartUU00 = rho_energy_val / expnu_val;
  initial_data->T4SphorCartUU01 = 0.0;
  initial_data->T4SphorCartUU02 = 0.0;
  initial_data->T4SphorCartUU03 = 0.0;
  initial_data->T4SphorCartUU11 = P_val / exp4phi_val;
  initial_data->T4SphorCartUU12 = 0.0;
  initial_data->T4SphorCartUU13 = 0.0;
  initial_data->T4SphorCartUU22 = P_val / (exp4phi_val * r_iso * r_iso);
  initial_data->T4SphorCartUU23 = 0.0;
  initial_data->T4SphorCartUU33 = P_val / (exp4phi_val * r_iso * r_iso * sin(theta) * sin(theta));
"""

    cfc.register_CFunction(
        subdirectory="TOVola",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )

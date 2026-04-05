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
/* Bisection index finder using binary search */
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
    }
  }
  fprintf(stderr, "INTERPOLATION BRACKETING ERROR: r_iso_min = %e ?<= r_iso = %.15e ?<= %e = r_iso_max\n", r_iso_arr[0], rr_iso,
          r_iso_arr[numpoints_arr - 1]);
  exit(EXIT_FAILURE);
}

/* Interpolation Function using Lagrange Polynomial */
static void TOVola_TOV_interpolate_1D(REAL rr_iso, const commondata_struct *restrict commondata,
                                      const int interpolation_stencil_size, const int numpoints_arr, const REAL *restrict r_Schw_arr,
                                      const REAL *restrict rho_energy_arr, const REAL *restrict rho_baryon_arr, const REAL *restrict P_arr,
                                      const REAL *restrict M_arr, const REAL *restrict expnu_arr, const REAL *restrict exp4phi_arr,
                                      const REAL *restrict r_iso_arr, REAL *restrict rho_energy, REAL *restrict rho_baryon, REAL *restrict P,
                                      REAL *restrict M, REAL *restrict expnu, REAL *restrict exp4phi) {
  const int R_idx = numpoints_arr - 1;
  const REAL M_star = M_arr[R_idx];
  const REAL r_iso_max_inside_star = r_iso_arr[R_idx];
  REAL r_Schw = 0.0;
  if (rr_iso < r_iso_max_inside_star) { // If we are INSIDE the star, we need to interpollate the data to the grid.
    // For this case, we know that for all our scalars, f(r) = f(-r)
    if (rr_iso < 0)
      rr_iso = -rr_iso;

    // First find the central interpolation stencil index:
    int idx_mid = TOVola_bisection_idx_finder(rr_iso, numpoints_arr, r_iso_arr);

    /* Use standard library functions instead of redefining macros */
    int idxmin = MAX(0, idx_mid - commondata->interpolation_stencil_size / 2 - 1);

    // -= Do not allow the interpolation stencil to cross the star's surface =-
    // max index is when idxmin + (commondata->interpolation_stencil_size-1) = R_idx
    //  -> idxmin at most can be R_idx - commondata->interpolation_stencil_size + 1
    idxmin = MIN(idxmin, R_idx - commondata->interpolation_stencil_size + 1);

    // Ensure that commondata->interpolation_stencil_size does not exceed the maximum
    if (commondata->interpolation_stencil_size > commondata->max_interpolation_stencil_size) {
      fprintf(stderr, "Interpolation stencil size exceeds maximum allowed.\n");
      exit(EXIT_FAILURE);
    }

    // Now perform the Lagrange polynomial interpolation:

    // First compute the interpolation coefficients:
    REAL r_iso_sample[commondata->max_interpolation_stencil_size];
    for (int i = idxmin; i < idxmin + commondata->interpolation_stencil_size; i++) {
      //if(i < 0 || i >= R_idx-1) { fprintf(stderr, "ERROR!\n"); exit(1); }
      r_iso_sample[i - idxmin] = r_iso_arr[i];
    }
    REAL l_i_of_r[commondata->max_interpolation_stencil_size];
    for (int i = 0; i < commondata->interpolation_stencil_size; i++) {
      REAL numer = 1.0;
      REAL denom = 1.0;
      for (int j = 0; j < commondata->interpolation_stencil_size; j++) {
        if (j != i) {
          numer *= (rr_iso - r_iso_sample[j]);
          denom *= (r_iso_sample[i] - r_iso_sample[j]);
        }
      }
      l_i_of_r[i] = numer / denom;
    }

    // Then perform the interpolation:
    *rho_energy = 0.0;
    *rho_baryon = 0.0;
    *P = 0.0;
    *M = 0.0;
    *expnu = 0.0;
    *exp4phi = 0.0;

    for (int i = idxmin; i < idxmin + commondata->interpolation_stencil_size; i++) {
      r_Schw += l_i_of_r[i - idxmin] * r_Schw_arr[i];
      *rho_energy += l_i_of_r[i - idxmin] * rho_energy_arr[i];
      *rho_baryon += l_i_of_r[i - idxmin] * rho_baryon_arr[i];
      *P += l_i_of_r[i - idxmin] * P_arr[i];
      *M += l_i_of_r[i - idxmin] * M_arr[i];
      *expnu += l_i_of_r[i - idxmin] * expnu_arr[i];
      *exp4phi += l_i_of_r[i - idxmin] * exp4phi_arr[i];
    }

  } else {
    // If we are OUTSIDE the star, the solution is just Schwarzschild.
    r_Schw = (rr_iso + M_star) + M_star * M_star / (4.0 * rr_iso); // Need to know what r_Schw is at our current grid location.
    *rho_energy = 0;
    *rho_baryon = 0;
    *P = 0;
    *M = M_star;
    *expnu = 1. - 2.0 * (M_star) / r_Schw;
    *exp4phi = (r_Schw * r_Schw) / (rr_iso * rr_iso);
  }
  //printf("%.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e hhhh\n", rr_iso, r_Schw, *rho_energy, *rho_baryon, *P, *M, *expnu, *exp4phi);
}
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
  TOVola_TOV_interpolate_1D(r_iso, commondata, commondata->max_interpolation_stencil_size, ID_persist->numpoints_arr, ID_persist->r_Schw_arr,
                            ID_persist->rho_energy_arr, ID_persist->rho_baryon_arr, ID_persist->P_arr, ID_persist->M_arr, ID_persist->expnu_arr,
                            ID_persist->exp4phi_arr, ID_persist->r_iso_arr, &rho_energy_val, &rho_baryon_val, &P_val, &M_val, &expnu_val,
                            &exp4phi_val);

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

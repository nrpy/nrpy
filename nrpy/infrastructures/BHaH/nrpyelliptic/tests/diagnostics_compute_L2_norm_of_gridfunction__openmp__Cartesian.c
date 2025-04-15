#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Kernel: compute_L2_norm_of_gridfunction_host.
 * Kernel to compute L2 quantities pointwise (summed).
 */
static void compute_L2_norm_of_gridfunction_host(const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1,
                                                 const REAL *restrict x2, const REAL *restrict in_gfs, REAL *restrict squared_sum_final,
                                                 REAL *restrict volume_sum_final, const REAL integration_radius, const int gf_index) {
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
  MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
  MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;

  // Load necessary parameters from params_struct
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;

  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum = 0.0;

#pragma omp parallel for reduction(+ : squared_sum, volume_sum)
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = x2[i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = x1[i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = x0[i0];

        /*
         *  Original SymPy expressions:
         *  "[const DOUBLE r = sqrt(xx0**2 + xx1**2 + xx2**2)]"
         *  "[const DOUBLE sqrtdetgamma = 1]"
         */
        const DOUBLE r = sqrt(((xx0) * (xx0)) + ((xx1) * (xx1)) + ((xx2) * (xx2)));
        const DOUBLE sqrtdetgamma = 1;

        if (r < integration_radius) {
          const DOUBLE gf_of_x = in_gfs[IDX4(gf_index, i0, i1, i2)];
          const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;

          squared_sum += gf_of_x * gf_of_x * dV;
          volume_sum += dV;
        } // END if(r < integration_radius)

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
  *squared_sum_final = squared_sum;
  *volume_sum_final = volume_sum;
} // END FUNCTION compute_L2_norm_of_gridfunction_host

/**
 * Compute l2-norm of a gridfunction assuming a single grid.
 */
void compute_L2_norm_of_gridfunction__rfm__Cartesian(commondata_struct *restrict commondata, params_struct *restrict params, REAL *restrict xx[3],
                                                     const REAL integration_radius, const int gf_index, REAL *l2norm, const REAL *restrict in_gfs) {
#include "../set_CodeParameters.h"

  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS_tot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;
  REAL *restrict x0 = xx[0];
  REAL *restrict x1 = xx[1];
  REAL *restrict x2 = xx[2];

  // Set summation variables to compute l2-norm
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum = 0.0;
  compute_L2_norm_of_gridfunction_host(params, x0, x1, x2, in_gfs, &squared_sum, &volume_sum, integration_radius, gf_index);

  // Compute and output the log of the l2-norm.
  *l2norm = log10(1e-16 + sqrt(squared_sum / volume_sum)); // 1e-16 + ... avoids log10(0)
} // END FUNCTION compute_L2_norm_of_gridfunction__rfm__Cartesian

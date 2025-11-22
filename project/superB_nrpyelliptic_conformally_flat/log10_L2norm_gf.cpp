#include "BHaH_defines.h"

/**
 * Compute l2-norm of a gridfunction assuming a single grid.
 */
void log10_L2norm_gf(commondata_struct *restrict commondata, griddata_struct *restrict griddata, const REAL integration_radius, const int gf_index,
                     const REAL *restrict in_gf, REAL localsums_for_residualH[2]) {
  // Unpack grid parameters assuming a single grid
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Define reference metric grid
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Set summation variables to compute l2-norm
  REAL squared_sum = 0.0;
  REAL volume_sum = 0.0;

#pragma omp parallel for reduction(+ : squared_sum, volume_sum)
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    MAYBE_UNUSED const REAL xx2 = xx[2][i2];
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      MAYBE_UNUSED const REAL xx1 = xx[1][i1];
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
        MAYBE_UNUSED const REAL xx0 = xx[0][i0];

        /*
         *  Original SymPy expressions:
         *  "[const REAL r = AMPL*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1/SINHW) - exp(-1/SINHW))]"
         *  "[const REAL sqrtdetgamma = AMPL**6*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)**2*(exp(xx0/SINHW) -
         * exp(-xx0/SINHW))**4*sin(xx1)**2/(exp(1/SINHW) - exp(-1/SINHW))**6]"
         */
        const REAL tmp0 = (1.0 / (SINHW));
        const REAL tmp1 = exp(tmp0) - exp(-tmp0);
        const REAL tmp3 = exp(tmp0 * xx0);
        const REAL tmp4 = exp(-tmp0 * xx0);
        const REAL tmp5 = tmp3 - tmp4;
        const REAL r = AMPL * tmp5 / tmp1;
        const REAL sqrtdetgamma = pow(AMPL, 6) * ((tmp5) * (tmp5) * (tmp5) * (tmp5)) * ((tmp0 * tmp3 + tmp0 * tmp4) * (tmp0 * tmp3 + tmp0 * tmp4)) *
                                  ((sin(xx1)) * (sin(xx1))) / pow(tmp1, 6);

        if (r < integration_radius) {
          const REAL gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
          const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
          squared_sum += gf_of_x * gf_of_x * dV;
          volume_sum += dV;
        } // END if(r < integration_radius)

      } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
    } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
  } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)

  localsums_for_residualH[0] = squared_sum;
  localsums_for_residualH[1] = volume_sum;
} // END FUNCTION log10_L2norm_gf

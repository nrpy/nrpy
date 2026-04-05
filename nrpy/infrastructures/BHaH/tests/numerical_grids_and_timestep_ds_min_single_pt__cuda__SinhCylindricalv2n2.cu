#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min.
 */
__host__ __device__ void ds_min_single_pt__rfm__SinhCylindricalv2n2(const params_struct *restrict params, const REAL xx0, const REAL xx1,
                                                                    const REAL xx2, REAL *restrict ds_min) {
  const REAL AMPLRHO = params->AMPLRHO;
  const REAL AMPLZ = params->AMPLZ;
  const REAL SINHWRHO = params->SINHWRHO;
  const REAL SINHWZ = params->SINHWZ;
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;
  const REAL rho_slope = params->rho_slope;
  const REAL z_slope = params->z_slope;
  /*
   *  Original SymPy expressions:
   *  "[const REAL ds0 = Abs(dxx0*(rho_slope + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO)/SINHWRHO +
   * exp(-xx0/SINHWRHO)/SINHWRHO)/(exp(1/SINHWRHO) - exp(-1/SINHWRHO)) + 2*xx0*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) -
   * exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) - exp(-1/SINHWRHO))))]"
   *  "[const REAL ds1 = Abs(dxx1*(rho_slope*xx0 + xx0**2*(AMPLRHO - rho_slope)*(exp(xx0/SINHWRHO) - exp(-xx0/SINHWRHO))/(exp(1/SINHWRHO) -
   * exp(-1/SINHWRHO))))]"
   *  "[const REAL ds2 = Abs(dxx2*(xx2**2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ)/SINHWZ + exp(-xx2/SINHWZ)/SINHWZ)/(exp(1/SINHWZ) - exp(-1/SINHWZ)) +
   * 2*xx2*(AMPLZ - z_slope)*(exp(xx2/SINHWZ) - exp(-xx2/SINHWZ))/(exp(1/SINHWZ) - exp(-1/SINHWZ)) + z_slope))]"
   */
  const REAL tmp0 = (1.0 / (SINHWRHO));
  const REAL tmp7 = (1.0 / (SINHWZ));
  const REAL tmp5 = (AMPLRHO - rho_slope) / (exp(tmp0) - exp(-tmp0));
  const REAL tmp11 = (AMPLZ - z_slope) / (exp(tmp7) - exp(-tmp7));
  const REAL tmp2 = exp(tmp0 * xx0);
  const REAL tmp3 = exp(-tmp0 * xx0);
  const REAL tmp6 = tmp5 * ((xx0) * (xx0));
  const REAL tmp9 = exp(tmp7 * xx2);
  const REAL tmp10 = exp(-tmp7 * xx2);
  const REAL tmp4 = tmp2 - tmp3;
  const REAL ds0 = fabs(dxx0 * (rho_slope + 2 * tmp4 * tmp5 * xx0 + tmp6 * (tmp0 * tmp2 + tmp0 * tmp3)));
  const REAL ds1 = fabs(dxx1 * (rho_slope * xx0 + tmp4 * tmp6));
  const REAL ds2 = fabs(dxx2 * (tmp11 * ((xx2) * (xx2)) * (tmp10 * tmp7 + tmp7 * tmp9) - 2 * tmp11 * xx2 * (tmp10 - tmp9) + z_slope));
  *ds_min = MIN(ds0, MIN(ds1, ds2));
} // END FUNCTION ds_min_single_pt__rfm__SinhCylindricalv2n2

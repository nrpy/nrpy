#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min.
 */
void ds_min_single_pt__rfm__SinhSymTP(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min) {
  const REAL AMAX = params->AMAX;
  const REAL SINHWAA = params->SINHWAA;
  const REAL bScale = params->bScale;
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;
  /*
   *  Original SymPy expressions:
   *  "[const REAL ds0 = Abs(AMAX*dxx0*(exp(xx0/SINHWAA)/SINHWAA + exp(-xx0/SINHWAA)/SINHWAA)*sqrt(AMAX**2*(exp(xx0/SINHWAA) -
   * exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2*sin(xx1)**2)/(exp(1/SINHWAA) -
   * exp(-1/SINHWAA)))/Abs(sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 + bScale**2))]"
   *  "[const REAL ds1 = Abs(dxx1*sqrt(AMAX**2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))**2/(exp(1/SINHWAA) - exp(-1/SINHWAA))**2 +
   * bScale**2*sin(xx1)**2))]"
   *  "[const REAL ds2 = Abs(AMAX*dxx2*(exp(xx0/SINHWAA) - exp(-xx0/SINHWAA))*sin(xx1)/(exp(1/SINHWAA) - exp(-1/SINHWAA)))]"
   */
  const REAL tmp1 = (1.0 / (SINHWAA));
  const REAL tmp8 = sin(xx1);
  const REAL tmp2 = exp(tmp1) - exp(-tmp1);
  const REAL tmp4 = exp(tmp1 * xx0);
  const REAL tmp5 = exp(-tmp1 * xx0);
  const REAL tmp10 = AMAX / tmp2;
  const REAL tmp6 = tmp4 - tmp5;
  const REAL tmp7 = ((AMAX) * (AMAX)) * ((tmp6) * (tmp6)) / ((tmp2) * (tmp2));
  const REAL tmp9 = sqrt(((bScale) * (bScale)) * ((tmp8) * (tmp8)) + tmp7);
  const REAL ds0 = fabs(dxx0 * tmp10 * tmp9 * (tmp1 * tmp4 + tmp1 * tmp5)) / fabs(sqrt(((bScale) * (bScale)) + tmp7));
  const REAL ds1 = fabs(dxx1 * tmp9);
  const REAL ds2 = fabs(dxx2 * tmp10 * tmp6 * tmp8);
  *ds_min = MIN(ds0, MIN(ds1, ds2));
} // END FUNCTION ds_min_single_pt__rfm__SinhSymTP

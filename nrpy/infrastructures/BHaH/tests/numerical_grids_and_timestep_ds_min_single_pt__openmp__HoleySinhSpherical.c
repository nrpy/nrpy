#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min.
 */
void ds_min_single_pt__rfm__HoleySinhSpherical(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2,
                                               REAL *restrict ds_min) {
  const REAL AMPL = params->AMPL;
  const REAL SINHW = params->SINHW;
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;
  /*
   *  Original SymPy expressions:
   *  "[const REAL ds0 = Abs(AMPL*dxx0*(exp(xx0/SINHW)/SINHW + exp(-xx0/SINHW)/SINHW)/(exp(1/SINHW) - exp(-1/SINHW)))]"
   *  "[const REAL ds1 = Abs(AMPL*dxx1*(exp(xx0/SINHW) - exp(-xx0/SINHW))/(exp(1/SINHW) - exp(-1/SINHW)))]"
   *  "[const REAL ds2 = Abs(AMPL*dxx2*(exp(xx0/SINHW) - exp(-xx0/SINHW))*sin(xx1)/(exp(1/SINHW) - exp(-1/SINHW)))]"
   */
  const REAL tmp0 = (1.0 / (SINHW));
  const REAL tmp4 = AMPL / (exp(tmp0) - exp(-tmp0));
  const REAL tmp2 = exp(tmp0 * xx0);
  const REAL tmp3 = exp(-tmp0 * xx0);
  const REAL tmp5 = tmp4 * (tmp2 - tmp3);
  const REAL ds0 = fabs(dxx0 * tmp4 * (tmp0 * tmp2 + tmp0 * tmp3));
  const REAL ds1 = fabs(dxx1 * tmp5);
  const REAL ds2 = fabs(dxx2 * tmp5 * sin(xx1));
  *ds_min = MIN(ds0, MIN(ds1, ds2));
} // END FUNCTION ds_min_single_pt__rfm__HoleySinhSpherical

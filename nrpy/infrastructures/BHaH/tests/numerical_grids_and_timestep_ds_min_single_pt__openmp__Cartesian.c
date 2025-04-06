#include "../BHaH_defines.h"
#include "../BHaH_function_prototypes.h"
/**
 * Examining all three directions at a given point on a numerical grid, find the minimum grid spacing ds_min.
 */
void ds_min_single_pt__rfm__Cartesian(const params_struct *restrict params, const REAL xx0, const REAL xx1, const REAL xx2, REAL *restrict ds_min) {
  const REAL dxx0 = params->dxx0;
  const REAL dxx1 = params->dxx1;
  const REAL dxx2 = params->dxx2;
  /*
   *  Original SymPy expressions:
   *  "[const REAL ds0 = Abs(dxx0)]"
   *  "[const REAL ds1 = Abs(dxx1)]"
   *  "[const REAL ds2 = Abs(dxx2)]"
   */
  const REAL ds0 = fabs(dxx0);
  const REAL ds1 = fabs(dxx1);
  const REAL ds2 = fabs(dxx2);
  *ds_min = MIN(ds0, MIN(ds1, ds2));
} // END FUNCTION ds_min_single_pt__rfm__Cartesian

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Compute quasiKinnersley tetrad for psi4, with use_metric_to_construct_unit_normal=False
 */
void psi4_tetrad(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL cf, const REAL hDD00,
                 const REAL hDD01, const REAL hDD02, const REAL hDD11, const REAL hDD12, const REAL hDD22, REAL *mre4U0, REAL *mre4U1, REAL *mre4U2,
                 REAL *mre4U3, REAL *mim4U0, REAL *mim4U1, REAL *mim4U2, REAL *mim4U3, REAL *n4U0, REAL *n4U1, REAL *n4U2, REAL *n4U3,
                 REAL *restrict xx[3], const int i0, const int i1, const int i2) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    psi4_tetrad__rfm__SinhSpherical(commondata, params, cf, hDD00, hDD01, hDD02, hDD11, hDD12, hDD22, mre4U0, mre4U1, mre4U2, mre4U3, mim4U0, mim4U1,
                                    mim4U2, mim4U3, n4U0, n4U1, n4U2, n4U3, xx, i0, i1, i2);
    break;
  default:
    fprintf(stderr, "ERROR in psi4_tetrad(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

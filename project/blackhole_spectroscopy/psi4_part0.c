#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Compute psi4 at all interior gridpoints, part 0
 */
void psi4_part0(const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs,
                REAL *restrict diagnostic_output_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    psi4_part0__rfm__SinhSpherical(commondata, params, xx, in_gfs, diagnostic_output_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in psi4_part0(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

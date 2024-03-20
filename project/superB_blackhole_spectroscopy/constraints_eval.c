#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Evaluate BSSN constraints.
 */
void constraints_eval(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                      const REAL *restrict in_gfs, const REAL *restrict auxevol_gfs, REAL *restrict diagnostic_output_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    constraints_eval__rfm__SinhSpherical(commondata, params, rfmstruct, in_gfs, auxevol_gfs, diagnostic_output_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in constraints_eval(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

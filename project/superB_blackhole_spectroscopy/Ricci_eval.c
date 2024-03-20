#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Set Ricci tensor.
 */
void Ricci_eval(const commondata_struct *restrict commondata, const params_struct *restrict params, const rfm_struct *restrict rfmstruct,
                const REAL *restrict in_gfs, REAL *restrict auxevol_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    Ricci_eval__rfm__SinhSpherical(commondata, params, rfmstruct, in_gfs, auxevol_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in Ricci_eval(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * Enforce det(gammabar) = det(gammahat) constraint. Required for strong hyperbolicity.
 */
void enforce_detgammabar_equals_detgammahat(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                            const rfm_struct *restrict rfmstruct, REAL *restrict in_gfs) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    enforce_detgammabar_equals_detgammahat__rfm__SinhSpherical(commondata, params, rfmstruct, in_gfs);
    break;
  default:
    fprintf(stderr, "ERROR in enforce_detgammabar_equals_detgammahat(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

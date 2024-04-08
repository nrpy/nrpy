#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/*
 * rfm_precompute_defines: reference metric precomputed lookup arrays: defines
 */
void rfm_precompute_defines(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct,
                            REAL *restrict xx[3]) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    rfm_precompute_defines__rfm__SinhSpherical(commondata, params, rfmstruct, xx);
    break;
  default:
    fprintf(stderr, "ERROR in rfm_precompute_defines(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
}

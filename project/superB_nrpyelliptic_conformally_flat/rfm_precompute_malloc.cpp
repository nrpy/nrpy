#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc(const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct) {
  switch (params->CoordSystem_hash) {
  case SINHSPHERICAL:
    rfm_precompute_malloc__rfm__SinhSpherical(commondata, params, rfmstruct);
    break;
  default:
    fprintf(stderr, "ERROR in rfm_precompute_malloc(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
    exit(1);
  }
} // END FUNCTION rfm_precompute_malloc

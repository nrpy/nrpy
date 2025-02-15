#include "../BHaH_defines.h"
/**
 * GPU Kernel: rfm_precompute_free__deallocate.
 * Kernel to deallocate rfmstruct arrays.
 */
static void rfm_precompute_free__deallocate(rfm_struct *restrict rfmstruct) {
  // Temporary parameters
  free(rfmstruct->f0_of_xx0);
} // END FUNCTION rfm_precompute_free__deallocate

/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__Cylindrical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                           rfm_struct *restrict rfmstruct) {
  rfm_precompute_free__deallocate(rfmstruct);
} // END FUNCTION rfm_precompute_free__rfm__Cylindrical

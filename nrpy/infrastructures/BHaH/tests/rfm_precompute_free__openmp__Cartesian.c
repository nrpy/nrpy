#include "../BHaH_defines.h"
/**
 * GPU Kernel: rfm_precompute_free__deallocate.
 * Kernel to deallocate rfmstruct arrays.
 */
static void rfm_precompute_free__deallocate(rfm_struct *restrict rfmstruct) {
  // Temporary parameters
} // END FUNCTION rfm_precompute_free__deallocate

/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__Cartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                         rfm_struct *restrict rfmstruct) {
  rfm_precompute_free__deallocate(rfmstruct);
} // END FUNCTION rfm_precompute_free__rfm__Cartesian

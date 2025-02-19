#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_free__deallocate_host.
 * Kernel to deallocate rfmstruct arrays.
 */
static void rfm_precompute_free__deallocate_host(rfm_struct *restrict rfmstruct) {
  // Temporary parameters
} // END FUNCTION rfm_precompute_free__deallocate_host

/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__Cartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                         rfm_struct *restrict rfmstruct) {
  rfm_precompute_free__deallocate_host(rfmstruct);
} // END FUNCTION rfm_precompute_free__rfm__Cartesian

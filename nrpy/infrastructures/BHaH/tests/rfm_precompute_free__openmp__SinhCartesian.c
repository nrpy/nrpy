#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_free__deallocate_host.
 * Kernel to deallocate rfmstruct arrays.
 */
static void rfm_precompute_free__deallocate_host(rfm_struct *restrict rfmstruct) {
  // Temporary parameters
  free(rfmstruct->f0_of_xx0);
  free(rfmstruct->f0_of_xx0__D0);
  free(rfmstruct->f0_of_xx0__DD00);
  free(rfmstruct->f1_of_xx1);
  free(rfmstruct->f1_of_xx1__D1);
  free(rfmstruct->f1_of_xx1__DD11);
  free(rfmstruct->f3_of_xx2);
  free(rfmstruct->f3_of_xx2__D2);
  free(rfmstruct->f3_of_xx2__DD22);
} // END FUNCTION rfm_precompute_free__deallocate_host

/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__SinhCartesian(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                             rfm_struct *restrict rfmstruct) {
  rfm_precompute_free__deallocate_host(rfmstruct);
} // END FUNCTION rfm_precompute_free__rfm__SinhCartesian

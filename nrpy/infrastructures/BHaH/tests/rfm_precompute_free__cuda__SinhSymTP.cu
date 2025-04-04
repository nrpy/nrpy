#include "../BHaH_defines.h"
/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__SinhSymTP(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                         rfm_struct *restrict rfmstruct) {
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__D0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DD00);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DDD000);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1__D1);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f1_of_xx1__DD11);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f2_of_xx0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f2_of_xx0__D0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f2_of_xx0__DD00);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f4_of_xx1);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f4_of_xx1__D1);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f4_of_xx1__DD11);
} // END FUNCTION rfm_precompute_free__rfm__SinhSymTP

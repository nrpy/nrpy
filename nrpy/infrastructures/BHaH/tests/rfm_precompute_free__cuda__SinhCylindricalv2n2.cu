#include "../BHaH_defines.h"
/**
 * rfm_precompute_free: reference metric precomputed lookup arrays: free
 */
void rfm_precompute_free__rfm__SinhCylindricalv2n2(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                   rfm_struct *restrict rfmstruct) {
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__D0);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DD00);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DDD000);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f3_of_xx2);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f3_of_xx2__D2);
  BHAH_FREE_DEVICE__PtrMember(rfmstruct, f3_of_xx2__DD22);
} // END FUNCTION rfm_precompute_free__rfm__SinhCylindricalv2n2

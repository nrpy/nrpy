#include "../BHaH_defines.h"
/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__SinhCylindricalv2n2(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                     rfm_struct *restrict rfmstruct) {
  BHAH_MALLOC__PtrMember(rfmstruct, f0_of_xx0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC__PtrMember(rfmstruct, f0_of_xx0__D0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC__PtrMember(rfmstruct, f0_of_xx0__DD00, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC__PtrMember(rfmstruct, f0_of_xx0__DDD000, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC__PtrMember(rfmstruct, f3_of_xx2, sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
  BHAH_MALLOC__PtrMember(rfmstruct, f3_of_xx2__D2, sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
  BHAH_MALLOC__PtrMember(rfmstruct, f3_of_xx2__DD22, sizeof(REAL) * params->Nxx_plus_2NGHOSTS2);
} // END FUNCTION rfm_precompute_malloc__rfm__SinhCylindricalv2n2

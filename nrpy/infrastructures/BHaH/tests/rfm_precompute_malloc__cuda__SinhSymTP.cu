#include "BHaH_defines.h"

/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__SinhSymTP(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                           rfm_struct *restrict rfmstruct) {
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f0_of_xx0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f0_of_xx0__D0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DD00, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f0_of_xx0__DDD000, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1__D1, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f1_of_xx1__DD11, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f2_of_xx0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f2_of_xx0__D0, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f2_of_xx0__DD00, sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f4_of_xx1, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f4_of_xx1__D1, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  BHAH_MALLOC_DEVICE__PtrMember(rfmstruct, f4_of_xx1__DD11, sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
} // END FUNCTION rfm_precompute_malloc__rfm__SinhSymTP

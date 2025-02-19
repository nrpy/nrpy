#include "../BHaH_defines.h"
/**
 * Kernel: rfm_precompute_malloc__allocate_host.
 * Kernel to allocate rfmstruct arrays.
 */
static void rfm_precompute_malloc__allocate_host(const params_struct *restrict params, rfm_struct *restrict rfmstruct) {
  // Temporary parameters
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
  rfmstruct->f0_of_xx0 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__D0 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DDD000 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS0);
  rfmstruct->f1_of_xx1 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__D1 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__DD11 = (REAL *)malloc(sizeof(REAL) * params->Nxx_plus_2NGHOSTS1);
} // END FUNCTION rfm_precompute_malloc__allocate_host

/**
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__HoleySinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                    rfm_struct *restrict rfmstruct) {
  rfm_precompute_malloc__allocate_host(params, rfmstruct);
} // END FUNCTION rfm_precompute_malloc__rfm__HoleySinhSpherical

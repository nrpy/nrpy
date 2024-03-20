#include "../BHaH_defines.h"
/*
 * rfm_precompute_malloc: reference metric precomputed lookup arrays: malloc
 */
void rfm_precompute_malloc__rfm__SinhSpherical(const commondata_struct *restrict commondata, const params_struct *restrict params,
                                               rfm_struct *restrict rfmstruct) {
#include "../set_CodeParameters.h"
  rfmstruct->f0_of_xx0 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__D0 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DD00 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f0_of_xx0__DDD000 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS0);
  rfmstruct->f1_of_xx1 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__D1 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
  rfmstruct->f1_of_xx1__DD11 = (REAL *)malloc(sizeof(REAL) * Nxx_plus_2NGHOSTS1);
}

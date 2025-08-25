#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate the gamma function using GSL.
 */
double complex SEOBNRv5_aligned_spin_gamma_wrapper(const REAL z_real, const REAL z_imag) {

  gsl_sf_result lnr, arg;
  int status = gsl_sf_lngamma_complex_e(z_real, z_imag, &lnr, &arg);
  int status_desired[1] = {GSL_SUCCESS};
  char lngamma_name[] = "gsl_sf_lngamma_complex_e";
  handle_gsl_return_status(status, status_desired, 1, lngamma_name);
  return cexp(lnr.val + I * arg.val);
} // END FUNCTION SEOBNRv5_aligned_spin_gamma_wrapper

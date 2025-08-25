#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate the SEOBNRv5 dissipative initial conditions.
 */
int SEOBNRv5_aligned_spin_initial_conditions_dissipative(commondata_struct *restrict commondata) {

  int status;
  int iter = 0;
  const int max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  REAL prstar;
  REAL x_lo = -3e-2;
  REAL x_hi = 0.0;
  REAL xtol = 1e-12;
  REAL rtol = 1e-10;
  gsl_function F;
  F.function = &SEOBNRv5_aligned_spin_radial_momentum_conditions;
  F.params = commondata;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    int fsolver_status[1] = {GSL_SUCCESS};
    char fsolver_name[] = "gsl_root_fsolver_iterate";
    handle_gsl_return_status(status, fsolver_status, 1, fsolver_name);
    prstar = gsl_root_fsolver_root(s);

    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, xtol, rtol);
    int test_interval_status[2] = {GSL_SUCCESS, GSL_CONTINUE};
    char root_test_name[] = "gsl_root_test_interval";
    handle_gsl_return_status(status, test_interval_status, 2, root_test_name);

  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);
  commondata->prstar = prstar;
  return status;
} // END FUNCTION SEOBNRv5_aligned_spin_initial_conditions_dissipative

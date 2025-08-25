#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_vector.h>

/**
 * Evaluate the SEOBNRv5 conservative initial conditions.
 */
int SEOBNRv5_aligned_spin_initial_conditions_conservative(commondata_struct *restrict commondata) {

  const size_t n = 2;
  gsl_multiroot_function f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit, n, commondata};
  REAL omega = commondata->initial_omega;
  REAL pphi = pow(omega, -1. / 3.);
  REAL r = pphi * pphi;
  const REAL x_guess[2] = {r, pphi};
  REAL *restrict x_result = malloc(2 * sizeof(REAL));
  const gsl_multiroot_fsolver_type *restrict T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *restrict s = gsl_multiroot_fsolver_alloc(T, n);
  gsl_vector *restrict x = gsl_vector_alloc(n);
  size_t i, iter = 0;
  int status;
  const int maxiter = 100;
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, i, x_guess[i]);
  }
  gsl_multiroot_fsolver_set(s, &f, x);
  do {
    iter++;
    status = gsl_multiroot_fsolver_iterate(s);
    int f_solver_status[1] = {GSL_SUCCESS};
    char fsolver_name[] = "gsl_multiroot_fsolver_iterate";
    handle_gsl_return_status(status, f_solver_status, 1, fsolver_name);
    status = gsl_multiroot_test_residual(s->f, 6e-12);
    int test_residual_status[2] = {GSL_SUCCESS, GSL_CONTINUE};
    char residual_name[] = "gsl_multiroot_test_residual";
    handle_gsl_return_status(status, test_residual_status, 2, residual_name);
  } while (status == GSL_CONTINUE && iter < maxiter);
  for (i = 0; i < n; i++) {
    x_result[i] = gsl_vector_get(s->x, i);
  }
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(x);
  commondata->r = x_result[0];
  commondata->pphi = x_result[1];
  free(x_result);
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_initial_conditions_conservative

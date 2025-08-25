#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Multidimensional root finder using GSL.
 */
void SEOBNRv5_aligned_multidimensional_root_wrapper(gsl_multiroot_function_fdf f, const REAL *restrict x_guess, const size_t n,
                                                    REAL *restrict x_result) {

  size_t i, iter = 0;
  const int maxiter = 100;
  int status;
  const gsl_multiroot_fdfsolver_type *restrict T = gsl_multiroot_fdfsolver_hybridsj;
  gsl_multiroot_fdfsolver *restrict s = gsl_multiroot_fdfsolver_alloc(T, n);
  gsl_vector *restrict x = gsl_vector_alloc(n);
  for (i = 0; i < n; i++) {
    gsl_vector_set(x, i, x_guess[i]);
  }
  gsl_multiroot_fdfsolver_set(s, &f, x);
  do {
    iter++;
    status = gsl_multiroot_fdfsolver_iterate(s);
    int fdf_solver_status[1] = {GSL_SUCCESS};
    char fdfsolver_name[] = "gsl_multiroot_fdfsolver_iterate";
    handle_gsl_return_status(status, fdf_solver_status, 1, fdfsolver_name);
    status = gsl_multiroot_test_residual(s->f, 6e-12);
    int test_residual_status[2] = {GSL_SUCCESS, GSL_CONTINUE};
    char residual_name[] = "gsl_multiroot_test_residual";
    handle_gsl_return_status(status, test_residual_status, 2, residual_name);
  } while (status == GSL_CONTINUE && iter < maxiter);
  for (i = 0; i < n; i++) {
    x_result[i] = gsl_vector_get(s->x, i);
  }
  gsl_multiroot_fdfsolver_free(s);
  gsl_vector_free(x);
} // END FUNCTION SEOBNRv5_aligned_multidimensional_root_wrapper

"""
Set up C function library for SEOBNR related GSL wrappers.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_multidimensional_root_wrapper() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for performing multidimensional root-finding using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Wrapper function for performing multidimensional root-finding using GSL when higher derivatives of the root function are available.

@param f - The GSL structure containing the root function and derivatives.
@param x_guess - The initial guess for the root.
@param n - The dimensionality of the root-finding problem.
@param x_result - The result of the root-finding.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_multidimensional_root_wrapper"
    params = "gsl_multiroot_function_fdf f, const REAL *restrict x_guess, const size_t n, REAL *restrict x_result"
    body = """
size_t i , iter = 0;
const int maxiter = 100;
int status;
const gsl_multiroot_fdfsolver_type *restrict T = gsl_multiroot_fdfsolver_hybridsj;
if (T == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_multidimensional_root_wrapper(), could not assign gsl_multiroot_fsolver_hybrids to T\\n");
  exit(1);
}
gsl_multiroot_fdfsolver *restrict s = gsl_multiroot_fdfsolver_alloc(T , n);
if (s == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_multidimensional_root_wrapper(), gsl_multiroot_fdfsolver_alloc failed to initialize\\n");
  exit(1);
}
gsl_vector *restrict x = gsl_vector_alloc(n);
if (x == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_multidimensional_root_wrapper(), gsl_vector_alloc failed to initialize\\n");
  exit(1);
}
for (i = 0; i < n; i++){
gsl_vector_set(x , i , x_guess[i]);
}
gsl_multiroot_fdfsolver_set(s , &f , x);
do {
  iter++;
  status = gsl_multiroot_fdfsolver_iterate (s);
  int fdf_solver_status[1] = {GSL_SUCCESS};
  char fdfsolver_name[] = "gsl_multiroot_fdfsolver_iterate";
  handle_gsl_return_status(status,fdf_solver_status,1,fdfsolver_name);
  status = gsl_multiroot_test_residual (s->f, 6e-12);
  int test_residual_status[2] = {GSL_SUCCESS,GSL_CONTINUE};
  char residual_name[] = "gsl_multiroot_test_residual";
  handle_gsl_return_status(status,test_residual_status,2,residual_name);
}
while(status == GSL_CONTINUE && iter < maxiter);
for (i = 0; i < n; i++){
x_result[i] = gsl_vector_get(s->x , i);
}
gsl_multiroot_fdfsolver_free (s);
gsl_vector_free (x);
"""
    cfc.register_CFunction(
        subdirectory="initial_conditions",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

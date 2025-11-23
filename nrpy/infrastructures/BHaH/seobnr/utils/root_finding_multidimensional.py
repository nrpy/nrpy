"""
Set up C function library for multidiminsional root finding using GSL, when the Jacobian of the root function is not specified.

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


def register_CFunction_root_finding_multidimensional() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating multidimensional root finding using GSL, when the Jacobian of the root function is not specified.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
"""
    desc = """
Evaluates the multidimensional root finding using GSL, when the Jacobian of the root function is not specified.
@params n - The dimensionality of the root-finding problem.
@params x_guess - The initial guess for the root.
@params f - The GSL structure containing the root function and derivatives.
@params x_result - The result of the root-finding.
"""
    cfunc_type = "void"
    name = "root_finding_multidimensional"
    params = "const size_t n, const REAL *restrict x_guess, gsl_multiroot_function *restrict f, REAL *restrict x_result"
    body = """
if (x_result == NULL){
  fprintf(stderr, "In root_finding_multidimensional(), x_result passed uninitialized\\n");
  exit(1);
}
const gsl_multiroot_fsolver_type *restrict T = gsl_multiroot_fsolver_hybrids;
if (T == NULL){
  fprintf(stderr,"Error: in root_finding_multidimensional(), could not assign gsl_multiroot_fsolver_hybrids to T\\n");
  exit(1);
}
gsl_multiroot_fsolver *restrict s = gsl_multiroot_fsolver_alloc(T , n);
if (s == NULL){
  fprintf(stderr,"Error: in root_finding_multidimensional(), gsl_multiroot_fsolver_alloc failed to initialize\\n");
  exit(1);
}
gsl_vector *restrict x = gsl_vector_alloc(n);
if (x == NULL){
  fprintf(stderr,"Error: in root_finding_multidimensional(), gsl_vector_alloc failed to initialize\\n");
  exit(1);
}
size_t i , iter = 0;
int status;
const int maxiter = 100;
for (i = 0; i < n; i++){
gsl_vector_set(x , i , x_guess[i]);
}
gsl_multiroot_fsolver_set(s , f , x);
do {
  iter++;
  status = gsl_multiroot_fsolver_iterate (s);
  int f_solver_status[1] = {GSL_SUCCESS};
  char fsolver_name[] = "gsl_multiroot_fsolver_iterate";
  handle_gsl_return_status(status,f_solver_status,1,fsolver_name);
  status = gsl_multiroot_test_residual (s->f, 6e-12);
  int test_residual_status[2] = {GSL_SUCCESS,GSL_CONTINUE};
  char residual_name[] = "gsl_multiroot_test_residual";
  handle_gsl_return_status(status,test_residual_status,2,residual_name);
}
while(status == GSL_CONTINUE && iter < maxiter);
for (i = 0; i < n; i++){
x_result[i] = gsl_vector_get(s->x , i);
}
gsl_multiroot_fsolver_free (s);
gsl_vector_free (x);
"""
    cfc.register_CFunction(
        subdirectory="utils",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

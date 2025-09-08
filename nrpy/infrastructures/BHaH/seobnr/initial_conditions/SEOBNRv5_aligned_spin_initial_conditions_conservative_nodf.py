"""
Set up C function library for SEOBNR initial conditions.

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


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative_nodf() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5's conservative initial conditions using GSL.

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
Evaluates the conservative initial conditions for the SEOBNRv5 ODE integration.
The conservative initial conditions are given by a circular orbit (p_{r_*} = 0) such that the 
the time derivative of the tortoise momentum is zero and the time derivative of the orbital phase equals the input orbital frequency.
@params commondata - The Common data structure containing the model parameters.
@returns - GSL_SUCCESS (0) upon success.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_initial_conditions_conservative"
    params = "commondata_struct *restrict commondata"
    body = """
const size_t n = 2;
gsl_multiroot_function f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit,
                                n,
                                commondata};
REAL omega = commondata->initial_omega;
REAL pphi = pow(omega,-1./3.);
REAL r = pphi*pphi;
const REAL x_guess[2] = {r,pphi};
REAL *restrict x_result = malloc(2*sizeof(REAL));
const gsl_multiroot_fsolver_type *restrict T = gsl_multiroot_fsolver_hybrids;
gsl_multiroot_fsolver *restrict s = gsl_multiroot_fsolver_alloc(T , n);
gsl_vector *restrict x = gsl_vector_alloc(n);
size_t i , iter = 0;
int status;
const int maxiter = 100;
for (i = 0; i < n; i++){
gsl_vector_set(x , i , x_guess[i]);
}
gsl_multiroot_fsolver_set(s , &f , x);
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
commondata->r = x_result[0];
commondata->pphi = x_result[1];
free(x_result);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="initial_conditions",
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

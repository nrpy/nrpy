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


def register_CFunction_handle_gsl_return_status() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for handling error statuses returned by GSL calls.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Handle GSL return status."""
    cfunc_type = "void"
    name = "handle_gsl_return_status"
    params = "int status, int status_desired[], int num_desired, const char *restrict function_name"
    body = """
int count = 0;
for (int i = 0; i < num_desired; i++){
  if (status == status_desired[i]){
    count++;
  }
}
if (count == 0){
  printf ("In function %s, gsl returned error: %s\\nAborted", function_name, gsl_strerror(status));
  exit(EXIT_FAILURE);
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper_complex_out() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the complex gamma function using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate the gamma function using GSL."""
    cfunc_type = "double complex"
    name = "SEOBNRv5_aligned_spin_gamma_wrapper"
    params = "const double complex z"
    body = """
const REAL z_real = (REAL) creal(z);
const REAL z_imag = (REAL) cimag(z);
gsl_sf_result lnr, arg;
int status = gsl_sf_lngamma_complex_e(z_real, z_imag, &lnr, &arg);
int status_desired[1] = {GSL_SUCCESS};
char lngamma_name[] = "gsl_sf_lngamma_complex_e";
handle_gsl_return_status(status,status_desired,1,lngamma_name);
double complex complex_gamma = cexp(lnr.val + I*arg.val);
return complex_gamma;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_gamma_wrapper() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the complex gamma function using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate the gamma function using GSL."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_gamma_wrapper"
    params = "const REAL z_real, const REAL z_imag, REAL *restrict gamma_z"
    body = """
gsl_sf_result lnr, arg;
int status = gsl_sf_lngamma_complex_e(z_real, z_imag, &lnr, &arg);
int status_desired[1] = {GSL_SUCCESS};
char lngamma_name[] = "gsl_sf_lngamma_complex_e";
handle_gsl_return_status(status,status_desired,1,lngamma_name);
const REAL gamma_amp = exp(lnr.val);
const REAL gamma_phase = arg.val;
gamma_z[0] =  gamma_amp*cos(gamma_phase);
gamma_z[1] =  gamma_amp*sin(gamma_phase);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


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
    desc = """Multidimensional root finder using GSL."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_multidimensional_root_wrapper"
    params = "gsl_multiroot_function_fdf f, REAL *restrict x_guess, const size_t n, REAL *restrict x_result"
    body = """
size_t i , iter = 0;
const int maxiter = 100;
int status;
const gsl_multiroot_fdfsolver_type *restrict T = gsl_multiroot_fdfsolver_hybridsj;
gsl_multiroot_fdfsolver *restrict s = gsl_multiroot_fdfsolver_alloc(T , n);
gsl_vector *restrict x = gsl_vector_alloc(n);
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
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

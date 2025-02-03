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
    cfunc_type = "double complex"
    name = "SEOBNRv5_aligned_spin_gamma_wrapper"
    params = "const REAL z_real, const REAL z_imag"
    body = """
gsl_sf_result lnr, arg;
int status = gsl_sf_lngamma_complex_e(z_real, z_imag, &lnr, &arg);
int status_desired[1] = {GSL_SUCCESS};
char lngamma_name[] = "gsl_sf_lngamma_complex_e";
handle_gsl_return_status(status,status_desired,1,lngamma_name);
return cexp(lnr.val + I*arg.val);
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


def register_CFunction_SEOBNRv5_aligned_spin_process_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for windowing and zero-padding the (2,2) time-domain mode of SEOBNRv5 for Fourier tranformation.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
    Window and zero-pad the SEOBNRv5 time-domain 22 mode for Fourier transformation.
    
    A Tukey window is applied to the time-domain 22 mode using the window function given by
    Bloomfield, P. Fourier Analysis of Time Series: An Introduction. New York: Wiley-Interscience, 2000.
    """
    cfunc_type = "void"
    prefunc = "#include<fftw3.h>"
    name = "SEOBNRv5_aligned_spin_process_waveform"
    params = "commondata_struct *restrict commondata"
    body = r"""
const size_t window_length = commondata->nsteps_IMR;
const size_t window_length_minus_1 = window_length - 1;

// For now, padding length is 2x10 but we can make this tuneable in the future.
const size_t padding_length_one_side = 10;
const size_t padding_length = 2 * padding_length_one_side;

// For now, cosine fraction is fixed but we can make this tuneable in the future.
const REAL cosine_fraction = .1;
const REAL cosine_fraction_over_2 = cosine_fraction * 0.5;

COMPLEX *restrict processed_waveform = (COMPLEX *restrict)malloc((padding_length + commondata->nsteps_IMR) * NUMMODES * sizeof(COMPLEX));
REAL w;
const size_t width = (size_t)floor(window_length_minus_1 * cosine_fraction_over_2);
const size_t padding_left = padding_length_one_side;
const size_t cosine_left = width + 1 + padding_length_one_side;
const size_t cosine_right = window_length_minus_1 - width + padding_length_one_side;
const size_t padding_right = window_length + padding_length_one_side;
size_t total_size = commondata->nsteps_IMR + padding_length;
size_t i;
for (i = 0; i < padding_left; i++) {
processed_waveform[IDX_WF(i, TIME)] = i * commondata->dT;
processed_waveform[IDX_WF(i, STRAIN)] = 0. + I * 0.;
}
for (i = padding_left; i < cosine_left; i++) {
w = 0.5 * (1 + cos(M_PI * (-1 + (i - padding_left) / cosine_fraction_over_2 / window_length_minus_1)));
processed_waveform[IDX_WF(i, TIME)] = (i) * commondata->dT;
processed_waveform[IDX_WF(i, STRAIN)] = commondata->waveform_IMR[IDX_WF(i - padding_left, STRAIN)] * (COMPLEX)w;
}
for (i = cosine_left; i < cosine_right; i++) {
processed_waveform[IDX_WF(i, TIME)] = (i) * commondata->dT;
processed_waveform[IDX_WF(i, STRAIN)] = commondata->waveform_IMR[IDX_WF(i - padding_left, STRAIN)];
}
for (i = cosine_right; i < padding_right; i++) {
w = 0.5 * (1 + cos(M_PI * (-1 / cosine_fraction_over_2 + 1 + (i - padding_left) / cosine_fraction_over_2 / window_length_minus_1)));
processed_waveform[IDX_WF(i, TIME)] = (i) * commondata->dT;
processed_waveform[IDX_WF(i, STRAIN)] = commondata->waveform_IMR[IDX_WF(i - padding_left, STRAIN)] * (COMPLEX)w;
}
for (i = window_length; i < total_size; i++) {
processed_waveform[IDX_WF(i, TIME)] = (i) * commondata->dT;
processed_waveform[IDX_WF(i, STRAIN)] = 0. + I * 0.;
}
commondata->waveform_IMR = (double complex *restrict)realloc(commondata->waveform_IMR,total_size*NUMMODES*sizeof(double complex));
commondata->nsteps_IMR = total_size;
memcpy(commondata->waveform_IMR,processed_waveform,total_size*NUMMODES*sizeof(double complex));
free(processed_waveform);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        prefunc=prefunc,
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
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_multidimensional_root_wrapper"
    params = "gsl_multiroot_function_fdf f, const REAL *restrict x_guess, const size_t n, REAL *restrict x_result"
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

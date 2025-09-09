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
    desc = """
Wrapper function for evaluating the complex gamma function using GSL's lngamma_complex_e.

@param z_real - The real part of the input complex number.
@param z_imag - The imaginary part of the input complex number.
@returns - The gamma function evaluated at z.
"""
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
        subdirectory="inspiral_waveform",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

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
    desc = """
Error handler for calls to GSL-dependent routines.

@param status - The return value of the GSL-dependent call.
@param status_desired - The desired return statuses of the GSL call.
@param num_desired - The number of desired return statuses.
@param function_name - The name of the GSL-dependent function that was called.
"""
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
        subdirectory="utils",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

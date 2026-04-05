"""
Set up C function library for SEOBNR and BOB attachment routines.

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


def register_CFunction_SEOBNRv5_aligned_spin_unwrap() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction that performs numpy.unwrap on a given array.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Unwraps an array of angles with a period of 2pi.

@param angles_in - Array of angles to unwrap.
@param angles_out - Array to store the unwrapped angles.
@param nsteps_arr - length of the angles_in array.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_unwrap"
    params = "REAL *restrict angles_in , REAL *restrict angles_out, size_t nsteps_arr"
    body = """
angles_out[0] = angles_in[0];
REAL diff;
for (size_t i = 1; i < nsteps_arr; i++){
  diff = angles_in[i] - angles_in[i-1];
  diff = fabs(diff) > M_PI ? (diff < - M_PI ? diff + 2 * M_PI : diff - 2 * M_PI) : diff;
  angles_out[i] = angles_out[i - 1] + diff;
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

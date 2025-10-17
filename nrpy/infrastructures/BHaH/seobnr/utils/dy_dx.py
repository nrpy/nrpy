"""
Register CFunction for computing the derivative of an array by finite difference.

Authors: Siddharth Mahesh
sm0193 at mix dot wvu dot edu
Zachariah B. Etienne
zachetie at gmail dot com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_dy_dx() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for computing the derivative of an array by finite difference.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
    Compute the derivative of an array by finite difference.

    @param y - Array to differentiate.
    @param x - Array of independent variables.
    @param dy_dx - Array to store the derivative.
    @param nsteps - Size of the array.
    """
    cfunc_type = "void"
    name = "dy_dx"
    params = (
        "REAL *restrict y, REAL *restrict x, REAL *restrict dy_dx, const size_t nsteps"
    )
    body = """
REAL coeffs[8] = {0.};
int indices[8] = {0,0,0,0,0,0,0,0};
int offset = -4;
const REAL h = x[1] - x[0];
for (size_t i = 0; i < nsteps; i++){
    dy_dx[i] = 0.;
    offset = 4 - (int)MIN(i,8);
    finite_difference_stencil(offset,coeffs,indices);
    for (int j = 0; j < 8; j++){
        dy_dx[i] += coeffs[j]*y[i+indices[j]]*h;
    }
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

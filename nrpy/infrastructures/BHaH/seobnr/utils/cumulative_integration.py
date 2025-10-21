"""
Register CFunction for computing the cumulative integral of an array by interpolated integration stencils.

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


def register_CFunction_cumulative_integration() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for computing the cumulative integral of an array by interpolated integration stencils.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
    Compute the cumulative integral of an array by interpolated integration stencils.
    This function calculates the integrals
          /x_i
    C_i = | f(x) dx
          /x_0
    for i = 1,2,...,nsteps. The cumulative integration is performed by adding up the individual integrals
          j = i - 1
          _
          \
    C_i =    I_j
          /_
          j = 1
    where I_j is the integral from x_{j-1} to x_j calculated by the integration stencil function.

    @param y - Array to integrate.
    @param x - Array of independent variables.
    @param C - Array to store the cumulative integrals.
    @param nsteps - Size of the array.
    """
    cfunc_type = "void"
    name = "cumulative_integration"
    params = "REAL *restrict y, REAL *restrict x, REAL *restrict C, const size_t nsteps"
    body = """
if (y == NULL){
    fprintf(stderr,"Error: in cumulative_integration, y is passed uninitialised\\n");
    exit(1);
}
if (x == NULL){
    fprintf(stderr,"Error: in cumulative_integration, x is passed uninitialised\\n");
    exit(1);
}
if (C == NULL){
    fprintf(stderr,"Error: in cumulative_integration, C is passed uninitialised\\n");
    exit(1);
}
REAL coeffs[8] = {0.};
int indices[8] = {0,0,0,0,0,0,0,0};
int offset = 3;
const REAL h = x[1] - x[0]; // assumes equally sampled, true for PA integration
REAL running_integral = 0.;
//starting time is zero by default
C[0] = running_integral;
//begin forward interpolated integrations
for (size_t i = 1; i < 4; i++){
  integration_stencil(offset,coeffs,indices);
  for (int j = 0; j < 8; j++){
    running_integral += coeffs[j]*y[i+indices[j] - 1];
  }
  C[i] = running_integral * h;
  offset--;
}
//begin central interpolated integrations
offset = 0;
for (size_t i = 4; i < nsteps - 3; i++){
  integration_stencil(offset,coeffs,indices);
  for (int j = 0; j < 8; j++){
    running_integral += coeffs[j]*y[i+indices[j] - 1];
  }
  C[i] = running_integral * h;
}
//begin backward interpolated integrations
for (size_t i = nsteps - 3; i < nsteps; i++){
  offset--;
  integration_stencil(offset,coeffs,indices);
  for (int j = 0; j < 8; j++){
    running_integral += coeffs[j]*y[i+indices[j] - 1];
  }
  C[i] = running_integral * h;
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

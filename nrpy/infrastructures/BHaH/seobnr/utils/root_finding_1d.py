"""
Set up C function library for 1-dimensional root finding using GSL.

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


def register_CFunction_root_finding_1d() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for implementing 1-dimensional root finding using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Implements 1-dimensional root finding using GSL.

@param x_low - Lower bound of the root finding interval.
@param x_high - Upper bound of the root finding interval.
@param F - GSL function to find the root of.
@returns x - The root of the function F.
"""
    cfunc_type = "REAL"
    name = "root_finding_1d"
    params = "const REAL x_low, const REAL x_high, gsl_function *restrict F"
    body = """
int status;
int iter = 0;
const int max_iter = 100;
const gsl_root_fsolver_type *T;
gsl_root_fsolver *s;
REAL x;
REAL x_lo = x_low;
REAL x_hi = x_high;
REAL xtol = 1e-12;
REAL rtol = 1e-10;
T = gsl_root_fsolver_brent;
s = gsl_root_fsolver_alloc(T);
status = gsl_root_fsolver_set(s, F, x_lo, x_hi);
int fsolver_set_status[1] = {GSL_SUCCESS};
char fsolver_set_name[] = "gsl_root_fsolver_set";
handle_gsl_return_status(status,fsolver_set_status,1,fsolver_set_name);

do {
  iter++;
  status = gsl_root_fsolver_iterate(s);
  int fsolver_iterate_status[1] = {GSL_SUCCESS};
  char fsolver_name[] = "gsl_root_fsolver_iterate";
  handle_gsl_return_status(status,fsolver_status,1,fsolver_name);
  x = gsl_root_fsolver_root(s);

  x_lo = gsl_root_fsolver_x_lower(s);
  x_hi = gsl_root_fsolver_x_upper(s);
  status = gsl_root_test_interval(x_lo, x_hi, xtol, rtol);
  int test_interval_status[2] = {GSL_SUCCESS,GSL_CONTINUE};
  char root_test_name[] = "gsl_root_test_interval";
  handle_gsl_return_status(status,test_interval_status,2,root_test_name);


}
while (status == GSL_CONTINUE && iter < max_iter);

gsl_root_fsolver_free (s);
return x;
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

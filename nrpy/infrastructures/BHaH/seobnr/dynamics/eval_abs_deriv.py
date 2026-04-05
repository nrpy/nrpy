"""
Register CFunction for evaluating the absolute value of the derivative of a spline at a given point.

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


def register_CFunction_eval_abs_deriv() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the absolute value of the derivative of a spline at a given point.
    Needed by SEOBNRv5_aligned_spin_iterative_refinement.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the absolute value of the derivative of a spline at a given point.

@param t - The point at which to evaluate the derivative.
@param params - The spline data.
@returns - The absolute value of the derivative of the spline at the given point.
"""
    cfunc_type = "double"
    name = "eval_abs_deriv"
    params = "double t, void *params"
    body = """
spline_data *sdata = (spline_data *)params;
return fabs(gsl_spline_eval_deriv(sdata->spline, t, sdata->acc));
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

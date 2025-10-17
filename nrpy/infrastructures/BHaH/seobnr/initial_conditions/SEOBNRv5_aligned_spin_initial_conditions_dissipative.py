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


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_dissipative() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5's dissipative initial conditions using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the dissipative initial conditions for the SEOBNRv5 ODE integration.

@params commondata - The Common data structure containing the model parameters.
@returns - GSL_SUCCESS (0) upon success.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_initial_conditions_dissipative"
    params = "commondata_struct *restrict commondata"
    body = """
REAL x_lo = -3e-2;
REAL x_hi = 0.0;
gsl_function F;
F.function = &SEOBNRv5_aligned_spin_radial_momentum_condition;
F.params = commondata;
commondata->prstar = root_finding_1d(x_lo, x_hi, &F);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="initial_conditions",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

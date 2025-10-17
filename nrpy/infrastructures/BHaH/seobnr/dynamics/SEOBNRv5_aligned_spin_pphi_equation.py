"""
Register CFunction for evaluating the angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

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


def register_CFunction_SEOBNRv5_pphi_equation() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluate the angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

@params x - The angular momentum.
@params params - The Common data structure containing the model parameters.
@returns - The angular momentum equation.
"""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_pphi_equation"
    params = "REAL x, void *restrict params"
    body = """
const REAL r = ((commondata_struct *restrict) params)->r;
const REAL prstar = ((commondata_struct *restrict) params)->prstar;
const REAL dprstar_dr = ((commondata_struct *restrict) params)->dprstar_dr;
REAL pphi = x;
REAL dynamics[4] = {r,0.,prstar,pphi};
REAL f[4] = {0.,0.,0.,0.};
SEOBNRv5_aligned_spin_right_hand_sides(0,dynamics,f,params);
REAL prstar_dot = f[2];
REAL rdot = f[0];
REAL pphi_equation = rdot*dprstar_dr - prstar_dot;
return pphi_equation;
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

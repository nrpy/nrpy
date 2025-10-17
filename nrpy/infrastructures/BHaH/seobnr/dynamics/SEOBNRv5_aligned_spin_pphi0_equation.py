"""
Register CFunction for evaluating the adiabatic angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_pphi0_equation() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the adiabatic angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluate the adiabatic angular momentum equation for the SEOBNRv5 post-adiabatic equations of motion.

@params x - The angular momentum.
@params params - The Common data structure containing the model parameters.
@returns - The angular momentum equation.
"""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_pphi0_equation"
    params = "REAL x, void *restrict params"
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    dHreal_dr_circ = Hq.dHreal_dr_circ
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
const REAL r = ((commondata_struct *restrict) params)->r;
REAL pphi = x;
"""
    body += ccg.c_codegen(
        dHreal_dr_circ,
        "REAL pphi0_equation",
        verbose=False,
        include_braces=False,
    )
    body += """
return pphi0_equation;
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

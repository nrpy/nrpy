"""
Register CFunction for evaluating the right hand sides for the SEOBNRv5 equations of motion.

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


def register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the right hand sides for the SEOBNRv5 equations of motion.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluate SEOBNRv5 Hamiltonian, flux, and needed derivatives to compute binary dynamics.

@param t - The current time.
@param y - Array of dynamical variables.
@param f - Array to store the right-hand sides.
@param params - Common data structure containing the model parameters.
@return - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_right_hand_sides"
    params = "REAL t, const REAL *restrict y, REAL *restrict f, void *restrict params"
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
const REAL r = y[0];
const REAL prstar = y[2];
const REAL pphi = y[3];
"""
    body += ccg.c_codegen(
        [
            Hq.Hreal,
            Hq.xi,
            Hq.dHreal_dr,
            Hq.dHreal_dprstar,
            Hq.dHreal_dpphi,
            Hq.dHreal_dpphi_circ,
        ],
        [
            "const REAL Hreal",
            "const REAL xi",
            "const REAL dHreal_dr",
            "const REAL dHreal_dprstar",
            "const REAL dHreal_dpphi",
            "const REAL dHreal_dpphi_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
REAL flux[2];
SEOBNRv5_aligned_spin_flux(y,Hreal,dHreal_dpphi,dHreal_dpphi_circ,flux,params);
f[0] = xi * dHreal_dprstar;
f[1] = dHreal_dpphi;
f[2] = -xi * dHreal_dr + flux[0];
f[3] = flux[1];
return GSL_SUCCESS;
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

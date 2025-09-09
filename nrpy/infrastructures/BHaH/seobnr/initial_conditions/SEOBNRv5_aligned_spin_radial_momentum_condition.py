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

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_radial_momentum_condition() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 initial radial tortoise momentum.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the SEOBNRv5 adiabatic radial momentum condition.

@params x - The radial tortoise momentum.
@params params - The Common data structure containing the model parameters.
@returns - The radial momentum condition.
"""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_radial_momentum_condition"
    params = "REAL x, void *restrict params"
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
const REAL r = ((commondata_struct *restrict) params)->r;
const REAL pphi = ((commondata_struct *restrict) params)->pphi;
REAL prstar = x;
"""
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body += ccg.c_codegen(
        [
            H.dHreal_dr_dr,
            H.dHreal_dr_dpphi,
            H.xi,
            H.dHreal_dprstar,
            H.Hreal,
            H.dHreal_dpphi,
            H.dHreal_dpphi_circ,
        ],
        [
            "const REAL dHreal_dr_dr",
            "const REAL dHreal_dr_dpphi",
            "const REAL xi",
            "const REAL dHreal_dprstar",
            "const REAL Hreal",
            "const REAL Omega",
            "const REAL Omega_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
const REAL dLdr = - dHreal_dr_dr / dHreal_dr_dpphi;
const REAL rdot_dyn = xi * dHreal_dprstar;
REAL flux[2];
const REAL y[4] = {r , 0. , prstar , pphi};
SEOBNRv5_aligned_spin_flux(y,Hreal,Omega,Omega_circ,flux,params);
const REAL dLdt = flux[1];
const REAL rdot_rad = dLdt / dLdr;
const REAL prstar_condition = rdot_dyn - rdot_rad;
return prstar_condition;
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

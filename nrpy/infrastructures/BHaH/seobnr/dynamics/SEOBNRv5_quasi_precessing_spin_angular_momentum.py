"""
Register CFunction for evaluating the angular momentum from the SEOBNRv5 spin dynamics.

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
import nrpy.equations.seobnr.SEOBNRv5_spin_evolution_equations as SEOBNRv5_spin_eqns
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_quasi_precessing_spin_angular_momentum() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the angular momentum from the SEOBNRv5 spin dynamics.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluate SEOBNRv5 angular momentum from the spin dynamics.

@param z - Array of spin variables.
@param L - Array to store the angular momentum.
@param params - Common data struct containing the model parameters.
@return - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_quasi_precessing_spin_angular_momentum"
    params = "const REAL *restrict z, REAL *restrict L, void *restrict params"
    spin_eqns = SEOBNRv5_spin_eqns.SEOBNRv5_spin_evolution_equations()
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL ln_x = z[LN_X];
const REAL ln_y = z[LN_Y];
const REAL ln_z = z[LN_Z];
const REAL chi1_x = z[CHI1_X];
const REAL chi1_y = z[CHI1_Y];
const REAL chi1_z = z[CHI1_Z];
const REAL chi2_x = z[CHI2_X];
const REAL chi2_y = z[CHI2_Y];
const REAL chi2_z = z[CHI2_Z];
const REAL omega = z[OMEGA_PN];
"""
    body += ccg.c_codegen(
        [
            spin_eqns.L_x,
            spin_eqns.L_y,
            spin_eqns.L_z,
        ],
        [
            "const REAL L_x",
            "const REAL L_y",
            "const REAL L_z",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
L[0] = L_x;
L[1] = L_y;
L[2] = L_z;
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

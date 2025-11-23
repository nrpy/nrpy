"""
Register CFunction for evaluating the right hand sides for the SEOBNRv5 spin evolution equations.

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


def register_CFunction_SEOBNRv5_quasi_precessing_spin_equations() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the right hand sides for the SEOBNRv5 spin evolution equations.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluate SEOBNRv5 Hamiltonian, flux, and needed derivatives to compute binary dynamics.

@param t - The current time.
@param z - Array of spin variables.
@param f - Array to store the right-hand sides.
@param params - Common data structure containing the model parameters.
@return - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "SEOBNRv5_quasi_precessing_spin_equations"
    params = "REAL t, const REAL *restrict z, REAL *restrict f, void *restrict params"
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
            spin_eqns.ln_dot_x,
            spin_eqns.ln_dot_y,
            spin_eqns.ln_dot_z,
            spin_eqns.chi1_dot_x,
            spin_eqns.chi1_dot_y,
            spin_eqns.chi1_dot_z,
            spin_eqns.chi2_dot_x,
            spin_eqns.chi2_dot_y,
            spin_eqns.chi2_dot_z,
            spin_eqns.omega_dot,
        ],
        [
            "const REAL ln_dot_x",
            "const REAL ln_dot_y",
            "const REAL ln_dot_z",
            "const REAL chi1_dot_x",
            "const REAL chi1_dot_y",
            "const REAL chi1_dot_z",
            "const REAL chi2_dot_x",
            "const REAL chi2_dot_y",
            "const REAL chi2_dot_z",
            "const REAL omega_dot",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
f[LN_X] = ln_dot_x;
f[LN_Y] = ln_dot_y;
f[LN_Z] = ln_dot_z;
f[CHI1_X] = chi1_dot_x;
f[CHI1_Y] = chi1_dot_y;
f[CHI1_Z] = chi1_dot_z;
f[CHI2_X] = chi2_dot_x;
f[CHI2_Y] = chi2_dot_y;
f[CHI2_Z] = chi2_dot_z;
f[OMEGA_PN] = omega_dot;
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

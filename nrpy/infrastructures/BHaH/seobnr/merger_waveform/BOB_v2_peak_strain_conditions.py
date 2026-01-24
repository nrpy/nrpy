"""
Set up C function library for the BOBv2 peak strain conditions.

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
import nrpy.equations.seobnr.BOB_v2_waveform_quantities_kankani_etal as BOB_v2_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_v2_peak_strain_conditions() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the BOBv2 peak strain condition.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the BOBv2 peak strain conditions, 
given peak news time t_p, and reference orbital frequency Omega_0.
The peak strain conditions are given by
d |h_BOB| |
_         |           = 0
dt        |_(t = t_0)

d phi_BOB |
_         |           = omega22NR
dt        |_(t = t_0)

where t_0 is the time at which the strain amplitude |h_BOB| is maximum.
The BOB waveform is computed using the inputs
Binary parameters - m1, m2, chi1, chi2
QNM parameters - omega_qnm, tau_qnm
Reference orbital frequency - Omega_0
Peak strain time - t_0
Peak news time - t_p

Since all parameters except t_p are known, the function evaluates
the derivative of the strain amplitude at t_0 to find the value of t_p
using a root finder.

@params x - The gsl_vector object containing the peak news time and reference orbital frequency.
@params params - The Common data structure containing the model parameters.
@params f - The gsl_vector object to store the peak strain conditions.
@returns - GSL_SUCCESS (0) as required by GSL.
"""
    cfunc_type = "int"
    name = "BOB_v2_peak_strain_conditions"
    #params = (
    #    "const gsl_vector *restrict x , void *restrict params , gsl_vector *restrict f"
    #)
    params = ("double t_val, void *restrict params")
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL t_0 = ((commondata_struct *restrict) params)->t_attach;
const REAL omega_qnm = ((commondata_struct *restrict) params)->omega_qnm;
const REAL tau_qnm = ((commondata_struct *restrict) params)->tau_qnm;
const REAL Omega_0 = gsl_vector_get(x , 1);
const REAL Mf = ((commondata_struct *restrict) params)->Mf;
const REAL a_f = ((commondata_struct *restrict) params)->a_f;

const REAL t_p = 0.0;
const REAL t = t_val;
"""
    wf = BOB_v2_wf.BOB_v2_waveform_quantities()
    body += ccg.c_codegen(
        [wf.t_p_condition, wf.Omega_0_condition],
        ["COMPLEX t_p_condition", "COMPLEX Omega_0_condition"],
        verbose=False,
        include_braces=False,
        fp_type="double complex",
        fp_type_alias="COMPLEX",
    )
    body += """
gsl_vector_set(f , 0 , creal(t_p_condition));
gsl_vector_set(f , 1 , creal(Omega_0_condition));
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="merger_waveform",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

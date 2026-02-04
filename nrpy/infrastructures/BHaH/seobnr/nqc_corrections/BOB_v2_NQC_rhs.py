"""
Register CFunction for BOBv2-informed NQC right hand sides.

Authors:
        Anuj Kankani
        aak00009 **at** mix **dot** wvu **dot** edu
        Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.BOB_v2_waveform_quantities_kankani_etal as BOB_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_v2_NQC_rhs() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_wf.BOB_v2_waveform_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    BOB_code = (
        ccg.c_codegen(
            [
                wf.h_t_attach,
                wf.hdot_t_attach,
                wf.hddot_t_attach,
                wf.w_t_attach,
                wf.wdot_t_attach,
            ],
            [
                "const REAL h_t_attach",
                "const REAL hdot_t_attach",
                "const REAL hddot_t_attach",
                "const REAL w_t_attach",
                "const REAL wdot_t_attach",
            ],
            verbose=False,
            include_braces=False,
        )
        .replace("REAL", "double complex")
        .replace("exp", "cexp")
        .replace("sqrt", "csqrt")
        .replace("pow", "cpow")
        .replace("fabs", "cabs")
        .replace("tanh", "ctanh")
        .replace("sinh", "csinh")
        .replace("cosh", "ccosh")
        .replace("actanh", "catanh")
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the BOBv2-informed Non Quasi-Circular (NQC) right-hand side terms.

@param commondata - Common data structure containing the model parameters.
@param amps - Array to store the calculated amplitudes.
@param omegas - Array to store the calculated angular frequencies.
"""
    cfunc_type = "void"
    name = "BOB_v2_NQC_rhs"
    params = "commondata_struct *restrict commondata , REAL *restrict amps , REAL *restrict omegas"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
const REAL M_f = commondata->M_f;
const REAL a_f = commondata->a_f;
const REAL t_p = commondata->t_p_BOB;
const REAL t_attach = commondata->t_attach;

//compute
"""
    body += BOB_code
    body += """
amps[0] = cabs(h_t_attach);
amps[1] = creal(hdot_t_attach);
amps[2] = creal(hddot_t_attach);
omegas[0] = creal(w_t_attach);
omegas[1] = creal(wdot_t_attach);
"""
    cfc.register_CFunction(
        subdirectory="nqc_corrections",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

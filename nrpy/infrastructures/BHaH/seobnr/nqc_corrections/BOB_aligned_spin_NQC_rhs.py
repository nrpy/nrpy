"""
Register CFunction for BOB-informed NQC right hand sides.

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
import nrpy.equations.seobnr.BOB_aligned_spin_waveform_quantities as BOB_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_aligned_spin_NQC_rhs() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_wf.BOB_aligned_spin_waveform_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    BOB_code = ccg.c_codegen(
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
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the BOB-informed Non Quasi-Circular (NQC) right-hand side terms.

@param commondata - Common data structure containing the model parameters.
@param amps - Array to store the calculated amplitudes.
@param omegas - Array to store the calculated angular frequencies.
"""
    cfunc_type = "void"
    name = "BOB_aligned_spin_NQC_rhs"
    params = "commondata_struct *restrict commondata , REAL *restrict amps , REAL *restrict omegas"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
//compute
"""
    body += BOB_code
    body += """
amps[0] = h_t_attach;
amps[1] = hdot_t_attach;
amps[2] = hddot_t_attach;
omegas[0] = w_t_attach;
omegas[1] = wdot_t_attach;
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

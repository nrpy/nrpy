"""
Register CFunction for BOB-informed merger waveform for single timestep.

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


def register_CFunction_BOB_aligned_spin_waveform() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the (2,2) mode of BOB.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_wf.BOB_aligned_spin_waveform_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    BOB_code = ccg.c_codegen(
        [wf.h, wf.phi],
        ["const REAL h", "const REAL phi"],
        verbose=False,
        include_braces=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the BOB (2,2) mode for a single timestep.

@param t - Time at which to evaluate the waveform.
@param commondata - Common data structure containing the model parameters.
@param waveform - Array to store the calculated waveform.
"""
    cfunc_type = "void"
    name = "BOB_aligned_spin_waveform"
    params = "const REAL t , commondata_struct *restrict commondata , REAL *restrict waveform"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
const REAL t_0 = commondata->t_attach;
//compute
"""
    body += BOB_code
    body += """
waveform[0] = h;
waveform[1] = phi;
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

"""
Set up C function library for native SEOBNRv5 merger-related routines.

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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_merger_quantities as SEOBNRv5_mr
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_merger_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of the SEOBNRv5 merger.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_mr.SEOBNRv5_aligned_spin_merger_quantities()
    SEOBNRv5_code = ccg.c_codegen(
        [wf.h, wf.phi],
        ["const REAL h", "const REAL phi"],
        verbose=False,
        include_braces=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the (2,2) mode of the native SEOBNRv5 merger-ringdown model at a single timestep.

@params t - Time at which to evaluate the waveform.
@params t_0 - Attachment time.
@params h_0 - Amplitude at attachment time.
@params hdot_0 - Amplitude derivative at attachment time.
@params phi_0 - Phase at attachment time.
@params phidot_0 - Angular frequency at attachment time.
@params commondata - Common data structure containing the model parameters.
@params waveform - Array to store the amplitude and phase of the waveform.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_merger_waveform"
    params = "const REAL t , const REAL t_0, const REAL h_0 , const REAL hdot_0 , const REAL phi_0 , const REAL phidot_0, commondata_struct *restrict commondata , REAL *restrict waveform"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
//compute
"""
    body += SEOBNRv5_code
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

"""
Set up C function library for the SEOBNR aligned spin expressions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import ast
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of SEOBNRv5 inspiral waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    hlms_dict = wf.strain()
    hlms = []
    hlms_labels = []
    khatm = []
    khatm_labels = []

    for key in hlms_dict.keys():
        mode = ast.literal_eval(key)
        l, m = mode
        hlms.append(hlms_dict[f"({l} , {m})"])
        hlms_labels.append(f"const double complex h{l}{m}")
        khatm.append(wf.khat[m])
        khatm_labels.append(f"const REAL khat{m}")
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    h_code = ccg.c_codegen(
        hlms,
        hlms_labels,
        verbose=False,
        include_braces=False,
        fp_type="double complex",
        fp_type_alias="COMPLEX",
    )
    khat_code = ccg.c_codegen(
        khatm,
        khatm_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="khat",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Calculates the (2,2) mode of the SEOBNRv5 inspiral waveform for a single timestep.

@param dynamics - Array of dynamical variables.
@param commondata - Common data structure containing the model parameters.
@return - The (2,2) mode of the SEOBNRv5 inspiral waveform.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_waveform"
    params = "REAL *restrict dynamics, commondata_struct *restrict commondata, double complex *inspiral_modes"
    body = """
COMPLEX gamma_22;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL phi = dynamics[PHI];
const REAL Hreal = dynamics[H];
const REAL Omega = dynamics[OMEGA];
const REAL Omega_circ = dynamics[OMEGA_CIRC];
//compute
"""
    body += khat_code
    for key in hlms_dict.keys():
        mode = ast.literal_eval(key)
        l, m = mode
        body += f"""
    gamma_{l}{m} = SEOBNRv5_aligned_spin_gamma_wrapper(1 + 1, -2.*khat{m})"""

    body += h_code
    body += """
inspiral_modes[STRAIN22 - 1] = h22
inspiral_modes[STRAIN21 - 1] = h21
inspiral_modes[STRAIN33 - 1] = h33
inspiral_modes[STRAIN32 - 1] = h32
inspiral_modes[STRAIN44 - 1] = h44
inspiral_modes[STRAIN43 - 1] = h43
inspiral_modes[STRAIN55 - 1] = h55
"""
    cfc.register_CFunction(
        subdirectory="inspiral_waveform",
        includes=includes,
        desc=desc,
        prefunc=prefunc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

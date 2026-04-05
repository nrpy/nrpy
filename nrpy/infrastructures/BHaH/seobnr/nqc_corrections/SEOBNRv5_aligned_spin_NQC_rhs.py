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


def register_CFunction_SEOBNRv5_aligned_spin_NQC_rhs() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for calculating the NQC amplitudes and phase from SEOBNRv5.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_mr.SEOBNRv5_aligned_spin_merger_quantities()
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    SEOBNRv5_code = ccg.c_codegen(
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
Calculate the SEOBNRv5 NR-informed right-hand sides for the Non Quasi-Circular (NQC) corrections.

@params commondata - Common data structure containing the model parameters.
@params amps - Array to store the amplitude and its higher derivatives.
@params omegas - Array to store the angular frequency and its derivative.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_NQC_rhs"
    params = "commondata_struct *restrict commondata , REAL *restrict amps , REAL *restrict omegas"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
//compute
"""
    body += SEOBNRv5_code
    body += """
amps[0] = h_t_attach;
amps[1] = hdot_t_attach;
amps[2] = hddot_t_attach;
omegas[0] = fabs(w_t_attach);
omegas[1] = fabs(wdot_t_attach);
commondata->nr_amp_1 = amps[0];
commondata->nr_amp_2 = amps[1];
commondata->nr_amp_3 = amps[2];
commondata->nr_omega_1 = omegas[0];
commondata->nr_omega_2 = omegas[1];
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

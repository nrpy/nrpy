"""
Set up C function library for the SEOBNR aligned spin expressions.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_constants as SEOBNRv5_const
import nrpy.helpers.parallel_codegen as pcg


def register_Cfunction_SEOBNRv5_aligned_spin_hNR_fits_at_t_attach() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for computing and applying "special" amplitude coefficients needed for (2,1), (4,3), and (5,5) modes.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()

    hNR: List[sp.Expr] = []
    hNR_labels: List[str] = []

    modes = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3), (5, 5)]

    for l, m in modes:
        hNR.append(cast(sp.Expr, const.hNR[f"({l} , {m})"]))
        hNR_labels.append(f"const REAL hNR{l}{m}")

    hNR_code = ccg.c_codegen(
        hNR,
        hNR_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="hNR",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes and applies the special amplitude coefficients to inspiral waveform modes (2,1), (4,3), and (5,5).

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_hNR_fits_at_t_attach"
    params = "commondata_struct *restrict commondata,  REAL *restrict dynamics, REAL *hNR"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL c_21 = commondata->c_21;
const REAL c_43 = commondata->c_43;
const REAL c_55 = commondata->c_55;
const REAL Omega = dynamics[OMEGA];
"""

    body += hNR_code
    body += """
hNR[HNR22] = hNR22;
hNR[HNR21] = hNR21;
hNR[HNR33] = hNR33;
hNR[HNR32] = hNR32;
hNR[HNR44] = hNR44;
hNR[HNR43] = hNR43;
hNR[HNR55] = hNR55;
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
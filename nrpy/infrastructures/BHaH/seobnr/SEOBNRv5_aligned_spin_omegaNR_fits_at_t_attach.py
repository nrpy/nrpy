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


def register_Cfunction_SEOBNRv5_aligned_spin_omegaNR_fits_at_t_attach() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register C function for computing the peak frequencies from numerical relativity.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    const = SEOBNRv5_const.SEOBNR_aligned_spin_constants()

    omegaNR: List[sp.Expr] = []
    omegaNR_labels: List[str] = []

    modes = [(2, 2), (2, 1), (3, 3), (3, 2), (4, 4), (4, 3), (5, 5)]

    for l, m in modes:
        omegaNR.append(cast(sp.Expr, const.omegaNR[f"({l} , {m})"]))
        omegaNR_labels.append(f"const REAL omegaNR{l}{m}")

    omegaNR_code = ccg.c_codegen(
        omegaNR,
        omegaNR_labels,
        verbose=False,
        include_braces=False,
        cse_varprefix="omegaNR",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes the peak frequencies from numerical relativity.

@param commondata - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_omegaNR_fits_at_t_attach"
    params = "commondata_struct *restrict commondata, REAL *omegaNR"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
"""

    body += omegaNR_code
    body += """
omegaNR[HNR22] = omegaNR22;
omegaNR[HNR21] = omegaNR21;
omegaNR[HNR33] = omegaNR33;
omegaNR[HNR32] = omegaNR32;
omegaNR[HNR44] = omegaNR44;
omegaNR[HNR43] = omegaNR43;
omegaNR[HNR55] = omegaNR55;
"""

    cfc.register_CFunction(
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

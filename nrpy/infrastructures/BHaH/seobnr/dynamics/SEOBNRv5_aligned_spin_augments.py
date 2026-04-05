"""
Register CFunction for evaluating SEOBNRv5 Hamiltonian and circular derivatives.

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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_augments() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian and circular derivatives.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """
Evaluates SEOBNRv5 Hamiltonian, instantaneous angular frequency, and circular angular frequency.
These quantities are augmented to the dynamical variables for the calculation of the inspiral modes.

@param commondata - The commondata struct containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_augments"
    params = "commondata_struct *restrict commondata"
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = ccg.c_codegen(
        [H.Hreal, H.dHreal_dpphi, H.dHreal_dpphi_circ],
        ["commondata->Hreal", "commondata->dHreal_dpphi", "commondata->Omega_circ"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return pcg.NRPyEnv()

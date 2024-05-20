"""
Set up C function library for SEOBNR inspiral integrations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par

# Needed during integration and derived from other quantities; do not set!
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "r",
        "m1",
        "m2",
        "a6",
        "dSO",
        "prstar",
        "pphi",
    ],
    [10, 0.5, 0.5, 0.0, 0.0, 10, 0],  # m1, m2, a6, dSO, prstar, pphi
    commondata=True,
    add_to_parfile=False,
)
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "dHreal_dr",
        "dHreal_dprstar",
        "dHreal_dpphi",
        "Hreal",
        "xi",
    ],
    commondata=True,
    add_to_parfile=False,
)

# Here for example:
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "initial_separation",  # ideally this is suffixed with the unit; e.g., initial_separation_km
        "mass_Msun",
        "mass_ratio",
        "chi1",
        "chi2",
    ],
    [20, 100, 0.5, 0.4, -0.3],
    commondata=True,
    add_to_parfile=True,
)


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian by itself.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate SEOBNRv5 Hamiltonian by itself."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_Hamiltonian"
    params = "commondata_struct *restrict commondata"
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = ccg.c_codegen(
        H.Hreal, "commondata->Hreal", verbose=False, include_braces=False
    )
    body += r"""
printf("Hreal = %.15e\n", commondata->Hreal);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_and_derivs() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian and derivatives.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate SEOBNRv5 Hamiltonian and needed derivatives to compute binary dynamics."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_and_derivs"
    params = "commondata_struct *restrict commondata"
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    dHreal_dr = sp.diff(Hq.Hreal, Hq.r)
    dHreal_dprstar = sp.diff(Hq.Hreal, Hq.prstar)
    dHreal_dpphi = sp.diff(Hq.Hreal, Hq.pphi)
    body = ccg.c_codegen(
        [
            dHreal_dr / Hq.nu,
            dHreal_dprstar / Hq.nu,
            dHreal_dpphi / Hq.nu,
            Hq.Hreal,
            Hq.xi,
        ],
        [
            "commondata->dHreal_dr",
            "commondata->dHreal_dprstar",
            "commondata->dHreal_dpphi",
            "commondata->Hreal",
            "commondata->xi",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
printf("dHreal_dr = %.15e\n", commondata->dHreal_dr);
printf("dHreal_dprstar = %.15e\n", commondata->dHreal_dprstar);
printf("dHreal_dpphi = %.15e\n", commondata->dHreal_dpphi);
printf("Hreal = %.15e\n", commondata->Hreal);
printf("xi = %.15e\n", commondata->xi);
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

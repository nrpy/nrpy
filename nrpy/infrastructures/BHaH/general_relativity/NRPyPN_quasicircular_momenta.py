"""
Set up C function for setting 3.5PN quasicircular momenta for binary black holes, using NRPyPN.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Union, cast
from types import FrameType as FT
from inspect import currentframe as cfr

from nrpypn.NRPyPN_shortcuts import (
    m1,
    m2,
    chi1U,
    chi2U,
    r,
    n12U,
    n21U,
    S1U,
    S2U,
    p1U,
    p2U,
)
from nrpypn.PN_p_t import PN_p_t
from nrpypn.PN_p_r import PN_p_r

import nrpy.params as par
import nrpy.c_function as cfc
import nrpy.c_codegen as ccg
import nrpy.helpers.parallel_codegen as pcg

par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "initial_sep",
        "mass_ratio",
        "bbhxy_BH_M_chix",
        "bbhxy_BH_M_chiy",
        "bbhxy_BH_M_chiz",
        "bbhxy_BH_m_chix",
        "bbhxy_BH_m_chiy",
        "bbhxy_BH_m_chiz",
    ],
    [10.0, 1.0, 0, 0, 0, 0, 0, 0],
    commondata=True,
)
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "mass_M",
        "mass_m",
    ],
    [0.5, 0.5],
    commondata=True,
    add_to_parfile=False,
)
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "initial_p_t",
        "initial_p_r",
    ],
    [-1, -1],
    commondata=True,
    add_to_parfile=True,
)


def register_CFunction_NRPyPN_quasicircular_momenta() -> Union[None, pcg.NRPyEnv_type]:
    """Register CFunction for setting quasicircular momenta using NRPyPN."""
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Compute quasicircular momenta using validated expressions from NRPyPN."""
    c_type = "void"
    name = "NRPyPN_quasicircular_momenta"
    params = "commondata_struct *restrict commondata"
    body = r"""// compute quasicircular parameters if commondata.p_t and commondata.p_r not
    // already set in the parfile (i.e., reset from their default values of -1.0).
if(commondata->initial_p_t == -1.0 && commondata->initial_p_r == -1.0) {
    // In NRPyPN, q = m2/m1 = mass_M/mass_m. So mass_m = object 1, and mass_M is object 2.
    const REAL q = commondata->mass_ratio;
    const REAL r = commondata->initial_sep;
    const REAL chi1U0 = commondata->bbhxy_BH_M_chix;
    const REAL chi1U1 = commondata->bbhxy_BH_M_chiy;
    const REAL chi1U2 = commondata->bbhxy_BH_M_chiz;
    const REAL chi2U0 = commondata->bbhxy_BH_m_chix;
    const REAL chi2U1 = commondata->bbhxy_BH_m_chiy;
    const REAL chi2U2 = commondata->bbhxy_BH_m_chiz;

    const REAL mass_M =   q / (1.0 + q);
    const REAL mass_m = 1.0 / (1.0 + q);
    // In NRPyPN, q = m2/m1 = mass_M/mass_m. So mass_m = object 1, and mass_M is object 2.
    const REAL m1 = mass_m;
    const REAL S1U0 = chi1U0 * mass_m*mass_m;
    const REAL S1U1 = chi1U1 * mass_m*mass_m;
    const REAL S1U2 = chi1U2 * mass_m*mass_m;
    const REAL m2 = mass_M;
    const REAL S2U0 = chi2U0 * mass_M*mass_M;
    const REAL S2U1 = chi2U1 * mass_M*mass_M;
    const REAL S2U2 = chi2U2 * mass_M*mass_M;

    REAL Pt, Pr;
"""
    # Compute p_t, the tangential component of momentum
    pt = PN_p_t(m1, m2, chi1U, chi2U, r)
    body += ccg.c_codegen(pt.p_t, "Pt", verbose=False)

    # Compute p_r, the radial component of momentum
    pr = PN_p_r(m1, m2, n12U, n21U, chi1U, chi2U, S1U, S2U, p1U, p2U, r)
    body += ccg.c_codegen(pr.p_r, "Pr", verbose=False)
    body += r"""
commondata->initial_p_t = Pt;
commondata->initial_p_r = Pr;
printf("p_t, p_r = %.15e %.15e\n", Pt, Pr);
}
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

"""
Register CFunction for evaluating the SEOBNRv5 aligned spin flux.

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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_flux() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the SEOBNRv5 aligned spin flux.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """
Evaluates the SEOBNRv5 aligned spin flux for the ODE right-hand sides.

@param y - Array of dynamical variables.
@param Hreal - The real part of the Hamiltonian.
@param Omega - The instantaneous angular frequency of the EOB perturber.
@param Omega_circ - The circular angular frequency of the EOB perturber.
@param f - Array to store the flux.
@param params - Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_flux"
    params = "const REAL *restrict y, const REAL Hreal, const REAL Omega, const REAL Omega_circ, REAL *restrict f, void *restrict params"
    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL c_N = ((commondata_struct *restrict) params)->c_N;
const REAL prstar = y[2];
const REAL pphi = y[3];
"""
    body += ccg.c_codegen(
        [flux],
        [
            "const REAL flux",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
const REAL flux_over_omega = flux / Omega;
f[0] = prstar * flux_over_omega / pphi;
f[1] = flux_over_omega;
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

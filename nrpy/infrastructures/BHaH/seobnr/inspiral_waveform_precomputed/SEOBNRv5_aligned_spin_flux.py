"""
Set up C function library for the SEOBNR aligned spin expressions when using pre-computed waveform coefficients.
The waveform coefficients entering the factorized resummed waveform expressions are functions of the masses and spins.
The individual spins do not evolve when the spins are aligned and we can compute and store the flux coefficients
at the start of the code to accelerate the computation of fluxes and waveforms.


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
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_precomputed_waveform_quantities as SEOBNRv5_precomp_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_aligned_spin_flux() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the SEOBNRv5 aligned spin flux when the waveform coefficients have been precomputed.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """
Evaluates the SEOBNRv5 aligned spin flux.

@param y - Array of dynamical variables.
@param Hreal - The real EOB Hamiltonian.
@param Omega - The instantaneous angular frequency of the EOB perturber.
@param Omega_circ - The circular angular frequency of the EOB perturber.
@param f - Array to store the flux.
@param params - Pointer to the commondata structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_flux"
    params = "const REAL *restrict y, const REAL Hreal, const REAL Omega, const REAL Omega_circ, REAL *restrict f, void *restrict params"
    wf = SEOBNRv5_precomp_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    body = """
commondata_struct *restrict commondata = (commondata_struct *restrict) params;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL prstar = y[2];
const REAL pphi = y[3];
"""
    flux_ccode = ccg.c_codegen(
        [flux],
        [
            "const REAL flux",
        ],
        verbose=False,
        include_braces=False,
    )
    body += flux_ccode

    body += """
const REAL flux_over_omega = flux / Omega;
f[0] = prstar * flux_over_omega / pphi;
f[1] = flux_over_omega;
"""
    cfc.register_CFunction(
        subdirectory="inspiral_waveform_precomputed",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

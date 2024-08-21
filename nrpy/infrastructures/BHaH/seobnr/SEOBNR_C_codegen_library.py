"""
Set up C function library for SEOBNR inspiral integrations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_Hamiltonian as SEOBNRv5_Ham
import nrpy.equations.seobnr.SEOBNRv5_aligned_spin_waveform_quantities as SEOBNRv5_wf
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par

# Needed during integration and derived from other quantities; do not set!
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "r",
        "phi",
        "m1",
        "m2",
        "a6",
        "dSO",
        "prstar",
        "pphi",
    ],
    [20, 0, 0.5, 0.5, 0.0, 0.0, 0.0, 3.3],  # r, phi, m1, m2, a6, dSO, prstar, pphi
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
        "dHreal_dr_dr",
        "dHreal_dr_dpphi",
        "dHreal_dr_circ",
        "dHreal_dpphi_circ",
        "Hreal",
        "xi",
        "flux",
    ],
    commondata=True,
    add_to_parfile=False,
)

# This is sufficient for initial conditions. Order is the same as pySEOBNR.
par.register_CodeParameters(
    "REAL",
    __name__,
    ["mass_ratio", "chi1", "chi2", "initial_omega"],
    [
        1,
        0.4,
        -0.3,
        0.01118,
    ],  # mass_ratio convention is m_greater/m_lesser, initial_omega chosen for r ~ 20M
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


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_coefficients() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the masses and Hamiltonian coefficients.
    The inputs needed to generate an EOB approximant are mass ratio, spins and
    an initial orbital frequency. Therefore, one needs to compute the individual
    masses from the mass ratio and the Hamiltonian coeffients which are a function
    of mass ratio and spins.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate the SEOBNRv5 Hamiltonian coefficients."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_coefficients"
    params = "commondata_struct *restrict commondata"
    H = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = """
const REAL q = commondata->mass_ratio;
commondata->m1 = q / (1.0 + q);
commondata->m2 = 1.0 / (1.0 + q);
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
"""
    body += ccg.c_codegen(
        [H.pyseobnr_a6, H.pyseobnr_dSO],
        ["commondata->a6", "commondata->dSO"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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
    dHreal_dr = Hq.dHreal_dr
    dHreal_dprstar = Hq.dHreal_dprstar
    dHreal_dpphi = Hq.dHreal_dpphi
    body = ccg.c_codegen(
        [
            dHreal_dr,
            dHreal_dprstar,
            dHreal_dpphi,
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


def register_CFunction_SEOBNRv5_aligned_spin_initial_conditions_conservative() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5's conservative initial conditions using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
"""
    desc = """Evaluate the SEOBNRv5 conservative initial conditions."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_initial_conditions_conservative"
    params = "commondata_struct *restrict commondata"
    body = """
const gsl_multiroot_fsolver_type *T;
gsl_multiroot_fsolver *s;
int status;
size_t i , iter = 0;
const size_t n = 2;
const int maxiter = 100;
gsl_multiroot_function f = {&SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit, n , commondata};
gsl_vector *x = gsl_vector_alloc(n);
REAL omega = commondata->initial_omega;
REAL pphi = pow(omega,-1./3.);
REAL r = pphi*pphi;
gsl_vector_set(x , 0 , r);
gsl_vector_set(x , 1 , pphi);
T = gsl_multiroot_fsolver_hybrids;
s = gsl_multiroot_fsolver_alloc(T , 2);
gsl_multiroot_fsolver_set(s , &f , x);
do {
  iter++;
  status = gsl_multiroot_fsolver_iterate (s);
  if (status){
    break;
  }
  status = gsl_multiroot_test_residual (s->f, 6e-12);
}
while(status == GSL_CONTINUE && iter < maxiter);
commondata->r = gsl_vector_get(s->x , 0);
commondata->pphi = gsl_vector_get(s->x, 1);
gsl_multiroot_fsolver_free (s);
gsl_vector_free (x);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating SEOBNRv5 Hamiltonian circular orbit conditions.
    The CFunction is constructed for use with gsl/multiroots.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
"""
    desc = """Evaluate SEOBNRv5 Hamiltonian's circular orbit conditions to compute conservative initial conditions."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_Hamiltonian_circular_orbit"
    params = "const gsl_vector * x , void *params , gsl_vector * f"  # This convention is from the gsl multi-dimensional root finder example
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    dHreal_dr_circ = Hq.dHreal_dr_circ
    dHreal_dpphi_circ = Hq.dHreal_dpphi_circ
    body = r"""
const REAL r = gsl_vector_get(x , 0);
const REAL pphi = gsl_vector_get(x , 1);
const REAL m1 = ((commondata_struct *) params)->m1;
const REAL m2 = ((commondata_struct *) params)->m2;
const REAL chi1 = ((commondata_struct *) params)->chi1;
const REAL chi2 = ((commondata_struct *) params)->chi2;
const REAL a6 = ((commondata_struct *) params)->a6;
const REAL dSO = ((commondata_struct *) params)->dSO;
const REAL initial_omega = ((commondata_struct *) params)->initial_omega;
"""
    body += ccg.c_codegen(
        [
            dHreal_dr_circ,
            dHreal_dpphi_circ,
        ],
        [
            "const REAL dHreal_dr_circ",
            "const REAL dHreal_dpphi_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
gsl_vector_set(f , 0 , dHreal_dr_circ);
gsl_vector_set(f , 1 , dHreal_dpphi_circ - initial_omega);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_flux() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating SEOBNRv5 factorized resummed flux.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate SEOBNRv5 flux."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_flux"
    params = "commondata_struct *restrict commondata"
    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    body = """
REAL Omega = commondata->dHreal_dpphi;
REAL Omega_circ = commondata->dHreal_dpphi_circ;
"""
    body += ccg.c_codegen(
        [
            flux,
        ],
        [
            "commondata->flux",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
printf("flux = %.15e\n", commondata->flux);
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

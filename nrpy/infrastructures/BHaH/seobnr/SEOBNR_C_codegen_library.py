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
        "t_stepback",
    ],
    [
        20,
        0,
        0.5,
        0.5,
        0.0,
        0.0,
        0.0,
        3.3,
        250.0,
    ],  # r, phi, m1, m2, a6, dSO, prstar, pphi
    commondata=True,
    add_to_parfile=False,
)

par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "Delta_t",
        "dT",
        "t_ISCO",
        "t_attach",
        "omega_qnm",
        "tau_qnm",
        "a_f",
        "M_f",
        "a_1_NQC",
        "a_2_NQC",
        "a_3_NQC",
        "b_1_NQC",
        "b_2_NQC",
        "r_stop",
        "r_ISCO",
    ],
    [
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ],
    commondata=True,
    add_to_parfile=False,
)


par.register_CodeParameters(
    "REAL *restrict",
    __name__,
    [
        "dynamics_low",
        "dynamics_fine",
        "waveform_low",
        "waveform_fine",
        "waveform_inspiral",
        "waveform_IMR",
    ],
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
)


par.register_CodeParameters(
    "size_t",
    __name__,
    ["nsteps_low", "nsteps_fine", "nsteps_inspiral", "nsteps_IMR"],
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
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
        "Omega_circ",
    ],
    commondata=True,
    add_to_parfile=False,
)

# This is sufficient for initial conditions. Order is the same as pySEOBNR.
par.register_CodeParameters(
    "REAL",
    __name__,
    ["mass_ratio", "chi1", "chi2", "initial_omega", "total_mass", "dt"],
    [
        1,
        0.4,
        -0.3,
        0.01118,
        50,
        2.4627455127717882e-05,
    ],  # mass_ratio convention is m_greater/m_lesser, initial_omega chosen for r ~ 20M
    commondata=True,
    add_to_parfile=True,
)


def register_CFunction_SEOBNRv5_aligned_spin_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for calculating the (2,2) mode of SEOBNRv5.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    hlms = wf.strain()
    h22 = hlms["(2 , 2)"]
    # We are going to be doing this twice;
    # once for the fine dynamics and once for the coarse.
    h22_code = (
        ccg.c_codegen(
            h22,
            ["const double complex h22"],
            verbose=False,
            include_braces=False,
        )
        .replace("REAL", "double complex")
        .replace("exp", "cexp")
    )
    khat2_code = ccg.c_codegen(
        [wf.khat[2]],
        [
            "const REAL khat2",
        ],
        verbose=False,
        include_braces=False,
        cse_varprefix="khat",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the SEOBNRv5 22 mode."""
    cfunc_type = "int"
    prefunc = "#include<complex.h>"
    name = "SEOBNRv5_aligned_spin_waveform"
    params = "REAL *restrict dynamics, REAL *restrict waveform, commondata_struct *restrict commondata"
    body = """
const int dyn_size = 8; //t,r,phi,prstar,pphi,Hreal, Omega, Omega_circ
int status = 0;
int i;
REAL gamma_real_22 , gamma_imag_22;
REAL gamma_22[2];
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL a6 = commondata->a6;
const REAL dSO = commondata->dSO;
const REAL t = dynamics[TIME];
const REAL r = dynamics[R];
const REAL phi = dynamics[PHI];
const REAL prstar = dynamics[PRSTAR];
const REAL pphi = dynamics[PPHI];
const REAL Hreal = dynamics[H];
const REAL Omega = dynamics[OMEGA];
const REAL Omega_circ = dynamics[OMEGA_CIRC];
//compute
"""
    body += khat2_code
    body += """
  status = SEOBNRv5_aligned_spin_gamma_wrapper(3.,-2.*khat2,gamma_22);
  gamma_real_22 = gamma_22[0];
  gamma_imag_22 = gamma_22[1];
"""
    body += h22_code
    body += """
waveform[0] = (REAL) creal(h22);
waveform[1] = (REAL) cimag(h22);
return GSL_SUCCESS;
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for computing the (2,2) mode of SEOBNRv5 given the dynamics.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the SEOBNRv5 22 mode."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_waveform_from_dynamics"
    params = "commondata_struct *restrict commondata"
    body = """
int status = 0;
int i;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL a6 = commondata->a6;
const REAL dSO = commondata->dSO;
REAL dynamics[NUMVARS] , waveform[2];
commondata->waveform_low = (REAL *)malloc(commondata->nsteps_low*NUMMODES*sizeof(REAL)); //t , h_+ , h_x
commondata->waveform_fine = (REAL *)malloc(commondata->nsteps_fine*NUMMODES*sizeof(REAL)); //t , h_+ , h_x

//low sampling
for (i = 0; i < commondata->nsteps_low; i++) {
  //assign
  dynamics[TIME] = commondata->dynamics_low[IDX(i,TIME)];
  dynamics[R] = commondata->dynamics_low[IDX(i,R)];
  dynamics[PHI] = commondata->dynamics_low[IDX(i,PHI)];
  dynamics[PRSTAR] = commondata->dynamics_low[IDX(i,PRSTAR)];
  dynamics[PPHI] = commondata->dynamics_low[IDX(i,PPHI)];
  dynamics[H] = commondata->dynamics_low[IDX(i,H)];
  dynamics[OMEGA] = commondata->dynamics_low[IDX(i,OMEGA)];
  dynamics[OMEGA_CIRC] = commondata->dynamics_low[IDX(i,OMEGA_CIRC)];
  
  //compute
  status = SEOBNRv5_aligned_spin_waveform(dynamics, waveform, commondata);
  //store
  commondata->waveform_low[IDX_WF(i,TIME)] = dynamics[TIME];
  commondata->waveform_low[IDX_WF(i,HPLUS)] = waveform[0];
  commondata->waveform_low[IDX_WF(i,HCROSS)] = -1.0*waveform[1]; // polarizations are described as h = h_+ - I*h_x
}
//high sampling
for (i = 0; i < commondata->nsteps_fine; i++) {
  //assign
  dynamics[TIME] = commondata->dynamics_fine[IDX(i,TIME)];
  dynamics[R] = commondata->dynamics_fine[IDX(i,R)];
  dynamics[PHI] = commondata->dynamics_fine[IDX(i,PHI)];
  dynamics[PRSTAR] = commondata->dynamics_fine[IDX(i,PRSTAR)];
  dynamics[PPHI] = commondata->dynamics_fine[IDX(i,PPHI)];
  dynamics[H] = commondata->dynamics_fine[IDX(i,H)];
  dynamics[OMEGA] = commondata->dynamics_fine[IDX(i,OMEGA)];
  dynamics[OMEGA_CIRC] = commondata->dynamics_fine[IDX(i,OMEGA_CIRC)];
  
  //compute
  status = SEOBNRv5_aligned_spin_waveform(dynamics, waveform, commondata);
  //store
  commondata->waveform_fine[IDX_WF(i,TIME)] = dynamics[TIME];
  commondata->waveform_fine[IDX_WF(i,HPLUS)] = waveform[0];
  commondata->waveform_fine[IDX_WF(i,HCROSS)] = -1.0*waveform[1]; // polarizations are described as h = h_+ - I*h_x
}
return GSL_SUCCESS;
"""
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


def register_CFunction_SEOBNRv5_aligned_spin_flux() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the SEOBNRv5 aligned spin flux.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    desc = """Evaluate the SEOBNRv5 aligned spin flux."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_flux"
    params = "const REAL *restrict y, const REAL Hreal, const REAL Omega, const REAL Omega_circ, REAL *restrict f, void *restrict params"
    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL r = y[0];
const REAL phi = y[1];
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
return GSL_SUCCESS;
"""
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


def register_CFunction_SEOBNRv5_aligned_spin_right_hand_sides() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the right hand sides for the SEOBNRv5 equations of motion.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate SEOBNRv5 Hamiltonian and needed derivatives to compute binary dynamics."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_right_hand_sides"
    params = "REAL t, const REAL *restrict y, REAL *restrict f, void *restrict params"
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    body = """
const REAL m1 = ((commondata_struct *restrict) params)->m1;
const REAL m2 = ((commondata_struct *restrict) params)->m2;
const REAL chi1 = ((commondata_struct *restrict) params)->chi1;
const REAL chi2 = ((commondata_struct *restrict) params)->chi2;
const REAL a6 = ((commondata_struct *restrict) params)->a6;
const REAL dSO = ((commondata_struct *restrict) params)->dSO;
const REAL r = y[0];
const REAL phi = y[1];
const REAL prstar = y[2];
const REAL pphi = y[3];
"""
    body += ccg.c_codegen(
        [
            Hq.Hreal,
            Hq.xi,
            Hq.dHreal_dr,
            Hq.dHreal_dprstar,
            Hq.dHreal_dpphi,
            Hq.dHreal_dpphi_circ,
        ],
        [
            "const REAL Hreal",
            "const REAL xi",
            "const REAL dHreal_dr",
            "const REAL dHreal_dprstar",
            "const REAL dHreal_dpphi",
            "const REAL dHreal_dpphi_circ",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
REAL flux[2];
SEOBNRv5_aligned_spin_flux(y,Hreal,dHreal_dpphi,dHreal_dpphi_circ,flux,params);
f[0] = xi * dHreal_dprstar;
f[1] = dHreal_dpphi;
f[2] = -xi * dHreal_dr + flux[0];
f[3] = flux[1];
return GSL_SUCCESS;
"""
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

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
    "REAL *restrict",
    __name__,
    [
        "dynamics_low",
        "dynamics_fine",
        "dynamics_combined",
        "waveform_low",
        "waveform_fine",
        "waveform_combined",
    ],
    commondata=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=False,
)

par.register_CodeParameters(
    "size_t",
    __name__,
    ["nsteps_low", "nsteps_fine", "nsteps_combined"],
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


def register_CFunction_SEOBNRv5_aligned_spin_waveform_from_dynamics() -> (
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
    # It makes sense to save the
    h22_code = ccg.c_codegen(
        [h22],
        [
            "h22",
        ],
        verbose=False,
        include_braces=False,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Calculate the SEOBNRv5 22 mode."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_waveform_from_dynamics"
    params = "commondata_struct *restrict commondata"
    body = """
const int dyn_size = 8; //t,r,phi,prstar,pphi,Hreal, Omega, Omega_circ
int status = 0;
int i;
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL a6 = commondata->a6;
const REAL dSO = commondata->dSO;
REAL t , r , phi , prstar, pphi , Hreal , Omega , Omega_circ;
commondata->waveform_low = (REAL *)calloc(commondata->nsteps_low*3,sizeof(REAL)); //t , h_+ , h_x
commondata->waveform_high = (REAL *)calloc(commondata->nsteps_high*3,sizeof(REAL)); //t , h_+ , h_x

//low sampling
for (i = 0; i < commondata->nsteps_low; i++) {
  //assign
  t = commondata->dynamics_low[8*i + 0];
  r = commondata->dynamics_low[8*i + 1];
  phi = commondata->dynamics_low[8*i + 2];
  prstar = commondata->dynamics_low[8*i + 3];
  pphi = commondata->dynamics_low[8*i + 4];
  Hreal = commondata->dynamics_low[8*i + 5];
  Omega = commondata->dynamics_low[8*i + 6];
  Omega_circ = commondata->dynamics_low[8*i + 7];
  
  //compute
"""
    body += h22_code
    body += """
  //store
  commondata->waveform_low[3*i + 0] = t;
  commondata->waveform_low[3*i + 1] = (REAL)creal(h22);
  commondata->waveform_low[3*i + 1] = -1*(REAL)cimag(h22); // polarizations are described as h = h_+ - I*h_x
}
"""
    body += """
//high sampling
for (i = 0; i < commondata->nsteps_high; i++) {
  //assign
  t = commondata->dynamics_high[8*i + 0];
  r = commondata->dynamics_high[8*i + 1];
  phi = commondata->dynamics_high[8*i + 2];
  prstar = commondata->dynamics_high[8*i + 3];
  pphi = commondata->dynamics_high[8*i + 4];
  Hreal = commondata->dynamics_high[8*i + 5];
  Omega = commondata->dynamics_high[8*i + 6];
  Omega_circ = commondata->dynamics_high[8*i + 7];
  
  //compute
"""
    body += h22_code
    body += """
  //store
  commondata->waveform_high[3*i + 0] = t;
  commondata->waveform_high[3*i + 1] = (REAL)creal(h22);
  commondata->waveform_high[3*i + 1] = -1*(REAL)cimag(h22); // polarizations are described as h = h_+ - I*h_x
}
"""
    body += """
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

    includes = ["BHaH_defines.h"]
    desc = """Evaluate SEOBNRv5 Hamiltonian and needed derivatives to compute binary dynamics."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_right_hand_sides"
    params = "REAL t, const REAL *restrict y, REAL *restrict f, void *restrict params"
    Hq = SEOBNRv5_Ham.SEOBNRv5_aligned_spin_Hamiltonian_quantities()
    wf = SEOBNRv5_wf.SEOBNRv5_aligned_spin_waveform_quantities()
    flux = wf.flux()
    flux = (
        flux.subs(wf.Hreal, Hq.Hreal)
        .subs(wf.Omega, Hq.dHreal_dpphi)
        .subs(wf.Omega_circ, Hq.dHreal_dpphi_circ)
    )
    r_dot = Hq.xi * Hq.dHreal_dprstar
    phi_dot = Hq.dHreal_dpphi
    pphi_dot = flux / Hq.nu
    prstar_dot = -Hq.xi * Hq.dHreal_dr + Hq.prstar * pphi_dot / Hq.pphi
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
        [r_dot, phi_dot, prstar_dot, pphi_dot],
        [
            "const REAL rdot",
            "const REAL phidot",
            "const REAL prstardot",
            "const REAL pphidot",
        ],
        verbose=False,
        include_braces=False,
    )
    body += r"""
f[0] = rdot;
f[1] = phidot;
f[2] = prstardot;
f[3] = pphidot;
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

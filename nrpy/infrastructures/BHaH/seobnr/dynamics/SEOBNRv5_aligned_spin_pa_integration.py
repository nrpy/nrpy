"""
Register CFunction for integrating the SEOBNRv5 post-adiabatic equations of motion using GSL.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par


def register_CFunction_SEOBNRv5_aligned_spin_pa_integration() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for integrating the SEOBNRv5 post-adiabatic equations of motion using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "dprstar_dr",
            "dpphi_dr",
        ],
        commondata=True,
        add_to_parfile=False,
    )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
"""
    desc = """
Integrates the SEOBNRv5 post-adiabatic equations of motion to obtain the dynamics of the EOB perturber
as well as the augmented dynamical quantities for generating the inspiral waveform.

@param commondata - Common data struct containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_pa_integration"
    params = "commondata_struct *restrict commondata"
    body = """
const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL nu = m1 * m2 / ((m1 + m2) * (m1 + m2));
const REAL chi_eff = (m1 * chi1 + m2 * chi2) / (m1 + m2);
const REAL r_final_prefactor = 2.7 + chi_eff * (1.0 - 4.0 * nu);
REAL final_r = fmax(10.0, r_final_prefactor * commondata->r_ISCO);
REAL initial_r = commondata->r;
// fiducial dr = .3
REAL dr = .3;
size_t nsteps = (size_t)((initial_r - final_r) / dr + 1);
// update nsteps and dr to ensure a minimum of 10 steps
nsteps = nsteps > 10 ? nsteps : 10;
dr = (initial_r - final_r) / (nsteps - 1);
REAL *restrict r = (REAL *)malloc(nsteps * sizeof(REAL));
if (r == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for r\\n");
  exit(1);
}
REAL *restrict pphi = (REAL *)malloc(nsteps * sizeof(REAL));
if (pphi == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for pphi\\n");
  exit(1);
}
REAL *restrict prstar = (REAL *)malloc(nsteps * sizeof(REAL));
if (prstar == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for prstar\\n");
  exit(1);
}
REAL *restrict dpphi_dr = (REAL *)malloc(nsteps * sizeof(REAL));
if (dpphi_dr == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for dpphi_dr\\n");
  exit(1);
}
REAL *restrict dprstar_dr = (REAL *)malloc(nsteps * sizeof(REAL));
if (dprstar_dr == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for dprstar_dr\\n");
  exit(1);
}
REAL *restrict dt_dr = (REAL *)malloc(nsteps * sizeof(REAL));
if (dt_dr == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for dt_dr\\n");
  exit(1);
}
REAL *restrict dphi_dr = (REAL *)malloc(nsteps * sizeof(REAL));
if (dphi_dr == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for dphi_dr\\n");
  exit(1);
}

gsl_function F_pphi0;
F_pphi0.function = &SEOBNRv5_aligned_spin_pphi0_equation;
F_pphi0.params = commondata;
REAL x_lo, x_hi;
gsl_function F_pphi;
F_pphi.function = &SEOBNRv5_aligned_spin_pphi_equation;
F_pphi.params = commondata;
gsl_function F_prstar;
F_prstar.function = &SEOBNRv5_aligned_spin_pr_equation;
F_prstar.params = commondata;

// perform zero PA integration
for (size_t i = 0; i < nsteps; i++) {
  r[i] = initial_r - dr * i;
  prstar[i] = 0.;
  dprstar_dr[i] = 0.;
  commondata->r = r[i];
  x_lo = sqrt(r[i]) * 0.6;
  x_hi = sqrt(r[i]) * 1.4;
  F_pphi0.params = commondata;
  pphi[i] = root_finding_1d(x_lo, x_hi, &F_pphi0);
}
// only compute dpphi_dr for zero PA integration
dy_dx(pphi, r, dpphi_dr, nsteps);
// perform PA integration
// fiducial PA order = 8
const int PA_ORDER = 8;
for (int j = 1; j <= PA_ORDER; j++) {
  for (size_t i = 0; i < nsteps; i++) {
    commondata->r = r[i];
    commondata->pphi = pphi[i];
    commondata->prstar = prstar[i];
    commondata->dpphi_dr = dpphi_dr[i];
    commondata->dprstar_dr = dprstar_dr[i];
    switch (j % 2) {
    case 0:
      // even order, compute pphi
      x_lo = pphi[i] * 0.95;
      x_hi = pphi[i] * 1.05;
      pphi[i] = root_finding_1d(x_lo, x_hi, &F_pphi);
      break;
    case 1:
      // odd order, compute prstar
      // note that prstar is always negative (inspiral)
      // Since upper and lower bounds are a 5% window of current prstar
      // so upper bound = 0.95*prstar and lower bound = 1.05*prstar
      // If prstar ~ 0, we can use the bounds from the ODE initial conditions
      x_hi = fabs(prstar[i]) > 1e-14 ? prstar[i] * 0.95 : 0.;
      x_lo = fabs(prstar[i]) > 1e-14 ? prstar[i] * 1.05 : -3e-2;
      prstar[i] = root_finding_1d(x_lo, x_hi, &F_prstar);
      break;
    default:
      break;
    }
  }
  // compute dprstar/dr and dpphi/dr
  dy_dx(prstar, r, dprstar_dr, nsteps);
  dy_dx(pphi, r, dpphi_dr, nsteps);
}
// compute derivatives of time and phase
for (size_t i = 0; i < nsteps; i++) {
  commondata->pphi = pphi[i];
  commondata->prstar = prstar[i];
  dt_dr[i] = SEOBNRv5_aligned_spin_t_equation(r[i], commondata);
  dphi_dr[i] = SEOBNRv5_aligned_spin_phi_equation(r[i], commondata);
}
// compute cumulative integrals
REAL *restrict t = (REAL *)malloc(nsteps * sizeof(REAL));
if (t == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for t\\n");
  exit(1);
}
REAL *restrict phi = (REAL *)malloc(nsteps * sizeof(REAL));
if (phi == NULL) {
  fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for phi\\n");
  exit(1);
}
cumulative_integration(dt_dr, r, t, nsteps);
cumulative_integration(dphi_dr, r, phi, nsteps);


// run the ODE integration with the final PA values as initial conditions
commondata->r = r[nsteps - 1];
commondata->phi = phi[nsteps - 1];
commondata->prstar = prstar[nsteps - 1];
commondata->pphi = pphi[nsteps - 1];
SEOBNRv5_aligned_spin_ode_integration(commondata);
// Store the dynamics array.
// In the case where the ODE integration is needed
// for less than 250 M, dynamics->low is never initialized
// Also note that ODE integrator contains the zeroth timestep 
// which is the last timestep of PA.
// We will not include the last PA timestep to ensure there is no
// double counting of trajectory points
if (commondata->dynamics_low == NULL) {
  commondata->nsteps_low = nsteps - 1;
  commondata->dynamics_low = (REAL *)malloc(commondata->nsteps_low * NUMVARS * sizeof(REAL));
  if (commondata->dynamics_low == NULL) {
    fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for commondata->dynamics_low\\n");
    exit(1);
  }
  for (size_t i = 0; i < commondata->nsteps_low; i++) {
    commondata->r = r[i];
    commondata->phi = phi[i];
    commondata->prstar = prstar[i];
    commondata->pphi = pphi[i];
    SEOBNRv5_aligned_spin_augments(commondata);
    commondata->dynamics_low[IDX(i, TIME)] = t[i];
    commondata->dynamics_low[IDX(i, R)] = r[i];
    commondata->dynamics_low[IDX(i, PHI)] = phi[i];
    commondata->dynamics_low[IDX(i, PRSTAR)] = prstar[i];
    commondata->dynamics_low[IDX(i, PPHI)] = pphi[i];
    commondata->dynamics_low[IDX(i, H)] = commondata->Hreal;
    commondata->dynamics_low[IDX(i, OMEGA)] = commondata->dHreal_dpphi;
    commondata->dynamics_low[IDX(i, OMEGA_CIRC)] = commondata->Omega_circ;
  }
} else {
  size_t nsteps_ode_dynamics_low = commondata->nsteps_low;
  REAL *ode_dynamics_low = (REAL *)malloc(nsteps_ode_dynamics_low * NUMVARS * sizeof(REAL));
  if (ode_dynamics_low == NULL) {
    fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, malloc failed for ode_dynamics_low\\n");
    exit(1);
  }
  memcpy(ode_dynamics_low, commondata->dynamics_low, nsteps_ode_dynamics_low * NUMVARS * sizeof(REAL));
  commondata->nsteps_low += nsteps - 1; // recall, we don't use the last timestep of the PA dynamics as it is the first timestep of ODE dynamics
  commondata->dynamics_low = (REAL *)realloc(commondata->dynamics_low, commondata->nsteps_low * NUMVARS * sizeof(REAL));
  if (commondata->dynamics_low == NULL) {
    fprintf(stderr, "Error: in SEOBNRv5_aligned_spin_pa_integration, realloc failed for commondata->dynamics_low\\n");
    exit(1);
  }
  for (size_t i = 0; i < nsteps - 1; i++) {
    commondata->r = r[i];
    commondata->phi = phi[i];
    commondata->prstar = prstar[i];
    commondata->pphi = pphi[i];
    SEOBNRv5_aligned_spin_augments(commondata);
    commondata->dynamics_low[IDX(i, TIME)] = t[i];
    commondata->dynamics_low[IDX(i, R)] = r[i];
    commondata->dynamics_low[IDX(i, PHI)] = phi[i];
    commondata->dynamics_low[IDX(i, PRSTAR)] = prstar[i];
    commondata->dynamics_low[IDX(i, PPHI)] = pphi[i];
    commondata->dynamics_low[IDX(i, H)] = commondata->Hreal;
    commondata->dynamics_low[IDX(i, OMEGA)] = commondata->dHreal_dpphi;
    commondata->dynamics_low[IDX(i, OMEGA_CIRC)] = commondata->Omega_circ;
  }
  memcpy(commondata->dynamics_low + (nsteps - 1) * NUMVARS, ode_dynamics_low, nsteps_ode_dynamics_low * NUMVARS * sizeof(REAL));
  // add the time offset to the dynamics_low
  for (size_t i = nsteps - 1; i < commondata->nsteps_low; i++) {
    commondata->dynamics_low[IDX(i, TIME)] += t[nsteps - 1];
  }
  free(ode_dynamics_low);
}
// add the time offset to the dynamics_fine.
// Note that this works regardless of the ODE low dynamics since the
// last PA timestep is always removed.
for (size_t i = 0; i < commondata->nsteps_fine; i++) {
  commondata->dynamics_fine[IDX(i, TIME)] += t[nsteps - 1];
}
// free memory
free(r);
free(pphi);
free(prstar);
free(dpphi_dr);
free(dprstar_dr);
free(dt_dr);
free(dphi_dr);
free(t);
free(phi);
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        prefunc=prefunc,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

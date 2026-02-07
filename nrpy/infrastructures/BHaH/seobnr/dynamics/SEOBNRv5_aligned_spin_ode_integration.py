"""
Register CFunction for integrating the SEOBNRv5 equations of motion using GSL.

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


def register_CFunction_SEOBNRv5_aligned_spin_ode_integration() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for integrating the SEOBNRv5 equations of motion using GSL.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
"""
    desc = """
Integrates the SEOBNRv5 equations of motion to obtain the dynamics of the EOB perturber
as well as the augmented dynamical quantities for generating the inspiral waveform.

@param commondata - Common data struct containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_ode_integration"
    params = "commondata_struct *restrict commondata"
    body = """
const gsl_odeiv2_step_type *restrict T = gsl_odeiv2_step_rk8pd;
gsl_odeiv2_step *restrict s = gsl_odeiv2_step_alloc(T, 4);
gsl_odeiv2_control *restrict c = gsl_odeiv2_control_standard_new(1e-12, 1e-11, 1.0, 1.0);
gsl_odeiv2_system sys = {SEOBNRv5_aligned_spin_right_hand_sides, NULL, 4, commondata};
gsl_odeiv2_evolve *restrict e = gsl_odeiv2_evolve_alloc (4);

REAL t = 0.0;
const REAL tmax = 2e9;
REAL y[4], dydt[4];
int status = 0;
int i;
int stop = 0;
y[0] = commondata->r;
y[1] = commondata->phi;
y[2] = commondata->prstar;
y[3] = commondata->pphi;
status = SEOBNRv5_aligned_spin_right_hand_sides(t, y, dydt, commondata);
int rhs_status[1] = {GSL_SUCCESS};
char rhs_name[] = "gsl_odeiv2_evolve_apply";
REAL h = 2.0 * M_PI / dydt[1] / 5.0;
size_t bufferlength = (size_t)(tmax / h); // runs up to 0.01x maximum time (we should not ideally run that long)
REAL *restrict dynamics_RK = (REAL *)malloc(bufferlength * (NUMVARS)*sizeof(REAL));
if (dynamics_RK == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for dynamics_RK\\n");
  exit(1);
}
size_t nsteps = 0;

// store
dynamics_RK[IDX(nsteps,TIME)] = t;
for (i = 1; i < 5; i++) {
  dynamics_RK[IDX(nsteps,i)] = y[i - 1];
}
SEOBNRv5_aligned_spin_augments(commondata);
dynamics_RK[IDX(nsteps,H)] = commondata->Hreal;
dynamics_RK[IDX(nsteps,OMEGA)] = commondata->dHreal_dpphi;
dynamics_RK[IDX(nsteps,OMEGA_CIRC)] = commondata->Omega_circ;
nsteps++;
REAL Omega_previous = commondata->initial_omega;

while (t < tmax && stop == 0) {
  // integrate
  status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, tmax, &h, y);
  handle_gsl_return_status(status,rhs_status,1,rhs_name);
  commondata->r = y[0];
  commondata->phi = y[1];
  commondata->prstar = y[2];
  commondata->pphi = y[3];
  SEOBNRv5_aligned_spin_augments(commondata);

  // buffercheck
  if (nsteps >= bufferlength) {
    bufferlength = 2 * bufferlength;
    dynamics_RK = (REAL *)realloc(dynamics_RK, bufferlength * (NUMVARS) * sizeof(REAL));
    if (dynamics_RK == NULL){
      fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), realloc() failed for dynamics_RK\\n");
      exit(1);
    }
  }

  // update
  // store
  dynamics_RK[IDX(nsteps,TIME)] = t;
  for (i = 1; i < 5; i++) {
    dynamics_RK[IDX(nsteps,i)] = y[i-1];
  }
  dynamics_RK[IDX(nsteps,H)] = commondata->Hreal;
  dynamics_RK[IDX(nsteps,OMEGA)] = commondata->dHreal_dpphi;
  dynamics_RK[IDX(nsteps,OMEGA_CIRC)] = commondata->Omega_circ;
  nsteps++;

  // stopcheck
  if (commondata->r < 6.0) {
    // decrease in frequency: Omega peak stop = index of omega
    if (commondata->dHreal_dpphi < Omega_previous) {
      stop = OMEGA;
      break;
    }
    SEOBNRv5_aligned_spin_right_hand_sides(t, y, dydt, commondata);
    // outspiral dPR: PR peak stop = index of pr
    if (dydt[2] > 0.0) {
      stop = PRSTAR;
      break;
    }
    // outspiral dR
    if (dydt[0] > 0.0) {
      break;
    }
    // Stopping radius
    if (commondata->r < commondata->r_stop){
      break;
    }
    // unphysical frequency
    if (commondata->r < 3.0 && commondata->Omega_circ > 1.0) {
      break;
    }
  }
  Omega_previous = commondata->dHreal_dpphi;
}

// free up gsl ode solver
gsl_odeiv2_control_free(c);
gsl_odeiv2_step_free(s);
gsl_odeiv2_evolve_free(e);

// save the times as a separate array
REAL *restrict times = malloc(nsteps * sizeof(REAL));
for (i = 0; i < nsteps; i++){
  times[i] = dynamics_RK[IDX(i,TIME)];
}

// find the coarse-fine separation index
REAL t_desired;
if (stop == OMEGA){
  t_desired = times[nsteps-1] - commondata->t_stepback - 50.;
} else{
  t_desired = times[nsteps-1] - commondata->t_stepback;
}
size_t coarse_fine_separation_idx = gsl_interp_bsearch(times,t_desired,0,nsteps-1);
REAL *times_fine_prelim = NULL;
REAL *dynamics_fine_prelim = NULL;
size_t nsteps_fine_prelim = 0;
// handle case where t_desired is less than the first timestep (only happens when PA integration is performed).
if (coarse_fine_separation_idx == 0){
  nsteps_fine_prelim = nsteps;
  dynamics_fine_prelim = (REAL *)malloc(NUMVARS * nsteps_fine_prelim * sizeof(REAL));
  if (dynamics_fine_prelim == NULL){
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for dynamics_fine_prelim\\n");
    exit(1);
  }
  times_fine_prelim = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (times_fine_prelim == NULL){
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for times_fine_prelim\\n");
    exit(1);
  }
  memcpy(dynamics_fine_prelim,dynamics_RK,nsteps_fine_prelim * NUMVARS * sizeof(REAL));
  memcpy(times_fine_prelim,times,nsteps_fine_prelim * sizeof(REAL));
} else{
  commondata->dynamics_low = (REAL *)malloc(NUMVARS * coarse_fine_separation_idx * sizeof(REAL));
  if (commondata->dynamics_low == NULL){
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for commondata->dynamics_low\\n");
    exit(1);
  }
  memcpy(commondata->dynamics_low,dynamics_RK,coarse_fine_separation_idx * NUMVARS * sizeof(REAL));
  commondata->nsteps_low = coarse_fine_separation_idx;

  nsteps_fine_prelim = nsteps - coarse_fine_separation_idx;
  dynamics_fine_prelim = (REAL *)malloc(NUMVARS * nsteps_fine_prelim * sizeof(REAL));
  if (dynamics_fine_prelim == NULL){
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for dynamics_fine_prelim\\n");
    exit(1);
  }
  times_fine_prelim = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (times_fine_prelim == NULL){
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for times_fine_prelim\\n");
    exit(1);
  }
  memcpy(dynamics_fine_prelim,dynamics_RK + coarse_fine_separation_idx * NUMVARS,nsteps_fine_prelim * NUMVARS * sizeof(REAL));
  memcpy(times_fine_prelim,times + coarse_fine_separation_idx,nsteps_fine_prelim * sizeof(REAL));
}
commondata->dynamics_raw = malloc(NUMVARS * nsteps * sizeof(REAL));
if (commondata->dynamics_raw == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for commondata->dynamics_raw\\n");
  exit(1);
}
commondata->nsteps_raw = nsteps;
memcpy(commondata->dynamics_raw,dynamics_RK,nsteps * NUMVARS * sizeof(REAL));
free(dynamics_RK);

// t_peak = dynamics[-1] if there is no peak

REAL t_peak = times_fine_prelim[nsteps_fine_prelim - 1];
// perform iterative refinement to find the true peak of the dynamics
if (stop == OMEGA) {
  REAL *t_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (t_values == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for t_values\\n");
    exit(1);
  }
  REAL *omega_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (omega_values == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for omega_values\\n");
    exit(1);
  }
for (i = 0; i < nsteps_fine_prelim; i++) {
    t_values[i] = dynamics_fine_prelim[IDX(i, TIME)];
    omega_values[i] = dynamics_fine_prelim[IDX(i, OMEGA)];
  }
  REAL dt = times_fine_prelim[nsteps_fine_prelim - 1] - times_fine_prelim[nsteps_fine_prelim - 2];
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  if (acc == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), gsl_interp_accel_alloc() failed\\n");
    exit(1);
  }
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  if (spline == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), gsl_spline_alloc() failed\\n");
    exit(1);
  }
  gsl_spline_init(spline, t_values, omega_values, nsteps_fine_prelim);
  spline_data sdata = {spline, acc};
  t_peak = SEOBNRv5_aligned_spin_iterative_refinement(&sdata, times_fine_prelim[0], times_fine_prelim[nsteps_fine_prelim - 1], 3, dt, false);
  free(t_values);
  free(omega_values);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
if (stop == PRSTAR) {
  REAL *t_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (t_values == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for t_values\\n");
    exit(1);
  }
  REAL *prstar_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  if (prstar_values == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), malloc() failed for prstar_values\\n");
    exit(1);
  }
  for (i = 0; i < nsteps_fine_prelim; i++) {
    t_values[i] = dynamics_fine_prelim[IDX(i, TIME)];
    prstar_values[i] = dynamics_fine_prelim[IDX(i, PRSTAR)];
  }
  REAL dt = times_fine_prelim[nsteps_fine_prelim - 1] - times_fine_prelim[nsteps_fine_prelim - 2];
  gsl_interp_accel *acc = gsl_interp_accel_alloc();
  if (acc == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_ode_integration(), gsl_interp_accel_alloc() failed\\n");
    exit(1);
  }
  gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  if (spline == NULL) {
    fprintf(stderr, "Error: in SEOBNRSEOBNRv5_aligned_spin_ode_integration(), gsl_spline_alloc() failed\\n");
    exit(1);
  }
  gsl_spline_init(spline, t_values, prstar_values, nsteps_fine_prelim);
  spline_data sdata = {spline, acc};
  t_peak = SEOBNRv5_aligned_spin_iterative_refinement(&sdata, times_fine_prelim[0], times_fine_prelim[nsteps_fine_prelim - 1], 3, dt, true);
  free(t_values);
  free(prstar_values);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}    
// interpolate the dynamics

SEOBNRv5_aligned_spin_interpolate_dynamics(commondata,dynamics_fine_prelim,nsteps_fine_prelim,t_peak,stop);
free(dynamics_fine_prelim);
free(times_fine_prelim);
free(times);;
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

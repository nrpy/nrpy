"""
Set up C function library for SEOBNR inspiral integrations.

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
    desc = """Evaluate SEOBNRv5 Hamiltonian and circular derivatives."""
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
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_SEOBNRv5_aligned_spin_argrelmin() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction that performs scipy.argrelmin with order = 3.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """C function to perform scipy.argrelmin for order = 3."""
    cfunc_type = "size_t"
    name = "SEOBNRv5_aligned_spin_argrelmin"
    params = "REAL *restrict arr , size_t nsteps_arr"
    body = """
size_t order = 3;
size_t minima_count_or_idx = 0;
// Loop through the array with bounds of `order` on each side
for (size_t i = order; i < nsteps_arr - order; i++) {
  int isMinimum = 1;

  // Check `order` elements to the left and right
  for (int j = 1; j <= order; j++) {
    if (arr[i] >= arr[i - j] || arr[i] >= arr[i + j]) {
      isMinimum = 0;
      break;
    }
  }
  if (isMinimum) {
    minima_count_or_idx = i;
    break;
  }
}
return minima_count_or_idx;
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


def register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the peak in frequency or momentum in SEOBNRv5.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate the peak in frequency or momentum in SEOBNRv5"""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_iterative_refinement"
    params = "gsl_spline *restrict spline, gsl_interp_accel *restrict acc, const REAL left, const REAL right, const int stop"
    body = """
REAL left_refined = left;
REAL right_refined = right;
REAL dt = 0.1 , result;
REAL *restrict abs_F_dots , *restrict Ts;
size_t nsteps_Ts , i;
size_t minimaIDX;
for(int iter = 0; iter < 2; iter++) {
  dt = dt * 0.1;
  nsteps_Ts = (size_t) ((right - left) / dt );
  Ts = (REAL *)malloc(nsteps_Ts * sizeof(REAL));
  abs_F_dots = (REAL *)malloc(nsteps_Ts * sizeof(REAL));
  for (i = 0 ; i < nsteps_Ts; i++){
    Ts[i] = left + i * dt;
    abs_F_dots[i] = fabs(gsl_spline_eval_deriv(spline,Ts[i],acc));
  }
  minimaIDX = SEOBNRv5_aligned_spin_argrelmin(abs_F_dots,nsteps_Ts);
  result = Ts[minimaIDX];
  free(Ts);
  free(abs_F_dots);
  if (minimaIDX > 0){
    left_refined = MAX(left_refined,result - 10.);
    right_refined = MIN(right_refined,result + 10.);
  }
  else{
    if (stop == PRSTAR){
      return right_refined;
    }
    else{
      result = (left_refined + right_refined) * 0.5;
      return result;
    }
  } 
}
return result;
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


def register_CFunction_SEOBNRv5_aligned_spin_intepolate_dynamics() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for interpolating the SEOBNRv5 dynamics.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate the peak in frequency or momentum in SEOBNRv5"""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_interpolate_dynamics"
    params = "commondata_struct *restrict commondata, REAL *restrict dynamics_fine_prelim, const size_t nsteps_fine_prelim, const REAL t_peak, const int stop"
    body = """
int i;
// Intepolate the high sampled dynamics for NQCs.
REAL time_start = dynamics_fine_prelim[IDX(0,TIME)];
REAL time_end = dynamics_fine_prelim[IDX(nsteps_fine_prelim-1,TIME)];
if (stop != 0){
  time_start = MAX(t_peak - commondata->t_stepback,dynamics_fine_prelim[IDX(0,TIME)]);
  time_end = MIN(t_peak , time_end);
}

REAL *restrict ts = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict rs = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict phis = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict prs = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
REAL *restrict pphis = (REAL *)malloc(nsteps_fine_prelim*sizeof(REAL));
for (i = 0; i < nsteps_fine_prelim; i++) {
  ts[i] = dynamics_fine_prelim[IDX(i,TIME)];
  rs[i] = dynamics_fine_prelim[IDX(i,R)];
  phis[i] = dynamics_fine_prelim[IDX(i,PHI)];
  prs[i] = dynamics_fine_prelim[IDX(i,PRSTAR)];
  pphis[i] = dynamics_fine_prelim[IDX(i,PPHI)];
}
gsl_interp_accel *restrict r_acc = gsl_interp_accel_alloc();
gsl_spline *restrict r_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(r_spline, ts, rs, nsteps_fine_prelim);
gsl_interp_accel *restrict phi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict phi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(phi_spline, ts, phis, nsteps_fine_prelim);
gsl_interp_accel *restrict pr_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pr_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(pr_spline, ts, prs, nsteps_fine_prelim);
gsl_interp_accel *restrict pphi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pphi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
gsl_spline_init(pphi_spline, ts, pphis, nsteps_fine_prelim);

const REAL dt = 0.1;
commondata->nsteps_fine = (size_t)(time_end - time_start) / dt;
commondata->dynamics_fine = (REAL *)malloc(NUMVARS * commondata->nsteps_fine * sizeof(REAL));
REAL t;
for (i = 0; i < commondata->nsteps_fine; i++) {
  t = time_start + i * dt;
  commondata->dynamics_fine[IDX(i , TIME)] = t;
  commondata->dynamics_fine[IDX(i , R)] = gsl_spline_eval(r_spline, t, r_acc);
  commondata->dynamics_fine[IDX(i , PHI)] = gsl_spline_eval(phi_spline, t, phi_acc);
  commondata->dynamics_fine[IDX(i , PRSTAR)] = gsl_spline_eval(pr_spline, t, pr_acc);
  commondata->dynamics_fine[IDX(i , PPHI)] = gsl_spline_eval(pphi_spline, t, pphi_acc);
  commondata->r = commondata->dynamics_fine[IDX(i , R)];
  commondata->phi = commondata->dynamics_fine[IDX(i , PHI)];
  commondata->prstar = commondata->dynamics_fine[IDX(i , PRSTAR)];
  commondata->pphi = commondata->dynamics_fine[IDX(i , PPHI)];
  SEOBNRv5_aligned_spin_augments(commondata);
  commondata->dynamics_fine[IDX(i , H)] = commondata->Hreal;
  commondata->dynamics_fine[IDX(i , OMEGA)] = commondata->dHreal_dpphi;
  commondata->dynamics_fine[IDX(i , OMEGA_CIRC)] = commondata->Omega_circ;
}

gsl_spline_free(r_spline);
gsl_interp_accel_free(r_acc);
gsl_spline_free(phi_spline);
gsl_interp_accel_free(phi_acc);
gsl_spline_free(pr_spline);
gsl_interp_accel_free(pr_acc);
gsl_spline_free(pphi_spline);
gsl_interp_accel_free(pphi_acc);

free(ts);
free(rs);
free(phis);
free(prs);
free(pphis);
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


def register_CFunction_SEOBNRv5_aligned_spin_ode_integration(
    perform_iterative_refinement: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for integrating the SEOBNRv5 equations of motion using GSL.

    :param perform_iterative_refinement: True/False if using iterative refinement to find the frequency peak.
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
    desc = """Integrate the SEOBNRv5 equations of motion."""
    cfunc_type = "int"
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
}
else{
  t_desired = times[nsteps-1] - commondata->t_stepback;
}
size_t coarse_fine_separation_idx = gsl_interp_bsearch(times,t_desired,0,nsteps-1);
commondata->dynamics_low = (REAL *)malloc(NUMVARS * coarse_fine_separation_idx * sizeof(REAL));
commondata->nsteps_low = coarse_fine_separation_idx;
for (i = 0; i < commondata->nsteps_low; i++){
  commondata->dynamics_low[IDX(i,TIME)] = dynamics_RK[IDX(i,TIME)];
  commondata->dynamics_low[IDX(i,R)] = dynamics_RK[IDX(i,R)];
  commondata->dynamics_low[IDX(i,PHI)] = dynamics_RK[IDX(i,PHI)];
  commondata->dynamics_low[IDX(i,PRSTAR)] = dynamics_RK[IDX(i,PRSTAR)];
  commondata->dynamics_low[IDX(i,PPHI)] = dynamics_RK[IDX(i,PPHI)];
  commondata->dynamics_low[IDX(i,H)] = dynamics_RK[IDX(i,H)];
  commondata->dynamics_low[IDX(i,OMEGA)] = dynamics_RK[IDX(i,OMEGA)];
  commondata->dynamics_low[IDX(i,OMEGA_CIRC)] = dynamics_RK[IDX(i,OMEGA_CIRC)];
}

size_t nsteps_fine_prelim = nsteps - coarse_fine_separation_idx;
REAL *restrict dynamics_fine_prelim = (REAL *)malloc(NUMVARS * nsteps_fine_prelim * sizeof(REAL));
REAL *restrict times_fine_prelim = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
for (i = 0; i < nsteps_fine_prelim; i++){
  dynamics_fine_prelim[IDX(i,TIME)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,TIME)];
  times_fine_prelim[i] = dynamics_fine_prelim[IDX(i,TIME)];
  dynamics_fine_prelim[IDX(i,R)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,R)];
  dynamics_fine_prelim[IDX(i,PHI)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,PHI)];
  dynamics_fine_prelim[IDX(i,PRSTAR)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,PRSTAR)];
  dynamics_fine_prelim[IDX(i,PPHI)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,PPHI)];
  dynamics_fine_prelim[IDX(i,H)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,H)];
  dynamics_fine_prelim[IDX(i,OMEGA)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,OMEGA)];
  dynamics_fine_prelim[IDX(i,OMEGA_CIRC)] = dynamics_RK[IDX(i + coarse_fine_separation_idx,OMEGA_CIRC)];
}

free(dynamics_RK);

// t_peak = dynamics[-1] if there is no peak

REAL t_peak = times_fine_prelim[nsteps_fine_prelim - 1];
"""
    if perform_iterative_refinement:
        body += """
// perform iterative refinement to find the true peak of the dynamics
if (stop != 0){
  REAL *restrict fpeak_fine = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  for (i = 0; i < nsteps_fine_prelim;i++){
    fpeak_fine[i] = dynamics_fine_prelim[IDX(i,stop)];
  }
  gsl_interp_accel *acc
      = gsl_interp_accel_alloc ();
  gsl_spline *spline
      = gsl_spline_alloc (gsl_interp_cspline, nsteps_fine_prelim);
  gsl_spline_init(spline,times_fine_prelim,fpeak_fine,nsteps_fine_prelim);
  REAL right = times_fine_prelim[nsteps_fine_prelim - 1];
  REAL left = times_fine_prelim[0];
  if (stop == PRSTAR){
    left = times_fine_prelim[nsteps_fine_prelim - 1] - 10.;
  }
  t_peak = SEOBNRv5_aligned_spin_iterative_refinement(spline,acc,left,right,stop);
  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  free(fpeak_fine);
}
"""
    body += """
// interpolate the dynamics

status = SEOBNRv5_aligned_spin_interpolate_dynamics(commondata,dynamics_fine_prelim,nsteps_fine_prelim,t_peak,stop);
free(dynamics_fine_prelim);
free(times_fine_prelim);
free(times);
return GSL_SUCCESS;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        prefunc=prefunc,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

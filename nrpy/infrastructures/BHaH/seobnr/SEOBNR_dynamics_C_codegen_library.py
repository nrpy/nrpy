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
        [H.Hreal, H.dHreal_dpphi_circ],
        ["commondata->Hreal", "commondata->Omega_circ"],
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


def register_CFunction_SEOBNRv5_aligned_spin_find_peak() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction taking a 1-D interpolant of Omega or pr and returning the location of its peak.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Find the peak of Omega or pr.
Used to refine the location of the maximum in frequency or momentum."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_find_peak"
    params = "gsl_spline *restrict spline, gsl_interp_accel *restrict acc, const REAL left, const REAL right, const REAL dt, size_t result_idx, REAL result"
    body = """
size_t n = (size_t) (right - left)/dt;
REAL x;
REAL dx[n];
for (size_t i = 0; i < n; i++){
  x = left + i*dt;
  dx[i] = gsl_spline_eval_deriv(spline, x, acc);
}
result_idx = gsl_interp_bsearch(dx, 0.0, 0, n);
result = dx[result_idx];
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
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_iterative_refinement"
    params = "gsl_spline *restrict spline, gsl_interp_accel *restrict acc, const int stop, const REAL left, const REAL right, REAL t_ext"
    body = """
REAL left_refined = left;
REAL right_refined = right;
REAL dt = 0.1;
int status;
int iter = 2;
size_t result_idx;
while (iter > 0) {
  dt = dt * 0.1;
  status = SEOBNRv5_aligned_spin_find_peak(spline, acc, left_refined, right_refined, dt, result_idx, t_ext);
  if (result_idx == 0 || result_idx == (size_t)(right_refined - left_refined) / dt) {
    // no actual extremum detected
    t_ext = (0.5 * (left_refined + right_refined)) * (1 - stop % 2) + (stop % 2) * (right_refined);
    break;
  }
  left_refined = 0.5 * (t_ext - 10. + left + fabs(t_ext - 10. - left));
  right_refined = 0.5 * (t_ext + 10. + right - fabs(t_ext + 10 - right));
  iter--;
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
    params = "commondata_struct *restrict commondata, REAL *restrict times, REAL *restrict dynamics_RK, const size_t nsteps, const REAL t_max, const REAL t_stepback"
    body = """
int i;

// Store the low sampled dynamics.
commondata->nsteps_low = gsl_interp_bsearch(times, t_stepback, 0, nsteps - 1);
commondata->dynamics_low = calloc(8 * commondata->nsteps_low, sizeof(REAL));
for (i = 0; i < commondata->nsteps_low; i++) {
  commondata->dynamics_low[8 * i + 0] = dynamics_RK[8 * i + 0];
  commondata->dynamics_low[8 * i + 1] = dynamics_RK[8 * i + 1];
  commondata->dynamics_low[8 * i + 2] = dynamics_RK[8 * i + 2];
  commondata->dynamics_low[8 * i + 3] = dynamics_RK[8 * i + 3];
  commondata->dynamics_low[8 * i + 4] = dynamics_RK[8 * i + 4];
  commondata->dynamics_low[8 * i + 5] = dynamics_RK[8 * i + 5];
  commondata->dynamics_low[8 * i + 6] = dynamics_RK[8 * i + 6];
  commondata->dynamics_low[8 * i + 7] = dynamics_RK[8 * i + 7];
}

// Intepolate the high sampled dynamics for NQCs.
const REAL dt = 0.1;
commondata->nsteps_fine = (size_t)(t_max - t_stepback) / dt;
size_t len_dynamics_fine = nsteps - commondata->nsteps_low;
REAL times_fine[commondata->nsteps_fine];
REAL ts[len_dynamics_fine], rs[len_dynamics_fine], phis[len_dynamics_fine], prs[len_dynamics_fine], pphis[len_dynamics_fine], Hs[len_dynamics_fine],
    Omegas[len_dynamics_fine], Omega_circs[len_dynamics_fine];
for (i = 0; i < commondata->nsteps_fine; i++) {
  times_fine[i] = t_stepback + i * dt;
}
for (i = 0; i < len_dynamics_fine; i++) {
  ts[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 0];
  rs[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 1];
  phis[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 2];
  prs[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 3];
  pphis[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 4];
  Hs[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 5];
  Omegas[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 6];
  Omega_circs[i] = dynamics_RK[8 * (i + commondata->nsteps_low) + 7];
}
gsl_interp_accel *restrict r_acc = gsl_interp_accel_alloc();
gsl_spline *restrict r_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(r_spline, ts, rs, len_dynamics_fine);
gsl_interp_accel *restrict phi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict phi_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(phi_spline, ts, phis, len_dynamics_fine);
gsl_interp_accel *restrict pr_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pr_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(pr_spline, ts, prs, len_dynamics_fine);
gsl_interp_accel *restrict pphi_acc = gsl_interp_accel_alloc();
gsl_spline *restrict pphi_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(pphi_spline, ts, pphis, len_dynamics_fine);
gsl_interp_accel *restrict H_acc = gsl_interp_accel_alloc();
gsl_spline *restrict H_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(H_spline, ts, Hs, len_dynamics_fine);
gsl_interp_accel *restrict Omega_acc = gsl_interp_accel_alloc();
gsl_spline *restrict Omega_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(Omega_spline, ts, Omegas, len_dynamics_fine);
gsl_interp_accel *restrict Omega_circ_acc = gsl_interp_accel_alloc();
gsl_spline *restrict Omega_circ_spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
gsl_spline_init(Omega_circ_spline, ts, Omega_circs, len_dynamics_fine);

commondata->dynamics_fine = (REAL *)calloc(8 * commondata->nsteps_fine, sizeof(REAL));
for (i = 0; i < commondata->nsteps_fine; i++) {
  commondata->dynamics_fine[8 * i + 0] = times_fine[i];
  commondata->dynamics_fine[8 * i + 1] = gsl_spline_eval(r_spline, times_fine[i], r_acc);
  commondata->dynamics_fine[8 * i + 2] = gsl_spline_eval(phi_spline, times_fine[i], phi_acc);
  commondata->dynamics_fine[8 * i + 3] = gsl_spline_eval(pr_spline, times_fine[i], pr_acc);
  commondata->dynamics_fine[8 * i + 4] = gsl_spline_eval(pphi_spline, times_fine[i], pphi_acc);
  commondata->dynamics_fine[8 * i + 5] = gsl_spline_eval(H_spline, times_fine[i], H_acc);
  commondata->dynamics_fine[8 * i + 6] = gsl_spline_eval(Omega_spline, times_fine[i], Omega_acc);
  commondata->dynamics_fine[8 * i + 7] = gsl_spline_eval(Omega_circ_spline, times_fine[i], Omega_circ_acc);
}

gsl_spline_free(r_spline);
gsl_interp_accel_free(r_acc);
gsl_spline_free(phi_spline);
gsl_interp_accel_free(phi_acc);
gsl_spline_free(pr_spline);
gsl_interp_accel_free(pr_acc);
gsl_spline_free(pphi_spline);
gsl_interp_accel_free(pphi_acc);
gsl_spline_free(H_spline);
gsl_interp_accel_free(H_acc);
gsl_spline_free(Omega_spline);
gsl_interp_accel_free(Omega_acc);
gsl_spline_free(Omega_circ_spline);
gsl_interp_accel_free(Omega_circ_acc);

// Populate the combined dynamics

commondata->nsteps_combined = commondata->nsteps_fine + commondata->nsteps_low;
commondata->dynamics_combined = calloc(8 * commondata->nsteps_combined, sizeof(REAL));
for (i = 0; i < 8 * commondata->nsteps_low; i++) {
  commondata->dynamics_combined[i] = commondata->dynamics_low[i];
}
for (i = 0; i < 8 * (commondata->nsteps_fine); i++) {
  commondata->dynamics_combined[8 * commondata->nsteps_low + i] = commondata->dynamics_fine[i];
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
    desc = """Integrate the SEOBNRv5 equations of motion."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_ode_integration"
    params = "commondata_struct *restrict commondata"
    body = """
const gsl_odeiv2_step_type *restrict T = gsl_odeiv2_step_rk8pd;

gsl_odeiv2_step *restrict s = gsl_odeiv2_step_alloc(T, 4);
gsl_odeiv2_control *restrict c = gsl_odeiv2_control_standard_new(1e-12, 1e-11, 1.0, 1.0);
gsl_odeiv2_system sys = {SEOBNRv5_aligned_spin_right_hand_sides, NULL, 4, commondata};

REAL t = 0.0;
REAL t_new;
const REAL t1 = 2e9;
const int dyn_size = 8; // t,r,phi,prstar,pphi,Hreal, Omega, Omega_circ
REAL y[4], yerr[4], dydt_in[4], dydt_out[4];
int status = 0;
int i;
int stop = 0;
y[0] = commondata->r;
y[1] = commondata->phi;
y[2] = commondata->prstar;
y[3] = commondata->pphi;
status = SEOBNRv5_aligned_spin_right_hand_sides(t, y, dydt_in, commondata);
SEOBNRv5_aligned_spin_augments(commondata);
REAL h = 2.0 * M_PI / dydt_in[1] / 5.0;
size_t bufferlength = (size_t)(t1 / h); // runs up to 0.01x maximum time (we should not ideally run that long)
REAL *restrict dynamics_RK = (REAL *)calloc(bufferlength * (dyn_size), sizeof(REAL));
size_t nsteps = 0;

// store
dynamics_RK[dyn_size * nsteps] = t;
for (i = 0; i < 4; i++) {
  dynamics_RK[dyn_size * nsteps + 1 + i] = y[i];
}
dynamics_RK[dyn_size * nsteps + 5] = commondata->Hreal;
dynamics_RK[dyn_size * nsteps + 6] = dydt_in[1];
dynamics_RK[dyn_size * nsteps + 7] = commondata->Omega_circ;
nsteps++;

while (stop == 0) {
  // integrate
  status = gsl_odeiv2_step_apply(s, t, h, y, yerr, dydt_in, dydt_out, &sys);
  if (status != GSL_SUCCESS) {
    printf("Error in step_apply!\\n");
    return status;
  }
  status = gsl_odeiv2_control_hadjust(c, s, y, yerr, dydt_out, &h);
  t_new = t + h;
  commondata->r = y[0];
  commondata->phi = y[1];
  commondata->prstar = y[2];
  commondata->pphi = y[3];
  SEOBNRv5_aligned_spin_augments(commondata);

  // buffercheck
  if (nsteps >= bufferlength) {
    bufferlength = 2 * bufferlength;
    dynamics_RK = (REAL *)realloc(dynamics_RK, bufferlength * (dyn_size) * sizeof(REAL));
  }

  // update
  t = t_new;
  memcpy(dydt_in, dydt_out, 4 * sizeof(REAL));
  // store
  dynamics_RK[dyn_size * nsteps] = t;
  for (i = 0; i < 4; i++) {
    dynamics_RK[dyn_size * nsteps + 1 + i] = y[i];
  }
  dynamics_RK[dyn_size * nsteps + 5] = commondata->Hreal;
  dynamics_RK[dyn_size * nsteps + 6] = dydt_out[1];
  dynamics_RK[dyn_size * nsteps + 7] = commondata->Omega_circ;
  nsteps++;

  // stopcheck
  if (commondata->r < 6.0) {
    // decrease in frequency: Omega peak stop = index of omega
    if (dydt_out[1] < dydt_in[1]) {
      stop = 6;
      break;
    }
    // outspiral dPR: PR peak stop = index of pr
    if (dydt_out[2] > 0.0) {
      stop = 3;
      break;
    }
    // outspiral dR
    if (dydt_out[0] > 0.0) {
      break;
    }
    // unphysical frequency
    if (commondata->r < 3.0 && commondata->Omega_circ > 1.0) {
      break;
    }
  }
}

// free up gsl ode solver
gsl_odeiv2_control_free(c);
gsl_odeiv2_step_free(s);

// High sampling.
// Get an estimate of the stepback time.

REAL times[nsteps];
for (i = 0; i < nsteps; i++) {
  times[i] = dynamics_RK[dyn_size * i + 0];
}
REAL t_stepback = times[nsteps - 1] - commondata->t_stepback;
size_t idx_stepback = gsl_interp_bsearch(times, t_stepback, 0, nsteps - 1);

// Get the time of peak Omega or pr if peak was detected.

REAL t_max = times[nsteps - 1];
if (stop != 0) {
  size_t len_dynamics_fine = nsteps - idx_stepback;
  REAL *restrict fpeak_fine = (REAL *)calloc(len_dynamics_fine, sizeof(REAL));
  REAL *restrict times_fine = (REAL *)calloc(len_dynamics_fine, sizeof(REAL));
  for (i = 0; i < len_dynamics_fine; i++) {
    times_fine[i] = times[i + idx_stepback];
    fpeak_fine[i] = dynamics_RK[8 * (i + idx_stepback) + stop];
  }
  gsl_interp_accel *restrict acc = gsl_interp_accel_alloc();
  REAL t_left = times_fine[0] * (1 - stop % 2) + (stop % 2) * (times[nsteps - 1] - 10); // t_stepback if stop == 6 and t_end - 10M if stop == 3
  gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, len_dynamics_fine);
  gsl_spline_init(spline, times_fine, fpeak_fine, len_dynamics_fine);
  status = SEOBNRv5_aligned_spin_iterative_refinement(spline, acc, stop, t_left, dynamics_RK[8 * (nsteps - 1) + 0], t_max);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}
t_stepback = t_max - commondata->t_stepback;
status = SEOBNRv5_aligned_spin_interpolate_dynamics(commondata, times, dynamics_RK, nsteps, t_max, t_stepback);
// free the dynamics_RK pointer as it is no longer needed
free(dynamics_RK);
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

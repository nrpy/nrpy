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
import nrpy.equations.seobnr.BOB_aligned_spin_waveform_quantities as BOB
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_NQC_corrections() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the Non Quasi-Circular (NQC) corrections for the SEOBNRv5 waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate SEOBNRv5 NQC corrections."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_NQC_corrections"
    params = "commondata_struct *restrict commondata"
    bob_wf = BOB.BOB_aligned_spin_waveform_quantities()

    body = """
REAL *restrict times = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict Q1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict Q2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict Q3 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict P1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict P2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict r = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict Omega = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict hamp = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict phase = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL *restrict phase_unwrapped = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
REAL radius, omega, prstar, hplus, hcross; 
size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  hplus = commondata->waveform_fine[IDX_WF(i,HPLUS)];
  hcross = commondata->waveform_fine[IDX_WF(i,HCROSS)];
  hamp[i] = sqrt(hplus * hplus + hcross * hcross);
  phase[i] = atan2(-hcross,hplus);
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  Q1[i] = hamp[i] * prstar * prstar / (r[i] * r[i] * Omega[i] * Omega[i]);
  Q2[i] = Q1[i] / r[i];
  Q3[i] = Q2[i] / sqrt(r[i]);
  P1[i] = -prstar / r[i] /Omega[i];
  P2[i] = -P1[i] * prstar * prstar;
}
REAL diff;
int wrap = 0;
phase_unwrapped[0] = phase[0];
for (i = 1; i < commondata->nsteps_fine; i++){
  diff = phase[i] - phase[i-1];
  wrap += (diff < -M_PI) - (diff > M_PI);
  phase_unwrapped[i] = phase[i] + wrap * 2 * M_PI;
}

// Find t_ISCO:

if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
}
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);  
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
  for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;

    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  }  
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];
  
  gsl_interp_accel_free(acc_r);
  gsl_spline_free(spline_r);
  free(t_zoom);
  free(minus_r_zoom);
}

REAL t_peak = commondata->t_ISCO - commondata->Delta_t;
if (t_peak > times[commondata->nsteps_fine - 1]){
  t_peak = times[commondata->nsteps_fine - 1];
}
size_t peak_idx = gsl_interp_bsearch(times, t_peak, 0, commondata->nsteps_fine);

size_t N = 5;
size_t left = (peak_idx - N + abs(peak_idx - N))/2;
size_t right = (peak_idx + N + commondata->nsteps_fine - abs(peak_idx + N - commondata->nsteps_fine))/2;
REAL Q_cropped1[right - left], Q_cropped2[right - left], Q_cropped3[right - left], P_cropped1[right - left], P_cropped2[right - left];
REAL t_cropped[right-left], phase_cropped[right - left], amp_cropped[right-left];
for (i = left; i < right; i++){
  t_cropped[i - left] = times[i];
  phase_cropped[i - left] = phase_unwrapped[i];
  amp_cropped[i - left] = hamp[i];
  Q_cropped1[i - left] = Q1[i];
  Q_cropped2[i - left] = Q2[i];
  Q_cropped3[i - left] = Q3[i];
  P_cropped1[i - left] = P1[i];
  P_cropped2[i - left] = P2[i];
}
gsl_interp_accel *restrict acc = gsl_interp_accel_alloc();
gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, right - left);

gsl_matrix *restrict Q = gsl_matrix_alloc (3, 3);
gsl_matrix *restrict P = gsl_matrix_alloc (2, 2);

gsl_spline_init(spline,t_cropped, Q_cropped1,right-left);
gsl_matrix_set(Q,0,0,gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(Q,1,0,gsl_spline_eval_deriv(spline, t_peak, acc));
gsl_matrix_set(Q,2,0,gsl_spline_eval_deriv2(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, Q_cropped2,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(Q,0,1,gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(Q,1,1,gsl_spline_eval_deriv(spline, t_peak, acc));
gsl_matrix_set(Q,2,1,gsl_spline_eval_deriv2(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, Q_cropped3,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(Q,0,2,gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(Q,1,2,gsl_spline_eval_deriv(spline, t_peak, acc));
gsl_matrix_set(Q,2,2,gsl_spline_eval_deriv2(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, P_cropped1,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(P,0,0,-gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(P,1,0,-gsl_spline_eval_deriv(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, P_cropped2,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(P,0,1,-gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(P,1,1,-gsl_spline_eval_deriv(spline, t_peak, acc));

gsl_vector *restrict A = gsl_vector_alloc(3);
gsl_vector *restrict O = gsl_vector_alloc(2);


gsl_spline_init(spline,t_cropped, amp_cropped,right-left);
gsl_interp_accel_reset(acc);
const REAL amp_insp = gsl_spline_eval(spline, t_peak, acc);
const REAL ampdot_insp = gsl_spline_eval_deriv(spline, t_peak, acc);
const REAL ampddot_insp = gsl_spline_eval_deriv2(spline, t_peak, acc);

gsl_spline_init(spline,t_cropped, phase_cropped,right-left);
gsl_interp_accel_reset(acc);
REAL omega_insp = gsl_spline_eval_deriv(spline, t_peak, acc);
REAL omegadot_insp =  gsl_spline_eval_deriv2(spline, t_peak, acc);

gsl_spline_free(spline);
gsl_interp_accel_free(acc);

if (omega_insp * omegadot_insp > 0.0){
  omega_insp = fabs(omega_insp);
  omegadot_insp = fabs(omegadot_insp);
}
else{
  omega_insp = fabs(omega_insp);
  omegadot_insp = -fabs(omegadot_insp);
}

const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
const REAL t = t_peak;
const REAL t_0 = t_peak;
"""
    body += ccg.c_codegen(
        [
            bob_wf.h_t_attach,
            bob_wf.hdot_t_attach,
            bob_wf.hddot_t_attach,
            bob_wf.w_t_attach,
            bob_wf.wdot_t_attach,
        ],
        [
            "const REAL h_t_attach",
            "const REAL hdot_t_attach",
            "const REAL hddot_t_attach",
            "const REAL w_t_attach",
            "const REAL wdot_t_attach",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
gsl_vector_set(A , 0 , h_t_attach - amp_insp);
gsl_vector_set(A , 1 , hdot_t_attach - ampdot_insp);
gsl_vector_set(A , 2 , hddot_t_attach - ampddot_insp);
gsl_vector_set(O , 0 , w_t_attach - omega_insp);
gsl_vector_set(O , 1 , wdot_t_attach - omegadot_insp);

gsl_vector *restrict a = gsl_vector_alloc (3);

int s;
gsl_permutation *restrict p_A = gsl_permutation_alloc (3);
gsl_linalg_LU_decomp(Q, p_A, &s);
gsl_linalg_LU_solve(Q, p_A, A, a);
gsl_permutation_free(p_A);

commondata->a_1_NQC = gsl_vector_get(a,0);
commondata->a_2_NQC = gsl_vector_get(a,1);
commondata->a_3_NQC = gsl_vector_get(a,2);

gsl_vector_free(a);
gsl_vector_free(A);
gsl_matrix_free(Q);

gsl_vector *restrict b = gsl_vector_alloc(2);
gsl_permutation * p_B = gsl_permutation_alloc(2);
gsl_linalg_LU_decomp(P, p_B, &s);
gsl_linalg_LU_solve (P, p_B, O, b);
gsl_permutation_free (p_B);

commondata->b_1_NQC = gsl_vector_get(b,0);
commondata->b_2_NQC = gsl_vector_get(b,1);

gsl_vector_free (b);
gsl_vector_free(O);
gsl_matrix_free(P);

free(times);
free(Q1);
free(Q2);
free(Q3);
free(P1);
free(P2);
free(r);
free(hamp);
free(phase);
free(phase_unwrapped);

// apply the nqc correction
commondata->waveform_inspiral = (REAL *)malloc(commondata->nsteps_inspiral*NUMMODES*sizeof(REAL));
REAL nqc_amp, nqc_phase, q1, q2, q3,p1, p2;
for (i = 0; i < commondata->nsteps_low; i++){
  prstar = commondata->dynamics_low[IDX(i,PRSTAR)];
  radius = commondata->dynamics_low[IDX(i,R)];
  omega = commondata->dynamics_low[IDX(i,OMEGA)];
  q1 = prstar * prstar / (radius * radius * omega * omega);
  q2 = q1 / radius;
  q3 = q2 / sqrt(radius);
  p1 = -prstar / radius /omega;
  p2 = -p1 * prstar * prstar;
  nqc_amp = 1 + commondata->a_1_NQC*q1 + commondata->a_2_NQC*q2 + commondata->a_3_NQC*q3;
  nqc_phase =  commondata->b_1_NQC*p1 + commondata->b_2_NQC*p2;
  commondata->waveform_inspiral[IDX_WF(i,TIME)] = commondata->dynamics_low[IDX(i,TIME)];
  commondata->waveform_inspiral[IDX_WF(i,HPLUS)] = nqc_amp * (commondata->waveform_low[IDX_WF(i,HPLUS)]*cos(nqc_phase) + commondata->waveform_low[IDX_WF(i,HCROSS)]*sin(nqc_phase));
  commondata->waveform_inspiral[IDX_WF(i,HCROSS)] = nqc_amp * (commondata->waveform_low[IDX_WF(i,HCROSS)]*cos(nqc_phase) - commondata->waveform_low[IDX_WF(i,HPLUS)]*sin(nqc_phase));
}
for (i = 0; i< commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  radius = commondata->dynamics_fine[IDX(i,R)];
  omega = commondata->dynamics_fine[IDX(i,OMEGA)];
  q1 = prstar * prstar / (radius * radius * omega * omega);
  q2 = q1 / radius;
  q3 = q2 / sqrt(radius);
  p1 = -prstar / radius /omega;
  p2 = -p1 * prstar * prstar;
  nqc_amp = 1 + commondata->a_1_NQC*q1 + commondata->a_2_NQC*q2 + commondata->a_3_NQC*q3;
  nqc_phase =  commondata->b_1_NQC*p1 + commondata->b_2_NQC*p2;
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,TIME)] = commondata->dynamics_fine[IDX(i,TIME)];
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,HPLUS)] = nqc_amp * (commondata->waveform_fine[IDX_WF(i,HPLUS)]*cos(nqc_phase) + commondata->waveform_fine[IDX_WF(i,HCROSS)]*sin(nqc_phase));
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,HCROSS)] = nqc_amp * (commondata->waveform_fine[IDX_WF(i,HCROSS)]*cos(nqc_phase) - commondata->waveform_fine[IDX_WF(i,HPLUS)]*sin(nqc_phase));
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


def register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the (2,2) IMR mode for the SEOBNRv5 waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Evaluate SEOBNRv5 (2,2) mode."""
    cfunc_type = "int"
    name = "SEOBNRv5_aligned_spin_IMR_waveform"
    params = "commondata_struct *restrict commondata"
    bob_wf = BOB.BOB_aligned_spin_waveform_quantities()
    body = """
int i;
REAL *restrict times_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_amp_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_phase_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict orbital_phase_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_phase_unwrapped_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_phase_unmodulated_old = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL hplus,hcross;
for (i = 0; i < commondata->nsteps_inspiral; i++){
  hplus = commondata->waveform_inspiral[IDX_WF(i,HPLUS)];
  hcross = commondata->waveform_inspiral[IDX_WF(i,HCROSS)];
  times_old[i] = commondata->waveform_inspiral[IDX_WF(i,TIME)];
  orbital_phase_old[i] = commondata->dynamics_inspiral[IDX(i,PHI)];
  h22_amp_old[i] = hplus*hplus + hcross*hcross;
  h22_phase_old[i] = atan2(-hcross,hplus);
}
int wrap = 0;
REAL diff;
h22_phase_unmodulated_old[0] = h22_phase_old[0] + 2 * orbital_phase_old[0];
h22_phase_unwrapped_old[0] = h22_phase_old[0];
for (i = 1; i < commondata->nsteps_inspiral; i++){
  diff = h22_phase_old[i] - h22_phase_old[i-1];
  wrap += (diff < -M_PI) - (diff > M_PI);
  h22_phase_unwrapped_old[i] = h22_phase_old[i] + wrap * 2 * M_PI;
  h22_phase_unmodulated_old[i] = h22_phase_unwrapped_old[i] + 2 * orbital_phase_old[i];
}
free(h22_phase_old);
gsl_interp_accel *restrict acc_amp = gsl_interp_accel_alloc();
gsl_spline *restrict spline_amp = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_inspiral);
gsl_spline_init(spline_amp,times_old,h22_amp_old, commondata->nsteps_inspiral);
gsl_interp_accel *restrict acc_phase = gsl_interp_accel_alloc();
gsl_spline *restrict spline_phase = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_inspiral);
gsl_spline_init(spline_phase,times_old,h22_phase_unmodulated_old, commondata->nsteps_inspiral);
gsl_interp_accel *restrict acc_phase_orb = gsl_interp_accel_alloc();
gsl_spline *restrict spline_phase_orb = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_inspiral);
gsl_spline_init(spline_phase_orb,times_old,orbital_phase_old, commondata->nsteps_inspiral);

const REAL dT = commondata->dt/(commondata->total_mass*4.925490947641266978197229498498379006e-6);
const size_t nsteps_insp_plunge = (size_t) commondata->dynamics_inspiral[IDX(commondata->nsteps_inspiral - 1,TIME)] / dT;
REAL *restrict times_new = malloc(nsteps_insp_plunge*sizeof(REAL));
REAL *restrict h22_amp_new = (REAL *)malloc(nsteps_insp_plunge*sizeof(REAL));
REAL *restrict h22_phase_new = (REAL *)malloc(nsteps_insp_plunge*sizeof(REAL));
for(i = 0; i < nsteps_insp_plunge; i++){
  times_new[i] = i * dT;
  h22_amp_new[i] = gsl_spline_eval(spline_amp,times_new[i],acc_amp);
  h22_phase_new[i] = gsl_spline_eval(spline_phase,times_new[i],acc_phase) - 2*gsl_spline_eval(spline_phase_orb,times_new[i],acc_phase_orb);
}

free(h22_amp_old);
free(times_old);
free(h22_phase_unwrapped_old);
free(h22_phase_unmodulated_old);
free(orbital_phase_old);
gsl_spline_free(spline_amp);
gsl_interp_accel_free(acc_amp);
gsl_spline_free(spline_phase);
gsl_interp_accel_free(acc_phase);
gsl_spline_free(spline_phase_orb);
gsl_interp_accel_free(acc_phase_orb);
REAL t_peak = commondata->t_ISCO - commondata->Delta_t;
size_t idx_attach;
if (t_peak > times_new[nsteps_insp_plunge - 1]){
  t_peak = times_new[nsteps_insp_plunge - 1];
  idx_attach = nsteps_insp_plunge - 2;
}
else{
  idx_attach = gsl_interp_bsearch(times_new, t_peak, 0, nsteps_insp_plunge);
  if (times_new[idx_attach] > t_peak){
    idx_attach--;
  }
}
const REAL t_attach = times_new[idx_attach];
const size_t ringdown_time = 15 * (size_t) (commondata->tau_qnm / dT);
REAL *restrict ringdown_waveform_amp_phase = (REAL *)malloc(NUMMODES*ringdown_time*sizeof(REAL));
REAL *restrict ringdown_phase_wrapped = (REAL *)malloc(ringdown_time*sizeof(REAL));

const REAL m1 = commondata->m1;
const REAL m2 = commondata->m2;
const REAL chi1 = commondata->chi1;
const REAL chi2 = commondata->chi2;
const REAL omega_qnm = commondata->omega_qnm;
const REAL tau_qnm = commondata->tau_qnm;
const REAL t_0 = t_peak;
REAL t, phase_before_sign;
REAL true_sign = copysign(1.0,h22_phase_new[idx_attach]);
const REAL phase_to_be_matched = h22_phase_new[idx_attach + 1];
for (i = 0; i < ringdown_time; i++){
  t = (i + 1) * dT + t_attach;
  ringdown_waveform_amp_phase[IDX_WF(i,TIME)] = t;
"""
    body += ccg.c_codegen(
        [
            bob_wf.h,
            bob_wf.phi,
        ],
        [
            "ringdown_waveform_amp_phase[IDX_WF(i,HAMP)]",
            "phase_before_sign",
        ],
        verbose=False,
        include_braces=False,
    )
    body += """
  ringdown_phase_wrapped[i] = true_sign * fabs(phase_before_sign);
  ringdown_phase_wrapped[i] += phase_to_be_matched - ringdown_phase_wrapped[0];
}
ringdown_waveform_amp_phase[IDX_WF(0,HPHASE)] = ringdown_phase_wrapped[0];
wrap = 0;
for (i = 1; i < ringdown_time; i++){
  diff = ringdown_phase_wrapped[i] - ringdown_phase_wrapped[i-1];
  wrap += (diff < -M_PI) - (diff > M_PI);
  ringdown_waveform_amp_phase[IDX_WF(0,HPHASE)] = ringdown_phase_wrapped[i] + wrap * 2 * M_PI;
}
free(ringdown_phase_wrapped);

commondata->nsteps_IMR = idx_attach + ringdown_time;
commondata->waveform_IMR = malloc(NUMMODES * commondata->nsteps_IMR*sizeof(REAL));
for (i = 0; i < idx_attach + 1; i++){
  commondata->waveform_IMR[IDX_WF(i,TIME)] = times_new[i];
  commondata->waveform_IMR[IDX_WF(i,HPLUS)] = h22_amp_new[i] * cos(h22_phase_new[i]);
  commondata->waveform_IMR[IDX_WF(i,HCROSS)] = h22_amp_new[i] * sin(-1. * h22_phase_new[i]);
}
free(h22_amp_new);
free(h22_phase_new);
free(times_new);

for(i = 0; i < ringdown_time; i++){
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_attach,TIME)] = ringdown_waveform_amp_phase[IDX_WF(i,TIME)];
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_attach,HPLUS)] = ringdown_waveform_amp_phase[IDX_WF(i,HAMP)] * cos(ringdown_waveform_amp_phase[IDX_WF(i,HPHASE)]);
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_attach,HCROSS)] = ringdown_waveform_amp_phase[IDX_WF(i,HAMP)] * sin(-1. * ringdown_waveform_amp_phase[IDX_WF(i,HPHASE)]);
}

free(ringdown_waveform_amp_phase);
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

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
REAL times[commondata->nsteps_fine], Q1[commondata->nsteps_fine], Q2[commondata->nsteps_fine], Q3[commondata->nsteps_fine], P1[commondata->nsteps_fine], P2[commondata->nsteps_fine];
REAL r[commondata->nsteps_fine], hamp[commondata->nsteps_fine], phase[commondata->nsteps_fine], phase_unwrapped[commondata->nsteps_fine];
REAL radius, prstar,Omega; 
size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  Omega = commondata->dynamics_fine[IDX(i,OMEGA)];
  hamp[i] = sqrt(commondata->waveform_fine[IDX_WF(i,HPLUS)]*commondata->waveform_fine[IDX_WF(i,HPLUS)] + commondata->waveform_fine[IDX_WF(i,HCROSS)]*commondata->waveform_fine[IDX_WF(i,HCROSS)] );
  phase[i] = atan2(-commondata->waveform_fine[IDX_WF(i,HCROSS)],commondata->waveform_fine[IDX_WF(i,HPLUS)]);
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  Q1[i] = hamp[i] * prstar * prstar / (r[i] * r[i] * Omega * Omega);
  Q2[i] = Q1[i] / r[i];
  Q3[i] = Q2[i] / sqrt(r[i]);
  P1[i] = prstar / Omega;
  P2[i] = P1[i] * prstar * prstar;
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

size_t ISCO_idx = gsl_interp_bsearch(r, commondata->r_ISCO, 0, commondata->nsteps_fine);
REAL t_peak = times[ISCO_idx] - commondata->Delta_t;
if (t_peak > times[commondata->nsteps_fine - 1]){
  t_peak = times[commondata->nsteps_fine - 1];
}
size_t peak_idx = gsl_interp_bsearch(times, t_peak, 0, commondata->nsteps_fine);


size_t N = 5;
size_t left = (peak_idx - N + abs(peak_idx - N))/2;
size_t right = (peak_idx + N + commondata->nsteps_fine - abs(peak_idx + N - commondata->nsteps_fine))/2;
REAL Q_cropped1[right - left], Q_cropped2[right - left], Q_cropped3[right - left], P_cropped1[right - left], P_cropped2[right - left];
REAL t_cropped[right-left], phase_cropped[right-left], amp_cropped[right-left];
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
gsl_matrix_set(P,0,0,gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(P,1,0,gsl_spline_eval_deriv(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, P_cropped2,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(P,0,1,gsl_spline_eval(spline, t_peak, acc));
gsl_matrix_set(P,1,1,gsl_spline_eval_deriv(spline, t_peak, acc));

gsl_vector *restrict A = gsl_vector_alloc(3);
gsl_vector *restrict O = gsl_vector_alloc(2);


gsl_spline_init(spline,t_cropped, amp_cropped,right-left);
gsl_interp_accel_reset(acc);
const REAL amp_insp = -1.*gsl_spline_eval(spline, t_peak, acc);
const REAL ampdot_insp = -1.*gsl_spline_eval_deriv(spline, t_peak, acc);
const REAL ampddot_insp = -1.*gsl_spline_eval_deriv2(spline, t_peak, acc);

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
gsl_vector_set(A , 2 , hddot_t_attach - ampdot_insp);
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

// apply the nqc correction
commondata->nsteps_combined = commondata->nsteps_low + commondata->nsteps_fine;
commondata->waveform_combined = (REAL *)calloc(commondata->nsteps_combined*NUMMODES, sizeof(REAL));
REAL nqc_amp, nqc_phase, q1, q2, q3,p1, p2;
for (i = 0; i < commondata->nsteps_low; i++){
  prstar = commondata->dynamics_low[IDX(i,PRSTAR)];
  radius = commondata->dynamics_low[IDX(i,R)];
  Omega = commondata->dynamics_low[IDX(i,OMEGA)];
  q1 = prstar * prstar / (radius * radius * Omega * Omega);
  q2 = q1 / radius;
  q3 = q2 / sqrt(radius);
  p1 = prstar / Omega;
  p2 = p1 * prstar * prstar;
  nqc_amp = 1 + commondata->a_1_NQC*q1 + commondata->a_2_NQC*q2 + commondata->a_3_NQC*q3;
  nqc_phase =  1 + commondata->b_1_NQC*p1 + commondata->b_2_NQC*p2;
  commondata->waveform_combined[IDX_WF(i,TIME)] = commondata->dynamics_low[IDX(i,TIME)];
  commondata->waveform_combined[IDX_WF(i,HPLUS)] = nqc_amp * (commondata->waveform_low[IDX_WF(i,HPLUS)]*cos(nqc_phase) + commondata->waveform_low[IDX_WF(i,HCROSS)]*sin(nqc_phase));
  commondata->waveform_combined[IDX_WF(i,HCROSS)] = nqc_amp * (commondata->waveform_low[IDX_WF(i,HCROSS)]*cos(nqc_phase) - commondata->waveform_low[IDX_WF(i,HPLUS)]*sin(nqc_phase));
}
for (i = 0; i< commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  radius = commondata->dynamics_fine[IDX(i,R)];
  Omega = commondata->dynamics_fine[IDX(i,OMEGA)];
  q1 = prstar * prstar / (radius * radius * Omega * Omega);
  q2 = q1 / radius;
  q3 = q2 / sqrt(radius);
  p1 = prstar / Omega;
  p2 = p1 * prstar * prstar;
  nqc_amp = 1 + commondata->a_1_NQC*q1 + commondata->a_2_NQC*q2 + commondata->a_3_NQC*q3;
  nqc_phase =  1 + commondata->b_1_NQC*p1 + commondata->b_2_NQC*p2;
  commondata->waveform_combined[IDX_WF(i+commondata->nsteps_low,TIME)] = commondata->dynamics_fine[IDX(i,TIME)];
  commondata->waveform_combined[IDX_WF(i+commondata->nsteps_low,HPLUS)] = nqc_amp * (commondata->waveform_fine[IDX_WF(i,HPLUS)]*cos(nqc_phase) + commondata->waveform_fine[IDX_WF(i,HCROSS)]*sin(nqc_phase));
  commondata->waveform_combined[IDX_WF(i+commondata->nsteps_low,HCROSS)] = nqc_amp * (commondata->waveform_fine[IDX_WF(i,HCROSS)]*cos(nqc_phase) - commondata->waveform_fine[IDX_WF(i,HPLUS)]*sin(nqc_phase));
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

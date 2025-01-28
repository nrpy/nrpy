"""
Set up C function library for SEOBNR and BOB attachment routines.

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


def register_CFunction_SEOBNRv5_aligned_spin_unwrap() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction that performs numpy.unwrap on a given array.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """C function to perform numpy.unwrap."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_unwrap"
    params = "REAL *restrict angles_in , REAL *restrict angles_out, size_t nsteps_arr"
    body = """
angles_out[0] = angles_in[0];
REAL diff;
for (size_t i = 1; i < nsteps_arr; i++){
  diff = angles_in[i] - angles_in[i-1];
  diff = fabs(diff) > M_PI ? (diff < - M_PI ? diff + 2 * M_PI : diff - 2 * M_PI) : diff;
  angles_out[i] = angles_out[i - 1] + diff;
}
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
REAL radius, omega, prstar; 
double complex h22;
size_t i;

for (i = 0; i < commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  h22 = commondata->waveform_fine[IDX_WF(i,STRAIN)];
  hamp[i] = cabs(h22);
  phase[i] = carg(h22);
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  Q1[i] = hamp[i] * prstar * prstar / (r[i] * r[i] * Omega[i] * Omega[i]);
  Q2[i] = Q1[i] / r[i];
  Q3[i] = Q2[i] / sqrt(r[i]);
  P1[i] = -prstar / r[i] /Omega[i];
  P2[i] = -P1[i] * prstar * prstar;
}
SEOBNRv5_aligned_spin_unwrap(phase,phase_unwrapped,commondata->nsteps_fine);
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
size_t peak_idx;
if (t_peak > times[commondata->nsteps_fine - 1]){
  t_peak = times[commondata->nsteps_fine - 1];
  peak_idx = commondata->nsteps_fine - 1;
}
else{
  peak_idx = gsl_interp_bsearch(times, t_peak, 0, commondata->nsteps_fine);
}
commondata->t_attach = t_peak;
size_t N = 5;
size_t left = MAX(peak_idx - N, 0);
size_t right = MIN(peak_idx + N , commondata->nsteps_fine);
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
gsl_matrix_set(P,0,0,-gsl_spline_eval_deriv(spline, t_peak, acc));
gsl_matrix_set(P,1,0,-gsl_spline_eval_deriv2(spline, t_peak, acc));

gsl_spline_init(spline,t_cropped, P_cropped2,right-left);
gsl_interp_accel_reset(acc);
gsl_matrix_set(P,0,1,-gsl_spline_eval_deriv(spline, t_peak, acc));
gsl_matrix_set(P,1,1,-gsl_spline_eval_deriv2(spline, t_peak, acc));

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

REAL omegas[2] , amps[3];
BOB_aligned_spin_NQC_rhs(commondata,amps,omegas);
gsl_vector_set(A , 0 , amps[0] - amp_insp);
gsl_vector_set(A , 1 , amps[1] - ampdot_insp);
gsl_vector_set(A , 2 , amps[2] - ampddot_insp);
gsl_vector_set(O , 0 , omegas[0] - omega_insp);
gsl_vector_set(O , 1 , omegas[1] - omegadot_insp);

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
free(Omega);
free(hamp);
free(phase);
free(phase_unwrapped);

// apply the nqc correction
commondata->nsteps_inspiral = commondata->nsteps_low + commondata->nsteps_fine;
commondata->waveform_inspiral = (double complex *)malloc(commondata->nsteps_inspiral*NUMMODES*sizeof(double complex));
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
  commondata->waveform_inspiral[IDX_WF(i,STRAIN)] = nqc_amp * commondata->waveform_low[IDX_WF(i,STRAIN)] *cexp(I * nqc_phase);
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
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,STRAIN)] = nqc_amp * commondata->waveform_fine[IDX_WF(i,STRAIN)] * cexp(I * nqc_phase);
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


def register_CFunction_SEOBNRv5_aligned_spin_interpolate_modes() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for interpolating the (2,2) inspiral mode for the SEOBNRv5 waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Interpolate the SEOBNRv5 (2,2) inspiral mode."""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_interpolate_modes"
    params = "commondata_struct *restrict commondata, const REAL dT"
    body = """
size_t i;
const size_t nsteps_inspiral_old = commondata->nsteps_inspiral;
const size_t nsteps_new = (size_t) (commondata->waveform_inspiral[IDX_WF(commondata->nsteps_inspiral - 1,TIME)] - commondata->waveform_inspiral[IDX_WF(0,TIME)]) / dT;
// build the complex inspiral modes
REAL *restrict times_old = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict orbital_phases = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict h22_nophase_real = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
REAL *restrict h22_nophase_imag = (REAL *)malloc(nsteps_inspiral_old * sizeof(REAL));
double complex h22_nophase , h22_rescaled;
for (i = 0; i < commondata->nsteps_low; i++){
  times_old[i] = commondata->waveform_inspiral[IDX_WF(i,TIME)];
  orbital_phases[i] = commondata->dynamics_low[IDX(i,PHI)];
  h22_nophase = cexp(2 * I * orbital_phases[i])*(commondata->waveform_inspiral[IDX_WF(i,STRAIN)]);
  h22_nophase_real[i] = creal(h22_nophase);
  h22_nophase_imag[i] = cimag(h22_nophase);
}
for (i = 0; i < commondata->nsteps_fine; i++){
  times_old[i + commondata->nsteps_low] = commondata->waveform_inspiral[IDX_WF(i + commondata->nsteps_low,TIME)];
  orbital_phases[i + commondata->nsteps_low] = commondata->dynamics_fine[IDX(i,PHI)];
  h22_nophase = cexp(2 * I * orbital_phases[i + commondata->nsteps_low])*(commondata->waveform_inspiral[IDX_WF(i + commondata->nsteps_low,STRAIN)]);
  h22_nophase_real[i + commondata->nsteps_low] = creal(h22_nophase);
  h22_nophase_imag[i + commondata->nsteps_low] = cimag(h22_nophase);
}

//interpolate and set the inspiral modes
REAL orbital_phase, h22_real , h22_imag, time;
const REAL tstart = commondata->waveform_inspiral[IDX_WF(0,TIME)];
gsl_interp_accel *restrict acc_real = gsl_interp_accel_alloc();
gsl_interp_accel *restrict acc_imag = gsl_interp_accel_alloc();
gsl_interp_accel *restrict acc = gsl_interp_accel_alloc();
gsl_spline *restrict spline_real = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline *restrict spline_imag = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_inspiral_old);
gsl_spline_init(spline,times_old,orbital_phases, nsteps_inspiral_old);
gsl_spline_init(spline_real,times_old,h22_nophase_real, nsteps_inspiral_old);
gsl_spline_init(spline_imag,times_old,h22_nophase_imag, nsteps_inspiral_old);
//realloc the amount of memory needed to store the interpolated modes
commondata->nsteps_inspiral = nsteps_new;
commondata->waveform_inspiral = (double complex *)realloc(commondata->waveform_inspiral,commondata->nsteps_inspiral * NUMMODES * sizeof(double complex));
for (i = 0; i < commondata->nsteps_inspiral; i++){
  time = tstart + i * dT;
  commondata->waveform_inspiral[IDX_WF(i,TIME)] = time;
  orbital_phase = gsl_spline_eval(spline,time,acc);
  h22_real = gsl_spline_eval(spline_real,time,acc_real);
  h22_imag = gsl_spline_eval(spline_imag,time,acc_imag);
  h22_rescaled = cexp(-2. * I * orbital_phase) * (h22_real + I * h22_imag);
  commondata->waveform_inspiral[IDX_WF(i,STRAIN)] = h22_rescaled;
}
gsl_interp_accel_free(acc);
gsl_interp_accel_free(acc_real);
gsl_interp_accel_free(acc_imag);
gsl_spline_free(spline_real);
gsl_spline_free(spline_imag);
gsl_spline_free(spline);
free(times_old);
free(orbital_phases);
free(h22_nophase_real);
free(h22_nophase_imag);
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
    body = """
size_t i;
const REAL dT = commondata->dt/(commondata->total_mass*4.925490947641266978197229498498379006e-6);
SEOBNRv5_aligned_spin_interpolate_modes(commondata , dT);
double complex h22;
REAL *restrict times_new = malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_amp_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_phase_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
REAL *restrict h22_wrapped_phase_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
for(i = 0; i < commondata->nsteps_inspiral; i++){
  times_new[i] = commondata->waveform_inspiral[IDX_WF(i,TIME)];
  h22 = commondata->waveform_inspiral[IDX_WF(i,STRAIN)];
  h22_amp_new[i] = cabs(h22);
  h22_wrapped_phase_new[i] = carg(h22);
}
SEOBNRv5_aligned_spin_unwrap(h22_wrapped_phase_new,h22_phase_new,commondata->nsteps_inspiral);
free(h22_wrapped_phase_new);
size_t idx_match = gsl_interp_bsearch(times_new,commondata->t_attach,0, commondata->nsteps_inspiral);
if (times_new[idx_match] > commondata->t_attach){
  idx_match--;
}
if (idx_match == commondata->nsteps_inspiral - 1){
  idx_match--;
}
const REAL t_match = times_new[idx_match];
const REAL phase_match = h22_phase_new[idx_match + 1];
const size_t nsteps_ringdown = 15 * (size_t) (commondata->tau_qnm / dT);
REAL *restrict ringdown_time = (REAL *)malloc(nsteps_ringdown*sizeof(REAL)); 
REAL *restrict ringdown_amp = (REAL *)malloc(nsteps_ringdown*sizeof(REAL));
REAL *restrict ringdown_phase = (REAL *)malloc(nsteps_ringdown*sizeof(REAL));
for(i = 0; i < nsteps_ringdown; i++){
  ringdown_time[i] = t_match + (i + 1) * dT;
}
BOB_aligned_spin_waveform_from_times(ringdown_time,ringdown_amp,ringdown_phase,nsteps_ringdown,commondata);
const REAL true_sign = copysign(phase_match,1.);
for(i = 0; i < nsteps_ringdown; i++){
  ringdown_phase[i] = true_sign*ringdown_phase[i] + phase_match;
}
commondata->nsteps_IMR = idx_match + 1 + nsteps_ringdown;
commondata->waveform_IMR = (double complex *)malloc(NUMMODES * commondata->nsteps_IMR*sizeof(double complex));
for (i = 0; i <= idx_match; i++){
  commondata->waveform_IMR[IDX_WF(i,TIME)] = times_new[i];
  h22 = h22_amp_new[i] * cexp(I * h22_phase_new[i]);
  commondata->waveform_IMR[IDX_WF(i,STRAIN)] = h22;
}
free(h22_amp_new);
free(h22_phase_new);
free(times_new);

for(i = 0; i < nsteps_ringdown; i++){
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_match,TIME)] = ringdown_time[i];
  h22 = ringdown_amp[i] * cexp(I * ringdown_phase[i]);
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_match,STRAIN)] = h22;
}
free(ringdown_time);
free(ringdown_phase);
free(ringdown_amp);
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

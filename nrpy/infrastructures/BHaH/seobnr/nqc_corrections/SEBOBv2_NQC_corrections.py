"""
Set up C function library for SEBOBv2 NQC attachment routines.

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


def register_CFunction_SEBOBv2_NQC_corrections() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the Non Quasi-Circular (NQC) corrections for the SEBOBv2 waveform.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Computes and applies the Non Quasi-Circular (NQC) corrections to the inspiral waveform.

@param commondata - Common data structure containing the model parameters.
@return - Returns GSL_SUCCESS (0) on success or a nonzero error code on failure.
"""
    cfunc_type = "void"
    name = "SEBOBv2_NQC_corrections"
    params = "commondata_struct *restrict commondata"
    body = """
REAL *restrict times = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (times == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for times\\n");
  exit(1);
}
REAL *restrict Q1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q1 == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for Q1\\n");
  exit(1);
}
REAL *restrict Q2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q2 == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for Q2\\n");
  exit(1);
}
REAL *restrict Q3 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q3 == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for Q3\\n");
  exit(1);
}
REAL *restrict P1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (P1 == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for P1\\n");
  exit(1);
}
REAL *restrict P2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (P2 == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for P2\\n");
  exit(1);
}
REAL *restrict r = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (r == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for r\\n");
  exit(1);
}
REAL *restrict Omega = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Omega == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for Omega\\n");
  exit(1);
}
REAL *restrict hamp = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (hamp == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for hamp\\n");
  exit(1);
}
REAL *restrict phase = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (phase == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for phase\\n");
  exit(1);
}
REAL *restrict phase_unwrapped = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (phase_unwrapped == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for phase_unwrapped\\n");
  exit(1);
}
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
  if (t_zoom == NULL){
    fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for t_zoom\\n");
    exit(1);
  }
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (minus_r_zoom == NULL){
    fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for times\\n");
    exit(1);
  }
  gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  if (acc_r == NULL){
    fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
    exit(1);
  }
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_r == NULL){
    fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_spline_alloc() failed to initialize\\n");
    exit(1);
  }
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
// if t_peak > the last point of the ODE trajector,
// use the second last point in time for t_peak, instead of the last point as in pySEOBNR.
// gsl's cubic spline uses natural boundary conditions that sets second derivatives to zero
// resulting in a singular matrix.
if (t_peak > times[commondata->nsteps_fine - 1]){
  t_peak = times[commondata->nsteps_fine - 2];
  peak_idx = commondata->nsteps_fine - 2;
}
else{
  peak_idx = gsl_interp_bsearch(times, t_peak, 0, commondata->nsteps_fine);
}
commondata->t_attach = t_peak;
// Compute t_p and Omega_0 in BOB
//not needed anymore
//BOB_v2_find_tp_Omega0(commondata);
BOB_v2_setup_peak_attachment(commondata);
// Compute Q_cropped, P_cropped, t_cropped, phase_cropped, amp_cropped
size_t N = 5;
size_t left = MAX(peak_idx, N) - N;
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
if (acc == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_interp_accel_alloc() failed to initialize\\n");
  exit(1);
}
gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, right - left);
if (spline == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_splin_alloc() failed to initialize\\n");
  exit(1);
}

gsl_matrix *restrict Q = gsl_matrix_alloc (3, 3);
if (Q == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_matrix_alloc() failed to initialize\\n");
  exit(1);
}
gsl_matrix *restrict P = gsl_matrix_alloc (2, 2);
if (P == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_matrix_alloc() failed to initialize\\n");
  exit(1);
}

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
if (A == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}
gsl_vector *restrict O = gsl_vector_alloc(2);
if (O == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}


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
BOB_v2_NQC_rhs(commondata,amps,omegas);
gsl_vector_set(A , 0 , amps[0] - amp_insp);
gsl_vector_set(A , 1 , amps[1] - ampdot_insp);
gsl_vector_set(A , 2 , amps[2] - ampddot_insp);
gsl_vector_set(O , 0 , omegas[0] - omega_insp);
gsl_vector_set(O , 1 , omegas[1] - omegadot_insp);

gsl_vector *restrict a = gsl_vector_alloc (3);
if (a == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}

int s;
gsl_permutation *restrict p_A = gsl_permutation_alloc (3);
if (p_A == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_permutation_alloc() failed to initialize\\n");
  exit(1);
}
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
if (b == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}
gsl_permutation * p_B = gsl_permutation_alloc(2);
if (p_B == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), gsl_permutation_alloc() failed to initialize\\n");
  exit(1);
}
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
if (commondata->waveform_inspiral == NULL){
  fprintf(stderr,"Error: in SEBOBv2_NQC_corrections(), malloc() failed for commondata->waveform_inspiral\\n");
  exit(1);
}
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
"""
    cfc.register_CFunction(
        subdirectory="nqc_corrections",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

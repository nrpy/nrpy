"""
Set up C function library for SEOBNR and BOB attachment routines.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Shreya Reddy
        redd7331 **at** vandals **dot** uidaho **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.equations.seobnr.BOB_aligned_spin_waveform_quantities_higher_modes as BOB_HM_wf
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_BOB_aligned_spin_NQC_corrections_higher_modes(
    use_numerical_relativity_nqc: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the Non Quasi-Circular (NQC) corrections for the SEOBNRv5 waveform.

    :param use_numerical_relativity_nqc: Flag to specify if NQCs are informed by NR fits or BOB
    :raises ValueError: if higher order NR-informed nqc corrections are called
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    wf = BOB_HM_wf.BOB_aligned_spin_waveform_quantities_higher_modes()
    modes = wf.modes
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Solves and applies BOB-informed higher-mode Non Quasi-Circular (NQC) corrections.

For each stored inspiral mode, this matches amplitude and waveform-frequency
data at the shared attachment time to BOB-informed targets and applies the
resulting NQC amplitude and phase corrections.

@param[in,out] commondata Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "BOB_aligned_spin_NQC_corrections_higher_modes"
    params = "commondata_struct *restrict commondata"
    body = """
REAL *restrict times = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (times == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for times\\n");
  exit(1);
} // END IF: times allocation failed
REAL *restrict Q1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q1 == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for Q1\\n");
  exit(1);
} // END IF: Q1 allocation failed
REAL *restrict Q2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q2 == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for Q2\\n");
  exit(1);
}  // END IF: Q2 allocation failed
REAL *restrict Q3 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Q3 == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for Q3\\n");
  exit(1);
}  // END IF: Q3 allocation failed
REAL *restrict P1 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (P1 == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for P1\\n");
  exit(1);
}  // END IF: P1 allocation failed
REAL *restrict P2 = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (P2 == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for P2\\n");
  exit(1);
}  // END IF: P2 allocation failed
REAL *restrict r = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (r == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for r\\n");
  exit(1);
}  // END IF: r allocation failed
REAL *restrict Omega = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (Omega == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for Omega\\n");
  exit(1);
}  // END IF: Omega allocation failed
REAL *restrict hamp = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (hamp == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for hamp\\n");
  exit(1);
}  // END IF: hamp allocation failed
REAL *restrict phase = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (phase == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for phase\\n");
  exit(1);
}  // END IF: phase allocation failed
REAL *restrict phase_unwrapped = (REAL *)malloc(commondata->nsteps_fine*sizeof(REAL));
if (phase_unwrapped == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for phase_unwrapped\\n");
  exit(1);
}  // END IF: phase_unwrapped allocation failed
REAL radius, omega, prstar;
size_t i;
"""
    body += """
for (i = 0; i < commondata->nsteps_fine; i++){
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  times[i] = commondata->dynamics_fine[IDX(i,TIME)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  P1[i] = -prstar / r[i] /Omega[i];
  P2[i] = -P1[i] * prstar * prstar;
} // END LOOP: for i over fine dynamics samples

// Step 2: Locate the ISCO crossing used to set the common attachment time.
if (commondata->r_ISCO < r[commondata->nsteps_fine - 1]){
  commondata->t_ISCO = times[commondata->nsteps_fine - 1];
} // END IF: final fine sample is inside r_ISCO
else{
  const REAL dt_ISCO = 0.001;
  const size_t N_zoom = (size_t) ((times[commondata->nsteps_fine - 1] - times[0]) / dt_ISCO);
  REAL *restrict t_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (t_zoom == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for t_zoom\\n");
  exit(1);
}  // END IF: t_zoom allocation failed
  REAL *restrict minus_r_zoom = (REAL *) malloc(N_zoom * sizeof(REAL));
  if (minus_r_zoom == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for minus_r_zoom\\n");
  exit(1);
}  // END IF: minus_r_zoom allocation failed
  gsl_interp_accel *restrict acc_r = gsl_interp_accel_alloc();
  if (acc_r == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for acc_r\\n");
  exit(1);
}  // END IF: acc_r allocation failed
  gsl_spline *restrict spline_r = gsl_spline_alloc(gsl_interp_cspline, commondata->nsteps_fine);
  if (spline_r == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for spline_r\\n");
  exit(1);
}  // END IF: spline_r allocation failed
  gsl_spline_init(spline_r,times,r,commondata->nsteps_fine);
  for (i = 0; i < N_zoom; i++){
    t_zoom[i] = times[0] + i * dt_ISCO;
    minus_r_zoom[i] = -1.0*gsl_spline_eval(spline_r,t_zoom[i],acc_r);
  } // END LOOP: for i over zoomed ISCO samples
  const size_t ISCO_zoom_idx = gsl_interp_bsearch(minus_r_zoom, -commondata->r_ISCO, 0 , N_zoom);
  commondata->t_ISCO = t_zoom[ISCO_zoom_idx];

  gsl_interp_accel_free(acc_r);
  gsl_spline_free(spline_r);
  free(t_zoom);
  free(minus_r_zoom);
} // END ELSE: interpolate ISCO crossing on zoomed grid

// Step 3: Select the attachment time and shared interpolation window.
REAL t_peak = commondata->t_ISCO - commondata->Delta_t;
size_t peak_idx;
// if t_peak > the last point of the ODE trajectory,
// use the second last point in time for t_peak, instead of the last point as in pySEOBNR.
// gsl's cubic spline uses natural boundary conditions that sets second derivatives to zero
// resulting in a singular matrix.
if (t_peak > times[commondata->nsteps_fine - 1]){
  t_peak = times[commondata->nsteps_fine - 2];
  peak_idx = commondata->nsteps_fine - 2;
} // END IF: attachment time is beyond final fine sample
else{
  peak_idx = gsl_interp_bsearch(times, t_peak, 0, commondata->nsteps_fine);
} // END ELSE: attachment time lies inside fine dynamics
commondata->t_attach = t_peak;
size_t N = 5;
size_t left = NRPYMAX(peak_idx, N) - N;
size_t right = NRPYMIN(peak_idx + N , commondata->nsteps_fine);
"""
    for lm in modes:
        body += f"""
double complex h{lm[0]}{lm[1]};

for (i = 0; i < commondata->nsteps_fine; i++){{
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  r[i] = commondata->dynamics_fine[IDX(i,R)];
  Omega[i] = commondata->dynamics_fine[IDX(i,OMEGA)];
  h{lm[0]}{lm[1]} = commondata->waveform_fine[IDX_WF(i,STRAIN{lm[0]}{lm[1]})];
  hamp[i] = cabs(h{lm[0]}{lm[1]});
  phase[i] = carg(h{lm[0]}{lm[1]});
  Q1[i] = hamp[i] * prstar * prstar / (r[i] * r[i] * Omega[i] * Omega[i]);
  Q2[i] = Q1[i] / r[i];
  Q3[i] = Q2[i] / sqrt(r[i]);
}} // END LOOP: for i over fine dynamics
SEOBNRv5_aligned_spin_unwrap(phase,phase_unwrapped,commondata->nsteps_fine);

REAL Q{lm[0]}{lm[1]}_cropped1[right - left], Q{lm[0]}{lm[1]}_cropped2[right - left], Q{lm[0]}{lm[1]}_cropped3[right - left], P{lm[0]}{lm[1]}_cropped1[right - left], P{lm[0]}{lm[1]}_cropped2[right - left];
REAL t_cropped{lm[0]}{lm[1]}[right-left], phase_cropped{lm[0]}{lm[1]}[right - left], amp_cropped{lm[0]}{lm[1]}[right-left];
for (i = left; i < right; i++){{
  t_cropped{lm[0]}{lm[1]}[i - left] = times[i];
  phase_cropped{lm[0]}{lm[1]}[i - left] = phase_unwrapped[i];
  amp_cropped{lm[0]}{lm[1]}[i - left] = hamp[i];
  Q{lm[0]}{lm[1]}_cropped1[i - left] = Q1[i];
  Q{lm[0]}{lm[1]}_cropped2[i - left] = Q2[i];
  Q{lm[0]}{lm[1]}_cropped3[i - left] = Q3[i];
  P{lm[0]}{lm[1]}_cropped1[i - left] = P1[i];
  P{lm[0]}{lm[1]}_cropped2[i - left] = P2[i];
}} // END LOOP: for i over cropped attachment samples
gsl_interp_accel *restrict acc{lm[0]}{lm[1]} = gsl_interp_accel_alloc();
if (acc{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_interp_accel_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: acc{lm[0]}{lm[1]} allocation failed
gsl_spline *restrict spline{lm[0]}{lm[1]} = gsl_spline_alloc(gsl_interp_cspline, right - left);
if (spline{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_spline_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: spline{lm[0]}{lm[1]} allocation failed

gsl_matrix *restrict Q{lm[0]}{lm[1]} = gsl_matrix_alloc (3, 3);
if (Q{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_matrix_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: Q{lm[0]}{lm[1]} allocation failed
gsl_matrix *restrict P{lm[0]}{lm[1]} = gsl_matrix_alloc (2, 2);
if (P{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_matrix_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: P{lm[0]}{lm[1]} allocation failed
gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, Q{lm[0]}{lm[1]}_cropped1,right-left);
gsl_matrix_set(Q{lm[0]}{lm[1]},0,0,gsl_spline_eval(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},1,0,gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},2,0,gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, Q{lm[0]}{lm[1]}_cropped2,right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
gsl_matrix_set(Q{lm[0]}{lm[1]},0,1,gsl_spline_eval(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},1,1,gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},2,1,gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, Q{lm[0]}{lm[1]}_cropped3,right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
gsl_matrix_set(Q{lm[0]}{lm[1]},0,2,gsl_spline_eval(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},1,2,gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(Q{lm[0]}{lm[1]},2,2,gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, P{lm[0]}{lm[1]}_cropped1,right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
gsl_matrix_set(P{lm[0]}{lm[1]},0,0,-gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(P{lm[0]}{lm[1]},1,0,-gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, P{lm[0]}{lm[1]}_cropped2,right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
gsl_matrix_set(P{lm[0]}{lm[1]},0,1,-gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));
gsl_matrix_set(P{lm[0]}{lm[1]},1,1,-gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]}));

gsl_vector *restrict A{lm[0]}{lm[1]} = gsl_vector_alloc(3);
if (A{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: A{lm[0]}{lm[1]} allocation failed
gsl_vector *restrict O{lm[0]}{lm[1]} = gsl_vector_alloc(2);
if (O{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: O{lm[0]}{lm[1]} allocation failed

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, amp_cropped{lm[0]}{lm[1]},right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
const REAL amp_insp{lm[0]}{lm[1]} = gsl_spline_eval(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]});
const REAL ampdot_insp{lm[0]}{lm[1]} = gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]});
const REAL ampddot_insp{lm[0]}{lm[1]} = gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]});

gsl_spline_init(spline{lm[0]}{lm[1]},t_cropped{lm[0]}{lm[1]}, phase_cropped{lm[0]}{lm[1]},right-left);
gsl_interp_accel_reset(acc{lm[0]}{lm[1]});
REAL omega_insp{lm[0]}{lm[1]} = gsl_spline_eval_deriv(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]});
REAL omegadot_insp{lm[0]}{lm[1]} =  gsl_spline_eval_deriv2(spline{lm[0]}{lm[1]}, t_peak, acc{lm[0]}{lm[1]});

gsl_spline_free(spline{lm[0]}{lm[1]});
gsl_interp_accel_free(acc{lm[0]}{lm[1]});

if (omega_insp{lm[0]}{lm[1]} * omegadot_insp{lm[0]}{lm[1]} > 0.0){{
  omega_insp{lm[0]}{lm[1]} = fabs(omega_insp{lm[0]}{lm[1]});
  omegadot_insp{lm[0]}{lm[1]} = fabs(omegadot_insp{lm[0]}{lm[1]});
}} // END IF: waveform frequency and derivative have matching signs
else{{
  omega_insp{lm[0]}{lm[1]} = fabs(omega_insp{lm[0]}{lm[1]});
  omegadot_insp{lm[0]}{lm[1]} = -fabs(omegadot_insp{lm[0]}{lm[1]});
}} // END ELSE: waveform frequency and derivative have opposing signs
"""
    if use_numerical_relativity_nqc:
        raise ValueError(
            "NR informed higher mode NQC corrections are not currently implemented"
        )

    body += """
REAL omegas[NUMMODES][2] , amps[NUMMODES][3];
REAL a_nqc[NUMMODES][3], b_nqc[NUMMODES][2];

BOB_aligned_spin_NQC_rhs_HM(commondata,amps,omegas);
"""

    for lm in modes:
        body += f"""
gsl_vector_set(A{lm[0]}{lm[1]} , 0 , amps[HNR{lm[0]}{lm[1]}][0] - amp_insp{lm[0]}{lm[1]});
gsl_vector_set(A{lm[0]}{lm[1]} , 1 , amps[HNR{lm[0]}{lm[1]}][1] - ampdot_insp{lm[0]}{lm[1]});
gsl_vector_set(A{lm[0]}{lm[1]} , 2 , amps[HNR{lm[0]}{lm[1]}][2] - ampddot_insp{lm[0]}{lm[1]});
gsl_vector_set(O{lm[0]}{lm[1]} , 0 , omegas[HNR{lm[0]}{lm[1]}][0] - omega_insp{lm[0]}{lm[1]});
gsl_vector_set(O{lm[0]}{lm[1]} , 1 , omegas[HNR{lm[0]}{lm[1]}][1] - omegadot_insp{lm[0]}{lm[1]});

gsl_vector *restrict a{lm[0]}{lm[1]} = gsl_vector_alloc (3);
if (a{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: a{lm[0]}{lm[1]} allocation failed
int s{lm[0]}{lm[1]};
gsl_permutation *restrict p_A{lm[0]}{lm[1]} = gsl_permutation_alloc (3);
if (p_A{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_permutation_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: p_A{lm[0]}{lm[1]} allocation failed
gsl_linalg_LU_decomp(Q{lm[0]}{lm[1]}, p_A{lm[0]}{lm[1]}, &s{lm[0]}{lm[1]});
gsl_linalg_LU_solve(Q{lm[0]}{lm[1]}, p_A{lm[0]}{lm[1]}, A{lm[0]}{lm[1]}, a{lm[0]}{lm[1]});
gsl_permutation_free(p_A{lm[0]}{lm[1]});

a_nqc[HNR{lm[0]}{lm[1]}][0] = gsl_vector_get(a{lm[0]}{lm[1]},0);
a_nqc[HNR{lm[0]}{lm[1]}][1] = gsl_vector_get(a{lm[0]}{lm[1]},1);
a_nqc[HNR{lm[0]}{lm[1]}][2] = gsl_vector_get(a{lm[0]}{lm[1]},2);

gsl_vector_free(a{lm[0]}{lm[1]});
gsl_vector_free(A{lm[0]}{lm[1]});


gsl_vector *restrict b{lm[0]}{lm[1]} = gsl_vector_alloc(2);
if (b{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_vector_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: b{lm[0]}{lm[1]} allocation failed
gsl_permutation * p_B{lm[0]}{lm[1]} = gsl_permutation_alloc(2);
if (p_B{lm[0]}{lm[1]} == NULL){{
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), gsl_permutation_alloc() failed to initialize\\n");
  exit(1);
}} // END IF: p_B{lm[0]}{lm[1]} allocation failed
gsl_linalg_LU_decomp(P{lm[0]}{lm[1]}, p_B{lm[0]}{lm[1]}, &s{lm[0]}{lm[1]});
gsl_linalg_LU_solve (P{lm[0]}{lm[1]}, p_B{lm[0]}{lm[1]}, O{lm[0]}{lm[1]}, b{lm[0]}{lm[1]});
gsl_permutation_free (p_B{lm[0]}{lm[1]});

b_nqc[HNR{lm[0]}{lm[1]}][0]= gsl_vector_get(b{lm[0]}{lm[1]},0);
b_nqc[HNR{lm[0]}{lm[1]}][1] = gsl_vector_get(b{lm[0]}{lm[1]},1);

gsl_vector_free(b{lm[0]}{lm[1]});
gsl_vector_free(O{lm[0]}{lm[1]});
gsl_matrix_free(P{lm[0]}{lm[1]});
gsl_matrix_free(Q{lm[0]}{lm[1]});


"""

    body += """
    
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
    

// Step 4: Apply the higher-mode NQC correction to low and fine waveform samples.
commondata->nsteps_inspiral = commondata->nsteps_low + commondata->nsteps_fine;
commondata->waveform_inspiral = (double complex *)malloc(commondata->nsteps_inspiral*NUMMODES*sizeof(double complex));
if (commondata->waveform_inspiral == NULL){
  fprintf(stderr,"Error: in BOB_aligned_spin_NQC_corrections_higher_modes(), malloc() failed for commondata->waveform_inspiral\\n");
  exit(1);
} // END IF: waveform_inspiral allocation failed
for (i = 0; i < commondata->nsteps_low; i++){
  commondata->waveform_inspiral[IDX_WF(i,TIME)] = commondata->dynamics_low[IDX(i,TIME)];
} // END LOOP: for i over low-resolution waveform times
for (i = 0; i< commondata->nsteps_fine; i++){
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,TIME)] = commondata->dynamics_fine[IDX(i,TIME)];
} // END LOOP: for i over fine-resolution waveform times
"""
    for lm in modes:
        body += f"""


REAL nqc_amp{lm[0]}{lm[1]}, nqc_phase{lm[0]}{lm[1]}, q1_{lm[0]}{lm[1]}, q2_{lm[0]}{lm[1]}, q3_{lm[0]}{lm[1]},p1_{lm[0]}{lm[1]}, p2_{lm[0]}{lm[1]};
for (i = 0; i < commondata->nsteps_low; i++){{
  prstar = commondata->dynamics_low[IDX(i,PRSTAR)];
  radius = commondata->dynamics_low[IDX(i,R)];
  omega = commondata->dynamics_low[IDX(i,OMEGA)];
  q1_{lm[0]}{lm[1]} = prstar * prstar / (radius * radius * omega * omega);
  q2_{lm[0]}{lm[1]} = q1_{lm[0]}{lm[1]} / radius;
  q3_{lm[0]}{lm[1]} = q2_{lm[0]}{lm[1]} / sqrt(radius);
  p1_{lm[0]}{lm[1]} = -prstar / radius /omega;
  p2_{lm[0]}{lm[1]} = -p1_{lm[0]}{lm[1]} * prstar * prstar;
  nqc_amp{lm[0]}{lm[1]} = 1 + a_nqc[HNR{lm[0]}{lm[1]}][0]*q1_{lm[0]}{lm[1]} + a_nqc[HNR{lm[0]}{lm[1]}][1]*q2_{lm[0]}{lm[1]} + a_nqc[HNR{lm[0]}{lm[1]}][2]*q3_{lm[0]}{lm[1]};
  nqc_phase{lm[0]}{lm[1]} =  b_nqc[HNR{lm[0]}{lm[1]}][0]*p1_{lm[0]}{lm[1]} + b_nqc[HNR{lm[0]}{lm[1]}][1]*p2_{lm[0]}{lm[1]};
  commondata->waveform_inspiral[IDX_WF(i,STRAIN{lm[0]}{lm[1]})] = nqc_amp{lm[0]}{lm[1]} * commondata->waveform_low[IDX_WF(i,STRAIN{lm[0]}{lm[1]})] *cexp(I * nqc_phase{lm[0]}{lm[1]});
}} // END LOOP: for i over low-resolution waveform samples
for (i = 0; i< commondata->nsteps_fine; i++){{
  prstar = commondata->dynamics_fine[IDX(i,PRSTAR)];
  radius = commondata->dynamics_fine[IDX(i,R)];
  omega = commondata->dynamics_fine[IDX(i,OMEGA)];
  q1_{lm[0]}{lm[1]} = prstar * prstar / (radius * radius * omega * omega);
  q2_{lm[0]}{lm[1]} = q1_{lm[0]}{lm[1]} / radius;
  q3_{lm[0]}{lm[1]} = q2_{lm[0]}{lm[1]} / sqrt(radius);
  p1_{lm[0]}{lm[1]} = -prstar / radius /omega;
  p2_{lm[0]}{lm[1]} = -p1_{lm[0]}{lm[1]} * prstar * prstar;
  nqc_amp{lm[0]}{lm[1]} = 1 + a_nqc[HNR{lm[0]}{lm[1]}][0]*q1_{lm[0]}{lm[1]} + a_nqc[HNR{lm[0]}{lm[1]}][1]*q2_{lm[0]}{lm[1]} + a_nqc[HNR{lm[0]}{lm[1]}][2]*q3_{lm[0]}{lm[1]};
  nqc_phase{lm[0]}{lm[1]} =  b_nqc[HNR{lm[0]}{lm[1]}][0]*p1_{lm[0]}{lm[1]} + b_nqc[HNR{lm[0]}{lm[1]}][1]*p2_{lm[0]}{lm[1]};
  commondata->waveform_inspiral[IDX_WF(i+commondata->nsteps_low,STRAIN{lm[0]}{lm[1]})] = nqc_amp{lm[0]}{lm[1]} * commondata->waveform_fine[IDX_WF(i,STRAIN{lm[0]}{lm[1]})] * cexp(I * nqc_phase{lm[0]}{lm[1]});
}} // END LOOP: for i over fine-resolution waveform samples

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

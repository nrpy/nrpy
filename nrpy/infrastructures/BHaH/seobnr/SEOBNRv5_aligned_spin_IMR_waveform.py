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


def register_CFunction_SEOBNRv5_aligned_spin_IMR_waveform(
    use_seobnrv5_merger_ringdown: bool,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for evaluating the (2,2) IMR mode for the SEOBNRv5 waveform.

    :param use_seobnrv5_merger_ringdown: Flag to specify whether or not to use the native SEOBNRv5 merger ringdown model instead of BOB.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Builds the aligned-spin (2,2) inspiral-merger-ringdown waveform.

The inspiral is interpolated through the attachment time. The ringdown segment
uses either the native SEOBNRv5 merger-ringdown model or the BOB model,
depending on how this function was registered.

@param[in,out] commondata Common data structure containing the model parameters.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_aligned_spin_IMR_waveform"
    params = "commondata_struct *restrict commondata"
    body = """
// Step 1: Interpolate inspiral modes and allocate amplitude and phase arrays.
size_t i;
const REAL dT = commondata->dt/(commondata->total_mass*4.925490947641266978197229498498379006e-6);
SEOBNRv5_aligned_spin_interpolate_modes(commondata , dT);
double complex h22;
REAL *restrict times_new = malloc(commondata->nsteps_inspiral*sizeof(REAL));
if (times_new == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for times_new\\n");
  exit(1);
}
REAL *restrict h22_amp_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
if (h22_amp_new == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for h22_amp_new\\n");
  exit(1);
}
REAL *restrict h22_phase_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
if (h22_phase_new == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for h22_phase_new\\n");
  exit(1);
}
REAL *restrict h22_wrapped_phase_new = (REAL *)malloc(commondata->nsteps_inspiral*sizeof(REAL));
if (h22_wrapped_phase_new == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for h22_wrapped_phase_new\\n");
  exit(1);
}
// Step 2: Extract interpolated inspiral amplitude and phase samples.
for(i = 0; i < commondata->nsteps_inspiral; i++){
  times_new[i] = commondata->waveform_inspiral[IDX_WF(i,TIME)];
  h22 = commondata->waveform_inspiral[IDX_WF(i,STRAIN)];
  h22_amp_new[i] = cabs(h22);
  h22_wrapped_phase_new[i] = carg(h22);
} // END LOOP: for i over inspiral samples
SEOBNRv5_aligned_spin_unwrap(h22_wrapped_phase_new,h22_phase_new,commondata->nsteps_inspiral);
free(h22_wrapped_phase_new);

// Step 3: Select matching index and allocate ringdown arrays.
size_t idx_match = gsl_interp_bsearch(times_new,commondata->t_attach,0, commondata->nsteps_inspiral);
if (times_new[idx_match] > commondata->t_attach){
  idx_match--;
}
if (idx_match == commondata->nsteps_inspiral - 1){
  idx_match--;
}
const REAL t_match = times_new[idx_match];
const size_t nsteps_ringdown = 15 * (size_t) (commondata->tau_qnm / dT);
REAL *restrict ringdown_time = (REAL *)malloc(nsteps_ringdown*sizeof(REAL));
if (ringdown_time == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for ringdown_time\\n");
  exit(1);
}
REAL *restrict ringdown_amp = (REAL *)malloc(nsteps_ringdown*sizeof(REAL));
if (ringdown_amp == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for ringdown_amp\\n");
  exit(1);
}
REAL *restrict ringdown_phase = (REAL *)malloc(nsteps_ringdown*sizeof(REAL));
if (ringdown_phase == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for ringdown_phase\\n");
  exit(1);
}
// Step 4: Build ringdown sample times.
for(i = 0; i < nsteps_ringdown; i++){
  ringdown_time[i] = t_match + (i + 1) * dT;
} // END LOOP: for i over ringdown sample times
"""
    if use_seobnrv5_merger_ringdown:
        body += """
// Step 5: Build native merger-ringdown initial data at the match time.
size_t left = NRPYMAX(5,idx_match) - 5;
size_t right = NRPYMIN(commondata->nsteps_inspiral,idx_match + 5);
REAL times_cropped[right-left] , amps_cropped[right-left] , phases_cropped[right-left];
for (i = left; i < right; i++){
  times_cropped[i - left] = times_new[i];
  amps_cropped[i - left] = h22_amp_new[i];
  phases_cropped[i - left] = h22_phase_new[i];
} // END LOOP: for i over cropped inspiral samples
gsl_interp_accel *restrict acc = gsl_interp_accel_alloc();
if (acc == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), gsl_interp_accel_alloc failed to initialize\\n");
  exit(1);
}
gsl_spline *restrict spline = gsl_spline_alloc(gsl_interp_cspline, right - left);
if (spline == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), gsl_spline_alloc failed to initialize\\n");
  exit(1);
}
gsl_spline_init(spline,times_cropped, amps_cropped,right-left);
const REAL h_0 = gsl_spline_eval(spline, t_match, acc);
const REAL hdot_0 = gsl_spline_eval_deriv(spline, t_match, acc);

gsl_spline_init(spline,times_cropped, phases_cropped,right-left);
gsl_interp_accel_reset(acc);
const REAL phi_0 = gsl_spline_eval(spline, t_match, acc);
const REAL phidot_0 = gsl_spline_eval_deriv(spline, t_match, acc);

gsl_spline_free(spline);
gsl_interp_accel_free(acc);


SEOBNRv5_aligned_spin_merger_waveform_from_times(ringdown_time,ringdown_amp,ringdown_phase,t_match,h_0,hdot_0,phi_0,phidot_0,nsteps_ringdown,commondata);
"""
    else:
        body += """
// Step 5: Build BOB ringdown samples and align their phase.
const REAL phase_match = h22_phase_new[idx_match + 1];
BOB_aligned_spin_waveform_from_times(ringdown_time,ringdown_amp,ringdown_phase,nsteps_ringdown,commondata);
const REAL true_sign = copysign(1.,phase_match);
for(i = 0; i < nsteps_ringdown; i++){
  ringdown_phase[i] = true_sign*ringdown_phase[i] + phase_match;
} // END LOOP: for i over BOB ringdown phases
"""
    body += """
// Step 6: Combine inspiral and ringdown into the IMR waveform.
commondata->nsteps_IMR = idx_match + 1 + nsteps_ringdown;
commondata->waveform_IMR = (double complex *)malloc(NUMMODES * commondata->nsteps_IMR*sizeof(double complex));
if (commondata->waveform_IMR == NULL){
  fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_IMR_waveform(), malloc() failed to for commondata->waveform_IMR\\n");
  exit(1);
}
for (i = 0; i <= idx_match; i++){
  commondata->waveform_IMR[IDX_WF(i,TIME)] = times_new[i] - commondata->t_attach;
  h22 = h22_amp_new[i] * cexp(I * h22_phase_new[i]);
  commondata->waveform_IMR[IDX_WF(i,STRAIN)] = h22;
} // END LOOP: for i over inspiral IMR samples
free(h22_amp_new);
free(h22_phase_new);
free(times_new);

for(i = 0; i < nsteps_ringdown; i++){
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_match,TIME)] = ringdown_time[i] - commondata->t_attach;
  h22 = ringdown_amp[i] * cexp(I * ringdown_phase[i]);
  commondata->waveform_IMR[IDX_WF(i + 1 + idx_match,STRAIN)] = h22;
} // END LOOP: for i over ringdown IMR samples
free(ringdown_time);
free(ringdown_phase);
free(ringdown_amp);
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
    return pcg.NRPyEnv()

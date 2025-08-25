#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate SEOBNRv5 (2,2) mode.
 */
int SEOBNRv5_aligned_spin_IMR_waveform(commondata_struct *restrict commondata) {

  size_t i;
  const REAL dT = commondata->dt / (commondata->total_mass * 4.925490947641266978197229498498379006e-6);
  SEOBNRv5_aligned_spin_interpolate_modes(commondata, dT);
  double complex h22;
  REAL *restrict times_new = malloc(commondata->nsteps_inspiral * sizeof(REAL));
  REAL *restrict h22_amp_new = (REAL *)malloc(commondata->nsteps_inspiral * sizeof(REAL));
  REAL *restrict h22_phase_new = (REAL *)malloc(commondata->nsteps_inspiral * sizeof(REAL));
  REAL *restrict h22_wrapped_phase_new = (REAL *)malloc(commondata->nsteps_inspiral * sizeof(REAL));
  for (i = 0; i < commondata->nsteps_inspiral; i++) {
    times_new[i] = commondata->waveform_inspiral[IDX_WF(i, TIME)];
    h22 = commondata->waveform_inspiral[IDX_WF(i, STRAIN)];
    h22_amp_new[i] = cabs(h22);
    h22_wrapped_phase_new[i] = carg(h22);
  }
  SEOBNRv5_aligned_spin_unwrap(h22_wrapped_phase_new, h22_phase_new, commondata->nsteps_inspiral);
  free(h22_wrapped_phase_new);
  size_t idx_match = gsl_interp_bsearch(times_new, commondata->t_attach, 0, commondata->nsteps_inspiral);
  if (times_new[idx_match] > commondata->t_attach) {
    idx_match--;
  }
  if (idx_match == commondata->nsteps_inspiral - 1) {
    idx_match--;
  }
  const REAL t_match = times_new[idx_match];
  const size_t nsteps_ringdown = 15 * (size_t)(commondata->tau_qnm / dT);
  REAL *restrict ringdown_time = (REAL *)malloc(nsteps_ringdown * sizeof(REAL));
  REAL *restrict ringdown_amp = (REAL *)malloc(nsteps_ringdown * sizeof(REAL));
  REAL *restrict ringdown_phase = (REAL *)malloc(nsteps_ringdown * sizeof(REAL));
  for (i = 0; i < nsteps_ringdown; i++) {
    ringdown_time[i] = t_match + (i + 1) * dT;
  }

  const REAL phase_match = h22_phase_new[idx_match + 1];
  BOB_aligned_spin_waveform_from_times(ringdown_time, ringdown_amp, ringdown_phase, nsteps_ringdown, commondata);
  const REAL true_sign = copysign(1., phase_match);
  for (i = 0; i < nsteps_ringdown; i++) {
    ringdown_phase[i] = true_sign * ringdown_phase[i] + phase_match;
  }

  commondata->nsteps_IMR = idx_match + 1 + nsteps_ringdown;
  commondata->waveform_IMR = (double complex *)malloc(NUMMODES * commondata->nsteps_IMR * sizeof(double complex));
  for (i = 0; i <= idx_match; i++) {
    commondata->waveform_IMR[IDX_WF(i, TIME)] = times_new[i] - commondata->t_attach;
    h22 = h22_amp_new[i] * cexp(I * h22_phase_new[i]);
    commondata->waveform_IMR[IDX_WF(i, STRAIN)] = h22;
  }
  free(h22_amp_new);
  free(h22_phase_new);
  free(times_new);

  for (i = 0; i < nsteps_ringdown; i++) {
    commondata->waveform_IMR[IDX_WF(i + 1 + idx_match, TIME)] = ringdown_time[i] - commondata->t_attach;
    h22 = ringdown_amp[i] * cexp(I * ringdown_phase[i]);
    commondata->waveform_IMR[IDX_WF(i + 1 + idx_match, STRAIN)] = h22;
  }
  free(ringdown_time);
  free(ringdown_phase);
  free(ringdown_amp);
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_IMR_waveform

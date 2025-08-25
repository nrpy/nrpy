#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Calculate the BOB 22 mode.
 */
void BOB_aligned_spin_waveform_from_times(REAL *restrict times, REAL *restrict amps, REAL *restrict phases, const size_t nsteps_BOB,
                                          commondata_struct *restrict commondata) {

  size_t i;
  REAL *restrict wrapped_phases = (REAL *)malloc(nsteps_BOB * sizeof(REAL));
  REAL waveform[2];
  for (i = 0; i < nsteps_BOB; i++) {
    // compute
    BOB_aligned_spin_waveform(times[i], commondata, waveform);
    // store
    amps[i] = waveform[0];
    wrapped_phases[i] = waveform[1];
  }
  // unwrap the phase
  SEOBNRv5_aligned_spin_unwrap(wrapped_phases, phases, nsteps_BOB);

  // shift and take absolute value of the phase
  const REAL pshift = fabs(phases[0]);
  for (i = 0; i < nsteps_BOB; i++) {
    phases[i] = fabs(phases[i]) - pshift;
  }
  free(wrapped_phases);
} // END FUNCTION BOB_aligned_spin_waveform_from_times

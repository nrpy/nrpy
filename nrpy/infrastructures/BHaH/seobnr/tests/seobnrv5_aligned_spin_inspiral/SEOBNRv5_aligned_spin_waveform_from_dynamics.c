#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Calculate the SEOBNRv5 22 mode.
 */
int SEOBNRv5_aligned_spin_waveform_from_dynamics(commondata_struct *restrict commondata) {

  int i;
  REAL dynamics[NUMVARS];
  commondata->waveform_low = (double complex *)malloc(commondata->nsteps_low * NUMMODES * sizeof(double complex));   // t , h_+ , h_x
  commondata->waveform_fine = (double complex *)malloc(commondata->nsteps_fine * NUMMODES * sizeof(double complex)); // t , h_+ , h_x

  // low sampling
  for (i = 0; i < commondata->nsteps_low; i++) {
    // assign
    dynamics[TIME] = commondata->dynamics_low[IDX(i, TIME)];
    dynamics[R] = commondata->dynamics_low[IDX(i, R)];
    dynamics[PHI] = commondata->dynamics_low[IDX(i, PHI)];
    dynamics[PRSTAR] = commondata->dynamics_low[IDX(i, PRSTAR)];
    dynamics[PPHI] = commondata->dynamics_low[IDX(i, PPHI)];
    dynamics[H] = commondata->dynamics_low[IDX(i, H)];
    dynamics[OMEGA] = commondata->dynamics_low[IDX(i, OMEGA)];
    dynamics[OMEGA_CIRC] = commondata->dynamics_low[IDX(i, OMEGA_CIRC)];

    // compute
    // store
    commondata->waveform_low[IDX_WF(i, TIME)] = dynamics[TIME];
    commondata->waveform_low[IDX_WF(i, STRAIN)] = SEOBNRv5_aligned_spin_waveform(dynamics, commondata);
  }
  // high sampling
  for (i = 0; i < commondata->nsteps_fine; i++) {
    // assign
    dynamics[TIME] = commondata->dynamics_fine[IDX(i, TIME)];
    dynamics[R] = commondata->dynamics_fine[IDX(i, R)];
    dynamics[PHI] = commondata->dynamics_fine[IDX(i, PHI)];
    dynamics[PRSTAR] = commondata->dynamics_fine[IDX(i, PRSTAR)];
    dynamics[PPHI] = commondata->dynamics_fine[IDX(i, PPHI)];
    dynamics[H] = commondata->dynamics_fine[IDX(i, H)];
    dynamics[OMEGA] = commondata->dynamics_fine[IDX(i, OMEGA)];
    dynamics[OMEGA_CIRC] = commondata->dynamics_fine[IDX(i, OMEGA_CIRC)];

    // compute
    // store
    commondata->waveform_fine[IDX_WF(i, TIME)] = dynamics[TIME];
    commondata->waveform_fine[IDX_WF(i, STRAIN)] = SEOBNRv5_aligned_spin_waveform(dynamics, commondata);
  }
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_waveform_from_dynamics

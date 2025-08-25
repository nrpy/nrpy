#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

/**
 * Integrate the SEOBNRv5 equations of motion.
 */
int SEOBNRv5_aligned_spin_ode_integration(commondata_struct *restrict commondata) {

  const gsl_odeiv2_step_type *restrict T = gsl_odeiv2_step_rk8pd;
  gsl_odeiv2_step *restrict s = gsl_odeiv2_step_alloc(T, 4);
  gsl_odeiv2_control *restrict c = gsl_odeiv2_control_standard_new(1e-12, 1e-11, 1.0, 1.0);
  gsl_odeiv2_system sys = {SEOBNRv5_aligned_spin_right_hand_sides, NULL, 4, commondata};
  gsl_odeiv2_evolve *restrict e = gsl_odeiv2_evolve_alloc(4);

  REAL t = 0.0;
  const REAL tmax = 2e9;
  REAL y[4], dydt[4];
  int status = 0;
  int i;
  int stop = 0;
  y[0] = commondata->r;
  y[1] = commondata->phi;
  y[2] = commondata->prstar;
  y[3] = commondata->pphi;
  status = SEOBNRv5_aligned_spin_right_hand_sides(t, y, dydt, commondata);
  int rhs_status[1] = {GSL_SUCCESS};
  char rhs_name[] = "gsl_odeiv2_evolve_apply";
  REAL h = 2.0 * M_PI / dydt[1] / 5.0;
  size_t bufferlength = (size_t)(tmax / h); // runs up to 0.01x maximum time (we should not ideally run that long)
  REAL *restrict dynamics_RK = (REAL *)malloc(bufferlength * (NUMVARS) * sizeof(REAL));
  size_t nsteps = 0;

  // store
  dynamics_RK[IDX(nsteps, TIME)] = t;
  for (i = 1; i < 5; i++) {
    dynamics_RK[IDX(nsteps, i)] = y[i - 1];
  }
  SEOBNRv5_aligned_spin_augments(commondata);
  dynamics_RK[IDX(nsteps, H)] = commondata->Hreal;
  dynamics_RK[IDX(nsteps, OMEGA)] = commondata->dHreal_dpphi;
  dynamics_RK[IDX(nsteps, OMEGA_CIRC)] = commondata->Omega_circ;
  nsteps++;
  REAL Omega_previous = commondata->initial_omega;

  while (t < tmax && stop == 0) {
    // integrate
    status = gsl_odeiv2_evolve_apply(e, c, s, &sys, &t, tmax, &h, y);
    handle_gsl_return_status(status, rhs_status, 1, rhs_name);
    commondata->r = y[0];
    commondata->phi = y[1];
    commondata->prstar = y[2];
    commondata->pphi = y[3];
    SEOBNRv5_aligned_spin_augments(commondata);

    // buffercheck
    if (nsteps >= bufferlength) {
      bufferlength = 2 * bufferlength;
      dynamics_RK = (REAL *)realloc(dynamics_RK, bufferlength * (NUMVARS) * sizeof(REAL));
    }

    // update
    // store
    dynamics_RK[IDX(nsteps, TIME)] = t;
    for (i = 1; i < 5; i++) {
      dynamics_RK[IDX(nsteps, i)] = y[i - 1];
    }
    dynamics_RK[IDX(nsteps, H)] = commondata->Hreal;
    dynamics_RK[IDX(nsteps, OMEGA)] = commondata->dHreal_dpphi;
    dynamics_RK[IDX(nsteps, OMEGA_CIRC)] = commondata->Omega_circ;
    nsteps++;

    // stopcheck
    if (commondata->r < 6.0) {
      // decrease in frequency: Omega peak stop = index of omega
      if (commondata->dHreal_dpphi < Omega_previous) {
        stop = OMEGA;
        break;
      }
      SEOBNRv5_aligned_spin_right_hand_sides(t, y, dydt, commondata);
      // outspiral dPR: PR peak stop = index of pr
      if (dydt[2] > 0.0) {
        stop = PRSTAR;
        break;
      }
      // outspiral dR
      if (dydt[0] > 0.0) {
        break;
      }
      // Stopping radius
      if (commondata->r < commondata->r_stop) {
        break;
      }
      // unphysical frequency
      if (commondata->r < 3.0 && commondata->Omega_circ > 1.0) {
        break;
      }
    }
    Omega_previous = commondata->dHreal_dpphi;
  }

  // free up gsl ode solver
  gsl_odeiv2_control_free(c);
  gsl_odeiv2_step_free(s);
  gsl_odeiv2_evolve_free(e);

  // save the times as a separate array
  REAL *restrict times = malloc(nsteps * sizeof(REAL));
  for (i = 0; i < nsteps; i++) {
    times[i] = dynamics_RK[IDX(i, TIME)];
  }

  // find the coarse-fine separation index
  REAL t_desired;
  if (stop == OMEGA) {
    t_desired = times[nsteps - 1] - commondata->t_stepback - 50.;
  } else {
    t_desired = times[nsteps - 1] - commondata->t_stepback;
  }
  size_t coarse_fine_separation_idx = gsl_interp_bsearch(times, t_desired, 0, nsteps - 1);
  commondata->dynamics_low = (REAL *)malloc(NUMVARS * coarse_fine_separation_idx * sizeof(REAL));
  commondata->nsteps_low = coarse_fine_separation_idx;
  for (i = 0; i < commondata->nsteps_low; i++) {
    commondata->dynamics_low[IDX(i, TIME)] = dynamics_RK[IDX(i, TIME)];
    commondata->dynamics_low[IDX(i, R)] = dynamics_RK[IDX(i, R)];
    commondata->dynamics_low[IDX(i, PHI)] = dynamics_RK[IDX(i, PHI)];
    commondata->dynamics_low[IDX(i, PRSTAR)] = dynamics_RK[IDX(i, PRSTAR)];
    commondata->dynamics_low[IDX(i, PPHI)] = dynamics_RK[IDX(i, PPHI)];
    commondata->dynamics_low[IDX(i, H)] = dynamics_RK[IDX(i, H)];
    commondata->dynamics_low[IDX(i, OMEGA)] = dynamics_RK[IDX(i, OMEGA)];
    commondata->dynamics_low[IDX(i, OMEGA_CIRC)] = dynamics_RK[IDX(i, OMEGA_CIRC)];
  }

  size_t nsteps_fine_prelim = nsteps - coarse_fine_separation_idx;
  REAL *restrict dynamics_fine_prelim = (REAL *)malloc(NUMVARS * nsteps_fine_prelim * sizeof(REAL));
  REAL *restrict times_fine_prelim = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  for (i = 0; i < nsteps_fine_prelim; i++) {
    dynamics_fine_prelim[IDX(i, TIME)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, TIME)];
    times_fine_prelim[i] = dynamics_fine_prelim[IDX(i, TIME)];
    dynamics_fine_prelim[IDX(i, R)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, R)];
    dynamics_fine_prelim[IDX(i, PHI)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, PHI)];
    dynamics_fine_prelim[IDX(i, PRSTAR)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, PRSTAR)];
    dynamics_fine_prelim[IDX(i, PPHI)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, PPHI)];
    dynamics_fine_prelim[IDX(i, H)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, H)];
    dynamics_fine_prelim[IDX(i, OMEGA)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, OMEGA)];
    dynamics_fine_prelim[IDX(i, OMEGA_CIRC)] = dynamics_RK[IDX(i + coarse_fine_separation_idx, OMEGA_CIRC)];
  }

  commondata->dynamics_raw = malloc(NUMVARS * nsteps * sizeof(REAL));
  commondata->nsteps_raw = nsteps;
  for (i = 0; i < nsteps; i++) {
    commondata->dynamics_raw[IDX(i, TIME)] = dynamics_RK[IDX(i, TIME)];
    commondata->dynamics_raw[IDX(i, R)] = dynamics_RK[IDX(i, R)];
    commondata->dynamics_raw[IDX(i, PHI)] = dynamics_RK[IDX(i, PHI)];
    commondata->dynamics_raw[IDX(i, PRSTAR)] = dynamics_RK[IDX(i, PRSTAR)];
    commondata->dynamics_raw[IDX(i, PPHI)] = dynamics_RK[IDX(i, PPHI)];
    commondata->dynamics_raw[IDX(i, H)] = dynamics_RK[IDX(i, H)];
    commondata->dynamics_raw[IDX(i, OMEGA)] = dynamics_RK[IDX(i, OMEGA)];
    commondata->dynamics_raw[IDX(i, OMEGA_CIRC)] = dynamics_RK[IDX(i, OMEGA_CIRC)];
  }

  free(dynamics_RK);

  // t_peak = dynamics[-1] if there is no peak

  REAL t_peak = times_fine_prelim[nsteps_fine_prelim - 1];

  // perform iterative refinement to find the true peak of the dynamics
  if (stop == OMEGA) {
    REAL *t_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
    REAL *omega_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
    for (i = 0; i < nsteps_fine_prelim; i++) {
      t_values[i] = dynamics_fine_prelim[IDX(i, TIME)];
      omega_values[i] = dynamics_fine_prelim[IDX(i, OMEGA)];
    }
    REAL dt = times_fine_prelim[nsteps_fine_prelim - 1] - times_fine_prelim[nsteps_fine_prelim - 2];
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
    gsl_spline_init(spline, t_values, omega_values, nsteps_fine_prelim);
    spline_data sdata = {spline, acc};
    t_peak = SEOBNRv5_aligned_spin_iterative_refinement(&sdata, times_fine_prelim[0], times_fine_prelim[nsteps_fine_prelim - 1], 3, dt, false);
    free(t_values);
    free(omega_values);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
  if (stop == PRSTAR) {
    REAL *t_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
    REAL *prstar_values = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
    for (i = 0; i < nsteps_fine_prelim; i++) {
      t_values[i] = dynamics_fine_prelim[IDX(i, TIME)];
      prstar_values[i] = dynamics_fine_prelim[IDX(i, PRSTAR)];
    }
    REAL dt = times_fine_prelim[nsteps_fine_prelim - 1] - times_fine_prelim[nsteps_fine_prelim - 2];
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
    gsl_spline_init(spline, t_values, prstar_values, nsteps_fine_prelim);
    spline_data sdata = {spline, acc};
    t_peak = SEOBNRv5_aligned_spin_iterative_refinement(&sdata, times_fine_prelim[0], times_fine_prelim[nsteps_fine_prelim - 1], 3, dt, true);
    free(t_values);
    free(prstar_values);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }

  // interpolate the dynamics

  status = SEOBNRv5_aligned_spin_interpolate_dynamics(commondata, dynamics_fine_prelim, nsteps_fine_prelim, t_peak, stop);
  free(dynamics_fine_prelim);
  free(times_fine_prelim);
  free(times);
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_ode_integration

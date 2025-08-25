#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate the peak in frequency or momentum in SEOBNRv5
 */
int SEOBNRv5_aligned_spin_interpolate_dynamics(commondata_struct *restrict commondata, REAL *restrict dynamics_fine_prelim,
                                               const size_t nsteps_fine_prelim, const REAL t_peak, const int stop) {

  int i;
  // Intepolate the high sampled dynamics for NQCs.
  REAL time_start = dynamics_fine_prelim[IDX(0, TIME)];
  REAL time_end = dynamics_fine_prelim[IDX(nsteps_fine_prelim - 1, TIME)];
  if (stop != 0) {
    time_start = MAX(t_peak - commondata->t_stepback, dynamics_fine_prelim[IDX(0, TIME)]);
    time_end = MIN(t_peak, time_end);
  }

  REAL *restrict ts = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  REAL *restrict rs = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  REAL *restrict phis = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  REAL *restrict prs = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  REAL *restrict pphis = (REAL *)malloc(nsteps_fine_prelim * sizeof(REAL));
  for (i = 0; i < nsteps_fine_prelim; i++) {
    ts[i] = dynamics_fine_prelim[IDX(i, TIME)];
    rs[i] = dynamics_fine_prelim[IDX(i, R)];
    phis[i] = dynamics_fine_prelim[IDX(i, PHI)];
    prs[i] = dynamics_fine_prelim[IDX(i, PRSTAR)];
    pphis[i] = dynamics_fine_prelim[IDX(i, PPHI)];
  }
  gsl_interp_accel *restrict r_acc = gsl_interp_accel_alloc();
  gsl_spline *restrict r_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  gsl_spline_init(r_spline, ts, rs, nsteps_fine_prelim);
  gsl_interp_accel *restrict phi_acc = gsl_interp_accel_alloc();
  gsl_spline *restrict phi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  gsl_spline_init(phi_spline, ts, phis, nsteps_fine_prelim);
  gsl_interp_accel *restrict pr_acc = gsl_interp_accel_alloc();
  gsl_spline *restrict pr_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  gsl_spline_init(pr_spline, ts, prs, nsteps_fine_prelim);
  gsl_interp_accel *restrict pphi_acc = gsl_interp_accel_alloc();
  gsl_spline *restrict pphi_spline = gsl_spline_alloc(gsl_interp_cspline, nsteps_fine_prelim);
  gsl_spline_init(pphi_spline, ts, pphis, nsteps_fine_prelim);

  const REAL dt = 0.1;
  commondata->nsteps_fine = (size_t)((time_end - time_start) / dt + 1);
  commondata->dynamics_fine = (REAL *)malloc(NUMVARS * commondata->nsteps_fine * sizeof(REAL));
  REAL t;
  for (i = 0; i < commondata->nsteps_fine; i++) {
    t = time_start + i * dt;
    commondata->dynamics_fine[IDX(i, TIME)] = t;
    commondata->dynamics_fine[IDX(i, R)] = gsl_spline_eval(r_spline, t, r_acc);
    commondata->dynamics_fine[IDX(i, PHI)] = gsl_spline_eval(phi_spline, t, phi_acc);
    commondata->dynamics_fine[IDX(i, PRSTAR)] = gsl_spline_eval(pr_spline, t, pr_acc);
    commondata->dynamics_fine[IDX(i, PPHI)] = gsl_spline_eval(pphi_spline, t, pphi_acc);
    commondata->r = commondata->dynamics_fine[IDX(i, R)];
    commondata->phi = commondata->dynamics_fine[IDX(i, PHI)];
    commondata->prstar = commondata->dynamics_fine[IDX(i, PRSTAR)];
    commondata->pphi = commondata->dynamics_fine[IDX(i, PPHI)];
    SEOBNRv5_aligned_spin_augments(commondata);
    commondata->dynamics_fine[IDX(i, H)] = commondata->Hreal;
    commondata->dynamics_fine[IDX(i, OMEGA)] = commondata->dHreal_dpphi;
    commondata->dynamics_fine[IDX(i, OMEGA_CIRC)] = commondata->Omega_circ;
  }

  gsl_spline_free(r_spline);
  gsl_interp_accel_free(r_acc);
  gsl_spline_free(phi_spline);
  gsl_interp_accel_free(phi_acc);
  gsl_spline_free(pr_spline);
  gsl_interp_accel_free(pr_acc);
  gsl_spline_free(pphi_spline);
  gsl_interp_accel_free(pphi_acc);

  free(ts);
  free(rs);
  free(phis);
  free(prs);
  free(pphis);
  return GSL_SUCCESS;
} // END FUNCTION SEOBNRv5_aligned_spin_interpolate_dynamics

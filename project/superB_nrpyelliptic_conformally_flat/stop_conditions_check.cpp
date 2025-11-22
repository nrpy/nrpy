#include "BHaH_defines.h"

/**
 * Evaluate stop conditions.
 */
void stop_conditions_check(commondata_struct *restrict commondata) {
  // Check if total number of iteration steps has been reached:
  //  EITHER the number of pseudo-time steps nn has reached nn_max,
  //  OR the log10(residual integral) < log10(residual integral tolerance).
  if (commondata->nn >= commondata->nn_max) {
    printf("\nMax iterations reached: iter=%d (max=%d), log10(res)=%.6e, log10(tol)=%.6e\n", commondata->nn, commondata->nn_max,
           commondata->log10_current_residual, commondata->log10_residual_tolerance);
    commondata->stop_relaxation = true;
  } else if (commondata->log10_current_residual < commondata->log10_residual_tolerance) {
    printf("\nConvergence reached: iter=%d (max=%d), log10(res)=%.6e, log10(tol)=%.6e\n", commondata->nn, commondata->nn_max,
           commondata->log10_current_residual, commondata->log10_residual_tolerance);
    commondata->stop_relaxation = true;
  } // END IF max iterations or convergence reached.
} // END FUNCTION stop_conditions_check

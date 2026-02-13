#include "BHaH_defines.h"

/**
 * @brief Computes the initial time-component of the 4-velocity (u^0)
 *
 *         Solves the quadratic Hamiltonian constraint equation:
 *             g_munu u^mu u^nu = -1
 *         for the positive root of u^0, given the spatial velocity components.
 *
 *         Input:
 *             metric: The metric tensor components at the current location.
 *             f[8]: The state vector (specifically spatial velocities f[5]..f[7]).
 *         Output:
 *             u0_out: The computed u^0 component.
 */
void u0_massive(const metric_struct *restrict metric, const double f[8], double *restrict u0_out) {
  // Unpack spatial velocity components from f[5]..f[7]
  // Note: f[4] is u^0 (which we are computing), so we skip it.
  const double uU1 = f[5];
  const double uU2 = f[6];
  const double uU3 = f[7];
  const REAL tmp0 = 2 * uU1;
  const REAL tmp2 = metric->g4DD01 * tmp0 + 2 * metric->g4DD02 * uU2 + 2 * metric->g4DD03 * uU3;
  *u0_out = (1.0 / 2.0) *
            (-tmp2 - sqrt(-4 * metric->g4DD00 *
                              (metric->g4DD11 * ((uU1) * (uU1)) + metric->g4DD12 * tmp0 * uU2 + metric->g4DD13 * tmp0 * uU3 +
                               metric->g4DD22 * ((uU2) * (uU2)) + 2 * metric->g4DD23 * uU2 * uU3 + metric->g4DD33 * ((uU3) * (uU3)) + 1) +
                          ((tmp2) * (tmp2)))) /
            metric->g4DD00;
} // END FUNCTION u0_massive

#include "BHaH_defines.h"

/**
 * @brief Computes the initial time-component of the 4-momentum (p^0).
 *
 *         Solves the quadratic Hamiltonian constraint equation:
 *             g_munu p^mu p^nu = 0
 *         for the negative root of p^0, given the spatial momentum components.
 *
 *         Input:
 *             metric: The metric tensor components at the current location.
 *             f[9]: The state vector (specifically spatial momentum f[5]..f[7]).
 *         Output:
 *             p0_out: The computed p^0 component.
 */
void p0_reverse(const metric_struct *restrict metric, const double f[9], double *restrict p0_out) {
  // Unpack spatial momentum components from f[5]..f[7]
  // Note: f[4] is p^0 (which we are computing), so we skip it.
  const double pU1 = f[5];
  const double pU2 = f[6];
  const double pU3 = f[7];
  const REAL tmp0 = 2 * pU1;
  const REAL tmp2 = metric->g4DD01 * tmp0 + 2 * metric->g4DD02 * pU2 + 2 * metric->g4DD03 * pU3;
  *p0_out = (1.0 / 2.0) *
            (-tmp2 + sqrt(-4 * metric->g4DD00 *
                              (metric->g4DD11 * ((pU1) * (pU1)) + metric->g4DD12 * pU2 * tmp0 + metric->g4DD13 * pU3 * tmp0 +
                               metric->g4DD22 * ((pU2) * (pU2)) + 2 * metric->g4DD23 * pU2 * pU3 + metric->g4DD33 * ((pU3) * (pU3))) +
                          ((tmp2) * (tmp2)))) /
            metric->g4DD00;
} // END FUNCTION p0_reverse

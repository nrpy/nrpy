#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * C function to perform numpy.unwrap.
 */
void SEOBNRv5_aligned_spin_unwrap(REAL *restrict angles_in, REAL *restrict angles_out, size_t nsteps_arr) {

  angles_out[0] = angles_in[0];
  REAL diff;
  for (size_t i = 1; i < nsteps_arr; i++) {
    diff = angles_in[i] - angles_in[i - 1];
    diff = fabs(diff) > M_PI ? (diff < -M_PI ? diff + 2 * M_PI : diff - 2 * M_PI) : diff;
    angles_out[i] = angles_out[i - 1] + diff;
  }
} // END FUNCTION SEOBNRv5_aligned_spin_unwrap

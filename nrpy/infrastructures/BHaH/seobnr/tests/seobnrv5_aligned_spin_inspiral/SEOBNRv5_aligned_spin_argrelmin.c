#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * C function to perform scipy.argrelmin for order = 3.
 */
size_t SEOBNRv5_aligned_spin_argrelmin(REAL *restrict arr, size_t nsteps_arr) {

  size_t order = 3;
  size_t minima_count_or_idx = 0;
  // Loop through the array with bounds of `order` on each side
  for (size_t i = order; i < nsteps_arr - order; i++) {
    int isMinimum = 1;

    // Check `order` elements to the left and right
    for (int j = 1; j <= order; j++) {
      if (arr[i] >= arr[i - j] || arr[i] >= arr[i + j]) {
        isMinimum = 0;
        break;
      }
    }
    if (isMinimum) {
      minima_count_or_idx = i;
      break;
    }
  }
  return minima_count_or_idx;
} // END FUNCTION SEOBNRv5_aligned_spin_argrelmin

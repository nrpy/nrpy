#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Find the local minimum index in an array.
 */
size_t find_local_minimum_index(REAL *restrict arr, size_t size, int order) {

  if (size < 2 * order + 1) {
    return -1; // Not enough points to apply the order
  }

  for (size_t i = 0; i < size; ++i) {
    bool is_min = true;
    for (int shift = 1; shift <= order; ++shift) {
      // clipped indexing: return 0 or size-1 if out of bounds
      size_t left_idx = MAX(i, shift) - shift;     // returns 0 if i < shift else i - shift
      size_t right_idx = MIN(i + shift, size - 1); // returns size-1 if i > size-1 else i + shift
      is_min = is_min && (arr[i] < arr[left_idx]) && (arr[i] < arr[right_idx]);
      if (is_min)
        return i;
    }
  }
  return -1; // No local minimum found
} // END FUNCTION find_local_minimum_index

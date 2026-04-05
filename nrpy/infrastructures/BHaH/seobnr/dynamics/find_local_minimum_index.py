"""
Register CFunction for finding the local minimum index in an array.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
        Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_find_local_minimum_index() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for finding the local minimum index in an array.
    Needed by SEOBNRv5_aligned_spin_iterative_refinement.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Finds the local minimum index in an array. 
Implements the logic of scipy.argrelmin with order = 3 and clipped indexing.
However, it returns the index of the first minimum from the left, 
as required by SEOBNRv5_aligned_spin_iterative_refinement.

@param arr - The array to search for the minimum in.
@param size - The size of the array.
@param order - The maximum array index shift to consider for the minimum.
@returns - The index of the local minimum in the array or -1 if no minimum is found.
"""
    cfunc_type = "size_t"
    name = "find_local_minimum_index"
    params = "REAL *restrict arr, size_t size, int order"
    body = """
  if (size < 2 * order + 1) {
    return -1; // Not enough points to apply the order
  }

  for (size_t i = 0; i < size; ++i) {
    bool is_min = true;
    for (int shift = 1; shift <= order; ++shift) {
      // clipped indexing: return 0 or size-1 if out of bounds
      size_t left_idx = MAX(i,shift) - shift; // returns 0 if i < shift else i - shift
      size_t right_idx = MIN(i+shift,size-1); // returns size-1 if i > size-1 else i + shift
      is_min = is_min && (arr[i] < arr[left_idx]) && (arr[i] < arr[right_idx]);
      if (is_min) return i;
    }
  }
  return -1; // No local minimum found
"""
    cfc.register_CFunction(
        subdirectory="dynamics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

"""
Register CFunction that performs scipy.argrelmin with order = 3.

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


def register_CFunction_SEOBNRv5_aligned_spin_argrelmin() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction that performs scipy.argrelmin with order = 3.
    Needed by SEOBNRv5_aligned_spin_iterative_refinement.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
C function to calculate the index of the first minimum in an array.
The logic follows that of scipy.argrelmin with order = 3,
however it outputs the index of the first relative minimum from the left.

@param arr - The array to find the minimum in.
@param nsteps_arr - The number of steps in the array.
@returns - The index of the first minimum in the array.
"""
    cfunc_type = "size_t"
    name = "SEOBNRv5_aligned_spin_argrelmin"
    params = "REAL *restrict arr , size_t nsteps_arr"
    body = """
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

"""
Register CFunction for evaluating the peak in frequency or momentum in SEOBNRv5.

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


def register_CFunction_SEOBNRv5_aligned_spin_iterative_refinement() -> (
    Union[None, pcg.NRPyEnv_type]
):
    """
    Register CFunction for evaluating the peak in frequency or momentum in SEOBNRv5.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
Evaluates the time at which the peak in orbital frequency or tortoise momentum occurs in the SEOBNRv5 dynamics.
This is done by iteratively interpolating and refining the interval in which the peak occurs.

@param sdata - Struct containing the spline of the tortoise momentum or orbital frequency.
@param initial_left - The left endpoint of the interval to search for the peak.
@param initial_right - The right endpoint of the interval to search for the peak.
@param levels - The number of iterations of refinement to perform.
@param dt_initial - The initial time step to use for interpolation.
@param pr - Boolean flag to determine whether to find the peak in tortoise momentum or orbital frequency.
@returns - The time at which the peak in orbital frequency or tortoise momentum occurs, or the midpoint of the interval if no peak is found.
"""
    cfunc_type = "REAL"
    name = "SEOBNRv5_aligned_spin_iterative_refinement"
    params = "spline_data *sdata, double initial_left, double initial_right, int levels, double dt_initial, bool pr"
    body = """
if (levels < 2) {
  printf("Error: levels must be greater than 1\\n");
  return (initial_left + initial_right) / 2.0;
}

double current_left = initial_left;
double current_right = initial_right;
double result = (current_left + current_right) / 2.0; // Default if nothing found
bool found_any_min = false;

for (int n = 1; n <= levels; ++n) {
  double dt = dt_initial / pow(10.0, n);
    
  // Generate t_fine points dynamically, similar to numpy.arange
  // Ensure at least 2 points for interval, more for derivative
  size_t num_t_fine = (size_t)floor((current_right - current_left) / dt) + 1;
  if (num_t_fine < 5) { // Need enough points for argrelmin order=3
    if (pr) return current_right;
    return (current_left + current_right) / 2.0;
  }

  double *t_fine = (double *)malloc(num_t_fine * sizeof(double));
  if (t_fine == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_iterative_refinement(), malloc() failed to for t_fine\\n");
    exit(1);
  }
  double *deriv_values = (double *)malloc(num_t_fine * sizeof(double));
  if (deriv_values == NULL) {
    fprintf(stderr,"Error: in SEOBNRv5_aligned_spin_iterative_refinement(), malloc() failed to for deriv_values\\n");
    exit(1);
  }

  for (size_t i = 0; i < num_t_fine; ++i) {
    t_fine[i] = current_left + i * dt;
    if (t_fine[i] > current_right) t_fine[i] = current_right; // Cap at right boundary
    deriv_values[i] = fabs(gsl_spline_eval_deriv(sdata->spline, t_fine[i], sdata->acc));
  }

  // Find local minimum using a simplified argrelmin logic
  int min_idx = find_local_minimum_index(deriv_values, num_t_fine, 3); // order=3 from Python

  if (min_idx != -1) {
    result = t_fine[min_idx];
    
    // Refine interval
    current_left = fmax(result - 10.0 * dt, initial_left);
    current_right = fmin(result + 10.0 * dt, initial_right);
    found_any_min = true;
  } else {
    // No local minimum found in this iteration
    free(t_fine);
    free(deriv_values);
    if (pr) return current_right;
    return (current_left + current_right) / 2.0;
  }
    
  free(t_fine);
  free(deriv_values);
}

if (!found_any_min) {
  if (pr) return initial_right;
  return (initial_left + initial_right) / 2.0;
}
return result;
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

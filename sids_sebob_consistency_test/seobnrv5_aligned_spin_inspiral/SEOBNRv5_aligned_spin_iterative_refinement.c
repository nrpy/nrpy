#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate the peak in frequency or momentum in SEOBNRv5
 */
REAL SEOBNRv5_aligned_spin_iterative_refinement(spline_data *sdata, double initial_left, double initial_right, int levels, double dt_initial,
                                                bool pr) {

  if (levels < 2) {
    printf("Error: levels must be greater than 1\n");
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
      if (pr)
        return current_right;
      return (current_left + current_right) / 2.0;
    }

    double *t_fine = (double *)malloc(num_t_fine * sizeof(double));
    double *deriv_values = (double *)malloc(num_t_fine * sizeof(double));

    if (!t_fine || !deriv_values) {
      printf("Memory allocation failed for t_fine or deriv_values.\n");
      if (t_fine)
        free(t_fine);
      if (deriv_values)
        free(deriv_values);
      return result; // Return last known good result or default
    }

    for (size_t i = 0; i < num_t_fine; ++i) {
      t_fine[i] = current_left + i * dt;
      if (t_fine[i] > current_right)
        t_fine[i] = current_right; // Cap at right boundary
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
      if (pr)
        return current_right;
      return (current_left + current_right) / 2.0;
    }

    free(t_fine);
    free(deriv_values);
  }

  if (!found_any_min) {
    if (pr)
      return initial_right;
    return (initial_left + initial_right) / 2.0;
  }
  return result;
} // END FUNCTION SEOBNRv5_aligned_spin_iterative_refinement

#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
/**
 * Evaluate the absolute derivative at a given point.
 */
double eval_abs_deriv(double t, void *params) {

  spline_data *sdata = (spline_data *)params;
  return fabs(gsl_spline_eval_deriv(sdata->spline, t, sdata->acc));
} // END FUNCTION eval_abs_deriv

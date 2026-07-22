#ifndef __SPLINE_STRUCT_H__
#define __SPLINE_STRUCT_H__

// Provide the shared spline-data container used by SEOBNR interpolation helpers.
// Translation units including this header must include the relevant GSL spline
// headers first; this header intentionally does not add those includes itself.

typedef struct {
  // Cubic spline allocated with gsl_spline_alloc and released with gsl_spline_free.
  gsl_spline *spline;

  // Interpolation accelerator allocated with gsl_interp_accel_alloc and released
  // with gsl_interp_accel_free.
  gsl_interp_accel *acc;
} spline_data;

#endif // __SPLINE_STRUCT_H__

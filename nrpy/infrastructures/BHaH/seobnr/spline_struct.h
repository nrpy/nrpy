#ifndef __SPLINE_STRUCT_H__
#define __SPLINE_STRUCT_H__

typedef struct {
  gsl_spline *spline;
  gsl_interp_accel *acc;
} spline_data;

#endif // __SPLINE_STRUCT_H__
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "gsl/gsl_errno.h"

/**
 * @brief GSL-compatible wrapper for photon geodesics in KerrSchild_Cartesian.
 *
 *         Unpacks the GSL 'params' void pointer into the BHaH 'commondata' struct,
 *         computes the local metric and connections, and calls the RHS calculation routine.
 *
 *         Input:
 *             t: Current value for affine parameter
 *             y[9]: Current state vector.
 *             params: Pointer to commondata_struct.
 *         Output:
 *             f[9]: Computed derivatives (RHS).
 */
int ode_gsl_wrapper_KerrSchild_Cartesian(double t, const double y[9], double f[9], void *params) {
  (void)t; // Mark affine parameter 't' as unused to avoid compiler warnings.

  // 1. Unpack parameters
  commondata_struct *commondata = (commondata_struct *)params;

  // 2. Declare geometric structs to hold intermediate results
  metric_struct metric;
  connection_struct conn;

  // 4. Compute Metric
  // Signature: (commondata, y, &metric)
  g4DD_metric_KerrSchild_Cartesian(commondata, y, &metric);

  // 3. Compute Connections (Christoffel Symbols)
  // Signature: (commondata, y, &conn)
  connections_KerrSchild_Cartesian(commondata, y, &conn);

  // 4. Compute Geodesic RHS
  // Signature: (y, metric, conn, rhs_out)
  calculate_ode_rhs(y, &metric, &conn, f);

  return GSL_SUCCESS;

} // END FUNCTION ode_gsl_wrapper_KerrSchild_Cartesian

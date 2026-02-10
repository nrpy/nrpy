#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "gsl/gsl_errno.h"

/**
 * @brief GSL-compatible wrapper for massive particle geodesics in KerrSchild_Cartesian.
 *
 * Unpacks the GSL 'params' void pointer into the BHaH 'commondata' struct, computes the local metric and connections, and calls the RHS calculation
 * routine.
 *
 * Input:
 *     t: Current proper time (unused in autonomous systems).
 *     y[8]: Current state vector.
 *     params: Pointer to commondata_struct.
 * Output:
 *     f[8]: Computed derivatives (RHS).
 */
int ode_gsl_wrapper_massive_KerrSchild_Cartesian(double t, const double y[8], double f[8], void *params) {
  (void)t; // Mark proper time 't' as unused to avoid compiler warnings.

  // 1. Unpack parameters
  commondata_struct *commondata = (commondata_struct *)params;

  // 2. Declare geometric structs to hold intermediate results
  connection_struct conn;

  // 3. Compute Connections (Christoffel Symbols)
  // Signature: (commondata, y, &conn)
  connections_KerrSchild_Cartesian(commondata, y, &conn);

  // 4. Compute Geodesic RHS
  // Signature: (y, &conn, f)
  calculate_ode_rhs_massive(y, &conn, f);

  return GSL_SUCCESS;

} // END FUNCTION ode_gsl_wrapper_massive_KerrSchild_Cartesian

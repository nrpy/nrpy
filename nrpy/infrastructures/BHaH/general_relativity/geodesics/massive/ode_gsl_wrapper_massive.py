"""
Register C wrapper function to interface with the GSL ODE solver.

This module registers the 'ode_gsl_wrapper_massive_{spacetime}' C function.
It acts as a dispatcher, bridge the GSL solver to the NRPy-generated SoA
(Structure of Arrays) kernels by treating a single GSL call as a batch of size 1.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def ode_gsl_wrapper_massive(spacetime_name: str) -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our SoA-compatible kernels with
    a batch size of 1, and returning the computed derivatives.

    :param spacetime_name: String used to identify the metric and connection functions.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]

    desc = f"""@brief GSL-compatible wrapper for massive particle geodesics in {spacetime_name}.

    Unpacks the GSL 'params' void pointer into the BHaH 'commondata' struct,
    computes the local metric and Christoffel symbols (connections) using
    SoA kernels (batch size 1), and calls the RHS calculation routine.

    Input:
        t: Current proper time (unused in autonomous systems).
        y[8]: Current state vector [t, x, y, z, u^t, u^x, u^y, u^z].
        params: Pointer to commondata_struct.
    Output:
        f[8]: Computed derivatives (RHS)."""

    cfunc_type = "int"
    name = f"ode_gsl_wrapper_massive_{spacetime_name}"
    params = "double t, const double y[8], double f[8], void *params"

    # Construct the body using the SoA-compatible signatures
    # We use batch_size=1 and batch_id=0 to adapt the SoA logic to a single ray
    body = f"""
    (void)t; // Mark proper time 't' as unused to avoid compiler warnings.

    // 1. Unpack parameters
    commondata_struct *commondata = (commondata_struct *)params;

    // 2. Declare geometric buffers to hold intermediate results (SoA style)
    // metric_local: 10 components (symmetric 4x4)
    // conn_local: 40 components (Gamma^alpha_{{beta gamma}}, symmetric in lower indices)
    double metric_local[10];
    double conn_local[40];

    // 3. Compute Metric and Connections
    // Signature: (commondata, y, out_array, batch_size, batch_id)
    g4DD_metric_{spacetime_name}(commondata, y, metric_local, 1, 0);
    connections_{spacetime_name}(commondata, y, conn_local, 1, 0);

    // 4. Compute Geodesic RHS
    // We map GSL's output array 'f' directly as the 'k_array' for RK stage 1.
    // Signature: (f_temp, metric, conn, k_array, bundle_capacity, stage, j)
    calculate_ode_rhs_massive(y, metric_local, conn_local, f, 1, 1, 0);

    return GSL_SUCCESS;
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
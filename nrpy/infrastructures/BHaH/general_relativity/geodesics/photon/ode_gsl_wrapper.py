"""
Register C wrapper function to interface with the GSL ODE solver.

This module registers the 'ode_gsl_wrapper_{spacetime}' C function.
It acts as a dispatcher, unpacking the generic `void *params` pointer into
NRPy-specific structures (commondata), computing the metric and Christoffel
symbols, and invoking the RHS calculation engine.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def ode_gsl_wrapper(spacetime_name: str) -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our project-specific dispatchers and
    engines, and returning the computed derivatives to the GSL solver.

    :param spacetime_name: String used to define metric in analytic_spacetimes.py.

    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]
    desc = f"""@brief GSL-compatible wrapper for photon geodesics in {spacetime_name}.

        Unpacks the GSL 'params' void pointer into the BHaH 'commondata' struct,
        computes the local metric and connections, and calls the RHS calculation routine.

        Input:
            t: Current value for affine parameter
            y[9]: Current state vector.
            params: Pointer to commondata_struct.
        Output:
            f[9]: Computed derivatives (RHS)."""

    cfunc_type = "int"
    name = f"ode_gsl_wrapper_{spacetime_name}"
    params = "double t, const double y[9], double f[9], void *params"

    # Construct the body with specific function calls utilizing the flat array SoA architecture
    body = f"""
    (void)t; // Mark affine parameter 't' as unused to avoid compiler warnings.

    // 1. Unpack parameters
    commondata_struct *commondata = (commondata_struct *)params;

    // 2. Declare flat geometric arrays to hold intermediate results
    double metric_g4DD[10];
    double conn_GammaUDD[40];

    // 3. Compute Metric for a single particle (batch_size=1, batch_id=0)
    g4DD_metric_{spacetime_name}(commondata, y, metric_g4DD, 1, 0);

    // 4. Compute Connections for a single particle
    connections_{spacetime_name}(commondata, y, conn_GammaUDD, 1, 0);

    // 5. Compute Geodesic RHS directly into the GSL f-vector
    // f serves as k_array. By setting stage=1, (stage-1)*9 evaluates to 0, matching the f index.
    calculate_ode_rhs(y, metric_g4DD, conn_GammaUDD, f, 1, 1, 0);

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

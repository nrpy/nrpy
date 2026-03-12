"""
Define C wrapper function to interface with the GSL ODE solver.

This module defines the 'ode_gsl_wrapper_massive_{spacetime}' C function.
It acts as a CPU-bound dispatcher, bridging the GSL solver to the NRPy-generated
thread-local kernels by treating a single GSL call as a localized evaluation.
Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc


def ode_gsl_wrapper_massive(spacetime_name: str) -> None:
    """
    Define the C wrapper function to interface with the GSL ODE solver.

    This function unpacks the generic void pointer to access simulation data,
    calls the thread-local geometric kernels, and returns the computed derivatives.

    :param spacetime_name: String used to identify the metric and connection functions.
    """
    # Define C-function arguments and metadata in Master Order
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]

    desc = rf"""@brief GSL-compatible wrapper for massive particle geodesics in {spacetime_name}.

    Unpacks the GSL parameters void pointer into the global commondata struct,
    computes the local metric $g_{{\mu\nu}}$ and Christoffel symbols $\Gamma^\alpha_{{\mu\nu}}$
    using thread-local kernels, and calls the RHS calculation routine.

    @param t Current proper time $\tau$ (unused in autonomous systems).
    @param y Current state vector $[t, x, y, z, u^t, u^x, u^y, u^z]$.
    @param f Computed derivatives (RHS) $dx^\mu/d\tau$ and $du^\mu/d\tau$.
    @param params Pointer to commondata_struct."""

    cfunc_type = "int"

    name = f"ode_gsl_wrapper_massive_{spacetime_name}"

    params = "double t, const double y[8], double f[8], void *params"

    body = rf"""  // --- GSL TO THREAD-LOCAL BRIDGE ---
    // Algorithmic Step: Cast the generic parameter pointer and prepare local buffers.
    // Hardware Justification: Thread-local execution ensures that the ODE integration remains in L1 cache.
    (void)t; // Mark proper time $\tau$ as unused to prevent compiler warnings.

    commondata_struct *commondata = (commondata_struct *)params; // Cast GSL void pointer to the global parameter struct.

    double metric_local[10]; // Buffer for the 10 unique symmetric metric components $g_{{\mu\nu}}$.
    double conn_local[40]; // Buffer for the 40 unique Christoffel symbols $\Gamma^\alpha_{{\mu\nu}}$.

    // --- GEOMETRIC EVALUATIONS ---
    // Algorithmic Step: Evaluate $g_{{\mu\nu}}$ and $\Gamma^\alpha_{{\mu\nu}}$ at the current coordinate position.
    // Hardware Justification: Inline geometric evaluations prevent the need for pre-computed global memory grids.
    g4DD_metric_{spacetime_name}(commondata, y, metric_local); // Evaluate the spacetime metric $g_{{\mu\nu}}$.
    connections_{spacetime_name}(commondata, y, conn_local); // Evaluate the connections $\Gamma^\alpha_{{\mu\nu}}$.

    // --- RIGHT-HAND SIDE COMPUTATION ---
    // Algorithmic Step: Compute $dx^\mu/d\tau$ and $du^\mu/d\tau$.
    // Hardware Justification: Reusing the GSL output array $f$ as the target buffer avoids an intermediate memory copy.
    calculate_ode_rhs_massive(y, metric_local, conn_local, f); // Evaluate the geodesic equations.

    return GSL_SUCCESS; // Signal successful evaluation to the GSL solver.
    """

    # Register the function, omitting empty kwargs to prevent bloat
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

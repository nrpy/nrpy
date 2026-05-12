# nrpy/infrastructures/BHaH/general_relativity/geodesics/massive/ode_gsl_wrapper_massive.py
r"""
Provides the C wrapper function generation for interfacing with the GSL ODE solver.

This module defines the Python function that constructs and registers the C wrapper
for massive particle geodesics. The generated C function acts as a dispatcher,
bridging the GSL solver to the NRPy-generated local kernels by treating a single GSL
call as a localized evaluation. Localized execution organizes memory access during
ODE integration. Inline geometric evaluations prevent the need for pre-computed
global memory grids. Reusing the GSL output array $f$ as the target buffer avoids an
intermediate memory copy.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def ode_gsl_wrapper_massive(spacetime_name: str) -> None:
    r"""
    Define the Python logic to register the C wrapper for the GSL ODE solver.

    The generated C function unpacks the void pointer to access simulation data, calls
    the thread-local geometric kernels, and returns the computed derivatives to the
    GSL solver.

    :param spacetime_name: String used to identify the metric and connection functions.
    """
    # Define C-function arguments and metadata
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]

    desc = rf""" GSL-compatible wrapper for massive particle geodesics in {spacetime_name}.

    Unpacks the GSL parameters void pointer into the global commondata struct,
    computes the local metric $g_{{\mu\nu}}$ and Christoffel symbols $\Gamma^\alpha_{{\mu\nu}}$
    using thread-local kernels, and calls the RHS calculation routine.
    Thread-local execution organizes memory access during ODE integration. Inline
    geometric evaluations prevent the need for pre-computed global memory grids.
    Reusing the GSL output array $f$ as the target buffer avoids an intermediate
    memory copy.

    @param t Current proper time $\tau$ (unused in autonomous systems).
    @param y Current state vector $[t, x, y, z, u^t, u^x, u^y, u^z]$.
    @param f Computed derivatives (RHS) $dx^\mu/d\tau$ and $du^\mu/d\tau$.
    @param params Pointer to commondata_struct."""

    cfunc_type = "int"

    name = f"ode_gsl_wrapper_massive_{spacetime_name}"

    params = "double t, const double y[8], double f[8], void *params"

    body = rf"""  //==========================================
    // GSL TO THREAD-LOCAL BRIDGE
    //==========================================
    // Cast the generic parameter pointer and prepare local buffers.
    (void)t; // Mark proper time $\tau$ as unused to prevent compiler warnings.

    commondata_struct *commondata = (commondata_struct *)params; // Cast GSL void pointer to the global parameter struct.

    double metric_local[10]; // Buffer for the 10 unique symmetric metric components $g_{{\mu\nu}}$.
    double conn_local[40]; // Buffer for the 40 unique Christoffel symbols $\Gamma^\alpha_{{\mu\nu}}$.

    //==========================================
    // GEOMETRIC EVALUATIONS
    //==========================================
    // Evaluate $g_{{\mu\nu}}$ and $\Gamma^\alpha_{{\mu\nu}}$ at the current coordinate position.
    g4DD_metric_{spacetime_name}(commondata, y, metric_local); // Evaluate the spacetime metric $g_{{\mu\nu}}$.
    connections_{spacetime_name}(commondata, y, conn_local); // Evaluate the connections $\Gamma^\alpha_{{\mu\nu}}$.

    //==========================================
    // RIGHT-HAND SIDE COMPUTATION
    //==========================================
    // Compute $dx^\mu/d\tau$ and $du^\mu/d\tau$.
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

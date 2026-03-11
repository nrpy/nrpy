"""
Register C dispatcher function for the internal RKF45 ODE solver.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def ode_wrapper(spacetime_name: str) -> None:
    """
    Register a C function to dispatch ODE RHS evaluations for a specific spacetime.

    This function uses the NRPy+ infrastructure to generate and register a C
    wrapper. The wrapper handles the local computation of the 4-metric and
    Christoffel symbols for a given spacetime, populates Structure-of-Array
    (SoA) buffers, and triggers the core ODE RHS calculation.

    :param spacetime_name: The name of the spacetime (e.g., 'KerrSchild_Cartesian'). This string is used to suffix the C function names for the metric and connection lookups.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"""@brief Internal dispatcher for photon geodesics in {spacetime_name}.

        Computes the local metric and connections directly into local flat arrays,
        writes them to the SoA batch buffers, then calls the RHS calculation routine.

        Input:
            f_in[9]: Current state vector.
            commondata: Pointer to the simulation's common data parameters.
        Output:
            k_array: Computed derivatives (RHS)."""

    cfunc_type = "void"
    name = f"ode_wrapper_{spacetime_name}"
    params = (
        "const double *restrict f_in, "
        "const commondata_struct *restrict commondata, "
        "double *restrict metric_g4DD, "
        "double *restrict conn_GammaUDD, "
        "double *restrict k_array, "
        "int batch_size, "
        "int batch_id"
    )

    body = f"""
    double local_metric[10];
    double local_conn[40];

    // Compute metric and connections directly into local flat arrays
    g4DD_metric_{spacetime_name}(commondata, f_in, local_metric);
    connections_{spacetime_name}(commondata, f_in, local_conn);

    // Populate the 1D batch arrays explicitly using the IDX_LOCAL macro
    for (int i = 0; i < 10; i++) {{
        metric_g4DD[IDX_LOCAL(i, batch_id, batch_size)] = local_metric[i];
    }}
    for (int i = 0; i < 40; i++) {{
        conn_GammaUDD[IDX_LOCAL(i, batch_id, batch_size)] = local_conn[i];
    }}

    // Call the RHS function using the populated batch arrays
    calculate_ode_rhs(f_in, metric_g4DD, conn_GammaUDD, k_array, batch_size, batch_id);
    """

    # Add OpenMP pragmas outside the function scope to ensure dual-architecture compilation
    prefunc = """
    #ifdef USE_GPU
    #pragma omp declare target
    #endif
    """

    postfunc = """
    #ifdef USE_GPU
    #pragma omp end declare target
    #endif
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        prefunc=prefunc,
        postfunc=postfunc,
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

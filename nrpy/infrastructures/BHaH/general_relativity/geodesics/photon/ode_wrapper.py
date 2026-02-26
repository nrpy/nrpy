"""
Register C dispatcher function for the internal RKF45 ODE solver.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc

def ode_wrapper(spacetime_name: str) -> None:
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
        postfunc=postfunc
    )

if __name__ == "__main__":
    import logging
    import os
    import sys
    import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

    sys.path.append(os.getcwd())
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestODEWrapper")

    SPACETIME = "KerrSchild_Cartesian"
    logger.info("Test: Generating Internal Wrapper C-code for %s...", SPACETIME)

    try:
        ode_wrapper(SPACETIME)
        func_name = f"ode_wrapper_{SPACETIME}"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(" -> Success! Wrote %s", filename)
            Bdefines_h.output_BHaH_defines_h(project_dir=".")
        else:
            raise RuntimeError(f"Function {func_name} not registered.")
    except (RuntimeError, OSError) as e:
        logger.error("Test failed: %s", e)
        sys.exit(1)

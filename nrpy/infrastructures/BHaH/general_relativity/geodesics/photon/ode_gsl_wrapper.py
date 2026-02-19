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

    # Construct the body with specific function calls
    body = f"""
    (void)t; // Mark affine parameter 't' as unused to avoid compiler warnings.
    
    // 1. Unpack parameters
    commondata_struct *commondata = (commondata_struct *)params;
    
    // 2. Declare geometric structs to hold intermediate results
    metric_struct metric;
    connection_struct conn;

    // 3. Compute Metric 
    // Signature: (commondata, y, &metric)
    g4DD_metric_{spacetime_name}(commondata, y, &metric);

    
    // 3. Compute Connections (Christoffel Symbols)
    // Signature: (commondata, y, &conn)
    connections_{spacetime_name}(commondata, y, &conn);
    
    // 5. Compute Geodesic RHS
    // Signature: (y, metric, conn, rhs_out)
    calculate_ode_rhs(y, &metric, &conn, f);
    
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
    import logging
    import os
    import sys

    import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestGSLWrapper")

    SPACETIME = "KerrSchild_Cartesian"

    logger.info("Test: Generating GSL Wrapper C-code for %s...", SPACETIME)

    try:
        # 1. Run the Generator
        ode_gsl_wrapper(SPACETIME)

        # 2. Output
        func_name = f"ode_gsl_wrapper_{SPACETIME}"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(" -> Success! Wrote %s", filename)

            # Also output defines to check struct registration
            Bdefines_h.output_BHaH_defines_h(project_dir=".")
            logger.info(" -> Updated BHaH_defines.h")
        else:
            raise RuntimeError(f"Function {func_name} not registered.")

    except (RuntimeError, OSError) as e:
        logger.error("Test failed: %s", e)
        sys.exit(1)

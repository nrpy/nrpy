"""
Generate the C wrapper function to interface with the GSL ODE solver.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.generic import clang_format
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h


def ode_gsl_wrapper_massive(spacetime_name: str) -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our project-specific dispatchers and
    engines, and returning the computed derivatives to the GSL solver.

    It also registers the `gsl_params` struct definition to BHaH_defines.h.

    :param spacetime_name: Name of the spacetime used to construct function calls.
    """
    # Register the gsl_params struct definition
    gsl_params_def = r"""
typedef struct __gsl_params__ {
    commondata_struct *restrict commondata;
    params_struct *restrict params;
    metric_params *restrict metric;
} gsl_params;
"""
    Bdefines_h.register_BHaH_defines("after_general", gsl_params_def)

    name = f"ode_gsl_wrapper_massive_{spacetime_name}"
    cfunc_type = "int"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]
    desc = f"GSL wrapper for the massive particle geodesic ODEs in {spacetime_name}."
    params = "double t, const double y[8], double f[8], void *params"

    # Construct the body with specific function calls
    body = f"""
    (void)t; // Mark proper time 't' as unused to avoid compiler warnings.
    
    // 1. Unpack generic parameters
    gsl_params *gsl_parameters = (gsl_params *)params;
    
    // 2. Declare geometric structs to hold intermediate results
    metric_struct g4DD;
    connection_struct conn;
    
    // 3. Compute Metric
    // Signature: (commondata, params, y, &g4DD)
    g4DD_metric_{spacetime_name}(gsl_parameters->commondata, gsl_parameters->params, y, &g4DD);
    
    // 4. Compute Connections (Christoffel Symbols)
    // Signature: (commondata, params, y, &conn)
    connections_{spacetime_name}(gsl_parameters->commondata, gsl_parameters->params, y, &conn);
    
    // 5. Compute Geodesic RHS
    // Signature: (y, &conn, f)
    calculate_ode_rhs_massive(y, &conn, f);
    
    return GSL_SUCCESS;
    """

    cfc.register_CFunction(
        name=name,
        cfunc_type=cfunc_type,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
    )


if __name__ == "__main__":
    import logging
    import sys
    import os

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestGSLWrapper")

    SPACETIME = "KerrSchild_Cartesian"

    logger.info("Test: Generating GSL Wrapper C-code for %s...", SPACETIME)

    try:
        # 1. Run the Generator
        ode_gsl_wrapper_massive(SPACETIME)

        # 2. Output
        func_name = f"ode_gsl_wrapper_massive_{SPACETIME}"
        if func_name in cfc.CFunction_dict:
            filename = f"{func_name}.c"
            with open(filename, "w", encoding="utf-8") as f:
                f.write(cfc.CFunction_dict[func_name].full_function)
            logger.info(f" -> Success! Wrote {filename}")

            # Also output defines to check struct registration
            Bdefines_h.output_BHaH_defines_h(project_dir=".")
            logger.info(" -> Updated BHaH_defines.h")
        else:
            raise RuntimeError(f"Function {func_name} not registered.")

    except Exception as e:
        logger.error(f"Test failed: {e}")
        sys.exit(1)

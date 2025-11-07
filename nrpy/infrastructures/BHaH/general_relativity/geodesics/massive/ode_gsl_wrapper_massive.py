"""
Generate the C wrapper function to interface with the GSL ODE solver.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def ode_gsl_wrapper_massive() -> None:
    """
    Generate the C wrapper function to interface with the GSL ODE solver.

    This function acts as a bridge, unpacking a generic `void*` pointer
    to access simulation data, calling our project-specific dispatchers and
    engines, and returning the computed derivatives to the GSL solver.
    """
    name = "ode_gsl_wrapper_massive"
    cfunc_type = "int"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "gsl/gsl_errno.h"]
    desc = "GSL wrapper for the massive particle geodesic ODEs."
    params = "double t, const double y[8], double f[8], void *params"
    body = r"""
    (void)t; // Mark proper time 't' as unused to avoid compiler warnings.
    gsl_params *gsl_parameters = (gsl_params *)params;
    metric_struct g4DD;
    connection_struct conn;
    g4DD_metric(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &g4DD);
    connections(gsl_parameters->commondata, gsl_parameters->params, gsl_parameters->metric, y, &conn);
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
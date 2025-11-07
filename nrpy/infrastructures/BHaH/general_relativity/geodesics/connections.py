"""
Generate the C dispatcher function for Christoffel symbol calculations.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def connections() -> None:
    """
    Generate and register the C dispatcher for Christoffel symbol calculations.

    Similar to g4DD_metric(), this function generates a C dispatcher that calls
    the appropriate worker function based on the runtime `metric->type`.
    """
    name = "connections"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "stdio.h",
        "stdlib.h",
    ]
    desc = "Dispatcher to compute Christoffel symbols for the chosen analytic metric."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  connection_struct *restrict conn"""
    body = r"""
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            con_kerr_schild(commondata, params, y_pos, conn);
            break;
        case Schwarzschild_Standard:
            con_schwarzschild_cartesian(commondata, params, y_pos, conn);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in connections().n", metric->type);
            exit(1);
    }
    """

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
    )
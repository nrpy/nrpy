"""
Generate the C dispatcher function for metric calculations.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def g4dd_metric() -> None:
    """
    Generate and register the C dispatcher for metric calculations.

    This function generates a C function `g4DD_metric()` that contains a switch
    statement. Based on the runtime `metric->type`, it calls the appropriate
    specialized worker function (e.g., `g4DD_kerr_schild`).
    """
    name = "g4DD_metric"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Dispatcher to compute the 4-metric g_munu for the chosen analytic metric."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric,
                  const double y[8],
                  metric_struct *restrict metric_out"""
    body = r"""
    const double y_pos[4] = {y[0], y[1], y[2], y[3]};
    switch(metric->type) {
        case Schwarzschild:
        case Kerr:
            g4DD_kerr_schild(commondata, params, y_pos, metric_out);
            break;
        case Schwarzschild_Standard:
            g4DD_schwarzschild_cartesian(commondata, params, y_pos, metric_out);
            break;
        default:
            fprintf(stderr, "Error: MetricType %d not supported in g4DD_metric().n", metric->type);
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
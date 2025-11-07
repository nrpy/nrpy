"""
Generates the C placeholder for the external numerical interpolation engine.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def placeholder_interpolation_engine() -> None:
    """
    Generate and register the high-fidelity C placeholder interpolation engine.

    This function generates the C code for `placeholder_interpolation_engine()`.
    This function serves as a stand-in for a future, high-performance numerical
    interpolation library. It mimics the required batch-processing API but
    computes the metric and Christoffel symbols by calling the existing, trusted
    analytic C worker functions. This allows for the complete validation of the
    numerical pipeline's control flow using a known-good analytic baseline.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""@brief Placeholder for the external batch-processing interpolation engine.

    ========================================================================
    ================== THIS IS A VALIDATION PLACEHOLDER ==================
    This function will be replaced by the high-performance numerical engine.
    ========================================================================

    It mimics the required API but computes the metric and Christoffel symbols
    by calling the low-level ANALYTIC WORKER functions directly. This provides
    a ground-truth analytic result for validating the numerical control flow.

    @param[in]  num_photons         The number of photons in the batch.
    @param[in]  requests            An array of photon_request_t structs.
    @param[out] metric_outputs      An array to be filled with metric components.
    @param[out] conn_outputs        An array to be filled with Christoffel symbols.
    @param[in]  commondata          Pointer to commondata struct with runtime parameters.
    @param[in]  params              Pointer to params struct (unused, for signature compatibility).
    @param[in]  metric              Pointer to the metric_params struct specifying the metric type.
    """
    name = "placeholder_interpolation_engine"
    params = """int num_photons,
                const photon_request_t requests[],
                metric_struct metric_outputs[],
                connection_struct conn_outputs[],
                const commondata_struct *restrict commondata,
                const params_struct *restrict params,
                const metric_params *restrict metric"""

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // This function loops through each request and computes the metric and
    // Christoffels individually by calling the high-level dispatchers.
    #pragma omp parallel for
    for (int i = 0; i < num_photons; ++i) {
        // The analytic dispatcher functions expect a 9-element state vector.
        // To call them safely, we create a temporary padded array on the stack.
        // We only need to fill the first 4 position components from the request.
        // The other 5 components (momentum, path length) are not used by these
        // specific dispatcher functions.
        double y_padded[9];
        y_padded[0] = requests[i].pos[0]; // t
        y_padded[1] = requests[i].pos[1]; // x
        y_padded[2] = requests[i].pos[2]; // y
        y_padded[3] = requests[i].pos[3]; // z

        // Now, call the high-level DISPATCHERS with the correctly-sized array.
        // These functions will internally use the 'metric->type' to call the
        // correct low-level worker (e.g., g4DD_kerr_schild).
        g4DD_metric(commondata, params, metric, y_padded, &metric_outputs[i]);
        connections(commondata, params, metric, y_padded, &conn_outputs[i]);
    }
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )
# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/interpolation_kernel.py
r"""
Generates the CUDA/OpenMP function for tensor interpolation.

This module registers a C function that runs the parallel evaluation of the
spacetime metric and Christoffel symbols for a photon chunk. It generates a CUDA
kernel or an OpenMP loop for the selected parallelization target. Photon states,
metric components, and connection components use flattened Structure-of-Arrays
layouts. Each parallel iteration copies one photon state into a local array, invokes
the metric and connection evaluators, and writes their tensor components to output
arrays.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def interpolation_kernel(spacetime_name: str) -> None:
    r"""
    Register the CUDA/OpenMP function for tensor interpolation.

    The generated parallel calculation reads the photon state vector $f^{\mu}$,
    evaluates the metric and connection components using the specified spacetime
    evaluators, and writes the resulting tensors to output arrays.

    :param spacetime_name: The string identifier for the target numerical spacetime.
    :raises ValueError: If the provided spacetime_name string is empty.
    """
    if not spacetime_name:
        raise ValueError("spacetime_name must contain a valid string identifier.")

    parallelization = par.parval_from_str("parallelization")

    metric_worker = f"g4DD_metric_{spacetime_name}"
    conn_worker = f"connections_{spacetime_name}"

    metric_c_code = cfc.CFunction_dict[metric_worker].full_function
    conn_c_code = cfc.CFunction_dict[conn_worker].full_function

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "double *restrict",
        "d_connection_bundle": "double *restrict",
        "chunk_size": "const long int",
    }

    arg_dict_host = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "double *restrict",
        "d_connection_bundle": "double *restrict",
        "chunk_size": "const long int",
    }

    # Pass commondata explicitly to the OpenMP function.
    if parallelization != "cuda":
        arg_dict_cuda["commondata"] = "const commondata_struct *restrict"
        arg_dict_host["commondata"] = "const commondata_struct *restrict"

    if parallelization == "cuda":
        loop_preamble = """
    //==========================================
    // CUDA THREAD IDENTIFICATION
    //==========================================
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Ignore CUDA threads beyond the active ray chunk.
    if (i >= chunk_size) return;
    """
        cd_ptr = "&d_commondata"
        loop_postamble = ""
    else:
        loop_preamble = """
    //==========================================
    // OPENMP PARALLEL LOOP
    //==========================================
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        cd_ptr = "commondata"
        loop_postamble = "    } // END LOOP: for i over chunk_size rays"

    core_math = rf"""
    //==========================================
    // MACRO DEFINITIONS FOR ARRAY ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state array using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric array.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_CONN maps a component to the flattened Christoffel connection array.
    #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    //==========================================
    // STATE UNPACKING
    //==========================================
    double f_local[9]; // Local array storing the 9-component state vector $f^{{\mu}}$.
    int comp; // Loop index for iterating over the tensor components.
    for (comp = 0; comp < 9; ++comp) {{
        // Load one photon's state-vector components into a local array.
        f_local[comp] = d_f_bundle[IDX_F(comp, i)]; // Component of the photon state vector $f^{{\mu}}$.
    }} // END LOOP: for comp over 9 state vector components

    //==========================================
    // METRIC TENSOR EVALUATION
    //==========================================
    double metric_local[10]; // Local array storing the 10 upper-triangular components of $g_{{\mu\nu}}$.

    // Evaluate the spacetime metric geometry.
    {metric_worker}({cd_ptr}, f_local, metric_local);
    //==========================================
    // METRIC OUTPUT
    //==========================================
    for (comp = 0; comp < 10; ++comp) {{
        // Write the computed metric components $g_{{\mu\nu}}$ to the output array.
        d_metric_bundle[IDX_METRIC(comp, i)] = metric_local[comp]; // Component of the spacetime metric $g_{{\mu\nu}}$.
    }} // END LOOP: for comp over 10 metric components

    //==========================================
    // CHRISTOFFEL CONNECTION EVALUATION
    //==========================================
    // Conditional logic skips connection calculation during the initialization phase if the pointer is NULL.
    if (d_connection_bundle != NULL) {{
        // Local array storing the 40 components of $\Gamma^{{\alpha}}_{{\beta\gamma}}$.
        double Gamma_local[40];

        // Evaluate the Christoffel symbols.
        {conn_worker}({cd_ptr}, f_local, Gamma_local);

        //==========================================
        // CONNECTION OUTPUT
        //==========================================
        for (comp = 0; comp < 40; ++comp) {{
            // Write the computed connection components $\Gamma^{{\alpha}}_{{\beta\gamma}}$ to the output array.
            d_connection_bundle[IDX_CONN(comp, i)] = Gamma_local[comp];
        }} // END LOOP: for comp over 40 connection components
    }} // END IF: d_connection_bundle is not NULL

    //==========================================
    // MACRO CLEANUP
    //==========================================
    // Undefine local indexing macros before the next generated function.
    #undef IDX_F
    #undef IDX_METRIC
    #undef IDX_CONN
"""

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    prefunc_kernel, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name=f"interpolation_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    prefunc = "\n\n".join([metric_c_code, conn_c_code, prefunc_kernel])

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = rf""" Evaluates the {spacetime_name} metric and connection for a photon chunk.

    @param commondata Pointer to spacetime parameters shared by all rays.
    @param d_f_bundle Pointer to the state vector array $f^{{\mu}}$ in memory.
    @param d_metric_bundle Pointer to the destination metric array $g_{{\mu\nu}}$ in memory.
    @param d_connection_bundle Pointer to the destination connection array $\Gamma^{{\alpha}}_{{\beta\gamma}}$ in memory.
    @param chunk_size The number of active rays in the current ray chunk.
    @param stream_idx Work-array index; CUDA uses the corresponding stream.
    """

    cfunc_type = "void"

    name = f"interpolation_kernel_{spacetime_name}"

    params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict d_f_bundle, "
        "double *restrict d_metric_bundle, "
        "double *restrict d_connection_bundle, "
        "const long int chunk_size,"
        "const int stream_idx"
    )

    body = launch_code

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
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

r"""
Provides the global memory kernel and orchestrator for the interpolation engine.

This module evaluates the spacetime metric $g_{\mu\nu}$ and Christoffel symbols
$\Gamma^{\alpha}_{\beta\gamma}$ for a batch of photons. It operates on Structure
of Arrays (SoA) bundles, ensuring memory coalescence and minimizing latency during
the integration step.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def interpolation_kernel(spacetime_name: str) -> None:
    r"""
    Register the global kernel for tensor interpolation.

    This kernel unpacks the photon state vector $f^{\mu}$ from global memory,
    evaluates the metric and connection components using the specified spacetime
    evaluators, and writes the resulting tensors back to memory bundles.

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

    # Pass commondata explicitly when not using CUDA's global memory
    if parallelization != "cuda":
        arg_dict_cuda["commondata"] = "const commondata_struct *restrict"
        arg_dict_host["commondata"] = "const commondata_struct *restrict"

    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $i$ represents the global thread index mapped to a specific photon ray.
    const long int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return;
    """
        cd_ptr = "&d_commondata"
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int i = 0; i < chunk_size; i++) {
    """
        cd_ptr = "commondata"
        loop_postamble = "    } // End OpenMP loop"

    core_math = rf"""
    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_CONN maps a component to the flattened Christoffel connection bundle.
    #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

// --- STATE UNPACKING ---
    double f_local[9]; // Local register array storing the 9-component state vector $f^{{\mu}}$.
    int comp; // Loop index for iterating over the tensor components.
    for (comp = 0; comp < 9; ++comp) {{
        // Load the state vector components from global memory into registers.
        f_local[comp] = d_f_bundle[IDX_F(comp, i)]; // Component of the photon state vector $f^{{\mu}}$.
    }}

    // --- METRIC TENSOR EVALUATION ---
    double metric_local[10]; // Local register array storing the 10 upper-triangular components of $g_{{\mu\nu}}$.

    // Evaluate the spacetime metric geometry.
    {metric_worker}({cd_ptr}, f_local, metric_local);
    // --- GLOBAL MEMORY WRITE (METRIC) ---
    for (comp = 0; comp < 10; ++comp) {{
        // Write the computed metric components $g_{{\mu\nu}}$ back to the global memory bundle.
        d_metric_bundle[IDX_METRIC(comp, i)] = metric_local[comp]; // Component of the spacetime metric $g_{{\mu\nu}}$.
    }}

    // --- CHRISTOFFEL CONNECTION EVALUATION ---
    // Conditional logic skips connection calculation during the initialization phase if the pointer is NULL.
    if (d_connection_bundle != NULL) {{
        // Local register array storing the 40 components of $\Gamma^{{\alpha}}_{{\beta\gamma}}$.
        double Gamma_local[40];

        // Evaluate the Christoffel symbols.
        {conn_worker}({cd_ptr}, f_local, Gamma_local);

        // --- GLOBAL MEMORY WRITE (CONNECTION) ---
        for (comp = 0; comp < 40; ++comp) {{
            // Write the computed connection components $\Gamma^{{\alpha}}_{{\beta\gamma}}$ to the global memory bundle.
            d_connection_bundle[IDX_CONN(comp, i)] = Gamma_local[comp];
        }}
    }}

    // --- MACRO CLEANUP ---
    // Undefine macros to ensure hermetic compilation and prevent redefinition errors.
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

    prefunc_kernel, launch_code = generate_kernel_and_launch_code(
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

    desc = rf"""@brief Orchestrates the memory kernel for the {spacetime_name} interpolation engine.

    @param d_f_bundle Pointer to the state vector bundle $f^{{\mu}}$ in memory.
    @param d_metric_bundle Pointer to the destination metric bundle $g_{{\mu\nu}}$ in memory.
    @param d_connection_bundle Pointer to the destination connection bundle $\Gamma^{{\alpha}}_{{\beta\gamma}}$ in memory.
    @param chunk_size The number of active rays in the current bundle batch.
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

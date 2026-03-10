r"""
Provides the CUDA kernel and host-side orchestrator for the interpolation engine.

This module provides the global memory kernel responsible for evaluating the spacetime
metric $g_{\mu\nu}$ and Christoffel symbols $\Gamma^{\alpha}_{\beta\gamma}$ for a batch
of photons. It operates on Structure of Arrays (SoA) bundles residing in VRAM, 
ensuring memory coalescence and minimizing latency during the integration step.

Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def interpolation_kernel(spacetime_name: str) -> None:
    r"""
    Register the global CUDA kernel for tensor interpolation.

    This kernel unpacks the photon state vector $f^{\mu}$ from global VRAM, 
    evaluates the metric and connection components using the specified spacetime 
    evaluators, and writes the resulting tensors back to VRAM bundles.

    :param spacetime_name: The string identifier for the target numerical spacetime.
    :raises ValueError: If the provided spacetime_name string is empty.
    """
    # Input validation.
    if not spacetime_name:
        raise ValueError("spacetime_name must contain a valid string identifier.")

    # Identify the external device functions based on the selected spacetime.
    metric_worker = f"g4DD_metric_{spacetime_name}"
    conn_worker = f"connections_{spacetime_name}"

    # Extract function bodies for inlining from the CFunction dictionary.
    # This ensures they are visible to the compiler within the same .cu file.
    metric_c_code = cfc.CFunction_dict[metric_worker].full_function
    conn_c_code = cfc.CFunction_dict[conn_worker].full_function

    # Define the argument dictionary for the CUDA kernel generation.
    # Note: commondata is removed here because it is accessed directly via constant memory.
    arg_dict = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "double *restrict",
        "d_connection_bundle": "double *restrict",
        "chunk_size": "const int"
    }

    # Define the GPU kernel body.
    kernel_body = fr"""
    // --- THREAD IDENTIFICATION ---
    // The identifier i represents the global thread index mapped to a specific photon ray.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Hardware Justification: Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_CONN maps a component to the flattened Christoffel connection bundle.
    #define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    // --- STATE UNPACKING ---
    // Local register array storing the 9-component state vector $f^{{\mu}}$.
    // Hardware Justification: Required to interface with external evaluators expecting contiguous memory.
    double f_local[9];
    for (int comp = 0; comp < 9; ++comp) {{
        f_local[comp] = ReadCUDA(&d_f_bundle[IDX_F(comp, i)]);
    }}

    // --- METRIC TENSOR EVALUATION ---
    // Local register array storing the 10 upper-triangular components of $g_{{\mu\nu}}$.
    double metric_local[10];
    
    // Evaluate the spacetime metric geometry.
    // Hardware Justification: Uses the globally available __constant__ d_commondata cache.
    {metric_worker}(&d_commondata, f_local, metric_local);

    // --- GLOBAL VRAM WRITE (METRIC) ---
    // Write the computed metric components to global VRAM.
    for (int comp = 0; comp < 10; ++comp) {{
        WriteCUDA(&d_metric_bundle[IDX_METRIC(comp, i)], metric_local[comp]);
    }}

    // --- CHRISTOFFEL CONNECTION EVALUATION ---
    // Conditional logic skips connection calculation during the initialization phase if the pointer is NULL.
    if (d_connection_bundle != NULL) {{
        // Local register array storing the 40 components of $\Gamma^{{\alpha}}_{{\beta\gamma}}$.
        double Gamma_local[40];
        
        // Evaluate the Christoffel symbols.
        {conn_worker}(&d_commondata, f_local, Gamma_local);

        // --- GLOBAL VRAM WRITE (CONNECTION) ---
        // Write the computed connection components to global VRAM.
        for (int comp = 0; comp < 40; ++comp) {{
            WriteCUDA(&d_connection_bundle[IDX_CONN(comp, i)], Gamma_local[comp]);
        }}
    }}

    // --- MACRO CLEANUP ---
    // Undefine macros to ensure hermetic compilation and prevent redefinition errors.
    #undef IDX_F
    #undef IDX_METRIC
    #undef IDX_CONN
    """

    # Generate the kernel and the C host wrapper.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx"
    }

    prefunc_kernel, body = generate_kernel_and_launch_code(
        kernel_name=f"interpolation_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
        thread_tiling_macro_suffix="RKF45"
    )

    # Consolidation: Prepend the worker functions to satisfy the inlining mandate.
    prefunc = "\n\n".join([metric_c_code, conn_c_code, prefunc_kernel])

    # Define arguments for C-Function registration strictly before the call.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "cuda_intrinsics.h"]
    
    desc = fr"""@brief Orchestrates the CUDA kernel for the {spacetime_name} interpolation engine.
    
    @param d_f_bundle Pointer to the state vector bundle $f^{{\mu}}$ in VRAM.
    @param d_metric_bundle Pointer to the destination metric bundle $g_{{\mu\nu}}$ in VRAM.
    @param d_connection_bundle Pointer to the destination connection bundle $\Gamma^{{\alpha}}_{{\beta\gamma}}$ in VRAM.
    @param chunk_size The number of active rays in the current bundle batch.
    """

    cfunc_type = "void"
    name = f"interpolation_kernel_{spacetime_name}"
    
    params = (
        "const double *restrict d_f_bundle, "
        "double *restrict d_metric_bundle, "
        "double *restrict d_connection_bundle, "
        "const int chunk_size,"
        "const int stream_idx"
    )

    # Register the complete C function using the canonical Master Order.
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body
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
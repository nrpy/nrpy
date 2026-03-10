# nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/normalization_constraint.py
"""
Evaluates the normalization constraint along trajectories using a streaming bundle architecture.

This module implements a streaming bundle architecture to compute the scalar invariant $C = g_{\\mu\\nu} v^\\mu v^\\nu$.
It strictly manages VRAM usage by processing particles in limited batches, using pinned memory transfers to bridge the Host-Device gap.
Author: Dalton J. Moone.
"""

from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
)
from nrpy.helpers.loop import loop
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

def normalization_constraint(norm_expr: sp.Expr, PARTICLE: str) -> None:
    """
    Evaluates the C function and native kernel for normalization constraint evaluation.

    :param norm_expr: The SymPy expression for the contraction $g_{\\mu\\nu} v^\\mu v^\\nu$.
    :param PARTICLE: The type of particle (e.g., photon or massive).
    :raises ValueError: If the particle type is not supported.
    """
    if PARTICLE == "massive":
        vec_desc = "4-velocity $u^\\mu$"
        expected_val = "-1.0"
        array_size = 8
    elif PARTICLE == "photon":
        vec_desc = "4-momentum $p^\\mu$"
        expected_val = "0.0"
        array_size = 9
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    norm_struct_def = """
    // --- NORMALIZATION CONSTRAINT STRUCTURE ---
    // Defines the physical normalization constraint evaluated along a trajectory.
    // Hardware Justification: This Structure of Arrays (AoS) definition is injected into the global BHaH_defines.h header to ensure uniform memory mapping across the Host orchestrator and VRAM kernels.
    typedef struct {
    double C;   // Scalar invariant $C = g_{\\mu\\nu} v^\\mu v^\\nu$.
    } normalization_constraint_t;
    """

    # Register the struct definition to the global header generation pipeline.
    Bdefines_h.register_BHaH_defines("normalization_constraint", norm_struct_def)

    # Generate the highly optimized math evaluation block.
    math_kernel = ccg.c_codegen(
        [norm_expr],
        ["d_norm_bundle[c].C"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    ).replace("SIMD", "CUDA")

    # Define the C-kernel name string.
    kernel_name = f"calculate_normalization_constraint_{PARTICLE}_kernel"

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_norm_bundle": "normalization_constraint_t *restrict",
        "current_chunk_size": "const long int",
    }

    # Dynamically generate the unpacking logic based on the specific vector coordinates.
    preamble_lines = [
        "    // --- COMPONENT HYDRATION ---",
        "    // Hardware Justification: Unpack 4-vector and metric components using strict SoA macros directly from the VRAM bundle."
    ]

    for i in range(4):
        preamble_lines.append(f"    const double vU{i} = d_f_bundle[IDX_LOCAL({i+4}, c, BUNDLE_CAPACITY)]; // 4-vector component $v^{i}$ evaluated from VRAM.")
        preamble_lines.append(f"    (void)vU{i}; // Suppress unused variable warning for $v^{i}$.")

    # Map the metric tensor components.
    k = 0
    for i in range(4):
        for j in range(i, 4):
            preamble_lines.append(f"    const double metric_g4DD{i}{j} = d_metric_bundle[IDX_LOCAL({k}, c, BUNDLE_CAPACITY)]; // Metric tensor component $g_{{{i}{j}}}$ evaluated from VRAM.")
            preamble_lines.append(f"    (void)metric_g4DD{i}{j}; // Suppress unused variable warning for $g_{{{i}{j}}}$.")
            k += 1

    preamble = "\n".join(preamble_lines)

    # Define the core __global__ kernel body executed on the device.
    kernel_body = f"""
    // --- MACRO DEFINITIONS ---
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    // Layout: [Component][RayID]
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(comp, ray_id, N) ((comp) * (N) + (ray_id))
    #endif

    // --- THREAD IDENTIFICATION ---
    // Local 1D thread mapping within the current VRAM bundle.
    // Thread ID maps to a unique trajectory index $c$.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Global thread evaluation index $c$.

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid bundle addresses.
    if (c >= current_chunk_size) return;

    {preamble}

    // --- CONSTRAINT EVALUATION ---
    // Evaluates the analytic SymPy expression for the normalization constraint.
    {math_kernel}
    """

    launch_dict = {
        "threads_per_block": [
            "BHAH_THREADS_IN_X_DIR_DEFAULT",
            "BHAH_THREADS_IN_Y_DIR_DEFAULT",
            "BHAH_THREADS_IN_Z_DIR_DEFAULT",
        ],
        "blocks_per_grid": [
            "(current_chunk_size + BHAH_THREADS_IN_X_DIR_DEFAULT - 1) / BHAH_THREADS_IN_X_DIR_DEFAULT",
            "1",
            "1",
        ],
    }

    # Generate the kernel definition and the internal launch string.
    prefunc, launch_body = generate_kernel_and_launch_code(
        kernel_name=kernel_name,
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        thread_tiling_macro_suffix="DEFAULT",
        cfunc_decorators="__global__",
    )

    # Host-side chunked execution loop.
    loop_body = f"""
    // --- BUNDLE SIZING ---
    // Variable $current\_chunk\_size$ defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY); // Safely bounded active chunk parameter $current\_chunk\_size$.

    // --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---
    // Transfer the state vector $f^\\mu$ for the current bundle to VRAM.
    // Hardware Justification: Component-wise transfers bounded by BUNDLE_CAPACITY prevent PCIe saturation.
    for(int m=0; m<{array_size}; m++) {{ // Loop over all {array_size} elements of the state vector $f^\\mu$.
        cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
                   all_photons->f + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice); // Transmits component vector array Host-to-Device.
    }}

    // Transfer the metric tensor $g_{{\\mu\\nu}}$ for the current bundle to VRAM.
    // Hardware Justification: Component-wise transfers bounded by BUNDLE_CAPACITY prevent PCIe saturation.
    for(int m=0; m<10; m++) {{ // Loop over all 10 independent components of the symmetric metric tensor $g_{{\\mu\\nu}}$.
        cudaMemcpy(d_metric_bundle + (m * BUNDLE_CAPACITY),
                   all_metrics + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice); // Transmits metric tensor array Host-to-Device.
    }}

    // --- KERNEL LAUNCH ---
    {launch_body}

    // --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---
    // Retrieve the calculated normalization constraints back to the Pinned memory array.
    cudaMemcpy(norm_result + start_idx, d_norm_bundle,
               sizeof(normalization_constraint_t) * current_chunk_size, cudaMemcpyDeviceToHost); // Transmits physical diagnostic records Device-to-Host.
    """

    host_loop = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        loop_body=loop_body,
    )

    # 7. Variable Definition (The Master Order)
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    
    desc = f"""@brief Computes the normalization constraint for a batch of trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors $f^\\mu$.
    @param all_metrics The master Structure of Arrays containing the symmetric metric tensor $g_{{\\mu\\nu}}$.
    @param num_rays The total number of particle trajectories.
    @param norm_result The array of diagnostic constraint structures to be populated.

    Detailed algorithm:
    1. Allocates VRAM staging buffers for state vectors $f^\\mu$, metrics $g_{{\\mu\\nu}}$, and quantity results.
    2. Iterates over the global dataset in chunks of BUNDLE_CAPACITY.
    3. Transfers state and metric data Host-to-Device, computes expressions, and transfers Device-to-Host.
    
    Expected Value: {expected_val} for {vec_desc}."""

    cfunc_type = "void"
    name = f"calculate_normalization_constraint_{PARTICLE}"
    params = (
        "const PhotonStateSoA *restrict all_photons, "
        "const double *restrict all_metrics, "
        "const long int num_rays, "
        "normalization_constraint_t *restrict norm_result"
    )
    include_CodeParameters_h = False

    body = f"""
    // --- VRAM STAGING ALLOCATION ---
    // Hardware Justification: Allocate buffers statically sized to BUNDLE_CAPACITY to fit within 10GB VRAM limits.
    double *d_f_bundle; // Device pointer for the bundled state vector $f^\\mu$.
    double *d_metric_bundle; // Device pointer for the bundled metric tensor $g_{{\\mu\\nu}}$.
    normalization_constraint_t *d_norm_bundle; // Device pointer for the bundled diagnostic outputs.

    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * {array_size} * BUNDLE_CAPACITY); // Allocates device target array for state vector.
    BHAH_MALLOC_DEVICE(d_metric_bundle, sizeof(double) * 10 * BUNDLE_CAPACITY); // Allocates device target array for metric tensor.
    BHAH_MALLOC_DEVICE(d_norm_bundle, sizeof(normalization_constraint_t) * BUNDLE_CAPACITY); // Allocates device return array.

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop}

    // --- VRAM CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle); // Releases memory for state vector $f^\\mu$.
    BHAH_FREE_DEVICE(d_metric_bundle); // Releases memory for metric tensor $g_{{\\mu\\nu}}$.
    BHAH_FREE_DEVICE(d_norm_bundle); // Releases memory for diagnostic outputs.
    """

    # 8. Function Registration
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
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
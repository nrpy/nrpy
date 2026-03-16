# nrpy/equations/general_relativity/geodesics/geodesic_diagnostics/normalization_constraint.py
r"""
Evaluate the normalization constraint along trajectories using a streaming bundle architecture.

This module implements a pure mathematical execution kernel to compute the scalar invariant $C = g_{\mu\nu} v^\mu v^\nu$.

Author: Dalton J. Moone.
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
)


def normalization_constraint(norm_expr: sp.Expr, PARTICLE: str) -> None:
    r"""
    Evaluate the C function and native kernel for normalization constraint evaluation.

    :param norm_expr: The SymPy expression for the contraction $g_{\mu\nu} v^\mu v^\nu$.
    :param PARTICLE: The type of particle (e.g., photon or massive).
    :raises ValueError: If the particle type is not supported.
    """
    if PARTICLE == "massive":
        vec_desc = r"4-velocity $u^\mu$"
        expected_val = "-1.0"
    elif PARTICLE == "photon":
        vec_desc = r"4-momentum $p^\mu$"
        expected_val = "0.0"
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    # --- ARCHITECTURE DETECTION ---
    parallelization = par.parval_from_str("parallelization")

    norm_struct_def = r"""
    // --- NORMALIZATION CONSTRAINT STRUCTURE ---
    // Defines the physical normalization constraint evaluated along a trajectory.
    // Hardware Justification: This Structure of Arrays (AoS) definition is injected into the global BHaH_defines.h header to ensure uniform memory mapping across the Host orchestrator and computational threads.
    typedef struct {
        double C;   // Scalar invariant $C = g_{\\mu\\nu} v^\\mu v^\\nu$.
    } normalization_constraint_t;
    """

    # Register the struct definition to the global header generation pipeline.
    Bdefines_h.register_BHaH_defines("normalization_constraint", norm_struct_def)

    # Define the highly optimized math evaluation block.
    math_kernel = ccg.c_codegen(
        [norm_expr],
        ["d_norm_bundle[c].C"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )
    if parallelization == "cuda":
        math_kernel = math_kernel.replace("SIMD", "CUDA")

    # Define the C-kernel name string.
    kernel_name = f"normalization_constraint_{PARTICLE}_kernel"

    arg_dict = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_norm_bundle": "normalization_constraint_t *restrict",
        "current_chunk_size": "const long int",
    }

    # Dynamically generate the unpacking logic based on the specific vector coordinates.
    preamble_lines = [
        "    // --- COMPONENT HYDRATION ---",
        "    // Hardware Justification: Unpack 4-vector and metric components using strict SoA macros directly from the memory bundle.",
    ]

    for i in range(4):
        preamble_lines.append(
            f"    const double vU{i} = d_f_bundle[IDX_LOCAL({i+4}, c, BUNDLE_CAPACITY)]; // 4-vector component $v^{i}$ evaluated from memory."
        )
        preamble_lines.append(
            f"    (void)vU{i}; // Suppress unused variable warning for $v^{i}$."
        )

    # Map the metric tensor components.
    k = 0
    for i in range(4):
        for j in range(i, 4):
            preamble_lines.append(
                f"    const double metric_g4DD{i}{j} = d_metric_bundle[IDX_LOCAL({k}, c, BUNDLE_CAPACITY)]; // Metric tensor component $g_{{{i}{j}}}$ evaluated from memory."
            )
            preamble_lines.append(
                f"    (void)metric_g4DD{i}{j}; // Suppress unused variable warning for $g_{{{i}{j}}}$."
            )
            k += 1

    preamble = "\n".join(preamble_lines)

    # --- THE KERNEL SANDWICH ---
    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $c$ maps directly to a unique particle index in the global VRAM batch.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Global thread evaluation index $c$.

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (c >= current_chunk_size) return;
    """
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute particle trajectories across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int c = 0; c < current_chunk_size; c++) {
    """
        loop_postamble = "    } // End OpenMP loop"

    core_math = f"""
    // --- MACRO DEFINITIONS ---
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    // Layout: [Component][RayID]
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(comp, ray_id, N) ((comp) * (N) + (ray_id))
    #endif

    {preamble}

    // --- CONSTRAINT EVALUATION ---
    // Evaluates the analytic SymPy expression for the normalization constraint.
    {math_kernel}
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

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
        "stream": "stream_idx",
    }

    # Generate the kernel definition and the internal launch string.
    prefunc, launch_body = generate_kernel_and_launch_code(
        kernel_name=kernel_name,
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        thread_tiling_macro_suffix="DEFAULT",
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    # Variables strictly ordered immediately prior to registration
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]

    desc = rf"""@brief Computes the normalization constraint for a batch of trajectories.
    @param d_f_bundle The device pointer array containing the state vectors $f^\mu$.
    @param d_metric_bundle The device pointer array containing the symmetric metric tensor $g_{{\mu\nu}}$.
    @param d_norm_bundle The device array of diagnostic constraint structures to be populated.
    @param current_chunk_size The dynamically sized operational boundary for the active chunk.

    Expected Value: {expected_val} for {vec_desc}."""

    cfunc_type = "void"

    name = f"normalization_constraint_{PARTICLE}"

    params = (
        "const double *restrict d_f_bundle, "
        "const double *restrict d_metric_bundle, "
        "normalization_constraint_t *restrict d_norm_bundle, "
        "const long int current_chunk_size, "
        "const int stream_idx"
    )

    include_CodeParameters_h = False

    body = launch_body

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

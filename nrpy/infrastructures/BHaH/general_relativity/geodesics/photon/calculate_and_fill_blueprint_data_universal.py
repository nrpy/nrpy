"""
Module implements a "Streaming Bundle" architecture to project escaped photon trajectories.

Module evaluates the spatial coordinates of escaped photons and projects them onto the
celestial sphere. It strictly manages memory usage by processing photons in fixed batches,
ensuring compliance with hardware memory limits across both CPU and GPU execution contexts.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.loop import loop
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code


def calculate_and_fill_blueprint_data_universal() -> None:
    """Compute the universal blueprint data for escaped photon trajectories."""
    kernel_name = "calculate_and_fill_blueprint_data_universal_kernel"
    parallelization = par.parval_from_str("parallelization")

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_status_bundle": "const termination_type_t *restrict",
        "d_result_bundle": "blueprint_data_t *restrict",
        "current_chunk_size": "const long int",
    }

    arg_dict_host = {
        "d_f_bundle": "const double *restrict",
        "d_status_bundle": "const termination_type_t *restrict",
        "d_result_bundle": "blueprint_data_t *restrict",
        "current_chunk_size": "const long int",
    }

    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // The identifier $c$ represents the global thread index mapped to a specific photon ray.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Thread index $c$ maps to photon ID.

    // Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (c >= current_chunk_size) return; // Fatal unrecoverable error: Out-of-bounds thread execution aborted.
    """
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int c = 0; c < current_chunk_size; c++) { // Thread index $c$ maps to photon ID.
    """
        loop_postamble = "    } // End OpenMP loop"

    core_math = r"""
    // --- MACRO DEFINITIONS ---
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    // Layout: [Component][RayID]
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(c, ray_id, N) ((c) * (N) + (ray_id))
    #endif

    // --- STATUS SYNCHRONIZATION ---
    // Synchronize final exit status directly into the persistent blueprint array.
    d_result_bundle[c].termination_type = d_status_bundle[c]; // Stores final termination state.

    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // === Termination Dispatch: Celestial Sphere (Escape) ===
    if (d_status_bundle[c] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Unpack final coordinates from the bundled state vector $f^mu$ in memory.
        const double x_final = d_f_bundle[IDX_LOCAL(1, c, BUNDLE_CAPACITY)]; // Photon $x$-coordinate.
        const double y_final = d_f_bundle[IDX_LOCAL(2, c, BUNDLE_CAPACITY)]; // Photon $y$-coordinate.
        const double z_final = d_f_bundle[IDX_LOCAL(3, c, BUNDLE_CAPACITY)]; // Photon $z$-coordinate.

        // Compute the final radial distance $r = \sqrt{x^2 + y^2 + z^2}$.
        const double r_final = sqrt(x_final*x_final + y_final*y_final + z_final*z_final); // Radial distance $r$.

        // --- CELESTIAL PROJECTION ---
        // Map the Cartesian escape coordinates to the celestial sphere.
        d_result_bundle[c].final_theta = acos(z_final / r_final); // Polar angle $\theta$ relative to the $z$-axis.
        d_result_bundle[c].final_phi = atan2(y_final, x_final);   // Azimuthal angle $\phi$ in the $x$-$y$ plane.
    }
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

    prefunc, launch_body = generate_kernel_and_launch_code(
        kernel_name=kernel_name,
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_host,
        parallelization=parallelization,
        launch_dict=launch_dict,
        thread_tiling_macro_suffix="DEFAULT",
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    # Determine correct memory transfer semantics based on the target architecture.
    if parallelization == "cuda":
        memcpy_status = "cudaMemcpy(d_status_bundle, all_photons->status + start_idx, sizeof(termination_type_t) * current_chunk_size, cudaMemcpyHostToDevice);"
        memcpy_f = "cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY), all_photons->f + (m * num_rays) + start_idx, sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice);"
        memcpy_result_in = "cudaMemcpy(d_result_bundle, result + start_idx, sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyHostToDevice);"
        memcpy_result_out = "cudaMemcpy(result + start_idx, d_result_bundle, sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyDeviceToHost);"
        transfer_comment_in = "// --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---"
        transfer_comment_out = "// --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---"
    else:
        memcpy_status = "memcpy(d_status_bundle, all_photons->status + start_idx, sizeof(termination_type_t) * current_chunk_size);"
        memcpy_f = "memcpy(d_f_bundle + (m * BUNDLE_CAPACITY), all_photons->f + (m * num_rays) + start_idx, sizeof(double) * current_chunk_size);"
        memcpy_result_in = "memcpy(d_result_bundle, result + start_idx, sizeof(blueprint_data_t) * current_chunk_size);"
        memcpy_result_out = "memcpy(result + start_idx, d_result_bundle, sizeof(blueprint_data_t) * current_chunk_size);"
        transfer_comment_in = "// --- MEMORY TRANSFER (HOST TO STAGING BUFFER) ---"
        transfer_comment_out = "// --- MEMORY TRANSFER (STAGING BUFFER TO HOST) ---"

    loop_body = f"""
    // Variable current_chunk_size defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY); // Active chunk size.

    {transfer_comment_in}
    // Transfer the status array for the current bundle to supply the kernel with termination states.
    {memcpy_status} // Transfer of status.

    // Transfer the 9-component state vector $f^mu$ for the current bundle for coordinate unpacking.
    for(int m=0; m<9; m++) {{ // Iterate over the 9 components of the $f^mu$ state vector.
        {memcpy_f} // Transfer of $f^mu$.
    }}

    // Pre-load existing results to prevent overwriting valid memory with garbage during the subsequent transfer.
    {memcpy_result_in} // Transfer of previous results.

    // --- KERNEL LAUNCH ---
    {launch_body}

    {transfer_comment_out}
    // Retrieve the calculated blueprint results for the current bundle to persist the final data.
    {memcpy_result_out} // Transfer of final results.
    """

    host_loop = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        loop_body=loop_body,
    )

    # 7. Variable Definition
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r"""@brief Evaluates the blueprint data for a batch of photon trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors.
    @param num_rays The total number of photon trajectories.
    @param result The array of blueprint data structures to be populated.
    @param stream_idx The stream index identifier for asynchronous scheduling.

    Detailed algorithm:
    1. Allocates staging buffers for state vectors, status, and results.
    2. Iterates over the global dataset in chunks of BUNDLE_CAPACITY.
    3. Transfers data, computes projections, and transfers results back.
    4. Evaluates final 3D positions onto a celestial sphere $(\theta, \phi)$ for escaped rays."""
    cfunc_type = "void"
    name = "calculate_and_fill_blueprint_data_universal"
    params = "const PhotonStateSoA *restrict all_photons, const long int num_rays, blueprint_data_t *restrict result, const int stream_idx"
    include_CodeParameters_h = False
    body = f"""
    // --- STAGING ALLOCATION ---
    // Pointers for the bundled data processing.
    double *d_f_bundle; // Buffer for state vector $f^mu$.
    termination_type_t *d_status_bundle; // Buffer for photon termination status.
    blueprint_data_t *d_result_bundle; // Buffer for calculated blueprint data.

    // Hardware Justification: Allocate buffers sized to BUNDLE_CAPACITY to fit within hardware memory limits.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^mu$ buffer.
    BHAH_MALLOC_DEVICE(d_status_bundle, sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate status buffer.
    BHAH_MALLOC_DEVICE(d_result_bundle, sizeof(blueprint_data_t) * BUNDLE_CAPACITY); // Allocate results buffer.

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop}

    // --- BUFFER CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle); // Free $f^mu$ buffer.
    BHAH_FREE_DEVICE(d_status_bundle); // Free status buffer.
    BHAH_FREE_DEVICE(d_result_bundle); // Free results buffer.
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

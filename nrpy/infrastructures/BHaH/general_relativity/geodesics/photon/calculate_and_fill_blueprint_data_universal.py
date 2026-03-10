"""
This module implements a "Streaming Bundle" architecture to project escaped photon trajectories
onto the celestial sphere. It strictly manages VRAM usage by processing photons in 32k batches,
using pinned memory transfers to bridge the Host-Device gap, ensuring compliance with the
RTX 3080 10GB memory limit.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code
from nrpy.helpers.loop import loop

def calculate_and_fill_blueprint_data_universal() -> None:
    """
    This function computes the universal blueprint data for escaped photon trajectories.

    :param None: This function requires no arguments.
    :raises Exception: Propagates nrpy core exceptions on generation failure.
    """
    kernel_name = "calculate_and_fill_blueprint_data_universal_kernel"

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_status_bundle": "const termination_type_t *restrict",
        "d_result_bundle": "blueprint_data_t *restrict",
        "current_chunk_size": "const long int",
    }

    kernel_body = r"""
    // --- MACRO DEFINITIONS ---
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    // Layout: [Component][RayID]
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(c, ray_id, N) ((c) * (N) + (ray_id))
    #endif

    // --- THREAD IDENTIFICATION ---
    // Local 1D thread mapping within the current VRAM bundle.
    // Thread ID maps to a unique photon index.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Thread index $c$ maps to photon ID.

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid bundle addresses.
    if (c >= current_chunk_size) return; // Fatal unrecoverable error: Out-of-bounds thread execution aborted.

    // --- STATUS SYNCHRONIZATION ---
    // Synchronize final exit status directly into the persistent blueprint array.
    d_result_bundle[c].termination_type = d_status_bundle[c]; // Stores final termination state.

    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // === Termination Dispatch: Celestial Sphere (Escape) ===
    if (d_status_bundle[c] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Unpack final coordinates from the bundled state vector $f^mu$ in VRAM.
        // Mapped to thread-local variables to minimize global memory reads.
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

    launch_dict = {
        "threads_per_block": [
            "BHAH_THREADS_IN_X_DIR_DEFAULT", 
            "BHAH_THREADS_IN_Y_DIR_DEFAULT", 
            "BHAH_THREADS_IN_Z_DIR_DEFAULT"
        ],
        "blocks_per_grid": [
            "(current_chunk_size + BHAH_THREADS_IN_X_DIR_DEFAULT - 1) / BHAH_THREADS_IN_X_DIR_DEFAULT", 
            "1", 
            "1"
        ],
        "stream": "stream_idx"
    }

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

    loop_body = f"""
    // Variable current_chunk_size defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY); // Active chunk size.

    // --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---
    // Transfer the status array for the current bundle from Host-to-Device to supply the kernel with termination states.
    cudaMemcpy(d_status_bundle, all_photons->status + start_idx,
               sizeof(termination_type_t) * current_chunk_size, cudaMemcpyHostToDevice); // Host-to-Device transfer of status.

    // Transfer the 9-component state vector $f^mu$ for the current bundle from Host-to-Device for coordinate unpacking.
    for(int m=0; m<9; m++) {{ // Iterate over the 9 components of the $f^mu$ state vector.
        cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
                   all_photons->f + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice); // Host-to-Device transfer of $f^mu$.
    }}

    // Pre-load existing results from Host-to-Device to prevent overwriting valid memory with garbage during the subsequent Device-to-Host transfer.
    cudaMemcpy(d_result_bundle, result + start_idx,
               sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyHostToDevice); // Host-to-Device transfer of previous results.

    // --- KERNEL LAUNCH ---
    {launch_body}

    // --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---
    // Retrieve the calculated blueprint results for the current bundle from Device-to-Host to persist the final data.
    cudaMemcpy(result + start_idx, d_result_bundle,
               sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyDeviceToHost); // Device-to-Host transfer of final results.
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
    prefunc = prefunc
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    desc = r"""@brief Evaluates the blueprint data for a batch of photon trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors.
    @param num_rays The total number of photon trajectories.
    @param result The array of blueprint data structures to be populated.

    Detailed algorithm:
    1. Allocates VRAM staging buffers for state vectors, status, and results.
    2. Iterates over the global dataset in chunks of BUNDLE_CAPACITY.
    3. Transfers data Host-to-Device, computes projections, and transfers Device-to-Host.
    4. Evaluates final 3D positions onto a celestial sphere $(\theta, \phi)$ for escaped rays."""
    cfunc_type = "void"
    name = "calculate_and_fill_blueprint_data_universal"
    params = "const PhotonStateSoA *restrict all_photons, const long int num_rays, blueprint_data_t *restrict result, const int stream_idx"
    include_CodeParameters_h = False
    body = f"""
    // --- VRAM STAGING ALLOCATION ---
    // Device pointers for the bundled data processing.
    double *d_f_bundle; // Device buffer for state vector $f^mu$.
    termination_type_t *d_status_bundle; // Device buffer for photon termination status.
    blueprint_data_t *d_result_bundle; // Device buffer for calculated blueprint data.

    // Hardware Justification: Allocate buffers sized to BUNDLE_CAPACITY to fit within VRAM.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocate $f^mu$ buffer.
    BHAH_MALLOC_DEVICE(d_status_bundle, sizeof(termination_type_t) * BUNDLE_CAPACITY); // Allocate status buffer.
    BHAH_MALLOC_DEVICE(d_result_bundle, sizeof(blueprint_data_t) * BUNDLE_CAPACITY); // Allocate results buffer.

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop}

    // --- VRAM CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle); // Free $f^mu$ buffer.
    BHAH_FREE_DEVICE(d_status_bundle); // Free status buffer.
    BHAH_FREE_DEVICE(d_result_bundle); // Free results buffer.
    """

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
"""
Generates the CUDA kernel and launch code to finalize photon blueprint data.

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
    Generate the C function and native kernel for universal blueprint data finalization.

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
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid bundle addresses.
    if (c >= current_chunk_size) return;

    // --- STATUS SYNCHRONIZATION ---
    // Synchronize final exit status directly into the persistent blueprint array.
    d_result_bundle[c].termination_type = d_status_bundle[c];

    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // === Termination Dispatch: Celestial Sphere (Escape) ===
    if (d_status_bundle[c] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Unpack final coordinates from the bundled state vector $f^\mu$ in VRAM.
        // Mapped to thread-local variables to minimize global memory reads.
        const double x_final = d_f_bundle[IDX_LOCAL(1, c, BUNDLE_CAPACITY)];
        const double y_final = d_f_bundle[IDX_LOCAL(2, c, BUNDLE_CAPACITY)];
        const double z_final = d_f_bundle[IDX_LOCAL(3, c, BUNDLE_CAPACITY)];

        // Compute the final radial distance $r = \sqrt{x^2 + y^2 + z^2}$.
        const double r_final = sqrt(x_final*x_final + y_final*y_final + z_final*z_final);

        // --- CELESTIAL PROJECTION ---
        // Map the Cartesian escape coordinates to the celestial sphere to prevent division by zero.
        if (r_final > 1e-9) {
            d_result_bundle[c].final_theta = acos(z_final / r_final); // Polar angle $\theta$ relative to the $z$-axis.
            d_result_bundle[c].final_phi = atan2(y_final, x_final);   // Azimuthal angle $\phi$ in the $xy$-plane.
        }
    }
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(current_chunk_size + 256 - 1) / 256", "1", "1"],
    }

    prefunc, launch_body = generate_kernel_and_launch_code(
        kernel_name=kernel_name,
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]

    desc = r"""@brief Finalizes the blueprint data for a batch of photon trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors.
    @param num_rays The total number of photon trajectories.
    @param result The array of blueprint data structures to be populated.

    Detailed algorithm:
    1. Allocates VRAM staging buffers for state vectors, status, and results.
    2. Iterates over the global dataset in chunks of BUNDLE_CAPACITY.
    3. Transfers data Host->Device, computes projections, and transfers Device->Host.
    4. Projects final 3D positions onto a celestial sphere $(\theta, \phi)$ for escaped rays."""

    cfunc_type = "void"
    name = "calculate_and_fill_blueprint_data_universal"
    params = "const PhotonStateSoA *restrict all_photons, const long int num_rays, blueprint_data_t *restrict result"

    # Python: Generate the loop body for the streaming bundle architecture.
    # Note: Includes CRITICAL FIX to pre-load results from Host to Device.
    loop_body = f"""
    // Variable current_chunk_size defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY);

    // --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---
    // Transfer the status array for the current bundle.
    cudaMemcpy(d_status_bundle, all_photons->status + start_idx,
               sizeof(termination_type_t) * current_chunk_size, cudaMemcpyHostToDevice);

    // Transfer the 9-component state vector $f^\mu$ for the current bundle.
    for(int m=0; m<9; m++) {{
        cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
                   all_photons->f + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice);
    }}

    // Pre-load existing results from Host to Device.
    // This prevents valid Window/Source hits from being overwritten with garbage by the D2H copy.
    cudaMemcpy(d_result_bundle, result + start_idx,
               sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyHostToDevice);

    // --- KERNEL LAUNCH ---
    {launch_body}

    // --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---
    // Retrieve the calculated blueprint results for the current bundle.
    cudaMemcpy(result + start_idx, d_result_bundle,
               sizeof(blueprint_data_t) * current_chunk_size, cudaMemcpyDeviceToHost);
    """

    host_loop = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        loop_body=loop_body,
    )

    body = f"""
    // --- VRAM STAGING ALLOCATION ---
    // Device pointers for the bundled data processing.
    double *d_f_bundle;
    termination_type_t *d_status_bundle;
    blueprint_data_t *d_result_bundle;

    // Hardware Justification: Allocate buffers sized to BUNDLE_CAPACITY to fit within 10GB VRAM.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_status_bundle, sizeof(termination_type_t) * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_result_bundle, sizeof(blueprint_data_t) * BUNDLE_CAPACITY);

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop}

    // --- VRAM CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_status_bundle);
    BHAH_FREE_DEVICE(d_result_bundle);
    """

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
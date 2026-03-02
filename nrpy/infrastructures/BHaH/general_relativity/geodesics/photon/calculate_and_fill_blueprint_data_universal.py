"""
Generates the CUDA kernel and launch code to finalize photon blueprint data.

This module extracts the final state from the Structure of Arrays (SoA)
and populates the results buffer, specifically projecting escaped trajectories
onto the celestial sphere. It couples the batch-processing physics engine 
with the final data extraction, ensuring all termination conditions are correctly logged.
Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def calculate_and_fill_blueprint_data_universal() -> None:
    """
    Generate the C function and native kernel for universal blueprint data finalization.

    :param None: This function requires no arguments.
    :raises Exception: Propagates nrpy core exceptions on generation failure.
    """
    kernel_name = "calculate_and_fill_blueprint_data_universal_kernel"
    
    arg_dict_cuda = {
        "all_photons": "const PhotonStateSoA *restrict",
        "num_rays": "const long int",
        "result": "blueprint_data_t *restrict",
    }

    kernel_body = r"""
    // --- THREAD IDENTIFICATION ---
    // Global 1D thread mapping for photon batch processing.
    const long int photon_idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid memory addresses.
    if (photon_idx >= num_rays) return;

    // --- STATUS SYNCHRONIZATION ---
    // Synchronize final exit status directly into the persistent blueprint array.
    result[photon_idx].termination_type = all_photons->status[photon_idx];

    // --- EVENT DETECTION & TERMINATION CHECKS ---
    // === Termination Dispatch: Celestial Sphere (Escape) ===
    if (all_photons->status[photon_idx] == TERMINATION_TYPE_CELESTIAL_SPHERE) {
        // Unpack final coordinates from the master state vector $f^\mu$.
        // Mapped to thread-local variables to minimize global memory reads.
        const double x_final = all_photons->f[IDX_GLOBAL(1, photon_idx, num_rays)];
        const double y_final = all_photons->f[IDX_GLOBAL(2, photon_idx, num_rays)];
        const double z_final = all_photons->f[IDX_GLOBAL(3, photon_idx, num_rays)];

        // Compute the final radial distance $r = \sqrt{x^2 + y^2 + z^2}$.
        const double r_final = sqrt(x_final*x_final + y_final*y_final + z_final*z_final);

        // --- CELESTIAL PROJECTION ---
        // Map the Cartesian escape coordinates to the celestial sphere to prevent division by zero.
        if (r_final > 1e-9) {
            result[photon_idx].final_theta = acos(z_final / r_final); // Polar angle $\theta$ relative to the $z$-axis.
            result[photon_idx].final_phi = atan2(y_final, x_final);   // Azimuthal angle $\phi$ in the $xy$-plane.
        }
    }
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(num_rays + 256 - 1) / 256", "1", "1"],
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
    1. Computes the global 1D thread index mapping.
    2. Synchronizes the termination status from the master SoA.
    3. For escaped photons, projects the final 3D position onto a celestial sphere $(\theta, \phi)$."""
    
    cfunc_type = "void"
    name = "calculate_and_fill_blueprint_data_universal"
    params = "const PhotonStateSoA *restrict all_photons, const long int num_rays, blueprint_data_t *restrict result"

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=launch_body,
    )
# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py
"""
Defines the C kernel and orchestrator for photon initialization.

This module provides a C function that initializes photon trajectories in Cartesian
coordinates. It allocates a staging buffer to process rays in batches, initializing
the spatial positions, spatial momenta, and adaptive step sizes. It sets the temporal
momentum to zero for downstream constraint solving. The implementation orchestrates
asynchronous data transfers and hardware synchronization.

Single coalesced memory writes prevent thread serialization and ensure aligned cache
access. An explicit hardware error synchronization trap prevents silent link-time
symbol failures caused by compiling with -rdc=true. Hydrating pinned memory via data
bus seeds the Time Slot Manager. Evaluating the initial side of the observer and
source planes natively prevents redundant device memory allocation and data transfers.
Thread identification boundaries prevent out-of-bounds access for threads exceeding
the active chunk. Parallelized batch processing distributes execution across threads.
Processing memory in static bundles protects hardware limits. A synchronization
transfer updates the master Structure of Arrays state.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.loop import loop

# Define the C-structs required for the simulation pipeline.
# These are registered here to ensure they appear in BHaH_defines.h before
# the initialization kernel is compiled.
batch_structs_c_code = r"""
    // Maximum number of photons processed per batch to fit within L1/L2 cache.
    #define BUNDLE_CAPACITY 524288

    // Defines the physical planes where a photon trajectory might terminate.
    typedef enum {
        WINDOW_EVENT, // Intersection with the observer's camera window.
        SOURCE_EVENT  // Intersection with the emission source plane.
    } event_type_t; // END ENUM: event_type_t

    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        TERMINATION_TYPE_CELESTIAL_SPHERE, // 0: Photon escaped to infinity.
        TERMINATION_TYPE_SOURCE_PLANE,     // 1: Photon successfully hit the source emission plane.
        FAILURE_PT_TOO_BIG,                // 2: Integration failed due to unbounded momentum $p_t$.
        FAILURE_RKF45_REJECTION_LIMIT,     // 3: Adaptive step-size rejected too many consecutive times.
        FAILURE_T_MAX_EXCEEDED,            // 4: Integration exceeded maximum allowable physical coordinate time $t$.
        FAILURE_SLOT_MANAGER_ERROR,        // 5: TimeSlotManager failed to allocate or retrieve the photon.
        TERMINATION_TYPE_FAILURE,          // 6: Generic unclassified numerical failure.
        ACTIVE,                            // 7: Photon is currently undergoing integration.
        REJECTED,                          // 8: Photon RKF45 step was rejected.
    } termination_type_t; // END ENUM: termination_type_t

    // Stores the final physical properties of a photon upon integration termination.
    typedef struct {
        termination_type_t termination_type; // The exit condition of the photon.
        double y_w; // Local $y$-coordinate intersection on the observer window.
        double z_w; // Local $z$-coordinate intersection on the observer window.
        double y_s; // Local $y$-coordinate intersection on the source plane.
        double z_s; // Local $z$-coordinate intersection on the source plane.
        double final_theta; // Final polar angle $\theta$ at termination.
        double final_phi;   // Final azimuthal angle $\phi$ at termination.
        double L_w; // Affine parameter $\lambda$ at the window intersection.
        double t_w; // Physical coordinate time $t$ at the window intersection.
        double L_s; // Affine parameter $\lambda$ at the source intersection.
        double t_s; // Physical coordinate time $t$ at the source intersection.
    } __attribute__((packed)) blueprint_data_t; // END STRUCT: blueprint_data_t

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; // Flattened state vector mapping 9 components $t, x, y, z, p_t, p_x, p_y, p_z, \text{aux}$.
        double *f_p; // State vector at the previous integration step.
        double *f_p_p; // State vector at two integration steps prior.
        double *affine_param; // Current affine parameter $\lambda$ for the trajectory.
        double *affine_param_p; // Affine parameter $\lambda$ at the previous step.
        double *affine_param_p_p; // Affine parameter $\lambda$ at two steps prior.
        double *h; // Current adaptive step size $h$ for the RKF45 integrator.
        termination_type_t *status; // Current physical/numerical status of the photon.
        int *rejection_retries; // Counter for consecutive RKF45 error tolerance rejections.

        // Event Detection State Flags (Persistence Layer for Batch C)
        bool *on_positive_side_of_window_prev; // True if photon was previously 'above' the window plane.
        bool *on_positive_side_of_source_prev; // True if photon was previously 'above' the source plane.

        bool *source_event_found; // Flag indicating a source plane intersection was detected.
        double *source_event_lambda; // Exact affine parameter $\lambda$ of the source intersection.
        double *source_event_f_intersect; // Interpolated 9-component state vector at the source intersection.

        bool *window_event_found; // Flag indicating an observer window intersection was detected.
        double *window_event_lambda; // Exact affine parameter $\lambda$ of the window intersection.
        double *window_event_f_intersect; // Interpolated 9-component state vector at the window intersection.
    } PhotonStateSoA; // END STRUCT: PhotonStateSoA
"""
Bdefines_h.register_BHaH_defines("photon_batch_structs", batch_structs_c_code)


def set_initial_conditions_kernel(spacetime_name: str) -> None:
    """
    Register the C function and device kernel for Cartesian photon initialization.

    :param spacetime_name: The specific metric or spacetime identifier (e.g., 'KerrSchild').
    """
    # Register necessary global parameters
    par.register_CodeParameter(
        "int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "window_center_x",
            "window_center_y",
            "window_center_z",
        ],
        [50.0, 0.0, 0.0],
        commondata=True,
        add_to_parfile=False,
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "original_window_center_x",
            "original_window_center_y",
            "original_window_center_z",
            "camera_pos_x",
            "camera_pos_y",
            "camera_pos_z",
            "window_up_vec_x",
            "window_up_vec_y",
            "window_up_vec_z",
            "window_width",
            "window_height",
            "t_start",
        ],
        [50.0, 0.0, 0.0, 51.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 100.0],
        commondata=True,
        add_to_parfile=True,
    )

    # Dynamic architecture detection.
    parallelization = par.parval_from_str("parallelization")
    cd_access = parallel_utils.get_commondata_access(parallelization)

    # Dictionary mapping for GPU/CPU kernel arguments.
    arg_dict = {
        "num_rays": "const long int",
        "d_f_bundle": "double *restrict",
        "d_h_bundle": "double *restrict",
        "nx_0": "const double",
        "nx_1": "const double",
        "nx_2": "const double",
        "ny_0": "const double",
        "ny_1": "const double",
        "ny_2": "const double",
        "start_idx": "const long int",
        "chunk_size": "const long int",
    }
    # Pass commondata explicitly when not using CUDA's global memory
    if parallelization != "cuda":
        arg_dict["commondata"] = "const commondata_struct *restrict"

    # ==========================================
    # ARCHITECTURE-SPECIFIC KERNEL PREAMBLE/POSTAMBLE
    # ==========================================
    if parallelization == "cuda":
        loop_preamble = r"""
    //==========================================
    // THREAD IDENTIFICATION & BOUNDARY CHECKS
    //==========================================
    // Thread ID maps to a unique photon index within the current bundle batch via the identifier $c$.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;

    if (c >= chunk_size) return;
    """
        loop_postamble = ""
    else:
        loop_preamble = r"""
    //==========================================
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    #pragma omp parallel for
    for (long int c = 0; c < chunk_size; c++) {
    """
        loop_postamble = "} // END LOOP: for c over chunk_size"

    # ==========================================
    # CORE MATH (Hardware Agnostic)
    # ==========================================
    core_math = r"""
    // The identifier $i$ represents the global ray index within the master $num\_rays$ SoA.
    const long int i = start_idx + c;

    //==========================================
    // MACRO DEFINITIONS FOR BUNDLE ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state bundle using SoA layout aligned to the active BUNDLE_CAPACITY.
    #define IDX_F(comp, ray_id) ((comp) * BUNDLE_CAPACITY + (ray_id))
    // IDX_H maps to the 1D adaptive step size bundle.
    #define IDX_H(ray_id) (ray_id)

    //==========================================
    // PIXEL MAPPING & GEOMETRY
    //==========================================
    // Vertical pixel coordinate index within the virtual observer's projection frame.
    const int row = i / {cd_access}scan_density;
    // Horizontal pixel coordinate index within the virtual observer's projection frame.
    const int col = i % {cd_access}scan_density;

    // Local physical distance $x_{pix}$ along the horizontal camera axis.
    const double x_pix = -{cd_access}window_width/2.0 + (col + 0.5) * ({cd_access}window_width / {cd_access}scan_density);
    // Local physical distance $y_{pix}$ along the vertical camera axis.
    const double y_pix = -{cd_access}window_height/2.0 + (row + 0.5) * ({cd_access}window_height / {cd_access}scan_density);

    // The array $target\_pos$ stores the global Cartesian intersection point $x^\mu$ on the projection window.
    const double target_pos[3] = {
        {cd_access}window_center_x + x_pix*nx_0 + y_pix*ny_0,
        {cd_access}window_center_y + x_pix*nx_1 + y_pix*ny_1,
        {cd_access}window_center_z + x_pix*nx_2 + y_pix*ny_2
    };

    //==========================================
    // INITIAL STATE POPULATION
    //==========================================
    // Write the starting position and spatial momentum explicitly to the VRAM bundle.
    d_f_bundle[IDX_F(0, c)] = {cd_access}t_start; // Coordinate time $t$
    d_f_bundle[IDX_F(1, c)] = {cd_access}camera_pos_x; // Spatial position $x$
    d_f_bundle[IDX_F(2, c)] = {cd_access}camera_pos_y; // Spatial position $y$
    d_f_bundle[IDX_F(3, c)] = {cd_access}camera_pos_z; // Spatial position $z$

    // Vector component $V^x$ for the unnormalized geometric trajectory.
    const double V_x = target_pos[0] - {cd_access}camera_pos_x;
    // Vector component $V^y$ for the unnormalized geometric trajectory.
    const double V_y = target_pos[1] - {cd_access}camera_pos_y;
    // Vector component $V^z$ for the unnormalized geometric trajectory.
    const double V_z = target_pos[2] - {cd_access}camera_pos_z;

    // Normalization ensures initial 4-momentum $p^\mu$ satisfies null trajectory constraints.
    const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z);

    // Explicitly set the temporal momentum $p^t = 0$ for downstream Hamiltonian constraint solving.
    d_f_bundle[IDX_F(4, c)] = 0.0;

    d_f_bundle[IDX_F(5, c)] = V_x * inv_mag_V; // Initial momentum component $p^x$
    d_f_bundle[IDX_F(6, c)] = V_y * inv_mag_V; // Initial momentum component $p^y$
    d_f_bundle[IDX_F(7, c)] = V_z * inv_mag_V; // Initial momentum component $p^z$

    // Explicitly set the initial path length to zero.
    d_f_bundle[IDX_F(8, c)] = 0.0;

    // Initialize the adaptive step size $h$ for the RKF45 integrator.
    d_h_bundle[IDX_H(c)] = {cd_access}numerical_initial_h;

    //==========================================
    // MACRO CLEANUP
    //==========================================
    #undef IDX_F
    #undef IDX_H
    """.replace("{cd_access}", cd_access)

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = (
        {
            "threads_per_block": ["BHAH_THREADS_IN_X_DIR_DEFAULT", "1", "1"],
            "blocks_per_grid": [
                "(chunk_size + BHAH_THREADS_IN_X_DIR_DEFAULT - 1) / BHAH_THREADS_IN_X_DIR_DEFAULT",
                "1",
                "1",
            ],
            "stream": "stream_idx",
        }
        if parallelization == "cuda"
        else None
    )

    kernel_prefunc, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name=f"set_initial_conditions_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    # ==========================================
    # MEMORY BRIDGE ARCHITECTURE
    # ==========================================
    if parallelization == "cuda":
        sync_and_transfer_code = r"""
        //==========================================
        // EXPLICIT HARDWARE ERROR SYNCHRONIZATION
        //==========================================
        #ifdef DEBUG
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            // Fatal unrecoverable error during kernel synchronization.
            printf("Init Kernel Failed on Batch starting at %ld: %s\\n", (long int)start_idx, cudaGetErrorString(err));
            exit(1);
        } // END IF: check for kernel launch errors
        #endif

        //==========================================
        // 9-STRIDED BRIDGE TRANSFER (DEVICE-TO-HOST)
        //==========================================
        // Transfer initialized state vectors $f^\mu$ from VRAM back to host RAM.
        for(int m=0; m<9; m++) {
            cudaMemcpy(all_photons->f + (m * num_rays) + start_idx,
                    d_f_bundle + (m * BUNDLE_CAPACITY),
                    sizeof(double) * chunk_size,
                    cudaMemcpyDeviceToHost);
        } // END LOOP: for m over 9-strided tensor transfer

        //==========================================
        // 1-STRIDED BRIDGE TRANSFER (DEVICE-TO-HOST)
        //==========================================
        // Transfer initialized adaptive step sizes $h$ from VRAM back to host RAM.
        cudaMemcpy(all_photons->h + start_idx,
                   d_h_bundle,
                   sizeof(double) * chunk_size,
                   cudaMemcpyDeviceToHost);
        """
    else:
        sync_and_transfer_code = r"""
        //==========================================
        // 9-STRIDED BRIDGE TRANSFER (HOST-TO-HOST)
        //==========================================
        // Transfer initialized state vectors $f^\mu$ from staging buffer to master SoA.
        for(int m=0; m<9; m++) {
            memcpy(all_photons->f + (m * num_rays) + start_idx,
                   d_f_bundle + (m * BUNDLE_CAPACITY),
                   sizeof(double) * chunk_size);
        } // END LOOP: for m over 9-strided tensor transfer

        //==========================================
        // 1-STRIDED BRIDGE TRANSFER (HOST-TO-HOST)
        //==========================================
        // Transfer initialized adaptive step sizes $h$ from staging buffer to master SoA.
        memcpy(all_photons->h + start_idx,
               d_h_bundle,
               sizeof(double) * chunk_size);
        """

    # Generate the host-side loop string to iterate over the dataset in bundles.
    loop_body = f"""
        // Variable chunk_size defines the active range for the current streaming bundle.
        const long int chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY);

        {launch_code}

        {sync_and_transfer_code}
    """

    host_loop_code = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        loop_body=loop_body,
    )

    body = r"""
    //==========================================
    // HOST-SIDE GEOMETRY SETUP
    //==========================================
    // Pre-calculate the projection plane basis vectors to save device registers.
    // Calculations use the static original window center to maintain consistent projection framing across all local tiles.
    const double cam_x = commondata->camera_pos_x; // The $x$-coordinate of the camera.
    const double cam_y = commondata->camera_pos_y; // The $y$-coordinate of the camera.
    const double cam_z = commondata->camera_pos_z; // The $z$-coordinate of the camera.

    const double owc_x = commondata->original_window_center_x; // Original $x$-coordinate of the master window center.
    const double owc_y = commondata->original_window_center_y; // Original $y$-coordinate of the master window center.
    const double owc_z = commondata->original_window_center_z; // Original $z$-coordinate of the master window center.

    const double wc_x = commondata->window_center_x; // The $x$-coordinate of the shifted tile window center.
    const double wc_y = commondata->window_center_y; // The $y$-coordinate of the shifted tile window center.
    const double wc_z = commondata->window_center_z; // The $z$-coordinate of the shifted tile window center.

    // Vector $n_z^i$ normal to the camera window representing the global line of sight.
    double n_z[3] = {owc_x - cam_x, owc_y - cam_y, owc_z - cam_z};
    // Magnitude of the normal vector $n_z^i$.
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]);
    // Normalize the $n_z^i$ vector.
    for(int j=0; j<3; j++) n_z[j] /= mag_n_z;

    // Geometric reference vector defining the upward orientation.
    const double guide_up[3] = {commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z};

    // Basis vector $n_x^i$ describing the horizontal camera axis.
    double n_x[3] = {guide_up[1]*n_z[2] - guide_up[2]*n_z[1],
                     guide_up[2]*n_z[0] - guide_up[0]*n_z[2],
                     guide_up[0]*n_z[1] - guide_up[1]*n_z[0]};
    // Magnitude of the horizontal basis vector $n_x^i$.
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);

    // Fallback logic prevents cross-product singularities at geometric poles for numerical stability.
    if (mag_n_x < 1e-9) {
        // Alternative upward reference vector to avoid cross-product singularities.
        double alternative_up[3] = {0.0, 1.0, 0.0};
        if (fabs(n_z[1]) > 0.999) { alternative_up[1] = 0.0; alternative_up[2] = 1.0; }
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    } // END IF: near-nadir fallback

    // Normalize the $n_x^i$ vector.
    for(int j=0; j<3; j++) n_x[j] /= mag_n_x;

    // Basis vector $n_y^i$ describing the vertical camera axis.
    double n_y[3] = {n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]};

    // Output the local tile components.
    window_center_out[0] = wc_x;
    window_center_out[1] = wc_y;
    window_center_out[2] = wc_z;
    for(int j=0; j<3; j++) {
        // Output the normalized global basis vectors $n_x^i, n_y^i, n_z^i$.
        n_x_out[j] = n_x[j]; n_y_out[j] = n_y[j]; n_z_out[j] = n_z[j];
    } // END LOOP: for j over basis vector components

    // Extracted basis components.
    const double nx_0 = n_x[0];
    const double nx_1 = n_x[1];
    const double nx_2 = n_x[2];
    const double ny_0 = n_y[0];
    const double ny_1 = n_y[1];
    const double ny_2 = n_y[2];

    //==========================================
    // DYNAMIC GEOMETRIC PLANE INITIALIZATION
    //==========================================
    // Evaluate the initial side of the observer and source planes natively on the CPU.
    // Evaluates against the original global window to correctly flag observer window crossings regardless of tile offset.
    const double val_window = n_z[0] * (cam_x - owc_x) +
                              n_z[1] * (cam_y - owc_y) +
                              n_z[2] * (cam_z - owc_z);

    // Evaluates strictly greater than zero per the physical observer constraints.
    const bool init_window_side = (val_window > 0.0);

    // Initial dot product to evaluate proximity to the source plane.
    const double val_source = commondata->source_plane_normal_x * (cam_x - commondata->source_plane_center_x) +
                              commondata->source_plane_normal_y * (cam_y - commondata->source_plane_center_y) +
                              commondata->source_plane_normal_z * (cam_z - commondata->source_plane_center_z);

    // Evaluates true if exactly zero per physical emission constraints.
    const bool init_source_side = (val_source >= 0.0);

    // Loop iterator traversing the entire global ray count to hydrate initial plane boundaries.
    for (long int plane_i = 0; plane_i < num_rays; ++plane_i) {
        all_photons->on_positive_side_of_window_prev[plane_i] = init_window_side;
        all_photons->on_positive_side_of_source_prev[plane_i] = init_source_side;
    } // END LOOP: for plane_i over initial plane boundaries

    //==========================================
    // VRAM STAGING ALLOCATION
    //==========================================
    // Device pointer for the chunked VRAM state staging buffer d_f_bundle.
    double *d_f_bundle;
    // Device pointer for the chunked VRAM step size staging buffer d_h_bundle.
    double *d_h_bundle;

    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_h_bundle, sizeof(double) * BUNDLE_CAPACITY);

    //==========================================
    // HOST-SIDE PAGINATION LOOP
    //==========================================
    {host_loop_code}

    //==========================================
    // VRAM DEALLOCATION
    //==========================================
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_h_bundle);
    """.replace("{host_loop_code}", str(host_loop_code))

    # Establish the final strings to satisfy the Translation Unit Inlining Mandate.
    prefunc = f"{kernel_prefunc}"

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r""" Initializes Cartesian starting conditions for photons.

    Detailed algorithm: Maps global thread IDs to pixel coordinates to calculate initial
    spatial coordinates $x, y, z$ and unnormalized momenta $p^x, p^y, p^z$. Explicitly enforces
    temporal momentum $p^t = 0$ and affine parameter $\lambda = 0$ for downstream constraint solvers.

    @param commondata Master configuration struct containing global parameters.
    @param num_rays Total number of photon trajectories in the simulation.
    @param all_photons Pointer to the master SoA state vector $f^\mu$ in host memory.
    @param window_center_out Output buffer for the geometric window center.
    @param n_x_out Output buffer for the $x$-axis basis vector $n_x^i$.
    @param n_y_out Output buffer for the $y$-axis basis vector $n_y^i$.
    @param n_z_out Output buffer for the $z$-axis basis vector $n_z^i$."""

    cfunc_type = "void"
    name = f"set_initial_conditions_kernel_{spacetime_name}"

    stream_param = "const int stream_idx" if parallelization == "cuda" else ""
    params = (
        "const commondata_struct *restrict commondata, "
        "long int num_rays, "
        "PhotonStateSoA *restrict all_photons, "
        "double window_center_out[3], "
        "double n_x_out[3], "
        "double n_y_out[3], "
        "double n_z_out[3]"
    )
    if stream_param:
        params += f", {stream_param}"

    include_CodeParameters_h = False

    # Register the complete C function
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

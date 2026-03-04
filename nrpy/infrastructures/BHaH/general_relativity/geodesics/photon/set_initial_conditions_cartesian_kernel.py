"""
Module for generating the native CUDA kernel and host-side orchestrator for photon initialization.

This module generates a C function that initializes photon trajectories in Cartesian coordinates
using a "Split-Pipeline" architecture. It allocates a VRAM staging buffer to process rays
in batches, initializing the spatial positions, spatial momenta, and adaptive step sizes,
while explicitly setting the temporal momentum to zero for downstream constraint solving.

Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
    get_commondata_access,
)
from nrpy.helpers.loop import loop

# Python: Define the C-structs required for the simulation pipeline.
# These are registered here to ensure they appear in BHaH_defines.h before
# the initialization kernel is compiled.
batch_structs_c_code = r"""
    #define BUNDLE_CAPACITY 32768 // Maximum number of photons processed per batch to fit within L1/L2 cache.

    // Defines the physical planes where a photon trajectory might terminate.
    typedef enum {
        WINDOW_EVENT, // Intersection with the observer's camera window.
        SOURCE_EVENT  // Intersection with the emission source plane.
    } event_type_t;

    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        TERMINATION_TYPE_CELESTIAL_SPHERE, // 0: Photon escaped to infinity (exceeded r_escape).
        TERMINATION_TYPE_SOURCE_PLANE,     // 1: Photon successfully hit the source emission plane.
        FAILURE_PT_TOO_BIG,                // 2: Integration failed due to unbounded momentum $p_t$.
        FAILURE_RKF45_REJECTION_LIMIT,     // 3: Adaptive step-size rejected too many consecutive times.
        FAILURE_T_MAX_EXCEEDED,            // 4: Integration exceeded maximum allowable physical time.
        FAILURE_SLOT_MANAGER_ERROR,        // 5: TimeSlotManager failed to allocate or retrieve the photon.
        TERMINATION_TYPE_FAILURE,          // 6: Generic unclassified numerical failure.
        ACTIVE                             // 7: Photon is currently undergoing integration.
    } termination_type_t;

    // Stores the final physical properties of a photon upon integration termination.
    typedef struct {
        termination_type_t termination_type; // The exit condition of the photon.
        double y_w; // Local y-coordinate intersection on the observer window.
        double z_w; // Local z-coordinate intersection on the observer window.
        double y_s; // Local y-coordinate intersection on the source plane.
        double z_s; // Local z-coordinate intersection on the source plane.
        double final_theta; // Final polar angle $\theta$ at termination.
        double final_phi;   // Final azimuthal angle $\phi$ at termination.
        double L_w; // Affine parameter $\lambda$ at the window intersection.
        double t_w; // Physical coordinate time $t$ at the window intersection.
        double L_s; // Affine parameter $\lambda$ at the source intersection.
        double t_s; // Physical coordinate time $t$ at the source intersection.
    } __attribute__((packed)) blueprint_data_t;

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; // Flattened state vector: 9 components (t, x, y, z, p_t, p_x, p_y, p_z, aux).
        double *f_p; // State vector at the previous integration step.
        double *f_p_p; // State vector at two integration steps prior.
        double *affine_param; // Current affine parameter $\lambda$ for the trajectory.
        double *affine_param_p; // Affine parameter $\lambda$ at the previous step.
        double *affine_param_p_p; // Affine parameter $\lambda$ at two steps prior.
        double *h; // Current adaptive step size for the RKF45 integrator.
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
    } PhotonStateSoA;
"""
Bdefines_h.register_BHaH_defines("photon_02_batch_structs", batch_structs_c_code)

def set_initial_conditions_kernel(spacetime_name: str) -> None:
    """
    Registers the C function and device kernel for Cartesian photon initialization.

    :param spacetime_name: The specific metric or spacetime identifier (e.g., 'KerrSchild').
    :raises Exception: Propagates NRPy+ core exceptions during AST generation.
    """
    # Python: Register necessary global parameters for the grid setup and numerical controls.
    par.register_CodeParameter("int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "t_start", 100.0, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "h_init", 0.1, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "window_width", 10.0, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "window_height", 10.0, commondata=True, add_to_parfile=True)

    # Python: Dictionary mapping for GPU kernel arguments.
    # Note: Global parameters like camera position are accessed directly via Constant Memory.
    arg_dict_cuda = {
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

    # Python: Retrieve the correct accessor for constant memory (e.g., "d_commondata.").
    cd_access = get_commondata_access("cuda")

    # Python: Define the GPU kernel body utilizing raw strings for LaTeX and strict variable documentation.
    kernel_body = rf"""
    // --- THREAD IDENTIFICATION & BOUNDARY CHECKS ---
    // The identifier $c$ represents the local thread index within the current bundle batch.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;

    // Hardware Justification: Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (c >= chunk_size) return;

    // The identifier $i$ represents the global ray index within the master $num\_rays$ SoA.
    const long int i = start_idx + c;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout aligned to the active chunk_size.
    #define IDX_F(comp, ray_id) ((comp) * chunk_size + (ray_id))
    // IDX_H maps to the 1D adaptive step size bundle.
    #define IDX_H(ray_id) (ray_id)

    // --- PIXEL MAPPING & GEOMETRY ---
    // Vertical pixel coordinate index within the virtual observer's projection frame.
    const int row = i / {cd_access}scan_density;
    // Horizontal pixel coordinate index within the virtual observer's projection frame.
    const int col = i % {cd_access}scan_density;

    // The variable $x_{{pix}}$ is the local physical distance along the horizontal camera axis.
    const double x_pix = -{cd_access}window_width/2.0 + (col + 0.5) * ({cd_access}window_width / {cd_access}scan_density);
    // The variable $y_{{pix}}$ is the local physical distance along the vertical camera axis.
    const double y_pix = -{cd_access}window_height/2.0 + (row + 0.5) * ({cd_access}window_height / {cd_access}scan_density);

    // The array $target\_pos$ stores the global Cartesian intersection point on the projection window $x^\mu$.
    const double target_pos[3] = {{
        {cd_access}window_center_x + x_pix*nx_0 + y_pix*ny_0,
        {cd_access}window_center_y + x_pix*nx_1 + y_pix*ny_1,
        {cd_access}window_center_z + x_pix*nx_2 + y_pix*ny_2
    }};

    // --- INITIAL STATE POPULATION ---
    // Algorithmic Step: Write the starting position and spatial momentum explicitly to the VRAM bundle.
    // Hardware Justification: Single coalesced writes via WriteCUDA prevent warp serialization on sm_86.
    
    WriteCUDA(&d_f_bundle[IDX_F(0, c)], {cd_access}t_start); // Coordinate time $t$
    WriteCUDA(&d_f_bundle[IDX_F(1, c)], {cd_access}camera_pos_x); // Spatial position $x$
    WriteCUDA(&d_f_bundle[IDX_F(2, c)], {cd_access}camera_pos_y); // Spatial position $y$
    WriteCUDA(&d_f_bundle[IDX_F(3, c)], {cd_access}camera_pos_z); // Spatial position $z$

    // Vector component $V^x$ for the unnormalized geometric trajectory.
    const double V_x = target_pos[0] - {cd_access}camera_pos_x;
    // Vector component $V^y$ for the unnormalized geometric trajectory.
    const double V_y = target_pos[1] - {cd_access}camera_pos_y;
    // Vector component $V^z$ for the unnormalized geometric trajectory.
    const double V_z = target_pos[2] - {cd_access}camera_pos_z;

    // Functional Justification: Normalization ensures initial 4-momentum $p^\mu$ satisfies null trajectory constraints.
    const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z);

    // Explicitly set the temporal momentum $p^t = 0$ for downstream Hamiltonian constraint solving.
    WriteCUDA(&d_f_bundle[IDX_F(4, c)], 0.0);

    WriteCUDA(&d_f_bundle[IDX_F(5, c)], V_x * inv_mag_V); // Initial momentum component $p^x$
    WriteCUDA(&d_f_bundle[IDX_F(6, c)], V_y * inv_mag_V); // Initial momentum component $p^y$
    WriteCUDA(&d_f_bundle[IDX_F(7, c)], V_z * inv_mag_V); // Initial momentum component $p^z$
    
    // Explicitly set the intitial path length to zero. 
    // See geodesics.py in /nrpy/equations/general_relativity/geodesics to define path length.
    WriteCUDA(&d_f_bundle[IDX_F(8, c)], 0.0);

    // Initialize the adaptive step size $h$ for the RKF45 integrator.
    WriteCUDA(&d_h_bundle[IDX_H(c)], {cd_access}h_init);

    // --- MACRO CLEANUP ---
    #undef IDX_F
    #undef IDX_H
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
    }

    kernel_prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name=f"set_initial_conditions_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    # Python: Generate the host-side loop string to iterate over the dataset in bundles.
    # We use 'start_idx' as the iterator and 'BUNDLE_CAPACITY' as the stride.
    loop_body = f"""
    // Variable chunk_size defines the active range for the current streaming bundle.
    const long int chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY);

    {launch_code}

    // --- EXPLICIT HARDWARE ERROR SYNCHRONIZATION ---
    // Algorithmic Step: Catch unrecoverable device faults immediately after kernel execution.
    // Fatal unrecoverable error exit if the device execution fails.
    #ifdef DEBUG
    cudaDeviceSynchronize();
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {{
        printf("Init Kernel Failed: %s\\n", cudaGetErrorString(err));
        exit(1);
    }}
    #endif

    // --- 9-STRIDED BRIDGE TRANSFER (DEVICE-TO-HOST) ---
    // Algorithmic Step: Transfer initialized state vectors $f^\mu$ from VRAM back to host RAM.
    // Hardware Justification: Pinned memory on the host is hydrated via PCIe to seed the Time Slot Manager.
    for(int m=0; m<9; m++) {{
        cudaMemcpy(all_photons->f + (m * num_rays) + start_idx,
                   d_f_bundle + (m * chunk_size),
                   sizeof(double) * chunk_size,
                   cudaMemcpyDeviceToHost);
    }}

    // --- 1-STRIDED BRIDGE TRANSFER (DEVICE-TO-HOST) ---
    // Algorithmic Step: Transfer initialized adaptive step sizes $h$ from VRAM back to host RAM.
    // Hardware Justification: Device-to-Host transfer required to synchronize the master SoA state.
    cudaMemcpy(all_photons->h + start_idx,
               d_h_bundle,
               sizeof(double) * chunk_size,
               cudaMemcpyDeviceToHost);
    """
    
    host_loop_code = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        loop_body=loop_body,
    )

    body = rf"""
    // --- HOST-SIDE GEOMETRY SETUP ---
    // Algorithmic Step: Pre-calculate the projection plane basis vectors to save device registers.
    const double cam_x = commondata->camera_pos_x;
    const double cam_y = commondata->camera_pos_y;
    const double cam_z = commondata->camera_pos_z;

    const double wc_x = commondata->window_center_x;
    const double wc_y = commondata->window_center_y;
    const double wc_z = commondata->window_center_z;

    // Vector $n_z$ normal to the camera window representing the line of sight.
    double n_z[3] = {{wc_x - cam_x, wc_y - cam_y, wc_z - cam_z}};
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]);
    for(int j=0; j<3; j++) n_z[j] /= mag_n_z;

    // Geometric reference vector defining the upward orientation.
    const double guide_up[3] = {{commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}};

    // Basis vector $n_x$ describing the horizontal camera axis.
    double n_x[3] = {{n_z[1]*guide_up[2] - n_z[2]*guide_up[1], n_z[2]*guide_up[0] - n_z[0]*guide_up[2], n_z[0]*guide_up[1] - n_z[1]*guide_up[0]}};
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);

    // Functional Justification: Fallback logic prevents cross-product singularities at geometric poles for numerical stability.
    if (mag_n_x < 1e-9) {{
        double alternative_up[3] = {{0.0, 1.0, 0.0}};
        if (fabs(n_z[1]) > 0.999) {{ alternative_up[1] = 0.0; alternative_up[2] = 1.0; }}
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    }}

    for(int j=0; j<3; j++) n_x[j] /= mag_n_x;

    // Basis vector $n_y$ describing the vertical camera axis.
    double n_y[3] = {{n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]}};

    window_center_out[0] = wc_x;
    window_center_out[1] = wc_y;
    window_center_out[2] = wc_z;
    for(int j=0; j<3; j++) {{
        n_x_out[j] = n_x[j]; n_y_out[j] = n_y[j]; n_z_out[j] = n_z[j];
    }}

    const double nx_0 = n_x[0];
    const double nx_1 = n_x[1];
    const double nx_2 = n_x[2];
    const double ny_0 = n_y[0];
    const double ny_1 = n_y[1];
    const double ny_2 = n_y[2];

    // --- VRAM STAGING ALLOCATION ---
    // Device pointer for the chunked VRAM state staging buffer $d\_f\_bundle$.
    double *d_f_bundle;
    // Device pointer for the chunked VRAM step size staging buffer $d\_h\_bundle$.
    double *d_h_bundle;
    
    // Hardware Justification: Memory is processed in bundles to protect the 10GB hardware limit.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_h_bundle, sizeof(double) * BUNDLE_CAPACITY);

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop_code}

    // --- VRAM DEALLOCATION ---
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_h_bundle);
    """

    # Python: Establish the final strings to satisfy the Translation Unit Inlining Mandate.
    prefunc = f"{kernel_prefunc}"

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "cuda_intrinsics.h",
        "<math.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]

    desc = r"""@brief Initializes Cartesian starting conditions for photons.
    
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
    
    params = (
        "const commondata_struct *restrict commondata, "
        "long int num_rays, "
        "PhotonStateSoA *restrict all_photons, "
        "double window_center_out[3], "
        "double n_x_out[3], "
        "double n_y_out[3], "
        "double n_z_out[3]"
    )

    # Python: Register the complete C function using the canonical Master Order.
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
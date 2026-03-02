"""
Module for generating Cartesian initial conditions for the photon integration pipeline.

This module defines the C function responsible for initializing the spatial coordinates
and spatial momentum of photons on a camera window. It utilizes a dedicated GPU kernel
to execute the initialization directly in VRAM, enforcing thread-local memory access
for the metric tensor to comply with hardware register limits.
Author: Dalton J. Moone.
"""
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def set_initial_conditions_cartesian(spacetime_name: str) -> None:
    """
    Register the C function and device kernel to set up Cartesian initial conditions for photons.

    :param spacetime_name: The specific metric or spacetime identifier (e.g., 'Kerr').
    :raises Exception: Propagates nrpy core exceptions on generation failure.
    """
    par.register_CodeParameter("int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "t_start", 100, commondata=True, add_to_parfile=True)

    arg_dict_cuda = {
        "commondata": "const commondata_struct *restrict",
        "num_rays": "const long int",
        "all_photons": "PhotonStateSoA *restrict",
        "cam_x": "const double", "cam_y": "const double", "cam_z": "const double",
        "wc_x": "const double", "wc_y": "const double", "wc_z": "const double",
        "nx_0": "const double", "nx_1": "const double", "nx_2": "const double",
        "ny_0": "const double", "ny_1": "const double", "ny_2": "const double",
        "start_idx": "const long int",
        "current_chunk_size": "const long int"
    }

    kernel_body = f"""
    // --- THREAD IDENTIFICATION ---
    // Calculate the local thread index within the current bundle batch.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid memory addresses.
    if (c >= current_chunk_size) return;

    // The global index of the photon currently being initialized.
    const long int i = start_idx + c;

    // --- THREAD-LOCAL TENSOR ALLOCATION ---
    // Thread-local array for the 9-element state vector $f^\mu$ and momentum $p^\mu$.
    // Keeps intermediate math within the 255 register limit per thread for sm_86.
    double f_local[9];
    // Thread-local array for the 10 independent metric components $g_{{\mu\nu}}$.
    double metric_local[10];

    // --- PIXEL MAPPING ---
    // The vertical pixel coordinate index derived from the 1D global ray identifier.
    const int row = i / commondata->scan_density;
    // The horizontal pixel coordinate index derived from the 1D global ray identifier.
    const int col = i % commondata->scan_density;

    // Local physical distance along the camera frame's horizontal and vertical axes.
    const double x_pix = -commondata->window_width/2.0 + (col + 0.5) * (commondata->window_width / commondata->scan_density);
    const double y_pix = -commondata->window_height/2.0 + (row + 0.5) * (commondata->window_height / commondata->scan_density);

    // Global Cartesian intersection point on the projection window $x^\mu$.
    const double target_pos[3] = {{
        wc_x + x_pix*nx_0 + y_pix*ny_0,
        wc_y + x_pix*nx_1 + y_pix*ny_1,
        wc_z + x_pix*nx_2 + y_pix*ny_2
    }};

    // --- INITIAL STATE POPULATION ---
    // Populate the thread-local state vector $f^\mu_{{local}}$ with starting coordinates.
    f_local[0] = commondata->t_start;
    f_local[1] = cam_x;
    f_local[2] = cam_y;
    f_local[3] = cam_z;

    // Unnormalized geometric trajectory vector connecting the camera to the target pixel $V^i$.
    const double V_x = target_pos[0] - cam_x;
    const double V_y = target_pos[1] - cam_y;
    const double V_z = target_pos[2] - cam_z;

    // Normalization scalar to generate a unit direction vector.
    const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z);

    f_local[5] = V_x * inv_mag_V;
    f_local[6] = V_y * inv_mag_V;
    f_local[7] = V_z * inv_mag_V;
    f_local[8] = 0.0; // Affine parameter $\lambda$

    // --- INTERPOLATION & METRIC EVALUATION ---
    // Extract local position array to execute the single-thread inline interpolation helper.
    double pos_local[4] = {{f_local[0], f_local[1], f_local[2], f_local[3]}};
    placeholder_interpolation_engine_{spacetime_name}(commondata, pos_local, metric_local, NULL);

    // --- HAMILTONIAN CONSTRAINT ---
    // Solve the negative root of the Hamiltonian constraint $p_\mu p^\mu = 0$ for initial $p^0$.
    p0_reverse(metric_local, f_local, &f_local[4]);

    // --- GLOBAL MEMORY WRITE ---
    // Writes the local state vector to the global SoA $f^\mu$.
    // Coalesced to global memory after local computation to prevent warp divergence.
    for(int m=0; m<9; m++) {{
        all_photons->f[IDX_GLOBAL(m, i, num_rays)] = f_local[m];
    }}
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(current_chunk_size + 256 - 1) / 256", "1", "1"],
    }

    prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name=f"set_initial_conditions_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "math.h",
        "stdio.h",
        "stdlib.h",
    ]

    desc = f"""@brief Initializes Cartesian starting conditions for photons in {spacetime_name}.
    @param commondata The master configuration struct.
    @param num_rays Total photons in the simulation batch.
    @param all_photons Pointer to the master SoA state vector.
    @param window_center_out Output buffer for geometric center.
    @param n_x_out Output buffer for the $x$-axis geometric basis vector.
    @param n_y_out Output buffer for the $y$-axis geometric basis vector.
    @param n_z_out Output buffer for the $z$-axis geometric basis vector.

    Detailed algorithm: Computes the geometric basis vectors of the camera window
    on the host, then launches a standalone device kernel in controlled batches to map 
    individual photon rays to pixel coordinates and populate the initial flattened Structure of Arrays
    (SoA) state vector. Enforces local register boundaries to prevent VRAM overflow."""

    cfunc_type = "void"
    name = f"set_initial_conditions_cartesian_{spacetime_name}"
    
    params = (
        "const commondata_struct *restrict commondata, "
        "long int num_rays, "
        "PhotonStateSoA *restrict all_photons, "
        "double window_center_out[3], "
        "double n_x_out[3], "
        "double n_y_out[3], "
        "double n_z_out[3]"
    )

    body = f"""
    // --- HOST-SIDE GEOMETRY SETUP ---
    // Spatial coordinates of the observer's camera in the global Cartesian frame $x^\mu_{{cam}}$.
    const double cam_x = commondata->camera_pos_x;
    const double cam_y = commondata->camera_pos_y;
    const double cam_z = commondata->camera_pos_z;

    // Geometric center of the projection window in the global Cartesian frame.
    const double wc_x = commondata->window_center_x;
    const double wc_y = commondata->window_center_y;
    const double wc_z = commondata->window_center_z;

    // Vector normal to the camera window (line of sight direction).
    double n_z[3] = {{wc_x - cam_x, wc_y - cam_y, wc_z - cam_z}};
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]);
    for(int i=0; i<3; i++) n_z[i] /= mag_n_z;

    const double guide_up[3] = {{commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}};
    
    // Orthonormal basis vector $n_x^i$ describing the horizontal axis of the camera frame.
    double n_x[3] = {{n_z[1]*guide_up[2] - n_z[2]*guide_up[1], n_z[2]*guide_up[0] - n_z[0]*guide_up[2], n_z[0]*guide_up[1] - n_z[1]*guide_up[0]}};
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);

    // Fallback mechanism to prevent cross-product singularities at the poles.
    if (mag_n_x < 1e-9) {{
        double alternative_up[3] = {{0.0, 1.0, 0.0}};
        if (fabs(n_z[1]) > 0.999) {{ alternative_up[1] = 0.0; alternative_up[2] = 1.0; }}
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    }}

    for(int i=0; i<3; i++) n_x[i] /= mag_n_x;

    // Orthonormal basis vector $n_y^i$ describing the vertical axis of the camera frame.
    double n_y[3] = {{n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]}};

    for(int i=0; i<3; i++) {{
        window_center_out[i] = (&wc_x)[i];
        n_x_out[i] = n_x[i]; n_y_out[i] = n_y[i]; n_z_out[i] = n_z[i];
    }}

    // Maps local basis variables to kernel parameters to prevent host-to-device memory copies.
    const double nx_0 = n_x[0]; const double nx_1 = n_x[1]; const double nx_2 = n_x[2];
    const double ny_0 = n_y[0]; const double ny_1 = n_y[1]; const double ny_2 = n_y[2];

    // --- BATCH PROCESSING & KERNEL LAUNCH ---
    // Launch the CUDA kernel to populate the initial state vector directly in VRAM.
    // BUNDLE_CAPACITY limits the streaming chunk size to prevent PCIe and VRAM overflow.
    for (long int start_idx = 0; start_idx < num_rays; start_idx += BUNDLE_CAPACITY) {{
        const long int current_chunk_size = (num_rays - start_idx < BUNDLE_CAPACITY) ? (num_rays - start_idx) : BUNDLE_CAPACITY;
        
        {launch_code}
    }}
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
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

    :param spacetime_name: The specific metric or spacetime identifier.
    :raises Exception: Propagates nrpy core exceptions on generation failure.
    """
    # Python: Register necessary global parameters for the grid setup.
    par.register_CodeParameter("int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "t_start", 100, commondata=True, add_to_parfile=True)

    # Python: Dictionary mapping for GPU kernel arguments.
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

    # Python: Define the GPU kernel body utilizing raw strings to protect LaTeX comments.
    kernel_body = rf"""
    // --- THREAD IDENTIFICATION ---
    /*
     * Algorithmic Step: Calculate the local thread index within the current bundle batch.
     * Architectural Justification: Maps execution tightly to the 1D streaming bundle grid,
     * ensuring each thread uniquely handles a single independent photon trajectory.
     */
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Local thread index within batch block.

    // --- BOUNDARY CHECK ---
    /*
     * Algorithmic Step: Ensure out-of-bounds threads do not access invalid memory addresses.
     * Architectural Justification: Truncates execution for boundary blocks to prevent warp divergence
     * and segmentation faults when the batch size is not perfectly divisible by the block size.
     */
    if (c >= current_chunk_size) return;

    const long int i = start_idx + c; // Global index of the photon being initialized.

    // --- THREAD-LOCAL TENSOR ALLOCATION ---
    /*
     * Algorithmic Step: Allocate purely thread-local arrays for physical state and metric data.
     * Architectural Justification: Restricting allocations to local arrays keeps intermediate math
     * strictly within the 255 register limit per thread for the sm_86 architecture, preventing VRAM spillage.
     */
    double f_local[9]; // Thread-local array for the 9-element state vector $f^\mu$ and momentum $p^\mu$.
    double metric_local[10]; // Thread-local array for the 10 independent metric components $g_{{\mu\nu}}$.

    // --- PIXEL MAPPING ---
    /*
     * Algorithmic Step: Map the 1D global ray identifier to a 2D camera pixel coordinate.
     * Architectural Justification: Integer arithmetic mapping provides deterministic sub-pixel
     * sampling without relying on floating-point accumulation drift.
     */
    const int row = i / commondata->scan_density; // Vertical pixel coordinate index.
    const int col = i % commondata->scan_density; // Horizontal pixel coordinate index.

    const double x_pix = -commondata->window_width/2.0 + (col + 0.5) * (commondata->window_width / commondata->scan_density); // Local physical distance $x$.
    const double y_pix = -commondata->window_height/2.0 + (row + 0.5) * (commondata->window_height / commondata->scan_density); // Local physical distance $y$.

    // Global Cartesian intersection point on the projection window $x^\mu$.
    const double target_pos[3] = {{
        wc_x + x_pix*nx_0 + y_pix*ny_0,
        wc_y + x_pix*nx_1 + y_pix*ny_1,
        wc_z + x_pix*nx_2 + y_pix*ny_2
    }};

    // --- INITIAL STATE POPULATION ---
    /*
     * Algorithmic Step: Populate the thread-local state vector $f^\mu_{{local}}$ with starting coordinates.
     * Architectural Justification: Immediate sequential assignment into registers avoids temporary VRAM reads.
     */
    f_local[0] = commondata->t_start;
    f_local[1] = cam_x;
    f_local[2] = cam_y;
    f_local[3] = cam_z;

    const double V_x = target_pos[0] - cam_x; // Unnormalized geometric trajectory vector $V^x$.
    const double V_y = target_pos[1] - cam_y; // Unnormalized geometric trajectory vector $V^y$.
    const double V_z = target_pos[2] - cam_z; // Unnormalized geometric trajectory vector $V^z$.

    const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z); // Normalization scalar.

    f_local[5] = V_x * inv_mag_V;
    f_local[6] = V_y * inv_mag_V;
    f_local[7] = V_z * inv_mag_V;
    f_local[8] = 0.0; // Initial Affine parameter $\lambda$.

    // --- INTERPOLATION & METRIC EVALUATION ---
    /*
     * Algorithmic Step: Evaluate the metric tensor at the photon's starting coordinate.
     * Architectural Justification: Utilizing the inline physics helper forces the metric calculations
     * directly into the calling kernel's thread-local registers, preventing external function call overhead.
     */
    double pos_local[4] = {{f_local[0], f_local[1], f_local[2], f_local[3]}};
    placeholder_interpolation_engine_{spacetime_name}(commondata, pos_local, metric_local, NULL);

    // --- HAMILTONIAN CONSTRAINT ---
    /*
     * Algorithmic Step: Solve the negative root of the Hamiltonian constraint $p_\mu p^\mu = 0$ for initial $p^0$.
     * Architectural Justification: The inline device call executes strictly within the SM core.
     */
    p0_reverse(metric_local, f_local, &f_local[4]);

    // --- GLOBAL MEMORY WRITE ---
    /*
     * Algorithmic Step: Write the completed local state vector to the global SoA $f^\mu$.
     * Architectural Justification: Executing a single coalesced 9-strided write to global memory at the conclusion
     * of the kernel minimizes VRAM transactions and prevents warp serialization.
     */
    for(int m=0; m<9; m++) {{
        all_photons->f[IDX_GLOBAL(m, i, num_rays)] = f_local[m];
    }}
    """

    # Python: Hardware threading configuration for the RTX 3080.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(current_chunk_size + 256 - 1) / 256", "1", "1"],
    }

    # Python: Generate the standalone CUDA kernel utilizing the BHaH infrastructure.
    kernel_prefunc, launch_code = generate_kernel_and_launch_code(
        kernel_name=f"set_initial_conditions_kernel_{spacetime_name}",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict_cuda,
        arg_dict_host=arg_dict_cuda,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
    )

    # Python: Extract physics helpers to force them into a single translation unit.
    metric_func = cfc.CFunction_dict[f"g4DD_metric_{spacetime_name}"].full_function
    conn_func = cfc.CFunction_dict[f"connections_{spacetime_name}"].full_function
    interp_func = cfc.CFunction_dict[f"placeholder_interpolation_engine_{spacetime_name}"].full_function
    p0_func = cfc.CFunction_dict["p0_reverse"].full_function

    # Python: Define the ordered variables for CFunction registration.
    prefunc = f"{metric_func}\n{conn_func}\n{interp_func}\n{p0_func}\n{kernel_prefunc}"

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]

    desc = rf"""@brief Initializes Cartesian starting conditions for photons in {spacetime_name}.
    
    Detailed algorithm: Computes the geometric basis vectors of the camera window
    on the host, then launches a standalone device kernel in controlled batches to map 
    individual photon rays to pixel coordinates and populate the initial flattened Structure of Arrays
    (SoA) state vector. Enforces local register boundaries to prevent VRAM overflow.

    @param commondata The master configuration struct containing global spacetime parameters.
    @param num_rays Total photons in the simulation batch.
    @param all_photons Pointer to the master SoA state vector mapped to VRAM.
    @param window_center_out Output buffer for geometric center mapped to CPU.
    @param n_x_out Output buffer for the $x$-axis geometric basis vector mapped to CPU.
    @param n_y_out Output buffer for the $y$-axis geometric basis vector mapped to CPU.
    @param n_z_out Output buffer for the $z$-axis geometric basis vector mapped to CPU."""

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

    include_CodeParameters_h = False

    body = rf"""
    // --- HOST-SIDE GEOMETRY SETUP ---
    /*
     * Algorithmic Step: Construct orthonormal basis vectors for the virtual observer camera.
     * Architectural Justification: Evaluating static camera geometry exclusively within the host CPU 
     * context prevents redundant calculations inside the massively parallel VRAM kernel.
     */
    const double cam_x = commondata->camera_pos_x; // Spatial coordinate $x$ of observer's camera $x^\mu_{{cam}}$.
    const double cam_y = commondata->camera_pos_y; // Spatial coordinate $y$ of observer's camera $x^\mu_{{cam}}$.
    const double cam_z = commondata->camera_pos_z; // Spatial coordinate $z$ of observer's camera $x^\mu_{{cam}}$.

    const double wc_x = commondata->window_center_x; // Geometric center coordinate $x$ of projection window.
    const double wc_y = commondata->window_center_y; // Geometric center coordinate $y$ of projection window.
    const double wc_z = commondata->window_center_z; // Geometric center coordinate $z$ of projection window.

    double n_z[3] = {{wc_x - cam_x, wc_y - cam_y, wc_z - cam_z}}; // Vector normal to the camera window.
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]); // Normalization factor for $z$-axis vector.
    for(int j=0; j<3; j++) n_z[j] /= mag_n_z;

    const double guide_up[3] = {{commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}}; // Geometric reference vector defining "up".
    
    double n_x[3] = {{n_z[1]*guide_up[2] - n_z[2]*guide_up[1], n_z[2]*guide_up[0] - n_z[0]*guide_up[2], n_z[0]*guide_up[1] - n_z[1]*guide_up[0]}}; // Orthonormal basis vector $n_x^i$.
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]); // Normalization factor for $x$-axis vector.

    // Fallback mechanism to prevent cross-product singularities at the geometric poles.
    if (mag_n_x < 1e-9) {{
        double alternative_up[3] = {{0.0, 1.0, 0.0}}; // Alternative "up" direction mapped to the $y$-axis.
        if (fabs(n_z[1]) > 0.999) {{ alternative_up[1] = 0.0; alternative_up[2] = 1.0; }}
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    }}

    for(int j=0; j<3; j++) n_x[j] /= mag_n_x;

    double n_y[3] = {{n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]}}; // Orthonormal basis vector $n_y^i$.

    window_center_out[0] = wc_x;
    window_center_out[1] = wc_y;
    window_center_out[2] = wc_z;
    for(int j=0; j<3; j++) {{
        n_x_out[j] = n_x[j]; n_y_out[j] = n_y[j]; n_z_out[j] = n_z[j];
    }}

    const double nx_0 = n_x[0]; // Local basis variable mapped for $x$-axis vector $0$.
    const double nx_1 = n_x[1]; // Local basis variable mapped for $x$-axis vector $1$.
    const double nx_2 = n_x[2]; // Local basis variable mapped for $x$-axis vector $2$.
    const double ny_0 = n_y[0]; // Local basis variable mapped for $y$-axis vector $0$.
    const double ny_1 = n_y[1]; // Local basis variable mapped for $y$-axis vector $1$.
    const double ny_2 = n_y[2]; // Local basis variable mapped for $y$-axis vector $2$.

    // --- BATCH PROCESSING & KERNEL LAUNCH ---
    /*
     * Algorithmic Step: Launch the CUDA kernel to populate the initial state vector directly in VRAM.
     * Architectural Justification: The host-side loop orchestrates execution within strict BUNDLE_CAPACITY limits
     * to prevent PCIe bottlenecks and avoid exceeding the 10GB VRAM ceiling of the target hardware.
     */
    for (long int start_idx = 0; start_idx < num_rays; start_idx += BUNDLE_CAPACITY) {{
        // Computes the strict size of the active processing block.
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
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
    )
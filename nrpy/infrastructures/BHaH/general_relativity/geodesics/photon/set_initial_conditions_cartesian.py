"""
This module generates the native CUDA kernel and host-side orchestrator for initializing photon trajectories.

The high-level architecture utilizes a chunked VRAM staging buffer to populate the initial state vector
$f^\mu$ directly within the device context. This design ensures that all initial geometric mappings and 
Hamiltonian constraint solutions occur within the hardware registers of the sm_86 architecture, 
protecting the 10GB VRAM limit by processing rays in controlled bundles. Data is bridged back to 
pinned host memory via a 9-strided asynchronous memory copy to maintain the Streaming Bundle pipeline.

Author: Dalton J. Moone.
"""

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def set_initial_conditions_cartesian(spacetime_name: str) -> None:
    """
    Registers the C function and device kernel for Cartesian photon initialization.

    :param spacetime_name: The specific metric or spacetime identifier (e.g., 'KerrSchild').
    :raises Exception: Propagates NRPy+ core exceptions during AST generation.
    """
    # Python: Register necessary global parameters for the grid setup.
    par.register_CodeParameter("int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "t_start", 100, commondata=True, add_to_parfile=True)

    # Python: Dictionary mapping for GPU kernel arguments using the VRAM staging pointer.
    arg_dict_cuda = {
        "commondata": "const commondata_struct *restrict",
        "num_rays": "const long int",
        "d_f_chunk": "double *restrict",
        "cam_x": "const double", "cam_y": "const double", "cam_z": "const double",
        "wc_x": "const double", "wc_y": "const double", "wc_z": "const double",
        "nx_0": "const double", "nx_1": "const double", "nx_2": "const double",
        "ny_0": "const double", "ny_1": "const double", "ny_2": "const double",
        "start_idx": "const long int",
        "current_chunk_size": "const long int"
    }

    # Python: Define the GPU kernel body utilizing raw strings for LaTeX and strict variable documentation.
    kernel_body = rf"""
    // --- THREAD IDENTIFICATION & BOUNDARY CHECKS ---
    // The identifier $c$ represents the local thread index within the current bundle batch.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; 

    // Architectural Justification: Truncation prevents out-of-bounds VRAM access for non-uniform batches.
    if (c >= current_chunk_size) return;

    // The identifier $i$ represents the global ray index within the master $num\_rays$ SoA.
    const long int i = start_idx + c; 

    // --- THREAD-LOCAL TENSOR ALLOCATION ---
    // Thread-local array for the 9-element state vector $f^\mu$ and 4-momentum $p^\mu$.
    double f_local[9]; 
    // Thread-local array for the 10 independent metric components $g_{{\mu\nu}}$.
    double metric_local[10]; 

    // --- PIXEL MAPPING & GEOMETRY ---
    // Vertical pixel coordinate index within the virtual observer's projection frame.
    const int row = i / commondata->scan_density; 
    // Horizontal pixel coordinate index within the virtual observer's projection frame.
    const int col = i % commondata->scan_density; 

    // The variable $x_{{pix}}$ is the local physical distance along the horizontal camera axis.
    const double x_pix = -commondata->window_width/2.0 + (col + 0.5) * (commondata->window_width / commondata->scan_density); 
    // The variable $y_{{pix}}$ is the local physical distance along the vertical camera axis.
    const double y_pix = -commondata->window_height/2.0 + (row + 0.5) * (commondata->window_height / commondata->scan_density); 

    // The array $target\_pos$ stores the global Cartesian intersection point on the projection window $x^\mu$.
    const double target_pos[3] = {{
        wc_x + x_pix*nx_0 + y_pix*ny_0,
        wc_y + x_pix*nx_1 + y_pix*ny_1,
        wc_z + x_pix*nx_2 + y_pix*ny_2
    }};

    // --- INITIAL STATE POPULATION ---
    f_local[0] = commondata->t_start; // Coordinate time $t$
    f_local[1] = cam_x; // Spatial position $x$
    f_local[2] = cam_y; // Spatial position $y$
    f_local[3] = cam_z; // Spatial position $z$

    // Vector component $V^x$ for the unnormalized geometric trajectory.
    const double V_x = target_pos[0] - cam_x; 
    // Vector component $V^y$ for the unnormalized geometric trajectory.
    const double V_y = target_pos[1] - cam_y; 
    // Vector component $V^z$ for the unnormalized geometric trajectory.
    const double V_z = target_pos[2] - cam_z; 

    // Functional Justification: Normalization ensures initial 4-momentum $p^\mu$ satisfies null trajectory constraints.
    const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z); 

    f_local[5] = V_x * inv_mag_V; // Initial momentum component $p^x$
    f_local[6] = V_y * inv_mag_V; // Initial momentum component $p^y$
    f_local[7] = V_z * inv_mag_V; // Initial momentum component $p^z$
    f_local[8] = 0.0; // Initial affine parameter $\lambda$

    // --- INTERPOLATION & METRIC EVALUATION ---
    // Thread-local temporary array for the 4-position $x^\mu$.
    double pos_local[4] = {{f_local[0], f_local[1], f_local[2], f_local[3]}};
    placeholder_interpolation_engine_{spacetime_name}(commondata, pos_local, metric_local, NULL);

    // --- HAMILTONIAN CONSTRAINT ---
    // Algorithmic Step: Solve the quadratic Hamiltonian constraint $p_\mu p^\mu = 0$ specifically for $p^0$.
    p0_reverse(metric_local, f_local, &f_local[4]);

    // --- GLOBAL VRAM CHUNK WRITE ---
    // Algorithmic Step: Populate the device staging buffer $d\_f\_chunk$ with the validated state vector.
    // Hardware Justification: Single coalesced writes prevent warp serialization on the sm_86 memory controller.
    for(int m=0; m<9; m++) {{
        d_f_chunk[IDX_LOCAL(m, c, BUNDLE_CAPACITY)] = f_local[m];
    }}
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(current_chunk_size + 256 - 1) / 256", "1", "1"],
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

    # Python: Extract bodies from the registration dictionary to satisfy the Translation Unit Inlining Mandate.
    metric_func = cfc.CFunction_dict[f"g4DD_metric_{spacetime_name}"].full_function
    conn_func = cfc.CFunction_dict[f"connections_{spacetime_name}"].full_function
    interp_func = cfc.CFunction_dict[f"placeholder_interpolation_engine_{spacetime_name}"].full_function
    p0_func = cfc.CFunction_dict["p0_reverse"].full_function

    prefunc = f"{metric_func}\n{conn_func}\n{interp_func}\n{p0_func}\n{kernel_prefunc}"

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>", "<stdio.h>", "<stdlib.h>"]

    desc = rf"""@brief Initializes Cartesian starting conditions for photons in {spacetime_name}.
    @param commondata Master configuration struct containing global parameters.
    @param num_rays Total number of photon trajectories in the simulation, $num\_rays$.
    @param all_photons Pointer to the master SoA state vector $f^\mu$ in host memory.
    @param window_center_out Output buffer for the geometric window center.
    @param n_x_out Output buffer for the $x$-axis basis vector $n_x^i$.
    @param n_y_out Output buffer for the $y$-axis basis vector $n_y^i$.
    @param n_z_out Output buffer for the $z$-axis basis vector $n_z^i$."""

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

    body = rf"""
    // --- HOST-SIDE GEOMETRY SETUP ---
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

    // Fallback logic prevents cross-product singularities at geometric poles for numerical stability.
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
    // Device pointer for the chunked VRAM staging buffer $d\_f\_chunk$.
    double *d_f_chunk;
    // Hardware Justification: Memory is processed in bundles to protect the 10GB hardware limit.
    BHAH_MALLOC_DEVICE(d_f_chunk, sizeof(double) * 9 * BUNDLE_CAPACITY);

    for (long int start_idx = 0; start_idx < num_rays; start_idx += BUNDLE_CAPACITY) {{
        // Variable current_chunk_size defines the active range for the current streaming bundle.
        const long int current_chunk_size = (num_rays - start_idx < BUNDLE_CAPACITY) ? (num_rays - start_idx) : BUNDLE_CAPACITY;
        
        {launch_code}

        // --- EXPLICIT HARDWARE ERROR SYNCHRONIZATION ---
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
        for(int m=0; m<9; m++) {{
            cudaMemcpy(all_photons->f + (m * num_rays) + start_idx, 
                       d_f_chunk + (m * BUNDLE_CAPACITY), 
                       sizeof(double) * current_chunk_size, 
                       cudaMemcpyDeviceToHost);
        }}
    }}
    BHAH_FREE_DEVICE(d_f_chunk);
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
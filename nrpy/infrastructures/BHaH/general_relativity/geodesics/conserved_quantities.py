"""
Evaluates the physical conserved quantities along photon trajectories using a streaming bundle architecture.

This module implements a streaming bundle architecture to compute physical 
conserved quantities along photon trajectories. It strictly manages memory usage 
by processing photons in limited batches, bridging the Host-Device gap, ensuring 
compliance with hardware execution limits for both CUDA and OpenMP.

Author: Dalton J. Moone.
"""

from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.parallelization.utilities import (
    generate_kernel_and_launch_code,
    get_commondata_access,
    get_params_access,
)
from nrpy.helpers.loop import loop
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    Geodesic_Diagnostics,
)
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

def conserved_quantities(spacetime_name: str, particle_type: str = "photon") -> None:
    """
    Evaluates the C function and native kernel for conserved quantity evaluation.

    :param spacetime_name: The target analytic or numerical spacetime.
    :param particle_type: The particle classification (default is photon).
    :raises ValueError: If the particle type is not supported.
    """
    if particle_type == "massive":
        array_size = 8
    elif particle_type == "photon":
        array_size = 9
    else:
        raise ValueError(f"Unsupported particle_type: {particle_type}")

    # 1. Architecture Detection
    parallelization = par.parval_from_str("parallelization")
    cd_access = get_commondata_access(parallelization)
    params_access = get_params_access(parallelization)

    config_key = f"{spacetime_name}_{particle_type}"
    diagnostics = Geodesic_Diagnostics[config_key]

    # Set up SymPy expressions and target C variables for the struct assignment.
    list_of_syms: List[sp.Expr] = []
    list_of_c_vars: List[str] = []

    if diagnostics.E_expr is not None:
        list_of_syms.append(diagnostics.E_expr)
        list_of_c_vars.append("d_cq_bundle[c].E")

    if diagnostics.L_exprs:
        list_of_syms.extend(diagnostics.L_exprs)
        list_of_c_vars.extend(["d_cq_bundle[c].Lx", "d_cq_bundle[c].Ly", "d_cq_bundle[c].Lz"])

    if diagnostics.Q_expr is not None:
        list_of_syms.append(diagnostics.Q_expr)
        list_of_c_vars.append("d_cq_bundle[c].Q")

    # Define the C-structure for tracking physical conserved quantities.
    cq_struct_def = """
    // --- CONSERVED QUANTITIES STRUCTURE ---
    // Defines the physical conserved quantities evaluated along a photon trajectory.
    // Hardware Justification: This Structure of Arrays (AoS) definition is injected into the global BHaH_defines.h header to ensure uniform memory mapping across the Host orchestrator and VRAM kernels.
    typedef struct {
    double E;   // Energy $E$ extracted from the temporal Killing vector.
    double Lx;  // Angular momentum projection $L_x$.
    double Ly;  // Angular momentum projection $L_y$.
    double Lz;  // Angular momentum projection $L_z$ extracted from the azimuthal Killing vector.
    double Q;   // Carter constant $Q$ separating the Hamilton-Jacobi equations.
    } conserved_quantities_t;
    """

    # Register the struct definition to the global header generation pipeline.
    Bdefines_h.register_BHaH_defines("conserved_quantities", cq_struct_def)

    # Generate the highly optimized math evaluation block.
    math_kernel = ccg.c_codegen(
        list_of_syms,
        list_of_c_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    kernel_name = f"calculate_conserved_quantities_{spacetime_name}_{particle_type}_kernel"

    arg_dict_cuda = {
        "d_f_bundle": "const double *restrict",
        "d_cq_bundle": "conserved_quantities_t *restrict",
        "current_chunk_size": "const long int",
    }
    
    # Pass commondata explicitly when not using CUDA's global memory
    if parallelization != "cuda":
        arg_dict_cuda["commondata"] = "const commondata_struct *restrict"

  
    arg_dict_host = arg_dict_cuda.copy()

    # Dynamically generate the unpacking logic based on the specific spacetime coordinates.
    preamble_lines = [
        "    // --- COORDINATE HYDRATION ---",
        "    // Hardware Justification: Unpack position and momentum using strict SoA macros directly from the bundle."
    ]

    # Map the analytic spatial coordinates (e.g., $x, y, z$ or $r, \theta, \phi$).
    for i, symbol in enumerate(diagnostics.xx):
        var_name = str(symbol)
        preamble_lines.append(f"    const double {var_name} = d_f_bundle[IDX_LOCAL({i}, c, BUNDLE_CAPACITY)]; // Spatial coordinate ${var_name}$ evaluated from memory.")
        preamble_lines.append(f"    (void){var_name}; // Suppress unused variable warning for ${var_name}$.")

    # Map the temporal and spatial momenta ($p_t, p_x, p_y, p_z$).
    for i in range(4):
        preamble_lines.append(f"    const double p{i} = d_f_bundle[IDX_LOCAL({i+4}, c, BUNDLE_CAPACITY)]; // Momentum component $p_{i}$ evaluated from memory.")
        preamble_lines.append(f"    (void)p{i}; // Suppress unused variable warning for $p_{i}$.")

    preamble_lines.append("")
    preamble_lines.append("    // --- CONSTANT MEMORY HYDRATION ---")
    preamble_lines.append(f"    // Hardware Justification: Extract spacetime physics constants natively from the {cd_access} cache.")
    
    # Algorithmic Step: Scrape the SymPy abstract syntax tree for free variables and map them to the memory cache.
    free_syms = set()
    for expr in list_of_syms:
        free_syms.update(expr.free_symbols)
        
    for sym in free_syms:
        sym_name = str(sym)
        # Verify the variable is a registered global parameter (e.g., $a_{spin}$, $M_{scale}$).
        if sym_name in par.glb_code_params_dict:
            preamble_lines.append(f"    const double {sym_name} = {cd_access}{sym_name}; // Parameter ${sym_name}$ extracted from constant cache.")

    preamble = "\n".join(preamble_lines)
    
    # 2. The Kernel "Sandwich"
    if parallelization == "cuda":
        loop_preamble = """
    // --- CUDA THREAD IDENTIFICATION ---
    // Local 1D thread mapping within the current VRAM bundle.
    // Thread ID maps to a unique photon index $c$.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x; // Global thread evaluation index $c$.

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid bundle addresses.
    if (c >= current_chunk_size) return;
"""
        loop_postamble = ""
    else:
        loop_preamble = """
    // --- OPENMP LOOP ARCHITECTURE ---
    // Distribute photon rays across available CPU threads for parallel evaluation.
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

    // --- DIAGNOSTIC EVALUATION ---
    // Evaluates the analytic SymPy expressions for the conserved quantities.
{math_kernel}

    // --- MACRO CLEANUP ---
    #undef IDX_LOCAL
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
    }

    # 3. Update the Kernel Launch
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

    # --- CONDITIONAL HOST ORCHESTRATION ---
    if parallelization == "cuda":
        host_to_device_transfer = r"""
    // --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---
    // Transfer the state vector $f^\mu$ for the current bundle to VRAM.
    for(int m=0; m<9; m++) { // Loop over all 9 elements of the state vector $f^\mu$.
        cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
                   all_photons->f + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice);
    }
    """
        device_to_host_transfer = r"""
    // --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---
    // Retrieve the calculated conserved quantities back to the Pinned memory array.
    cudaMemcpy(cq_result + start_idx, d_cq_bundle,
               sizeof(conserved_quantities_t) * current_chunk_size, cudaMemcpyDeviceToHost);
    """
    else:
        host_to_device_transfer = r"""
    // --- SYNCHRONOUS MEMORY TRANSFER (HOST TO HOST CHUNK) ---
    // Populating the bundle array for localized CPU chunk processing.
    for(int m=0; m<9; m++) { // Loop over all 9 elements of the state vector $f^\mu$.
        memcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
               all_photons->f + (m * num_rays) + start_idx,
               sizeof(double) * current_chunk_size);
    }
    """
        device_to_host_transfer = r"""
    // --- SYNCHRONOUS MEMORY TRANSFER (HOST CHUNK TO HOST) ---
    // Retrieve the calculated conserved quantities back to the main memory array.
    memcpy(cq_result + start_idx, d_cq_bundle,
           sizeof(conserved_quantities_t) * current_chunk_size);
    """

    # Host-side chunked execution loop.
    loop_body = f"""
    // --- BUNDLE SIZING ---
    // Variable $current_chunk_size$ defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY); // Safely bounded active chunk parameter $current_chunk_size$.
    {host_to_device_transfer}
        // --- KERNEL LAUNCH ---
    {launch_body}
    {device_to_host_transfer}
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
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")
    
    desc = r"""@brief Computes conserved quantities for a batch of trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors $f^\mu$.
    @param num_rays The total number of photon trajectories.
    @param cq_result The array of diagnostic structures to be populated.

    Detailed algorithm:
    1. Allocates execution staging buffers for state vectors $f^\mu$ and quantity results.
    2. Iterates over the global dataset in chunks of $BUNDLE\_CAPACITY$.
    3. Prepares bundle memory, computes expressions natively on the targeted hardware, and unloads results."""

    cfunc_type = "void"
    name = f"calculate_conserved_quantities_universal_{spacetime_name}_{particle_type}"
    params = (
        "const commondata_struct *restrict commondata, "
        "const PhotonStateSoA *restrict all_photons, "
        "const long int num_rays, "
        "conserved_quantities_t *restrict cq_result"
    )
    include_CodeParameters_h = False

    body = fr"""
    // --- MEMORY STAGING ALLOCATION ---
    // Hardware Justification: Allocate buffers statically sized to $BUNDLE\_CAPACITY$. Maps to VRAM or RAM based on platform.
    double *d_f_bundle; // Pointer for the bundled state vector $f^\mu$.
    conserved_quantities_t *d_cq_bundle; // Pointer for the bundled diagnostic outputs.

    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY); // Allocates device/host target array.
    BHAH_MALLOC_DEVICE(d_cq_bundle, sizeof(conserved_quantities_t) * BUNDLE_CAPACITY); // Allocates device/host return array.

    // --- PAGINATION LOOP ---
    {host_loop}

    // --- MEMORY CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle); // Releases memory for state vector $f^\mu$.
    BHAH_FREE_DEVICE(d_cq_bundle); // Releases memory for diagnostic outputs.
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
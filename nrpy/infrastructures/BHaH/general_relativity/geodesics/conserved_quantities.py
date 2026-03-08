"""
Generates the CUDA kernel and launch code to evaluate conserved quantities.

This module implements a streaming bundle architecture to compute physical 
conserved quantities along photon trajectories. It strictly manages VRAM usage 
by processing photons in limited batches, using pinned memory transfers to bridge 
the Host-Device gap, ensuring compliance with the hardware memory limits.

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
)
from nrpy.helpers.loop import loop
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    Geodesic_Diagnostics,
)
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h

def conserved_quantities(spacetime_name: str, particle_type: str = "photon") -> None:
    """
    Generate the C function and native kernel for conserved quantity evaluation.

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
    # Hardware Justification: This Structure of Arrays (AoS) definition is injected into the global 
    # BHaH_defines.h header to ensure uniform memory mapping across the Host orchestrator and VRAM kernels.
    cq_struct_def = """
    // Defines the physical conserved quantities evaluated along a photon trajectory.
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


   # Fetch the architecture-specific constant memory prefix (e.g., "d_commondata.").
    cd_access = get_commondata_access("cuda")

    # Dynamically generate the unpacking logic based on the specific spacetime coordinates.
    preamble_lines = [
        "    // --- COORDINATE HYDRATION ---",
        "    // Hardware Justification: Unpack position and momentum using strict SoA macros directly from the VRAM bundle."
    ]

    # Map the analytic spatial coordinates (e.g., $x, y, z$ or $r, \theta, \phi$).
    for i, symbol in enumerate(diagnostics.xx):
        var_name = str(symbol)
        preamble_lines.append(f"    const double {var_name} = d_f_bundle[IDX_LOCAL({i}, c, BUNDLE_CAPACITY)];")
        preamble_lines.append(f"    (void){var_name};")

    # Map the temporal and spatial momenta ($p_t, p_x, p_y, p_z$).
    for i in range(4):
        preamble_lines.append(f"    const double p{i} = d_f_bundle[IDX_LOCAL({i+4}, c, BUNDLE_CAPACITY)];")
        preamble_lines.append(f"    (void)p{i};")

    preamble_lines.append("")
    preamble_lines.append("    // --- CONSTANT MEMORY HYDRATION ---")
    preamble_lines.append("    // Hardware Justification: Extract spacetime physics constants natively from the 64KB d_commondata cache.")
    
    # Algorithmic Step: Scrape the SymPy abstract syntax tree for free variables and map them to the GPU constant memory.
    free_syms = set()
    for expr in list_of_syms:
        free_syms.update(expr.free_symbols)
        
    for sym in free_syms:
        sym_name = str(sym)
        # Verify the variable is a registered global parameter (e.g., $a_{spin}$, $M_{scale}$).
        if sym_name in par.glb_code_params_dict:
            preamble_lines.append(f"    const double {sym_name} = {cd_access}{sym_name};")

    preamble = "\n".join(preamble_lines)
    

# Define the core __global__ kernel body executed on the device.
    kernel_body = f"""
    // --- MACRO DEFINITIONS ---
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    // Layout: [Component][RayID]
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(comp, ray_id, N) ((comp) * (N) + (ray_id))
    #endif

    // --- THREAD IDENTIFICATION ---
    // Local 1D thread mapping within the current VRAM bundle.
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;

    // --- BOUNDARY CHECK ---
    // Ensure out-of-bounds threads do not access invalid bundle addresses.
    if (c >= current_chunk_size) return;

    {preamble}

    // --- DIAGNOSTIC EVALUATION ---
    // Evaluates the analytic SymPy expressions for the conserved quantities.
    {math_kernel}
    """



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

    # Generate the kernel definition and the internal launch string.
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

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    
    desc = f"""@brief Computes conserved quantities for a batch of trajectories.
    @param all_photons The master Structure of Arrays containing the state vectors.
    @param num_rays The total number of photon trajectories.
    @param cq_result The array of diagnostic structures to be populated.

    Detailed algorithm:
    1. Allocates VRAM staging buffers for state vectors and quantity results.
    2. Iterates over the global dataset in chunks of BUNDLE_CAPACITY.
    3. Transfers state data Host->Device, computes expressions, and transfers Device->Host."""

    cfunc_type = "void"
    name = f"calculate_conserved_quantities_universal_{spacetime_name}_{particle_type}"
    params = (
        "const PhotonStateSoA *restrict all_photons, "
        "const long int num_rays, "
        "conserved_quantities_t *restrict cq_result"
    )

    # Host-side chunked execution loop.
    loop_body = f"""
    // Variable current_chunk_size defines the active range for the current streaming bundle.
    const long int current_chunk_size = MIN(num_rays - start_idx, BUNDLE_CAPACITY);

    // --- ASYNC MEMORY TRANSFER (HOST TO DEVICE) ---
    // Transfer the 9-component state vector $f^\mu$ for the current bundle to VRAM.
    // Hardware Justification: Component-wise transfers bounded by BUNDLE_CAPACITY prevent PCIe saturation.
    for(int m=0; m<9; m++) {{
        cudaMemcpy(d_f_bundle + (m * BUNDLE_CAPACITY),
                   all_photons->f + (m * num_rays) + start_idx,
                   sizeof(double) * current_chunk_size, cudaMemcpyHostToDevice);
    }}

    // --- KERNEL LAUNCH ---
    {launch_body}

    // --- ASYNC MEMORY TRANSFER (DEVICE TO HOST) ---
    // Retrieve the calculated conserved quantities back to the Pinned memory array.
    cudaMemcpy(cq_result + start_idx, d_cq_bundle,
               sizeof(conserved_quantities_t) * current_chunk_size, cudaMemcpyDeviceToHost);
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
    // Device pointer for the bundled state vector $f^\mu$.
    double *d_f_bundle;
    // Device pointer for the bundled diagnostic outputs.
    conserved_quantities_t *d_cq_bundle;

    // Hardware Justification: Allocate buffers statically sized to BUNDLE_CAPACITY to fit within 10GB VRAM limits.
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_cq_bundle, sizeof(conserved_quantities_t) * BUNDLE_CAPACITY);

    // --- HOST-SIDE PAGINATION LOOP ---
    {host_loop}

    // --- VRAM CLEANUP ---
    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_cq_bundle);
    """

    # Enforce strict canonical top-to-bottom sequence for registration without bloating postfunc.
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
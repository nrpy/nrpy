"""
Generates the CUDA kernel and host-side orchestrator for computing the temporal momentum component.

This module evaluates the quadratic Hamiltonian constraint to enforce physical null
trajectories for a batch of photons. It operates explicitly on VRAM bundles to resolve
the negative root of the constraint equation for $p^t$.

Author: Dalton J. Moone.
"""
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.helpers.parallelization.utilities import generate_kernel_and_launch_code

def p0_reverse_kernel(p0_expr: sp.Expr) -> None:
    """
    Register the global CUDA kernel for the initial temporal momentum calculation.

    The generated kernel loads the photon's spatial momenta and metric components 
    from global VRAM, executes the algebraically generated Hamiltonian constraint solver, 
    and writes the resulting temporal momentum $p^0$ back to the state vector bundle.

    :param p0_expr: The SymPy expression representing the negative root of the Hamiltonian constraint.
    :raises ValueError: If the symbolic expression fails to parse during code generation.
    """
    # Python: Input validation.
    if not p0_expr:
        raise ValueError("p0_expr must contain a valid symbolic expression.")

    # Python: Define the argument dictionary for the CUDA kernel generation.
    arg_dict = {
        "d_f_bundle": "double *restrict",
        "d_metric_bundle": "const double *restrict",
        "chunk_size": "const int"
    }

    # Python: Generate the raw C math string from the SymPy expression.
    # The output variable is named p0_val to act as a local register target.
    body_math = ccg.c_codegen(
        [p0_expr], ["p0_val"], enable_cse=True, enable_simd=True, verbose=False, include_braces=False
    )
    
    # Python: Translate SIMD macro signatures to native CUDA hardware intrinsics.
    body_math = body_math.replace("SIMD", "CUDA")

    # Python: Map the 2D symmetric metric components to explicit local register loads.
    metric_loads = []
    k = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            metric_loads.append(f"const double {comp_name} = ReadCUDA(&d_metric_bundle[IDX_METRIC({k}, i)]);")
            k += 1
    metric_load_str = "\n    ".join(metric_loads)

    # Python: Define the GPU kernel body.
    kernel_body = f"""
    // --- THREAD IDENTIFICATION ---
    // The identifier i represents the global thread index mapped to a specific photon ray.
    const int i = blockIdx.x * blockDim.x + threadIdx.x;

    // Hardware Justification: Guard prevents out-of-bounds VRAM access for threads exceeding the active chunk.
    if (i >= chunk_size) return;

    // --- MACRO DEFINITIONS FOR BUNDLE ACCESS ---
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))


    // --- SPATIAL MOMENTUM UNPACKING ---
    // Load contravariant spatial momentum components $p^x$, $p^y$, $p^z$ from VRAM.
    const double pU1 = ReadCUDA(&d_f_bundle[IDX_F(5, i)]);
    const double pU2 = ReadCUDA(&d_f_bundle[IDX_F(6, i)]);
    const double pU3 = ReadCUDA(&d_f_bundle[IDX_F(7, i)]);

    // --- METRIC TENSOR UNPACKING ---
    // Load the 10 independent metric components $g_{{\mu\\nu}}$ from VRAM into explicitly named registers.
    // Hardware Justification: These register names match the symbols expected by the generated C code.
    {metric_load_str}

    // --- HAMILTONIAN CONSTRAINT ROOT FINDING ---
    // Evaluate the algebraically generated solution for the temporal momentum $p^0$.
    // Local register storing the final evaluated negative root.
    double p0_val = 0.0;
    {body_math}

    // --- GLOBAL VRAM WRITE ---
    // Write the resulting temporal momentum back to component index 4 of the state bundle.
    WriteCUDA(&d_f_bundle[IDX_F(4, i)], p0_val);

    // --- MACRO CLEANUP ---
    #undef IDX_F
    #undef IDX_METRIC
    """

    # Python: Generate the kernel and the C host wrapper.
    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
    }

    prefunc, body = generate_kernel_and_launch_code(
        kernel_name="p0_reverse_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization="cuda",
        launch_dict=launch_dict,
        cfunc_decorators="__global__",
        thread_tiling_macro_suffix="RKF45"
    )

    # Python: Define arguments for C-Function registration strictly before the call.
    includes = ["BHaH_defines.h", "cuda_intrinsics.h"]
    
    desc = r"""@brief Orchestrates the CUDA kernel for the initial temporal momentum calculation.
    
    @param d_f_bundle Pointer to the state vector bundle $f^{\mu}$ in VRAM.
    @param d_metric_bundle Pointer to the pre-calculated metric bundle $g_{\mu\nu}$ in VRAM.
    @param chunk_size The number of active rays in the current bundle batch.
    """

    cfunc_type = "void"
    name = "p0_reverse_kernel"
    
    params = (
        "double *restrict d_f_bundle, "
        "const double *restrict d_metric_bundle, "
        "const int chunk_size"
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
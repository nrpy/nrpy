"""
Generates the C function for computing the photon geodesic ordinary differential equation right-hand sides.

This module writes the highly optimized C code required to evaluate the spatial and
temporal derivatives during the Runge-Kutta-Fehlberg (RKF45) integration step. The
generated C code maps local, intermediate tensor math (such as the metric $g_{\mu\nu}$
and Christoffel symbols $\Gamma^{\alpha}_{\beta\gamma}$) directly to hardware registers.
This explicit localization avoids excessive register spilling and global memory
accesses, maximizing throughput for High-Performance Computing (HPC) environments
executing highly parallel ray-tracing simulations on the NVIDIA sm_86 architecture.

Author: Dalton J. Moone.
"""

from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def calculate_ode_rhs(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    """
    Generate and register the C function to compute the ODE right-hand side.

    The generated C code unpacks spacetime coordinates and physical momenta from
    a thread-local state vector, mapping them to explicit hardware instructions. 
    It utilizes thread-local metric and connection arrays to formulate the 
    geodesic integration step without global memory indexing, staying within 
    the strict 255-register limit of the RTX 3080.

    :param geodesic_rhs_expressions: The mathematical right-hand side evaluations 
                                     representing the geodesic equations.
    :param coordinate_symbols: The spatial and temporal coordinate variables in order.
    :raises ValueError: If the provided geodesic expression list is empty.
    """
    if not geodesic_rhs_expressions:
        raise ValueError("geodesic_rhs_expressions must contain at least one mathematical expression.")

    # Include headers for hardware-agnostic definitions and optimized CUDA intrinsics.
    includes = ["BHaH_defines.h", "cuda_intrinsics.h"]

    # 1. Define function description with Doxygen and LaTeX
    desc = """@brief Computes the right-hand side of the photon geodesic ODE system.
    @param f_local Thread-local state vector containing coordinates $x^{\\mu}$ and momenta $p^{\\mu}$.
    @param metric_local Thread-local array containing the 10 upper-triangular components of $g_{\\mu\\nu}$.
    @param Gamma_local Thread-local array containing the 40 Christoffel symbols $\\Gamma^{\\alpha}_{\\mu\\nu}$.
    @param k_local Thread-local output array storing the evaluated derivatives $dk / d\\lambda$ for the current stage."""

    # 2. Define function type and signature
    # BHAH_HD_INLINE maps to __device__ __inline__ on GPU, forcing logic into registers.
    cfunc_type = "BHAH_HD_INLINE void"
    name = "calculate_ode_rhs"
    params = """const double *restrict f_local,
                  const double *restrict metric_local,
                  const double *restrict Gamma_local,
                  double *restrict k_local"""

    include_CodeParameters_h = False

    # 3. Generate the Body
    used_symbol_names = {
        str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols
    }

    preamble_lines = [
        "// --- STATE VECTOR & COORDINATE UNPACKING ---",
        "/* Unpack spacetime coordinates $x^{\\mu}$ from the local state vector. */",
        "/* Local register mapping prevents VRAM latency and avoids dead-code elimination. */"
    ]

    # Map the first 4 components as spacetime coordinates (t, x, y, z)
    for i, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(
                f"const double {str(sym)} = f_local[{i}]; // Spacetime coordinate representing ${str(sym)}$."
            )

    preamble_lines.extend([
        "\n  // --- MOMENTUM UNPACKING ---",
        "/* Unpack contravariant four-momenta $p^{\\mu}$ from the state vector indices 4-7. */"
    ])
    for i in range(4):
        if f"pU{i}" in used_symbol_names:
            preamble_lines.append(
                f"const double pU{i} = f_local[{i+4}]; // Contravariant momentum component $p^{{{i}}}$."
            )

    preamble_lines.extend([
        "\n  // --- METRIC TENSOR MAPPING ---",
        "/* Unpack the symmetric covariant metric $g_{\\mu\\nu}$ from thread-local memory. */",
        "/* Stored in registers to minimize the footprint of the fused integration kernel. */"
    ])
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"const double {comp_name} = metric_local[{curr_idx}]; // Covariant metric component $g_{{{m}{n}}}$."
                )
            curr_idx += 1

    preamble_lines.extend([
        "\n  // --- CHRISTOFFEL CONNECTION MAPPING ---",
        "/* Unpack Christoffel symbols $\\Gamma^{\\alpha}_{\\mu\\nu}$ from thread-local memory. */",
        "/* This architectural step occurs here to ensure coefficients are available for immediate RHS evaluation. */"
    ])
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(
                        f"const double {comp_name} = Gamma_local[{curr_idx}]; // Christoffel symbol $\\Gamma^{{{a}}}_{{{m}{n}}}$."
                    )
                curr_idx += 1

    preamble_lines.extend([
        "\n  // --- GEODESIC RHS EVALUATION ---",
        "/* Compute the derivatives $dx^{\\mu}/d\\lambda$ and $dp^{\\mu}/d\\lambda$ using hardware intrinsics. */",
        "/* enable_intrinsics=True ensures the use of FusedMulAddCUDA for peak throughput. */"
    ])

    body = "\n  ".join(preamble_lines) + "\n\n"

    # Map the output directly to the local derivative array indices (0-8)
    k_array_outputs = [f"k_local[{i}]" for i in range(9)]

    # Python: Generate the C code using the SIMD backend for vectorized mathematical operations
    raw_c_code = ccg.c_codegen(
        geodesic_rhs_expressions,
        k_array_outputs,
        enable_cse=True,
        enable_simd=True,
        include_braces=False,
        verbose=False,
    )

    # Python: Translate SIMD macro signatures to native CUDA hardware intrinsics
    body += raw_c_code.replace("SIMD", "CUDA")

    # 5. Register the function adhering to the canonical Master Order
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body
    )
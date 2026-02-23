"""
Register C function for computing photon geodesic ODE right-hand sides.
Optimized with a Preamble to minimize register pressure for GPU architectures.
"""
from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

def calculate_ode_rhs(
    geodesic_rhs_expressions: List[sp.Expr], coordinate_symbols: List[sp.Symbol]
) -> None:
    includes = ["BHaH_defines.h"]
    desc = """@brief Computes the 9 derivatives required for the geodesic ODE system."""
    name = "calculate_ode_rhs"

    c_params_str = """
                  const double *restrict f_temp,
                  const double *restrict metric_g4DD,
                  const double *restrict conn_GammaUDD,
                  double *restrict k_array,
                  const int bundle_capacity,
                  const int stage,
                  const long int j"""

    # --- Step 1: Identify used symbols ---
    # This prevents the preamble from being bloated with unused memory reads.
    used_symbol_names = {str(sym) for expr in geodesic_rhs_expressions for sym in expr.free_symbols}

    preamble_lines = [
        "#ifndef IDX_LOCAL",
        "#define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))",
        "#endif\n",
        "// Unpack position and momentum from f_temp"
    ]

    # --- Step 2: Build the Preamble ---
    
    # 2.a. Unpack State Vector (Coordinates and Momenta)
    for i, sym in enumerate(coordinate_symbols):
        if str(sym) in used_symbol_names:
            preamble_lines.append(f"const double {str(sym)} = f_temp[IDX_LOCAL({i}, j, bundle_capacity)];")

    for i in range(4):
        if f"pU{i}" in used_symbol_names:
            preamble_lines.append(f"const double pU{i} = f_temp[IDX_LOCAL({i+4}, j, bundle_capacity)];")

    # 2.b. Unpack Metric Components
    preamble_lines.append("\n  // Unpack metric components")
    curr_idx = 0
    for m in range(4):
        for n in range(m, 4):
            comp_name = f"metric_g4DD{m}{n}"
            if comp_name in used_symbol_names:
                preamble_lines.append(f"const double {comp_name} = metric_g4DD[IDX_LOCAL({curr_idx}, j, bundle_capacity)];")
            curr_idx += 1

    # 2.c. Unpack Christoffel Symbols
    preamble_lines.append("\n  // Unpack Christoffel symbols")
    curr_idx = 0
    for a in range(4):
        for m in range(4):
            for n in range(m, 4):
                comp_name = f"conn_Gamma4UDD{a}{m}{n}"
                if comp_name in used_symbol_names:
                    preamble_lines.append(f"const double {comp_name} = conn_GammaUDD[IDX_LOCAL({curr_idx}, j, bundle_capacity)];")
                curr_idx += 1

    # Combine preamble into a single string block
    body = "\n  ".join(preamble_lines) + "\n\n"

    # --- Step 3: Define Output Mapping ---
    k_array_outputs = [
        f"k_array[IDX_LOCAL((stage - 1) * 9 + {i}, j, bundle_capacity)]" for i in range(9)
    ]

    # --- Step 4: Generate the Kernel ---
    # Note: 'read_symbol_replacement_dict' is removed to avoid the TypeError.
    # The generator will now look for the local variables defined in the preamble.
    body += ccg.c_codegen(
        geodesic_rhs_expressions,
        k_array_outputs,
        enable_cse=True,
        include_braces=False,
        verbose=False,
    )

    # --- Step 5: Register ---
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=c_params_str,
        include_CodeParameters_h=False,
        body=body,
    )
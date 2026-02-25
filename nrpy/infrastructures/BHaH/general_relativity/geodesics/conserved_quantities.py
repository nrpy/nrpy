"""
Register C function for computing conserved quantities along geodesics.

Author: Dalton J. Moone
"""

from typing import List
import sympy as sp
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    Geodesic_Diagnostics,
)

def conserved_quantities(spacetime_name: str, particle_type: str = "massive") -> None:
    if particle_type == "massive":
        array_size = 8
        vec_desc = "4-velocity u^mu"
    elif particle_type == "photon":
        array_size = 9
        vec_desc = "4-momentum p^mu"
    else:
        raise ValueError(f"Unsupported particle_type: {particle_type}")

    config_key = f"{spacetime_name}_{particle_type}"
    diagnostics = Geodesic_Diagnostics[config_key]

    list_of_syms: List[sp.Expr] = []
    list_of_c_vars: List[str] = []

    c_params = (
        "const commondata_struct *restrict commondata, "
        "const double *restrict all_photons_f, "
        "long int num_rays, long int photon_idx"
    )

    if diagnostics.E_expr is not None:
        list_of_syms.append(diagnostics.E_expr)
        list_of_c_vars.append("*E")
        c_params += ", double *restrict E"

    if diagnostics.L_exprs:
        list_of_syms.extend(diagnostics.L_exprs)
        list_of_c_vars.extend(["*Lx", "*Ly", "*Lz"])
        c_params += ", double *restrict Lx, double *restrict Ly, double *restrict Lz"

    if diagnostics.Q_expr is not None:
        list_of_syms.append(diagnostics.Q_expr)
        list_of_c_vars.append("*Q")
        c_params += ", double *restrict Q"

    preamble_lines = ["// Unpack position and momentum using strict SoA macros"]
    
    # Added (void) casts to each unpacked variable to suppress unused variable warnings
    for i, symbol in enumerate(diagnostics.xx):
        var_name = str(symbol)
        preamble_lines.append(f"const double {var_name} = all_photons_f[IDX_GLOBAL({i}, photon_idx, num_rays)];")
        preamble_lines.append(f"(void){var_name};")
        
    for i in range(4):
        preamble_lines.append(f"const double p{i} = all_photons_f[IDX_GLOBAL({i+4}, photon_idx, num_rays)];")
        preamble_lines.append(f"(void)p{i};")

    preamble = "\n    ".join(preamble_lines)

    includes = ["BHaH_defines.h"]
    desc = f"""@brief Computes conserved quantities for {spacetime_name} ({particle_type}).
    Expects flattened SoA master buffer and local batch metric buffer."""
    name = f"conserved_quantities_{spacetime_name}_{particle_type}"
    params = c_params

    kernel = ccg.c_codegen(
        list_of_syms,
        list_of_c_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

   # Wrap for GPU-side diagnostic checks
    body = f"""
        #ifdef USE_GPU
        #pragma omp declare target
        #endif
        {preamble}

        {kernel}
        #ifdef USE_GPU
        #pragma omp end declare target
        #endif
        """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    print(f"    ... {name}() registration complete.")

if __name__ == "__main__":
    import logging
    import os
    import sys
    sys.path.append(os.getcwd())
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestConservedGen")
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "massive"
    logger.info("Test: Generating Conserved Quantities C-code...")
    try:
        conserved_quantities(SPACETIME, PARTICLE)
        cfunc_name = f"conserved_quantities_{SPACETIME}_{PARTICLE}"
        if cfunc_name not in cfc.CFunction_dict:
            raise RuntimeError(f"FAIL: '{cfunc_name}' was not registered.")
        cfunc = cfc.CFunction_dict[cfunc_name]
        filename = f"{cfunc_name}.c"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(cfunc.full_function)
    except Exception as e:
        logger.error(" -> FAIL: Test failed with error: %s", e)
        sys.exit(1)
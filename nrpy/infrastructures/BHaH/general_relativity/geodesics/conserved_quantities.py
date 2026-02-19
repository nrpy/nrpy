"""
Register C function for computing conserved quantities along geodesics.

This module registers the 'conserved_quantities_{spacetime}_{particle_type}' C function.
It uses a C-code preamble to unpack the input state vector f[N_state] into local variables
matching the symbolic names (e.g., xx0, p0) used in the physics modules. For massive
particles N_state = 8 (position + 4-velocity, assuming m=1); for photon particles
N_state = 9 (position + 4-momentum + path-length diagnostic).

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
from typing import List

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.equations.general_relativity.geodesics.geodesic_diagnostics.conserved_quantities import (
    Geodesic_Diagnostics,
)


def conserved_quantities(spacetime_name: str, particle_type: str = "massive") -> None:
    """
    Generate and register the C function to compute conserved quantities.

    This function instantiates the GeodesicDiagnostics class, identifies available
    conserved quantities (E, L, Q), and generates a C function.
    Assumes mass of all massive particles is m=1.

    Particle Support:
    - Massive: f[8], v^mu = u^mu (4-velocity).
    - Photon: f[9], v^mu = p^mu (4-momentum).

    :param spacetime_name: Name of the spacetime (e.g., "KerrSchild_Cartesian").
    :param particle_type: "massive" or "photon".
    :raises ValueError: If PARTICLE is not "massive" or "photon".
    """
    # Step 1: Specific setup based on particle type
    if particle_type == "massive":
        array_size = 8
        vec_desc = "4-velocity u^mu"
    elif particle_type == "photon":
        array_size = 9
        vec_desc = "4-momentum p^mu"
    else:
        raise ValueError(f"Unsupported particle_type: {particle_type}")

    # Step 2: Acquire Symbolic Diagnostics
    config_key = f"{spacetime_name}_{particle_type}"
    diagnostics = Geodesic_Diagnostics[config_key]

    # Step 3: Determine Active Quantities and Build Lists
    list_of_syms: List[sp.Expr] = []
    list_of_c_vars: List[str] = []

    # Function signature: standard generic interface
    # f[0..3]: Coordinates (t, x, y, z or similar)
    # f[4..7]: 4-Momentum (p0, p1, p2, p3)
    # f[8]   : Proper Length L (Photons only)
    c_params = (
        f"const commondata_struct *restrict commondata, const double f[{array_size}]"
    )

    # 3.a: Energy
    if diagnostics.E_expr is not None:
        list_of_syms.append(diagnostics.E_expr)
        list_of_c_vars.append("*E")
        c_params += ", double *restrict E"

    # 3.b: Angular Momentum (Vector)
    if diagnostics.L_exprs:
        list_of_syms.extend(diagnostics.L_exprs)
        list_of_c_vars.extend(["*Lx", "*Ly", "*Lz"])
        c_params += ", double *restrict Lx, double *restrict Ly, double *restrict Lz"

    # 3.c: Carter Constant
    if diagnostics.Q_expr is not None:
        list_of_syms.append(diagnostics.Q_expr)
        list_of_c_vars.append("*Q")
        c_params += ", double *restrict Q"

    # Step 4: Generate the Dynamic Preamble
    # We unpack f[] into local 'const double' variables that match the SymPy symbols.

    preamble_lines = ["// Unpack position coordinates from f[0]..f[3]"]

    # 4.a: Unpack Coordinates (xx0, xx1, etc.)
    # The diagnostics object holds the specific symbols used (e.g., xx0 vs t).
    for i, symbol in enumerate(diagnostics.xx):
        preamble_lines.append(f"const double {str(symbol)} = f[{i}];")

    preamble_lines.append("")
    preamble_lines.append("// Unpack 4-vector components from f[4]..f[7]")
    preamble_lines.append("// For massive particles this is the 4-velocity u^mu (m=1),")
    preamble_lines.append("// for photon particles this is the 4-momentum p^mu.")

    # 4.b: Unpack components (p0, p1, p2, p3)
    # The symbolic diagnostics module uses the names p0..p3 for the 4-vector;
    # here p0..p3 represent u^mu for massive particles and p^mu for photons.
    for i in range(4):
        # SymPy symbol is "p0", maps to f[4]
        preamble_lines.append(f"const double p{i} = f[{i+4}];")

    preamble = "\n    ".join(preamble_lines)

    # Step 5: Define Metadata
    includes = ["BHaH_defines.h"]
    desc = f"""@brief Computes conserved quantities for {spacetime_name} ({particle_type}).
            
        Computed: {', '.join([v.strip('*') for v in list_of_c_vars])}.

        Expects state vector f[{array_size}] containing position and {vec_desc}.

        Note: Input vector components are mapped to local variables p0..p3 to match 
        the symbolic variable names used in the equation generation module."""
    name = f"conserved_quantities_{spacetime_name}_{particle_type}"
    params = c_params

    # Step 6: Generate C Body

    # Generate the computation kernel using the ORIGINAL symbols (no substitution).
    # The C compiler will resolve "p0" because we defined it in the preamble.
    kernel = ccg.c_codegen(
        list_of_syms,
        list_of_c_vars,
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Combine Preamble + Kernel
    body = preamble + "\n\n" + kernel

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

    # Ensure local modules can be imported
    sys.path.append(os.getcwd())

    # Configure logging
    logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger("TestConservedGen")

    # Test Configuration
    SPACETIME = "KerrSchild_Cartesian"
    PARTICLE = "massive"

    logger.info(
        "Test: Generating Conserved Quantities C-code for %s_%s...", SPACETIME, PARTICLE
    )

    try:
        # 1. Run the Generator
        conserved_quantities(SPACETIME, PARTICLE)

        # 2. Validation
        cfunc_name = f"conserved_quantities_{SPACETIME}_{PARTICLE}"

        # Check Registration
        if cfunc_name not in cfc.CFunction_dict:
            raise RuntimeError(f"FAIL: '{cfunc_name}' was not registered.")

        cfunc = cfc.CFunction_dict[cfunc_name]
        logger.info(" -> PASS: '%s' registered.", cfunc_name)

        # 3. Output to file
        filename = f"{cfunc_name}.c"
        with open(filename, "w", encoding="utf-8") as f:
            f.write(cfunc.full_function)
        logger.info("    ... Wrote %s", filename)

    except Exception as e:  # pylint: disable=broad-exception-caught
        logger.error(" -> FAIL: Test failed with error: %s", e)
        import traceback

        traceback.print_exc()
        sys.exit(1)

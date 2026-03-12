"""
Register the C function for computing the initial time component of 4-velocity.

This module registers the 'u0_massive' C function. It enforces
the 4-velocity normalization constraint for massive particles (u.u = -1) by solving
the quadratic Hamiltonian constraint for the time component u^0.

It utilizes the flattened SoA architecture via local and global batch indexing.

Author: Dalton J. Moone
"""

import sys

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def u0_massive(u0_expr: sp.Expr) -> None:
    """
    Generate and register the C function to compute u^0 for a massive particle.
    Utilizes the flattened SoA architecture via local and global batch indexing.

    :param u0_expr: The SymPy expression for u^0.
    """
    # Step 1: Define C function metadata
    includes = ["BHaH_defines.h"]
    desc = """@brief Computes the initial time-component of the 4-velocity (u^0).

    Solves the quadratic Hamiltonian constraint equation:
        g_munu u^mu u^nu = -1
    for the positive root of u^0 using the SoA architecture.

    Input:
        metric_g4DD: Flattened array of metric components.
        all_photons_f: Global state vector array.
        num_rays: Total number of rays for global indexing.
        batch_size: Size of the SIMD/GPU batch.
        particle_idx: Global index of the particle.
        batch_id: Local index within the current batch.
    Output:
        u0_out: The computed u^0 component, also synced to all_photons_f[4]."""

    name = "u0_massive"

    params = (
        "const double *restrict metric_g4DD, "
        "double *restrict all_photons_f, "
        "const long int num_rays, "
        "const int batch_size, "
        "const long int particle_idx, "
        "const int batch_id, "
        "double *restrict u0_out"
    )

    print(f" -> Generating C worker function: {name} ...")

    # Step 2: Map symmetric metric components to flat array indices
    metric_map = {}
    k = 0
    for i in range(4):
        for j in range(i, 4):
            metric_map[f"metric_g4DD{i}{j}"] = (
                f"metric_g4DD[IDX_LOCAL({k}, batch_id, batch_size)]"
            )
            k += 1

    # Step 3: Generate the Math Body (using CSE)
    body_math = ccg.c_codegen(
        [u0_expr], ["*u0_out"], enable_cse=True, verbose=False, include_braces=False
    )

    # Replace symbolic metric names with array indexing
    for old, new in metric_map.items():
        body_math = body_math.replace(old, new)

    # Step 4: Generate Preamble and Postamble
    # Unpack spatial velocity components from global state vector (indices 5, 6, 7)
    preamble = (
        "// Unpack spatial velocity components from global state vector (indices 5, 6, 7)\n"
        "  // Note: index 4 is u^0 (which we are computing)\n"
        "  const double uU1 = all_photons_f[IDX_GLOBAL(5, particle_idx, num_rays)];\n"
        "  const double uU2 = all_photons_f[IDX_GLOBAL(6, particle_idx, num_rays)];\n"
        "  const double uU3 = all_photons_f[IDX_GLOBAL(7, particle_idx, num_rays)];\n"
    )

    # Sync result back to the global state vector (index 4)
    postamble = "\n  // Sync result back to the global state vector (index 4)\n"
    postamble += "  all_photons_f[IDX_GLOBAL(4, particle_idx, num_rays)] = *u0_out;"

    body = f"{preamble}\n{body_math}\n{postamble}"

    # Step 5: Define GPU Offloading Pragmas
    prefunc = "#ifdef USE_GPU\n#pragma omp declare target\n#endif"
    postfunc = "#ifdef USE_GPU\n#pragma omp end declare target\n#endif"

    # Step 6: Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
        postfunc=postfunc,
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

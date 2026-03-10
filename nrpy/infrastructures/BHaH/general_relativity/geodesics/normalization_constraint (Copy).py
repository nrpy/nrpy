"""
Register C function for computing the normalization constraint of the state vector.

This module registers the 'normalization_constraint_{PARTICLE}' C function.
It computes the scalar invariant C = g_munu v^mu v^nu, where v^mu is the
4-vector (4-velocity or 4-momentum) contained in the global state vector f.

Particle Support:
- Massive: f[8], v^mu = u^mu (4-velocity). Expected C = -1.
- Photon: f[9], v^mu = p^mu (4-momentum). Expected C = 0.

It utilizes the flattened SoA architecture via local and global batch indexing.

Author: Dalton J. Moone
"""

import logging
import sys

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def normalization_constraint(norm_expr: sp.Expr, PARTICLE: str) -> None:
    """
    Generate and register the C function to compute the normalization constraint.

    :param norm_expr: The SymPy expression for the contraction g_munu v^mu v^nu.
    :param PARTICLE: The type of particle ("massive" or "photon").
                          Determines array size and naming convention.
    :raises ValueError: If PARTICLE is not "massive" or "photon".
    """
    # Step 1: Specific setup based on particle type
    if PARTICLE == "massive":
        vec_desc = "4-velocity u^mu"
        expected_val = "-1.0"
    elif PARTICLE == "photon":
        vec_desc = "4-momentum p^mu"
        expected_val = "0.0"
    else:
        raise ValueError(f"Unsupported PARTICLE: {PARTICLE}")

    # Step 2: Define C function metadata
    includes = ["BHaH_defines.h"]
    name = f"normalization_constraint_{PARTICLE}"

    desc = f"""@brief Computes the normalization constraint of the 4-vector.

        Evaluates the scalar invariant:
            C = g_munu v^mu v^nu
        where v^mu corresponds to the {vec_desc} stored in the global SoA array.

        Expected Value: {expected_val}"""

    params = (
        "const double *restrict metric_g4DD, "
        "const double *restrict all_photons_f, "
        "const long int num_rays, "
        "const int batch_size, "
        "const long int photon_idx, "
        "const int batch_id, "
        "double *restrict norm_out"
    )

    # Step 3: Generate C body
    print(f" -> Generating C worker function: {name}...")

    # 3a. Map metric symbolic names to 1D SoA local flat array
    metric_map = {}
    k = 0
    for i in range(4):
        for j in range(i, 4):
            old_symbol = f"metric_g4DD{i}{j}"
            new_access = f"metric_g4DD[IDX_LOCAL({k}, batch_id, batch_size)]"
            metric_map[old_symbol] = new_access
            k += 1

    # 3b. Generate the Math Body (using CSE)
    body_math = ccg.c_codegen(
        [norm_expr],
        ["*norm_out"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    # Apply the string replacements for the metric tensor components
    for old, new in metric_map.items():
        body_math = body_math.replace(old, new)

    # 3c. Generate the Dynamic Preamble for Global SoA extraction
    preamble = f"""
    // Unpack {vec_desc} components from global state vector
    const double vU0 = all_photons_f[IDX_GLOBAL(4, photon_idx, num_rays)];
    const double vU1 = all_photons_f[IDX_GLOBAL(5, photon_idx, num_rays)];
    const double vU2 = all_photons_f[IDX_GLOBAL(6, photon_idx, num_rays)];
    const double vU3 = all_photons_f[IDX_GLOBAL(7, photon_idx, num_rays)];
    """

    # Combine preamble and math
    full_body = f"{preamble}\n{body_math}"

    # Step 4: Define GPU OpenMP wrappers
    prefunc = """
#ifdef USE_GPU
#pragma omp declare target
#endif
"""

    postfunc = """
#ifdef USE_GPU
#pragma omp end declare target
#endif
"""

    # Step 5: Register the C function
    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=full_body,
        postfunc=postfunc,
    )
    print(f"    ... {name}() registration complete.")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")
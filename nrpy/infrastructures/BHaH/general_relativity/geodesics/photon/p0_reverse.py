"""
Module for generating the C code that computes the temporal momentum component (p^0).

This constraint solver forces the four-momentum vector to remain null
by solving the quadratic Hamiltonian constraint for the negative root.
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc


def p0_reverse(p0_expr: sp.Expr) -> None:
    """
    Generate and registers the C function to compute p^0 for a photon.

    Utilizes the flattened SoA architecture via local and global batch indexing.

    :param p0_expr: The SymPy expression representing the negative root of the Hamiltonian constraint for p^0.

    >>> import sympy as sp
    >>> p0_reverse(sp.sympify("0.0"))
    """
    includes = ["BHaH_defines.h"]

    desc = """@brief Computes the initial temporal component of the 4-momentum (p^0).
    Detailed algorithm: Solves the quadratic Hamiltonian constraint p_mu p^mu = 0
    specifically for the negative root to enforce physical null photon trajectories
    within the SoA architecture."""

    name = "p0_reverse"

    params = (
        "const double *restrict metric_g4DD, "
        "double *restrict all_photons_f, "
        "const long int num_rays, "
        "const int batch_size, "
        "const long int photon_idx, "
        "const int batch_id, "
        "double *restrict p0_out"
    )

    metric_map = {}
    k = 0
    for i in range(4):
        for j in range(i, 4):
            old_symbol = f"metric_g4DD{i}{j}"
            new_access = f"metric_g4DD[IDX_LOCAL({k}, batch_id, batch_size)]"
            metric_map[old_symbol] = new_access
            k += 1

    body_math = ccg.c_codegen(
        [p0_expr], ["*p0_out"], enable_cse=True, verbose=False, include_braces=False
    )

    for old, new in metric_map.items():
        body_math = body_math.replace(old, new)

    preamble = """
    // Contravariant spatial momentum component p^x extracted from global state vector.
    const double pU1 = all_photons_f[IDX_GLOBAL(5, photon_idx, num_rays)];
    // Contravariant spatial momentum component p^y extracted from global state vector.
    const double pU2 = all_photons_f[IDX_GLOBAL(6, photon_idx, num_rays)];
    // Contravariant spatial momentum component p^z extracted from global state vector.
    const double pU3 = all_photons_f[IDX_GLOBAL(7, photon_idx, num_rays)];
    """

    postamble = "\n    all_photons_f[IDX_GLOBAL(4, photon_idx, num_rays)] = *p0_out;"

    body = f"""
    {preamble}
    {body_math}
    {postamble}
    """

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

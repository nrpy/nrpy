# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/normalization_constraint_photon_normalized.py
r"""
Defines the normalized photon constraint evaluated along batched trajectories.

This module implements a CPU/OpenMP execution kernel to compute the scalar
diagnostic ``C = \gamma^{ij} \Pi_i \Pi_j`` for normalized photon trajectories.
It mirrors the batch-oriented photon infrastructure pattern used elsewhere in
this directory: the generated helper reads the normalized photon state from a
flattened Structure of Arrays (SoA) bundle, unpacks the local metric data
needed by the symbolic expression, and writes one scalar diagnostic structure
per active ray.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par


def normalization_constraint_photon_normalized(norm_expr: sp.Expr) -> None:
    r"""
    Register the CPU/OpenMP normalized photon constraint kernel.

    :param norm_expr: The SymPy expression for the scalar diagnostic
        ``\gamma^{ij} \Pi_i \Pi_j``.
    :raises ValueError: If ``norm_expr`` is not a SymPy expression.
    :raises ValueError: If the configured parallelization mode is not OpenMP.

    Doctests:
    >>> import os
    >>> import nrpy.c_function as cfc
    >>> import nrpy.equations.general_relativity.geodesics.geodesics as geo
    >>> os.environ["XDG_CACHE_HOME"] = "/tmp"
    >>> cfc.CFunction_dict.clear()
    >>> generic_geodesic_equations = geo.GeodesicEquations.__new__(geo.GeodesicEquations)
    >>> norm_expr = generic_geodesic_equations.normalization_constraint_photon_normalized()
    >>> normalization_constraint_photon_normalized(norm_expr)
    >>> generated = cfc.CFunction_dict["normalization_constraint_photon_normalized"].full_function
    >>> "#pragma omp parallel for" in generated
    True
    >>> "d_norm_bundle[c].C" in generated
    True
    """
    if not isinstance(norm_expr, sp.Expr):
        raise ValueError("norm_expr must be a valid SymPy expression.")

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "normalization_constraint_photon_normalized currently supports only "
            "parallelization='openmp'."
        )

    norm_struct_def = r"""
    //==========================================
    // NORMALIZATION CONSTRAINT STRUCTURE
    //==========================================
    // Defines the physical normalization constraint evaluated along a trajectory.
    typedef struct {
        double C;   // Scalar diagnostic $C = \gamma^{ij} \Pi_i \Pi_j$.
    } normalization_constraint_t; // END STRUCT: normalization_constraint_t
    """
    Bdefines_h.register_BHaH_defines("normalization_constraint", norm_struct_def)

    used_symbol_names = {str(sym) for sym in norm_expr.free_symbols}

    math_kernel = ccg.c_codegen(
        [norm_expr],
        ["d_norm_bundle[c].C"],
        enable_cse=True,
        verbose=False,
        include_braces=False,
    )

    kernel_name = "normalization_constraint_photon_normalized_kernel"
    arg_dict = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_norm_bundle": "normalization_constraint_t *restrict",
        "current_chunk_size": "const long int",
    }

    preamble_lines = [
        "    //==========================================",
        "    // COMPONENT HYDRATION",
        "    //==========================================",
    ]

    if "u" in used_symbol_names:
        preamble_lines.append(
            "    const double u = d_f_bundle[IDX_LOCAL(4, c, BUNDLE_CAPACITY)]; // Normalized photon energy variable $u$."
        )
        preamble_lines.append(
            "    (void)u; // Suppress unused variable warnings when symbolic simplification removes $u$."
        )

    for i in range(3):
        symbol_name = f"PiD{i}"
        if symbol_name in used_symbol_names:
            preamble_lines.append(
                f"    const double {symbol_name} = d_f_bundle[IDX_LOCAL({i + 5}, c, BUNDLE_CAPACITY)]; // Normalized covariant momentum component $\\Pi_{{{i}}}$."
            )
            preamble_lines.append(
                f"    (void){symbol_name}; // Suppress unused variable warnings for $\\Pi_{{{i}}}$."
            )

    metric_idx = 0
    for i in range(4):
        for j in range(i, 4):
            comp_name = f"metric_g4DD{i}{j}"
            if comp_name in used_symbol_names:
                preamble_lines.append(
                    f"    const double {comp_name} = d_metric_bundle[IDX_LOCAL({metric_idx}, c, BUNDLE_CAPACITY)]; // Metric tensor component $g_{{{i}{j}}}$."
                )
                preamble_lines.append(
                    f"    (void){comp_name}; // Suppress unused variable warnings for $g_{{{i}{j}}}$."
                )
            metric_idx += 1

    preamble = "\n".join(preamble_lines)

    loop_preamble = """
    //==========================================
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    // Distribute photon trajectories across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(long int c = 0; c < current_chunk_size; c++) {
    """
    loop_postamble = "    } // END LOOP: for c over current_chunk_size"

    core_math = f"""
    //==========================================
    // MACRO DEFINITIONS
    //==========================================
    // IDX_LOCAL maps a component to the flattened state bundle using SoA layout.
    #ifndef IDX_LOCAL
    #define IDX_LOCAL(comp, ray_id, N) ((comp) * (N) + (ray_id))
    #endif

{preamble}

    //==========================================
    // CONSTRAINT EVALUATION
    //==========================================
    // Evaluates the analytic SymPy expression for the normalized photon constraint.
    {math_kernel}
    """

    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    prefunc, launch_body = parallel_utils.generate_kernel_and_launch_code(
        kernel_name=kernel_name,
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=None,
        thread_tiling_macro_suffix="DEFAULT",
        cfunc_decorators="",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]

    desc = r""" Computes the normalized photon constraint for a batch of trajectories.
    @param d_f_bundle The pointer array containing the normalized photon state vectors.
    @param d_metric_bundle The pointer array containing the symmetric metric tensor $g_{\mu\nu}$.
    @param d_norm_bundle The array of diagnostic constraint structures to be populated.
    @param current_chunk_size The dynamically sized operational boundary for the active chunk.

    Expected Value: 1.0 for normalized photon trajectories."""

    cfunc_type = "void"
    name = "normalization_constraint_photon_normalized"
    params = (
        "const double *restrict d_f_bundle, "
        "const double *restrict d_metric_bundle, "
        "normalization_constraint_t *restrict d_norm_bundle, "
        "const long int current_chunk_size, "
        "const int stream_idx"
    )
    include_CodeParameters_h = False
    body = launch_body

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

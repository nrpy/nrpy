r"""
Defines the CPU/OpenMP helper that converts direct photon momentum to normalized variables.

This module registers a host-only C helper used during numerical photon
initialization after the temporal momentum component ``p^0`` has been recovered.
It processes a batch of photons in parallel, reads the direct four-momentum
state ``(p^0, p^1, p^2, p^3)`` together with the local four-metric ``g_{\mu\nu}``,
and overwrites the momentum slots with the normalized variables
``u = \ln|\alpha p^0|`` and ``\Pi_i = p_i / (\alpha p^0)`` needed by the
coordinate-time evolution system.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

from typing import List

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def photon_momentum_to_normalized_kernel(
    u_expr: sp.Expr, PiD_exprs: List[sp.Expr]
) -> None:
    r"""
    Register the CPU/OpenMP helper that converts photon momentum to normalized form.

    This helper is intended for the numerical photon initialization path after
    ``p0_reverse_kernel()`` has solved for the signed temporal momentum
    component. It operates on a batch of photons and overwrites the state bundle
    slots storing the direct momentum variables with the normalized quantities
    used by the Bohn-style coordinate-time formulation.

    :param u_expr: SymPy expression for ``u = \ln|\alpha p^0|``.
    :param PiD_exprs: Three SymPy expressions for the covariant normalized
        spatial momentum components ``\Pi_i``.
    :raises ValueError: If ``u_expr`` is empty.
    :raises ValueError: If ``PiD_exprs`` does not contain exactly three expressions.
    :raises ValueError: If the configured parallelization mode is not OpenMP.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.equations.general_relativity.geodesics.geodesics as geo
    >>> cfc.CFunction_dict.clear()
    >>> u_expr, PiD_exprs = geo.GeodesicEquations.photon_momentum_to_normalized_quantities()
    >>> photon_momentum_to_normalized_kernel(u_expr, PiD_exprs)
    >>> generated = cfc.CFunction_dict["photon_momentum_to_normalized_kernel"].full_function
    >>> "#pragma omp parallel for" in generated
    True
    >>> "d_f_bundle[IDX_F(4, i)] = u_out;" in generated
    True
    """
    if not u_expr:
        raise ValueError("u_expr must contain a valid symbolic expression.")
    if len(PiD_exprs) != 3:
        raise ValueError(
            "PiD_exprs must contain exactly three symbolic expressions for Pi_i."
        )

    parallelization = par.parval_from_str("parallelization")
    if parallelization != "openmp":
        raise ValueError(
            "photon_momentum_to_normalized_kernel currently supports only "
            "parallelization='openmp'."
        )

    # Keep the argument metadata in the standard kernel-generator layout, but
    # reject non-OpenMP builds above so this remains a CPU-only helper.
    arg_dict = {
        "d_f_bundle": "double *restrict",
        "d_metric_bundle": "const double *restrict",
        "chunk_size": "const int",
    }

    # Step 1: Generate the raw C math for the conversion
    #
    #   alpha = 1 / sqrt(-g^{00})
    #   p_i = g_{{i0}} p^0 + g_{{ij}} p^j
    #   u = ln|alpha p^0|
    #   Pi_i = p_i / (alpha p^0)
    #
    # from the symbolic expressions supplied by geodesics.py.
    output_names = ["u_out", "PiD0_out", "PiD1_out", "PiD2_out"]
    body_math = ccg.c_codegen(
        [u_expr] + PiD_exprs,
        output_names,
        enable_cse=True,
        enable_simd=False,
        verbose=False,
        include_braces=False,
    )

    # Step 2: Load the metric tensor components into explicit local registers.
    metric_loads = []
    metric_idx = 0
    for mu in range(4):
        for nu in range(mu, 4):
            metric_loads.append(
                f"// Covariant metric component $g_{{{mu}{nu}}}$.\n    const double metric_g4DD{mu}{nu} = d_metric_bundle[IDX_METRIC({metric_idx}, i)];"
            )
            metric_idx += 1
    metric_load_str = "\n    ".join(metric_loads)

    # Step 3: Build the kernel body for the CPU/OpenMP batch conversion.
    loop_preamble = """
    //==========================================
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    // Distribute photon rays across available CPU threads for parallel evaluation.
    #pragma omp parallel for
    for(int i = 0; i < chunk_size; i++) {
    """
    loop_postamble = "    } // END LOOP: for i over chunk_size rays"

    core_math = f"""
    //==========================================
    // MACRO DEFINITIONS FOR BUNDLE ACCESS
    //==========================================
    // IDX_F maps a component to the flattened state bundle using SoA layout.
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    // IDX_METRIC maps a component to the flattened symmetric metric bundle.
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    //==========================================
    // DIRECT PHOTON MOMENTUM UNPACKING
    //==========================================
    // Load the direct four-momentum components from the state bundle.
    const double pU0 = d_f_bundle[IDX_F(4, i)];
    const double pU1 = d_f_bundle[IDX_F(5, i)];
    const double pU2 = d_f_bundle[IDX_F(6, i)];
    const double pU3 = d_f_bundle[IDX_F(7, i)];

    //==========================================
    // METRIC TENSOR UNPACKING
    //==========================================
    // Load the 10 independent metric components $g_{{\\mu\\nu}}$ into explicitly named registers.
    {metric_load_str}

    //==========================================
    // NORMALIZED PHOTON CONVERSION
    //==========================================
    // Evaluate the one-time initialization conversion
    //
    //   alpha = 1 / sqrt(-g^{00})
    //   p_i = g_{{i0}} p^0 + g_{{ij}} p^j
    //   u = ln|alpha p^0|
    //   Pi_i = p_i / (alpha p^0)
    //
    // from the local four-metric and the direct four-momentum.
    double u_out = 0.0;
    double PiD0_out = 0.0;
    double PiD1_out = 0.0;
    double PiD2_out = 0.0;
    {body_math}

    //==========================================
    // GLOBAL MEMORY WRITE
    //==========================================
    // Overwrite only the legacy momentum slots with the normalized state
    // variables. The affine-parameter slot f[8] is left untouched here.
    d_f_bundle[IDX_F(4, i)] = u_out;
    d_f_bundle[IDX_F(5, i)] = PiD0_out;
    d_f_bundle[IDX_F(6, i)] = PiD1_out;
    d_f_bundle[IDX_F(7, i)] = PiD2_out;

    //==========================================
    // MACRO CLEANUP
    //==========================================
    #undef IDX_F
    #undef IDX_METRIC
    """
    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }

    prefunc, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="photon_momentum_to_normalized_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="",
        thread_tiling_macro_suffix="RKF45",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]

    desc = r""" Convert direct photon four-momentum to normalized photon variables.

    @param d_f_bundle Pointer to the photon state vector bundle in memory.
    @param d_metric_bundle Pointer to the pre-calculated metric bundle $g_{\mu\nu}$ in memory.
    @param chunk_size The number of active rays in the current bundle batch.
    @param stream_idx Unused placeholder kept only for interface compatibility with other photon helpers.
    """

    cfunc_type = "void"
    name = "photon_momentum_to_normalized_kernel"
    params = (
        "double *restrict d_f_bundle, "
        "const double *restrict d_metric_bundle, "
        "const int chunk_size, "
        "const int stream_idx"
    )

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=launch_code,
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

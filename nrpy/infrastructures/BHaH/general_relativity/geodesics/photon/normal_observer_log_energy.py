r"""
Emit the direct-EOM helper for the common normal-observer log-energy measure.

The helper reads direct ``p^0`` and the covariant four-metric ``g4DD`` and
computes ``ln|alpha p^0|``. The symbolic expression supplied by
``geodesics.py`` constructs the needed inverse-metric component internally.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par


def normal_observer_log_energy(log_energy_expr: sp.Expr) -> None:
    r"""
    Register the direct-EOM log-energy helper for CPU or CUDA execution.

    :param log_energy_expr: Symbolic expression for ``ln|alpha p^0|``.
    :raises ValueError: If the expression is empty or the parallelization is unsupported.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.equations.general_relativity.geodesics.geodesics as geo
    >>> cfc.CFunction_dict.clear()
    >>> u_expr, _ = geo.GeodesicEquations.photon_momentum_to_normalized_quantities()
    >>> normal_observer_log_energy(u_expr)
    >>> generated = cfc.CFunction_dict["normal_observer_log_energy"].full_function
    >>> "#pragma omp parallel for" in generated
    True
    >>> "WriteCUDA(&d_log_energy_bundle[i], log_energy_out);" in generated
    True
    >>> "metric_g4DD00" in generated
    True
    """
    if not log_energy_expr:
        raise ValueError("log_energy_expr must contain a valid symbolic expression.")

    parallelization = par.parval_from_str("parallelization")
    if parallelization not in ("openmp", "cuda"):
        raise ValueError(
            "normal_observer_log_energy supports only "
            "parallelization='openmp' or 'cuda'."
        )

    arg_dict = {
        "d_f_bundle": "const double *restrict",
        "d_metric_bundle": "const double *restrict",
        "d_log_energy_bundle": "double *restrict",
        "chunk_size": "const int",
    }

    body_math = ccg.c_codegen(
        [log_energy_expr],
        ["log_energy_out"],
        enable_cse=True,
        enable_simd=False,
        verbose=False,
        include_braces=False,
    )

    metric_loads = []
    metric_idx = 0
    for mu in range(4):
        for nu in range(mu, 4):
            metric_loads.append(
                f"// Covariant metric component $g_{{{mu}{nu}}}$.\n"
                f"    const double metric_g4DD{mu}{nu} = "
                f"ReadCUDA(&d_metric_bundle[IDX_METRIC({metric_idx}, i)]);"
            )
            metric_idx += 1
    metric_load_str = "\n    ".join(metric_loads)

    if parallelization == "cuda":
        loop_preamble = """
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= chunk_size) {
        return;
    } // END IF: thread is outside active ray chunk
        """
        loop_open = ""
        loop_close = ""
        loop_postamble = ""
    else:
        loop_preamble = """
    #pragma omp parallel for
    for (int i = 0; i < chunk_size; ++i) {
        """
        loop_open = ""
        loop_close = "    }"
        loop_postamble = ""

    kernel_body = f"""
    #define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
    #define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

    {loop_preamble}
    {loop_open}
        const double pU0 = ReadCUDA(&d_f_bundle[IDX_F(4, i)]);
        {metric_load_str}

        double log_energy_out = 0.0;
        {body_math}
        WriteCUDA(&d_log_energy_bundle[i], log_energy_out);
    {loop_close} // END LOOP: for i over chunk_size rays
    {loop_postamble}

    #undef IDX_F
    #undef IDX_METRIC
    """

    launch_dict = {
        "threads_per_block": ["256", "1", "1"],
        "blocks_per_grid": ["(chunk_size + 256 - 1) / 256", "1", "1"],
        "stream": "stream_idx",
    }
    prefunc, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="normal_observer_log_energy",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
        thread_tiling_macro_suffix="RKF45",
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h"]
    if parallelization == "cuda":
        includes.extend(["cuda_intrinsics.h", "BHaH_device_defines.h"])

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=r"""Compute the direct-EOM normal-observer log-energy measure.

        @param[in] d_f_bundle Direct photon state bundle containing $p^0$ in
        component $f^4$.
        @param[in] d_metric_bundle Covariant four-metric bundle $g_{\mu\nu}$.
        @param[out] d_log_energy_bundle Per-ray $\ln|\alpha p^0|$ values.
        @param chunk_size Number of active rays in the bundle.
        @param stream_idx CPU compatibility argument.
        """,
        cfunc_type="void",
        name="normal_observer_log_energy",
        params=(
            "const double *restrict d_f_bundle, "
            "const double *restrict d_metric_bundle, "
            "double *restrict d_log_energy_bundle, "
            "const int chunk_size, "
            "const int stream_idx"
        ),
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
    print(f"Doctest passed: All {results.attempted} test(s) passed")

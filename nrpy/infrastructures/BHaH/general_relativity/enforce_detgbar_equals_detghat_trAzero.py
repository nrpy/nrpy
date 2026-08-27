"""
Generate C function enforcing conformal metric determinant and A trace constraints.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_algebraic_constraints import (
    BSSN_algebraic_constraints,
)
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH


def register_CFunction_enforce_detgbar_equals_detghat_trAzero(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register combined determinant and trace-free conformal-tensor projection.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None in registration phase, otherwise the updated NRPy environment.

    Doctests:
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> par.set_parval_from_str("enable_parallel_codegen", False)
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> import contextlib
    >>> import io
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _ = register_CFunction_enforce_detgbar_equals_detghat_trAzero(
    ...         "Spherical", False, False, 1
    ...     )
    >>> print(
    ...     cfc.CFunction_dict[
    ...         "enforce_detgbar_equals_detghat_trAzero__rfm__Spherical"
    ...     ].full_function
    ... )
    #include "BHaH_defines.h"
    <BLANKLINE>
    /**
     * Kernel: enforce_detgbar_equals_detghat_trAzero_host.
     * Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0.
     */
    static void enforce_detgbar_equals_detghat_trAzero_host(const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1,
                                                            const REAL *restrict x2, REAL *restrict in_gfs, const REAL *restrict auxevol_gfs) {
      MAYBE_UNUSED const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
      MAYBE_UNUSED const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
      MAYBE_UNUSED const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;
    <BLANKLINE>
      MAYBE_UNUSED const REAL invdxx0 = params->invdxx0;
      MAYBE_UNUSED const REAL invdxx1 = params->invdxx1;
      MAYBE_UNUSED const REAL invdxx2 = params->invdxx2;
    <BLANKLINE>
    #pragma omp parallel for
      for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
        for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
          for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    <BLANKLINE>
            /*
             * NRPy-Generated GF Access/FD Code, Step 1 of 2:
             * Read gridfunction(s) from main memory and compute FD stencils as needed.
             */
            const REAL aDD00 = in_gfs[IDX4(ADD00GF, i0, i1, i2)];
            const REAL aDD01 = in_gfs[IDX4(ADD01GF, i0, i1, i2)];
            const REAL aDD02 = in_gfs[IDX4(ADD02GF, i0, i1, i2)];
            const REAL aDD11 = in_gfs[IDX4(ADD11GF, i0, i1, i2)];
            const REAL aDD12 = in_gfs[IDX4(ADD12GF, i0, i1, i2)];
            const REAL aDD22 = in_gfs[IDX4(ADD22GF, i0, i1, i2)];
            const REAL hDD00 = in_gfs[IDX4(HDD00GF, i0, i1, i2)];
            const REAL hDD01 = in_gfs[IDX4(HDD01GF, i0, i1, i2)];
            const REAL hDD02 = in_gfs[IDX4(HDD02GF, i0, i1, i2)];
            const REAL hDD11 = in_gfs[IDX4(HDD11GF, i0, i1, i2)];
            const REAL hDD12 = in_gfs[IDX4(HDD12GF, i0, i1, i2)];
            const REAL hDD22 = in_gfs[IDX4(HDD22GF, i0, i1, i2)];
    <BLANKLINE>
            /*
             * NRPy-Generated GF Access/FD Code, Step 2 of 2:
             * Evaluate SymPy expressions and write to main memory.
             */
            const REAL FDPart3tmp0 = hDD00 + 1;
            const REAL FDPart3tmp2 = hDD22 + 1;
            const REAL FDPart3tmp4 = hDD11 + 1;
            const REAL FDPart3tmp6 = FDPart3tmp0 * FDPart3tmp2 * FDPart3tmp4 - FDPart3tmp0 * ((hDD12) * (hDD12)) - FDPart3tmp2 * ((hDD01) * (hDD01)) -
                                     FDPart3tmp4 * ((hDD02) * (hDD02)) + 2 * hDD01 * hDD02 * hDD12;
            const REAL FDPart3tmp7 = (1.0 / cbrt(FDPart3tmp6));
            const REAL FDPart3tmp8 = (1.0 / (FDPart3tmp6));
            const REAL FDPart3tmp9 = 2 * FDPart3tmp8;
            const REAL FDPart3tmp10 =
                FDPart3tmp8 * aDD00 * (FDPart3tmp2 * FDPart3tmp4 - ((hDD12) * (hDD12))) +
                FDPart3tmp8 * aDD11 * (FDPart3tmp0 * FDPart3tmp2 - ((hDD02) * (hDD02))) +
                FDPart3tmp8 * aDD22 * (FDPart3tmp0 * FDPart3tmp4 - ((hDD01) * (hDD01))) + FDPart3tmp9 * aDD01 * (-FDPart3tmp2 * hDD01 + hDD02 * hDD12) +
                FDPart3tmp9 * aDD02 * (-FDPart3tmp4 * hDD02 + hDD01 * hDD12) + FDPart3tmp9 * aDD12 * (-FDPart3tmp0 * hDD12 + hDD01 * hDD02);
            const REAL FDPart3tmp11 = (1.0 / 3.0) * FDPart3tmp10;
            in_gfs[IDX4(HDD00GF, i0, i1, i2)] = FDPart3tmp0 * FDPart3tmp7 - 1;
            in_gfs[IDX4(HDD01GF, i0, i1, i2)] = FDPart3tmp7 * hDD01;
            in_gfs[IDX4(HDD02GF, i0, i1, i2)] = FDPart3tmp7 * hDD02;
            in_gfs[IDX4(HDD11GF, i0, i1, i2)] = FDPart3tmp4 * FDPart3tmp7 - 1;
            in_gfs[IDX4(HDD12GF, i0, i1, i2)] = FDPart3tmp7 * hDD12;
            in_gfs[IDX4(HDD22GF, i0, i1, i2)] = FDPart3tmp2 * FDPart3tmp7 - 1;
            in_gfs[IDX4(ADD00GF, i0, i1, i2)] = -FDPart3tmp10 * ((1.0 / 3.0) * hDD00 + 1.0 / 3.0) + aDD00;
            in_gfs[IDX4(ADD01GF, i0, i1, i2)] = -FDPart3tmp11 * hDD01 + aDD01;
            in_gfs[IDX4(ADD02GF, i0, i1, i2)] = -FDPart3tmp11 * hDD02 + aDD02;
            in_gfs[IDX4(ADD11GF, i0, i1, i2)] = -FDPart3tmp10 * ((1.0 / 3.0) * hDD11 + 1.0 / 3.0) + aDD11;
            in_gfs[IDX4(ADD12GF, i0, i1, i2)] = -FDPart3tmp11 * hDD12 + aDD12;
            in_gfs[IDX4(ADD22GF, i0, i1, i2)] = -FDPart3tmp10 * ((1.0 / 3.0) * hDD22 + 1.0 / 3.0) + aDD22;
    <BLANKLINE>
          } // END LOOP: for i0 over [0, Nxx_plus_2NGHOSTS0)
        } // END LOOP: for i1 over [0, Nxx_plus_2NGHOSTS1)
      } // END LOOP: for i2 over [0, Nxx_plus_2NGHOSTS2)
    } // END FUNCTION: enforce_detgbar_equals_detghat_trAzero_host
    <BLANKLINE>
    /**
     * Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0.
     */
    void enforce_detgbar_equals_detghat_trAzero__rfm__Spherical(const params_struct *restrict params, const REAL *restrict x0, const REAL *restrict x1,
                                                                const REAL *restrict x2, REAL *restrict in_gfs, const REAL *restrict auxevol_gfs) {
      enforce_detgbar_equals_detghat_trAzero_host(params, x0, x1, x2, in_gfs, auxevol_gfs);
    } // END FUNCTION: enforce_detgbar_equals_detghat_trAzero__rfm__Spherical
    <BLANKLINE>
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")

    includes = ["BHaH_defines.h"]
    desc = r"""Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0."""
    cfunc_type = "void"
    name = "enforce_detgbar_equals_detghat_trAzero"
    arg_dict_cuda = {
        "in_gfs": "REAL *restrict",
        "auxevol_gfs": "const REAL *restrict",
    }
    if enable_rfm_precompute:
        arg_dict_cuda = {
            "rfmstruct": "const rfm_struct *restrict",
            **arg_dict_cuda,
        }
    else:
        arg_dict_cuda = {
            "x0": "const REAL *restrict",
            "x1": "const REAL *restrict",
            "x2": "const REAL *restrict",
            **arg_dict_cuda,
        }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    params = ",".join([f"{v} {k}" for k, v in arg_dict_host.items()])

    hprimeDD, aprimeDD = BSSN_algebraic_constraints(CoordSystem, enable_rfm_precompute)
    access_gfs = [
        gri.BHaHGridFunction.access_gf(
            f"{basename}{i}{j}", 0, 0, 0, gf_array_name="in_gfs"
        )
        for basename in ("hDD", "aDD")
        for i in range(3)
        for j in range(i, 3)
    ]
    projected_exprs = [
        expressions[i][j]
        for expressions in (hprimeDD, aprimeDD)
        for i in range(3)
        for j in range(i, 3)
    ]

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    kernel_body = BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            projected_exprs,
            access_gfs,
            enable_fd_codegen=True,
            enable_simd=False,
            enable_fd_functions=enable_fd_functions,
        ),
        loop_region="all points",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=False,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )
    loop_params = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=False
    )
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(projected_exprs)
    params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=False,
        var_access=parallel_utils.get_params_access(parallelization),
    )
    kernel_body = f"{loop_params}\n{params_definitions}\n{kernel_body}"

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=f"static {cfunc_type}",
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="DETGTRAZERO",
    )

    prefunc = ""
    if parallelization == "cuda" and enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc(
            cfunc_decorators="__device__ "
        ).replace("SIMD", "CUDA")
    elif enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc()

    cfc.register_CFunction(
        include_CodeParameters_h=False,
        prefunc=prefunc + kernel,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=launch_body,
        enable_simd=False,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

"""
Generate C functions to store the first derivatives of hDD in curvilinear coordinates.

Ricci needs every mixed second derivative of hDD. NRPy builds a mixed stencil as the
tensor product of two centered first-derivative stencils, so storing partial_k h_ij once
per point lets Ricci_eval rebuild each mixed second derivative as a single first
derivative of the stored gridfunction, and read partial_k h_ij itself instead of
recomputing it. The gridfunction is named hDDdD after the reference-metric convention
ghatDDdD. It lives in the SCRATCH group: no array is allocated for it, and the caller supplies
storage that is dead between hDDdD_eval and Ricci_eval within one right-hand-side
evaluation, such as the Method of Lines buffer the right-hand sides are about to overwrite.

A stored derivative must be valid wherever it is differentiated again, so each
direction is produced over the interior grown by fd_order/2 points in every direction in
which Ricci_eval differentiates it: partial_0 h in i1 and i2, partial_1 h in i2, and
partial_2 h is read pointwise only. A direction that a symmetry axis has zeroed is
neither differentiated nor grown.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Dict, List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.infrastructures import BHaH


def register_hDDdD_gridfunctions() -> None:
    """
    Register the SCRATCH gridfunctions hDDdD if they are not already registered.

    Both the producer and its consumers call this, because parallel code generation runs
    each registration in its own worker.
    """
    if "hDDdD000" not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions_for_single_rankN(
            "hDDdD",
            rank=3,
            symmetry="sym01",
            group="SCRATCH",
        )


def hDDdD_substitutions() -> Dict[sp.Symbol, sp.Symbol]:
    """
    Map hDD derivative symbols onto reads of the stored hDDdD gridfunctions.

    Each first derivative becomes the stored gridfunction, and each mixed second
    derivative becomes a first derivative of the stored gridfunction, canonicalized to
    k < l so that h_{ij,kl} and h_{ij,lk} remain one shared value. Unmixed second
    derivatives are not tensor products of first-derivative stencils and are left alone.
    A derivative that a symmetry axis has already zeroed is not a symbol, so it is
    absent from the map and stays zero; hDDdD_eval stores nothing along such an axis.

    :return: Substitution dictionary for sympy xreplace.

    DocTests:
    >>> subs = hDDdD_substitutions()
    >>> subs[sp.Symbol("hDD_dD012")], subs[sp.Symbol("hDD_dDD0112")]
    (hDDdD012, hDDdD_dD0112)
    >>> len(subs)  # 18 first derivatives and 18 mixed second derivatives of six components
    36
    >>> sp.Symbol("hDD_dDD0111") in subs
    False
    """
    hDD_dD = ixp.declarerank3("hDD_dD", symmetry="sym01")
    hDD_dDD = ixp.declarerank4("hDD_dDD", symmetry="sym01_sym23")
    hDDdD = ixp.declarerank3("hDDdD", symmetry="sym01")
    hDDdD_dD = ixp.declarerank4("hDDdD_dD", symmetry="sym01")
    subs: Dict[sp.Symbol, sp.Symbol] = {}
    for i in range(3):
        for j in range(3):
            for k in range(3):
                if hDD_dD[i][j][k] != 0:
                    subs[hDD_dD[i][j][k]] = hDDdD[i][j][k]
                for l in range(3):
                    if k != l and hDD_dDD[i][j][k][l] != 0:
                        subs[hDD_dDD[i][j][k][l]] = hDDdD_dD[i][j][min(k, l)][max(k, l)]
    return subs


def register_CFunction_hDDdD_eval(
    CoordSystem: str,
    enable_intrinsics: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that stores the first derivatives of hDD.

    :param CoordSystem: The coordinate system to be used.
    :param enable_intrinsics: Whether to enable SIMD/CUDA intrinsics.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    is_cuda = par.parval_from_str("parallelization") == "cuda"
    # Registers the evolved BSSN gridfunctions this kernel differentiates. Parallel code
    # generation runs each registration in its own worker, so the dependency is stated
    # here, exactly as Ricci_eval states it.
    _ = BSSN_quantities[CoordSystem + "_rfm_precompute"]
    register_hDDdD_gridfunctions()
    hDD_dD = ixp.declarerank3("hDD_dD", symmetry="sym01")
    hDD_dDD = ixp.declarerank4("hDD_dDD", symmetry="sym01_sym23")
    # Each stored direction must be valid wherever Ricci_eval differentiates it: the region is the
    # interior grown by the stencil radius in every transverse direction that some non-zero mixed
    # second derivative differentiates it in. Directions zeroed by a symmetry axis need no halo.
    loop_regions = []
    for direction in range(3):
        grown = [
            f"i{transverse}"
            for transverse in range(direction + 1, 3)
            if any(
                hDD_dDD[i][j][direction][transverse] != 0
                for i in range(3)
                for j in range(i, 3)
            )
        ]
        loop_regions += [
            "interior plus stencil halo in " + " ".join(grown) if grown else "interior"
        ]

    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += [
            str(
                Path("intrinsics") / "cuda_intrinsics.h"
                if is_cuda
                else Path("intrinsics") / "simd_intrinsics.h"
            )
        ]
    desc = r"""Set the first derivatives of hDD."""
    cfunc_type = "void"
    name = "hDDdD_eval"
    arg_dict_cuda = {
        "in_gfs": "const REAL *restrict",
        "scratch_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    params = ",".join([f"{v} {k}" for k, v in arg_dict_host.items()])

    # c_codegen() clears the finite-difference helper registry at the start of every call,
    # so each direction's helpers are collected before the next call and emitted once.
    collected_FDFunctions = {}
    prefunc = ""
    launch_body = ""
    for direction, loop_region in enumerate(loop_regions):
        exprs = []
        access_gfs: List[str] = []
        for i in range(3):
            for j in range(i, 3):
                # Nothing to store along a symmetry axis; hDDdD_substitutions() reads
                # nothing there either.
                if hDD_dD[i][j][direction] == 0:
                    continue
                exprs += [hDD_dD[i][j][direction]]
                access_gfs += [
                    f"scratch_gfs[IDX4(HDDDD{i}{j}{direction}GF, i0, i1, i2)]"
                ]
        if not exprs:
            continue

        point_body = ccg.c_codegen(
            exprs,
            access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            enable_fd_functions=enable_fd_functions,
            rational_const_alias=("static constexpr" if is_cuda else "static const"),
        ).replace("SIMD", "CUDA" if is_cuda else "SIMD")
        collected_FDFunctions.update(fin.FDFunctions_dict)

        kernel_body = BHaH.simple_loop.simple_loop(
            loop_body=point_body,
            loop_region=loop_region,
            enable_intrinsics=enable_intrinsics,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=False,
            read_xxs=False,
            OMP_collapse=OMP_collapse,
        )
        loop_params = parallel_utils.get_loop_parameters(
            "cuda" if is_cuda else "openmp", enable_intrinsics=enable_intrinsics
        )
        param_symbols, _ = get_params_commondata_symbols_from_expr_list(exprs)
        params_definitions = generate_definition_header(
            param_symbols,
            enable_intrinsics=enable_intrinsics,
            var_access=parallel_utils.get_params_access(
                "cuda" if is_cuda else "openmp"
            ),
        )
        kernel_body = f"{loop_params}\n{params_definitions}\n{kernel_body}"

        kernel, direction_launch_body = parallel_utils.generate_kernel_and_launch_code(
            f"{name}_dirn{direction}",
            kernel_body.replace("SIMD", "CUDA" if is_cuda else "SIMD"),
            arg_dict_cuda,
            arg_dict_host,
            parallelization="cuda" if is_cuda else "openmp",
            comments=f"{desc} Direction {direction}.",
            cfunc_type=f"static {cfunc_type}",
            # One function launches every direction, so each launch block needs its own
            # scope for the thread-count locals it declares.
            launchblock_with_braces=True,
            launch_dict={
                **BHaH.parallelization.cuda_utilities.default_launch_dictionary,
                "threads_per_block": ["64", "1", "1"],
            },
            thread_tiling_macro_suffix="HDDDD_EVAL",
        )
        prefunc += kernel
        launch_body += direction_launch_body

    if enable_fd_functions:
        fin.FDFunctions_dict.clear()
        fin.FDFunctions_dict.update(collected_FDFunctions)
        prefunc = (
            fin.construct_FD_functions_prefunc(
                cfunc_decorators="__device__ " if is_cuda else ""
            ).replace("SIMD", "CUDA" if is_cuda else "SIMD")
            + prefunc
        )

    cfc.register_CFunction(
        include_CodeParameters_h=False,
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=launch_body,
        enable_simd=enable_intrinsics,
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

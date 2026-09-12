# nrpy/infrastructures/BHaH/general_relativity/cfdD_alphadD_vetUdD_eval.py
"""
Generate a C function that stores the first derivatives of cf, alpha and vetU.

The BSSN right-hand sides need every mixed second derivative of cf, alpha and the three
components of vetU: fifteen of them, each the tensor product of two centered
first-derivative stencils, which at eighth order is 64 nonzero terms costing 74
double-precision operations. Storing partial_0 f and partial_1 f once per point lets
rhs_eval rebuild each mixed second derivative as a single nine-point first derivative of
the stored gridfunction, and read partial_0 f and partial_1 f themselves instead of
recomputing them. The gridfunctions are named cfdD, alphadD and vetUdD after the
reference-metric convention ghatDDdD.

The FD layer selects stored derivatives while leaving mathematical expressions unchanged.
cfdD_alphadD_vetUdD_gridfunction_expressions() supplies the producer expressions;
register_cfdD_alphadD_vetUdD_gridfunctions() registers their output gridfunctions, and
register_CFunction_cfdD_alphadD_vetUdD_eval() emits the C function that writes them.
rhs_eval selects this storage with enable_cfdD_alphadD_vetUdD_gridfunctions=True.
The BHaH BSSN examples switch the scheme on with
enable_cfdD_alphadD_vetUdD_gridfunctions_for_GPU, off by default, and then call
cfdD_alphadD_vetUdD_eval before rhs_eval within each Method of Lines substep.

Unlike hDDdD, which Ricci_eval consumes (hDDdD_eval.py), these gridfunctions cannot live in
the SCRATCH group: rhs_eval reads them with a stencil while writing its own output to the
Method of Lines buffer, so a pointwise store there would overwrite a neighbor's stencil
point. They are ordinary AUXEVOL gridfunctions, which rhs_eval already receives, and so
they cost memory.

Only the directions that a mixed second derivative differentiates are stored. With the
mixed index pair canonicalized to j < k, partial_0 f is differentiated in i1 and i2 and
partial_1 f in i2, so each stored direction is produced over the interior grown by
fd_order/2 points in those directions. A component that a symmetry axis has zeroed is
neither stored nor grown.

Ten gridfunctions rather than five: storing partial_0 alone would rebuild ten of the
fifteen mixed derivatives and remove 10*74 + 5*9 - 10*9 = 695 double-precision operations
per point for five gridfunctions, while storing partial_1 as well rebuilds all fifteen and
removes 15*74 + 10*9 - 15*9 = 1065 for ten. The second five gridfunctions therefore buy 370
operations per point, 7% of the 5166 that rhs_eval evaluates, for the second half of the
memory. The ten-gridfunction variant is the one measured in
nrpy/examples/blackhole_spectroscopy.py.

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


def register_cfdD_alphadD_vetUdD_gridfunctions() -> List[str]:
    """
    Register the AUXEVOL gridfunctions cfdD, alphadD and vetUdD if not already registered.

    Both the producer and its consumer call this, because parallel code generation runs each
    registration in its own worker.

    The components are registered by name, so each carries rank 0 and therefore scalar parity.
    That is correct only because nothing ever applies boundary conditions to them: the producer
    computes every point rhs_eval reads, ghost zones included, so these gridfunctions must never
    be added to an inner-boundary synchronization list. Declaring them rank 1 would not fix the
    parity either, since vetUdD is a mixed rank-2 object.

    :return: Stored gridfunction names that RHS kernels may select during FD lowering.
    """
    names = list(cfdD_alphadD_vetUdD_gridfunction_expressions())
    if names and names[0] not in gri.glb_gridfcs_dict:
        _ = gri.register_gridfunctions(names, group="AUXEVOL", is_basename=False)
    return names


def cfdD_alphadD_vetUdD_gridfunction_expressions() -> Dict[str, sp.Expr]:
    """
    Construct the first derivatives written to the AUXEVOL gridfunctions.

    Store only lower directions of nonzero mixed second derivatives. FD lowering
    selects these names explicitly; unmixed and upwind derivatives retain their
    original finite-difference stencils.

    :return: Stored gridfunction names mapped to the expressions that produce them.

    Doctests:
    >>> saved_symmetry_axes = par.parval_from_str("symmetry_axes")
    >>> try:
    ...     par.set_parval_from_str("symmetry_axes", "")
    ...     full = cfdD_alphadD_vetUdD_gridfunction_expressions()
    ...     par.set_parval_from_str("symmetry_axes", "1")
    ...     axis1 = cfdD_alphadD_vetUdD_gridfunction_expressions()
    ...     par.set_parval_from_str("symmetry_axes", "02")
    ...     axes02 = cfdD_alphadD_vetUdD_gridfunction_expressions()
    ... finally:
    ...     par.set_parval_from_str("symmetry_axes", saved_symmetry_axes)
    >>> full["cfdD0"], full["vetUdD21"], len(full)
    (cf_dD0, vetU_dD21, 10)
    >>> sorted(axis1), axes02
    (['alphadD0', 'cfdD0', 'vetUdD00', 'vetUdD10', 'vetUdD20'], {})
    """
    cf_dDD = ixp.declarerank2("cf_dDD", symmetry="sym01")
    vetU_dD = ixp.declarerank2("vetU_dD")
    fields = [(scalar, ixp.declarerank1(f"{scalar}_dD")) for scalar in ("cf", "alpha")]
    expressions: Dict[str, sp.Expr] = {}
    for basename, first in fields:
        for j in range(3):
            if any(cf_dDD[j][k] != 0 for k in range(j + 1, 3)):
                expressions[f"{basename}dD{j}"] = first[j]
    for i in range(3):
        for j in range(3):
            if any(cf_dDD[j][k] != 0 for k in range(j + 1, 3)):
                expressions[f"vetUdD{i}{j}"] = vetU_dD[i][j]
    return expressions


def register_CFunction_cfdD_alphadD_vetUdD_eval(
    CoordSystem: str,
    enable_intrinsics: bool,
    enable_fd_functions: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that stores the first derivatives of cf, alpha and vetU.

    The emitted cfdD_alphadD_vetUdD_eval(params, in_gfs, auxevol_gfs) launches one kernel
    per stored direction, each over the interior grown by the stencil radius in every
    direction rhs_eval differentiates that stored derivative in. It must run after the
    evolved gridfunctions in in_gfs are current and before rhs_eval reads auxevol_gfs
    within the same right-hand-side evaluation, and it pairs with
    register_CFunction_rhs_eval(..., enable_cfdD_alphadD_vetUdD_gridfunctions=True), which
    selects stored derivatives during finite-difference lowering.

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
    # generation runs each registration in its own worker, so the dependency is stated here,
    # exactly as rhs_eval states it.
    _ = BSSN_quantities[CoordSystem + "_rfm_precompute"]

    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += [
            str(
                Path("intrinsics") / "cuda_intrinsics.h"
                if is_cuda
                else Path("intrinsics") / "simd_intrinsics.h"
            )
        ]
    desc = r"""Set the first derivatives of cf, alpha and vetU."""
    cfunc_type = "void"
    name = "cfdD_alphadD_vetUdD_eval"
    arg_dict_cuda = {
        "in_gfs": "const REAL *restrict",
        "auxevol_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    params = ",".join([f"{v} {k}" for k, v in arg_dict_host.items()])

    register_cfdD_alphadD_vetUdD_gridfunctions()
    expressions = cfdD_alphadD_vetUdD_gridfunction_expressions()
    stored_directions = sorted({int(name[-1]) for name in expressions})
    cf_dDD = ixp.declarerank2("cf_dDD", symmetry="sym01")

    # c_codegen() clears the finite-difference helper registry at the start of every call,
    # so each direction's helpers are collected before the next call and emitted once.
    collected_FDFunctions = {}
    prefunc = ""
    launch_body = ""
    for direction in stored_directions:
        # The stored direction must be valid wherever rhs_eval differentiates it: the interior
        # grown by the stencil radius in every transverse direction that some non-zero mixed
        # second derivative differentiates it in.
        grown = [
            f"i{transverse}"
            for transverse in range(direction + 1, 3)
            if cf_dDD[direction][transverse] != 0
        ]
        loop_region = (
            "interior plus stencil halo in " + " ".join(grown) if grown else "interior"
        )

        exprs = [
            expr for gf, expr in expressions.items() if gf.endswith(str(direction))
        ]
        access_gfs = [
            f"auxevol_gfs[IDX4({str(gf).upper()}GF, i0, i1, i2)]"
            for gf in expressions
            if gf.endswith(str(direction))
        ]

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
            launchblock_with_braces=False,
            launch_dict={
                **BHaH.parallelization.cuda_utilities.default_launch_dictionary,
                "threads_per_block": ["64", "1", "1"],
            },
            thread_tiling_macro_suffix="CFDD_ALPHADD_VETUDD_EVAL",
        )
        prefunc += kernel
        # Keep CUDA launch locals in independent, semantically marked scopes.
        if is_cuda:
            launch_body += (
                f"{{\n{direction_launch_body}\n}} "
                f"// END BLOCK: Launch {name} direction {direction}\n"
            )
        else:
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

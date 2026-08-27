"""
Generate C function enforcing conformal metric determinant and A trace constraints.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Union, cast

import sympy as sp
from mpmath import mpc, mpf  # type: ignore

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
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
    validate_expressions: bool = False,
) -> Union[None, Dict[str, Union[mpc, mpf]], pcg.NRPyEnv_type]:
    """
    Register combined determinant and trace-free conformal-tensor projection.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param validate_expressions: Whether to validate generated SymPy expressions against trusted values.

    :return: None in registration phase, processed expressions when validating, otherwise the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    Bq = BSSN_quantities[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]
    is_general_rfm = CoordSystem.startswith("GeneralRFM")

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

    hprimeDD = ixp.zerorank2()
    aprimeDD = ixp.zerorank2()
    if is_general_rfm:
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        _, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
        nrpyAbs = sp.Function("nrpyAbs")
        q = (nrpyAbs(rfm.detgammahat) / detgammabar) ** sp.Rational(1, 3)
        gammabarprimeDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammabarprimeDD[i][j] = q * Bq.gammabarDD[i][j]
        gammabarprimeUU, _ = ixp.symm_matrix_inverter3x3(gammabarprimeDD)
        trAbarprime = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                trAbarprime += gammabarprimeUU[i][j] * Bq.AbarDD[i][j]
        for i in range(3):
            for j in range(3):
                hprimeDD[i][j] = (gammabarprimeDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[
                    i
                ][j]
                aprimeDD[i][j] = (
                    Bq.AbarDD[i][j]
                    - sp.Rational(1, 3) * gammabarprimeDD[i][j] * trAbarprime
                ) / rfm.ReDD[i][j]
    else:
        # Every standard RFM is orthogonal. With diagonal scale-factor matrix D,
        # gammabar = D H D, gammahat = D I D, and Abar = D a D, where H = I+h.
        # Thus det(gammabar)/det(gammahat) = det(H), so q = det(H)^(-1/3)
        # and h' = q H - I. Also gammabar'^(-1) = q^(-1) D^(-1) H^(-1)
        # D^(-1), hence tr(gammabar'^(-1) Abar) = q^(-1) (H^(-1):a).
        # The q factors cancel in gammabar' tr/3, leaving a' = a-H(H^(-1):a)/3.
        HDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                HDD[i][j] = Bq.hDD[i][j] + sp.KroneckerDelta(i, j)
        HUU, detH = ixp.symm_matrix_inverter3x3(HDD)
        q = detH ** sp.Rational(-1, 3)
        trA = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                trA += HUU[i][j] * Bq.aDD[i][j]
        for i in range(3):
            for j in range(3):
                hprimeDD[i][j] = q * HDD[i][j] - sp.KroneckerDelta(i, j)
                aprimeDD[i][j] = Bq.aDD[i][j] - sp.Rational(1, 3) * HDD[i][j] * trA

    access_gfs: List[str] = []
    projected_exprs_dict: Dict[str, sp.Expr] = {}
    for basename, expressions in (("hDD", hprimeDD), ("aDD", aprimeDD)):
        for i in range(3):
            for j in range(i, 3):
                varname = f"{basename}{i}{j}"
                access_gfs.append(
                    gri.BHaHGridFunction.access_gf(
                        varname, 0, 0, 0, gf_array_name="in_gfs"
                    )
                )
                projected_exprs_dict[varname] = expressions[i][j]

    if validate_expressions:
        validation_exprs_dict = {
            varname: expr.subs(sp.Function("nrpyAbs"), sp.Abs)
            for varname, expr in projected_exprs_dict.items()
        }
        return ve.process_dictionary_of_expressions(
            validation_exprs_dict, fixed_mpfs_for_free_symbols=True
        )

    projected_exprs = list(projected_exprs_dict.values())

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
    import os

    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("enable_parallel_codegen", False)
    for validation_coord, validation_precompute in (
        ("Cartesian", False),
        ("SinhSpherical", True),
        ("GeneralRFM", True),
    ):
        results_dict = register_CFunction_enforce_detgbar_equals_detghat_trAzero(
            validation_coord,
            validation_precompute,
            False,
            1,
            validate_expressions=True,
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_coord}"
            + ("_rfm_precompute" if validation_precompute else ""),
            cast(Dict[str, Union[mpc, mpf]], results_dict),
        )

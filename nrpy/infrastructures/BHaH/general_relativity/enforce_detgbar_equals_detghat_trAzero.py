"""
Generate C function enforcing conformal metric determinant and A trace constraints.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
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
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register combined determinant and trace-free conformal-tensor projection.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
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
    projected_exprs: List[sp.Expr] = []
    for basename, expressions in (("hDD", hprimeDD), ("aDD", aprimeDD)):
        for i in range(3):
            for j in range(i, 3):
                access_gfs.append(
                    gri.BHaHGridFunction.access_gf(
                        f"{basename}{i}{j}", 0, 0, 0, gf_array_name="in_gfs"
                    )
                )
                projected_exprs.append(expressions[i][j])

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
    import nrpy.validate_expressions.validate_expressions as ve

    validation_HDD = [
        [sp.Integer(4), sp.Integer(1), sp.Rational(1, 2)],
        [sp.Integer(1), sp.Integer(3), sp.Rational(1, 4)],
        [sp.Rational(1, 2), sp.Rational(1, 4), sp.Integer(2)],
    ]
    validation_aDD = [
        [sp.Integer(2), sp.Integer(3), sp.Integer(5)],
        [sp.Integer(3), sp.Integer(7), sp.Integer(11)],
        [sp.Integer(5), sp.Integer(11), sp.Integer(13)],
    ]
    validation_HUU, validation_detH = ixp.symm_matrix_inverter3x3(validation_HDD)
    validation_q = validation_detH ** sp.Rational(-1, 3)
    validation_HprimeDD = [
        [validation_q * validation_HDD[i][j] for j in range(3)] for i in range(3)
    ]
    validation_trA = sum(
        validation_HUU[i][j] * validation_aDD[i][j] for i in range(3) for j in range(3)
    )
    validation_HprimeUU, validation_detHprime = ixp.symm_matrix_inverter3x3(
        validation_HprimeDD
    )
    validation_aprimeDD = [
        [
            validation_aDD[i][j]
            - sp.Rational(1, 3) * validation_HDD[i][j] * validation_trA
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_trAprime = sum(
        validation_HprimeUU[i][j] * validation_aprimeDD[i][j]
        for i in range(3)
        for j in range(3)
    )
    ve.assert_equal(validation_detHprime, sp.sympify(1), suppress_message=True)
    ve.assert_equal(validation_trAprime, sp.sympify(0), suppress_message=True)

    validation_scaleD = [sp.Integer(2), sp.Integer(3), sp.Integer(5)]
    validation_gammahatDD = [
        [
            sp.KroneckerDelta(i, j) * validation_scaleD[i] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_gammabarDD = [
        [
            validation_scaleD[i] * validation_HDD[i][j] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_AbarDD = [
        [
            validation_scaleD[i] * validation_aDD[i][j] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    _, validation_detgammahat = ixp.symm_matrix_inverter3x3(validation_gammahatDD)
    _, validation_detgammabar = ixp.symm_matrix_inverter3x3(validation_gammabarDD)
    validation_full_q = (
        validation_detgammahat / validation_detgammabar
    ) ** sp.Rational(1, 3)
    validation_gammabarprimeDD = [
        [validation_full_q * validation_gammabarDD[i][j] for j in range(3)]
        for i in range(3)
    ]
    validation_gammabarprimeUU, _ = ixp.symm_matrix_inverter3x3(
        validation_gammabarprimeDD
    )
    validation_full_trA = sum(
        validation_gammabarprimeUU[i][j] * validation_AbarDD[i][j]
        for i in range(3)
        for j in range(3)
    )
    validation_full_hprimeDD = [
        [
            (validation_gammabarprimeDD[i][j] - validation_gammahatDD[i][j])
            / (validation_scaleD[i] * validation_scaleD[j])
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_full_aprimeDD = [
        [
            (
                validation_AbarDD[i][j]
                - sp.Rational(1, 3)
                * validation_gammabarprimeDD[i][j]
                * validation_full_trA
            )
            / (validation_scaleD[i] * validation_scaleD[j])
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_reduced_hprimeDD = [
        [validation_HprimeDD[i][j] - sp.KroneckerDelta(i, j) for j in range(3)]
        for i in range(3)
    ]
    ve.assert_equal(
        {
            "hprimeDD": validation_full_hprimeDD,
            "aprimeDD": validation_full_aprimeDD,
        },
        {
            "hprimeDD": validation_reduced_hprimeDD,
            "aprimeDD": validation_aprimeDD,
        },
        suppress_message=True,
    )

    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("enable_parallel_codegen", False)
    for validation_coord, validation_precompute in (
        ("Cartesian", False),
        ("SinhSpherical", True),
        ("GeneralRFM", True),
    ):
        cfc.CFunction_dict.clear()
        register_CFunction_enforce_detgbar_equals_detghat_trAzero(
            validation_coord, validation_precompute, False, 1
        )
        validation_cfunc = next(iter(cfc.CFunction_dict.values()))
        generated = (validation_cfunc.prefunc or "") + (validation_cfunc.body or "")
        if validation_coord.startswith("GeneralRFM"):
            for validation_i, validation_j in (
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 1),
                (1, 2),
                (2, 2),
            ):
                if f"const REAL ghatDD{validation_i}{validation_j} =" not in generated:
                    raise AssertionError(
                        "GeneralRFM projection must read all six ghat components"
                    )
        elif (
            "ghatDD" in generated
            or "rfmstruct->" in generated
            or any(f"const REAL xx{axis} =" in generated for axis in range(3))
        ):
            raise AssertionError(
                "standard projection must not evaluate coordinate or reference-metric data"
            )
        if validation_precompute:
            if "const rfm_struct *restrict rfmstruct" not in validation_cfunc.params:
                raise AssertionError("precomputed-RFM function signature changed")
        elif "const REAL *restrict x0" not in validation_cfunc.params:
            raise AssertionError("analytic-RFM function signature changed")
        if generated.count("Step 1 of 2") != 1 or generated.count("Step 2 of 2") != 1:
            raise AssertionError(
                "projection must contain one generated twelve-output block"
            )
        if generated.count("END LOOP: for i2 over") != 1:
            raise AssertionError("projection must contain one all-points loop nest")
        write_section = generated.split("Step 2 of 2", maxsplit=1)[1]
        if write_section.count("in_gfs[IDX4(") != 12:
            raise AssertionError("projection must store twelve outputs in one block")
        first_store = generated.index("in_gfs[IDX4(", generated.index("Step 2 of 2"))
        for validation_basename in ("aDD", "hDD"):
            for validation_i, validation_j in (
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 1),
                (1, 2),
                (2, 2),
            ):
                load = f"const REAL {validation_basename}{validation_i}{validation_j} ="
                if load not in generated or generated.index(load) > first_store:
                    raise AssertionError(
                        "projection must load all raw inputs before stores"
                    )

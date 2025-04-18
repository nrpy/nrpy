"""
Generate C code for computing the RHS of the BSSN equations in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict as ODict
from inspect import currentframe as cfr
from pathlib import Path
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
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def register_CFunction_rhs_eval(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_RbarDD_gridfunctions: bool,
    enable_T4munu: bool,
    enable_intrinsics: bool,
    enable_fd_functions: bool,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    enable_KreissOliger_dissipation: bool,
    KreissOliger_strength_gauge: float = 0.3,
    KreissOliger_strength_nongauge: float = 0.3,
    enable_CAKO: bool = False,
    enable_CAHD: bool = False,
    enable_SSL: bool = False,
    OMP_collapse: int = 1,
    validate_expressions: bool = False,
) -> Union[None, Dict[str, Union[mpf, mpc]], pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_RbarDD_gridfunctions: Whether to enable RbarDD gridfunctions.
    :param enable_T4munu: Whether to enable T4munu (stress-energy terms).
    :param enable_intrinsics: Whether to enable SIMD (Single Instruction, Multiple Data).
    :param enable_fd_functions: Whether to enable finite difference functions.
    :param LapseEvolutionOption: Lapse evolution equation choice.
    :param ShiftEvolutionOption: Lapse evolution equation choice.
    :param enable_KreissOliger_dissipation: Whether to enable Kreiss-Oliger dissipation.
    :param KreissOliger_strength_gauge: Gauge strength for Kreiss-Oliger dissipation.
    :param KreissOliger_strength_nongauge: Non-gauge strength for Kreiss-Oliger dissipation.
    :param enable_CAKO: Whether to enable curvature-aware Kreiss-Oliger dissipation (multiply strength by W).
    :param enable_CAHD: Whether to enable curvature-aware Hamiltonian-constraint damping.
    :param enable_SSL: Whether to enable slow-start lapse.
    :param OMP_collapse: Degree of OpenMP loop collapsing.
    :param validate_expressions: Whether to validate generated sympy expressions against trusted values.

    :raises ValueError: If EvolvedConformalFactor_cf not set to a supported value: {phi, chi, W}.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h"]
    if enable_intrinsics:
        includes += [
            str(
                Path("intrinsics") / "cuda_intrinsics.h"
                if parallelization == "cuda"
                else Path("intrinsics") / "simd_intrinsics.h"
            )
        ]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    cfunc_type = "void"
    name = "rhs_eval"
    arg_dict_cuda = {
        "auxevol_gfs": "const REAL *restrict",
        "in_gfs": "const REAL *restrict",
        "rhs_gfs": "REAL *restrict",
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
        "commondata": "const commondata_struct *restrict",
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
    }
    params = ",".join([f"{v} {k}" for k, v in arg_dict_host.items()])

    # Populate BSSN rhs variables
    rhs = BSSN_RHSs[
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        + ("_T4munu" if enable_T4munu else "")
    ]
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_T4munu=enable_T4munu,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
    )
    rhs.BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] = alpha_rhs
    for i in range(3):
        rhs.BSSN_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] = vet_rhsU[i]
        rhs.BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] = bet_rhsU[i]

    # local_BSSN_RHSs_varname_to_expr_dict is modified below if e.g., we add KO terms;
    #    DO NOT MODIFY rhs.BSSN_RHSs_varname_to_expr_dict!
    local_BSSN_RHSs_varname_to_expr_dict = rhs.BSSN_RHSs_varname_to_expr_dict.copy()
    local_BSSN_RHSs_varname_to_expr_dict = ODict(
        sorted(local_BSSN_RHSs_varname_to_expr_dict.items())
    )

    # Define conformal factor W.
    Bq = BSSN_quantities[
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
    ]
    EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    if EvolvedConformalFactor_cf == "W":
        W = Bq.cf
    elif EvolvedConformalFactor_cf == "chi":
        W = sp.sqrt(Bq.cf)
    elif EvolvedConformalFactor_cf == "phi":
        W = sp.exp(-2 * Bq.cf)
    else:
        raise ValueError(
            "Error: only EvolvedConformalFactor_cf = (W or chi or phi) supported."
        )

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    if enable_KreissOliger_dissipation:
        diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
            "REAL",
            __name__,
            ["KreissOliger_strength_gauge", "KreissOliger_strength_nongauge"],
            [KreissOliger_strength_gauge, KreissOliger_strength_nongauge],
            commondata=True,
        )

        # vvv BEGIN CAKO vvv
        if enable_CAKO:
            diss_strength_gauge *= W
            diss_strength_nongauge *= W
        # ^^^ END CAKO ^^^

        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        alpha_dKOD = ixp.declarerank1("alpha_dKOD")
        cf_dKOD = ixp.declarerank1("cf_dKOD")
        trK_dKOD = ixp.declarerank1("trK_dKOD")
        betU_dKOD = ixp.declarerank2("betU_dKOD", symmetry="nosym")
        vetU_dKOD = ixp.declarerank2("vetU_dKOD", symmetry="nosym")
        lambdaU_dKOD = ixp.declarerank2("lambdaU_dKOD", symmetry="nosym")
        aDD_dKOD = ixp.declarerank3("aDD_dKOD", symmetry="sym01")
        hDD_dKOD = ixp.declarerank3("hDD_dKOD", symmetry="sym01")
        for k in range(3):
            local_BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] += (
                diss_strength_gauge * alpha_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            local_BSSN_RHSs_varname_to_expr_dict["cf_rhs"] += (
                diss_strength_nongauge * cf_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            local_BSSN_RHSs_varname_to_expr_dict["trK_rhs"] += (
                diss_strength_nongauge * trK_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftEvolutionOption:
                    local_BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] += (
                        diss_strength_gauge * betU_dKOD[i][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                local_BSSN_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] += (
                    diss_strength_gauge * vetU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                local_BSSN_RHSs_varname_to_expr_dict[f"lambda_rhsU{i}"] += (
                    diss_strength_nongauge * lambdaU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(i, 3):
                    local_BSSN_RHSs_varname_to_expr_dict[f"a_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * aDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    local_BSSN_RHSs_varname_to_expr_dict[f"h_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * hDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    # vvv BEGIN CAHD vvv
    if enable_CAHD:
        Bcon = BSSN_constraints[
            CoordSystem
            + ("_rfm_precompute" if enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
            + ("_T4munu" if enable_T4munu else "")
        ]
        if "cahdprefactor" not in gri.glb_gridfcs_dict:
            _ = gri.register_gridfunctions(
                "cahdprefactor",
                group="AUXEVOL",
                gf_array_name="auxevol_gfs",
            )
        _C_CAHD = par.register_CodeParameter(
            "REAL", __name__, "C_CAHD", 0.15, commondata=True, add_to_parfile=True
        )
        # Initialize CAHD_term assuming phi is the evolved conformal factor. CFL_FACTOR is defined in MoL.
        # CAHD_term = -C_CAHD * (sp.symbols("CFL_FACTOR") * sp.symbols("dsmin")) * Bcon.H
        # -> cahdprefactor = C_CAHD * sp.symbols("CFL_FACTOR") * sp.symbols("dsmin")
        CAHD_term = -1 * sp.symbols("cahdprefactor") * Bcon.H
        if EvolvedConformalFactor_cf == "phi":
            pass  # CAHD_term already assumes phi is the evolved conformal factor.
        elif EvolvedConformalFactor_cf == "W":
            # \partial_t W = \partial_t e^{-2 phi} = -2 W \partial_t phi
            CAHD_term *= -2 * Bq.cf
        elif EvolvedConformalFactor_cf == "chi":
            # \partial_t chi = \partial_t e^{-4 phi} = -4 chi \partial_t phi
            CAHD_term *= -4 * Bq.cf
        else:
            raise ValueError(
                "Error: only EvolvedConformalFactor_cf = (W or chi or phi) supported."
            )
        local_BSSN_RHSs_varname_to_expr_dict["cf_rhs"] += CAHD_term
    # ^^^ END CAHD ^^^

    # vvv BEGIN SSL vvv
    if enable_SSL:
        SSL_Gaussian_prefactor = par.register_CodeParameter(
            "REAL",
            __name__,
            "SSL_Gaussian_prefactor",
            1.0,
            commondata=True,
            add_to_parfile=False,
        )
        _SSL_h, _SSL_sigma = par.register_CodeParameters(
            "REAL",
            __name__,
            ["SSL_h", "SSL_sigma"],
            [0.6, 20.0],
            commondata=True,
            add_to_parfile=True,
        )
        local_BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] -= (
            W * SSL_Gaussian_prefactor * (Bq.alpha - W)
        )
    # ^^^ END SSL ^^^

    BSSN_RHSs_access_gf: List[str] = []
    for var in local_BSSN_RHSs_varname_to_expr_dict.keys():
        BSSN_RHSs_access_gf += [
            gri.BHaHGridFunction.access_gf(
                var.replace("_rhs", ""),
                0,
                0,
                0,
                gf_array_name="rhs_gfs",
            )
        ]
    # Set up upwind control vector (betaU)
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = ixp.zerorank1()
    vetU = ixp.declarerank1("vetU")
    for i in range(3):
        # self.lambda_rhsU[i] = self.Lambdabar_rhsU[i] / rfm.ReU[i]
        betaU[i] = vetU[i] * rfm.ReU[i]

    # Perform validation of BSSN_RHSs against trusted version.
    if validate_expressions:
        return ve.process_dictionary_of_expressions(
            local_BSSN_RHSs_varname_to_expr_dict, fixed_mpfs_for_free_symbols=True
        )
    # ve.compare_or_generate_trusted_results(
    #     os.path.abspath(__file__),
    #     os.getcwd(),
    #     # File basename. If this is set to "trusted_module_test1", then
    #     #   trusted results_dict will be stored in tests/trusted_module_test1.py
    #     f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolutionOption}_{ShiftEvolutionOption}_{CoordSystem}_T4munu{enable_T4munu}_KO{enable_KreissOliger_dissipation}",
    #     cast(Dict[str, Union[mpf, mpc]], results_dict),
    # )

    expr_list = list(local_BSSN_RHSs_varname_to_expr_dict.values())

    # Find symbols stored in params
    param_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=[f"xx{j}" for j in range(3)]
    )

    arg_dict_cuda = {
        **arg_dict_cuda,
        **{
            f"NOCUDA{k}" if enable_intrinsics else k: "const REAL"
            for k in commondata_symbols
        },
    }

    arg_dict_host = {
        "params": "const params_struct *restrict",
        **{k.replace("CUDA", "SIMD"): v for k, v in arg_dict_cuda.items()},
    }

    kernel_body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            expr_list,
            BSSN_RHSs_access_gf,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            upwind_control_vec=betaU,
            enable_fd_functions=enable_fd_functions,
            rational_const_alias=(
                "static constexpr" if parallelization == "cuda" else "static const"
            ),
        ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    loop_params = parallel_utils.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )
    if enable_intrinsics:
        for symbol in commondata_symbols:
            loop_params += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {symbol} = ConstSIMD(NOSIMD{symbol});\n"

    params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=enable_intrinsics,
        var_access=parallel_utils.get_params_access(parallelization),
    )
    kernel_body = f"{loop_params}\n{params_definitions}\n{kernel_body}"

    kernel, launch_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body.replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="BSSN_RHS",
    )

    for symbol in commondata_symbols:
        tmp_sym = (
            f"NOCUDA{symbol}"
            if parallelization == "cuda" and enable_intrinsics
            else (
                f"NOSIMD{symbol}"
                if enable_intrinsics
                else symbol if enable_intrinsics else symbol
            )
        )
        launch_body = launch_body.replace(tmp_sym, f"commondata->{symbol}")

    prefunc = ""
    if parallelization == "cuda" and enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc(
            cfunc_decorators="__device__ "
        ).replace("SIMD", "CUDA")
    elif enable_fd_functions:
        prefunc = fin.construct_FD_functions_prefunc()

    cfc.register_CFunction(
        include_CodeParameters_h=False,
        includes=includes,
        prefunc=prefunc + kernel,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=launch_body,
        enable_simd=enable_intrinsics,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


if __name__ == "__main__":
    import os

    Coord = "SinhSpherical"
    LapseEvolOption = "OnePlusLog"
    ShiftEvolOption = "GammaDriving2ndOrder_Covariant"
    for Rbar_gfs in [True, False]:
        for T4munu_enable in [True, False]:
            for enable_Improvements in [True, False]:
                results_dict = register_CFunction_rhs_eval(
                    CoordSystem=Coord,
                    enable_rfm_precompute=True,
                    enable_RbarDD_gridfunctions=Rbar_gfs,
                    enable_T4munu=T4munu_enable,
                    enable_intrinsics=False,
                    enable_fd_functions=False,
                    enable_KreissOliger_dissipation=True,
                    LapseEvolutionOption=LapseEvolOption,
                    ShiftEvolutionOption=ShiftEvolOption,
                    enable_CAKO=enable_Improvements,
                    enable_CAHD=enable_Improvements,
                    enable_SSL=enable_Improvements,
                    validate_expressions=True,
                )
                ve.compare_or_generate_trusted_results(
                    os.path.abspath(__file__),
                    os.getcwd(),
                    # File basename. If this is set to "trusted_module_test1", then
                    #   trusted results_dict will be stored in tests/trusted_module_test1.py
                    f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolOption}_{ShiftEvolOption}_{Coord}_Rbargfs{Rbar_gfs}_T4munu{T4munu_enable}_Improvements{enable_Improvements}",
                    cast(Dict[str, Union[mpf, mpc]], results_dict),
                )

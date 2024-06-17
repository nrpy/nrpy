"""
Generates function to compute the right-hand sides of the BSSN evolution equations.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

import re
from collections import OrderedDict as ODict
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Any, Dict, List, Union, cast

import sympy as sp
from mpmath import mpc, mpf  # type: ignore

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.CarpetX.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.infrastructures.CarpetX.CarpetX_include_header import define_standard_includes


def register_CFunction_rhs_eval(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_T4munu: bool,
    enable_simd: bool,
    fd_order: int,
    LapseEvolutionOption: str,
    ShiftEvolutionOption: str,
    enable_KreissOliger_dissipation: bool,
    KreissOliger_strength_mult_by_W: bool = False,
    # when mult by W, strength_gauge=0.99 & strength_nongauge=0.3 is best.
    KreissOliger_strength_gauge: float = 0.1,
    KreissOliger_strength_nongauge: float = 0.1,
    fp_type: str = "double",
    validate_expressions: bool = False,
) -> Union[None, Dict[str, Union[mpf, mpc]], pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_T4munu: Whether to include the stress-energy tensor. Defaults to False.
    :param enable_simd: Whether to enable SIMD (Single Instruction, Multiple Data).
    :param fd_order: Order of finite difference method
    :param LapseEvolutionOption: Lapse evolution equation choice.
    :param ShiftEvolutionOption: Lapse evolution equation choice.
    :param enable_KreissOliger_dissipation: Whether to enable Kreiss-Oliger dissipation.
    :param KreissOliger_strength_mult_by_W: Whether to multiply Kreiss-Oliger strength by W.
    :param KreissOliger_strength_gauge: Gauge strength for Kreiss-Oliger dissipation.
    :param KreissOliger_strength_nongauge: Non-gauge strength for Kreiss-Oliger dissipation.
    :param fp_type: Floating point type, e.g., "double".
    :param validate_expressions: Whether to validate generated sympy expressions against trusted values.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    old_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    enable_RbarDD_gridfunctions = True

    includes = define_standard_includes()
    if enable_simd:
        includes += [("./simd/simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    name = f"{thorn_name}_rhs_eval_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTSX_{name};
"""
    if enable_simd:
        body += f"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));
  const CCTK_REAL *param_PI CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("PI", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY PI CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_PI);
  const CCTK_REAL *param_eta CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("eta", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY eta CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_eta);
  const CCTK_REAL *param_diss_strength CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("diss_strength", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY diss_strength CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_diss_strength);

"""
    else:
        body += """  const CCTK_REAL invdxx0 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdxx1 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdxx2 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(2);
  DECLARE_CCTK_PARAMETERS;

  #define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0

"""

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

    rhs.BSSN_RHSs_varname_to_expr_dict = ODict(
        sorted(rhs.BSSN_RHSs_varname_to_expr_dict.items())
    )

    # Add Kreiss-Oliger dissipation to the BSSN RHSs:
    if enable_KreissOliger_dissipation:
        diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
            "CCTK_REAL",
            __name__,
            ["diss_strength", "diss_strength"],
            [KreissOliger_strength_gauge, KreissOliger_strength_nongauge],
            commondata=True,
        )

        if KreissOliger_strength_mult_by_W:
            Bq = BSSN_quantities[
                CoordSystem
                + ("_rfm_precompute" if enable_rfm_precompute else "")
                + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
            ]
            EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
            if EvolvedConformalFactor_cf == "W":
                diss_strength_gauge *= Bq.cf
                diss_strength_nongauge *= Bq.cf
            elif EvolvedConformalFactor_cf == "chi":
                diss_strength_gauge *= sp.sqrt(Bq.cf)
                diss_strength_nongauge *= sp.sqrt(Bq.cf)
            elif EvolvedConformalFactor_cf == "phi":
                diss_strength_gauge *= sp.exp(-2 * Bq.cf)
                diss_strength_nongauge *= sp.exp(-2 * Bq.cf)

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
            rhs.BSSN_RHSs_varname_to_expr_dict["alpha_rhs"] += (
                diss_strength_gauge * alpha_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.BSSN_RHSs_varname_to_expr_dict["cf_rhs"] += (
                diss_strength_nongauge * cf_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            rhs.BSSN_RHSs_varname_to_expr_dict["trK_rhs"] += (
                diss_strength_nongauge * trK_dKOD[k] * rfm.ReU[k]
            )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
            for i in range(3):
                if "2ndOrder" in ShiftEvolutionOption:
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] += (
                        diss_strength_gauge * betU_dKOD[i][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.BSSN_RHSs_varname_to_expr_dict[f"vet_rhsU{i}"] += (
                    diss_strength_gauge * vetU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                rhs.BSSN_RHSs_varname_to_expr_dict[f"lambda_rhsU{i}"] += (
                    diss_strength_nongauge * lambdaU_dKOD[i][k] * rfm.ReU[k]
                )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                for j in range(i, 3):
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"a_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * aDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
                    rhs.BSSN_RHSs_varname_to_expr_dict[f"h_rhsDD{i}{j}"] += (
                        diss_strength_nongauge * hDD_dKOD[i][j][k] * rfm.ReU[k]
                    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]

    BSSN_RHSs_access_gf: List[str] = []
    pattern = re.compile(r"([a-zA-Z_]+)_rhs([a-zA-Z_]+[0-2]+)")
    for var in rhs.BSSN_RHSs_varname_to_expr_dict.keys():
        patmatch = pattern.match(var)
        if patmatch:
            var_name = patmatch.group(1) + patmatch.group(2) + "_rhs"
        else:
            var_name = var
        BSSN_RHSs_access_gf += [gri.CarpetXGridFunction.access_gf(gf_name=var_name)]

    # Set up upwind control vector (betaU)
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = ixp.zerorank1()
    vetU = ixp.declarerank1("vetU")
    for i in range(3):
        betaU[i] = vetU[i] * rfm.ReU[i]
    if validate_expressions:
        return ve.process_dictionary_of_expressions(
            rhs.BSSN_RHSs_varname_to_expr_dict, fixed_mpfs_for_free_symbols=True
        )
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list(rhs.BSSN_RHSs_varname_to_expr_dict.values()),
            BSSN_RHSs_access_gf,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            upwind_control_vec=betaU,
            enable_fd_functions=True,
            enable_GoldenKernels=True,
            fp_type=fp_type,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in ODESolvers_RHS as {thorn_name}_RHS after {thorn_name}_Ricci
  {{
    LANG: C
    READS:  evol_variables(everywhere),
            auxevol_variables(interior)
    WRITES: evol_variables_rhs(interior)
  }} "Evaluate BSSN RHSs, at finite-differencing order {fd_order}"
}}
"""

    params = ["eta", "diss_strength", "FD_order"]
    if thorn_name == "Baikal":
        params += ["PI"]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        prefunc=fin.construct_FD_functions_prefunc(),
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("ODESolvers_RHS", schedule)],
        ET_current_thorn_CodeParams_used=params,
    )

    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


if __name__ == "__main__":
    import os

    Coord = "Cartesian"
    LapseEvolOption = "OnePlusLog"
    ShiftEvolOption = "GammaDriving2ndOrder_Covariant"
    for T4munu_enable in [True, False]:
        results_dict = register_CFunction_rhs_eval(
            thorn_name="dummy_thorn_name",
            CoordSystem=Coord,
            enable_rfm_precompute=False,
            enable_T4munu=T4munu_enable,
            enable_simd=False,
            fd_order=4,  # unused for this validation.
            enable_KreissOliger_dissipation=True,
            LapseEvolutionOption=LapseEvolOption,
            ShiftEvolutionOption=ShiftEvolOption,
            validate_expressions=True,
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolOption}_{ShiftEvolOption}_{Coord}_T4munu{T4munu_enable}",
            cast(dict[Any, Any | Any], results_dict),
        )

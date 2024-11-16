"""
Generates function to compute the right-hand sides of the BSSN evolution equations.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

import os
import re
from collections import OrderedDict as ODict
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
import nrpy.indexedexp as ixp
import nrpy.infrastructures.CarpetX.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
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
    KreissOliger_strength_gauge: float = 0.1,
    KreissOliger_strength_nongauge: float = 0.1,
    enable_CAKO: bool = False,
    enable_CAHD: bool = False,
    enable_SSL: bool = False,
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
    :param KreissOliger_strength_gauge: Gauge strength for Kreiss-Oliger dissipation.
    :param KreissOliger_strength_nongauge: Non-gauge strength for Kreiss-Oliger dissipation.
    :param enable_CAKO: Whether to enable curvature-aware Kreiss-Oliger dissipation (multiply strength by W).
    :param enable_CAHD: Whether to enable curvature-aware Hamiltonian-constraint damping.
    :param enable_SSL: Whether to enable slow-start lapse.
    :param validate_expressions: Whether to validate generated sympy expressions against trusted values.

    :raises ValueError: If EvolvedConformalFactor_cf not set to a supported value: {phi, chi, W}.

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
        includes += [("./intrinsics/simd_intrinsics.h")]
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
"""
        if enable_CAKO:
            body += f"""
  const CCTK_REAL *param_diss_strength_gauge CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("diss_strength_gauge", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY diss_strength_gauge CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_diss_strength_gauge);
  const CCTK_REAL *param_diss_strength_nongauge CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("diss_strength_nongauge", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY diss_strength_nongauge CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_diss_strength_nongauge);
"""
        else:
            body += f"""
  const CCTK_REAL *param_diss_strength CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("diss_strength", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY diss_strength CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_diss_strength);
"""
        if enable_CAHD:
            body += f"""
  const CCTK_REAL *param_C_CAHD CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("C_CAHD", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY C_CAHD CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_C_CAHD);
  // cahdprefactor = C_CAHD * sp.symbols("CFL_FACTOR") * sp.symbols("dsmin")
  const REAL_SIMD_ARRAY cahdprefactor CCTK_ATTRIBUTE_UNUSED = MulSIMD(C_CAHD, MulSIMD(,));
"""
        if enable_SSL:
            body += f"""
    const CCTK_REAL *SSL_h CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("SSL_h", "{thorn_name}", NULL);
    const CCTK_REAL *SSL_sigma CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("SSL_sigma", "{thorn_name}", NULL);

    const CCTK_REAL noSIMD_SSL_Gaussian_prefactor CCTK_ATTRIBUTE_UNUSED = *SSL_h * exp(-CCTK_TIME * CCTK_TIME / (2 * (*SSL_sigma) * (*SSL_sigma)));
    const REAL_SIMD_ARRAY SSL_Gaussian_prefactor CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*noSIMD_SSL_Gaussian_prefactor);
"""
    else:
        body += """  const CCTK_REAL invdxx0 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(0);
  const CCTK_REAL invdxx1 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(1);
  const CCTK_REAL invdxx2 CCTK_ATTRIBUTE_UNUSED = 1.0/CCTK_DELTA_SPACE(2);
  DECLARE_CCTK_PARAMETERS;

  #define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
"""
        if enable_CAHD:
            body += """
  // cahdprefactor = C_CAHD * sp.symbols("CFL_FACTOR") * sp.symbols("dsmin")
  const CCTK_REAL cahdprefactor = C_CAHD * ;
"""
        if enable_SSL:
            body += """
  const CCTK_REAL SSL_Gaussian_prefactor CCTK_ATTRIBUTE_UNUSED = SSL_h * exp(-CCTK_TIME * CCTK_TIME / (2 * (SSL_sigma) * (SSL_sigma)));
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
        if enable_CAKO:
            diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
                "CCTK_REAL",
                __name__,
                ["diss_strength_gauge", "diss_strength_nongauge"],
                [KreissOliger_strength_gauge, KreissOliger_strength_nongauge],
                commondata=True,
            )
        else:
            diss_strength_gauge, diss_strength_nongauge = par.register_CodeParameters(
                "CCTK_REAL",
                __name__,
                ["diss_strength", "diss_strength"],
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
    pattern = re.compile(r"([a-zA-Z_]+)_rhs([a-zA-Z_]+[0-2]+)")
    for var in local_BSSN_RHSs_varname_to_expr_dict.keys():
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

    # Perform validation of BSSN_RHSs against trusted version.
    results_dictionary = ve.process_dictionary_of_expressions(
        local_BSSN_RHSs_varname_to_expr_dict, fixed_mpfs_for_free_symbols=True
    )
    if validate_expressions:
        return results_dictionary
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        # File basename. If this is set to "trusted_module_test1", then
        #   trusted results_dict will be stored in tests/trusted_module_test1.py
        f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolutionOption}_{ShiftEvolutionOption}_{CoordSystem}_T4munu{enable_T4munu}_KO{enable_KreissOliger_dissipation}_improvements{enable_SSL}",
        results_dictionary,
    )

    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list(local_BSSN_RHSs_varname_to_expr_dict.values()),
            BSSN_RHSs_access_gf,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            upwind_control_vec=betaU,
            enable_fd_functions=True,
            enable_GoldenKernels=True,
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

    params = ["eta", "FD_order"]
    if enable_CAKO:
        params += ["diss_strength_gauge", "diss_strength_nongauge"]
    else:
        params += ["diss_strength"]
    if enable_SSL:
        params += ["SSL_h", "SSL_sigma"]
    if enable_CAHD:
        params += [
            "C_CAHD",
            "CFL_FACTOR__ignore_repeats_Carpet_timeref_factors",
        ]
    if thorn_name == "Baikal":
        params += ["PI"]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        prefunc=fin.construct_FD_functions_prefunc().replace(
            "NO_INLINE", "CCTK_ATTRIBUTE_NOINLINE"
        ),  # This prevents a hang when compiling higher-order FD kernels with certain versions of GCC. I'd prefer not adjusting construct_FD_functions_prefunc() for just this infrastructure.
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("ODESolvers_RHS", schedule)],
        ET_current_thorn_CodeParams_used=params,
    )

    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


if __name__ == "__main__":
    Coord = "Cartesian"
    LapseEvolOption = "OnePlusLog"
    ShiftEvolOption = "GammaDriving2ndOrder_Covariant"
    for T4munu_enable in [True, False]:
        for improvements_enable in [True, False]:
            results_dict = register_CFunction_rhs_eval(
                thorn_name="dummy_thorn_name",
                CoordSystem=Coord,
                enable_rfm_precompute=False,
                enable_T4munu=T4munu_enable,
                enable_simd=False,
                fd_order=4,  # unused for this validation.
                LapseEvolutionOption=LapseEvolOption,
                ShiftEvolutionOption=ShiftEvolOption,
                enable_KreissOliger_dissipation=True,
                enable_CAKO=improvements_enable,
                enable_CAHD=improvements_enable,
                enable_SSL=improvements_enable,
                validate_expressions=True,
            )
            ve.compare_or_generate_trusted_results(
                os.path.abspath(__file__),
                os.getcwd(),
                # File basename. If this is set to "trusted_module_test1", then
                #   trusted results_dict will be stored in tests/trusted_module_test1.py
                f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolOption}_{ShiftEvolOption}_{Coord}_T4munu{T4munu_enable}_improvements{improvements_enable}",
                cast(Dict[str, Union[mpf, mpc]], results_dict),
            )

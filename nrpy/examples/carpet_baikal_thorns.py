"""
Generates Einstein Toolkit thorns for solving the wave equation on Cartesian AMR grids with Carpet.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Notes from development:
#
# BSSN_to_ADM:
# - I don't understand the purpose of the second scheduling in Baikal
#   but it is necessary for tests to pass.
# MoL_PostStep:
# - scheduling adjusted to be deterministic w.ur.t other functions

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import re
from collections import OrderedDict as ODict
from typing import Union, cast, List
from pathlib import Path
from inspect import currentframe as cfr
from types import FrameType as FT
import sympy
import shutil
import os

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers import simd

import nrpy.infrastructures.ETLegacy.simple_loop as lp
from nrpy.infrastructures.ETLegacy import boundary_conditions
from nrpy.infrastructures.ETLegacy import CodeParameters
from nrpy.infrastructures.ETLegacy import make_code_defn
from nrpy.infrastructures.ETLegacy import MoL_registration
from nrpy.infrastructures.ETLegacy import Symmetry_registration
from nrpy.infrastructures.ETLegacy import zero_rhss
from nrpy.infrastructures.ETLegacy import schedule_ccl
from nrpy.infrastructures.ETLegacy import interface_ccl
from nrpy.infrastructures.ETLegacy import param_ccl

import nrpy.equations.general_relativity.g4munu_conversions as g4conv
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs
from nrpy.equations.general_relativity.BSSN_gauge_RHSs import BSSN_gauge_RHSs
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support

par.set_parval_from_str("Infrastructure", "ETLegacy")

# Code-generation-time parameters:
project_name = "et_baikal"
thorn_names = ["BaikalVacuum", "Baikal"]
enable_rfm_precompute = False
MoL_method = "RK4"
enable_simd = True
enable_KreissOliger_dissipation = True
parallel_codegen_enable = True
CoordSystem = "Cartesian"
coord_name = ["x", "y", "z"]
coord_name_4D = ["t", "x", "y", "z"]
OMP_collapse = 1

register_MU_gridfunctions = True

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("register_MU_gridfunctions", register_MU_gridfunctions)

standard_ET_includes = ["math.h", "cctk.h", "cctk_Arguments.h", "cctk_Parameters.h"]

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

##############################################
#  Declare the ADM_to_BSSN_order_* functions #
##############################################

def register_CFunction_ADM_to_BSSN(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param fd_order: Order of finite difference method

    :return: A string representing the full C function.
    """

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    default_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)


    desc = """Converting from ADM to BSSN quantities is required in the Einstein Toolkit,
as initial data are given in terms of ADM quantities, and {thorn_name} evolves the BSSN quantities."""
    name = f"{thorn_name}_ADM_to_BSSN_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""
    lapse_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="alp", use_GF_suffix=False)
    lapse2_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="alpha")
    loop_body = lapse2_gf_access + " = " + lapse_gf_access + ";\n"

    for i in range(3):
        shift_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="beta"+coord_name[i], use_GF_suffix=False)
        loop_body += f"const CCTK_REAL local_betaU{i} = {shift_gf_access};\n"

    for i in range(3):
        dtshift_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="dtbeta"+coord_name[i], use_GF_suffix=False)
        loop_body += f"const CCTK_REAL local_BU{i} = {dtshift_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            gamma_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="g"+coord_name[i]+coord_name[j], use_GF_suffix=False)
            loop_body += f"const CCTK_REAL local_gDD{i}{j} = {gamma_gf_access};\n"

    for i in range(3):
        for j in range(i, 3):
            curv_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="k"+coord_name[i]+coord_name[j], use_GF_suffix=False)
            loop_body += f"const CCTK_REAL local_kDD{i}{j} = {curv_gf_access};\n"

    gammaCartDD = ixp.declarerank2("local_gDD", symmetry="sym01")
    KCartDD = ixp.declarerank2("local_kDD", symmetry="sym01")

    betaU = ixp.declarerank1("local_betaU")
    BU = ixp.declarerank1("local_BU")
    adm2bssn = ADM_to_BSSN(gammaCartDD, KCartDD, betaU, BU, "Cartesian")

    cf_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="cf")
    trK_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="trK")
    list_of_output_exprs = [adm2bssn.cf, adm2bssn.trK]
    list_of_output_varnames = [cf_gf_access, cf_gf_access]
    for i in range(3):
        vetU_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"vetU{i}")
        list_of_output_exprs += [adm2bssn.vetU[i]]
        list_of_output_varnames += [vetU_gf_access]
    for i in range(3):
        betU_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"betU{i}")
        list_of_output_exprs += [adm2bssn.betU[i]]
        list_of_output_varnames += [betU_gf_access]
    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"hDD{i}{j}")
            list_of_output_exprs += [adm2bssn.hDD[i][j]]
            list_of_output_varnames += [hDD_gf_access]
    for i in range(3):
        for j in range(i, 3):
            aDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"aDD{i}{j}")
            list_of_output_exprs += [adm2bssn.aDD[i][j]]
            list_of_output_varnames += [aDD_gf_access]

    loop_body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )
    loop_body = loop_body.rstrip()

    body += lp.simple_loop(
        loop_body=loop_body,
        enable_simd=False,
        loop_region="all points",
        enable_OpenMP=True,
    )

    body += "\n"

    Bq = BSSN_quantities[CoordSystem]
    gammabarUU = Bq.gammabarUU
    GammabarUDD = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    # (from the reference metric):
    rfm = refmetric.reference_metric[CoordSystem]
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (
                    GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                )

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    body += lp.simple_loop(
        ccg.c_codegen(
            lambdaU,
            [
                gri.ETLegacyGridFunction.access_gf("lambdaU0"),
                gri.ETLegacyGridFunction.access_gf("lambdaU1"),
                gri.ETLegacyGridFunction.access_gf("lambdaU2"),
            ],
            verbose=False,
            include_braces=False,
            enable_fd_codegen=True,
        ),
        loop_region="interior",
    )

    body += f"""
  ExtrapolateGammas(cctkGH, lambdaU0GF);
  ExtrapolateGammas(cctkGH, lambdaU1GF);
  ExtrapolateGammas(cctkGH, lambdaU2GF);
"""

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME at CCTK_INITIAL after ADMBase_PostInitial
  {{
    LANG: C
    READS:  ADMBase::metric,
            ADMBase::shift,
            ADMBase::curv,
            ADMBase::dtshift,
            ADMBase::lapse
    WRITES: evol_variables
    SYNC: evol_variables
  }} "Convert initial data into BSSN variables"
}}
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=standard_ET_includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("CCTK_INITIAL", schedule)],
    )
    par.set_parval_from_str("fd_order", default_fd_order)
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

#############################################
#  Declare the Ricci_eval_order_* functions #
#############################################

def register_CFunction_Ricci_eval(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int,
    OMP_collapse: int = 1,
    ) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param fd_order: Order of finite difference method
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    default_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    par.set_parval_from_str("enable_RbarDD_gridfunctions", False)

    includes = standard_ET_includes
    if enable_simd:
        includes += [("./SIMD/simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    name = f"{thorn_name}_Ricci_eval_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
"""
    if enable_simd:
        body +=f"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));

"""
    else:
        body +=f"""  DECLARE_CCTK_PARAMETERS;

"""

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Populate Ricci tensor
    Ricci_access_gfs: List[str] = []
    for var in Bq.Ricci_varnames:
        Ricci_access_gfs += [
            gri.ETLegacyGridFunction.access_gf(gf_name=var)
        ]
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            Bq.Ricci_exprs,
            Ricci_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in MoL_CalcRHS as {thorn_name}_Ricci before {thorn_name}_RHS
  {{
    LANG: C
    READS:  hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
            lambdaU0GF, lambdaU1GF, lambdaU2GF
    WRITES: RbarDD00GF, RbarDD01GF, RbarDD02GF, RbarDD11GF, RbarDD12GF, RbarDD22GF
  }} "Compute Ricci tensor, needed for BSSN RHSs, at finite-differencing order {fd_order}"
}}
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_CalcRHS", schedule)],
    )
    par.set_parval_from_str("fd_order", default_fd_order)
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

###########################################
#  Declare the rhs_eval_order_* functions #
###########################################

def register_CFunction_rhs_eval(
    thorn_name: str,
    CoordSystem: str,
    enable_T4munu: bool,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int,
    enable_KreissOliger_dissipation: bool,
    KreissOliger_strength_mult_by_W: bool = False,
    # when mult by W, strength_gauge=0.99 & strength_nongauge=0.3 is best.
    KreissOliger_strength_gauge: float = 0.3,
    KreissOliger_strength_nongauge: float = 0.3,
    OMP_collapse: int = 1,
    ) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the BSSN equations.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_T4munu: Whether to include the stress-energy tensor. Defaults to False.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD (Single Instruction, Multiple Data).
    :param fd_order: Order of finite difference method
    :param enable_KreissOliger_dissipation: Whether or not to enable Kreiss-Oliger dissipation.
    :param KreissOliger_strength_mult_by_W: Whether to multiply Kreiss-Oliger strength by W.
    :param KreissOliger_strength_gauge: Gauge strength for Kreiss-Oliger dissipation.
    :param KreissOliger_strength_nongauge: Non-gauge strength for Kreiss-Oliger dissipation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    default_fd_order = par.parval_from_str("fd_order")
    default_enable_T4munu = par.parval_from_str("enable_T4munu")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    par.set_parval_from_str("enable_T4munu", enable_T4munu)
    par.set_parval_from_str("enable_RbarDD_gridfunctions", True)

    includes = standard_ET_includes
    if enable_simd:
        includes += [("./SIMD/simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    name = f"{thorn_name}_rhs_eval_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
"""
    if enable_simd:
        body +=f"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));
  const CCTK_REAL *param_eta CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("eta", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY eta CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_eta);
  const CCTK_REAL *param_diss_strength CCTK_ATTRIBUTE_UNUSED = CCTK_ParameterGet("diss_strength", "{thorn_name}", NULL);
  const REAL_SIMD_ARRAY diss_strength CCTK_ATTRIBUTE_UNUSED = ConstSIMD(*param_diss_strength);

"""
    else:
        body +=f"""  DECLARE_CCTK_PARAMETERS;

"""

    rhs = BSSN_RHSs[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    alpha_rhs, vet_rhsU, bet_rhsU = BSSN_gauge_RHSs(
        CoordSystem,
        enable_rfm_precompute,
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
                CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else CoordSystem
            ]
            EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
            if EvolvedConformalFactor_cf == "W":
                diss_strength_gauge *= Bq.cf
                diss_strength_nongauge *= Bq.cf
            elif EvolvedConformalFactor_cf == "chi":
                diss_strength_gauge *= sympy.sqrt(Bq.cf)
                diss_strength_nongauge *= sympy.sqrt(Bq.cf)
            elif EvolvedConformalFactor_cf == "phi":
                diss_strength_gauge *= sympy.exp(-2 * Bq.cf)
                diss_strength_nongauge *= sympy.exp(-2 * Bq.cf)

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
                #if "2ndOrder" in ShiftEvolutionOption:
                #    rhs.BSSN_RHSs_varname_to_expr_dict[f"bet_rhsU{i}"] += (
                #        diss_strength_gauge * betU_dKOD[i][k] * rfm.ReU[k]
                #    )  # ReU[k] = 1/scalefactor_orthog_funcform[k]
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
    for var in rhs.BSSN_RHSs_varname_to_expr_dict.keys():
        pattern = re.compile(r'([a-zA-Z_]+)_rhs([a-zA-Z_]+[0-2]+)')
        patmatch = pattern.match(var)
        if patmatch:
            var_name = patmatch.group(1) + patmatch.group(2) + "_rhs"
        else:
            var_name = var
        BSSN_RHSs_access_gf += [
            gri.ETLegacyGridFunction.access_gf(gf_name=var)
        ]

    # Set up upwind control vector (betaU)
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    betaU = ixp.zerorank1()
    vetU = ixp.declarerank1("vetU")
    for i in range(3):
        betaU[i] = vetU[i] * rfm.ReU[i]
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list(rhs.BSSN_RHSs_varname_to_expr_dict.values()),
            BSSN_RHSs_access_gf,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            upwind_control_vec=betaU,
            enable_fd_functions=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in MoL_CalcRHS as {thorn_name}_RHS after {thorn_name}_Ricci
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
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_CalcRHS", schedule)],
        ET_current_thorn_CodeParams_used=params
    )
    par.set_parval_from_str("fd_order", default_fd_order)
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

######################################
#  Declare the BSSN_to_ADM function  #
######################################

def register_CFunction_BSSN_to_ADM(
    thorn_name: str,
    CoordSystem: str,
    OMP_collapse: int = 1,
    ) -> Union[None, pcg.NRPyEnv_type]:
    """
    Convert BSSN variables in the Cartesian basis to ADM variables in the Cartesian basis.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = standard_ET_includes
    if enable_simd:
        includes += [("./SIMD/simd_intrinsics.h")]
    desc = r"""Set RHSs for the BSSN evolution equations."""
    name = f"{thorn_name}_BSSN_to_ADM"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    bssn2adm = BSSN_to_ADM(
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    )

    # lapse, shift, and dtshift are just straight copies
    lapse_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="alp", use_GF_suffix=False)
    lapse2_gf_access  = gri.ETLegacyGridFunction.access_gf(gf_name="alpha")

    loop_body = lapse_gf_access + " = " + lapse2_gf_access + ";\n"

    for i in range(3):
        bssn_shift_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="vetU"+str(i))
        shift_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="beta"+coord_name[i], use_GF_suffix=False)
        loop_body += f"{shift_gf_access} = {bssn_shift_gf_access};\n"

    for i in range(3):
        bssn_dtshift_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="bet"+str(i), use_GF_suffix=False)
        dtshift_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="dtbeta"+coord_name[i], use_GF_suffix=False)
        loop_body += f"{dtshift_gf_access} = {bssn_dtshift_gf_access};\n"

    list_of_output_exprs = []
    list_of_output_varnames = []

    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="hDD"+str(i)+str(j))
            loop_body += f"const CCTK_REAL hDD{i}{j} = {hDD_gf_access};\n"

            gamma_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="g"+coord_name[i]+coord_name[j], use_GF_suffix=False)
            list_of_output_exprs += [bssn2adm.gammaDD[i][j]]
            list_of_output_varnames += [gamma_gf_access]

    for i in range(3):
        for j in range(i, 3):
            aDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="k"+str(i)+str(j))
            loop_body += f"const CCTK_REAL aDD{i}{j} = {aDD_gf_access};\n"

            curv_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="k"+coord_name[i]+coord_name[j], use_GF_suffix=False)
            list_of_output_exprs += [bssn2adm.KDD[i][j]]
            list_of_output_varnames += [curv_gf_access]

    loop_body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
        enable_simd=False,
    )
    loop_body = loop_body.rstrip()

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule_poststep = ("MoL_PostStep",
        f"""
schedule FUNC_NAME in MoL_PostStep after {thorn_name}_enforce_detgammahat_constraint before ADMBase_SetADMVars
{{
  LANG: C
  READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          vetU0GF, vetU1GF, vetU2GF, betU0GF, betU1GF, betU2GF,
          cfGF, trKGF, alphaGF
  WRITES: ADMBase::metric(everywhere),
          ADMBase::shift(everywhere),
          ADMBase::curv(everywhere),
          ADMBase::dtshift(everywhere),
          ADMBase::lapse(everywhere)
}} "Perform BSSN-to-ADM conversion. Useful for diagnostics."
""")

    schedule_pseudoevol = ("MoL_PseudoEvolution",
        f"""
schedule FUNC_NAME in MoL_PseudoEvolution after {thorn_name}_auxgfs_ApplyBCs
{{
  LANG: C
  READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          vetU0GF, vetU1GF, vetU2GF, betU0GF, betU1GF, betU2GF,
          cfGF, trKGF, alphaGF
  WRITES: ADMBase::metric(everywhere),
          ADMBase::shift(everywhere),
          ADMBase::curv(everywhere),
          ADMBase::dtshift(everywhere),
          ADMBase::lapse(everywhere)
}} "Perform BSSN-to-ADM conversion in MoL_PseudoEvolution. Needed for proper HydroBase integration."
""")

    ET_schedule_bins_entries = [schedule_poststep]
    if thorn_name == "Baikal":
        ET_schedule_bins_entries += [schedule_pseudoevol]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

########################################################
#  Declare the enforce_detgammahat_constraint function #
########################################################
def register_CFunction_enforce_detgammahat_constraint(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    desc = r"""Enforce det(gammabar) = det(gammahat) constraint. Required for strong hyperbolicity."""
    name = f"{thorn_name}_enforce_detgammahat_constraint"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    # First define the Kronecker delta:
    KroneckerDeltaDD = ixp.zerorank2()
    for i in range(3):
        KroneckerDeltaDD[i][i] = sympy.sympify(1)

    # The detgammabar in BSSN_RHSs is set to detgammahat when BSSN_RHSs::detgbarOverdetghat_equals_one=True (default),
    #    so we manually compute it here:
    dummygammabarUU, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)

    # Next apply the constraint enforcement equation above.
    nrpyAbs = sympy.Function("nrpyAbs")
    hprimeDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            hprimeDD[i][j] = (nrpyAbs(rfm.detgammahat) / detgammabar) ** (
                sympy.Rational(1, 3)
            ) * (KroneckerDeltaDD[i][j] + Bq.hDD[i][j]) - KroneckerDeltaDD[i][j]

    hDD_access_gfs: List[str] = []
    hprimeDD_expr_list: List[sympy.Expr] = []
    for i in range(3):
        for j in range(i, 3):
            hDD_access_gfs += [
                gri.ETLegacyGridFunction.access_gf(gf_name=f"hDD{i}{j}")
            ]
            hprimeDD_expr_list += [hprimeDD[i][j]]

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            hprimeDD_expr_list,
            hDD_access_gfs,
            enable_simd=False,
            automatically_read_gf_data_from_memory=True
        ),
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
schedule FUNC_NAME in MoL_PostStep
{{
  LANG: C
  READS:  hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF
  WRITES: hDD00GF(everywhere), hDD01GF(everywhere), hDD02GF(everywhere),
          hDD11GF(everywhere), hDD12GF(everywhere), hDD22GF(everywhere)
}} "Enforce detgammabar = detgammahat (= 1 in Cartesian)"
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=standard_ET_includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_PostStep", schedule)],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

########################################################
#  Declare the T4DD_to_T4UU function #
########################################################
def register_CFunction_T4DD_to_T4UU(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    desc = r"""Compute T4UU from T4DD (provided by TmunuBase),
using BSSN quantities as inputs for the 4D raising operation

WARNING: Do not enable SIMD here, as it is not guaranteed that
         cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2] is a multiple of
         SIMD_width!"""
    name = f"{thorn_name}_T4DD_to_T4UU"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    g4UU = g4conv.BSSN_to_g4UU(CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute)

    T4DD = ixp.declarerank2("T4DD", dimension=4, symmetry="sym01")
    T4UU = ixp.zerorank2(dimension=4)
    for mu in range(4):
        for nu in range(4):
            T4UU[mu][nu] = sum([
                               sum([T4DD[alpha][beta] * g4UU[mu][alpha] * g4UU[nu][beta]
                                    for alpha in range(4)])
                                for beta in range(4)])

    T4DD_access_gfs: List[str] = []
    T4UU_expr_list: List[sympy.Expr] = []
    for i in range(4):
        for j in range(i, 4):
            T4DD_access_gfs += [
                gri.ETLegacyGridFunction.access_gf(
                    gf_name=f"T4UU{i}{j}")
            ]
            T4UU_expr_list += [T4UU[i][j]]

    loop_body = ""
    lapse_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"alpha")
    loop_body += f"const CCTK_REAL alpha = {lapse_gf_access};\n"
    cf_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"cf")
    loop_body += f"const CCTK_REAL cf = {cf_gf_access};\n"
    for i in range(3):
        vet_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"vetU{i}")
        loop_body += f"const CCTK_REAL vetU{i} = {vet_gf_access};\n"
    for i in range(3):
        for j in range(i, 3):
            hDD_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name=f"hDD{i}{j}")
            loop_body += f"const CCTK_REAL hDD{i}{j} = {hDD_gf_access};\n"
    for i in range(4):
        for j in range(i, 4):
            Tmunu_gf_access = gri.ETLegacyGridFunction.access_gf(gf_name="eT"+coord_name_4D[i]+coord_name_4D[j], use_GF_suffix=False)
            loop_body += f"const CCTK_REAL T4DD{i}{j} = {Tmunu_gf_access};\n"

    loop_body += ccg.c_codegen(
        T4UU_expr_list,
        T4DD_access_gfs,
        enable_simd=False,
    )

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule1 = f"""
schedule FUNC_NAME in MoL_CalcRHS before {thorn_name}_RHS
{{
  LANG: C
  READS:  TmunuBase::stress_energy_scalar,
          TmunuBase::stress_energy_vector,
          TmunuBase::stress_energy_tensor,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          alphaGF, cfGF, vetU0GF, vetU1GF, vetU2GF
  WRITES: T4UU00GF(everywhere), T4UU01GF(everywhere), T4UU02GF(everywhere), T4UU03GF(everywhere),
          T4UU11GF(everywhere), T4UU12GF(everywhere), T4UU13GF(everywhere),
          T4UU22GF(everywhere), T4UU23GF(everywhere), T4UU33GF(everywhere)
}} "Compute T4UU from T4DD (provided in eT?? from TmunuBase), needed for BSSN RHSs"
"""
    schedule2 = f"""
schedule FUNC_NAME in MoL_PseudoEvolution before {thorn_name}_BSSN_constraints
{{
  LANG: C
  READS:  TmunuBase::stress_energy_scalar,
          TmunuBase::stress_energy_vector,
          TmunuBase::stress_energy_tensor,
          hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          alphaGF, cfGF, vetU0GF, vetU1GF, vetU2GF
  WRITES: T4UU00GF(everywhere), T4UU01GF(everywhere), T4UU02GF(everywhere), T4UU03GF(everywhere),
          T4UU11GF(everywhere), T4UU12GF(everywhere), T4UU13GF(everywhere),
          T4UU22GF(everywhere), T4UU23GF(everywhere), T4UU33GF(everywhere)
}} "Compute T4UU from T4DD (provided in eT?? from TmunuBase), needed for BSSN constraints"
"""
    schedule_RHS = ("MoL_CalcRHS", schedule1)
    schedule_constraints = ("MoL_PseudoEvolution", schedule2)

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=standard_ET_includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[schedule_RHS, schedule_constraints],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

#########################################
#  Declare the floor_the_lapse function #
#########################################
def register_CFunction_floor_the_lapse(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the function that enforces the det(gammabar) = det(gammahat) constraint.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """

    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = """Apply floor to the lapse."""
    name = f"{thorn_name}_floor_the_lapse"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

#ifndef MAX
#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#endif

"""
    lapse_access_gfs = [
        gri.ETLegacyGridFunction.access_gf(gf_name=f"alpha")
            ]
    loop_body = lapse_access_gfs[0] + " = MAX(" + lapse_access_gfs[0] + ", lapse_floor);"

    body += lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
schedule FUNC_NAME in MoL_PostStep before {thorn_name}_enforce_detgammahat_constraint
{{
  LANG: C
  READS:  alphaGF(everywhere)
  WRITES: alphaGF(everywhere)
}} "Set lapse = max(lapse_floor, lapse)"
"""

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=standard_ET_includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_PostStep", schedule)],
        ET_current_thorn_CodeParams_used=["lapse_floor"]
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

###################################################
#  Declare the BSSN_constraints_order_* functions #
###################################################

def register_CFunction_constraints(
    thorn_name: str,
    enable_T4munu: bool,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    fd_order: int,
    OMP_collapse: int = 1,
    ) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_T4munu: Whether to include the stress-energy tensor.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD instructions.
    :param fd_order: Order of finite difference method
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    default_fd_order = par.parval_from_str("fd_order")
    default_enable_T4munu = par.parval_from_str("enable_T4munu")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    par.set_parval_from_str("enable_T4munu", enable_T4munu)

    includes = standard_ET_includes
    if enable_simd:
        includes += [("./SIMD/simd_intrinsics.h")]
    desc = r"""Evaluate BSSN constraints."""
    name = f"{thorn_name}_BSSN_constraints_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
"""
    if enable_simd:
        body +=f"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));

"""
    else:
        body +=f"""  DECLARE_CCTK_PARAMETERS;

"""

    Bcon = BSSN_constraints[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    list_of_output_exprs = [Bcon.H]
    Constraints_access_gfs: List[str] = [
        gri.ETLegacyGridFunction.access_gf(gf_name="H")]

    for index in range(3):
        list_of_output_exprs += [Bcon.MU[index]]
        Constraints_access_gfs += [
            gri.ETLegacyGridFunction.access_gf(gf_name="MU"+str(index))
        ]
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list_of_output_exprs,
            Constraints_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
if(FD_order == {fd_order}) {{
  schedule FUNC_NAME in MoL_PseudoEvolution as {thorn_name}_BSSN_constraints after BaikalVacuum_ApplyBCs
  {{
    LANG: C
    READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
            hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
            trKGF, cfGF, lambdaU0GF, lambdaU1GF, lambdaU2GF"""
    if enable_T4munu:
        schedule += f""",
            alphaGF, vetU0GF, vetU1GF, vetU2GF,
            T4UU00GF, T4UU01GF, T4UU02GF, T4UU03GF"""
    schedule += f"""
    WRITES: aux_variables
  }} "Compute BSSN (Hamiltonian and momentum) constraints, at finite-differencing order {fd_order}"
}}
"""

    params = None
    if thorn_name == "Baikal":
        params = ["PI"]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_PseudoEvolution", schedule)],
        ET_current_thorn_CodeParams_used=params
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

#########################################
#  Declare the RegisterSlicing function #
#########################################

def register_CFunction_slicing(
    thorn_name: str,
    ) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the slicing registration function.

    :param thorn_name: The Einstein Toolkit thorn name.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = f"""Register slicing condition for NRPy+-generated thorn {thorn_name}."""
    name = f"{thorn_name}_RegisterSlicing"
    body = f"""
    Einstein_RegisterSlicing ("{thorn_name}");
  return 0;"""

    schedule = """
schedule FUNC_NAME at STARTUP
{
  LANG: C
  OPTIONS: meta
} "Register 3+1 slicing condition"
"""


    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=["Slicing.h", "cctk.h"],
        desc=desc,
        c_type="int",
        name=name,
        params="",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("STARTUP", schedule)],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

###################################
#  Generate functions in parallel #
###################################

for evol_thorn_name in thorn_names:
    lapse_floor = par.register_CodeParameter(
        "CCTK_REAL",
        __name__,
        "lapse_floor",
        "1e-15",
        add_to_glb_code_params_dict=True,
    )
    fd_order_param = par.register_CodeParameter(
        "CCTK_INT",
        __name__,
        "FD_order",
        "4",
        add_to_glb_code_params_dict=True,
    )

    enable_T4munu = False
    fd_order_list = [4,6,8]
    if evol_thorn_name == "Baikal":
        enable_T4munu = True
        fd_order_list = [2,4]
        register_CFunction_T4DD_to_T4UU(thorn_name=evol_thorn_name, CoordSystem=CoordSystem,
                                        enable_rfm_precompute=False)

    for fd_order in fd_order_list:
        par.set_parval_from_str("fd_order", fd_order)
        register_CFunction_ADM_to_BSSN(thorn_name=evol_thorn_name, CoordSystem=CoordSystem, fd_order=fd_order,
                                       enable_rfm_precompute=enable_rfm_precompute, enable_simd=enable_simd)
        register_CFunction_Ricci_eval(thorn_name=evol_thorn_name, fd_order=fd_order, CoordSystem="Cartesian",
                                      enable_rfm_precompute=False, enable_simd=enable_simd)
        register_CFunction_rhs_eval(thorn_name=evol_thorn_name, enable_T4munu=enable_T4munu, fd_order=fd_order, CoordSystem="Cartesian",
                                    enable_rfm_precompute=False, enable_simd=enable_simd,
                                    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation)
        register_CFunction_constraints(thorn_name=evol_thorn_name, enable_T4munu=enable_T4munu, fd_order=fd_order, CoordSystem="Cartesian",
                                       enable_rfm_precompute=False, enable_simd=enable_simd)

    register_CFunction_slicing(thorn_name=evol_thorn_name)
    register_CFunction_BSSN_to_ADM(thorn_name=evol_thorn_name, CoordSystem="Cartesian")
    register_CFunction_floor_the_lapse(thorn_name=evol_thorn_name, CoordSystem="Cartesian",
                                       enable_rfm_precompute=False)
    register_CFunction_enforce_detgammahat_constraint(thorn_name=evol_thorn_name, CoordSystem="Cartesian",
                                                      enable_rfm_precompute=False)

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

for evol_thorn_name in thorn_names:
    ########################
    # STEP 2: Register functions that depend on all gridfunctions & CodeParameters having been set:
   
    Symmetry_registration.register_CFunction_Symmetry_registration_oldCartGrid3D(
        thorn_name=evol_thorn_name
    )
    boundary_conditions.register_CFunctions(thorn_name=evol_thorn_name)
    zero_rhss.register_CFunction_zero_rhss(thorn_name=evol_thorn_name)
    MoL_registration.register_CFunction_MoL_registration(thorn_name=evol_thorn_name)
    
    ########################
    # STEP 3: All functions have been registered at this point. Time to output the CCL files & thorns!
    
    CParams_registered_to_params_ccl: List[str] = []
    
    # CCL files: evol_thorn
    schedule_ccl.construct_schedule_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
        STORAGE="""
STORAGE: evol_variables[3]     # Evolution variables
STORAGE: evol_variables_rhs[1] # Variables storing right-hand-sides
STORAGE: auxevol_variables[1]  # Single-timelevel storage of variables needed for evolutions.
STORAGE: aux_variables[3]      # Diagnostics variables""",
    )
    interface_ccl.construct_interface_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
        inherits="ADMBase Boundary Grid TmunuBase",
        USES_INCLUDEs="""USES INCLUDE: Symmetry.h
USES INCLUDE: Boundary.h
USES INCLUDE: Slicing.h

# Needed to convert ADM initial data into BSSN initial data (gamma extrapolation)
CCTK_INT FUNCTION ExtrapolateGammas(CCTK_POINTER_TO_CONST IN cctkGH, CCTK_REAL ARRAY INOUT var)
REQUIRES FUNCTION ExtrapolateGammas""",
        is_evol_thorn=True,
        enable_NewRad=True,
    )
    
    params_str = f"""
shares: ADMBase

EXTENDS CCTK_KEYWORD evolution_method "evolution_method"
{{
  "{evol_thorn_name}" :: ""
}}

EXTENDS CCTK_KEYWORD lapse_evolution_method "lapse_evolution_method"
{{
  "{evol_thorn_name}" :: ""
}}

EXTENDS CCTK_KEYWORD shift_evolution_method "shift_evolution_method"
{{
  "{evol_thorn_name}" :: ""
}}

EXTENDS CCTK_KEYWORD dtshift_evolution_method "dtshift_evolution_method"
{{
  "{evol_thorn_name}" :: ""
}}

EXTENDS CCTK_KEYWORD dtlapse_evolution_method "dtlapse_evolution_method"
{{
  "{evol_thorn_name}" :: ""
}}"""
    CParams_registered_to_params_ccl += param_ccl.construct_param_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
        shares_extends_str=params_str,
    )
    
    make_code_defn.output_CFunctions_and_construct_make_code_defn(
        project_dir=project_dir, thorn_name=evol_thorn_name
    )

    simd.copy_simd_intrinsics_h(
        project_dir=str(Path(project_dir) / evol_thorn_name / "src")
    )

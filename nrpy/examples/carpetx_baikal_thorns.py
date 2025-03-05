"""
Generates Einstein Toolkit thorns for evolving the BSSN equations on Cartesian AMR grids with CarpetX.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

import os
from pathlib import Path

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
from typing import List

import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.CarpetX import (
    boundary_conditions,
    configuration_ccl,
    general_relativity,
    interface_ccl,
    make_code_defn,
    param_ccl,
    schedule_ccl,
    zero_rhss,
)

# Code-generation-time parameters:
project_name = "et_baikalx"
thorn_names = [
    "BaikalX",
    "BaikalVacuumX",
]  # note that this ordering matters due to how ccl files generate
enable_rfm_precompute = False
enable_simd = False
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_NoCovariant"
enable_KreissOliger_dissipation = True
parallel_codegen_enable = True
CoordSystem = "Cartesian"
register_MU_gridfunctions = True

project_dir = os.path.join("project", project_name)

par.set_parval_from_str("Infrastructure", "CarpetX")
par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("register_MU_gridfunctions", register_MU_gridfunctions)

########################
# STEP 1: register functions and parameters
########################

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
    fd_order_list = [4, 6, 8]
    if evol_thorn_name == "BaikalX":
        enable_T4munu = True
        fd_order_list = [2, 4]
        general_relativity.T4DD_to_T4UU.register_CFunction_T4DD_to_T4UU(
            thorn_name=evol_thorn_name,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=False,
        )

    for fd_order in fd_order_list:
        par.set_parval_from_str("fd_order", fd_order)
        general_relativity.ADM_to_BSSN.register_CFunction_ADM_to_BSSN(
            thorn_name=evol_thorn_name,
            CoordSystem=CoordSystem,
            fd_order=fd_order,
        )
        general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
            thorn_name=evol_thorn_name,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_simd=enable_simd,
            fd_order=fd_order,
        )
        general_relativity.rhs_eval.register_CFunction_rhs_eval(
            thorn_name=evol_thorn_name,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_T4munu=enable_T4munu,
            enable_simd=enable_simd,
            fd_order=fd_order,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
        general_relativity.BSSN_constraints.register_CFunction_BSSN_constraints(
            thorn_name=evol_thorn_name,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_T4munu=enable_T4munu,
            enable_simd=enable_simd,
            fd_order=fd_order,
        )

    general_relativity.BSSN_to_ADM.register_CFunction_BSSN_to_ADM(
        thorn_name=evol_thorn_name, CoordSystem="Cartesian"
    )
    general_relativity.floor_the_lapse.register_CFunction_floor_the_lapse(
        thorn_name=evol_thorn_name
    )
    general_relativity.enforce_detgammahat_constraint.register_CFunction_enforce_detgammahat_constraint(
        thorn_name=evol_thorn_name, CoordSystem="Cartesian", enable_rfm_precompute=False
    )

########################
#  STEP 2: Generate functions in parallel
########################

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

for evol_thorn_name in thorn_names:
    if evol_thorn_name == "BaikalVacuumX":
        for i in range(4):
            for j in range(i, 4):
                gri.glb_gridfcs_dict.pop("T4UU" + str(i) + str(j))

    ########################
    # STEP 3: Register functions that depend on all gridfunctions & CodeParameters having been set
    ########################

    boundary_conditions.register_CFunctions(thorn_name=evol_thorn_name)
    zero_rhss.register_CFunction_zero_rhss(thorn_name=evol_thorn_name)

    ########################
    # STEP 4: Generate ccl files for this thorn
    ########################

    CParams_registered_to_params_ccl: List[str] = []

    # CCL files
    # only 1 timelevel of storage (no subcycling)
    schedule_ccl.construct_schedule_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
        STORAGE="""
STORAGE: evol_variables[1]     # Evolution variables
STORAGE: evol_variables_rhs[1] # Variables storing right-hand-sides
STORAGE: auxevol_variables[1]  # Single-timelevel storage of variables needed for evolutions.
STORAGE: aux_variables[1]      # Diagnostics variables""",
    )

    configuration_ccl.construct_configuration_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
    )

    inherits = "ADMBaseX"
    if evol_thorn_name == "Baikal":
        inherits += " TmunuBaseX"
    interface_ccl.construct_interface_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
        inherits=inherits,
        USES_INCLUDEs="USES INCLUDE: loop_device.hxx",
        is_evol_thorn=True,
        enable_NewRad=True,
    )

    CParams_registered_to_params_ccl += param_ccl.construct_param_ccl(
        project_dir=project_dir,
        thorn_name=evol_thorn_name,
    )

    make_code_defn.output_CFunctions_and_construct_make_code_defn(
        project_dir=project_dir, thorn_name=evol_thorn_name
    )

    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=str(Path(project_dir) / evol_thorn_name / "src"),
        subdirectory="simd",
    )

    simd_name = str(
        Path(project_dir) / evol_thorn_name / "src" / "simd" / "simd_intrinsics.h"
    )
    with open(simd_name, "r", encoding="utf-8") as file:
        contents = file.read()

    new_contents = contents.replace("REAL_SIMD_ARRAY REAL", "REAL_SIMD_ARRAY CCTK_REAL")

    with open(simd_name, "w", encoding="utf-8") as file:
        file.write(new_contents)

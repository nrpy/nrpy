"""
Generates Einstein Toolkit thorns for evolving the BSSN equations on Cartesian AMR grids with CarpetX.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
from typing import List
from pathlib import Path
import os

import nrpy.grid as gri
import nrpy.params as par
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers import simd

from nrpy.infrastructures.CarpetX import boundary_conditions
from nrpy.infrastructures.CarpetX import make_code_defn
from nrpy.infrastructures.CarpetX import zero_rhss
from nrpy.infrastructures.CarpetX import schedule_ccl
from nrpy.infrastructures.CarpetX import interface_ccl
from nrpy.infrastructures.CarpetX import param_ccl
from nrpy.infrastructures.CarpetX import configuration_ccl

# All needed functions can be imported from the CarpetX infrastructure
from nrpy.infrastructures.CarpetX.general_relativity.ADM_to_BSSN import (
    register_CFunction_ADM_to_BSSN,
)
from nrpy.infrastructures.CarpetX.general_relativity.BSSN_to_ADM import (
    register_CFunction_BSSN_to_ADM,
)
from nrpy.infrastructures.CarpetX.general_relativity.enforce_detgammahat_constraint import (
    register_CFunction_enforce_detgammahat_constraint,
)
from nrpy.infrastructures.CarpetX.general_relativity.floor_the_lapse import (
    register_CFunction_floor_the_lapse,
)
from nrpy.infrastructures.CarpetX.general_relativity.BSSN_constraints import (
    register_CFunction_BSSN_constraints,
)
from nrpy.infrastructures.CarpetX.general_relativity.Ricci_eval import (
    register_CFunction_Ricci_eval,
)
from nrpy.infrastructures.CarpetX.general_relativity.rhs_eval import (
    register_CFunction_rhs_eval,
)
from nrpy.infrastructures.CarpetX.general_relativity.T4DD_to_T4UU import (
    register_CFunction_T4DD_to_T4UU,
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
        register_CFunction_T4DD_to_T4UU(
            thorn_name=evol_thorn_name,
            CoordSystem=CoordSystem,
            enable_rfm_precompute=False,
        )

    for fd_order in fd_order_list:
        par.set_parval_from_str("fd_order", fd_order)
        register_CFunction_ADM_to_BSSN(
            thorn_name=evol_thorn_name,
            CoordSystem=CoordSystem,
            fd_order=fd_order,
        )
        register_CFunction_Ricci_eval(
            thorn_name=evol_thorn_name,
            fd_order=fd_order,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_simd=enable_simd,
        )
        register_CFunction_rhs_eval(
            thorn_name=evol_thorn_name,
            enable_T4munu=enable_T4munu,
            fd_order=fd_order,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_simd=enable_simd,
            LapseEvolutionOption=LapseEvolutionOption,
            ShiftEvolutionOption=ShiftEvolutionOption,
            enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        )
        register_CFunction_BSSN_constraints(
            thorn_name=evol_thorn_name,
            enable_T4munu=enable_T4munu,
            fd_order=fd_order,
            CoordSystem="Cartesian",
            enable_rfm_precompute=False,
            enable_simd=enable_simd,
        )

    register_CFunction_BSSN_to_ADM(thorn_name=evol_thorn_name, CoordSystem="Cartesian")
    register_CFunction_floor_the_lapse(thorn_name=evol_thorn_name)
    register_CFunction_enforce_detgammahat_constraint(
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

    simd.copy_simd_intrinsics_h(
        project_dir=str(Path(project_dir) / evol_thorn_name / "src")
    )

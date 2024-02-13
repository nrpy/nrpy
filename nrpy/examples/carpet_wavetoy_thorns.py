"""
Generates Einstein Toolkit thorns for solving the wave equation on Cartesian AMR grids with Carpet.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
from typing import Union, cast, List
from pathlib import Path
from inspect import currentframe as cfr
from types import FrameType as FT
import os

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.helpers.parallel_codegen as pcg
from nrpy.helpers import simd

from nrpy.equations.wave_equation.WaveEquation_Solutions_InitialData import (
    WaveEquation_solution_Cartesian,
)
from nrpy.equations.wave_equation.WaveEquation_RHSs import WaveEquation_RHSs
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


par.set_parval_from_str("Infrastructure", "ETLegacy")

# Code-generation-time parameters:
project_name = "et_wavetoy"
ID_thorn_name = "IDWaveToyNRPy"
diag_thorn_name = "diagWaveToyNRPy"
evol_thorn_name = "WaveToyNRPy"
WaveType = "SphericalGaussian"
default_sigma = 3.0
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size
default_diagnostics_output_every = 0.5
default_checkpoint_every = 50.0
enable_rfm_precompute = False
MoL_method = "RK4"
fd_order = 8
enable_simd = True
enable_KreissOliger_dissipation = False
parallel_codegen_enable = True
CoordSystem = "Cartesian"
OMP_collapse = 1

project_dir = os.path.join("project", project_name)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
standard_ET_includes = ["math.h", "cctk.h", "cctk_Arguments.h", "cctk_Parameters.h"]


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_exact_solution_all_points(
    thorn_name: str = "",
    in_WaveType: str = "SphericalGaussian",
    in_default_sigma: float = 3.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a C function for the exact solution at a single point.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param in_WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param in_default_sigma: The default value for the Gaussian width (sigma).
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = standard_ET_includes

    # Populate uu_exactsoln, vv_exactsoln
    exactsoln = WaveEquation_solution_Cartesian(
        WaveType=in_WaveType,
        default_sigma=in_default_sigma,
    )

    prefunc = r"""
// Exact solution at a single point.
static void WaveToy_exact_solution_single_point(const CCTK_REAL time, const CCTK_REAL xCart0, const CCTK_REAL xCart1, const CCTK_REAL xCart2,
    CCTK_REAL *restrict exact_soln_UUGF, CCTK_REAL *restrict exact_soln_VVGF) {
DECLARE_CCTK_PARAMETERS;
"""
    prefunc += ccg.c_codegen(
        [exactsoln.uu_exactsoln, exactsoln.vv_exactsoln],
        ["*exact_soln_UUGF", "*exact_soln_VVGF"],
        verbose=False,
        include_braces=False,
    )
    prefunc += "}\n"

    prefunc += r"""
// Exact solution at a single point around the origin.
static void WaveToy_exact_solution_single_point_r0(const CCTK_REAL time, const CCTK_REAL xCart0, const CCTK_REAL xCart1, const CCTK_REAL xCart2,
    CCTK_REAL *restrict exact_soln_UUGF, CCTK_REAL *restrict exact_soln_VVGF) {
DECLARE_CCTK_PARAMETERS;
"""
    prefunc += ccg.c_codegen(
        [exactsoln.uu_exactsoln_r0, exactsoln.vv_exactsoln_r0],
        ["*exact_soln_UUGF", "*exact_soln_VVGF"],
        verbose=False,
        include_braces=False,
    )
    prefunc += "}\n"

    gri.register_gridfunctions(["uu_exact", "vv_exact"], group="AUX")
    desc = r"""Set the exact solution at all grid points."""
    c_type = "void"
    name = f"{thorn_name}_exact_solution_all_points"
    params = "CCTK_ARGUMENTS"
    body = f"DECLARE_CCTK_ARGUMENTS_{name};\n"

    x_gf_access = gri.ETLegacyGridFunction.access_gf("x", use_GF_suffix=False)
    y_gf_access = gri.ETLegacyGridFunction.access_gf("y", use_GF_suffix=False)
    z_gf_access = gri.ETLegacyGridFunction.access_gf("z", use_GF_suffix=False)
    uuGF = "uu"
    vvGF = "vv"
    if thorn_name == diag_thorn_name:
        uuGF = "uu_exact"
        vvGF = "vv_exact"

    uu_exact_gf_access = gri.ETLegacyGridFunction.access_gf(uuGF)
    vv_exact_gf_access = gri.ETLegacyGridFunction.access_gf(vvGF)

    body += lp.simple_loop(
        f"""CCTK_REAL local_x = {x_gf_access};
  CCTK_REAL local_y = {y_gf_access};
  CCTK_REAL local_z = {z_gf_access};
  if(local_x*local_x + local_y*local_y + local_z*local_z > 1e-20) {{
  WaveToy_exact_solution_single_point(cctk_time, local_x, local_y,
                                      local_z, &{uu_exact_gf_access}, &{vv_exact_gf_access});
}} else {{
  WaveToy_exact_solution_single_point_r0(cctk_time, local_x, local_y,
                                      local_z, &{uu_exact_gf_access}, &{vv_exact_gf_access});
}}
""",
        loop_region="all points",
    )

    schedule_bin = "CCTK_INITIAL"
    if thorn_name == diag_thorn_name:
        schedule_bin = "CCTK_ANALYSIS"
    ET_schedule_bin_entry = (
        schedule_bin,
        f"""
schedule FUNC_NAME IN {schedule_bin}
{{
  LANG: C
  READS: Grid::x(Everywhere)
  READS: Grid::y(Everywhere)
  READS: Grid::z(Everywhere)
  WRITES: {evol_thorn_name}::{uuGF}GF(Everywhere)
  WRITES: {evol_thorn_name}::{vvGF}GF(Everywhere)
}} "Set up metric fields for binary black hole initial data"
""",
    )
    ET_current_thorn_CodeParams_used = None
    ET_other_thorn_CodeParams_used = None
    if thorn_name == ID_thorn_name:
        ET_current_thorn_CodeParams_used = ["sigma", "wavespeed"]
    if thorn_name == diag_thorn_name:
        ET_other_thorn_CodeParams_used = ["sigma", "wavespeed"]

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[ET_schedule_bin_entry],
        ET_current_thorn_CodeParams_used=ET_current_thorn_CodeParams_used,
        ET_other_thorn_CodeParams_used=ET_other_thorn_CodeParams_used,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_rhs_eval(thorn_name: str) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the wave equation with specific parameters.

    :param thorn_name: The name of the thorn for which the right-hand side evaluation function is being registered.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = standard_ET_includes
    if enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
    desc = r"""Set RHSs for wave equation."""
    c_type = "void"
    name = f"{thorn_name}_rhs_eval"
    params = "CCTK_ARGUMENTS"
    # Populate uu_rhs, vv_rhs
    rhs = WaveEquation_RHSs()

    if enable_KreissOliger_dissipation:
        diss_strength = par.register_CodeParameter(
            "REAL",
            __name__,
            "KreissOliger_diss_strength",
            0.9,
            commondata=True,
        )
        uu_dKOD = ixp.declarerank1("uu_dKOD")
        vv_dKOD = ixp.declarerank1("vv_dKOD")
        for k in range(3):
            rhs.uu_rhs += diss_strength * uu_dKOD[k]
            rhs.vv_rhs += diss_strength * vv_dKOD[k]
        print(
            "WARNING: Kreiss-Oliger has been enabled.\n"
            "  In your parfile, you will need to increase by 1 the following parameters:\n"
            "   driver::ghost_size and all the\n"
            "   CoordBase::boundary_size_?_*"
        )

    body = f"DECLARE_CCTK_ARGUMENTS_{name};\n"
    body += CodeParameters.read_CodeParameters(
        list_of_tuples__thorn_CodeParameter=[(ID_thorn_name, "wavespeed")],
        enable_simd=enable_simd,
        declare_invdxxs=True,
    )
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.uu_rhs, rhs.vv_rhs],
            [
                gri.ETLegacyGridFunction.access_gf("uu_rhs"),
                gri.ETLegacyGridFunction.access_gf("vv_rhs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
    )
    ET_schedule_bin_entry = (
        "MoL_CalcRHS",
        """
schedule FUNC_NAME in MoL_CalcRHS as rhs_eval
{
  LANG: C
  READS: evol_variables(everywhere)
  WRITES: evol_variables_rhs(interior)
} "MoL: Evaluate WaveToy RHSs"
""",
    )

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[ET_schedule_bin_entry],
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Parallel codegen:
for thorn in [ID_thorn_name, diag_thorn_name]:
    register_CFunction_exact_solution_all_points(
        thorn_name=thorn, in_WaveType=WaveType, in_default_sigma=default_sigma
    )
register_CFunction_rhs_eval(thorn_name=evol_thorn_name)

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

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
STORAGE: aux_variables[3]      # Diagnostics variables
""",
)
interface_ccl.construct_interface_ccl(
    project_dir=project_dir,
    thorn_name=evol_thorn_name,
    inherits="Boundary Grid MethodofLines",
    USES_INCLUDEs="""USES INCLUDE: Symmetry.h
USES INCLUDE: Boundary.h
""",
    is_evol_thorn=True,
    enable_NewRad=False,
)
CParams_registered_to_params_ccl += param_ccl.construct_param_ccl(
    project_dir=project_dir,
    thorn_name=evol_thorn_name,
    shares_extends_str=f"""shares: {ID_thorn_name}
USES CCTK_REAL wavespeed""",
)

# CCL files: ID_thorn
schedule_ccl.construct_schedule_ccl(
    project_dir=project_dir,
    thorn_name=ID_thorn_name,
    STORAGE="""
STORAGE: evol_variables[3] # Evolution variables
""",
)
interface_ccl.construct_interface_ccl(
    project_dir=project_dir,
    thorn_name=ID_thorn_name,
    inherits=f"""Grid {evol_thorn_name} # {evol_thorn_name} provides all gridfunctions.""",
    USES_INCLUDEs="",
    is_evol_thorn=False,
    enable_NewRad=False,
)
CParams_registered_to_params_ccl += param_ccl.construct_param_ccl(
    project_dir=project_dir,
    thorn_name=ID_thorn_name,
    shares_extends_str="",
)

# CCL files: diagnostics thorn
schedule_ccl.construct_schedule_ccl(
    project_dir=project_dir,
    thorn_name=diag_thorn_name,
    STORAGE="""
STORAGE: aux_variables[3] # Diagnostics variables
""",
)
interface_ccl.construct_interface_ccl(
    project_dir=project_dir,
    thorn_name=diag_thorn_name,
    inherits=f"""Grid {evol_thorn_name} # {evol_thorn_name} provides all gridfunctions.""",
    USES_INCLUDEs="",
    is_evol_thorn=False,
    enable_NewRad=False,
)
CParams_registered_to_params_ccl += param_ccl.construct_param_ccl(
    project_dir=project_dir,
    thorn_name=diag_thorn_name,
    # FIXME: the following line generation could be automated:
    shares_extends_str=f"""
shares: {ID_thorn_name}
USES CCTK_REAL sigma
USES CCTK_REAL wavespeed
""",
)


for thorn in [evol_thorn_name, ID_thorn_name, diag_thorn_name]:
    make_code_defn.output_CFunctions_and_construct_make_code_defn(
        project_dir=project_dir, thorn_name=thorn
    )
simd.copy_simd_intrinsics_h(
    project_dir=str(Path(project_dir) / evol_thorn_name / "src")
)

"""
Spinning black hole example.

Specifically, evolve a spinning black hole in the UIUC initial data slicing,
  Liu, Etienne, & Shapiro (2009) https://arxiv.org/pdf/1001.4077.pdf;
  this example sets up a complete C code for solving the GR field
  equations in curvilinear coordinates on a cell-centered grid,
  using a reference metric approach.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import os

import nrpy.params as par
from nrpy.helpers import simd
import nrpy.helpers.parallel_codegen as pcg

from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
from nrpy.infrastructures.BHaH import rfm_precompute
from nrpy.infrastructures.BHaH import rfm_wrapper_functions
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numericalgrids
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCartxx
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "spinning_blackhole"
CoordSystem = "Cartesian"
IDtype = "UIUCBlackHole"
IDCoordSystem = "Spherical"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 1.0
grid_physical_size = 7.5
diagnostics_output_every = 0.25
t_final = 1.0 * grid_physical_size
Nxx_dict = {
    "Spherical": [72, 12, 2],
    "SinhSpherical": [72, 12, 2],
    "Cartesian": [64, 64, 64],
}
default_BH_mass = 1.0
default_BH_spin = +0.8
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
enable_simd = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
enable_fd_functions = True
enable_KreissOliger_dissipation = False
boundary_conditions_desc = "outgoing radiation"

OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    OMP_collapse = 2  # about 2x faster
if "Cylindrical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "1")
    OMP_collapse = 2  # might be slightly faster
project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("enable_RbarDD_gridfunctions", separate_Ricci_and_BSSN_RHS)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
BCl.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    ID_persist_struct_str="",
)

numericalgrids.register_CFunctions(
    list_of_CoordSystems=[CoordSystem],
    grid_physical_size=grid_physical_size,
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
BCl.register_CFunction_diagnostics(
    list_of_CoordSystems=[CoordSystem],
    default_diagnostics_out_every=diagnostics_output_every,
    enable_psi4_diagnostics=False,
    grid_center_filename_tuple=("out0d-conv_factor%.2f.txt", "convergence_factor"),
    axis_filename_tuple=(
        "out1d-AXIS-conv_factor%.2f-t%08.2f.txt",
        "convergence_factor, time",
    ),
    plane_filename_tuple=(
        "out2d-PLANE-conv_factor%.2f-t%08.2f.txt",
        "convergence_factor, time",
    ),
    out_quantities_dict="default",
)

if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(
        list_of_CoordSystems=[CoordSystem]
    )
BCl.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_Ricci_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_constraints(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

if __name__ == "__main__":
    pcg.do_parallel_codegen()

cbc.CurviBoundaryConditions_register_C_functions(
    list_of_CoordSystems=[CoordSystem], radiation_BC_fd_order=radiation_BC_fd_order
)
rhs_string = """
Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")

MoL.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(commondata, params, rfmstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
xxCartxx.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
xxCartxx.register_CFunction_xx_to_Cart(CoordSystem)
progress.register_CFunction_progress_indicator()
rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)
par.adjust_CodeParam_default("M", default_BH_mass)
par.adjust_CodeParam_default("chi", default_BH_spin)

CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_simd=enable_simd,
    fin_NGHOSTS_add_one_for_upwinding=True,
)

main.register_CFunction_main_c(
    initial_data_desc=IDtype,
    MoL_method=MoL_method,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    boundary_conditions_desc=boundary_conditions_desc,
)

if enable_simd:
    simd.copy_simd_intrinsics_h(project_dir=project_dir)
Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    compiler_opt_option="default",
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

# print(cfc.CFunction_dict["initial_data"].full_function)
# print(cfc.CFunction_dict["rhs_eval"].full_function)
# print(cfc.CFunction_dict["apply_bcs"].full_function)
# print(cfc.CFunction_dict["parameter_file_read_and_parse"].full_function)
# print(cfc.CFunction_dict["main"].full_function)

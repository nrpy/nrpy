"""
Sets up a complete C code project for solving the wave equation in curvilinear coordinates on a cell-centered grid, using a reference metric.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil

import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    CodeParameters,
    Makefile_helpers,
    checkpointing,
    cmdline_input_and_parfiles,
    griddata_commondata,
    main_c,
    numerical_grids_and_timestep,
    rfm_precompute,
    rfm_wrapper_functions,
    wave_equation,
    xx_tofrom_Cart,
)
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "curviwavetoy"
WaveType = "SphericalGaussian"
default_sigma = 3.0
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size
default_diagnostics_output_every = 0.5
default_checkpoint_every = 50.0
CoordSystem = "Spherical"
# symmetry_axes will be set on any i such that Nxx[i] = 2 below.
Nxx_dict = {
    "Spherical": [64, 2, 2],
    "SinhSpherical": [64, 2, 2],
    "Cartesian": [64, 64, 64],
    "SinhCartesian": [64, 64, 64],
}
OMP_collapse = 1
if (
    "Spherical" in CoordSystem
    and WaveType == "SphericalGaussian"
    and Nxx_dict["Spherical"][1] == Nxx_dict["Spherical"][2] == 2
):
    par.set_parval_from_str("symmetry_axes", "12")
    OMP_collapse = 2  # about 2x faster

enable_rfm_precompute = False
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 2
enable_simd = False
enable_KreissOliger_dissipation = False
parallel_codegen_enable = True
boundary_conditions_desc = "outgoing radiation"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
wave_equation.initial_data_exact_soln.register_CFunction_exact_solution_single_Cartesian_point(
    WaveType=WaveType, default_sigma=default_sigma
)
wave_equation.initial_data_exact_soln.register_CFunction_initial_data(
    enable_checkpointing=True, OMP_collapse=OMP_collapse
)
numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

wave_equation.diagnostics.register_CFunction_diagnostics(
    set_of_CoordSystems={CoordSystem},
    default_diagnostics_out_every=default_diagnostics_output_every,
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
    rfm_precompute.register_CFunctions_rfm_precompute(set_of_CoordSystems={CoordSystem})
wave_equation.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    OMP_collapse=OMP_collapse,
)


if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

cbc.CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems={CoordSystem}, radiation_BC_fd_order=radiation_BC_fd_order
)
rhs_string = """rhs_eval(commondata, params, rfmstruct,  RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
MoL_register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
checkpointing.register_CFunctions(default_checkpoint_every=default_checkpoint_every)

progress.register_CFunction_progress_indicator()
rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
CodeParameters.write_CodeParameters_h_files(project_dir=project_dir)
CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
cmdline_input_and_parfiles.generate_default_parfile(
    project_dir=project_dir, project_name=project_name
)
cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_simd,
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=enable_KreissOliger_dissipation,
)
main_c.register_CFunction_main_c(
    initial_data_desc=WaveType,
    pre_MoL_step_forward_in_time="write_checkpoint(&commondata, griddata);\n",
    MoL_method=MoL_method,
    boundary_conditions_desc=boundary_conditions_desc,
)
griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

if enable_simd:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

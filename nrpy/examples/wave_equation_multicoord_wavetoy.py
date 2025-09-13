"""
Sets up a complete C code project for solving the wave equation in curvilinear coordinates on a cell-centered grid, using a reference metric.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import argparse
import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

parser = argparse.ArgumentParser(
    description="NRPyElliptic Solver for Conformally Flat BBH initial data"
)
parser.add_argument(
    "--parallelization",
    type=str,
    help="Parallelization strategy to use (e.g. openmp, cuda).",
    default="openmp",
)
parser.add_argument(
    "--floating_point_precision",
    type=str,
    help="Floating point precision (e.g. float, double).",
    default="double",
)
parser.add_argument(
    "--disable_intrinsics",
    action="store_true",
    help="Flag to disable hardware intrinsics",
    default=False,
)
parser.add_argument(
    "--disable_rfm_precompute",
    action="store_true",
    help="Flag to disable RFM precomputation.",
    default=False,
)
args = parser.parse_args()

# Code-generation-time parameters:
fp_type = args.floating_point_precision.lower()
parallelization = args.parallelization.lower()
enable_intrinsics = not args.disable_intrinsics
enable_rfm_precompute = not args.disable_rfm_precompute

if parallelization not in ["openmp", "cuda"]:
    raise ValueError(
        f"Invalid parallelization strategy: {parallelization}. "
        "Choose 'openmp' or 'cuda'."
    )

par.set_parval_from_str("Infrastructure", "BHaH")
par.set_parval_from_str("parallelization", parallelization)
par.set_parval_from_str("fp_type", fp_type)

# Code-generation-time parameters:
project_name = "multicoords_curviwavetoy"
WaveType = "SphericalGaussian"
default_sigma = 3.0
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size
default_diagnostics_output_every = 0.2
default_checkpoint_every = 2.0
set_of_CoordSystems = {"Spherical", "SinhSpherical", "Cartesian", "SinhCartesian"}
# Must set this for all multi-patch/multi-coordinate systems. Otherwise only one CoordSystem will be registered in params.
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", "All")
list_of_grid_physical_sizes = []
for CoordSystem in set_of_CoordSystems:
    list_of_grid_physical_sizes.append(grid_physical_size)
NUMGRIDS = len(set_of_CoordSystems)
num_cuda_streams = NUMGRIDS

Nxx_dict = {
    "Spherical": [64, 2, 2],
    "SinhSpherical": [64, 2, 2],
    "Cartesian": [64, 64, 64],
    "SinhCartesian": [64, 64, 64],
}
OMP_collapse = 1
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 2
enable_KreissOliger_dissipation = False
parallel_codegen_enable = True
boundary_conditions_desc = "outgoing radiation"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.adjust_CodeParam_default("NUMGRIDS", NUMGRIDS)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

BHaH.wave_equation.initial_data_exact_soln.register_CFunction_initial_data_exact(
    OMP_collapse=OMP_collapse, WaveType=WaveType, default_sigma=default_sigma
)
BHaH.wave_equation.initial_data_exact_soln.register_CFunction_initial_data(
    enable_checkpointing=True,
)

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=list_of_grid_physical_sizes,
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

for CoordSystem in set_of_CoordSystems:
    BHaH.wave_equation.rhs_eval.register_CFunction_rhs_eval(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_intrinsics=enable_intrinsics,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        OMP_collapse=OMP_collapse,
    )
    BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

BHaH.wave_equation.diagnostics.register_CFunction_diagnostics(
    set_of_CoordSystems=set_of_CoordSystems,
    default_diagnostics_out_every=default_diagnostics_output_every,
    grid_center_filename_tuple=(
        "out0d-%s-conv_factor-%.2f.txt",
        "CoordSystemName, convergence_factor",
    ),
    axis_filename_tuple=(
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    plane_filename_tuple=(
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    out_quantities_dict="default",
)

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems=set_of_CoordSystems
    )

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
)
rhs_string = """rhs_eval(commondata, params, rfmstruct,  RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
BHaH.MoLtimestepping.register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
BHaH.checkpointing.register_CFunctions(
    default_checkpoint_every=default_checkpoint_every
)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
par.adjust_CodeParam_default("t_final", t_final)
BHaH.CodeParameters.write_CodeParameters_h_files(project_dir=project_dir)
BHaH.CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
BHaH.cmdline_input_and_parfiles.generate_default_parfile(
    project_dir=project_dir, project_name=project_name
)
BHaH.cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
gpu_defines_filename = BHaH.BHaH_device_defines_h.output_device_headers(
    project_dir, num_streams=num_cuda_streams
)
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_intrinsics,
    intrinsics_header_lst=(
        ["cuda_intrinsics.h"] if parallelization == "cuda" else ["simd_intrinsics.h"]
    ),
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=enable_KreissOliger_dissipation,
    DOUBLE_means="double" if fp_type == "float" else "REAL",
    restrict_pointer_type="*" if parallelization == "cuda" else "*restrict",
    supplemental_defines_dict=(
        {
            "C++/CUDA safe restrict": "#define restrict __restrict__",
            "GPU Header": f'#include "{gpu_defines_filename}"',
        }
        if parallelization == "cuda"
        else {}
    ),
)

compute_griddata = "griddata_device" if parallelization in ["cuda"] else "griddata"
write_checkpoint_call = f"write_checkpoint(&commondata, {compute_griddata});\n".replace(
    compute_griddata,
    (
        f"griddata_host, {compute_griddata}"
        if parallelization in ["cuda"]
        else compute_griddata
    ),
)

BHaH.main_c.register_CFunction_main_c(
    initial_data_desc=WaveType,
    MoL_method=MoL_method,
    pre_MoL_step_forward_in_time=write_checkpoint_call,
    boundary_conditions_desc=boundary_conditions_desc,
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

if enable_intrinsics:
    copy_files(
        package="nrpy.helpers",
        filenames_list=(
            ["cuda_intrinsics.h"]
            if parallelization == "cuda"
            else ["simd_intrinsics.h"]
        ),
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

cuda_makefiles_options = (
    {
        "CC": "nvcc",
        "src_code_file_ext": "cu",
        "compiler_opt_option": "nvcc",
    }
    if parallelization == "cuda"
    else {}
)
if parallelization == "cuda":
    BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=project_name,
        CC="nvcc",
        src_code_file_ext="cu",
        compiler_opt_option="nvcc",
    )
else:
    BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=project_name,
    )
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

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
    "--floating_point_precision",
    type=str,
    help="Floating point precision (e.g. float, double).",
    default="double",
)
parser.add_argument(
    "--cuda",
    action="store_true",
    help="Use CUDA parallelization.",
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
enable_intrinsics = not args.disable_intrinsics
enable_rfm_precompute = not args.disable_rfm_precompute

# Default to openmp; override with cuda if --cuda is set
parallelization = "cuda" if args.cuda else "openmp"
if parallelization not in ["openmp", "cuda"]:
    raise ValueError(
        f"Invalid parallelization strategy: {parallelization}. "
        "Choose 'openmp' or 'cuda'."
    )

par.set_parval_from_str("Infrastructure", "BHaH")
par.set_parval_from_str("parallelization", parallelization)
par.set_parval_from_str("fp_type", fp_type)

# Code-generation-time parameters:
project_name = "wave_equation_curvilinear"
WaveType = "SphericalGaussian"
default_sigma = 3.0
default_k0, default_k1, default_k2 = (1.0, 1.0, 1.0)
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size
default_diagnostics_output_every = 0.5
default_checkpoint_every = 50.0
CoordSystem = "SinhCylindrical"
set_of_CoordSystems = {CoordSystem}
list_of_grid_physical_sizes = []
for CoordSystem in set_of_CoordSystems:
    list_of_grid_physical_sizes.append(grid_physical_size)
NUMGRIDS = len(set_of_CoordSystems)
num_cuda_streams = NUMGRIDS
# symmetry_axes will be set on any i such that Nxx[i] = 2 below.
Nxx_dict = {
    "Spherical": [64, 2, 2],
    "SinhSpherical": [64, 2, 2],
    "SinhCylindrical": [64, 2, 64],
    "Cartesian": [64, 64, 64],
    "SinhCartesian": [64, 64, 64],
}
OMP_collapse = 1
if (
    "Spherical" in CoordSystem
    and WaveType == "SphericalGaussian"
    and Nxx_dict[CoordSystem][1] == Nxx_dict[CoordSystem][2] == 2
):
    par.set_parval_from_str("symmetry_axes", "12")
    OMP_collapse = 2  # about 2x faster
if (
    "Cylindrical" in CoordSystem
    and WaveType == "SphericalGaussian"
    and Nxx_dict[CoordSystem][1] == 2
):
    par.set_parval_from_str("symmetry_axes", "1")

MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 2
enable_KreissOliger_dissipation = False
enable_parallel_codegen = True
boundary_conditions_desc = "outgoing radiation"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.adjust_CodeParam_default("NUMGRIDS", NUMGRIDS)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

BHaH.wave_equation.initial_data_exact_soln.register_CFunction_initial_data_exact(
    OMP_collapse=OMP_collapse,
    WaveType=WaveType,
    default_sigma=default_sigma,
    default_k0=default_k0,
    default_k1=default_k1,
    default_k2=default_k2,
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

# Diagnostics C code registration
BHaH.diagnostics.diagnostics.register_all_diagnostics(
    set_of_CoordSystems=set_of_CoordSystems,
    project_dir=project_dir,
    default_diagnostics_out_every=default_diagnostics_output_every,
    enable_nearest_diagnostics=True,
    enable_interp_diagnostics=False,
    enable_volume_integration_diagnostics=True,
    enable_free_auxevol=False,
)
BHaH.wave_equation.diagnostics_volume_integration.register_CFunction_diagnostics_volume_integration()
BHaH.wave_equation.diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
    WaveType=WaveType,
    default_k0=default_k0,
    default_k1=default_k1,
    default_k2=default_k2,
    default_sigma=default_sigma,
)
BHaH.wave_equation.diagnostics_nearest.register_CFunction_diagnostics_nearest()

if __name__ == "__main__" and enable_parallel_codegen:
    pcg.do_parallel_codegen()

#########################################################
# STEP 3 (post parallel codegen): Generate header files,
#         register remaining C functions and command-line
#         parameters, set up boundary conditions, and
#         create a Makefile for this project.
#         Project is output to project/[project_name]/

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

par.adjust_CodeParam_default("t_final", t_final)

BHaH.diagnostics.diagnostic_gfs_h_create.diagnostics_gfs_h_create(
    project_dir=project_dir,
)

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
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=enable_KreissOliger_dissipation,
    DOUBLE_means="double" if fp_type == "float" else "REAL",
    restrict_pointer_type="*" if parallelization == "cuda" else "*restrict",
)

BHaH.main_c.register_CFunction_main_c(
    initial_data_desc=WaveType,
    MoL_method=MoL_method,
    pre_MoL_step_forward_in_time=(
        f"write_checkpoint(&commondata, "
        f"{'griddata_host, griddata_device' if parallelization == 'cuda' else 'griddata'});\n"
    ),
    boundary_conditions_desc=boundary_conditions_desc,
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

# SIMD intrinsics needed for 3D interpolation, constraints evaluation, etc.
intrinsics_file_list = ["simd_intrinsics.h"]
if parallelization == "cuda":
    # CUDA intrinsics needed for CUDA-enabled projects.
    intrinsics_file_list += ["cuda_intrinsics.h"]
copy_files(
    package="nrpy.helpers",
    filenames_list=intrinsics_file_list,
    project_dir=project_dir,
    subdirectory="intrinsics",
)

BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    compiler_opt_option=("nvcc" if parallelization == "cuda" else "default"),
    CC=("nvcc" if parallelization == "cuda" else "autodetect"),
    src_code_file_ext=("cu" if parallelization == "cuda" else "c"),
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

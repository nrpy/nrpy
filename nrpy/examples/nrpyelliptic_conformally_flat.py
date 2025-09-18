"""
Sets up a complete C code project for solving the hyperbolic relaxation equation in curvilinear coordinates on a cell-centered grid, using a reference metric.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
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

project_name = "nrpyelliptic_conformally_flat"
par.set_parval_from_str("fp_type", fp_type)
par.set_parval_from_str("parallelization", parallelization)
par.set_parval_from_str("Infrastructure", "BHaH")

grid_physical_size = 1.0e6
t_final = grid_physical_size  # This parameter is effectively not used in NRPyElliptic
nn_max = 10000  # Sets the maximum number of relaxation steps


def get_log10_residual_tolerance(fp_type_str: str = "double") -> float:
    """
    Determine the residual tolerance based on the fp_precision.

    :param fp_type_str: string representing the floating point type.
    :return: float of the residual tolerance based on fp_type.
    :raises ValueError: If the input fp_type_str branch is not defined.
    """
    res: float = -1.0
    if fp_type_str == "double":
        res = -15.8
    elif fp_type_str == "float":
        res = -10.0
    else:
        raise ValueError(f"residual tolerence not defined for {fp_type_str} precision")
    return res


# Set tolerance for log10(residual) to stop relaxation
log10_residual_tolerance = get_log10_residual_tolerance(fp_type_str=fp_type)
default_diagnostics_output_every = 100
default_checkpoint_every = 50.0
eta_damping = 11.0
MINIMUM_GLOBAL_WAVESPEED = 0.7
CFL_FACTOR = 1.0  # NRPyElliptic wave speed prescription assumes this parameter is ALWAYS set to 1
# CoordSystem = "SinhSpherical"
CoordSystem = "SinhSymTP"
# CoordSystem = "SymTP"
Nxx_dict = {
    "SymTP": [128, 128, 16],
    "SinhSymTP": [128, 128, 16],
    "SinhCylindricalv2": [128, 16, 256],
    "SinhSpherical": [128, 64, 16],
}
# Set parameters specific to SinhSymTP coordinates
AMAX = grid_physical_size
bScale = 5.0
SINHWAA = 0.07
# Set parameters specific to SinhCylindricalv2 coordinates
AMPLRHO = grid_physical_size
AMPLZ = grid_physical_size
SINHWRHO = 0.04
SINHWZ = 0.04
const_drho = 2.0e-5
const_dz = 2.0e-5
# Set parameters specific to SinhSpherical coordinates
AMPL = grid_physical_size
SINHW = 0.06

OMP_collapse = 1
enable_checkpointing = True
MoL_method = "RK4"
fd_order = 10
radiation_BC_fd_order = 6
parallel_codegen_enable = True
boundary_conditions_desc = "outgoing radiation"
set_of_CoordSystems = {CoordSystem}
NUMGRIDS = len(set_of_CoordSystems)
num_cuda_streams = NUMGRIDS
par.adjust_CodeParam_default("NUMGRIDS", NUMGRIDS)
# fmt: off
initial_data_type = "gw150914"  # choices are: "gw150914", "axisymmetric", and "single_puncture"

q = 36.0 / 29.0
Pr = -0.00084541526517121  # Radial linear momentum
Pphi = 0.09530152296974252  # Azimuthal linear momentum
S0_y_dimless = 0.31
S1_y_dimless = -0.46
m0_adm = q / (1.0 + q)
m1_adm = 1.0 / (1.0 + q)

gw150914_params = {
    "zPunc": 5.0,
    "q": q,
    "bare_mass_0": 0.51841993533587039,
    "bare_mass_1": 0.39193567996522616,
    "Pr": Pr,
    "Pphi": Pphi,
    "S0_y_dimless": S0_y_dimless,
    "S1_y_dimless": S1_y_dimless,
    "m0_adm": m0_adm,
    "m1_adm": m1_adm,
    "S0_y": S0_y_dimless * (m0_adm ** 2),
    "S1_y": S1_y_dimless * (m1_adm ** 2),
    "P0_x": Pphi,
    "P0_z": Pr,
    "P1_x": -Pphi,
    "P1_z": -Pr,
}

axisymmetric_params = {
    "zPunc": 5.0,
    "bare_mass_0": 0.5,
    "bare_mass_1": 0.5,
    "S0_z": +0.2,
    "S1_z": -0.2,
}

single_puncture_params = {
    "zPunc": 0.0,
    "bare_mass_0": 0.5,
    "S0_z": 0.2,
}
# fmt: on

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

# Generate functions to set initial guess
BHaH.nrpyelliptic.initial_data.register_CFunction_initial_guess_single_point()
BHaH.nrpyelliptic.initial_data.register_CFunction_initial_guess_all_points(
    OMP_collapse=OMP_collapse,
    enable_checkpointing=enable_checkpointing,
)

# Generate function that calls functions to set variable wavespeed and all other AUXEVOL gridfunctions
for CoordSystem in set_of_CoordSystems:
    BHaH.nrpyelliptic.auxevol_gfs_set_to_constant.register_CFunction_auxevol_gfs_set_to_constant(
        CoordSystem, OMP_collapse=OMP_collapse, enable_intrinsics=enable_intrinsics
    )

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size for c in set_of_CoordSystems],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

BHaH.nrpyelliptic.diagnostics.register_CFunction_diagnostics(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    default_diagnostics_out_every=default_diagnostics_output_every,
)

if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(set_of_CoordSystems)

# Generate function to compute RHSs
BHaH.nrpyelliptic.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_intrinsics=enable_intrinsics,
    OMP_collapse=OMP_collapse,
)

# Generate function to compute residuals
BHaH.nrpyelliptic.residual_H_compute_all_points.register_CFunction_residual_H_compute_all_points(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_intrinsics=enable_intrinsics,
    OMP_collapse=OMP_collapse,
)

# Generate diagnostics functions
BHaH.nrpyelliptic.log10_L2norm_gf.register_CFunction_log10_L2norm_gf(
    CoordSystem=CoordSystem
)

# Register function to check for stop conditions
BHaH.nrpyelliptic.stop_conditions_check.register_CFunction_stop_conditions_check()

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
)
wave_speed_assignment = (
    "cudaMemcpy(&wavespeed_at_outer_boundary, &auxevol_gfs[IDX4P(params, VARIABLE_WAVESPEEDGF, params->Nxx_plus_2NGHOSTS0 - NGHOSTS - 1, NGHOSTS, params->Nxx_plus_2NGHOSTS2 / 2)], sizeof(REAL), cudaMemcpyDeviceToHost);"
    if parallelization == "cuda"
    else "wavespeed_at_outer_boundary = auxevol_gfs[IDX4P(params, VARIABLE_WAVESPEEDGF, params->Nxx_plus_2NGHOSTS0 - NGHOSTS - 1, NGHOSTS, params->Nxx_plus_2NGHOSTS2 / 2)];"
)
rhs_string = f"""rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0) {{
  REAL wavespeed_at_outer_boundary;
  {wave_speed_assignment}
  const REAL custom_gridfunctions_wavespeed[2] = {{wavespeed_at_outer_boundary, wavespeed_at_outer_boundary}};
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata->xx,
                                     custom_gridfunctions_wavespeed, gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);
}}"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
BHaH.MoLtimestepping.register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
    enable_intrinsics=False,  # SIMD intrinsics in MoL is not properly supported -- MoL update loops are not properly bounds checked.
    rational_const_alias="static constexpr" if parallelization == "cuda" else "const",
)
BHaH.checkpointing.register_CFunctions(
    default_checkpoint_every=default_checkpoint_every
)

# Define string with print statement for progress indicator
progress_str = r"""
  fprintf(stderr, "nn / nn_max = %d / %d ; log10(residual) / log10(residual_target) =  %.4f / %.4f \r",
    commondata->nn,
    commondata->nn_max,
    commondata->log10_current_residual,
    commondata->log10_residual_tolerance);
  fflush(stderr); // Flush the stderr buffer
"""
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator(
    progress_str=progress_str, compute_ETA=False
)
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

# Update parameters needed for hyperbolic relaxation method
par.adjust_CodeParam_default("t_final", t_final)
par.adjust_CodeParam_default("eta_damping", eta_damping)
par.adjust_CodeParam_default("MINIMUM_GLOBAL_WAVESPEED", MINIMUM_GLOBAL_WAVESPEED)
par.adjust_CodeParam_default("CFL_FACTOR", CFL_FACTOR)
par.adjust_CodeParam_default("nn_max", nn_max)
par.adjust_CodeParam_default("log10_residual_tolerance", log10_residual_tolerance)

# Update parameters specific to the coordinate system
if CoordSystem == "SinhSymTP":
    par.adjust_CodeParam_default("AMAX", AMAX)
    par.adjust_CodeParam_default("bScale", bScale)
    par.adjust_CodeParam_default("SINHWAA", SINHWAA)

if CoordSystem == "SinhCylindricalv2":
    par.adjust_CodeParam_default("AMPLRHO", AMPLRHO)
    par.adjust_CodeParam_default("AMPLZ", AMPLZ)
    par.adjust_CodeParam_default("SINHWRHO", SINHWRHO)
    par.adjust_CodeParam_default("SINHWZ", SINHWZ)
    par.adjust_CodeParam_default("const_drho", const_drho)
    par.adjust_CodeParam_default("const_dz", const_dz)

if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("AMPL", AMPL)
    par.adjust_CodeParam_default("SINHW", SINHW)

# Update parameters specific to initial data type
if initial_data_type == "gw150914":
    for param, value in gw150914_params.items():
        if param in [
            "zPunc",
            "bare_mass_0",
            "bare_mass_1",
            "S0_y",
            "S1_y",
            "P0_x",
            "P0_z",
            "P1_x",
            "P1_z",
        ]:
            par.adjust_CodeParam_default(param, value)

if initial_data_type == "single_puncture":
    for param, value in single_puncture_params.items():
        if param in [
            "zPunc",
            "bare_mass_0",
            "S0_z",
        ]:
            par.adjust_CodeParam_default(param, value)

if initial_data_type == "axisymmetric":
    for param, value in axisymmetric_params.items():
        if param in [
            "zPunc",
            "bare_mass_0",
            "bare_mass_1",
            "S0_z",
            "S1_z",
        ]:
            par.adjust_CodeParam_default(param, value)

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/

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
    DOUBLE_means="double" if fp_type == "float" else "REAL",
    restrict_pointer_type="*" if parallelization == "cuda" else "*restrict",
    supplemental_defines_dict=(
        {
            "ADDITIONAL GPU DIAGNOSTICS": """
#define L2_DVGF 0
#define L2_SQUARED_DVGF 1
""",
            "ADDITIONAL HOST DIAGNOSTICS": """
#define HOST_RESIDUAL_HGF 0
#define HOST_UUGF 1
#define NUM_HOST_DIAG 2
""",
            "C++/CUDA safe restrict": "#define restrict __restrict__",
            "GPU Header": f'#include "{gpu_defines_filename}"',
        }
        if parallelization == "cuda"
        else {}
    ),
)

# Set griddata struct used for calculations to griddata_device for certain parallelizations
compute_griddata = "griddata_device" if parallelization in ["cuda"] else "griddata"

# Define {pre,post}_MoL_step_forward_in_time string for main function
write_checkpoint_call = (
    f"write_checkpoint(&commondata, "
    f"{'griddata_host, griddata_device' if parallelization == 'cuda' else 'griddata'});\n"
)
pre_MoL_step_forward_in_time = write_checkpoint_call
post_MoL_step_forward_in_time = rf"""    stop_conditions_check(&commondata, {compute_griddata});
    if (commondata.stop_relaxation) {{
      // Force a checkpoint when stop condition is reached.
      commondata.checkpoint_every = 1e-4*commondata.dt;
      {write_checkpoint_call}
      break;
    }}
"""
post_non_y_n_auxevol_mallocs = (
    "for (int grid = 0; grid < commondata.NUMGRIDS; grid++)\n"
    f"auxevol_gfs_set_to_constant(&commondata, &{compute_griddata}[grid].params, {compute_griddata}[grid].xx, &{compute_griddata}[grid].gridfuncs);\n"
)
BHaH.main_c.register_CFunction_main_c(
    MoL_method=MoL_method,
    initial_data_desc="",
    boundary_conditions_desc=boundary_conditions_desc,
    post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
    pre_MoL_step_forward_in_time=pre_MoL_step_forward_in_time,
    post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
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

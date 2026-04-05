"""
Sets up a complete C code project for solving the hyperbolic relaxation equation in curvilinear coordinates on a cell-centered grid, using a reference metric.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Nishita Jadoo; njadoo **at** uidaho **dot* edu

Note: This is the superB version.
"""

import argparse
import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
from pathlib import Path

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH, superB

parser = argparse.ArgumentParser(
    description="NRPyElliptic Solver for Conformally Flat BBH initial data"
)
parser.add_argument(
    "--floating_point_precision",
    type=str,
    help="Floating point precision (e.g. float, double).",
    default="double",
)
args = parser.parse_args()

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
fp_type = args.floating_point_precision.lower()
project_name = "superB_nrpyelliptic_conformally_flat"
enable_simd_intrinsics = True
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
        res = -11.2
    elif fp_type_str == "float":
        res = -6.5
    else:
        raise ValueError(f"residual tolerance not defined for {fp_type_str} precision")
    return res


# Set tolerance for log10(residual) to stop relaxation
log10_residual_tolerance = get_log10_residual_tolerance(fp_type_str=fp_type)
default_diagnostics_output_every = 8e-2
default_checkpoint_every = 50.0
eta_damping = 11.0
MINIMUM_GLOBAL_WAVESPEED = 0.7
CFL_FACTOR = 1.0  # NRPyElliptic wave speed prescription assumes this parameter is ALWAYS set to 1
CoordSystem = "SinhSpherical"
set_of_CoordSystems = {CoordSystem}
Nxx_dict = {
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
enable_checkpointing = False
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 10
radiation_BC_fd_order = 6
enable_simd = True
enable_parallel_codegen = True
outer_bcs_type = "radiation"
boundary_conditions_desc = "outgoing radiation"
# fmt: off
initial_data_type = "gw150914"  # choices are: "gw150914", "axisymmetric", and "single_puncture"
# Number of chares, Nchare0, Nchare1, and Nchare2, in each direction,
# should be chosen such that Nxx0/Nchare0, Nxx1/Nchare1, Nxx2/Nchare2 are integers greater than NGHOSTS,
# NGHOSTS is fd_order/2
if "Spherical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 16)
    par.adjust_CodeParam_default("Nchare1", 2)
    par.adjust_CodeParam_default("Nchare2", 2)
if "Cylindrical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 16)
    par.adjust_CodeParam_default("Nchare1", 2)
    par.adjust_CodeParam_default("Nchare2", 32)

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

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]


# Generate functions to set initial guess
BHaH.nrpyelliptic.initial_data.register_CFunction_initial_guess_single_point()
BHaH.nrpyelliptic.initial_data.register_CFunction_initial_guess_all_points(
    OMP_collapse=OMP_collapse, enable_checkpointing=enable_checkpointing
)

# Generate function that calls functions to set variable wavespeed and all other AUXEVOL gridfunctions
for CoordSystem in set_of_CoordSystems:
    BHaH.nrpyelliptic.auxevol_gfs_set_to_constant.register_CFunction_auxevol_gfs_set_to_constant(
        CoordSystem, OMP_collapse=OMP_collapse
    )

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

superB.numerical_grids.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

superB.diagnostics.diagnostics.register_all_diagnostics(
    set_of_CoordSystems=set_of_CoordSystems,
    project_dir=project_dir,
    default_diagnostics_out_every=default_diagnostics_output_every,
    enable_nearest_diagnostics=True,
    enable_interp_diagnostics=False,
    enable_volume_integration_diagnostics=True,
    enable_free_auxevol=False,
)
BHaH.nrpyelliptic.diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics=False
)
superB.nrpyelliptic.diagnostics_nearest.register_CFunction_diagnostics_nearest(
    CoordSystem
)

if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems=set_of_CoordSystems
    )

# Generate function to compute RHSs
BHaH.nrpyelliptic.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_intrinsics=enable_simd,
    OMP_collapse=OMP_collapse,
)

# Generate function to compute residuals
BHaH.nrpyelliptic.residual_H_compute_all_points.register_CFunction_residual_H_compute_all_points(
    CoordSystem=CoordSystem,
    enable_simd_intrinsics=enable_simd_intrinsics,
    OMP_collapse=OMP_collapse,
)

# Register function to check for stop conditions
BHaH.nrpyelliptic.stop_conditions_check.register_CFunction_stop_conditions_check()

if __name__ == "__main__" and enable_parallel_codegen:
    pcg.do_parallel_codegen()

superB.chare_communication_maps.chare_comm_register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems
)
superB.CurviBoundaryConditions.CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
    set_parity_on_aux=True,
)

rhs_string = "rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);"
if outer_bcs_type == "radiation":
    rhs_string += """
const REAL wavespeed_at_outer_boundary = griddata[grid].params.wavespeed_at_outer_boundary;
const REAL custom_gridfunctions_wavespeed[2] = {wavespeed_at_outer_boundary, wavespeed_at_outer_boundary};
apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata->xx,
                                    custom_gridfunctions_wavespeed, gridfunctions_f_infinity,
                                    RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")

post_rhs_bcs_str = ""
if outer_bcs_type != "radiation":
    post_rhs_bcs_str += """
apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);"""

superB.MoL.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_bcs_str=post_rhs_bcs_str,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
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


BHaH.diagnostics.diagnostic_gfs_h_create.diagnostics_gfs_h_create(
    project_dir=project_dir,
)

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

# Define post_MoL_step_forward_in_time string for main function
post_MoL_step_forward_in_time = r"""
serial {
  stop_conditions_check(&commondata);
  if (commondata.stop_relaxation) {
    mainProxy.done();
  }
}
"""

superB.superB.superB_pup.register_CFunction_superB_pup_routines(
    MoL_method=MoL_method,
)
copy_files(
    package="nrpy.infrastructures.superB.superB",
    filenames_list=["superB.h", "superB_pup_function_prototypes.h"],
    project_dir=project_dir,
    subdirectory="superB",
)

superB.main_chare.output_commondata_object_h_and_main_h_cpp_ci(
    project_dir=project_dir,
)
superB.timestepping_chare.output_timestepping_h_cpp_ci_register_CFunctions(
    project_dir=project_dir,
    MoL_method=MoL_method,
    outer_bcs_type=outer_bcs_type,
    enable_rfm_precompute=enable_rfm_precompute,
    nrpyelliptic_project=True,
    post_non_y_n_auxevol_mallocs="auxevol_gfs_set_to_constant(&commondata, &griddata_chare[grid].params, griddata_chare[grid].xx, &griddata_chare[grid].gridfuncs);\n",
    post_MoL_step_forward_in_time=post_MoL_step_forward_in_time,
)

BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

BHaH.BHaH_defines_h.output_BHaH_defines_h(
    additional_includes=[str(Path("superB") / Path("superB.h"))],
    project_dir=project_dir,
    enable_rfm_precompute=True,
)

if enable_simd:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

superB.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    compiler_opt_option="default",
    addl_CFLAGS=["-fpermissive "],
    addl_libraries=["-module CkIO"],
    CC="charmc",
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./charmrun +p4 ./{project_name} to run with 4 processors, for example."
)
print(f"    Parameter file can be found in {project_name}.par")

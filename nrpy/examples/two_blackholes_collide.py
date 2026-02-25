"""
Two black holes collide example.

Specifically, evolve Brill-Lindquist initial data forward in time;
this example sets up a complete C code for solving the GR field
  equations in curvilinear coordinates on a cell-centered grid,
  using a reference metric approach.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import argparse
import os
import shutil
import subprocess

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

parser = argparse.ArgumentParser(
    description="NRPyElliptic Solver for Conformally Flat BBH initial data"
)
parser.add_argument(
    "--cuda",
    action="store_true",
    help="Use CUDA parallelization.",
)
parser.add_argument(
    "--floating_point_precision",
    type=str,
    help="Floating point precision (e.g. float, double).",
    default="double",
)
args = parser.parse_args()

# Code-generation-time parameters:
fp_type = args.floating_point_precision.lower()
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
project_name = "two_blackholes_collide"
CoordSystem = "GeneralRFM_fisheyeN1"
IDtype = "BrillLindquist"
IDCoordSystem = "Cartesian"
num_fisheye_transitions = (
    int(CoordSystem.replace("GeneralRFM_fisheyeN", ""))
    if CoordSystem.startswith("GeneralRFM_fisheyeN")
    else None
)
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
    "GeneralRFM_fisheyeN1": [128, 128, 128],
}
default_BH1_mass = default_BH2_mass = 0.5
default_BH1_z_posn = +0.5
default_BH2_z_posn = -0.5
# Fisheye parameters
fisheye_param_defaults: dict[str, float] = {}
if num_fisheye_transitions is not None:
    for i in range(num_fisheye_transitions + 1):
        fisheye_param_defaults[f"fisheye_a{i}"] = float(2**i)
    fisheye_param_defaults["fisheye_phys_L"] = grid_physical_size
if num_fisheye_transitions == 1:
    fisheye_param_defaults.update(
        {
            "fisheye_a0": 1.0,
            "fisheye_a1": 2.0,
            "fisheye_phys_L": grid_physical_size,
            "fisheye_phys_r_trans1": 8.0,
            "fisheye_phys_w_trans1": 2.0,
        }
    )
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
separate_Ricci_and_BSSN_RHS = True
enable_parallel_codegen = True
enable_rfm_precompute = True  # WIP: Will remove; for ease of maintenance we are no longer supporting disabled
enable_intrinsics = True  # WIP: Will remove; for ease of maintenance we are no longer supporting disabled
enable_fd_functions = True
enable_KreissOliger_dissipation = False
enable_CAKO = True
boundary_conditions_desc = "outgoing radiation"

set_of_CoordSystems = {CoordSystem}
NUMGRIDS = len(set_of_CoordSystems)
num_cuda_streams = NUMGRIDS
par.adjust_CodeParam_default("NUMGRIDS", NUMGRIDS)
enable_bhahaha = parallelization == "openmp"

BHaHAHA_subdir = "BHaHAHA"
if fd_order != 6:
    BHaHAHA_subdir = f"BHaHAHA-{fd_order}o"

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

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
if enable_bhahaha:
    try:
        # Attempt to run as a script path
        subprocess.run(
            [
                "python",
                "nrpy/examples/bhahaha.py",
                "--fdorder",
                str(fd_order),
                "--outrootdir",
                project_dir,
            ],
            check=True,
        )
    except subprocess.CalledProcessError:
        # If it fails (e.g., from a pip install), try running as a module
        subprocess.run(
            [
                "python",
                "-m",
                "nrpy.examples.bhahaha",
                "--fdorder",
                str(fd_order),
                "--outrootdir",
                project_dir,
            ],
            check=True,
        )
    BHaH.BHaHAHA.BHaH_implementation.register_CFunction_bhahaha_find_horizons(
        CoordSystem=CoordSystem, max_horizons=3
    )
    BHaH.BHaHAHA.interpolation_3d_general__uniform_src_grid.register_CFunction_interpolation_3d_general__uniform_src_grid(
        enable_simd=enable_intrinsics, project_dir=project_dir
    )

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

BHaH.general_relativity.initial_data.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    ID_persist_struct_str="",
)

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size for _ in set_of_CoordSystems],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
BHaH.diagnostics.diagnostics.register_all_diagnostics(
    project_dir=project_dir,
    set_of_CoordSystems=set_of_CoordSystems,
    default_diagnostics_out_every=diagnostics_output_every,
    enable_nearest_diagnostics=True,
    enable_interp_diagnostics=False,
    enable_volume_integration_diagnostics=True,
    enable_free_auxevol=False,
    enable_psi4_diagnostics=False,
    enable_bhahaha=enable_bhahaha,
)
BHaH.general_relativity.diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics=False, enable_psi4=False
)
BHaH.general_relativity.diagnostics_nearest.register_CFunction_diagnostics_nearest()
BHaH.general_relativity.diagnostics_volume_integration.register_CFunction_diagnostics_volume_integration()
if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems=set_of_CoordSystems,
    )
BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=False,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    enable_CAKO=enable_CAKO,
    OMP_collapse=OMP_collapse,
)
if separate_Ricci_and_BSSN_RHS:
    BHaH.general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
        CoordSystem=CoordSystem,
        enable_intrinsics=enable_intrinsics,
        enable_fd_functions=enable_fd_functions,
        OMP_collapse=OMP_collapse,
    )
    if parallelization == "cuda":
        BHaH.general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
            CoordSystem=CoordSystem,
            enable_intrinsics=enable_intrinsics,
            enable_fd_functions=enable_fd_functions,
            OMP_collapse=OMP_collapse,
            host_only_version=True,
        )
BHaH.general_relativity.enforce_detgammabar_equals_detgammahat.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BHaH.general_relativity.constraints_eval.register_CFunction_constraints_eval(
    CoordSystem=CoordSystem,
    enable_T4munu=False,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

if __name__ == "__main__":
    pcg.do_parallel_codegen()

if num_fisheye_transitions is not None:
    BHaH.fisheye.phys_params_to_fisheye.register_CFunction_fisheye_params_from_physical_N(
        num_transitions=num_fisheye_transitions
    )

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
)
rhs_string = ""
if separate_Ricci_and_BSSN_RHS:
    rhs_string += "Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);"
rhs_string += """
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
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
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(params, rfmstruct, RK_OUTPUT_GFS, auxevol_gfs);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
    rational_const_alias=(
        "static constexpr" if parallelization == "cuda" else "static const"
    ),
)
BHaH.xx_tofrom_Cart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
par.adjust_CodeParam_default("t_final", t_final)
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)
par.adjust_CodeParam_default("BH1_mass", default_BH1_mass)
par.adjust_CodeParam_default("BH2_mass", default_BH2_mass)
par.adjust_CodeParam_default("BH1_posn_z", default_BH1_z_posn)
par.adjust_CodeParam_default("BH2_posn_z", default_BH2_z_posn)
if num_fisheye_transitions is not None:
    for parname, value in fisheye_param_defaults.items():
        par.adjust_CodeParam_default(parname, value)
if enable_bhahaha:
    # Set BHaHAHA defaults to reasonable values.
    par.adjust_CodeParam_default(
        "bah_initial_grid_z_center", [default_BH1_z_posn, default_BH2_z_posn, 0.0]
    )
    par.adjust_CodeParam_default("bah_Nr_interp_max", 40)
    par.adjust_CodeParam_default(
        "bah_M_scale",
        [default_BH1_mass, default_BH2_mass, default_BH1_mass + default_BH2_mass],
    )
    par.adjust_CodeParam_default(
        "bah_max_search_radius",
        [
            0.6 * default_BH1_mass,
            0.6 * default_BH2_mass,
            1.9 * (default_BH1_mass + default_BH2_mass),
        ],
    )
    par.adjust_CodeParam_default(
        "bah_Theta_L2_times_M_tolerance",
        [2e-4, 2e-4, 2e-4],
    )
    par.adjust_CodeParam_default(
        "bah_Nphi_array_multigrid",
        [2, 2, 2],
    )
    par.adjust_CodeParam_default(
        "bah_Nr_interp_max",
        40,
    )
    par.adjust_CodeParam_default(
        "bah_enable_BBH_mode",
        1,
    )

BHaH.diagnostics.diagnostic_gfs_h_create.diagnostics_gfs_h_create(
    project_dir=project_dir
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
    additional_includes=(
        [os.path.join(BHaHAHA_subdir, "BHaHAHA.h")] if enable_bhahaha else []
    ),
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
    DOUBLE_means="double" if fp_type == "float" else "REAL",
    restrict_pointer_type="*" if parallelization == "cuda" else "*restrict",
)

# Set griddata struct used for calculations to griddata_device for certain parallelizations
compute_griddata = "griddata_device" if parallelization == "cuda" else "griddata"
post_params_struct_set_to_default = ""
if num_fisheye_transitions is not None:
    post_params_struct_set_to_default = BHaH.fisheye.phys_params_to_fisheye.build_post_params_struct_set_to_default_hook(
        num_transitions=num_fisheye_transitions,
        compute_griddata=compute_griddata,
    )

BHaH.main_c.register_CFunction_main_c(
    initial_data_desc=IDtype,
    MoL_method=MoL_method,
    boundary_conditions_desc=boundary_conditions_desc,
    post_params_struct_set_to_default=post_params_struct_set_to_default,
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    enable_bhahaha=enable_bhahaha,
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
    addl_dirs_to_make=([BHaHAHA_subdir] if enable_bhahaha else []),
    CC=("nvcc" if parallelization == "cuda" else "autodetect"),
    src_code_file_ext=("cu" if parallelization == "cuda" else "c"),
    addl_libraries=(
        [f"-L{BHaHAHA_subdir}", f"-l{BHaHAHA_subdir}"] if enable_bhahaha else []
    ),
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

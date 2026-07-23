"""
Two black holes collide example.

Specifically, evolve Brill-Lindquist initial data forward in time;
this example sets up a complete C code for solving the GR field
  equations in curvilinear coordinates on a cell-centered grid,
  using a reference metric approach.

Optionally, this script can generate metric spacetime data for raytracing
workflows and write a helper script that combines raytracing time slices.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1: Import needed Python modules, then set codegen and compile-time
# parameters.
import argparse
import json
import os
import shutil
import subprocess
from typing import Dict

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

SUPPORTED_RAYTRACING_COORD_SYSTEMS = ("SinhCylindricalv2n2",)
DEFAULT_RAYTRACING_NXX = [162, 2, 256]
DEFAULT_RAYTRACING_DOMAIN = [30.0, 0.075, 0.05, 1.0, 4.0]
DEFAULT_RAYTRACING_T_FINAL = 30.0
DEFAULT_RAYTRACING_DIAGNOSTICS_OUTPUT_EVERY = 0.4

parser = argparse.ArgumentParser(
    description="""Generate a BHaH two-black-hole collision example with optional raytracing spacetime outputs."""
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
parser.add_argument(
    "--raytracing-time",
    nargs="*",
    default=None,
    type=float,
    metavar="VALUE",
    help="""Enable metric raytracing spacetime outputs. Provide T_FINAL and DIAGNOSTICS_OUTPUT_EVERY, or omit both to use the tuned 30.0/0.4 defaults. Currently supported only for OpenMP builds with double precision.""",
)
parser.add_argument(
    "--raytracing-data-mode",
    choices=("g4DD", "g4DD_d0", "GammaUDD", "all"),
    default="g4DD",
    help="""Raytracing data written by the spacetime evolution. 'all' writes all three method-specific datasets in one evolution.""",
)
parser.add_argument(
    "--raytracing-static-christoffels",
    action="store_true",
    help="""Use static-spacetime Christoffels for the output slice nearest to T_FINAL when the selected mode includes GammaUDD.""",
)
parser.add_argument(
    "--raytracing-domain",
    nargs="+",
    default=None,
    type=float,
    help="""Override the SinhCylindricalv2n2 domain parameters GRID_PHYSICAL_SIZE SINHWRHO SINHWZ RHO_SLOPE Z_SLOPE. Omit to use 30.0 0.075 0.05 1.0 4.0.""",
)
parser.add_argument(
    "--raytracing-coord-system",
    type=str,
    choices=SUPPORTED_RAYTRACING_COORD_SYSTEMS,
    default=None,
    help="""Override the generated project's raytracing coordinate system when --raytracing-time is enabled.""",
)
parser.add_argument(
    "--raytracing-Nxx",
    dest="raytracing_nxx",
    nargs=3,
    type=int,
    default=None,
    metavar=("NXX0", "NXX1", "NXX2"),
    help="""Override the tuned SinhCylindricalv2n2 base Nxx, e.g. --raytracing-Nxx 162 2 256. The reduced angular direction is NXX1 == 2.""",
)
parser.add_argument(
    "--raytracing-bhs",
    nargs=4,
    type=float,
    default=None,
    metavar=("Z_1", "Z_2", "M_1", "M_2"),
    help="""Override the generated project's black-hole z positions and masses, e.g. --raytracing-bhs 0.5 -0.5 0.5 0.5.""",
)
args = parser.parse_args()


# Step 1.a: Validate command-line inputs and select build settings.
fp_type = args.floating_point_precision.lower()
# Step 1.b: Default to OpenMP; override with CUDA if --cuda is set.
parallelization = "cuda" if args.cuda else "openmp"
if parallelization not in ["openmp", "cuda"]:
    raise ValueError(
        f"""Invalid parallelization strategy: {parallelization}. Choose 'openmp' or 'cuda'."""
    )
# Step 1.c: Raytracing spacetime output requires OpenMP, double precision,
# and consistent raytracing-specific option combinations.
enable_raytracing_data_output = args.raytracing_time is not None
raytracing_data_mode = args.raytracing_data_mode

if fp_type not in ("float", "double"):
    raise ValueError("--floating_point_precision must be either 'float' or 'double'.")
if enable_raytracing_data_output and fp_type != "double":
    raise ValueError(
        "--raytracing-time currently requires --floating_point_precision double."
    )

if args.cuda and enable_raytracing_data_output:
    raise ValueError("--raytracing-time is currently supported only for OpenMP builds.")
if not enable_raytracing_data_output and args.raytracing_data_mode != "g4DD":
    raise ValueError("--raytracing-data-mode requires --raytracing-time.")
if not enable_raytracing_data_output and args.raytracing_static_christoffels:
    raise ValueError("--raytracing-static-christoffels requires --raytracing-time.")
if args.raytracing_static_christoffels and raytracing_data_mode not in (
    "GammaUDD",
    "all",
):
    raise ValueError(
        "--raytracing-static-christoffels requires "
        "--raytracing-data-mode GammaUDD or all."
    )
if args.raytracing_coord_system is not None and not enable_raytracing_data_output:
    raise ValueError("--raytracing-coord-system requires --raytracing-time.")
if args.raytracing_domain is not None and not enable_raytracing_data_output:
    raise ValueError("--raytracing-domain requires --raytracing-time.")
if args.raytracing_nxx is not None and not enable_raytracing_data_output:
    raise ValueError("--raytracing-Nxx requires --raytracing-time.")
par.set_parval_from_str("Infrastructure", "BHaH")
par.set_parval_from_str("parallelization", parallelization)
par.set_parval_from_str("fp_type", fp_type)

# Step 1.b: Set code-generation-time parameters.
project_name = "two_blackholes_collide"
CoordSystem = (
    args.raytracing_coord_system
    if args.raytracing_coord_system is not None
    else ("SinhCylindricalv2n2" if enable_raytracing_data_output else "Spherical")
)
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
sinh_width = None
sinh_width_rho = None
sinh_width_z = None
rho_slope = None
z_slope = None
if enable_raytracing_data_output:
    if len(args.raytracing_time) == 0:
        t_final = DEFAULT_RAYTRACING_T_FINAL
        diagnostics_output_every = DEFAULT_RAYTRACING_DIAGNOSTICS_OUTPUT_EVERY
    elif len(args.raytracing_time) == 2:
        t_final = args.raytracing_time[0]
        diagnostics_output_every = args.raytracing_time[1]
    else:
        raise ValueError(
            "--raytracing-time accepts either no values or "
            "T_FINAL DIAGNOSTICS_OUTPUT_EVERY."
        )
    domain = (
        list(args.raytracing_domain)
        if args.raytracing_domain is not None
        else list(DEFAULT_RAYTRACING_DOMAIN)
    )

    if t_final <= 0.0:
        raise ValueError("--raytracing-time T_FINAL must be positive.")
    if diagnostics_output_every <= 0.0:
        raise ValueError("--raytracing-time DIAGNOSTICS_OUTPUT_EVERY must be positive.")

    if len(domain) != 5:
        raise ValueError(
            """SinhCylindricalv2n2 expects --raytracing-domain GRID_PHYSICAL_SIZE SINHWRHO SINHWZ RHO_SLOPE Z_SLOPE."""
        )
    (
        grid_physical_size,
        sinh_width_rho,
        sinh_width_z,
        rho_slope,
        z_slope,
    ) = domain

    if grid_physical_size <= 0.0:
        raise ValueError("--raytracing-domain GRID_PHYSICAL_SIZE must be positive.")
    if sinh_width_rho is not None and sinh_width_rho <= 0.0:
        raise ValueError("--raytracing-domain SINHWRHO must be positive.")
    if sinh_width_z is not None and sinh_width_z <= 0.0:
        raise ValueError("--raytracing-domain SINHWZ must be positive.")
    if rho_slope is not None and rho_slope <= 0.0:
        raise ValueError("--raytracing-domain RHO_SLOPE must be positive.")
    if z_slope is not None and z_slope <= 0.0:
        raise ValueError("--raytracing-domain Z_SLOPE must be positive.")
else:
    t_final = 7.5
    grid_physical_size = 7.5
    diagnostics_output_every = 0.25

Nxx_dict = {
    "Spherical": [72, 12, 2],
    "SinhSpherical": [72, 12, 2],
    "SinhCylindricalv2n2": list(DEFAULT_RAYTRACING_NXX),
    "Cartesian": [64, 64, 64],
    "GeneralRFM_fisheyeN1": [128, 128, 128],
    "GeneralRFM_fisheyeN2": [128, 128, 128],
}
if CoordSystem not in Nxx_dict:
    raise ValueError(
        f"Unsupported coordinate system '{CoordSystem}'. "
        f"Choose one of {sorted(Nxx_dict)}."
    )
if enable_raytracing_data_output:
    nxx_override = (
        list(args.raytracing_nxx)
        if args.raytracing_nxx is not None
        else list(DEFAULT_RAYTRACING_NXX)
    )
    if any(nxx_value <= 0 for nxx_value in nxx_override):
        raise ValueError("--raytracing-Nxx values must all be positive.")
    Nxx_dict[CoordSystem] = nxx_override
if args.raytracing_bhs is None:
    default_BH1_z_posn = +0.5
    default_BH2_z_posn = -0.5
    default_BH1_mass = 0.5
    default_BH2_mass = 0.5
else:
    default_BH1_z_posn = args.raytracing_bhs[0]
    default_BH2_z_posn = args.raytracing_bhs[1]
    default_BH1_mass = args.raytracing_bhs[2]
    default_BH2_mass = args.raytracing_bhs[3]
    if default_BH1_mass <= 0.0 or default_BH2_mass <= 0.0:
        raise ValueError("--raytracing-bhs masses M_1 and M_2 must be positive.")
# Step 1.d: Derive fisheye parameter defaults from grid_physical_size when a
# GeneralRFM_fisheyeN* coordinate system is selected.
fisheye_param_defaults: Dict[str, float] = {}
if num_fisheye_transitions == 1:
    fisheye_param_defaults = {
        "fisheye_phys_a0": 1.0,
        "fisheye_phys_a1": 1.5,
        "fisheye_phys_L": grid_physical_size,
        "fisheye_phys_r_trans1": 2.0,
        "fisheye_phys_w_trans1": 1.0,
    }
elif num_fisheye_transitions == 2:
    fisheye_param_defaults = {
        "fisheye_phys_a0": 1.0,
        "fisheye_phys_a1": 2.0,
        "fisheye_phys_a2": 4.0,
        "fisheye_phys_L": grid_physical_size,
        "fisheye_phys_r_trans1": 1.5,
        "fisheye_phys_w_trans1": 0.8,
        "fisheye_phys_r_trans2": 4.5,
        "fisheye_phys_w_trans2": 1.5,
    }
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
separate_Ricci_and_BSSN_RHS = True
enable_parallel_codegen = True
enable_rfm_precompute = True
enable_intrinsics = True
enable_fd_functions = True
enable_KreissOliger_dissipation = False
enable_CAKO = True
boundary_conditions_desc = "outgoing radiation"

set_of_CoordSystems = {CoordSystem}
num_cuda_streams = 1
enable_bhahaha = parallelization == "openmp"
if enable_bhahaha and fp_type != "double":
    raise ValueError(
        "BHaHAHA integration currently requires --floating_point_precision double."
    )

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

# Step 1.e: Clean the project directory before regenerating the example.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

# Step 2: Declare core C functions and register each to
# cfc.CFunction_dict["function_name"].
if enable_bhahaha:
    try:
        # Step 2.a: Attempt to run the BHaHAHA generator as a script path.
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
        # Step 2.b: Fall back to the module entry point when the script path
        # is unavailable.
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
        max_horizons=3
    )
    BHaH.BHaHAHA.interpolation_3d_general__uniform_src_grid.register_CFunction_interpolation_3d_general__uniform_src_grid(
        enable_simd=enable_intrinsics, project_dir=project_dir
    )

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

BHaH.general_relativity.initial_data.register_CFunction_initial_data(
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    set_of_CoordSystems=set_of_CoordSystems,
    ID_persist_struct_str="",
)

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size for _ in set_of_CoordSystems],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
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
BHaH.diagnostics.diagnostics.register_all_diagnostics(
    project_dir=project_dir,
    set_of_CoordSystems=set_of_CoordSystems,
    default_diagnostics_out_every=diagnostics_output_every,
    enable_nearest_diagnostics=True,
    enable_interp_diagnostics=False,
    enable_volume_integration_diagnostics=True,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_free_auxevol=False,
    enable_psi4_diagnostics=False,
    enable_bhahaha=enable_bhahaha,
    enable_raytracing_data_output=enable_raytracing_data_output,
    raytracing_data_mode=raytracing_data_mode,
    enable_static_christoffels=args.raytracing_static_christoffels,
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
BHaH.xx_tofrom_Cart.register_CFunction_Cart_to_xx_and_nearest_i0i1i2_assume_valid(
    CoordSystem
)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.general_relativity.basis_transforms.register_all.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
)
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

# Step 3: Generate header files, register C functions and command-line
# parameters, set up boundary conditions, and create a Makefile in
# project/[project_name]/.
par.adjust_CodeParam_default("t_final", t_final)
if CoordSystem == "SinhSpherical":
    assert sinh_width is not None
    par.adjust_CodeParam_default("SINHW", sinh_width)
if CoordSystem == "SinhCylindrical":
    assert sinh_width_rho is not None
    assert sinh_width_z is not None
    par.adjust_CodeParam_default("SINHWRHO", sinh_width_rho)
    par.adjust_CodeParam_default("SINHWZ", sinh_width_z)
if CoordSystem == "SinhCylindricalv2n2":
    assert sinh_width_rho is not None
    assert sinh_width_z is not None
    assert rho_slope is not None
    assert z_slope is not None
    par.adjust_CodeParam_default("SINHWRHO", sinh_width_rho)
    par.adjust_CodeParam_default("SINHWZ", sinh_width_z)
    par.adjust_CodeParam_default("rho_slope", rho_slope)
    par.adjust_CodeParam_default("z_slope", z_slope)
par.adjust_CodeParam_default("eta", GammaDriving_eta)
par.adjust_CodeParam_default("BH1_mass", default_BH1_mass)
par.adjust_CodeParam_default("BH2_mass", default_BH2_mass)
par.adjust_CodeParam_default("BH1_posn_z", default_BH1_z_posn)
par.adjust_CodeParam_default("BH2_posn_z", default_BH2_z_posn)
if num_fisheye_transitions is not None:
    for parname, value in fisheye_param_defaults.items():
        par.adjust_CodeParam_default(parname, value)
if enable_bhahaha:
    # Step 3.a: Set BHaHAHA defaults to reasonable values.
    par.adjust_CodeParam_default(
        "bah_initial_grid_z_center", [default_BH1_z_posn, default_BH2_z_posn, 0.0]
    )
    par.adjust_CodeParam_default("bah_Nr_interp_max", 150)
    par.adjust_CodeParam_default(
        "bah_M_scale",
        [default_BH1_mass, default_BH2_mass, default_BH1_mass + default_BH2_mass],
    )
    par.adjust_CodeParam_default(
        "bah_max_search_radius",
        [
            1.0 * default_BH1_mass,
            1.0 * default_BH2_mass,
            2.4 * (default_BH1_mass + default_BH2_mass),
        ],
    )
    par.adjust_CodeParam_default(
        "bah_Theta_L2_times_M_tolerance",
        [2e-4, 2e-4, 2e-4],
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

# Step 3.b: Select the griddata handle for the active parallelization mode.
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

# Step 3.c: Copy SIMD intrinsics needed for interpolation and diagnostics.
intrinsics_file_list = ["simd_intrinsics.h"]
if parallelization == "cuda":
    # Step 3.c.i: Copy CUDA intrinsics for CUDA-enabled projects.
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

if enable_raytracing_data_output:
    # Step 3.d: Copy the raytracing time-slice combiner into the generated
    # project.
    copy_files(
        package="nrpy.infrastructures.BHaH.diagnostics",
        filenames_list=["combine_raytracing_time_slices.py"],
        project_dir=project_dir,
        subdirectory="",
    )

    # Step 3.d.i: Build a reproducible combined-output filename from the
    # raytracing parameters.
    def _format_output_name_value(output_value: float) -> str:
        return str(output_value).replace("-", "neg").replace(".", "p")

    t_final_name = _format_output_name_value(t_final)
    grid_physical_size_name = _format_output_name_value(grid_physical_size)
    diagnostics_output_every_name = _format_output_name_value(diagnostics_output_every)
    bh1_z_name = _format_output_name_value(default_BH1_z_posn)
    bh2_z_name = _format_output_name_value(default_BH2_z_posn)
    bh1_mass_name = _format_output_name_value(default_BH1_mass)
    bh2_mass_name = _format_output_name_value(default_BH2_mass)
    nxx0, nxx1, nxx2 = Nxx_dict[CoordSystem]

    assert sinh_width_rho is not None
    assert sinh_width_z is not None
    assert rho_slope is not None
    assert z_slope is not None
    sinh_width_rho_name = _format_output_name_value(sinh_width_rho)
    sinh_width_z_name = _format_output_name_value(sinh_width_z)
    rho_slope_name = _format_output_name_value(rho_slope)
    z_slope_name = _format_output_name_value(z_slope)
    domain_suffix = (
        f"sinhwrho_{sinh_width_rho_name}_sinhwz_{sinh_width_z_name}_"
        f"rho-slope_{rho_slope_name}_z-slope_{z_slope_name}_"
    )

    raytracing_combined_bin_stem = f"{project_name}_{t_final_name}_{grid_physical_size_name}_{diagnostics_output_every_name}_z1_{bh1_z_name}_z2_{bh2_z_name}_M1_{bh1_mass_name}_M2_{bh2_mass_name}_{CoordSystem}_{domain_suffix}{nxx0}_{nxx1}_{nxx2}"
    raytracing_modes = (
        ("g4DD", "g4DD_d0", "GammaUDD")
        if raytracing_data_mode == "all"
        else (raytracing_data_mode,)
    )
    raytracing_combined_bin_names = {
        mode: f"{raytracing_combined_bin_stem}_{mode}.bin" for mode in raytracing_modes
    }
    combined_output_paths = {
        mode: os.path.join("..", "raytracing_data", filename)
        for mode, filename in raytracing_combined_bin_names.items()
    }
    raytracing_run_metadata_name = "raytracing_run_metadata.json"
    raytracing_run_metadata_path = os.path.join(
        project_dir, raytracing_run_metadata_name
    )
    raytracing_run_metadata = {
        "black_holes": {
            "BH1_mass": default_BH1_mass,
            "BH1_posn_z": default_BH1_z_posn,
            "BH2_mass": default_BH2_mass,
            "BH2_posn_z": default_BH2_z_posn,
        },
        "raytracing_data_mode": raytracing_data_mode,
        "datasets": raytracing_combined_bin_names,
        "raytracing_static_christoffels": args.raytracing_static_christoffels,
        "generator_script": "two_blackholes_collide.py",
        "grid_physical_size": grid_physical_size,
        "project_name": project_name,
        "raytracing_Nxx": [nxx0, nxx1, nxx2],
        "raytracing_coord_system": CoordSystem,
        "raytracing_domain": domain,
        "raytracing_time": {
            "diagnostics_output_every": diagnostics_output_every,
            "t_final": t_final,
        },
        "slice_filename_pattern": "raytracing_data_t????????.bin",
    }
    if raytracing_data_mode == "all":
        raytracing_run_metadata["generated_datasets"] = list(raytracing_modes)
    else:
        raytracing_run_metadata["combined_output_filename"] = (
            raytracing_combined_bin_names[raytracing_data_mode]
        )
    with open(raytracing_run_metadata_path, "w", encoding="utf-8") as metadata_file:
        json.dump(
            raytracing_run_metadata,
            metadata_file,
            sort_keys=True,
            indent=2,
        )
        metadata_file.write("\n")

    # Step 3.d.ii: Generate a helper script that builds, runs, and combines
    # time slices.
    script_path = os.path.join(project_dir, "run_raytracing_data_pipeline.sh")
    combine_commands = "\n".join(f"""echo "Combining {mode} raytracing time slices..."
python3 combine_raytracing_time_slices.py \\
  --input-dir "raytracing_slices/{mode}" \\
  --pattern "raytracing_data_t????????.bin" \\
  --run-metadata "{raytracing_run_metadata_name}" \\
  --output "{combined_output_paths[mode]}" \\
  --force""" for mode in raytracing_modes)
    combined_output_messages = "\n".join(
        f'echo "{combined_output_paths[mode]}"' for mode in raytracing_modes
    )
    script_text = rf"""#!/usr/bin/env bash
set -euo pipefail

echo "Building {project_name}..."
make

echo "Running ./{project_name}..."
./{project_name}

{combine_commands}

echo "Combined raytracing data written to:"
{combined_output_messages}
"""
    with open(script_path, "w", encoding="utf-8") as script_file:
        script_file.write(script_text)
    os.chmod(script_path, 0o755)

    # Step 3.d.iii: Report the generated run instructions and output paths.
    print(
        f"""Finished! Now go into project/{project_name} and run `./run_raytracing_data_pipeline.sh`."""
    )
    for mode in raytracing_modes:
        print(
            f"    Single-time-slice {mode} files will remain in "
            f"project/{project_name}/raytracing_slices/{mode}/."
        )
        print(
            f"    Combined {mode} output will be written to "
            f"project/raytracing_data/{raytracing_combined_bin_names[mode]}"
        )
else:
    print(
        f"""Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."""
    )
print(f"    Parameter file can be found in {project_name}.par")

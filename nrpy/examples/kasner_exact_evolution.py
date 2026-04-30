"""
Kasner exact-solution example using the BHaH infrastructure.

This generates a complete C project for evolving the vacuum Kasner spacetime
with BSSN.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* com
"""

#########################################################
# Step 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import argparse
import os
import shutil
from typing import Dict

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

parser = argparse.ArgumentParser(description="Generate a BSSN Kasner benchmark project")
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
project_name = "kasner_exact_evolution"
CoordSystem = "GeneralRFM_fisheyeN1"
IDtype = "Kasner"
IDCoordSystem = "Cartesian"
num_fisheye_transitions = (
    int(CoordSystem.replace("GeneralRFM_fisheyeN", ""))
    if CoordSystem.startswith("GeneralRFM_fisheyeN")
    else None
)
LapseEvolutionOption = "Frozen"
ShiftEvolutionOption = "Frozen"
grid_physical_size = 12.0
diagnostics_output_every = 0.05
t_final = 7.0
Nxx_dict = {
    "Cartesian": [32, 32, 32],
    "GeneralRFM_fisheyeN1": [128, 128, 128],
}
# Fisheye parameters
fisheye_param_defaults: Dict[str, float] = {}
if num_fisheye_transitions == 1:
    fisheye_param_defaults = {
        "fisheye_phys_a0": 1.0,
        "fisheye_phys_a1": 2.0,
        "fisheye_phys_L": grid_physical_size,
        "fisheye_phys_r_trans1": 3.0,
        "fisheye_phys_w_trans1": 1.0,
    }
elif num_fisheye_transitions == 2:
    fisheye_param_defaults = {
        "fisheye_phys_a0": 1.0,
        "fisheye_phys_a1": 2.0,
        "fisheye_phys_a2": 3.0,
        "fisheye_phys_L": grid_physical_size,
        "fisheye_phys_r_trans1": 3.0,
        "fisheye_phys_w_trans1": 1.0,
        "fisheye_phys_r_trans2": 5.0,
        "fisheye_phys_w_trans2": 1.0,
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
enable_CAKO = False
boundary_conditions_desc = "extrapolation"
outer_bcs_type = "extrapolation"

set_of_CoordSystems = {CoordSystem}
OMP_collapse = 1
if CoordSystem not in Nxx_dict:
    raise ValueError(f"CoordSystem = {CoordSystem} not supported by Nxx_dict.")
par.adjust_CodeParam_default("NUMGRIDS", len(set_of_CoordSystems))

project_dir = os.path.join("project", project_name)
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.set_parval_from_str("EvolvedConformalFactor_cf", "W")

if parallelization == "cuda":
    BHaH.parallelization.cuda_utilities.register_CFunctions_HostDevice__operations()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_minimum()
    BHaH.parallelization.cuda_utilities.register_CFunction_find_global_sum()

if num_fisheye_transitions is not None:
    BHaH.fisheye.phys_params_to_fisheye.register_CFunction_fisheye_params_from_physical_N(
        num_transitions=num_fisheye_transitions
    )
    for parname, value in fisheye_param_defaults.items():
        par.adjust_CodeParam_default(parname, value)

kasner_param_guard = r"""
const REAL kasner_sum_p = commondata->KASNER_p1 + commondata->KASNER_p2 + commondata->KASNER_p3;
const REAL kasner_sum_p2 = commondata->KASNER_p1 * commondata->KASNER_p1
                         + commondata->KASNER_p2 * commondata->KASNER_p2
                         + commondata->KASNER_p3 * commondata->KASNER_p3;
const REAL kasner_constraint_tol =
    (sizeof(REAL) == sizeof(float)) ? 1.0e-6 : 1.0e-14;
const REAL kasner_abs_sum_err = (kasner_sum_p > 1.0) ? (kasner_sum_p - 1.0) : (1.0 - kasner_sum_p);
const REAL kasner_abs_sum2_err = (kasner_sum_p2 > 1.0) ? (kasner_sum_p2 - 1.0) : (1.0 - kasner_sum_p2);
if (kasner_abs_sum_err > kasner_constraint_tol || kasner_abs_sum2_err > kasner_constraint_tol) {
  fprintf(stderr, "Error: KASNER_p1, KASNER_p2, KASNER_p3 must satisfy p1+p2+p3=1 and p1^2+p2^2+p3^2=1.\\n");
  exit(1);
}
"""

BHaH.general_relativity.initial_data.register_CFunction_initial_data(
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    set_of_CoordSystems=set_of_CoordSystems,
    ID_persist_struct_str="",
    populate_ID_persist_struct_str=kasner_param_guard,
)

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size],
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
)
BHaH.general_relativity.Kasner.diagnostics.register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics=False,
    enable_psi4=False,
)
BHaH.general_relativity.Kasner.diagnostics.register_CFunction_diagnostics_nearest()
BHaH.general_relativity.diagnostics_volume_integration.register_CFunction_diagnostics_volume_integration()
if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems=set_of_CoordSystems,
    )
# Use separate Ricci only where that path is actually supported.
use_separate_ricci = separate_Ricci_and_BSSN_RHS and not (
    parallelization == "cuda" and CoordSystem.startswith("GeneralRFM")
)

BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=use_separate_ricci,
    enable_T4munu=False,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    enable_CAKO=enable_CAKO,
    OMP_collapse=OMP_collapse,
)
if use_separate_ricci:
    # GeneralRFM + CUDA does not support device-side Ricci_eval kernels.
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

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
)

par.adjust_CodeParam_default("outer_bc_type", outer_bcs_type)
rhs_string = ""
if use_separate_ricci:
    rhs_string += "Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);"
rhs_string += """
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);
"""
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

par.adjust_CodeParam_default("t_final", t_final)
BHaH.xx_tofrom_Cart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

if num_fisheye_transitions is not None:
    for parname, value in fisheye_param_defaults.items():
        par.adjust_CodeParam_default(parname, value)

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
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=[],
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
    DOUBLE_means="double" if fp_type == "float" else "REAL",
    restrict_pointer_type="*" if parallelization == "cuda" else "*restrict",
)

# Set griddata struct used for calculations to griddata_device for certain parallelizations.
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
)

# SIMD intrinsics needed for 3D interpolation, constraints evaluation, etc.
intrinsics_file_list = ["simd_intrinsics.h"]
if parallelization == "cuda":
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
    addl_dirs_to_make=[],
    CC=("nvcc" if parallelization == "cuda" else "autodetect"),
    src_code_file_ext=("cu" if parallelization == "cuda" else "c"),
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

"""
Set up BHaHAHA: the BlackHoles@Home Apparent Horizon Algorithm, based on NRPyElliptic, by Thiago Assumpção.

Handles the setup of the BlackHoles@Home Apparent Horizon Algorithm (BHaHAHA) for numerical relativity simulations,
including the registration of core functions, boundary conditions, and the generation of Makefile for compilation.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
"""

import argparse
import copy
import os
import pkgutil
import shutil
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

parser = argparse.ArgumentParser(
    description="BHaHAHA, the BlackHoles@Home Apparent Horizon Algorithm"
)
parser.add_argument(
    "--fdorder",
    type=int,
    help="Finite-difference order; interpolation order is one less than this. Default=6. For precision horizon finding, set to 8.",
    default=6,
)
parser.add_argument(
    "--outrootdir",
    type=str,
    help="Output root directory. Default = project/",
    default="project",
)
parser.add_argument(
    "--cpp",
    action="store_true",
    help="Generate C++-compatible BHaHAHA",
)
parser.add_argument(
    "--no-openmp",
    action="store_true",
    help="Disable OpenMP flags",
)
args = parser.parse_args()
fd_order = args.fdorder
outrootdir = args.outrootdir
use_cpp = args.cpp
use_openmp = not args.no_openmp  # default: True if flag omitted

par.set_parval_from_str("Infrastructure", "BHaH")

#########################################################
# STEP 1: Code-generation-time parameters
project_name = "BHaHAHA"
if fd_order != 6:
    project_name = f"BHaHAHA-{fd_order}o"
# These parameters must be set by the calling code.
CoordSystem = "Spherical"

OMP_collapse = 1
enable_rfm_precompute = True
enable_fd_functions = False
MoL_method = "SSPRK33"
enable_simd = False  # Does not work, as the SIMD vector would be in the radial direction, where only a single point is sampled.
parallel_codegen_enable = True
boundary_conditions_desc = "inner boundaries only"

#########################################################
# STEP 2: Project directory setup
project_dir = os.path.join(outrootdir, project_name)

# Clean the project directory if it exists
shutil.rmtree(project_dir, ignore_errors=True)

#########################################################
# STEP 3: Define the damping parameter, min global wavespeed, and residual stop condition
par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

#########################################################
# STEP 4: Declare core C functions & register each to cfc.CFunction_dict["function_name"]
BHaH.BHaHAHA.find_horizon.register_CFunction_find_horizon()
BHaH.BHaHAHA.poisoning_set_inputs.register_CFunction_poisoning_set_inputs()
BHaH.BHaHAHA.poisoning_check_inputs.register_CFunction_poisoning_check_inputs()
BHaH.BHaHAHA.over_relaxation.register_CFunction_over_relaxation()
BHaH.BHaHAHA.radial_grid_cell_centered_set_up.register_CFunction_radial_grid_cell_centered_set_up()
BHaH.BHaHAHA.xyz_center_r_minmax.register_CFunction_xyz_center_r_minmax()
BHaH.BHaHAHA.quadratic_extrapolation.register_CFunction_quadratic_extrapolation()
BHaH.BHaHAHA.numgrid__external_input_set_up.register_CFunction_numgrid__external_input_set_up()
BHaH.BHaHAHA.numgrid__interp_src_set_up.register_CFunction_numgrid__interp_src_set_up()
BHaH.BHaHAHA.numgrid__evol_set_up.register_CFunction_numgrid__evol_set_up()

BHaH.BHaHAHA.initial_data.register_CFunction_initial_data()
BHaH.BHaHAHA.interpolation_1d_radial_spokes_on_3d_src_grid.register_CFunction_interpolation_1d_radial_spokes_on_3d_src_grid(
    enable_simd=enable_simd, project_dir=project_dir
)
BHaH.BHaHAHA.interpolation_2d_external_input_to_interp_src_grid.register_CFunction_interpolation_2d_external_input_to_interp_src_grid()
BHaH.BHaHAHA.interpolation_2d_general__uniform_src_grid.register_CFunction_interpolation_2d_general__uniform_src_grid(
    enable_simd=enable_simd, project_dir=project_dir
)
BHaH.BHaHAHA.hDD_dD_and_W_dD_in_interp_src_grid_interior.register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior()
BHaH.BHaHAHA.error_message.register_CFunction_error_message()

# Register numerical grid functions
BHaH.BHaHAHA.cfl_limited_timestep_based_on_h_equals_r.register_CFunction_cfl_limited_timestep_based_on_h_equals_r()

BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

BHaH.BHaHAHA.diagnostics.register_CFunction_diagnostics()
BHaH.BHaHAHA.diagnostics_area_centroid_and_Theta_norms.register_CFunction_diagnostics_area_centroid_and_Theta_norms()
BHaH.BHaHAHA.diagnostics_file_output.register_CFunction_diagnostics_file_output()
BHaH.BHaHAHA.diagnostics_integration_weights.register_CFunction_diagnostics_integration_weights()
BHaH.BHaHAHA.diagnostics_min_max_mean_radii_wrt_centroid.register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid()
BHaH.BHaHAHA.diagnostics_proper_circumferences.register_CFunction_diagnostics_proper_circumferences()

if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems={CoordSystem},
    )

BHaH.BHaHAHA.rhs_eval_KO_apply.register_CFunction_rhs_eval(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
)
BHaH.BHaHAHA.rhs_eval_KO_apply.register_CFunction_KO_apply(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
)

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

#########################################################
# STEP 5: Register C function to set up the boundary condition struct
BHaH.BHaHAHA.apply_bcs_inner_only.register_CFunction_apply_bcs_inner_only()
# cbc.register_CFunction_apply_bcs_inner_only()
BHaH.BHaHAHA.apply_bcs_r_maxmin_partial_r_hDD_upwinding.register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding(
    upwinding_fd_order=fd_order
)
BHaH.BHaHAHA.bcstruct_set_up.register_CFunction_bcstruct_set_up(CoordSystem="Spherical")

rhs_string = """
bah_rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
// Theta was just evaluated at Nxx1*Nxx2 gridpoints; update counter:
commondata->bhahaha_diagnostics->Theta_eval_points_counter += params->Nxx1 * params->Nxx2;
if(commondata->KO_diss_strength > 0.0)
  bah_KO_apply(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
BHaH.MoLtimestepping.register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="bah_apply_bcs_inner_only(commondata, params, bcstruct, RK_OUTPUT_GFS);",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)

#########################################################
# STEP 6: Remove all rfm wrapper functions and disable output to Spherical/ subdirectory.
# Create a shallow copy of the original dictionary to iterate over (to avoid modifying while iterating)
original_dict_items = list(cfc.CFunction_dict.items())
for name, CFunction_obj in original_dict_items:
    # Remove the suffix from the name
    newname = name.replace("__rfm__Spherical", "")
    if newname != name:
        # Deep copy the CFunction object to avoid modifying the original directly
        cfc.CFunction_dict[newname] = copy.deepcopy(CFunction_obj)
        # Modify the copied object's name and regenerate the function
        cfc.CFunction_dict[newname].name = newname
        cfc.CFunction_dict[newname].CoordSystem_for_wrapper_func = ""
        cfc.CFunction_dict[newname].subdirectory = os.path.join(".")
        _, _, _ = cfc.CFunction_dict[newname].generate_full_function()
        # Remove the old key from the dictionary
        del cfc.CFunction_dict[name]

#########################################################
# STEP 7: Update parameters needed for hyperbolic relaxation method
par.adjust_CodeParam_default(
    "t_final", 1e30
)  # Ensure that pseudo-time limit is never hit; instead we limit by iteration counts.

#########################################################
# STEP 8: Generate header files, register C functions, set up boundary conditions, and create a Makefile
BHaH.CodeParameters.write_CodeParameters_h_files(project_dir=project_dir)
BHaH.CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_simd,
    define_no_simd_UPWIND_ALG=False,
)

#########################################################
# STEP 9: Copy files and construct Makefile
if enable_simd:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=f"lib{project_name}",
    lib_function_prefix="bah_",
    create_lib=True,
    static_lib=True,
    use_openmp=use_openmp,
)

# Append latest error codes & error message function prototype to BHaHAHA.h
# Load the header file using pkgutil
data_bytes = pkgutil.get_data("nrpy.infrastructures.BHaH.BHaHAHA", "BHaHAHA_header.h")
if data_bytes is None:
    raise FileNotFoundError("BHaHAHA_header.h not found via pkgutil.get_data")

BHaHAHA_h = data_bytes.decode("utf-8")

# Append finite-difference ghostzones setting
BHaHAHA_h += f"""
//===============================================
// Set the number of (finite-difference) ghostzones in BHaHAHA
//===============================================
#define BHAHAHA_NGHOSTS {int(par.parval_from_str("finite_difference::fd_order") / 2)}
"""

# Append the error codes
BHaHAHA_h += """
//===============================================
// BHaHAHA error handling
//===============================================
// Error codes, set in error_message.py
typedef enum {
"""
for item in BHaH.BHaHAHA.error_message.error_code_msg_tuples_list:
    BHaHAHA_h += f"  {item[0]},\n"
BHaHAHA_h += "} bhahaha_error_codes;\n"

BHaHAHA_h += """
// Function: bah_error_message
// Interprets bah_find_horizon() error codes & returns a useful string.
const char *bah_error_message(const bhahaha_error_codes error_code);
//===============================================

#endif // BHAHAHA_HEADER_H
"""

if use_cpp:
    # C++ compatibility: extern "C" and restrict mapping
    cpp_compatibility_preamble = (
        "#ifdef __cplusplus\n"
        'extern "C" {\n'
        "#endif\n\n"
        "#define restrict __restrict__\n\n"
    )
    cpp_compatibility_epilogue = "\n#ifdef __cplusplus\n" "}\n" "#endif"
    BHaHAHA_h = cpp_compatibility_preamble + BHaHAHA_h + cpp_compatibility_epilogue

    # Convert fixed-size parameter to pointer
    BHaHAHA_h = BHaHAHA_h.replace("REAL radii[Nr_interp_max]", "REAL radii[]")


# Write the updated content to the output file
with Path(project_dir, "BHaHAHA.h").open("w", encoding="utf-8") as output_file:
    output_file.write(BHaHAHA_h)

#########################################################
# STEP 10: Final message
print(
    f"Finished! Now go into ./{outrootdir}/{project_name} and type `make` to build BHaHAHA.\n"
    "For help linking to your NR code, start by reading BHaHAHA.h"
)

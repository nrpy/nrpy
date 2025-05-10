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
import shutil
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile

# import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import (
    rfm_precompute,
    xx_tofrom_Cart,
)
from nrpy.infrastructures.BHaH.BHaHAHA import (
    apply_bcs_inner_only,
    apply_bcs_r_maxmin_partial_r_hDD_upwinding,
    bcstruct_set_up,
    cfl_limited_timestep_based_on_h_equals_r,
    diagnostics,
    diagnostics_area_centroid_and_Theta_norms,
    diagnostics_file_output,
    diagnostics_integration_weights,
    diagnostics_min_max_mean_radii_wrt_centroid,
    diagnostics_proper_circumferences,
    error_message,
    find_horizon,
    hDD_dD_and_W_dD_in_interp_src_grid_interior,
    initial_data,
    interpolation_1d_radial_spokes_on_3d_src_grid,
    interpolation_2d_external_input_to_interp_src_grid,
    interpolation_2d_general__uniform_src_grid,
    numgrid__evol_set_up,
    numgrid__external_input_set_up,
    numgrid__interp_src_set_up,
    over_relaxation,
    poisoning_check_inputs,
    poisoning_set_inputs,
    quadratic_extrapolation,
    radial_grid_cell_centered_set_up,
    rhs_eval_KO_apply,
    xyz_center_r_minmax,
)
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all

parser = argparse.ArgumentParser(
    description="BHaHAHA, the BlackHoles@Home Apparent Horizon Algorithm"
)
parser.add_argument(
    "--fdorder",
    type=int,
    help="Finite-difference order; interpolation order is one less than this. Default=6. For precision horizon finding, set to 8.",
    default=6,
)
args = parser.parse_args()
fd_order = args.fdorder

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
project_dir = os.path.join("project", project_name)

# Clean the project directory if it exists
shutil.rmtree(project_dir, ignore_errors=True)

#########################################################
# STEP 3: Define the damping parameter, min global wavespeed, and residual stop condition
par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

#########################################################
# STEP 4: Declare core C functions & register each to cfc.CFunction_dict["function_name"]
find_horizon.register_CFunction_find_horizon()
poisoning_set_inputs.register_CFunction_poisoning_set_inputs()
poisoning_check_inputs.register_CFunction_poisoning_check_inputs()
over_relaxation.register_CFunction_over_relaxation()
radial_grid_cell_centered_set_up.register_CFunction_radial_grid_cell_centered_set_up()
xyz_center_r_minmax.register_CFunction_xyz_center_r_minmax()
quadratic_extrapolation.register_CFunction_quadratic_extrapolation()
numgrid__external_input_set_up.register_CFunction_numgrid__external_input_set_up()
numgrid__interp_src_set_up.register_CFunction_numgrid__interp_src_set_up()
numgrid__evol_set_up.register_CFunction_numgrid__evol_set_up()

initial_data.register_CFunction_initial_data()
interpolation_1d_radial_spokes_on_3d_src_grid.register_CFunction_interpolation_1d_radial_spokes_on_3d_src_grid(
    enable_simd=enable_simd, project_dir=project_dir
)
interpolation_2d_external_input_to_interp_src_grid.register_CFunction_interpolation_2d_external_input_to_interp_src_grid()
interpolation_2d_general__uniform_src_grid.register_CFunction_interpolation_2d_general__uniform_src_grid(
    enable_simd=enable_simd, project_dir=project_dir
)
hDD_dD_and_W_dD_in_interp_src_grid_interior.register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior()
error_message.register_CFunction_error_message()

# Register numerical grid functions
cfl_limited_timestep_based_on_h_equals_r.register_CFunction_cfl_limited_timestep_based_on_h_equals_r()

xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem=CoordSystem)

diagnostics.register_CFunction_diagnostics()
diagnostics_area_centroid_and_Theta_norms.register_CFunction_diagnostics_area_centroid_and_Theta_norms()
diagnostics_file_output.register_CFunction_diagnostics_file_output()
diagnostics_integration_weights.register_CFunction_diagnostics_integration_weights()
diagnostics_min_max_mean_radii_wrt_centroid.register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid()
diagnostics_proper_circumferences.register_CFunction_diagnostics_proper_circumferences()

if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems={CoordSystem},
    )

rhs_eval_KO_apply.register_CFunction_rhs_eval(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
)
rhs_eval_KO_apply.register_CFunction_KO_apply(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
)

if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

#########################################################
# STEP 5: Register C function to set up the boundary condition struct
apply_bcs_inner_only.register_CFunction_apply_bcs_inner_only()
# cbc.register_CFunction_apply_bcs_inner_only()
apply_bcs_r_maxmin_partial_r_hDD_upwinding.register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding(
    upwinding_fd_order=fd_order
)
bcstruct_set_up.register_CFunction_bcstruct_set_up(CoordSystem="Spherical")

rhs_string = """
bah_rhs_eval(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
// Theta was just evaluated at Nxx1*Nxx2 gridpoints; update counter:
commondata->bhahaha_diagnostics->Theta_eval_points_counter += params->Nxx1 * params->Nxx2;
if(commondata->KO_diss_strength > 0.0)
  bah_KO_apply(commondata, params, rfmstruct,  auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
MoL_register_all.register_CFunctions(
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
CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
Bdefines_h.output_BHaH_defines_h(
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

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=f"lib{project_name}",
    lib_function_prefix="bah_",
    create_lib=True,
    static_lib=True,
)

# Append latest error codes & error message function prototype to BHaHAHA.h
# Read the contents of the original file
with Path("nrpy/infrastructures/BHaH/BHaHAHA/BHaHAHA_header.h").open(
    "r",
    encoding="utf-8",
) as input_file:
    BHaHAHA_h = input_file.read()
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
for item in error_message.error_code_msg_tuples_list:
    BHaHAHA_h += f"  {item[0]},\n"
BHaHAHA_h += "} bhahaha_error_codes;\n"
BHaHAHA_h += """
// Function: bah_error_message
// Interprets bah_find_horizon() error codes & returns a useful string.
const char *bah_error_message(const bhahaha_error_codes error_code);
//===============================================

#endif // BHAHAHA_HEADER_H
"""
# Write the updated content to the output file
with Path(project_dir, "BHaHAHA.h").open("w", encoding="utf-8") as output_file:
    output_file.write(BHaHAHA_h)

#########################################################
# STEP 10: Final message
print(
    f"Finished! Now go into ./project/{project_name} and type `make` to build BHaHAHA.\n"
    "For help linking to your NR code, start by reading BHaHAHA.h"
)

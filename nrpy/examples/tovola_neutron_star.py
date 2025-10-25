"""
Neutron star initial data example.

This example sets up a complete C code for setting up and
  validating initial data for a neutron star, using the
  TOVola initial data solver, written by David Boyer.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "tovola_neutron_star"
CoordSystem = "SinhSpherical"
IDtype = "TOVola_interp"
IDCoordSystem = "Spherical"

grid_physical_size = 10.0
sinh_width = 0.0785
t_final = 1.0e-10
diagnostics_output_every = 0.5
default_checkpoint_every = 2.0
Nxx_dict = {
    "SinhSpherical": [64, 16, 2],
}
enable_rfm_precompute = True
fd_order = 4
radiation_BC_fd_order = 4
enable_intrinsics = True
enable_parallel_codegen = True
enable_fd_functions = True
separately_compute_Ricci = False

OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    OMP_collapse = 2  # about 2x faster
    if CoordSystem == "SinhSpherical":
        sinh_width = 0.2
if "Cylindrical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "1")
    OMP_collapse = 2  # might be slightly faster

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

# Set NRPy parameters that steer the code generation
par.set_parval_from_str("enable_parallel_codegen", enable_parallel_codegen)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

BHaH.general_relativity.TOVola.TOVola_interp.register_CFunction_TOVola_interp()
BHaH.general_relativity.TOVola.TOVola_solve.register_CFunction_TOVola_solve()
BHaH.general_relativity.BSSN.initial_data.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    enable_checkpointing=True,
    ID_persist_struct_str=BHaH.general_relativity.TOVola.ID_persist_struct.ID_persist_str(),
    populate_ID_persist_struct_str=r"""
TOVola_solve(commondata, &ID_persist);
""",
    free_ID_persist_struct_str=r"""
{
  free(ID_persist.r_Schw_arr);
  free(ID_persist.rho_energy_arr);
  free(ID_persist.rho_baryon_arr);
  free(ID_persist.P_arr);
  free(ID_persist.M_arr);
  free(ID_persist.expnu_arr);
  free(ID_persist.exp4phi_arr);
  free(ID_persist.r_iso_arr);
}
""",
    enable_T4munu=True,
)
BHaH.general_relativity.BSSN.diagnostics.register_CFunction_diagnostics(
    set_of_CoordSystems={CoordSystem},
    default_diagnostics_out_every=diagnostics_output_every,
    enable_psi4_diagnostics=False,
    use_Ricci_eval_func=False,
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
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems={CoordSystem}
    )
BHaH.general_relativity.BSSN.constraints.register_CFunction_constraints_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=False,
    enable_T4munu=True,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

if __name__ == "__main__":
    pcg.do_parallel_codegen()
# Does not need to be parallelized.
BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems={CoordSystem}, radiation_BC_fd_order=radiation_BC_fd_order
)

rhs_string = ""
BHaH.MoLtimestepping.register_all.register_CFunctions(
    MoL_method="RK4",
    rhs_string=rhs_string,
    post_rhs_string="",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)

BHaH.xx_tofrom_Cart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
BHaH.checkpointing.register_CFunctions(
    default_checkpoint_every=default_checkpoint_every
)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

# Reset CodeParameter defaults according to variables set above.
# Coord system parameters
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", sinh_width)
par.adjust_CodeParam_default("t_final", t_final)

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
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
)
post_non_y_n_auxevol_mallocs = ""
BHaH.main_c.register_CFunction_main_c(
    MoL_method="",
    initial_data_desc=IDtype,
    set_initial_data_after_auxevol_malloc=True,
    boundary_conditions_desc="No BCs",
    post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
    pre_MoL_step_forward_in_time="",
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

if enable_intrinsics:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    compiler_opt_option="default",
    addl_CFLAGS=["$(shell gsl-config --cflags)"],
    addl_libraries=["$(shell gsl-config --libs)"],
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

# print(cfc.CFunction_dict["initial_data"].full_function)
# print(cfc.CFunction_dict["rhs_eval"].full_function)
# print(cfc.CFunction_dict["apply_bcs"].full_function)
# print(cfc.CFunction_dict["parameter_file_read_and_parse"].full_function)
# print(cfc.CFunction_dict["main"].full_function)

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
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.checkpointing as chkpt
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
import nrpy.infrastructures.BHaH.general_relativity.TOVola.ID_persist_struct as IDps
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_interp as TOVinterp
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_solve as TOVsolve
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numericalgrids
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCartxx
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import (
    griddata_commondata,
    rfm_precompute,
    rfm_wrapper_functions,
)
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "tovola_neutron_star"
CoordSystem = "SinhSpherical"
IDtype = "TOVola_interp"
IDCoordSystem = "Spherical"

grid_physical_size = 10.0
sinh_width = 0.0785
t_final = 0.0001
diagnostics_output_every = 0.5
default_checkpoint_every = 2.0
Nxx_dict = {
    "SinhSpherical": [64, 16, 2],
}
enable_rfm_precompute = True
fd_order = 4
radiation_BC_fd_order = 4
enable_simd = True
parallel_codegen_enable = True
enable_fd_functions = True
separately_compute_Ricci = False

OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    par.adjust_CodeParam_default("CFL_FACTOR", 1.0)
    OMP_collapse = 2  # about 2x faster
    if CoordSystem == "SinhSpherical":
        sinh_width = 0.2
if "Cylindrical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "1")
    par.adjust_CodeParam_default("CFL_FACTOR", 1.0)
    OMP_collapse = 2  # might be slightly faster

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

# Set NRPy parameters that steer the code generation
par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

TOVinterp.register_CFunction_TOVola_interp()
TOVsolve.register_CFunction_TOVola_solve()
BCl.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    enable_checkpointing=True,
    ID_persist_struct_str=IDps.ID_persist_str(),
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
BCl.register_CFunction_diagnostics(
    list_of_CoordSystems=[CoordSystem],
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
    rfm_precompute.register_CFunctions_rfm_precompute(
        list_of_CoordSystems=[CoordSystem]
    )
BCl.register_CFunction_constraints(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=False,
    enable_T4munu=True,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

if __name__ == "__main__":
    pcg.do_parallel_codegen()
# Does not need to be parallelized.
numericalgrids.register_CFunctions(
    list_of_CoordSystems=[CoordSystem],
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

cbc.CurviBoundaryConditions_register_C_functions(
    list_of_CoordSystems=[CoordSystem], radiation_BC_fd_order=radiation_BC_fd_order
)

rhs_string = ""
MoL_register_all.register_CFunctions(
    MoL_method="RK4",
    rhs_string=rhs_string,
    post_rhs_string="",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)

xxCartxx.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
xxCartxx.register_CFunction_xx_to_Cart(CoordSystem)
chkpt.register_CFunctions(default_checkpoint_every=default_checkpoint_every)
progress.register_CFunction_progress_indicator()
rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

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
CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_simd,
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
)
post_non_y_n_auxevol_mallocs = ""
main.register_CFunction_main_c(
    MoL_method="",
    initial_data_desc=IDtype,
    set_initial_data_after_auxevol_malloc=True,
    boundary_conditions_desc="No BCs",
    post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
    pre_MoL_step_forward_in_time="",
)
griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute, enable_CurviBCs=True
)

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

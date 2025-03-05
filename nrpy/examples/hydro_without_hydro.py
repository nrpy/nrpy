"""
Sets up a complete C code project for evolving a TOV solution in a dynamical spacetime with a static fluid ("hydro without hydro").

Author: Leonardo Rosa Werneck
        wernecklr **at** gmail **dot* com
"""

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import os
import shutil

import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
import nrpy.infrastructures.BHaH.general_relativity.TOVola.ID_persist_struct as IDps
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_interp as TOVinterp
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_solve as TOVsolve
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import (
    BHaH_defines_h,
    CodeParameters,
    Makefile_helpers,
    checkpointing,
    cmdline_input_and_parfiles,
    griddata_commondata,
    main_c,
    numerical_grids_and_timestep,
    rfm_precompute,
    rfm_wrapper_functions,
    xx_tofrom_Cart,
)
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "hydro_without_hydro"
IDtype = "TOVola_interp"
CoordSystem = "Spherical"
LapseEvolutionOption = "HarmonicSlicing"
ShiftEvolutionOption = "Frozen"
fp_type = "double"
grid_physical_size = 7.5
diagnostics_output_every = 0.25
default_checkpoint_every = 2.0
t_final = 1.0 * grid_physical_size
Nxx_dict = {
    "Spherical": [72, 12, 2],
    # "SinhSpherical": [72, 12, 2],
    # "Cartesian": [64, 64, 64],
}
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
enable_simd = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
enable_fd_functions = True
enable_KreissOliger_dissipation = False
enable_CAKO = True
enable_T4munu = True
boundary_conditions_desc = "outgoing radiation"

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

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("fp_type", fp_type)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
TOVinterp.register_CFunction_TOVola_interp()
TOVsolve.register_CFunction_TOVola_solve()
BCl.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem="Spherical",
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

numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

BCl.register_CFunction_diagnostics(
    set_of_CoordSystems={CoordSystem},
    default_diagnostics_out_every=diagnostics_output_every,
    enable_psi4_diagnostics=False,
    grid_center_filename_tuple=("out0d-conv_factor%.2f.txt", "convergence_factor"),
    axis_filename_tuple=(
        "out1d-AXIS-conv_factor%.2f-t%08.4f.txt",
        "convergence_factor, time",
    ),
    plane_filename_tuple=(
        "out2d-PLANE-conv_factor%.2f-t%08.4f.txt",
        "convergence_factor, time",
    ),
    out_quantities_dict="default",
)
if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems={CoordSystem},
    )
BCl.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=enable_T4munu,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    enable_CAKO=enable_CAKO,
    OMP_collapse=OMP_collapse,
)
if separate_Ricci_and_BSSN_RHS:
    BCl.register_CFunction_Ricci_eval(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_simd=enable_simd,
        enable_fd_functions=enable_fd_functions,
        OMP_collapse=OMP_collapse,
    )
BCl.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_constraints(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=enable_T4munu,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

if __name__ == "__main__":
    pcg.do_parallel_codegen()

cbc.CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems={CoordSystem},
    radiation_BC_fd_order=radiation_BC_fd_order,
)
rhs_string = ""
if separate_Ricci_and_BSSN_RHS:
    rhs_string += (
        "Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);"
    )
rhs_string += """
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")

MoL_register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(commondata, params, rfmstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
xx_tofrom_Cart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
checkpointing.register_CFunctions(default_checkpoint_every=default_checkpoint_every)
progress.register_CFunction_progress_indicator()
rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)

CodeParameters.write_CodeParameters_h_files(project_dir=project_dir)
CodeParameters.register_CFunctions_params_commondata_struct_set_to_default()
cmdline_input_and_parfiles.generate_default_parfile(
    project_dir=project_dir, project_name=project_name
)
cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_simd,
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
)

main_c.register_CFunction_main_c(
    initial_data_desc=IDtype,
    MoL_method=MoL_method,
    boundary_conditions_desc=boundary_conditions_desc,
    set_initial_data_after_auxevol_malloc=True,
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
Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
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

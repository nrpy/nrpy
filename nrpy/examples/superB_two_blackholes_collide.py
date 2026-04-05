"""
Two black holes collide example.

Specifically, evolve Brill-Lindquist initial data forward in time;
this example sets up a complete C code for solving the GR field
  equations in curvilinear coordinates on a cell-centered grid,
  using a reference metric approach.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu

Note: This is the superB version.

"""

import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import subprocess
from pathlib import Path

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH, superB

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "superB_two_blackholes_collide"
CoordSystem = "Spherical"
set_of_CoordSystems = {CoordSystem}
IDtype = "BrillLindquist"
IDCoordSystem = "Cartesian"
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
}
default_BH1_mass = default_BH2_mass = 0.5
default_BH1_z_posn = +0.5
default_BH2_z_posn = -0.5
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
enable_intrinsics = True
separate_Ricci_and_BSSN_RHS = True
enable_parallel_codegen = True
enable_fd_functions = True
enable_KreissOliger_dissipation = False
enable_CAKO = True
enable_BHaHAHA = False
outer_bcs_type = "radiation"
boundary_conditions_desc = "outgoing radiation"
# Number of chares, Nchare0, Nchare1, and Nchare2, in each direction,
# should be chosen such that Nxx0/Nchare0, Nxx1/Nchare1, Nxx2/Nchare2 are integers greater than NGHOSTS,
# NGHOSTS is fd_order/2
if "Spherical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 18)
    par.adjust_CodeParam_default("Nchare1", 2)
    par.adjust_CodeParam_default("Nchare2", 1)
if "Cylindrical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 18)
    par.adjust_CodeParam_default("Nchare1", 2)
    par.adjust_CodeParam_default("Nchare2", 1)
if "Cartesian" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 16)
    par.adjust_CodeParam_default("Nchare1", 16)
    par.adjust_CodeParam_default("Nchare2", 16)

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

par.adjust_CodeParam_default("t_final", t_final)

if enable_BHaHAHA:
    #########################################################
    # STEP 2: Declare core C functions & register each to
    #         cfc.CFunction_dict["function_name"]
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
                "--cpp",
                "--no-openmp",
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
                "--cpp",
                "--no-openmp",
            ],
            check=True,
        )
    from nrpy.infrastructures.BHaH.BHaHAHA import (
        interpolation_3d_general__uniform_src_grid,
    )
    from nrpy.infrastructures.superB import (
        BHaH_implementation,
    )

    BHaH_implementation.register_CFunction_bhahaha_find_horizons(
        CoordSystem=CoordSystem, max_horizons=3
    )
    interpolation_3d_general__uniform_src_grid.register_CFunction_interpolation_3d_general__uniform_src_grid(
        enable_simd=enable_intrinsics,
        project_dir=project_dir,
        use_cpp=True,
    )
    superB.interpolator3d_chare.output_interpolator3d_h_cpp_ci(project_dir=project_dir)
    superB.horizon_finder_chare.output_horizon_finder_h_cpp_ci(project_dir=project_dir)

#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
superB.initial_data.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    ID_persist_struct_str="",
)

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
superB.numerical_grids.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
superB.diagnostics.diagnostics.register_all_diagnostics(
    set_of_CoordSystems=set_of_CoordSystems,
    project_dir=project_dir,
    default_diagnostics_out_every=diagnostics_output_every,
    enable_nearest_diagnostics=True,
    enable_interp_diagnostics=False,
    enable_volume_integration_diagnostics=True,
    enable_free_auxevol=False,
)
BHaH.general_relativity.diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics=False, enable_psi4=False
)
superB.general_relativity.diagnostics_nearest.register_CFunction_diagnostics_nearest(
    CoordSystem
)

if enable_rfm_precompute:
    BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
        set_of_CoordSystems={CoordSystem}
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
BHaH.general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
    CoordSystem=CoordSystem,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
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

superB.chare_communication_maps.chare_comm_register_C_functions(
    set_of_CoordSystems={CoordSystem}
)
superB.CurviBoundaryConditions.CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems={CoordSystem},
    radiation_BC_fd_order=radiation_BC_fd_order,
    set_parity_on_aux=True,
)
rhs_string = ""
if separate_Ricci_and_BSSN_RHS:
    rhs_string += "Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);"
if outer_bcs_type == "radiation":
    rhs_string += """
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                   gridfunctions_wavespeed,gridfunctions_f_infinity,
                                   RK_INPUT_GFS, RK_OUTPUT_GFS);"""
else:
    rhs_string += """
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);"""

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
    post_rhs_string="""enforce_detgammabar_equals_detgammahat(params, rfmstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
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
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)
par.adjust_CodeParam_default("BH1_mass", default_BH1_mass)
par.adjust_CodeParam_default("BH2_mass", default_BH2_mass)
par.adjust_CodeParam_default("BH1_posn_z", default_BH1_z_posn)
par.adjust_CodeParam_default("BH2_posn_z", default_BH2_z_posn)
if enable_BHaHAHA:
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
            1.1 * (default_BH1_mass + default_BH2_mass),
        ],
    )
    par.adjust_CodeParam_default("bah_verbosity_level", 0)

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
    enable_BHaHAHA=enable_BHaHAHA,
)
superB.timestepping_chare.output_timestepping_h_cpp_ci_register_CFunctions(
    project_dir=project_dir,
    MoL_method=MoL_method,
    outer_bcs_type=outer_bcs_type,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_BHaHAHA=enable_BHaHAHA,
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
BHaH.BHaH_defines_h.output_BHaH_defines_h(
    additional_includes=[
        str(Path("superB") / "superB.h"),
        *([os.path.join(BHaHAHA_subdir, "BHaHAHA.h")] if enable_BHaHAHA else []),
    ],
    project_dir=project_dir,
    enable_rfm_precompute=enable_rfm_precompute,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
)

if enable_intrinsics:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

superB.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    addl_dirs_to_make=[*([BHaHAHA_subdir] if enable_BHaHAHA else [])],
    exec_or_library_name=project_name,
    compiler_opt_option="default",
    addl_CFLAGS=["-fpermissive "],
    addl_libraries=[
        "-module CkIO",
        *([f"-L{BHaHAHA_subdir}", f"-l{BHaHAHA_subdir}"] if enable_BHaHAHA else []),
    ],
    CC="charmc",
    enable_BHaHAHA=enable_BHaHAHA,
)

print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./charmrun +p4 ./{project_name} to run with 4 processors, for example."
)
print(f"    Parameter file can be found in {project_name}.par")

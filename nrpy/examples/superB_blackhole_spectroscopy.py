"""

Black hole spectroscopy example.

Specifically, evolve Brill-Lindquist initial data forward in time
  and monitors the ringing of the merged black hole over time via
  psi4.
This example sets up a complete C code for solving the GR field
  equations in curvilinear coordinates on a cell-centered grid,
  using a reference metric approach.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu

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
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.BHaHAHA.interpolation_2d_general__uniform_src_grid as interpolation2d
import nrpy.infrastructures.BHaH.checkpointing as chkpt
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numericalgrids
import nrpy.infrastructures.BHaH.special_functions.spin_weight_minus2_spherical_harmonics as swm2sh
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCartxx
import nrpy.infrastructures.superB.chare_communication_maps as charecomm
import nrpy.infrastructures.superB.CurviBoundaryConditions as superBcbc
import nrpy.infrastructures.superB.diagnostics as superBdiagnostics
import nrpy.infrastructures.superB.initial_data as superBinitialdata
import nrpy.infrastructures.superB.main_chare as superBmain
import nrpy.infrastructures.superB.Makefile_helpers as superBMakefile
import nrpy.infrastructures.superB.MoL as superBMoL
import nrpy.infrastructures.superB.numerical_grids as superBnumericalgrids
import nrpy.infrastructures.superB.superB.superB_pup as superBpup
import nrpy.infrastructures.superB.timestepping_chare as superBtimestepping
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import rfm_precompute, rfm_wrapper_functions
from nrpy.infrastructures.BHaH.general_relativity import (
    BSSN,
    PSI4,
    NRPyPN_quasicircular_momenta,
    TwoPunctures,
)

parser = argparse.ArgumentParser()
parser.add_argument(
    "--paper", action="store_true", help="use the paper-version parameters"
)
args = parser.parse_args()
paper = args.paper

par.set_parval_from_str("Infrastructure", "BHaH")


# Code-generation-time parameters:
project_name = "superB_blackhole_spectroscopy"
CoordSystem = "SinhCylindrical"
IDtype = "TP_Interp"
IDCoordSystem = "Cartesian"

initial_sep = 0.5 if not paper else 10.0
mass_ratio = 1.0  # must be >= 1.0. Will need higher resolution for > 1.0.
BH_m_chix = 0.0  # dimensionless spin parameter for less-massive BH
BH_M_chix = 0.0  # dimensionless spin parameter for more-massive BH
initial_p_r = 0.0  # want this to be <= 0.0. 0.0 -> fall from rest, < 0.0 -> boosted toward each other.
TP_npoints_A = 48
TP_npoints_B = 48
TP_npoints_phi = 4

enable_KreissOliger_dissipation = True
enable_CAKO = True
enable_CAHD = False
enable_SSL = True
KreissOliger_strength_gauge = 0.99
KreissOliger_strength_nongauge = 0.3
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 2.0
grid_physical_size = 300.0
diagnostics_output_every = 0.5
enable_charm_checkpointing = True
default_checkpoint_every = 20.0
t_final = 1.5 * grid_physical_size
swm2sh_maximum_l_mode_generated = 8
swm2sh_maximum_l_mode_to_compute = 2 if not paper else 8
if paper:
    list_of_psi4_extraction_radii = [80.0, 160.0]
    num_psi4_extraction_radii = len(list_of_psi4_extraction_radii)
Nxx_dict = {
    "SinhSpherical": [800, 16, 2],
    "SinhCylindrical": [400, 2, 1200] if not paper else [800, 2, 2400],
}
default_BH1_mass = default_BH2_mass = 0.5
default_BH1_z_posn = +0.25 if not paper else +5.0
default_BH2_z_posn = -0.25 if not paper else -5.0
enable_rfm_precompute = True
MoL_method = "RK4" if not paper else "SSPRK33"
fd_order = 8
radiation_BC_fd_order = 4 if not paper else 8
enable_intrinsics = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
enable_fd_functions = True
outer_bcs_type = "radiation"
boundary_conditions_desc = "outgoing radiation"
# Number of chares, Nchare0, Nchare1, and Nchare2, in each direction,
# should be chosen such that Nxx0/Nchare0, Nxx1/Nchare1, Nxx2/Nchare2 are integers greater than NGHOSTS,
# NGHOSTS is fd_order/2
if "Spherical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 10)
    par.adjust_CodeParam_default("Nchare1", 2)
    par.adjust_CodeParam_default("Nchare2", 1)
if "Cylindrical" in CoordSystem:
    par.adjust_CodeParam_default("Nchare0", 4)
    par.adjust_CodeParam_default("Nchare1", 1)
    par.adjust_CodeParam_default("Nchare2", 4)

OMP_collapse = 1
sinh_width = 0.2
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    par.adjust_CodeParam_default("CFL_FACTOR", 1.0)
    OMP_collapse = 2  # about 2x faster
if "Cylindrical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "1")
    par.adjust_CodeParam_default("CFL_FACTOR", 0.5)
    OMP_collapse = 2  # might be slightly faster

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

# Set NRPy parameters that steer the code generation
par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.set_parval_from_str(
    "swm2sh_maximum_l_mode_generated", swm2sh_maximum_l_mode_generated
)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
NRPyPN_quasicircular_momenta.register_CFunction_NRPyPN_quasicircular_momenta()
TwoPunctures.TwoPunctures_lib.register_C_functions()
superBinitialdata.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    enable_checkpointing=False,
    ID_persist_struct_str=TwoPunctures.ID_persist_struct.ID_persist_str(),
    populate_ID_persist_struct_str=r"""
initialize_ID_persist_struct(commondata, &ID_persist);
TP_solve(&ID_persist);
""",
    free_ID_persist_struct_str=r"""
{
  extern void free_derivs (derivs * v, int n);  // <- Needed to free memory allocated by TwoPunctures.
  // <- Free memory allocated within ID_persist.
  // Now that we're finished with par.v and par.cf_v (needed in setting up ID, we can free up memory for TwoPunctures' grids...
  free_derivs (&ID_persist.v,    ID_persist.npoints_A * ID_persist.npoints_B * ID_persist.npoints_phi);
  free_derivs (&ID_persist.cf_v, ID_persist.npoints_A * ID_persist.npoints_B * ID_persist.npoints_phi);
}
""",
)
interpolation2d.register_CFunction_interpolation_2d_general__uniform_src_grid(
    enable_simd=enable_intrinsics, project_dir=project_dir
)
superBdiagnostics.register_CFunction_diagnostics(
    set_of_CoordSystems={CoordSystem},
    default_diagnostics_out_every=diagnostics_output_every,
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
    enable_psi4_diagnostics=True,
    enable_L2norm_BSSN_constraints_diagnostics=True,
)
if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(set_of_CoordSystems={CoordSystem})
BSSN.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=False,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    KreissOliger_strength_gauge=KreissOliger_strength_gauge,
    KreissOliger_strength_nongauge=KreissOliger_strength_nongauge,
    enable_CAKO=enable_CAKO,
    enable_CAHD=enable_CAHD,
    enable_SSL=enable_SSL,
    OMP_collapse=OMP_collapse,
)
if enable_CAHD:
    BSSN.cahdprefactor_gf.register_CFunction_cahdprefactor_auxevol_gridfunction(
        {CoordSystem}
    )
if separate_Ricci_and_BSSN_RHS:
    BSSN.Ricci_eval.register_CFunction_Ricci_eval(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_intrinsics=enable_intrinsics,
        enable_fd_functions=enable_fd_functions,
        OMP_collapse=OMP_collapse,
    )
BSSN.enforce_detgammabar_equals_detgammahat.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BSSN.constraints.register_CFunction_constraints(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=False,
    enable_intrinsics=enable_intrinsics,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

PSI4.compute_psi4.register_CFunction_psi4(
    CoordSystem=CoordSystem,
    OMP_collapse=OMP_collapse,
    enable_fd_functions=enable_fd_functions,
)
swm2sh.register_CFunction_spin_weight_minus2_sph_harmonics()

if __name__ == "__main__":
    pcg.do_parallel_codegen()
# Does not need to be parallelized.
superBdiagnostics.register_CFunction_psi4_spinweightm2_decomposition_on_sphlike_grids()
superBdiagnostics.register_CFunction_psi4_spinweightm2_decomposition_on_cylindlike_grids()
superBdiagnostics.register_CFunction_psi4_spinweightm2_decomposition_file_write()

numericalgrids.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)
superBnumericalgrids.register_CFunctions(
    set_of_CoordSystems={CoordSystem},
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    enable_psi4_diagnostics=True,
)
charecomm.chare_comm_register_C_functions(set_of_CoordSystems={CoordSystem})
superBcbc.CurviBoundaryConditions_register_C_functions(
    set_of_CoordSystems={CoordSystem},
    radiation_BC_fd_order=radiation_BC_fd_order,
    set_parity_on_aux=True,
)

rhs_string = ""
if enable_SSL:
    rhs_string += """
// Set SSL strength (SSL_Gaussian_prefactor):
commondata->SSL_Gaussian_prefactor = commondata->SSL_h * exp(-commondata->time * commondata->time / (2 * commondata->SSL_sigma * commondata->SSL_sigma));
"""
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

superBMoL.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_bcs_str=post_rhs_bcs_str,
    post_rhs_string="""enforce_detgammabar_equals_detgammahat(params, rfmstruct, RK_OUTPUT_GFS);""",
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
if CoordSystem == "SinhCylindrical":
    par.adjust_CodeParam_default("AMPLRHO", grid_physical_size)
    par.adjust_CodeParam_default("AMPLZ", grid_physical_size)
    par.adjust_CodeParam_default("SINHWRHO", sinh_width)
    par.adjust_CodeParam_default("SINHWZ", sinh_width)

par.adjust_CodeParam_default("t_final", t_final)
# Initial data parameters
par.adjust_CodeParam_default("initial_sep", initial_sep)
par.adjust_CodeParam_default("mass_ratio", mass_ratio)
par.adjust_CodeParam_default("bbhxy_BH_m_chix", BH_m_chix)
par.adjust_CodeParam_default("bbhxy_BH_M_chix", BH_M_chix)
par.adjust_CodeParam_default("initial_p_t", 0.0)
par.adjust_CodeParam_default("initial_p_r", initial_p_r)
par.adjust_CodeParam_default("TP_npoints_A", TP_npoints_A)
par.adjust_CodeParam_default("TP_npoints_B", TP_npoints_B)
par.adjust_CodeParam_default("TP_npoints_phi", TP_npoints_phi)
par.adjust_CodeParam_default("TP_bare_mass_m", 1.0 / (1.0 + mass_ratio))
par.adjust_CodeParam_default("TP_bare_mass_M", mass_ratio / (1.0 + mass_ratio))
# Evolution / diagnostics parameters
par.adjust_CodeParam_default("eta", GammaDriving_eta)
par.adjust_CodeParam_default(
    "swm2sh_maximum_l_mode_to_compute", swm2sh_maximum_l_mode_to_compute
)
if paper:
    par.adjust_CodeParam_default("num_psi4_extraction_radii", num_psi4_extraction_radii)
    par.adjust_CodeParam_default(
        "list_of_psi4_extraction_radii",
        list_of_psi4_extraction_radii,
        new_cparam_type=f"REAL[{num_psi4_extraction_radii}]",
    )

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

copy_files(
    package="nrpy.infrastructures.BHaH.general_relativity.TwoPunctures",
    filenames_list=["TwoPunctures.h", "TP_utilities.h"],
    project_dir=project_dir,
    subdirectory="TwoPunctures",
)

superBmain.output_commondata_object_h_and_main_h_cpp_ci(
    project_dir=project_dir,
    enable_charm_checkpointing=enable_charm_checkpointing,
)

post_non_y_n_auxevol_mallocs = ""
if enable_CAHD:
    post_non_y_n_auxevol_mallocs = """for(int grid=0; grid<commondata.NUMGRIDS; grid++) {
    cahdprefactor_auxevol_gridfunction(&commondata, &griddata_chare[grid].params, griddata_chare[grid].xx,  griddata_chare[grid].gridfuncs.auxevol_gfs);
}\n"""

superBtimestepping.output_timestepping_h_cpp_ci_register_CFunctions(
    project_dir=project_dir,
    MoL_method=MoL_method,
    post_non_y_n_auxevol_mallocs=post_non_y_n_auxevol_mallocs,
    enable_rfm_precompute=enable_rfm_precompute,
    outer_bcs_type=outer_bcs_type,
    enable_psi4_diagnostics=True,
    enable_charm_checkpointing=enable_charm_checkpointing,
    enable_L2norm_BSSN_constraints_diagnostics=True,
)

superBpup.register_CFunction_superB_pup_routines(
    MoL_method=MoL_method,
    enable_psi4_diagnostics=True,
)
copy_files(
    package="nrpy.infrastructures.superB.superB",
    filenames_list=["superB.h", "superB_pup_function_prototypes.h"],
    project_dir=project_dir,
    subdirectory="superB",
)

Bdefines_h.output_BHaH_defines_h(
    additional_includes=[
        str(Path("TwoPunctures") / Path("TwoPunctures.h")),
        str(Path("superB") / Path("superB.h")),
    ],
    project_dir=project_dir,
    enable_intrinsics=enable_intrinsics,
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

superBMakefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    compiler_opt_option="default",
    addl_CFLAGS=["$(shell gsl-config --cflags)", "-fpermissive "],
    addl_libraries=[
        "$(shell gsl-config --libs)",
        "-module CkIO",
    ],
    CC="charmc",
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./charmrun +p4 ./{project_name} to run with 4 processors, for example."
)
print(f"To restart from checkpoint, run ./charmrun +p4 ./{project_name} +restart log")
print(f"    Parameter file can be found in {project_name}.par")

# print(cfc.CFunction_dict["initial_data"].full_function)
# print(cfc.CFunction_dict["rhs_eval"].full_function)
# print(cfc.CFunction_dict["apply_bcs"].full_function)
# print(cfc.CFunction_dict["parameter_file_read_and_parse"].full_function)
# print(cfc.CFunction_dict["main"].full_function)

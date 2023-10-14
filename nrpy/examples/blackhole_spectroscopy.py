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
"""
#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import os
from pathlib import Path

import nrpy.params as par
from nrpy.helpers import simd
import nrpy.helpers.parallel_codegen as pcg

from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
from nrpy.infrastructures.BHaH import rfm_precompute
from nrpy.infrastructures.BHaH import rfm_wrapper_functions
import nrpy.infrastructures.BHaH.general_relativity.TwoPunctures.TwoPunctures_lib as TPl
import nrpy.infrastructures.BHaH.general_relativity.TwoPunctures.ID_persist_struct as IDps
import nrpy.infrastructures.BHaH.general_relativity.NRPyPN_quasicircular_momenta as NRPyPNqm
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.checkpointing as chkpt
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numericalgrids
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCartxx
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
import nrpy.infrastructures.BHaH.special_functions.spin_weight_minus2_spherical_harmonics as swm2sh

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "blackhole_spectroscopy"
CoordSystem = "SinhSpherical"
IDtype = "TP_Interp"
IDCoordSystem = "Cartesian"

initial_sep = 0.5
mass_ratio = 1.0  # must be >= 1.0. Will need higher resolution for > 1.0.
BH_m_chix = 0.0  # dimensionless spin parameter for less-massive BH
BH_M_chix = 0.0  # dimensionless spin parameter for more-massive BH
initial_p_r = 0.0  # want this to be <= 0.0. 0.0 -> fall from rest, < 0.0 -> boosted toward each other.
TP_npoints_A = 48
TP_npoints_B = 48
TP_npoints_phi = 4

enable_KreissOliger_dissipation = True
KreissOliger_strength_mult_by_W = True
KreissOliger_strength_gauge = 0.99
KreissOliger_strength_nongauge = 0.1
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 2.0
grid_physical_size = 300.0
diagnostics_output_every = 0.5
default_checkpoint_every = 2.0
t_final = 1.5 * grid_physical_size
swm2sh_maximum_l_mode_generated = 8
swm2sh_maximum_l_mode_to_compute = 2  # for consistency with NRPy 1.0 version.
Nxx_dict = {
    "SinhSpherical": [800, 16, 2],
}
default_BH1_mass = default_BH2_mass = 0.5
default_BH1_z_posn = +0.25
default_BH2_z_posn = -0.25
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 8
radiation_BC_fd_order = 4
enable_simd = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
enable_fd_functions = True
boundary_conditions_desc = "outgoing radiation"

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
par.set_parval_from_str("enable_RbarDD_gridfunctions", separate_Ricci_and_BSSN_RHS)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.set_parval_from_str(
    "swm2sh_maximum_l_mode_generated", swm2sh_maximum_l_mode_generated
)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
NRPyPNqm.register_CFunction_NRPyPN_quasicircular_momenta()
TPl.register_C_functions()
BCl.register_CFunction_initial_data(
    CoordSystem=CoordSystem,
    IDtype=IDtype,
    IDCoordSystem=IDCoordSystem,
    enable_checkpointing=True,
    ID_persist_struct_str=IDps.ID_persist_str(),
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
BCl.register_CFunction_diagnostics(
    list_of_CoordSystems=[CoordSystem],
    default_diagnostics_out_every=diagnostics_output_every,
    enable_psi4_diagnostics=True,
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
BCl.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    KreissOliger_strength_mult_by_W=KreissOliger_strength_mult_by_W,
    KreissOliger_strength_gauge=KreissOliger_strength_gauge,
    KreissOliger_strength_nongauge=KreissOliger_strength_nongauge,
    OMP_collapse=OMP_collapse,
)
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
    enable_simd=enable_simd,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
swm2sh.register_CFunction_spin_weight_minus2_sph_harmonics()

for which_part in range(3):
    BCl.register_CFunction_psi4_part(
        CoordSystem=CoordSystem,
        which_part=which_part,
        enable_fd_functions=enable_fd_functions,
        OMP_collapse=OMP_collapse,
        output_empty_function=False,
    )
BCl.register_CFunction_psi4_tetrad(
    CoordSystem=CoordSystem,
    output_empty_function=False,
)
BCl.register_CFunction_psi4_spinweightm2_decomposition_on_sphlike_grids()

if __name__ == "__main__":
    pcg.do_parallel_codegen()

numericalgrids.register_CFunctions(
    list_of_CoordSystems=[CoordSystem],
    grid_physical_size=grid_physical_size,
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

cbc.CurviBoundaryConditions_register_C_functions(
    list_of_CoordSystems=[CoordSystem], radiation_BC_fd_order=radiation_BC_fd_order
)

rhs_string = """
Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")

MoL.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(commondata, params, rfmstruct, RK_OUTPUT_GFS);""",
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
TPl.copy_TwoPunctures_header_files(TwoPunctures_Path=Path(project_dir) / "TwoPunctures")
Bdefines_h.output_BHaH_defines_h(
    additional_includes=[str(Path("TwoPunctures") / Path("TwoPunctures.h"))],
    project_dir=project_dir,
    fin_NGHOSTS_add_one_for_upwinding=True,
    enable_simd=enable_simd,
)
main.register_CFunction_main_c(
    initial_data_desc=IDtype,
    pre_MoL_step_forward_in_time="write_checkpoint(&commondata, griddata);\n",
    MoL_method=MoL_method,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    boundary_conditions_desc=boundary_conditions_desc,
)

if enable_simd:
    simd.copy_simd_intrinsics_h(project_dir=project_dir)

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

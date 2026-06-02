"""
TOV example using GRoovy.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

import multiprocessing
import os
import shutil
import subprocess

import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures import BHaH
from nrpy.infrastructures.BHaH.GRoovy import diagnostics_nearest

#########################################################
# Step 1: Import needed Python modules, then set codegen
#         and compile-time parameters.

par.set_parval_from_str("Infrastructure", "BHaH")
par.set_parval_from_str("parallelization", "openmp")
par.set_parval_from_str("fp_type", "double")

# Code-generation-time parameters:
CoordSystem = "Spherical"  # Some options: "Spherical", "SinhSpherical", "Cylindrical", "Cartesian", "SinhCartesian"
project_name = "groovy_TOV_BSSN__" + CoordSystem
fp_type = "double"
IDType = "TOV"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
evolving_spacetime = True
evolving_temperature = False
grid_physical_size = 20.0
diagnostics_output_every = 0.25
default_checkpoint_every = 2.0
t_final = 2300.0
CFL_FACTOR = 0.2
# symmetry_axes will be set on any i such that Nxx[i] = 2 below.
Nxx_dict = {
    "Spherical": [100, 2, 2],
    "SinhSpherical": [64, 2, 2],
    "Cylindrical": [142, 2, 284],
    "Cartesian": [64, 64, 64],
    "SinhCartesian": [64, 64, 64],
}

enable_rfm_precompute = False
enable_intrinsics = False
MoL_method = "RK4"
# PPM requires at least 3 ghost zones
fd_order = 4
radiation_BC_fd_order = 4
separate_Ricci_and_BSSN_RHS = True
enable_parallel_codegen = True
enable_fd_functions = True

enable_KreissOliger_dissipation = True
eta = 0.2
KreissOliger_strength_gauge = 0.2
KreissOliger_strength_nongauge = 0.2

enable_CAKO = False
enable_CAHD = False
enable_SSL = False
boundary_conditions_desc = "outgoing radiation"
set_of_CoordSystems = {CoordSystem}

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

#########################################################
# Step 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]

BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=True,
    enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
    enable_T4munu=True,
    enable_intrinsics=True,
    enable_fd_functions=enable_fd_functions,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
    KreissOliger_strength_gauge=KreissOliger_strength_gauge,
    KreissOliger_strength_nongauge=KreissOliger_strength_nongauge,
    OMP_collapse=OMP_collapse,
)
if separate_Ricci_and_BSSN_RHS:
    BHaH.general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
        CoordSystem=CoordSystem,
        enable_intrinsics=True,
        enable_fd_functions=enable_fd_functions,
        OMP_collapse=OMP_collapse,
    )
BHaH.general_relativity.enforce_detgammabar_equals_detgammahat.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=True,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)
BHaH.general_relativity.constraints_eval.register_CFunction_constraints_eval(
    CoordSystem=CoordSystem,
    enable_T4munu=True,
    enable_fd_functions=enable_fd_functions,
    OMP_collapse=OMP_collapse,
)

# Register all GRHD grid functions needed for evolution
BHaH.GRoovy.register_all_grhd_gridfunctions.register_all_grhd_gridfunctions(
    CoordSystem=CoordSystem, evolving_spacetime=evolving_spacetime
)

rho_baryon_central = 1.28e-3
neos = 1
Gamma_poly_tab = 2
# rho_poly_tab = 0.0
K_poly_tab0 = 100

grhayl_setup_str = rf"""

//========================================================
// Beginning of Initialize GRHayL
//========================================================

  // This section sets up the initial parameters that would normally
  // be provided by the simulation.
  const int backup_routine[3] = {{Palenzuela1D,Font1D,None}};
  const bool calc_prims_guess = true; // GRHayL's primitive guesses are better to use than the previous time step, at least for hybrid eos
  const double Psi6threshold = 1e100;

  commondata->eos.eos_type = ghl_eos_hybrid;
  const int neos = {neos};
  const double W_max = 10.0;
  const double rho_b_min = {rho_baryon_central * 1e-13}; //*1e-9
  const double rho_b_max = 1e300;
  const double Gamma_th = 2.0;
  const double rho_ppoly[{max(neos - 1, 1)}] = {{0.0}};
  const double Gamma_ppoly[{neos}] = {{2.0}};
  const double k_ppoly0 = {K_poly_tab0};
  const double Lorentz_damping_factor = 0./0.0;

  // Here, we initialize the structs that are (usually) static during
  // a simulation.

  ghl_initialize_params(Noble2D,
                      backup_routine,
                      false, // evolve entropy
                      false, // evolve temperature
                      calc_prims_guess,
                      Psi6threshold,
                      W_max,
                      Lorentz_damping_factor,
                      &commondata->ghl_params);

  ghl_initialize_hybrid_eos_functions_and_params(
                                             rho_b_min, rho_b_min, rho_b_max,
                                             neos, rho_ppoly, Gamma_ppoly,
                                             k_ppoly0, Gamma_th, &commondata->eos);

//========================================================
// End of Initialize GRHayL
//========================================================
"""

BHaH.GRoovy.hybrid_EoS_TOV_initial_data.register_CFunction_hybrid_EoS_TOV_initial_data(
    grhayl_setup_str=grhayl_setup_str,
    CoordSystem=CoordSystem,
    OMP_collapse=OMP_collapse,
)

BHaH.general_relativity.TOVola.TOVola_interp.register_CFunction_TOVola_interp()
BHaH.general_relativity.TOVola.TOVola_solve.register_CFunction_TOVola_solve()

BHaH.general_relativity.ADM_Initial_Data_Reader__BSSN_Converter.register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
    CoordSystem=CoordSystem,
    enable_T4munu=True,
    enable_fd_functions=False,
    ID_persist_struct_str=BHaH.general_relativity.TOVola.ID_persist_struct.ID_persist_str(),
)

# Set parameters for TOV ID solve by TOVola
par.adjust_CodeParam_default("initial_central_density", rho_baryon_central)
par.adjust_CodeParam_default("poly_eos_Gamma", Gamma_poly_tab)
par.adjust_CodeParam_default("poly_eos_K", K_poly_tab0)

# grhayl calls currently incompatible with SIMD
BHaH.GRoovy.calculate_all_source_terms.register_CFunction_calculate_all_source_terms(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_intrinsics=False,
    OMP_collapse=OMP_collapse,
    enable_GoldenKernels=True,
    evolving_temperature=evolving_temperature,
)

for flux_dirn in range(3):
    BHaH.GRoovy.calculate_HLL_flux_dirn_i.register_CFunction_calculate_HLL_flux_dirn_i(
        flux_dirn=flux_dirn,
        enable_intrinsics=False,
        enable_GoldenKernels=True,
        evolving_temperature=evolving_temperature,
    )

BHaH.GRoovy.calculate_flux_divergences.register_CFunction_calculate_flux_divergences(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    OMP_collapse=OMP_collapse,
    enable_GoldenKernels=True,
    evolving_temperature=evolving_temperature,
)

BHaH.GRoovy.compute_up_index_velocity_time_component_pointwise.register_CFunction_compute_up_index_velocity_time_component_pointwise(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_GoldenKernels=True,
)

BHaH.GRoovy.reconstruction_loop.register_CFunction_reconstruction_loop_wenoz(
    CoordSystem=CoordSystem,
)

BHaH.GRoovy.interpolate_metric_gfs_to_cell_faces.register_CFunction_interpolate_metric_gfs_to_cell_faces(
    evolving_spacetime=evolving_spacetime,
)

BHaH.GRoovy.primitives_to_conservatives_routine.register_CFunction_primitives_to_conservatives_routine(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_intrinsics=False,
    OMP_collapse=OMP_collapse,
    enable_GoldenKernels=True,
)

BHaH.GRoovy.conservatives_to_primitives_routine.register_CFunction_conservatives_to_primitives_routine(
    CoordSystem=CoordSystem,
)

BHaH.GRoovy.basis_transform_Cartesian_to_rfm_basis.register_CFunction_basis_transform_Cartesian_to_rfm_basis(
    CoordSystem=CoordSystem,
)

BHaH.GRoovy.basis_transform_rfm_basis_to_Cartesian.register_CFunction_basis_transform_rfm_basis_to_Cartesian(
    CoordSystem=CoordSystem,
)

BHaH.GRoovy.basis_transform_rfm_basis_to_Cartesian__read_cons_only.register_CFunction_basis_transform_rfm_basis_to_Cartesian__read_cons_only(
    CoordSystem=CoordSystem,
)

BHaH.GRoovy.apply_copy_and_outflow_bcs.register_CFunction_apply_copy_and_outflow_bcs(
    CoordSystem=CoordSystem,
    enable_GoldenKernels=True,
    evolving_spacetime=evolving_spacetime,
)

BHaH.GRoovy.grhd_rhs_eval.register_CFunction_grhd_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
)

BHaH.GRoovy.compute_stress_energy_tensor.register_CFunction_compute_stress_energy_tensor(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_intrinsics=enable_intrinsics,
    OMP_collapse=OMP_collapse,
    enable_GoldenKernels=True,
    evolving_temperature=evolving_temperature,
)

# grvy_aux.register_CFunction_full_diagnostics(
#     set_of_CoordSystems=[CoordSystem],
#     enable_RbarDD_gridfunctions=True,
#     default_diagnostics_out_every=default_diagnostics_output_every,
#     grid_center_filename_tuple=("out0d-conv_factor%.2f.txt", "convergence_factor"),
#     axis_filename_tuple=(
#         "out1d-AXIS-conv_factor%.2f-t%08.2f.txt",
#         "convergence_factor, time",
#     ),
#     plane_filename_tuple=(
#         "out2d-PLANE-conv_factor%.2f-t%08.2f.txt",
#         "convergence_factor, time",
#     ),
#     out_quantities_dict="default",
# )

BHaH.numerical_grids_and_timestep.register_CFunctions(
    set_of_CoordSystems=set_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=True,
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
BHaH.general_relativity.diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics=False,
    enable_psi4=False,
    enable_T4munu=True,
)
diagnostics_nearest.register_CFunction_diagnostics_nearest(
    set_of_CoordSystems=set_of_CoordSystems,
    evolving_temperature=evolving_temperature,
)
BHaH.general_relativity.diagnostics_volume_integration.register_CFunction_diagnostics_volume_integration()
BHaH.rfm_precompute.register_CFunctions_rfm_precompute(
    set_of_CoordSystems=set_of_CoordSystems,
)


if __name__ == "__main__" and enable_parallel_codegen:
    pcg.do_parallel_codegen()

BHaH.CurviBoundaryConditions.register_all.register_C_functions(
    set_of_CoordSystems=set_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
    set_parity_on_auxevol=True,
)

rhs_string = r"""
REAL *restrict xx[3]; 
for(int ww=0;ww<3;ww++) 
    xx[ww] = griddata[grid].xx[ww];

Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);

grhd_rhs_eval(commondata, params, &commondata->ghl_params, &commondata->eos, griddata[grid].xx, RK_INPUT_GFS, auxevol_gfs, RK_OUTPUT_GFS);

apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata[grid].xx,
                                   gridfunctions_wavespeed,gridfunctions_f_infinity,
                                   RK_INPUT_GFS, RK_OUTPUT_GFS);
"""

post_RHS_string = r"""
enforce_detgammabar_equals_detgammahat(params, rfmstruct, RK_OUTPUT_GFS, auxevol_gfs);
conservatives_to_primitives_routine(commondata, params, &commondata->ghl_params, &commondata->eos, xx, RK_OUTPUT_GFS, auxevol_gfs);
apply_copy_and_outflow_bcs(commondata, params, &commondata->ghl_params, bcstruct, xx, RK_OUTPUT_GFS, auxevol_gfs);
compute_stress_energy_tensor(commondata, params, xx, &commondata->eos, RK_OUTPUT_GFS, auxevol_gfs);
"""

BHaH.MoLtimestepping.register_all.register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string=post_RHS_string,
    enable_rfm_precompute=True,
    enable_curviBCs=True,
)
BHaH.xx_tofrom_Cart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
BHaH.xx_tofrom_Cart.register_CFunction_xx_to_Cart(CoordSystem)
BHaH.checkpointing.register_CFunctions(
    default_checkpoint_every=default_checkpoint_every
)
BHaH.diagnostics.progress_indicator.register_CFunction_progress_indicator()
BHaH.rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# Step 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
par.adjust_CodeParam_default("t_final", t_final)
par.adjust_CodeParam_default("CFL_FACTOR", CFL_FACTOR)
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)

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
# gpu_defines_filename = BHaH.BHaH_device_defines_h.output_device_headers(
#     project_dir, num_streams=num_cuda_streams
# )

BHaH.griddata_commondata.register_griddata_commondata(
    "GRoovy",
    "ghl_eos_parameters eos",
    "GRHayL's equation of state (eos) struct, used globally",
    is_commondata=True,
)

BHaH.griddata_commondata.register_griddata_commondata(
    "GRoovy",
    "ghl_parameters ghl_params",
    "GRHayL's parameter struct, used globally",
    is_commondata=True,
)

additional_includes = [
    "ghl.h",
    "ghl_reconstruction.h",
    "ghl_atmosphere.h",
    "ghl_con2prim.h",
]

copy_files(
    package="nrpy.helpers",
    filenames_list=["simd_intrinsics.h"],
    project_dir=project_dir,
    subdirectory="intrinsics",
)

BHaH.BHaH_defines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    additional_includes=additional_includes,
    enable_rfm_precompute=True,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
    DOUBLE_means="double" if fp_type == "float" else "REAL",
)

post_initial = r"""
for (int grid = 0; grid < commondata.NUMGRIDS; grid++) {
  const params_struct *restrict params = &griddata[grid].params;
  const bc_struct *restrict bcstruct = &griddata[grid].bcstruct;
  const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;

  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

  apply_bcs_outerextrap_and_inner(&commondata, params, bcstruct, in_gfs);

  enforce_detgammabar_equals_detgammahat(params, rfmstruct, in_gfs, auxevol_gfs);

  primitives_to_conservatives_routine(&commondata, params, xx, &commondata.eos, auxevol_gfs, in_gfs);
  conservatives_to_primitives_routine(&commondata, params, &commondata.ghl_params, &commondata.eos, xx, in_gfs, auxevol_gfs);

  compute_stress_energy_tensor(&commondata, params, xx, &commondata.eos, in_gfs, auxevol_gfs);
} // END LOOP: for grid over all numerical grids
"""
BHaH.main_c.register_CFunction_main_c(
    initial_data_desc=IDType,
    MoL_method=MoL_method,
    set_initial_data_after_auxevol_malloc=True,
    # pre_MoL_step_forward_in_time="write_checkpoint(&commondata, griddata);\n",
    boundary_conditions_desc=boundary_conditions_desc,
    post_non_y_n_auxevol_mallocs=post_initial,
)
BHaH.griddata_commondata.register_CFunction_griddata_free(
    enable_rfm_precompute=True, enable_CurviBCs=True
)

# Set GRHayL directory
ghl_INC_FLAG = "GRHayL/include/ghl/"
ghl_LIB_FLAG = "-LGRHayL/lib"
BHaH.Makefile_helpers.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name=project_name,
    # compiler_opt_option="fast",
    compiler_opt_option="debug",
    include_dirs=[ghl_INC_FLAG],
    addl_libraries=[ghl_LIB_FLAG + " -lghl", "$(shell gsl-config --libs)"],
    addl_CFLAGS=["$(shell gsl-config --cflags)"],
)
if __name__ == "__main__":
    print("Cloning and compiling GRHayL...")

    # Define the repository URL and directory
    repo_url = "https://github.com/GRHayL/GRHayL.git"
    repo_dir = "GRHayL"
    codes_root_dir = f"project/{project_name}"  # Replace with the actual directory

    # Change to the codes root directory
    os.chdir(codes_root_dir)

    # Clone the repository
    result = subprocess.run(
        ["git", "clone", repo_url], capture_output=True, text=True, check=True
    )
    print(result.stdout)

    # Change to the repository directory
    os.chdir(repo_dir)

    # Configure the build
    configure_options = ["./configure", "--prefix=./", "--buildtype=opt"]
    # For Terrence's laptop, uncomment the following lines and comment out the above line
    # configure_options = ['./configure', '--hdf5inc', '/usr/include/hdf5/mpich',
    #                       '--hdf5lib', '/usr/lib/x86_64-linux-gnu/hdf5/mpich/',
    #                       '--prefix=./', '--buildtype=opt']
    result = subprocess.run(
        configure_options, capture_output=True, text=True, check=True
    )
    print(result.stdout)

    # Get the number of CPU cores
    cpus = str(multiprocessing.cpu_count())

    # Build and install
    result = subprocess.run(
        ["make", "-j" + cpus], capture_output=True, text=True, check=True
    )
    print(result.stdout)
    result = subprocess.run(
        ["make", "install"], capture_output=True, text=True, check=True
    )
    print(result.stdout)

    # Change to the project directory
    os.chdir("../")

    print(
        f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
    )
    print(
        "    GRHayL must already be available under the generated project or provided by the build environment."
    )
    print(f"    Parameter file can be found in {project_name}.par")

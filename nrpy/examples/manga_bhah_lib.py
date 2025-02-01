"""
Generate bhah_lib, a NRPy2 library used by the unstructured mesh hydrodynamics code MANGA.

Author: Leonardo Rosa Werneck
        wernecklr **at** gmail **dot* com
"""

import os

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil

import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.bhah_lib as lib
import nrpy.infrastructures.BHaH.checkpointing as chkpt
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
import nrpy.infrastructures.BHaH.general_relativity.TOVola.ID_persist_struct as IDps
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_interp as TOVinterp
import nrpy.infrastructures.BHaH.general_relativity.TOVola.TOVola_solve as TOVsolve
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numericalgrids
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCartxx
import nrpy.params as par
from nrpy.helpers.generic import copy_files
from nrpy.infrastructures.BHaH import rfm_precompute, rfm_wrapper_functions
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all

par.set_parval_from_str("Infrastructure", "BHaH")

# Basic project setup
project_name = "bhah_lib"
list_of_CoordSystems = [
    "Spherical",
    # "SinhCylindrical",
]
NUMGRIDS = len(list_of_CoordSystems)
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 1.0
grid_physical_size = 1.6
OMP_collapse = 1
enable_simd = True
Nxx_dict = {
    "Spherical": [96, 16, 2],
    # "SinhCylindrical": [144, 8, 24],
}

# BSSN evolution parameters
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 4
enable_simd = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
enable_fd_functions = True
enable_KreissOliger_dissipation = False
boundary_conditions_desc = "outgoing radiation"
diagnostics_output_every = 100
default_checkpoint_every = 2.0
enable_T4munu = True
IDtype = "TOVola_interp"
IDCoordSystem = "Spherical"
project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.adjust_CodeParam_default("NUMGRIDS", NUMGRIDS)
chkpt.register_CFunctions(default_checkpoint_every=default_checkpoint_every)
progress.register_CFunction_progress_indicator()

numericalgrids.register_CFunctions(
    list_of_CoordSystems=list_of_CoordSystems,
    list_of_grid_physical_sizes=[grid_physical_size],
    Nxx_dict=Nxx_dict,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
)

BCl.register_CFunction_diagnostics(
    list_of_CoordSystems=list_of_CoordSystems,
    default_diagnostics_out_every=diagnostics_output_every,
    out_quantities_dict={
        (
            "REAL",
            "log10HL",
        ): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))",
        ("REAL", "alphaL"): "y_n_gfs[IDX4pt(ALPHAGF, idx3)]",
        ("REAL", "hxx"): "y_n_gfs[IDX4pt(HDD00GF, idx3)]",
        ("REAL", "hxy"): "y_n_gfs[IDX4pt(HDD01GF, idx3)]",
        ("REAL", "hxz"): "y_n_gfs[IDX4pt(HDD02GF, idx3)]",
        ("REAL", "hyy"): "y_n_gfs[IDX4pt(HDD11GF, idx3)]",
        ("REAL", "hyz"): "y_n_gfs[IDX4pt(HDD12GF, idx3)]",
        ("REAL", "hzz"): "y_n_gfs[IDX4pt(HDD22GF, idx3)]",
        ("REAL", "T4UU00"): "auxevol_gfs[IDX4pt(T4UU00GF, idx3)]",
    },
    enable_progress_indicator=True,
)

# cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
#     project_name=project_name, cmdline_inputs=["CFL_factor"]
# )
TOVinterp.register_CFunction_TOVola_interp()
TOVsolve.register_CFunction_TOVola_solve()

for CoordSystem in list_of_CoordSystems:
    par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
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

    BCl.register_CFunction_rhs_eval(
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        enable_RbarDD_gridfunctions=separate_Ricci_and_BSSN_RHS,
        enable_simd=enable_simd,
        enable_T4munu=enable_T4munu,
        enable_fd_functions=enable_fd_functions,
        enable_KreissOliger_dissipation=enable_KreissOliger_dissipation,
        LapseEvolutionOption=LapseEvolutionOption,
        ShiftEvolutionOption=ShiftEvolutionOption,
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

if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(
        list_of_CoordSystems=list_of_CoordSystems
    )

if __name__ == "__main__":
    pcg.do_parallel_codegen()

cbc.CurviBoundaryConditions_register_C_functions(
    list_of_CoordSystems=list_of_CoordSystems,
    radiation_BC_fd_order=radiation_BC_fd_order,
)
rhs_string = """
Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata->xx,
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

for CoordSystem in list_of_CoordSystems:
    xxCartxx.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
    xxCartxx.register_CFunction_xx_to_Cart(CoordSystem)
rfm_wrapper_functions.register_CFunctions_CoordSystem_wrapper_funcs()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
if "SinhSpherical" in list_of_CoordSystems:
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)

CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_intrinsics=enable_simd,
    fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
    supplemental_defines_dict={
        "BHaH Lib": """
typedef struct BHaH_struct {
  commondata_struct *commondata;
  griddata_struct *griddata;
} BHaH_struct;"""
    },
)

if enable_simd:
    copy_files(
        package="nrpy.helpers",
        filenames_list=["simd_intrinsics.h"],
        project_dir=project_dir,
        subdirectory="intrinsics",
    )

lib.register_CFunctions_bhah_lib()

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir,
    project_name=project_name,
    exec_or_library_name="lib" + project_name,
    compiler_opt_option="default",
    create_lib=True,
    addl_CFLAGS=["$(shell gsl-config --cflags)"],
    addl_libraries=["$(shell gsl-config --libs)"],
)
print(
    f"Finished! Now go into project/{project_name} and type `make` to build, then ./{project_name} to run."
)
print(f"    Parameter file can be found in {project_name}.par")

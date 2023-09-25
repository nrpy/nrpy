"""
Sets up a complete C code project for solving the wave equation
  in curvilinear coordinates on a cell-centered grid,
  using a reference metric.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import os
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import cast, Union

import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.grid as gri
import nrpy.c_codegen as ccg
from nrpy.helpers import simd
import nrpy.helpers.parallel_codegen as pcg

from nrpy.equations.wave_equation.WaveEquationCurvilinear_RHSs import (
    WaveEquationCurvilinear_RHSs,
)
from nrpy.equations.wave_equation.InitialData import InitialData
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
from nrpy.infrastructures.BHaH import rfm_precompute
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numgrids
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "curviwavetoy"
WaveType = "SphericalGaussian"
default_sigma = 3.0
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size
CoordSystem = "Spherical"
Nxx_dict = {
    "Spherical": [64, 2, 2],
    "SinhSpherical": [64, 2, 2],
    "Cartesian": [64, 64, 64],
}
OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "12")
    OMP_collapse = 1  # about 2x faster
enable_rfm_precompute = True
MoL_method = "RK4"
fd_order = 4
radiation_BC_fd_order = 2
enable_simd = True
parallel_codegen_enable = True
boundary_conditions_desc = "outgoing radiation"

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_exact_solution_single_point(
    in_CoordSystem: str,
    in_WaveType: str = "SphericalGaussian",
    in_default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Registers the C function for the exact solution at a single point.

    :param in_CoordSystem: The coordinate system to use in setting up initial data.
    :param in_WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param in_default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Populate uu_ID, vv_ID
    ID = InitialData(
        CoordSystem=in_CoordSystem,
        WaveType=in_WaveType,
        default_sigma=in_default_sigma,
        default_k0=default_k0,
        default_k1=default_k1,
        default_k2=default_k2,
    )

    includes = ["BHaH_defines.h"]

    desc = r"""Exact solution at a single point."""
    c_type = "void"
    name = "exact_solution_single_point"
    params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xx0, const REAL xx1, const REAL xx2,  REAL *restrict exact_soln_UUGF, REAL *restrict exact_soln_VVGF
"""
    body = ccg.c_codegen(
        [ID.uu_ID, ID.vv_ID],
        ["*exact_soln_UUGF", "*exact_soln_VVGF"],
        verbose=False,
        include_braces=False,
    )
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initial_data() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial data function for the wave equation with specific parameters.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    c_type = "void"
    name = "initial_data"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
    vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
    # Unpack griddata
    body = r"""
for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
"""
    body += lp.simple_loop(
        loop_body="exact_solution_single_point(commondata, params, xx0,xx1,xx2,"
        f"&{uu_gf_memaccess},"
        f"&{vv_gf_memaccess});",
        read_xxs=True,
        loop_region="all points",
        OMP_collapse=OMP_collapse,
    )
    body += "}\n"
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


_ = par.CodeParameter(
    "REAL", __name__, "diagnostics_output_every", 0.2, commondata=True
)


def register_CFunction_diagnostics() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side evaluation function for the wave equation with specific parameters.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Diagnostics."""
    c_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
// Explanation of the if() below:
// Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
// Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
// Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
if(fabs(round(currtime / outevery) * outevery - currtime) < 0.5*currdt) {
for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
  // Unpack griddata struct:
  const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
  const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // 0D output
  {
    char filename[256];
    sprintf(filename, "out0d-conv_factor%.2f.txt", convergence_factor);
    FILE *outfile;
    if (nn == 0)
      outfile = fopen(filename, "w");
    else
      outfile = fopen(filename, "a");
    if (outfile == NULL) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }

    const int center_of_grid_idx = IDX3(NGHOSTS, Nxx_plus_2NGHOSTS1 / 2, Nxx_plus_2NGHOSTS2 / 2);

    const REAL num_soln_at_center_UUGF = y_n_gfs[IDX4pt(UUGF, center_of_grid_idx)];
    const REAL num_soln_at_center_VVGF = y_n_gfs[IDX4pt(VVGF, center_of_grid_idx)];
    REAL exact_soln_at_center_UUGF, exact_soln_at_center_VVGF;
    exact_solution_single_point(commondata, params, xx[0][NGHOSTS], xx[1][Nxx_plus_2NGHOSTS1 / 2], xx[2][Nxx_plus_2NGHOSTS2 / 2],
                                &exact_soln_at_center_UUGF, &exact_soln_at_center_VVGF);

    fprintf(outfile, "%e %e %e %e %e\n", time,
            fabs(fabs(num_soln_at_center_UUGF - exact_soln_at_center_UUGF) / exact_soln_at_center_UUGF),
            fabs(fabs(num_soln_at_center_VVGF - exact_soln_at_center_VVGF) / (1e-16 + exact_soln_at_center_VVGF)),
            num_soln_at_center_UUGF, exact_soln_at_center_UUGF);

    fclose(outfile);
  }


  // 1D output
  {
    char filename[256];
    sprintf(filename, "out1d-conv_factor%.2f-t%.2f.txt", convergence_factor, time);
    FILE *outfile;
    outfile = fopen(filename, "w");
    if (outfile == NULL) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }

    for(int i0=NGHOSTS;i0<Nxx_plus_2NGHOSTS0-NGHOSTS;i0++) {
      const int idx3 = IDX3(i0, Nxx_plus_2NGHOSTS1 / 2, Nxx_plus_2NGHOSTS2 / 2);

      const REAL num_soln_UUGF = y_n_gfs[IDX4pt(UUGF, idx3)];
      const REAL num_soln_VVGF = y_n_gfs[IDX4pt(VVGF, idx3)];
      REAL exact_soln_UUGF, exact_soln_VVGF;
      exact_solution_single_point(commondata, params, xx[0][i0], xx[1][Nxx_plus_2NGHOSTS1 / 2], xx[2][Nxx_plus_2NGHOSTS2 / 2],
                                  &exact_soln_UUGF, &exact_soln_VVGF);

      fprintf(outfile, "%e %e %e %e %e\n", xx[0][i0],
              fabs(fabs(num_soln_UUGF - exact_soln_UUGF) / exact_soln_UUGF),
              fabs(fabs(num_soln_VVGF - exact_soln_VVGF) / (1e-16 + exact_soln_VVGF)),
              num_soln_UUGF, exact_soln_UUGF);
    }
    fclose(outfile);
  }
}
}
progress_indicator(commondata, griddata);
if(commondata->time + commondata->dt > commondata->t_final) printf("\n");
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_rhs_eval(enable_rfm_pre: bool) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the wave equation.

    This function sets the right-hand side of the wave equation according to the
    selected coordinate system and specified parameters.

    :param enable_rfm_pre: Whether or not to enable reference metric precomputation.
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
    desc = r"""Set RHSs for wave equation."""
    c_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_pre:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate uu_rhs, vv_rhs
    rhs = WaveEquationCurvilinear_RHSs(CoordSystem, enable_rfm_pre)
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.uu_rhs, rhs.vv_rhs],
            [
                gri.BHaHGridFunction.access_gf("uu", gf_array_name="rhs_gfs"),
                gri.BHaHGridFunction.access_gf("vv", gf_array_name="rhs_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_pre,
        read_xxs=not enable_rfm_pre,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


register_CFunction_exact_solution_single_point(
    in_CoordSystem=CoordSystem, in_WaveType=WaveType, in_default_sigma=default_sigma
)
register_CFunction_initial_data()
numgrids.register_CFunction_numerical_grids_and_timestep_setup(
    CoordSystem, grid_physical_size, Nxx_dict
)
register_CFunction_diagnostics()
if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(CoordSystem)
register_CFunction_rhs_eval(enable_rfm_precompute)


if __name__ == "__main__" and parallel_codegen_enable:
    pcg.do_parallel_codegen()

cbc.CurviBoundaryConditions_register_C_functions(
    CoordSystem, radiation_BC_fd_order=radiation_BC_fd_order
)
rhs_string = """rhs_eval(commondata, params, rfmstruct,  RK_INPUT_GFS, RK_OUTPUT_GFS);
if (strncmp(commondata->outer_bc_type, "radiation", 50) == 0)
  apply_bcs_outerradiation_and_inner(commondata, params, bcstruct, griddata->xx,
                                     gridfunctions_wavespeed,gridfunctions_f_infinity,
                                     RK_INPUT_GFS, RK_OUTPUT_GFS);"""
if not enable_rfm_precompute:
    rhs_string = rhs_string.replace("rfmstruct", "xx")
MoL.MoL_register_CFunctions(
    MoL_method=MoL_method,
    rhs_string=rhs_string,
    post_rhs_string="""if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
progress.register_CFunction_progress_indicator()

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
    enable_simd=enable_simd,
    MoL_method=MoL_method,
    CoordSystem=CoordSystem,
)
main.register_CFunction_main_c(
    initial_data_desc=WaveType,
    MoL_method=MoL_method,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_CurviBCs=True,
    boundary_conditions_desc=boundary_conditions_desc,
)

if enable_simd:
    simd.copy_simd_intrinsics_h(project_dir=project_dir)

Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir=project_dir, project_name=project_name, exec_name=project_name
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

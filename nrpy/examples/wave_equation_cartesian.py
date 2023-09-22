"""
Sets up a complete C code project for solving the wave equation
  in Cartesian coordinates, on a cell-centered Cartesian
  grid.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
#########################################################
# STEP 1: Import needed Python modules, then set codegen
#         and compile-time parameters.
import shutil
import os

import nrpy.params as par
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.c_codegen as ccg
from nrpy.helpers import simd

from nrpy.equations.wave_equation.WaveEquation_RHSs import WaveEquation_RHSs
from nrpy.equations.wave_equation.InitialData import InitialData
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress

par.set_parval_from_str("Infrastructure", "BHaH")

# Code-generation-time parameters:
project_name = "wavetoy"
WaveType = "SphericalGaussian"
default_sigma = 2.0
MoL_method = "RK4"
fd_order = 4
enable_simd = True
grid_physical_size = 10.0
t_final = 0.8 * grid_physical_size

project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("fd_order", fd_order)
par.adjust_CodeParam_default("t_final", t_final)


# fmt: off
for i in range(3):
    _ = par.CodeParameter("int", __name__, f"Nxx_plus_2NGHOSTS{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("int", __name__, f"Nxx{i}", 64)
    _ = par.CodeParameter("REAL", __name__, f"xxmin{i}", -grid_physical_size)
    _ = par.CodeParameter("REAL", __name__, f"xxmax{i}", grid_physical_size)
    _ = par.CodeParameter("REAL", __name__, f"invdxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
    _ = par.CodeParameter("REAL", __name__, f"dxx{i}", add_to_parfile=False, add_to_set_CodeParameters_h=True)
_ = par.CodeParameter("REAL", __name__, "convergence_factor", 1.0, commondata=True)
# fmt: on


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_numerical_grids_and_timestep_setup() -> None:
    """
    Register a C function to set up numerical grids and time steps for simulations.

    This function calls the cfc.register_CFunction method to register a C function
    responsible for the setup of numerical grids and time steps. It defines essential parameters,
    including grid dimensions, time steps, and other related quantities.

    :return: None
    """
    includes = ["BHaH_defines.h"]
    desc = r"""Set up cell-centered Cartesian grids."""
    c_type = "void"
    name = "numerical_grids_and_timestep_setup"
    params = "commondata_struct *restrict commondata, params_struct *restrict params, griddata_struct *restrict griddata"
    body = r"""
params->Nxx0 *= convergence_factor;
params->Nxx1 *= convergence_factor;
params->Nxx2 *= convergence_factor;

params->Nxx_plus_2NGHOSTS0 = params->Nxx0 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS1 = params->Nxx1 + 2*NGHOSTS;
params->Nxx_plus_2NGHOSTS2 = params->Nxx2 + 2*NGHOSTS;

params->dxx0 = (xxmax0 - xxmin0) / ((REAL)params->Nxx0);
params->dxx1 = (xxmax1 - xxmin1) / ((REAL)params->Nxx1);
params->dxx2 = (xxmax2 - xxmin2) / ((REAL)params->Nxx2);

params->invdxx0 = ((REAL)params->Nxx0) / (xxmax0 - xxmin0);
params->invdxx1 = ((REAL)params->Nxx1) / (xxmax1 - xxmin1);
params->invdxx2 = ((REAL)params->Nxx2) / (xxmax2 - xxmin2);

// Time parameters:
commondata->nn = 0;
commondata->nn_0 = 0;
commondata->t_0 = 0.0;
commondata->dt = CFL_FACTOR * MIN(params->dxx0,MIN(params->dxx1,params->dxx2)); // CFL condition
commondata->time = 0.0;

// Set up cell-centered Cartesian coordinate grid, centered at the origin.
griddata->xx[0] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS0);
griddata->xx[1] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS1);
griddata->xx[2] = (REAL *restrict)malloc(sizeof(REAL)*params->Nxx_plus_2NGHOSTS2);
for(int j=0;j<params->Nxx_plus_2NGHOSTS0;j++) griddata->xx[0][j] = xxmin0 + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx0;
for(int j=0;j<params->Nxx_plus_2NGHOSTS1;j++) griddata->xx[1][j] = xxmin1 + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx1;
for(int j=0;j<params->Nxx_plus_2NGHOSTS2;j++) griddata->xx[2][j] = xxmin2 + ((REAL)(j-NGHOSTS) + (1.0/2.0))*params->dxx2;
"""
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_exact_solution_single_point(
    in_WaveType: str = "SphericalGaussian",
    in_default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> None:
    """
    Registers the C function for the exact solution at a single point.

    :param in_WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param in_default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction
    """

    # Populate uu_ID, vv_ID
    ID = InitialData(
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


def register_CFunction_initial_data() -> None:
    """
    Register the initial data function for the wave equation with specific parameters.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    c_type = "void"
    name = "initial_data"
    params = "const commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    uu_gf_obj = gri.glb_gridfcs_dict["uu"]
    vv_gf_obj = gri.glb_gridfcs_dict["vv"]
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
        f"&{uu_gf_obj.read_gf_from_memory_Ccode_onept()},"
        f"&{vv_gf_obj.read_gf_from_memory_Ccode_onept()});",
        read_xxs=True,
        loop_region="all points",
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


_ = par.CodeParameter(
    "REAL", __name__, "diagnostics_output_every", 0.2, commondata=True
)


def register_CFunction_diagnostics() -> None:
    """
    Register the right-hand side evaluation function for the wave equation with specific parameters.
    """
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

    const int i0_center = Nxx_plus_2NGHOSTS0 / 2;
    const int i1_center = Nxx_plus_2NGHOSTS1 / 2;
    const int i2_center = Nxx_plus_2NGHOSTS2 / 2;
    const int center_of_grid_idx = IDX3(i0_center, i1_center, i2_center);

    const REAL num_soln_at_center_UUGF = y_n_gfs[IDX4pt(UUGF, center_of_grid_idx)];
    const REAL num_soln_at_center_VVGF = y_n_gfs[IDX4pt(VVGF, center_of_grid_idx)];
    REAL exact_soln_at_center_UUGF, exact_soln_at_center_VVGF;
    exact_solution_single_point(commondata, params, xx[0][i0_center], xx[1][i1_center], xx[2][i2_center],
                                &exact_soln_at_center_UUGF, &exact_soln_at_center_VVGF);

    fprintf(outfile, "%e %e %e %e %e\n", time,
            fabs(fabs(num_soln_at_center_UUGF - exact_soln_at_center_UUGF) / exact_soln_at_center_UUGF),
            fabs(fabs(num_soln_at_center_VVGF - exact_soln_at_center_VVGF) / (1e-16 + exact_soln_at_center_VVGF)),
            num_soln_at_center_UUGF, exact_soln_at_center_UUGF);

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


def register_CFunction_rhs_eval() -> None:
    """
    Register the right-hand side evaluation function for the wave equation with specific parameters.
    """
    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [os.path.join("simd", "simd_intrinsics.h")]
    desc = r"""Set RHSs for wave equation."""
    c_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    # Populate uu_rhs, vv_rhs
    rhs = WaveEquation_RHSs()
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


def register_CFunction_apply_bcs() -> None:
    """
    Register a C function to apply boundary conditions to gridfunctions.

    This function registers a C function responsible for applying spatial boundary
    conditions to the scalar wave grid functions. Quadratic polynomial extrapolation
    is used to update boundary conditions on all six faces of the 3D grid cube.

    :return: None
    """
    includes = ["BHaH_defines.h"]
    desc = """Apply (quadratic extrapolation) spatial boundary conditions to the scalar wave gridfunctions.
BCs are applied to all six boundary faces of the cube, filling in the innermost
ghost zone first, and moving outward."""
    c_type = "void"
    name = "apply_bcs"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params,REAL *restrict gfs"
    prefunc = r"""
// Declare boundary condition FACE_UPDATE macro,
//          which updates a single face of the 3D grid cube
//          using quadratic polynomial extrapolation.
const int MAXFACE = -1;
const int NUL     = +0;
const int MINFACE = +1;
#define FACE_UPDATE(which_gf, i0min,i0max, i1min,i1max, i2min,i2max, FACEX0,FACEX1,FACEX2) \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++) { \
        gfs[IDX4(which_gf,i0,i1,i2)] =                                  \
          +3.0*gfs[IDX4(which_gf,i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2)]  \
          -3.0*gfs[IDX4(which_gf,i0+2*FACEX0,i1+2*FACEX1,i2+2*FACEX2)]  \
          +1.0*gfs[IDX4(which_gf,i0+3*FACEX0,i1+3*FACEX1,i2+3*FACEX2)]; \
      }
"""
    body = r"""
#pragma omp parallel for
    for(int which_gf=0;which_gf<NUM_EVOL_GFS;which_gf++) {
#include "set_CodeParameters.h"
      int imin[3] = { NGHOSTS, NGHOSTS, NGHOSTS };
      int imax[3] = { Nxx_plus_2NGHOSTS0-NGHOSTS, Nxx_plus_2NGHOSTS1-NGHOSTS, Nxx_plus_2NGHOSTS2-NGHOSTS };
      for(int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
        // After updating each face, adjust imin[] and imax[]
        //   to reflect the newly-updated face extents.
        FACE_UPDATE(which_gf, imin[0]-1,imin[0], imin[1],imax[1], imin[2],imax[2], MINFACE,NUL,NUL); imin[0]--;
        FACE_UPDATE(which_gf, imax[0],imax[0]+1, imin[1],imax[1], imin[2],imax[2], MAXFACE,NUL,NUL); imax[0]++;

        FACE_UPDATE(which_gf, imin[0],imax[0], imin[1]-1,imin[1], imin[2],imax[2], NUL,MINFACE,NUL); imin[1]--;
        FACE_UPDATE(which_gf, imin[0],imax[0], imax[1],imax[1]+1, imin[2],imax[2], NUL,MAXFACE,NUL); imax[1]++;

        FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imin[2]-1,imin[2], NUL,NUL,MINFACE); imin[2]--;
        FACE_UPDATE(which_gf, imin[0],imax[0], imin[1],imax[1], imax[2],imax[2]+1, NUL,NUL,MAXFACE); imax[2]++;
    }
  }
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


register_CFunction_exact_solution_single_point()
register_CFunction_initial_data()
register_CFunction_numerical_grids_and_timestep_setup()
register_CFunction_diagnostics()
register_CFunction_rhs_eval()
register_CFunction_apply_bcs()
MoL.MoL_register_CFunctions(
    MoL_method=MoL_method,
    rhs_string="rhs_eval(commondata, params,  RK_INPUT_GFS, RK_OUTPUT_GFS);",
    post_rhs_string="apply_bcs(commondata, params,  RK_OUTPUT_GFS);",
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
    project_dir=project_dir, enable_simd=enable_simd, MoL_method=MoL_method
)
main.register_CFunction_main_c(
    MoL_method=MoL_method,
    initial_data_desc=WaveType,
    boundary_conditions_desc="Quadratic extrapolation, manually defined",
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

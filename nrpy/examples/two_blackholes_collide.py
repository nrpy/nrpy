"""
Sets up a complete C code for solving the general relativistic
  field equations in curvilinear coordinates on a cell-centered
  grid, using a reference metric.

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
import time

import nrpy.params as par
import nrpy.c_function as cfc
from nrpy.helpers import simd
import nrpy.helpers.parallel_codegen as pcg

from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
from nrpy.infrastructures.BHaH import rfm_precompute
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.infrastructures.BHaH.CodeParameters as CPs
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.infrastructures.BHaH.main_c as main
import nrpy.infrastructures.BHaH.diagnostics.progress_indicator as progress
import nrpy.infrastructures.BHaH.Makefile_helpers as Makefile
import nrpy.infrastructures.BHaH.cmdline_input_and_parfiles as cmdpar
import nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions as cbc
import nrpy.infrastructures.BHaH.numerical_grids_and_timestep as numgrids
import nrpy.infrastructures.BHaH.xx_tofrom_Cart as xxCart
import nrpy.infrastructures.BHaH.general_relativity.BSSN_C_codegen_library as BCl
import nrpy.infrastructures.BHaH.special_functions.spin_weight_minus2_spherical_harmonics as swm2sh

par.set_parval_from_str("Infrastructure", "BHaH")
start_time = time.time()

# Code-generation-time parameters:
project_name = "two_blackholes_collide"
CoordSystem = "Spherical"
IDtype = "BrillLindquist"
IDCoordSystem = "Cartesian"
LapseEvolutionOption = "OnePlusLog"
ShiftEvolutionOption = "GammaDriving2ndOrder_Covariant"
GammaDriving_eta = 1.0
grid_physical_size = 7.5
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
fd_order = 8
radiation_BC_fd_order = 8
enable_simd = True
separate_Ricci_and_BSSN_RHS = True
parallel_codegen_enable = True
boundary_conditions_desc = "outgoing radiation"

OMP_collapse = 1
if "Spherical" in CoordSystem:
    par.set_parval_from_str("symmetry_axes", "2")
    OMP_collapse = 2  # about 2x faster
project_dir = os.path.join("project", project_name)

# First clean the project directory, if it exists.
shutil.rmtree(project_dir, ignore_errors=True)

par.set_parval_from_str("parallel_codegen_enable", parallel_codegen_enable)
par.set_parval_from_str("fd_order", fd_order)
par.set_parval_from_str("LeaveRicciSymbolic", separate_Ricci_and_BSSN_RHS)
par.set_parval_from_str("CoordSystem_to_register_CodeParameters", CoordSystem)

# In the following function we overwrite t_final (registered as a CodeParameter in MoL),
#    so let's remove it from the parfile
par.adjust_CodeParam_default("t_final", t_final)


#########################################################
# STEP 2: Declare core C functions & register each to
#         cfc.CFunction_dict["function_name"]
def register_CFunction_diagnostics(
    in_CoordSystem: str, plane: str = "yz"
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for simulation diagnostics.

    :param in_CoordSystem: Specifies the coordinate system for the diagnostics.
    :param plane: The default plane for diagnostics; defaults to "yz".
    :return: None
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    _ = par.CodeParameter(
        "REAL", __name__, "diagnostics_output_every", 0.25, commondata=True
    )

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
    REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
    const params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

    // Constraint output
    {
      Ricci_eval(commondata, params, &griddata->rfmstruct, y_n_gfs, auxevol_gfs);
      constraints_eval(commondata, params, &griddata->rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }
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

      const int r_mid_idx = IDX3(Nxx_plus_2NGHOSTS1 / 2, Nxx_plus_2NGHOSTS1 / 2, Nxx_plus_2NGHOSTS2 / 2);

      const REAL H_at_center = diagnostic_output_gfs[IDX4pt(HGF, r_mid_idx)];
      const REAL M2_at_center = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, r_mid_idx)];
      const REAL alpha_at_center = y_n_gfs[IDX4pt(ALPHAGF, r_mid_idx)];
      const REAL trK_at_center = y_n_gfs[IDX4pt(TRKGF, r_mid_idx)];

      fprintf(outfile, "%e %e %e %e %e\n", time, log10(fabs(H_at_center + 1e-16)), log10(fabs(M2_at_center + 1e-16)), alpha_at_center, trK_at_center);

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
      for (int i1 = Nxx_plus_2NGHOSTS1 - NGHOSTS - 1; i1 >= NGHOSTS; i1 -= Nxx1-1) {
        int i0_start, i0_end, i0_step;

        if (i1 == (Nxx_plus_2NGHOSTS1 - NGHOSTS - 1)) {
          i0_start = Nxx_plus_2NGHOSTS0 - NGHOSTS - 1;
          i0_end = NGHOSTS - 1;
          i0_step = -1;
        } else if (i1 == NGHOSTS) {
          i0_start = NGHOSTS;
          i0_end = Nxx_plus_2NGHOSTS0 - NGHOSTS;
          i0_step = 1;
        } else {
          continue; // Skip this iteration if i1 is not one of the special values
        }

        for (int i0 = i0_start; i0 != i0_end; i0 += i0_step) {
          const int i2 = Nxx_plus_2NGHOSTS2 / 2;
          const int idx3 = IDX3(i0, i1, i2);
          REAL xCart[3];
          xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

          const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
          const REAL trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
          const REAL HL = diagnostic_output_gfs[IDX4pt(HGF, idx3)];
          const REAL M2L = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)];

          fprintf(outfile, "%e %e %e %e %e\n", xCart[2], log10(fabs(HL + 1e-16)), log10(fabs(M2L + 1e-16)), alphaL, trKL);
        }
      }
      fclose(outfile);
    }
"""
    body += rf"""
    // 2D output:
    {{
      char filename[256];
      sprintf(filename, "out2d-{plane}_plane-conv_factor%.2f-t%.2f.txt", convergence_factor, time);
      FILE *outfile;
      outfile = fopen(filename, "w");
      if (outfile == NULL) {{
        fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
        exit(1);
      }}
"""
    body += lp.simple_loop_2D(
        loop_body=r"""
const int idx3 = IDX3(i0, i1, i2);
REAL xCart[3];
xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

const REAL alphaL = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
const REAL trKL = y_n_gfs[IDX4pt(TRKGF, idx3)];
const REAL HL = diagnostic_output_gfs[IDX4pt(HGF, idx3)];
const REAL M2L = diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)];

fprintf(outfile, "%e %e %e %e %e %e %e\n", xCart[0], xCart[1], xCart[2], log10(fabs(HL + 1e-16)), log10(fabs(M2L + 1e-16)), alphaL, trKL);
""",
        CoordSystem=in_CoordSystem,
        plane=plane,
    )
    body += r"""
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


BCl.register_CFunction_initial_data(
    CoordSystem=CoordSystem, IDtype=IDtype, IDCoordSystem=IDCoordSystem
)
numgrids.register_CFunction_numerical_grids_and_timestep_setup(
    CoordSystem, grid_physical_size, Nxx_dict
)
register_CFunction_diagnostics(in_CoordSystem=CoordSystem)
if enable_rfm_precompute:
    rfm_precompute.register_CFunctions_rfm_precompute(CoordSystem)
BCl.register_CFunction_rhs_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    LapseEvolutionOption=LapseEvolutionOption,
    ShiftEvolutionOption=ShiftEvolutionOption,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_Ricci_eval(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_enforce_detgammabar_equals_detgammahat(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    OMP_collapse=OMP_collapse,
)
BCl.register_CFunction_constraints(
    CoordSystem=CoordSystem,
    enable_rfm_precompute=enable_rfm_precompute,
    enable_simd=enable_simd,
    OMP_collapse=OMP_collapse,
)
swm2sh.register_CFunction_spin_weight_minus2_sph_harmonics(maximum_l=8)
print(f"Section 1 finished at {time.time() - start_time:.4f} seconds")

if __name__ == "__main__":
    pcg.do_parallel_codegen()

print(f"Section 2 finished at {time.time() - start_time:.4f} seconds")

cbc.CurviBoundaryConditions_register_C_functions(
    CoordSystem, radiation_BC_fd_order=radiation_BC_fd_order
)
print(f"Section 3 finished at {time.time() - start_time:.4f} seconds")
rhs_string = """
Ricci_eval(commondata, params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);
rhs_eval(commondata, params, rfmstruct, auxevol_gfs, RK_INPUT_GFS, RK_OUTPUT_GFS);
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
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
  enforce_detgammabar_equals_detgammahat(commondata, params, rfmstruct, RK_OUTPUT_GFS);""",
    enable_rfm_precompute=enable_rfm_precompute,
    enable_curviBCs=True,
)
print(f"Section 4 finished at {time.time() - start_time:.4f} seconds")
xxCart.register_CFunction__Cart_to_xx_and_nearest_i0i1i2(CoordSystem)
xxCart.register_CFunction_xx_to_Cart(CoordSystem)
progress.register_CFunction_progress_indicator()

#########################################################
# STEP 3: Generate header files, register C functions and
#         command line parameters, set up boundary conditions,
#         and create a Makefile for this project.
#         Project is output to project/[project_name]/
if CoordSystem == "SinhSpherical":
    par.adjust_CodeParam_default("SINHW", 0.4)
par.adjust_CodeParam_default("eta", GammaDriving_eta)

CPs.write_CodeParameters_h_files(project_dir=project_dir)
CPs.register_CFunctions_params_commondata_struct_set_to_default()
cmdpar.generate_default_parfile(project_dir=project_dir, project_name=project_name)
cmdpar.register_CFunction_cmdline_input_and_parfile_parser(
    project_name=project_name, cmdline_inputs=["convergence_factor"]
)
Bdefines_h.output_BHaH_defines_h(
    project_dir=project_dir,
    enable_simd=enable_simd,
    fin_NGHOSTS_add_one_for_upwinding=True,
    MoL_method=MoL_method,
    CoordSystem=CoordSystem,
)
main.register_CFunction_main_c(
    initial_data_desc=IDtype,
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

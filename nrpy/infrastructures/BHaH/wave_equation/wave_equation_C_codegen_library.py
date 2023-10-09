"""
Library of C functions for solving the wave equation in curvilinear coordinates, using a reference-metric formalism.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Union, cast
from inspect import currentframe as cfr
from types import FrameType as FT
from pathlib import Path

import nrpy.grid as gri
import nrpy.params as par
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc

import nrpy.helpers.parallel_codegen as pcg

from nrpy.equations.wave_equation.WaveEquationCurvilinear_RHSs import (
    WaveEquationCurvilinear_RHSs,
)
from nrpy.equations.wave_equation.InitialData import InitialData
import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_exact_solution_single_point(
    CoordSystem: str,
    WaveType: str = "SphericalGaussian",
    default_sigma: float = 3.0,
    default_k0: float = 1.0,
    default_k1: float = 1.0,
    default_k2: float = 1.0,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for the exact solution at a single point.

    :param CoordSystem: The coordinate system to use in setting up initial data.
    :param WaveType: The type of wave: SphericalGaussian or PlaneWave
    :param default_sigma: The default value for the Gaussian width (sigma).
    :param default_k0: The default value for the plane wave wavenumber k in the x-direction.
    :param default_k1: The default value for the plane wave wavenumber k in the y-direction.
    :param default_k2: The default value for the plane wave wavenumber k in the z-direction.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Populate uu_ID, vv_ID
    ID = InitialData(
        CoordSystem=CoordSystem,
        WaveType=WaveType,
        default_sigma=default_sigma,
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
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_initial_data(
    OMP_collapse: int, enable_checkpointing: bool = False
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the initial data function for the wave equation with specific parameters.

    :param enable_checkpointing: Attempt to read from a checkpoint file before generating initial data.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""Set initial data to params.time==0 corresponds to the initial data."""
    c_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    uu_gf_memaccess = gri.BHaHGridFunction.access_gf("uu")
    vv_gf_memaccess = gri.BHaHGridFunction.access_gf("vv")
    body = ""
    if enable_checkpointing:
        body += """// Attempt to read checkpoint file. If it doesn't exist, then continue. Otherwise return.
if( read_checkpoint(commondata, griddata) ) return;
"""
    body += r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
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


def register_CFunction_diagnostics(
    default_diagnostics_out_every: float,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the simulation diagnostics function for the wave equation with specific parameters.

    :param default_diagnostics_out_every: The default frequency for diagnostics output.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _ = par.CodeParameter(
        "REAL",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = r"""Diagnostics."""
    c_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""  const REAL currtime = commondata->time, currdt = commondata->dt, outevery = commondata->diagnostics_output_every;
  // Explanation of the if() below:
  // Step 1: round(currtime / outevery) rounds to the nearest integer multiple of currtime.
  // Step 2: Multiplying by outevery yields the exact time we should output again, t_out.
  // Step 3: If fabs(t_out - currtime) < 0.5 * currdt, then currtime is as close to t_out as possible!
  if (fabs(round(currtime / outevery) * outevery - currtime) < 0.5 * currdt) {
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
        sprintf(filename, "out0d-%s-conv_factor%.2f.txt", CoordSystemName, convergence_factor);
        FILE *outfile;
        if (nn == 0)
          outfile = fopen(filename, "w");
        else
          outfile = fopen(filename, "a");
        if (outfile == NULL) {
          fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
          exit(1);
        }
        const int i0_mid = Nxx_plus_2NGHOSTS0 / 2;
        const int i1_mid = Nxx_plus_2NGHOSTS1 / 2;
        const int i2_mid = Nxx_plus_2NGHOSTS2 / 2;
        const int center_of_grid_idx = IDX3(i0_mid, i1_mid, i2_mid);

        const REAL num_soln_at_center_UUGF = y_n_gfs[IDX4pt(UUGF, center_of_grid_idx)];
        const REAL num_soln_at_center_VVGF = y_n_gfs[IDX4pt(VVGF, center_of_grid_idx)];
        REAL exact_soln_at_center_UUGF, exact_soln_at_center_VVGF;
        exact_solution_single_point(commondata, params, xx[0][i0_mid], xx[1][i1_mid], xx[2][i2_mid], &exact_soln_at_center_UUGF, &exact_soln_at_center_VVGF);

        fprintf(outfile, "%e %e %e %e %e\n", time, fabs(fabs(num_soln_at_center_UUGF - exact_soln_at_center_UUGF) / exact_soln_at_center_UUGF),
                fabs(fabs(num_soln_at_center_VVGF - exact_soln_at_center_VVGF) / (1e-16 + exact_soln_at_center_VVGF)), num_soln_at_center_UUGF, exact_soln_at_center_UUGF);

        fclose(outfile);
      }

      // 1D output
      {
        char filename[256];
        sprintf(filename, "out1d-%s-conv_factor%.2f-t%.2f.txt", CoordSystemName, convergence_factor, time);
        FILE *outfile;
        outfile = fopen(filename, "w");
        if (outfile == NULL) {
          fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
          exit(1);
        }

        for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
          const int i1 = Nxx_plus_2NGHOSTS1 / 2;
          const int i2 = Nxx_plus_2NGHOSTS2 / 2;
          const int idx3 = IDX3(i0, i1, i2);

          const REAL num_soln_UUGF = y_n_gfs[IDX4pt(UUGF, idx3)];
          const REAL num_soln_VVGF = y_n_gfs[IDX4pt(VVGF, idx3)];
          REAL exact_soln_UUGF, exact_soln_VVGF, xCart[3], rr;
          exact_solution_single_point(commondata, params, xx[0][i0], xx[1][i1], xx[2][i2], &exact_soln_UUGF, &exact_soln_VVGF);
          xx_to_Cart(commondata, params, xx,i0,i1,i2, xCart);
          rr = sqrt(xCart[0]*xCart[0] + xCart[1]*xCart[1] + xCart[2]*xCart[2]);

          fprintf(outfile, "%e %e %e %e %e\n", rr, fabs(fabs(num_soln_UUGF - exact_soln_UUGF) / exact_soln_UUGF), fabs(fabs(num_soln_VVGF - exact_soln_VVGF) / (1e-16 + exact_soln_VVGF)),
                  num_soln_UUGF, exact_soln_UUGF);
        }
        fclose(outfile);
      }
    }
  }
  progress_indicator(commondata, griddata);
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
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


def register_CFunction_rhs_eval(
    CoordSystem: str,
    WaveType: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the right-hand side (RHS) evaluation function for the wave equation.

    This function sets the right-hand side of the wave equation according to the
    selected coordinate system and specified parameters.

    :param CoordSystem: The coordinate system.
    :param WaveType: The type of wave.
    :param enable_rfm_precompute: Whether or not to enable reference metric precomputation.
    :param enable_simd: Whether or not to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    OMP_collapse = 1
    if "Spherical" in CoordSystem and WaveType == "SphericalGaussian":
        par.set_parval_from_str("symmetry_axes", "12")
        OMP_collapse = 2  # about 2x faster

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = r"""Set RHSs for wave equation."""
    c_type = "void"
    name = "rhs_eval"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], const REAL *restrict in_gfs, REAL *restrict rhs_gfs"
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate uu_rhs, vv_rhs
    rhs = WaveEquationCurvilinear_RHSs(CoordSystem, enable_rfm_precompute)
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
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        c_type=c_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

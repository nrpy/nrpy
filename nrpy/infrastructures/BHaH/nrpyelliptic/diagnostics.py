"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)


# Define function to compute the l^2 of a gridfunction
def register_CFunction_compute_L2_norm_of_gridfunction(
    CoordSystem: str,
) -> None:
    """
    Register function to compute l2-norm of a gridfunction assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.
    """
    includes = ["BHaH_defines.h"]
    desc = "Compute l2-norm of a gridfunction assuming a single grid."
    cfunc_type = "REAL"
    name = "compute_L2_norm_of_gridfunction"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL integration_radius, const int gf_index, const REAL *restrict in_gf"""

    rfm = refmetric.reference_metric[CoordSystem]

    fp_type = par.parval_from_str("fp_type")
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"
    loop_body = ccg.c_codegen(
        [
            rfm.xxSph[0],
            rfm.detgammahat,
        ],
        [
            "const DOUBLE r",
            "const DOUBLE sqrtdetgamma",
        ],
        include_braces=False,
        fp_type_alias=fp_type_alias,
    )

    loop_body += r"""
if(r < integration_radius) {
  const DOUBLE gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
  const DOUBLE dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r < integration_radius)
"""
    body = r"""
  // Unpack grid parameters assuming a single grid
  const int grid = 0;
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Define reference metric grid
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Set summation variables to compute l2-norm
  DOUBLE squared_sum = 0.0;
  DOUBLE volume_sum  = 0.0;
"""

    body += lp.simple_loop(
        loop_body="\n" + loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)",
    )

    body += r"""
  // Compute and output the log of the l2-norm.
  return log10(1e-16 + sqrt(squared_sum / volume_sum));  // 1e-16 + ... avoids log10(0)
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=body,
    )


# Define function to compute residual the solution
def register_CFunction_compute_residual_all_points(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_simd: bool,
    OMP_collapse: int,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the residual evaluation function.

    This function sets the residual of the Hamiltonian constraint in the hyperbolic
    relaxation equation according to the selected coordinate system and specified
    parameters.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_simd: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    if enable_simd:
        includes += [str(Path("intrinsics") / "simd_intrinsics.h")]
    desc = r"""Compute residual of the Hamiltonian constraint for the hyperbolic relaxation equation."""
    cfunc_type = "void"
    name = "compute_residual_all_points"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                REAL *restrict xx[3], const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict aux_gfs"""
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )
    # Populate residual_H
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute)
    body = lp.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.residual],
            [
                gri.BHaHGridFunction.access_gf("residual_H", gf_array_name="aux_gfs"),
            ],
            enable_fd_codegen=True,
            enable_simd=enable_simd,
        ),
        loop_region="interior",
        enable_intrinsics=enable_simd,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        body=body,
        enable_simd=enable_simd,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define diagnostics function
def register_CFunction_diagnostics(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    default_diagnostics_out_every: int,
    enable_progress_indicator: bool = False,
    axis_filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-n-%08d.txt",
        "nn",
    ),
    plane_filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-n-%08d.txt",
        "nn",
    ),
    out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for simulation diagnostics.

    :param CoordSystem: Coordinate system used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param enable_progress_indicator: Whether to enable the progress indicator.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    _ = par.CodeParameter(
        "int",
        __name__,
        "diagnostics_output_every",
        default_diagnostics_out_every,
        commondata=True,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Diagnostics."
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "numUU"): "y_n_gfs[IDX4pt(UUGF, idx3)]",
            ("REAL", "log10ResidualH"): "log10(fabs(diagnostic_output_gfs[IDX4pt(RESIDUAL_HGF, idx3)] + 1e-16))",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on
    for axis in ["y", "z"]:
        out012d.register_CFunction_diagnostics_nearest_1d_axis(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            axis=axis,
            filename_tuple=axis_filename_tuple,
        )
    for plane in ["xy", "yz"]:
        out012d.register_CFunction_diagnostics_nearest_2d_plane(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            plane=plane,
            filename_tuple=plane_filename_tuple,
        )

    body = r"""  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
  const int grid = 0;

  // Set gridfunctions aliases
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
"""
    if enable_rfm_precompute:
        body += r"""  // Set rfm_struct
  const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""
    else:
        body += r"""
  // Set xx
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  compute_residual_all_points(commondata, params, xx, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""

    body += r"""
  // Set integration radius for l2-norm computation
  const REAL integration_radius = 1000;

  // Compute l2-norm of Hamiltonian constraint violation
  const REAL residual_H = compute_L2_norm_of_gridfunction(commondata, griddata, integration_radius, RESIDUAL_HGF, diagnostic_output_gfs);

  // Update residual to be used in stop condition
  commondata->log10_current_residual = residual_H;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (nn == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", nn, time, residual_H);
    fclose(outfile);
  }

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;
  if (n_step % outevery == 0) {
    // Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
        xx[ww] = griddata[grid].xx[ww];

    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &griddata[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &griddata[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &griddata[grid].gridfuncs);
  }
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


# Define function to evaluate stop conditions
def register_CFunction_check_stop_conditions() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to evaluate stop conditions.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h"]
    desc = "Evaluate stop conditions."
    cfunc_type = "void"
    name = "check_stop_conditions"
    params = (
        """commondata_struct *restrict commondata, griddata_struct *restrict griddata"""
    )

    # Register a parameter to stop hyperbolic relaxation
    _stop_relaxation = par.register_CodeParameter(
        "bool", __name__, "stop_relaxation", False, commondata=True
    )

    # Register parameter that sets the total number of relaxation steps
    _nn_max = par.register_CodeParameter("int", __name__, "nn_max", 0, commondata=True)

    # Register parameter that sets the tolerance for log of residual
    _log10_residual_tolerance = par.register_CodeParameter(
        "REAL", __name__, "log10_residual_tolerance", -15.8, commondata=True
    )

    # Register parameter that sets log of residual to be updated at every time step
    _log10_current_residual = par.register_CodeParameter(
        "REAL", __name__, "log10_current_residual", 1.0, commondata=True
    )

    body = r"""  // Since this version of NRPyElliptic is unigrid, we simply set the grid index to 0
  const int grid = 0;

  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"

  // Check if total number of iteration steps has been reached
  if ((nn >= nn_max) || (log10_current_residual < log10_residual_tolerance)){
    printf("\nExiting main loop after %8d iterations\n", nn);
    printf("The tolerance for the logarithmic residual is %.8e\n", log10_residual_tolerance);
    printf("Exiting relaxation with logarithmic residual of %.8e\n", log10_current_residual);
    commondata->stop_relaxation = true;
  }
"""

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func="",
        name=name,
        params=params,
        include_CodeParameters_h=False,  # set_CodeParameters.h is manually included after the declaration of params_struct *restrict params
        body=body,
    )
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

"""
C function for hyperbolic relaxation diagnostics.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.infrastructures import BHaH


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

    Doctest:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name="diagnostics"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    _ = register_CFunction_diagnostics("Cartesian", True, 100)
    ...    generated_str = cfc.CFunction_dict[f'{name}'].full_function
    ...    validation_desc = f"{name}__{parallelization}".replace(" ", "_")
    ...    validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
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

    parallelization = par.parval_from_str("parallelization")
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Diagnostics."
    cfunc_type = "void"
    name = "diagnostics"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    params += (
        ", griddata_struct *restrict griddata_host"
        if parallelization in ["cuda"]
        else ""
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
        BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices.register_CFunction_diagnostics_nearest_1d_axis(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            axis=axis,
            filename_tuple=axis_filename_tuple,
        )
    for plane in ["xy", "yz"]:
        BHaH.diagnostics.output_0d_1d_2d_nearest_gridpoint_slices.register_CFunction_diagnostics_nearest_2d_plane(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            plane=plane,
            filename_tuple=plane_filename_tuple,
        )

    body = r"""  // Output progress to stderr
  progress_indicator(commondata, griddata);

  // Grid data output
  const int n_step = commondata->nn, outevery = commondata->diagnostics_output_every;
  const REAL time = commondata->time;

  REAL global_norm = -1e9;
  for(int grid = 0; grid < commondata->NUMGRIDS; ++grid) {

  // Set gridfunctions aliases
  REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata[grid].gridfuncs.diagnostic_output_gfs;
  MAYBE_UNUSED REAL *restrict diagnostic_output_gfs2 = griddata[grid].gridfuncs.diagnostic_output_gfs2;
  // Set params
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
"""
    if parallelization in ["cuda"]:
        body += r"""
  REAL *restrict host_y_n_gfs = griddata_host[grid].gridfuncs.y_n_gfs;
  REAL *restrict host_diag_gfs = griddata_host[grid].gridfuncs.diagnostic_output_gfs;
  if (n_step % outevery == 0) {
    size_t streamid = params->grid_idx % NUM_STREAMS;
    cpyDevicetoHost__gf(commondata, params, host_y_n_gfs, y_n_gfs, UUGF, UUGF, streamid);
  }
"""
    if enable_rfm_precompute:
        body += r"""  // Set rfm_struct
  const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  residual_H_compute_all_points(commondata, params, rfmstruct, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""
    else:
        body += r"""
  // Set xx
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];

  // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
  residual_H_compute_all_points(commondata, params, xx, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);
"""

    if parallelization in ["cuda"]:
        body += r"""if (n_step % outevery == 0) {
    BHAH_DEVICE_SYNC();
    size_t streamid = params->grid_idx % NUM_STREAMS;
    cpyDevicetoHost__gf(commondata, params, host_diag_gfs, diagnostic_output_gfs, RESIDUAL_HGF, RESIDUAL_HGF, streamid);
  }"""

    body += r"""
  // Set integration radius for l2-norm computation
  const REAL integration_radius = 1000;

  // Compute l2-norm of Hamiltonian constraint violation
  REAL residual_H;
  log10_L2norm_gf(commondata, &griddata[grid].params, griddata[grid].xx, integration_radius, RESIDUAL_HGF, &residual_H, diagnostic_output_gfs);
  global_norm = MAX(global_norm, residual_H);
  } // END for(grid=0; grid<commondata->NUMGRIDS; ++grid)

  // Update residual to be used in stop condition
  commondata->log10_current_residual = global_norm;

  // Output l2-norm of Hamiltonian constraint violation to file
  {
    char filename[256];
    sprintf(filename, "residual_l2_norm.txt");
    FILE *outfile = (n_step == 0) ? fopen(filename, "w") : fopen(filename, "a");
    if (!outfile) {
      fprintf(stderr, "Error: Cannot open file %s for writing.\n", filename);
      exit(1);
    }
    fprintf(outfile, "%6d %10.4e %.17e\n", n_step, time, global_norm);
    fclose(outfile);
  }

  if (n_step % outevery == 0) {
    // Only consider a single grid for now.
    const int grid = 0;
    params_struct *restrict params = &griddata[grid].params;
""".replace(
        "diagnostic_output_gfs)",
        (
            "diagnostic_output_gfs, diagnostic_output_gfs2)"
            if parallelization in ["cuda"]
            else "diagnostic_output_gfs)"
        ),
    )

    host_griddata = "griddata_host" if parallelization in ["cuda"] else "griddata"
    body += rf"""// Set reference metric grid xx
    REAL *restrict xx[3];
    for (int ww = 0; ww < 3; ww++)
        xx[ww] = {host_griddata}[grid].xx[ww];
"""
    if parallelization in ["cuda"]:
        body += r"""
    // Sync with device when relevant
    BHAH_DEVICE_SYNC();
"""
    body += rf"""
    // 1D output
    diagnostics_nearest_1d_y_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_1d_z_axis(commondata, params, xx, &{host_griddata}[grid].gridfuncs);

    // 2D output
    diagnostics_nearest_2d_xy_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
    diagnostics_nearest_2d_yz_plane(commondata, params, xx, &{host_griddata}[grid].gridfuncs);
  }}
"""
    if enable_progress_indicator:
        body += "progress_indicator(commondata, griddata);"
    body += r"""
  if (commondata->time + commondata->dt > commondata->t_final)
    printf("\n");
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

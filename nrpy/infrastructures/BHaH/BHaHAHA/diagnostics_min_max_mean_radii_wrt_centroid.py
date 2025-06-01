"""
Register the C function for computing minimum, maximum, and mean radius *from the coordinate centroid* of the horizon.
Needs: coordinate centroid of the horizon.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures.BHaH.BHaHAHA import area


def register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for computing minimum, maximum, and mean radius *from the coordinate centroid* of the horizon.
    Needs: coordinate centroid of the horizon.

    :param enable_fd_functions: Whether to enable finite difference functions, defaults to True.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    DocTests:
    >>> import nrpy.grid as gri
    >>> _ = gri.register_gridfunctions("hh")[0]
    >>> env = register_CFunction_diagnostics_min_max_mean_radii_wrt_centroid()
    Setting up reference_metric[Spherical]...
    Setting up ExpansionFunctionThetaClass[Spherical]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "BHaHAHA apparent horizon diagnostics: Compute Theta L2 and Linfinity norms."
    cfunc_type = "void"
    name = "diagnostics_min_max_mean_radii_wrt_centroid"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  const int grid = 0;
  bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;
  const params_struct *restrict params = &griddata[grid].params;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++)
    xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

  // Set integration weights.
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

  // Compute radii quantities. Mean radius is area-weighted, since physical gridspacing is quite uneven.
  REAL sum_curr_area = 0.0;
  REAL sum_mean_radius = 0.0;
  REAL min_radius_squared = +1e30;
  REAL max_radius_squared = -1e30;
#pragma omp parallel
{
#pragma omp for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
    const MAYBE_UNUSED REAL xx2 = xx[2][i2];
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
        const MAYBE_UNUSED REAL xx1 = xx[1][i1];
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
"""
    body += (
        ccg.c_codegen(
            area.area3(),
            "const REAL area_element",
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )
        + """
#pragma omp critical
          {
            sum_curr_area += area_element * weight1 * weight2;
            const REAL r = in_gfs[IDX4(HHGF, NGHOSTS, i1, i2)];
            const REAL theta = commondata->interp_src_r_theta_phi[1][i1];
            const REAL phi = commondata->interp_src_r_theta_phi[2][i2];
            const REAL xx = r * sin(theta) * cos(phi);
            const REAL yy = r * sin(theta) * sin(phi);
            const REAL zz = r * cos(theta);
            // Radius as measured from AH centroid:
            REAL radius_squared = ((xx - bhahaha_diags->x_centroid_wrt_coord_origin) * (xx - bhahaha_diags->x_centroid_wrt_coord_origin) + //
                                   (yy - bhahaha_diags->y_centroid_wrt_coord_origin) * (yy - bhahaha_diags->y_centroid_wrt_coord_origin) + //
                                   (zz - bhahaha_diags->z_centroid_wrt_coord_origin) * (zz - bhahaha_diags->z_centroid_wrt_coord_origin));
            sum_mean_radius += sqrt(radius_squared) * area_element * weight1 * weight2;
            if (radius_squared < min_radius_squared)
              min_radius_squared = radius_squared;
            if (radius_squared > max_radius_squared)
              max_radius_squared = radius_squared;
          } // END OMP CRITICAL
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END OMP PARALLEL

  // Store commondata diagnostic parameters.
  const REAL curr_area = sum_curr_area * params->dxx1 * params->dxx2;
  bhahaha_diags->min_coord_radius_wrt_centroid = sqrt(min_radius_squared);
  bhahaha_diags->max_coord_radius_wrt_centroid = sqrt(max_radius_squared);
  bhahaha_diags->mean_coord_radius_wrt_centroid = sum_mean_radius * params->dxx1 * params->dxx2 / curr_area;
"""
    )
    cfc.register_CFunction(
        subdirectory="",
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

    results = doctest.testmod()

    if results.failed > 0:
        raise RuntimeError(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)"
        )
    print(f"Doctest passed: All {results.attempted} test(s) passed")

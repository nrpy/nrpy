"""
Register the C function for BHaHAHA apparent horizon diagnostics.
Compute area, centroid location, and Theta (L2 and Linfinity) norms.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.general_relativity.bhahaha.ExpansionFunctionTheta import (
    ExpansionFunctionTheta,
)
from nrpy.infrastructures import BHaH


def register_CFunction_diagnostics_area_centroid_and_Theta_norms(
    CoordSystem: str = "Spherical",
    enable_rfm_precompute: bool = False,
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for BHaHAHA apparent horizon diagnostics: compute area, centroid location, and Theta (L2 and Linfinity) norms.

    :param CoordSystem: The coordinate system to use, defaults to "Spherical".
    :param enable_rfm_precompute: Whether to enable RFM precompute, defaults to True.
    :param enable_fd_functions: Whether to enable finite difference functions, defaults to True.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    DocTests:
    >>> import nrpy.grid as gri
    >>> _ = gri.register_gridfunctions("hh")[0]
    >>> env = register_CFunction_diagnostics_area_centroid_and_Theta_norms()
    Setting up ExpansionFunctionThetaClass[Spherical]...
    Setting up reference_metric[Spherical]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # This is only an internal diagnostic; the diagnostic interesting to users is wrt the
    #    horizon centroid. This diagnostic in general will exhibit a kink the second time
    #    a horizon is found, as the initial guess for the horizon centroid might be a little
    #    bit off from the guess.
    _ = par.register_CodeParameters(
        "REAL",
        __name__,
        ["min_radius_wrt_grid_center", "max_radius_wrt_grid_center"],
        -1.0,
        commondata=True,
        add_to_parfile=False,
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "BHaHAHA apparent horizon diagnostics: compute area, centroid location, and Theta (L2 and Linfinity) norms."
    cfunc_type = "void"
    name = "diagnostics_area_centroid_and_Theta_norms"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    Th = ExpansionFunctionTheta[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]
    body = r"""
  const int grid=0;
  const params_struct *restrict params = &griddata[grid].params;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict xx[3];
  for(int ww=0;ww<3;ww++) xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

  // Set integration weights.
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights, &weight_stencil_size);

  // Compute Linfinity and L2 norms of Theta
  REAL min_radius = 1e10, max_radius = -1e10;
  REAL sum_Theta_squared_for_L2_norm = 0.0;
  REAL sum_curr_area = 0;
  REAL sum_x_centroid = 0, sum_y_centroid = 0, sum_z_centroid = 0;
  REAL max_Theta_squared_for_Linf_norm = -1e30;
#pragma omp parallel
{
#pragma omp for
  for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
    const REAL weight2 = weights[(i2 - NGHOSTS) % weight_stencil_size];
    const REAL xx2 = xx[2][i2];
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const REAL weight1 = weights[(i1 - NGHOSTS) % weight_stencil_size];
        const REAL xx1 = xx[1][i1];
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
"""
    body += (
        ccg.c_codegen(
            [Th.Theta, BHaH.BHaHAHA.area.area3()],
            ["const REAL Theta", "const REAL area_element"],
            enable_fd_codegen=True,
            enable_fd_functions=enable_fd_functions,
        )
        + """
#pragma omp critical
          {
            sum_curr_area += area_element * weight1 * weight2;
            sum_Theta_squared_for_L2_norm += Theta * Theta * area_element * weight1 * weight2;
            {
              const REAL tmp0 = hh * sin(xx1);
              sum_x_centroid += (Cart_originx + tmp0 * cos(xx2)) * area_element * weight1 * weight2;
              sum_y_centroid += (Cart_originy + tmp0 * sin(xx2)) * area_element * weight1 * weight2;
              sum_z_centroid += (Cart_originz + hh * cos(xx1)) * area_element * weight1 * weight2;
            } // END centroid sums
            if (Theta * Theta > max_Theta_squared_for_Linf_norm)
              max_Theta_squared_for_Linf_norm = Theta * Theta;
            if (hh > max_radius)
              max_radius = hh;
            if (hh < min_radius)
              min_radius = hh;
          } // END OMP CRITICAL
        } // END LOOP over i0
      } // END LOOP over i1
    } // END LOOP over i2
  } // END OMP PARALLEL

  // Store diagnostics in commondata->bhahaha_diagnostics struct.
  {

    // {min,max}_radius_wrt_grid_center are strictly internal diagnostics; these will
    //   exhibit a kink the second time a horizon is found, as the initial guess for
    //   the horizon centroid might be a little bit off from the guess.
    // Users should view the radii wrt the AH centroid location, which will be smooth
    //   over time.
    commondata->min_radius_wrt_grid_center = min_radius;
    commondata->max_radius_wrt_grid_center = max_radius;

    bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;

    // Theta was just evaluated at Nxx1*Nxx2 gridpoints; update counter:
    bhahaha_diags->Theta_eval_points_counter += params->Nxx1 * params->Nxx2;

    // Update area, compute irreducible mass.
    bhahaha_diags->area = sum_curr_area * params->dxx1 * params->dxx2;

    // Compute area-weighted norms of Theta, normalized properly by M_scale, as [Theta] ~ 1/length.
    const REAL M_scale = commondata->bhahaha_params_and_data->M_scale;
    bhahaha_diags->Theta_Linf_times_M = M_scale * sqrt(max_Theta_squared_for_Linf_norm);
    bhahaha_diags->Theta_L2_times_M = M_scale * sqrt(sum_Theta_squared_for_L2_norm * params->dxx1 * params->dxx2 / bhahaha_diags->area);

    bhahaha_diags->x_centroid_wrt_coord_origin = sum_x_centroid * params->dxx1 * params->dxx2 / bhahaha_diags->area;
    bhahaha_diags->y_centroid_wrt_coord_origin = sum_y_centroid * params->dxx1 * params->dxx2 / bhahaha_diags->area;
    bhahaha_diags->z_centroid_wrt_coord_origin = sum_z_centroid * params->dxx1 * params->dxx2 / bhahaha_diags->area;
  } // END store diagnostics in commondata->bhahaha_diagnostics struct.
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

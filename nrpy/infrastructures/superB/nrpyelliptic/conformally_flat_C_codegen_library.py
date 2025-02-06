"""
Library of C functions for solving the hyperbolic relaxation equation in curvilinear coordinates, using a reference-metric formalism.

Authors: Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Nishita Jadoo; njadoo **at** uidaho **dot* edu

"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, Tuple, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.infrastructures.superB.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.infrastructures.BHaH import griddata_commondata


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
    cfunc_type = "void"
    name = "compute_L2_norm_of_gridfunction"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL integration_radius, const int gf_index, const REAL *restrict in_gf, REAL localsums_for_residualH[2]"""

    rfm = refmetric.reference_metric[CoordSystem]

    loop_body = ccg.c_codegen(
        [
            rfm.xxSph[0],
            rfm.detgammahat,
        ],
        [
            "const REAL r",
            "const REAL sqrtdetgamma",
        ],
        include_braces=False,
    )

    loop_body += r"""
if(r < integration_radius) {
  const REAL gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
  const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
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
  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;

"""

    body += lp.simple_loop(
        loop_body="\n" + loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)",
    )

    body += r"""
  localsums_for_residualH[0] = squared_sum;
	localsums_for_residualH[1] = volume_sum;
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


# Define function to compute the l^2 of a gridfunction between 2 radii
def register_CFunction_compute_L2_norm_of_gridfunction_between_r1_r2(
    CoordSystem: str,
) -> None:
    """
    Register function to compute l2-norm of a gridfunction between 2 radii assuming a single grid.

    Note that parallel codegen is disabled for this function, as it sometimes causes a
    multiprocess race condition on Python 3.6.7

    :param CoordSystem: the rfm coordinate system.
    """
    includes = ["BHaH_defines.h"]
    desc = "Compute l2-norm of a gridfunction between 2 radii  assuming a single grid."
    cfunc_type = "void"
    name = "compute_L2_norm_of_gridfunction_between_r1_r2"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL integration_radius1, const REAL integration_radius2, const int gf_index, const REAL *restrict in_gf, REAL localsums[2]"""

    rfm = refmetric.reference_metric[CoordSystem]

    loop_body = ccg.c_codegen(
        [
            rfm.xxSph[0],
            rfm.detgammahat,
        ],
        [
            "const REAL r",
            "const REAL sqrtdetgamma",
        ],
        include_braces=False,
    )

    loop_body += r"""
if(r > integration_radius1 && r < integration_radius2) {
  const REAL gf_of_x = in_gf[IDX4(gf_index, i0, i1, i2)];
  const REAL dV = sqrtdetgamma * dxx0 * dxx1 * dxx2;
  squared_sum += gf_of_x * gf_of_x * dV;
  volume_sum  += dV;
} // END if(r > integration_radius1 && r < integration_radius2)
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
  REAL squared_sum = 0.0;
  REAL volume_sum  = 0.0;

"""

    body += lp.simple_loop(
        loop_body="\n" + loop_body,
        read_xxs=True,
        loop_region="interior",
        OMP_custom_pragma=r"#pragma omp parallel for reduction(+:squared_sum,volume_sum)",
    )

    body += r"""
  localsums[0] = squared_sum;
	localsums[1] = volume_sum;
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


# Define diagnostics function
def register_CFunction_diagnostics(
    CoordSystem: str,
    default_diagnostics_out_every: int,
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
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
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
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata_chare, griddata_struct *restrict griddata, Ck::IO::Session token, const int which_output, const int grid, const int chare_index[3], REAL localsums_for_residualH[2]"

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
        )
        out012d.register_CFunction_diagnostics_set_up_nearest_1d_axis(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=axis_filename_tuple,
            axis=axis,
        )
    for plane in ["xy", "yz"]:
        out012d.register_CFunction_diagnostics_nearest_2d_plane(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            plane=plane,
        )
        out012d.register_CFunction_diagnostics_set_up_nearest_2d_plane(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=plane_filename_tuple,
            plane=plane,
        )

    body = r"""  // Unpack griddata and griddata_chare
  const int nn = commondata->nn;
  const params_struct *restrict params = &griddata[grid].params;
  const params_struct *restrict params_chare = &griddata_chare[grid].params;
  const rfm_struct *restrict rfmstruct_chare = griddata_chare[grid].rfmstruct;
  const charecomm_struct *restrict charecommstruct = &griddata_chare[grid].charecommstruct;
  const REAL *restrict y_n_gfs = griddata_chare[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata_chare[grid].gridfuncs.auxevol_gfs;
  REAL *restrict diagnostic_output_gfs = griddata_chare[grid].gridfuncs.diagnostic_output_gfs;
  REAL *restrict xx_chare[3];
  {
    for (int ww = 0; ww < 3; ww++)
      xx_chare[ww] = griddata_chare[grid].xx[ww];
  }
  REAL *restrict xx[3];
  {
    for (int ww = 0; ww < 3; ww++)
      xx[ww] = griddata[grid].xx[ww];
  }
  const int Nchare0 = commondata->Nchare0;
  const int Nchare1 = commondata->Nchare1;
  const int Nchare2 = commondata->Nchare2;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  if (which_output == OUTPUT_RESIDUAL) {
    // Compute Hamiltonian constraint violation and store it at diagnostic_output_gfs
    compute_residual_all_points(commondata, params_chare, rfmstruct_chare, auxevol_gfs, y_n_gfs, diagnostic_output_gfs);

    // Set integration radius for l2-norm computation
    const REAL integration_radius = 1000;

    // Compute local sums for l2-norm of Hamiltonian constraint violation
    compute_L2_norm_of_gridfunction(commondata, griddata_chare, integration_radius, RESIDUAL_HGF, diagnostic_output_gfs, localsums_for_residualH);

"""
    body += r"""
} else {
    const int num_diagnostic_1d_y_pts = griddata_chare[grid].diagnosticstruct.num_diagnostic_1d_y_pts;
    const int num_diagnostic_1d_z_pts = griddata_chare[grid].diagnosticstruct.num_diagnostic_1d_z_pts;
    const int num_diagnostic_2d_xy_pts = griddata_chare[grid].diagnosticstruct.num_diagnostic_2d_xy_pts;
    const int num_diagnostic_2d_yz_pts = griddata_chare[grid].diagnosticstruct.num_diagnostic_2d_yz_pts;


    const bool write_diagnostics = (which_output == OUTPUT_0D) ||
              (num_diagnostic_1d_y_pts > 0) ||
              (num_diagnostic_1d_z_pts > 0) ||
              (num_diagnostic_2d_xy_pts > 0) ||
              (num_diagnostic_2d_yz_pts > 0);

    if (write_diagnostics) {
      // 1D and 2D outputs
      if (which_output == OUTPUT_1D_Y) {
        if (num_diagnostic_1d_y_pts > 0) {
        diagnostics_nearest_1d_y_axis(commondata, params_chare, xx_chare, &griddata_chare[grid].gridfuncs, &griddata_chare[grid].diagnosticstruct, token);
        }
      } else if (which_output == OUTPUT_1D_Z) {
        if (num_diagnostic_1d_z_pts > 0) {
        diagnostics_nearest_1d_z_axis(commondata, params_chare, xx_chare, &griddata_chare[grid].gridfuncs, &griddata_chare[grid].diagnosticstruct, token);
        }
      } else if (which_output == OUTPUT_2D_XY) {
        if (num_diagnostic_2d_xy_pts > 0) {
        diagnostics_nearest_2d_xy_plane(commondata, params_chare, xx_chare, &griddata_chare[grid].gridfuncs, &griddata_chare[grid].diagnosticstruct, token);
        }
      } else if (which_output == OUTPUT_2D_YZ) {
        if (num_diagnostic_2d_yz_pts > 0) {
        diagnostics_nearest_2d_yz_plane(commondata, params_chare, xx_chare, &griddata_chare[grid].gridfuncs, &griddata_chare[grid].diagnosticstruct, token);
        }
      }
    }
  }
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

    # Register diagnostic_struct's contribution to griddata_struct:
    griddata_commondata.register_griddata_commondata(
        __name__,
        "diagnostic_struct diagnosticstruct",
        "store indices of 1d and 2d diagnostic points, the offset in the output file, etc",
    )

    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())

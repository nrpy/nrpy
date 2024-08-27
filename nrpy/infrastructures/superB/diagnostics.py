"""
C functions for diagnostics for the superB infrastructure.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.superB.output_0d_1d_2d_nearest_gridpoint_slices as out012d
import nrpy.params as par
from nrpy.infrastructures.BHaH import griddata_commondata


def register_CFunction_psi4_diagnostics_set_up() -> Union[None, pcg.NRPyEnv_type]:
    r"""
    Register C function for setting up diagnostic struct.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = "Setup psi4_diagnostics."
    cfunc_type = "void"
    name = "psi4_diagnostics_set_up"
    params = "commondata_struct *restrict commondata, const params_struct *restrict params, const params_struct *restrict params_chare, const charecomm_struct *restrict charecommstruct, REAL *restrict xx[3], const int chare_index[3], diagnostic_struct *restrict diagnosticstruct"

    body = r"""

  const int Nchare0 = commondata->Nchare0;
  const int Nchare1 = commondata->Nchare1;
  const int Nchare2 = commondata->Nchare2;
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  // Do psi4 output, but only if the grid is spherical-like.
  if (strstr(params_chare->CoordSystemName, "Spherical") != NULL) {

    // Adjusted to match Tutorial-Start_to_Finish-BSSNCurvilinear-Two_BHs_Collide-Psi4.ipynb
    const int psi4_spinweightm2_sph_harmonics_max_l = 2;
#define num_of_R_exts 24
    const REAL list_of_R_exts[num_of_R_exts] = {10.0, 20.0, 21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0,  30.0,
                                                31.0, 32.0, 33.0, 35.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 150.0};
    // Find which R_exts lie in the chares's grid
    int count_num_of_R_exts_chare = 0;
    for (int which_R_ext = 0; which_R_ext < num_of_R_exts; which_R_ext++) {
      const REAL R_ext = list_of_R_exts[which_R_ext];
      const REAL xCart_R_ext[3] = {R_ext, 0.0, 0.0};
      int Cart_to_i0i1i2[3];
      REAL closest_xx[3];
      Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart_R_ext, closest_xx, Cart_to_i0i1i2);

      // To lie in this section only chare index 0 needs to be same as this chare
      const int globalidx3 = IDX3GENERAL(Cart_to_i0i1i2[0], Cart_to_i0i1i2[1], Cart_to_i0i1i2[2], Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1);
      const int chareidx3_of_pt = charecommstruct->globalidx3pt_to_chareidx3[globalidx3];
      int chare_index0_pt;
      int chare_index1_pt;
      int chare_index2_pt;
      REVERSE_IDX3GENERAL(chareidx3_of_pt, Nchare0, Nchare1, chare_index0_pt, chare_index1_pt, chare_index2_pt);
      if (chare_index0_pt == chare_index[0]) {
        // count num R_exts_chare
        count_num_of_R_exts_chare++;
      }
    }

    // Declare list_of_R_exts_chare with the correct size after counting
    diagnosticstruct->num_of_R_exts_chare = count_num_of_R_exts_chare;
    diagnosticstruct->list_of_R_exts_chare = (REAL *restrict)malloc(sizeof(REAL) * count_num_of_R_exts_chare);
    diagnosticstruct->psi4_spinweightm2_sph_harmonics_max_l = psi4_spinweightm2_sph_harmonics_max_l;

    count_num_of_R_exts_chare = 0; // Reset the counter
    for (int which_R_ext = 0; which_R_ext < num_of_R_exts; which_R_ext++) {
      const REAL R_ext = list_of_R_exts[which_R_ext];
      const REAL xCart_R_ext[3] = {R_ext, 0.0, 0.0};
      int Cart_to_i0i1i2[3];
      REAL closest_xx[3];
      Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart_R_ext, closest_xx, Cart_to_i0i1i2);
      const int globalidx3 = IDX3GENERAL(Cart_to_i0i1i2[0], Cart_to_i0i1i2[1], Cart_to_i0i1i2[2], Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1);

      // To lie in this section only chare index 0 needs to be same as this chare
      const int chareidx3_of_pt = charecommstruct->globalidx3pt_to_chareidx3[globalidx3];
      int chare_index0_pt;
      int chare_index1_pt;
      int chare_index2_pt;
      REVERSE_IDX3GENERAL(chareidx3_of_pt, Nchare0, Nchare1, chare_index0_pt, chare_index1_pt, chare_index2_pt);
      if (chare_index0_pt == chare_index[0]) {
        // Set list_of_R_exts_chare
        diagnosticstruct->list_of_R_exts_chare[count_num_of_R_exts_chare] = R_ext;
        count_num_of_R_exts_chare++;
      }
    }

    // Allocate memory localsums_for_psi4_decomp
    const int size_psi4_decomp =
        diagnosticstruct->num_of_R_exts_chare * (psi4_spinweightm2_sph_harmonics_max_l - 1) * ((2 * psi4_spinweightm2_sph_harmonics_max_l) + 1) * 2;
    diagnosticstruct->localsums_for_psi4_decomp = (REAL *restrict)malloc(sizeof(REAL) * size_psi4_decomp);
    diagnosticstruct->globalsums_for_psi4_decomp = (REAL *restrict)malloc(sizeof(REAL) * size_psi4_decomp);
    diagnosticstruct->length_localsums_for_psi4_decomp = size_psi4_decomp;
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
    return cast(pcg.NRPyEnv_type, pcg.NRPyEnv())


def register_CFunction_diagnostics(
    list_of_CoordSystems: List[str],
    default_diagnostics_out_every: float,
    grid_center_filename_tuple: Tuple[str, str] = (
        "out0d-conv_factor%.2f.txt",
        "convergence_factor",
    ),
    axis_filename_tuple: Tuple[str, str] = (
        "out1d-AXIS-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    plane_filename_tuple: Tuple[str, str] = (
        "out2d-PLANE-%s-conv_factor%.2f-t%08.2f.txt",
        "CoordSystemName, convergence_factor, time",
    ),
    out_quantities_dict: Union[str, Dict[Tuple[str, str], str]] = "default",
    enable_BSSN_diagnostics: bool = True,
    enable_psi4_diagnostics: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function for simulation diagnostics.

    :param list_of_CoordSystems: Lists of unique CoordSystems used.
    :param default_diagnostics_out_every: Specifies the default diagnostics output frequency.
    :param grid_center_filename_tuple: Tuple containing filename and variables for grid center output.
    :param axis_filename_tuple: Tuple containing filename and variables for axis output.
    :param plane_filename_tuple: Tuple containing filename and variables for plane output.
    :param out_quantities_dict: Dictionary or string specifying output quantities.
    :param enable_BSSN_diagnostics: Whether or not to enable BSSN constraint violation diagnostics.
    :param enable_psi4_diagnostics: Whether or not to enable psi4 diagnostics.

    :return: None if in registration phase, else the updated NRPy environment.
    :raises TypeError: If `out_quantities_dict` is not a dictionary and not set to "default".
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

    # fmt: off
    if out_quantities_dict == "default":
        out_quantities_dict = {
            ("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))",
            ("REAL", "log10sqrtM2L"): "log10(sqrt(diagnostic_output_gfs[IDX4pt(MSQUAREDGF, idx3)]) + 1e-16)",
            ("REAL", "cfL"): "y_n_gfs[IDX4pt(CFGF, idx3)]",
            ("REAL", "alphaL"): "y_n_gfs[IDX4pt(ALPHAGF, idx3)]",
            ("REAL", "trKL"): "y_n_gfs[IDX4pt(TRKGF, idx3)]",
        }
    if not isinstance(out_quantities_dict, dict):
        raise TypeError(f"out_quantities_dict was initialized to {out_quantities_dict}, which is not a dictionary!")
    # fmt: on

    for CoordSystem in list_of_CoordSystems:
        out012d.register_CFunction_diagnostics_nearest_grid_center(
            CoordSystem=CoordSystem,
            out_quantities_dict=out_quantities_dict,
            filename_tuple=grid_center_filename_tuple,
        )
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

    if enable_psi4_diagnostics:
        register_CFunction_psi4_diagnostics_set_up()

    desc = r"""Diagnostics."""
    cfunc_type = "void"
    name = "diagnostics"
    params = "commondata_struct *restrict commondata, griddata_struct *restrict griddata_chare, griddata_struct *restrict griddata, Ck::IO::Session token, const int which_output, const int grid, const int chare_index[3]"
    body = r"""

  // Unpack griddata and griddata_chare
const params_struct *restrict params = &griddata[grid].params;
const params_struct *restrict params_chare = &griddata_chare[grid].params;
const charecomm_struct *restrict charecommstruct = &griddata_chare[grid].charecommstruct;
diagnostic_struct *restrict diagnosticstruct = &griddata_chare[grid].diagnosticstruct;
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
"""
    if enable_psi4_diagnostics:
        body += r"""
if (which_output == OUTPUT_PSI4) {
   // Do psi4 output, but only if the grid is spherical-like.
  if (strstr(params_chare->CoordSystemName, "Spherical") != NULL) {
    if (diagnosticstruct->num_of_R_exts_chare > 0) {
    // Set psi4.
    psi4_part0(commondata, params_chare, xx_chare, y_n_gfs, diagnostic_output_gfs);
    psi4_part1(commondata, params_chare, xx_chare, y_n_gfs, diagnostic_output_gfs);
    psi4_part2(commondata, params_chare, xx_chare, y_n_gfs, diagnostic_output_gfs);
    // Decompose psi4 into spin-weight -2  spherical harmonics & output to files.
    psi4_spinweightm2_decomposition_on_sphlike_grids(
            commondata, params, params_chare, diagnostic_output_gfs, diagnosticstruct->list_of_R_exts_chare, diagnosticstruct->num_of_R_exts_chare,
            diagnosticstruct->psi4_spinweightm2_sph_harmonics_max_l, xx, xx_chare, chare_index, diagnosticstruct->localsums_for_psi4_decomp, diagnosticstruct->length_localsums_for_psi4_decomp);
    }
  }
} else {
"""
    body += r"""

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
"""
    if enable_BSSN_diagnostics:
        body += r"""

    // Constraint output
    {
      Ricci_eval(commondata, params_chare, &griddata_chare[grid].rfmstruct, y_n_gfs, auxevol_gfs);
      constraints_eval(commondata, params_chare, &griddata_chare[grid].rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_output_gfs);
    }
"""

    body += r"""

    // // 0D, 1D and 2D outputs
    if (which_output == OUTPUT_0D) {
        diagnostics_nearest_grid_center(commondata, params_chare, &griddata_chare[grid].gridfuncs);
    } else if (which_output == OUTPUT_1D_Y) {
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
"""
    if enable_psi4_diagnostics:
        body += r"""
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


def register_CFunction_psi4_spinweightm2_decomposition_on_sphlike_grids() -> None:
    """Register C function for decomposing psi4 into spin-weighted spherical harmonics."""

    register_CFunction_psi4_spinweightm2_decomposition_file_write()

    prefunc = r"""
static void lowlevel_decompose_psi4_into_swm2_modes(const int Nxx_plus_2NGHOSTS1, const int Nxx_plus_2NGHOSTS2, const REAL dxx1, const REAL dxx2,
                                                    const int swm2sh_maximum_l_mode_to_compute, const REAL curr_time, const REAL R_ext,
                                                    const REAL *restrict th_array, const REAL *restrict sinth_array, const REAL *restrict ph_array,
                                                    const REAL *restrict psi4r_at_R_ext, const REAL *restrict psi4i_at_R_ext,
                                                    REAL *restrict localsums_for_psi4_decomp, const int which_R_ext) {

  int which_l = 0;
  int which_m = 0;
  for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
    which_m = 0;
    for (int m = -l; m <= l; m++) {
      // Parallelize the integration loop:
      REAL psi4r_l_m = 0.0;
      REAL psi4i_l_m = 0.0;
#pragma omp parallel for reduction(+ : psi4r_l_m, psi4i_l_m)
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1 - 2 * NGHOSTS; i1++) {
        const REAL th = th_array[i1];
        const REAL sinth = sinth_array[i1];
        for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS; i2++) {
          const REAL ph = ph_array[i2];
          // Construct integrand for psi4 spin-weight s=-2 spherical harmonic
          REAL ReY_sm2_l_m, ImY_sm2_l_m;
          spin_weight_minus2_sph_harmonics(l, m, th, ph, &ReY_sm2_l_m, &ImY_sm2_l_m);

          const int idx2d = i1 * (Nxx_plus_2NGHOSTS2 - 2 * NGHOSTS) + i2;
          const REAL a = psi4r_at_R_ext[idx2d];
          const REAL b = psi4i_at_R_ext[idx2d];
          const REAL c = ReY_sm2_l_m;
          const REAL d = ImY_sm2_l_m;
          psi4r_l_m += (a * c + b * d) * dxx2 * sinth * dxx1;
          psi4i_l_m += (b * c - a * d) * dxx2 * sinth * dxx1;
        }
      }
      localsums_for_psi4_decomp[IDX4PSI4(which_R_ext, which_l, which_m, 0, swm2sh_maximum_l_mode_to_compute - 1,
                                               (2 * swm2sh_maximum_l_mode_to_compute) + 1, 2)] = psi4r_l_m;
      localsums_for_psi4_decomp[IDX4PSI4(which_R_ext, which_l, which_m, 1, swm2sh_maximum_l_mode_to_compute - 1,
                                               (2 * swm2sh_maximum_l_mode_to_compute) + 1, 2)] = psi4i_l_m;
      which_m++;
    }
    which_l++;
  }
}
"""

    desc = "Decompose psi4 across all l,m modes from l=2 up to and including L_MAX (global variable)"
    name = "psi4_spinweightm2_decomposition_on_sphlike_grids"
    params = r"""const commondata_struct *restrict commondata, const params_struct *restrict params,
                                                      const params_struct *restrict params_chare, REAL *restrict diagnostic_output_gfs,
                                                      const REAL *restrict list_of_R_exts, const int num_of_R_exts,
                                                      const int psi4_spinweightm2_sph_harmonics_max_l, REAL *restrict xx[3], REAL *restrict xx_chare[3], const int chare_index[3],
                                                      REAL *restrict localsums_for_psi4_decomp, const int length_localsums_for_psi4_decomp"""
    body = r"""const int Nxx_plus_2NGHOSTS0chare = params_chare->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1chare = params_chare->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2chare = params_chare->Nxx_plus_2NGHOSTS2;
  const int Nxx0chare = params_chare->Nxx0;
  const int Nxx1chare = params_chare->Nxx1;
  const int Nxx2chare = params_chare->Nxx2;

  // Step 1: Allocate memory for 2D arrays used to store psi4, theta, sin(theta), and phi.
  const int sizeof_2Darray = sizeof(REAL) * (Nxx_plus_2NGHOSTS1chare - 2 * NGHOSTS) * (Nxx_plus_2NGHOSTS2chare - 2 * NGHOSTS);
  REAL *restrict psi4r_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  REAL *restrict psi4i_at_R_ext = (REAL *restrict)malloc(sizeof_2Darray);
  //         ... also store theta, sin(theta), and phi to corresponding 1D arrays.
  REAL *restrict sinth_array = (REAL *restrict)malloc(sizeof(REAL) * (Nxx_plus_2NGHOSTS1chare - 2 * NGHOSTS));
  REAL *restrict th_array = (REAL *restrict)malloc(sizeof(REAL) * (Nxx_plus_2NGHOSTS1chare - 2 * NGHOSTS));
  REAL *restrict ph_array = (REAL *restrict)malloc(sizeof(REAL) * (Nxx_plus_2NGHOSTS2chare - 2 * NGHOSTS));

  const int NinterpGHOSTS = MIN(2, NGHOSTS - 1);
  const int N0 = 2 * NinterpGHOSTS; // Interp stencil is 2*NinterpGHOSTS+1 in size;
  //                                 reaches NinterpGHOSTS to the left & right of
  //                                 central point.
  const REAL pow_dxx0__N0 = pow(params->dxx0, N0);

  // Reset to zero
  for (int i = 0; i < length_localsums_for_psi4_decomp; ++i) {
    localsums_for_psi4_decomp[i] = 0.0;
  }

  // Step 2: Loop over all extraction indices:
  for (int which_R_ext = 0; which_R_ext < num_of_R_exts; which_R_ext++) {
    // Step 2.a: Set the extraction radius R_ext based on the radial index R_ext_idx
    const REAL R_ext = list_of_R_exts[which_R_ext];
    const REAL xCart_R_ext[3] = {R_ext, 0.0, 0.0}; // just put a point on the x-axis.

    int Cart_to_i0i1i2[3];
    REAL closest_xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart_R_ext, closest_xx, Cart_to_i0i1i2);

    const int closest_i0 = Cart_to_i0i1i2[0];

    // We want a src grid point inside the source grid (duh) with
    //  mask=+0, as all mask=+0 points will have at least
    //  NGHOSTS>=NinterpGHOSTS of filled neighbor pts.
    // Is grid interior of global grid, not chare's local grid ...
    if (IS_IN_GRID_INTERIOR(Cart_to_i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {

      // Step 2.a.i: Set radial interpolation coefficients for r=R_ext.
      REAL l0i__times__w0i_inv[2 * NinterpGHOSTS + 1];
      {
        for (int i = 0; i <= N0; i++) {
          REAL prod_numer_i = 1.0;
          int prod_denom_i = 1;
          for (int l = 0; l < i; l++) {
            prod_denom_i *= i - l;
            prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0 / 2 + l];
          }
          for (int l = i + 1; l <= N0; l++) {
            prod_denom_i *= i - l;
            prod_numer_i *= closest_xx[0] - xx[0][closest_i0 - N0 / 2 + l];
          }
          l0i__times__w0i_inv[i] = prod_numer_i / (pow_dxx0__N0 * (REAL)prod_denom_i);
        }
      }

      // Step 2.b: Compute psi_4 at this extraction radius and store to a local 2D array.
#pragma omp parallel for
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1chare - NGHOSTS; i1++) {
        th_array[i1 - NGHOSTS] = xx_chare[1][i1];
        sinth_array[i1 - NGHOSTS] = sin(xx_chare[1][i1]);
        for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2chare - NGHOSTS; i2++) {
          ph_array[i2 - NGHOSTS] = xx_chare[2][i2];

          REAL sum_psi4r = 0;
          REAL sum_psi4i = 0;
          // Perform radial interpolation to get psi4 at desired extraction radius R_ext.
          for (int i = 0; i <= N0; i++) {
            // psi4r and psi4i in fixed frame have been stored to diagnostic_output_gfs.
            //  Here we interpolate to specific radius.
            // Convert to local chare index for accessing diagnostic_output_gfs
            const int locali0 = MAP_GLOBAL_TO_LOCAL_IDX0(chare_index[0], i + closest_i0 - N0 / 2, Nxx0chare);
            const int locali1 = i1;
            const int locali2 = i2;

            sum_psi4r += (diagnostic_output_gfs[IDX4GENERAL(PSI4_PART0REGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)] +
                          diagnostic_output_gfs[IDX4GENERAL(PSI4_PART1REGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)] +
                          diagnostic_output_gfs[IDX4GENERAL(PSI4_PART2REGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)]) *
                         l0i__times__w0i_inv[i];
            sum_psi4i += (diagnostic_output_gfs[IDX4GENERAL(PSI4_PART0IMGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)] +
                          diagnostic_output_gfs[IDX4GENERAL(PSI4_PART1IMGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)] +
                          diagnostic_output_gfs[IDX4GENERAL(PSI4_PART2IMGF, locali0, locali1, locali2, Nxx_plus_2NGHOSTS0chare,
                                                            Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare)]) *
                         l0i__times__w0i_inv[i];

          }
          // Store result to "2D" array (actually 1D array with 2D storage):
          const int idx2d = (i1 - NGHOSTS) * (Nxx_plus_2NGHOSTS2chare - 2 * NGHOSTS) + (i2 - NGHOSTS);
          psi4r_at_R_ext[idx2d] = sum_psi4r;
          psi4i_at_R_ext[idx2d] = sum_psi4i;
        }
      }

      // Step 3: Perform integrations across all l,m modes from l=2 up to and including L_MAX (global variable):
      lowlevel_decompose_psi4_into_swm2_modes(Nxx_plus_2NGHOSTS1chare, Nxx_plus_2NGHOSTS2chare, dxx1, dxx2, swm2sh_maximum_l_mode_to_compute, time, R_ext,
                                              th_array, sinth_array, ph_array, psi4r_at_R_ext, psi4i_at_R_ext, localsums_for_psi4_decomp,
                                              which_R_ext);
    }
  }
  // Step 4: Free all allocated memory:
  free(psi4r_at_R_ext);
  free(psi4i_at_R_ext);
  free(sinth_array);
  free(th_array);
  free(ph_array);
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def register_CFunction_psi4_spinweightm2_decomposition_file_write() -> None:
    """Register C function for file write of psi4 decomposition into spin-weighted spherical harmonics."""
    desc = "File write of psi4 decomposition into spin-weighted spherical harmonics"
    name = "psi4_spinweightm2_decomposition_file_write"
    params = r"""const commondata_struct *restrict commondata, diagnostic_struct *restrict diagnosticstruct"""
    body = r"""
const REAL curr_time = commondata->time;
  // Unpack diagnosticptoffset struct:
  const int swm2sh_maximum_l_mode_to_compute = diagnosticstruct->psi4_spinweightm2_sph_harmonics_max_l;
  const int num_of_R_exts = diagnosticstruct->num_of_R_exts_chare;
  const REAL *restrict list_of_R_exts = diagnosticstruct->list_of_R_exts_chare;
  const REAL *restrict globalsums_for_psi4_decomp = diagnosticstruct->globalsums_for_psi4_decomp;

  int which_l = 0;
  int which_m = 0;
  for (int which_R_ext = 0; which_R_ext < num_of_R_exts; which_R_ext++) {
    const REAL R_ext = list_of_R_exts[which_R_ext];
    char filename[100];
    FILE *outpsi4_l_m;
    // Output header at t=0:
    if (curr_time == 0) {
      for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
        sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
        outpsi4_l_m = fopen(filename, "w");
        fprintf(outpsi4_l_m, "# column 1: t-R_ext = [retarded time]\n");
        int col = 2;
        for (int m = -l; m <= l; m++) {
          fprintf(outpsi4_l_m, "# column %d: Re(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
          col++;
          fprintf(outpsi4_l_m, "# column %d: Im(psi4_{l=%d,m=%d}) * R_ext\n", col, l, m);
          col++;
        }
        fclose(outpsi4_l_m);
      }
    }
    // Output one file per l mode; each column represents a unique complex component of l,m
    which_l = 0;
    for (int l = 2; l <= swm2sh_maximum_l_mode_to_compute; l++) {
      sprintf(filename, "Rpsi4_l%d-r%06.1f.txt", l, (double)R_ext);
      outpsi4_l_m = fopen(filename, "a");
      char oneline[10000];
      sprintf(oneline, "%e", (double)(curr_time - R_ext));
      which_m = 0;
      for (int m = -l; m <= l; m++) {
        REAL psi4r_l_m = globalsums_for_psi4_decomp[IDX4PSI4(which_R_ext, which_l, which_m, 0, swm2sh_maximum_l_mode_to_compute - 1, (2 * swm2sh_maximum_l_mode_to_compute) + 1, 2)];
        REAL psi4i_l_m = globalsums_for_psi4_decomp[IDX4PSI4(which_R_ext, which_l, which_m, 1, swm2sh_maximum_l_mode_to_compute - 1, (2 * swm2sh_maximum_l_mode_to_compute) + 1, 2)];
        sprintf(oneline + strlen(oneline), " %.15e %.15e", (double)(R_ext * psi4r_l_m), (double)(R_ext * psi4i_l_m));
        which_m++;
      }
      fprintf(outpsi4_l_m, "%s\n", oneline);
      fclose(outpsi4_l_m);
      which_l++;
    }
  }
"""

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        desc=desc,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )

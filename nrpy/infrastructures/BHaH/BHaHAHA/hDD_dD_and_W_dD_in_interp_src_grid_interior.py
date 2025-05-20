"""
Register C function for computing h_{ij,k} and W_dD within the source grid interior.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp


def register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior() -> None:
    """
    Register the C function 'hDD_dD_and_W_dD_in_src_grid_interior' to compute h_{ij,k}.

    This function sets up the necessary parameters and code body to compute the derivatives
    hDD_dD and W_dD in the source grid interior using OpenMP for parallelization.

    >>> register_CFunction_hDD_dD_and_W_dD_in_interp_src_grid_interior()
    """
    includes = ["BHaH_defines.h"]
    desc = "Compute h_{ij,k}."
    cfunc_type = "void"
    name = "hDD_dD_and_W_dD_in_interp_src_grid_interior"
    params = "commondata_struct *restrict commondata"
    body = r"""
  int i0_min_shift = 0;
  if (commondata->bhahaha_params_and_data->r_min_external_input == 0)
    i0_min_shift = NGHOSTS;

  const int Nxx_plus_2NGHOSTS0 = commondata->interp_src_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->interp_src_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->interp_src_Nxx_plus_2NGHOSTS2;

  const REAL invdxx0 = commondata->interp_src_invdxx0;
  const REAL invdxx1 = commondata->interp_src_invdxx1;
  const REAL invdxx2 = commondata->interp_src_invdxx2;

  // PART 1 OF 2: Compute angular derivatives h_{ij,k} and W_{,k} (k = 1, 2) at all active radial points.
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
      for (int i0 = i0_min_shift; i0 < Nxx_plus_2NGHOSTS0; i0++) {
"""
    # Calling hDD_dD this ensures that c_codegen sets it as a derivative.
    hDD_dD = ixp.declarerank3("hDD_dD", symmetry="sym01")
    hDD_dD_and_W_dD_ang_derivs_expr_list = []
    hDD_dD_and_W_dD_ang_derivs_name_list = []
    for k in range(1, 3):  # only angular derivatives
        for i in range(3):
            for j in range(i, 3):
                hDD_dD_and_W_dD_ang_derivs_expr_list += [hDD_dD[i][j][k]]
                hDD_dD_and_W_dD_ang_derivs_name_list += [
                    f"commondata->interp_src_gfs[IDX4(SRC_PARTIAL_D_HDD{k}{i}{j}GF, i0, i1, i2)]"
                ]
    W_dD = ixp.declarerank1("WW_dD")
    for k in range(1, 3):  # only angular derivatives
        hDD_dD_and_W_dD_ang_derivs_expr_list += [W_dD[k]]
        hDD_dD_and_W_dD_ang_derivs_name_list += [
            f"commondata->interp_src_gfs[IDX4(SRC_PARTIAL_D_WW{k}GF, i0, i1, i2)]"
        ]
    body += ccg.c_codegen(
        hDD_dD_and_W_dD_ang_derivs_expr_list,
        hDD_dD_and_W_dD_ang_derivs_name_list,
        enable_fd_codegen=True,
    ).replace("auxevol_gfs[IDX4(", "commondata->interp_src_gfs[IDX4(SRC_")
    body += r"""
      } // END LOOP over non-inner-boundary points

  // PART 2 OF 2: Compute radial derivatives h_{ij,k} and W_{,k} (k = 0) at ALL interior points.
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
      for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
"""
    hDD_dD_and_W_dD_rad_derivs_expr_list = []
    hDD_dD_and_W_dD_rad_derivs_name_list = []
    for k in range(0, 1):  # only radial derivatives.
        for i in range(3):
            for j in range(i, 3):
                hDD_dD_and_W_dD_rad_derivs_expr_list += [hDD_dD[i][j][k]]
                hDD_dD_and_W_dD_rad_derivs_name_list += [
                    f"commondata->interp_src_gfs[IDX4(SRC_PARTIAL_D_HDD{k}{i}{j}GF, i0, i1, i2)]"
                ]
    W_dD = ixp.declarerank1("WW_dD")
    for k in range(0, 1):  # only radial derivatives.
        hDD_dD_and_W_dD_rad_derivs_expr_list += [W_dD[k]]
        hDD_dD_and_W_dD_rad_derivs_name_list += [
            f"commondata->interp_src_gfs[IDX4(SRC_PARTIAL_D_WW{k}GF, i0, i1, i2)]"
        ]
    body += ccg.c_codegen(
        hDD_dD_and_W_dD_rad_derivs_expr_list,
        hDD_dD_and_W_dD_rad_derivs_name_list,
        enable_fd_codegen=True,
    ).replace("auxevol_gfs[IDX4(", "commondata->interp_src_gfs[IDX4(SRC_")

    body += r"""
} // END LOOP over ALL interior points
"""

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,  # params not passed to function
        body=body,
    )

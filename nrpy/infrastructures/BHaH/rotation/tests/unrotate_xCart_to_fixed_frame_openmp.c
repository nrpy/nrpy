#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Map a Cartesian vector from rotating-frame coordinates to fixed-frame coordinates.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - This routine applies x_fixed = R * x_rot directly.
 *
 * @param[in] commondata Commondata structure with cumulative regrid basis vectors.
 * @param[in,out] xCart Cartesian point/vector to be unrotated in place.
 */
void unrotate_xCart_to_fixed_frame(const commondata_struct commondata, REAL xCart[3]) {
  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  if (so3_validate_and_optionally_fix_hats(xhatU, yhatU, zhatU, 0) != 0) {
    fprintf(stderr,
            "ERROR in %s: cumulative hats are invalid at the consuming boundary. "
            "Repair must be applied at the producing/update boundary.\n",
            __func__);
    exit(1);
  }

  REAL R[3][3];
  build_R_from_cumulative_hats(xhatU, yhatU, zhatU, R);
  so3_apply_R_to_vector(R, xCart, xCart);
} // END FUNCTION unrotate_xCart_to_fixed_frame

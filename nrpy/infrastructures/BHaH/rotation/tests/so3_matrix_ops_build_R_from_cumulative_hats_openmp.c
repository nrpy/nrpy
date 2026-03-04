#include "BHaH_defines.h"

/**
 * @brief Build SO(3) rotation matrix from cumulative rotating-frame hats.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - Columns are hats in fixed basis: R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat.
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout is row-major indexing: R[i][j] is row i, column j.
 *
 * @param[in] xhatU Rotating-frame xhat in fixed coordinates.
 * @param[in] yhatU Rotating-frame yhat in fixed coordinates.
 * @param[in] zhatU Rotating-frame zhat in fixed coordinates.
 * @param[out] R Rotation matrix mapping rotating -> fixed components.
 */
void build_R_from_cumulative_hats(const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]) {
  // R[i][j] means row i, column j. Fill R by columns from hats.
  R[0][0] = xhatU[0];
  R[1][0] = xhatU[1];
  R[2][0] = xhatU[2];

  R[0][1] = yhatU[0];
  R[1][1] = yhatU[1];
  R[2][1] = yhatU[2];

  R[0][2] = zhatU[0];
  R[1][2] = zhatU[1];
  R[2][2] = zhatU[2];
} // END FUNCTION build_R_from_cumulative_hats

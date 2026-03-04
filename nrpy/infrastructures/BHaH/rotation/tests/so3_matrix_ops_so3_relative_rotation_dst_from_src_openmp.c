#include "BHaH_defines.h"

/**
 * @brief Compute relative basis transform DeltaR_dst_from_src = R_dst^T R_src.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 *
 * Given two absolute rotating->fixed rotations, @p R_src and @p R_dst, this
 * returns the matrix that maps components in the src rotating basis to
 * components in the dst rotating basis:
 * v_dst = DeltaR_dst_from_src * v_src.
 *
 * @param[in] R_dst Absolute destination rotation (rotating->fixed).
 * @param[in] R_src Absolute source rotation (rotating->fixed).
 * @param[out] DeltaR_dst_from_src Relative transform from src rotating basis to dst rotating basis.
 */
void so3_relative_rotation_dst_from_src(const REAL R_dst[3][3], const REAL R_src[3][3], REAL DeltaR_dst_from_src[3][3]) {
  // DeltaR_dst_from_src[i][j] = sum_k (R_dst^T)[i][k] R_src[k][j] = sum_k R_dst[k][i] R_src[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += R_dst[k][i] * R_src[k][j];
      }
      DeltaR_dst_from_src[i][j] = accum;
    }
  }
} // END FUNCTION so3_relative_rotation_dst_from_src

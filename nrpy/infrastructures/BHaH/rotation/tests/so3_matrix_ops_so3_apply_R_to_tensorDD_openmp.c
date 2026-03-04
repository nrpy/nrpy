#include "BHaH_defines.h"

/**
 * @brief Apply R to a rank-2 tensor: T_fixed = R T_rot R^T.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout uses R[i][j] = row i, column j.
 *
 * @param[in] R Rotation matrix.
 * @param[in] T_rot Tensor components in rotating basis.
 * @param[out] T_fixed Tensor components in fixed basis.
 */
void so3_apply_R_to_tensorDD(const REAL R[3][3], const REAL T_rot[3][3], REAL T_fixed[3][3]) {
  REAL tmp[3][3];

  // tmp = R * T_rot: tmp[i][j] = sum_k R[i][k] * T_rot[k][j].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += R[i][k] * T_rot[k][j];
      }
      tmp[i][j] = accum;
    }
  }

  // T_fixed = tmp * R^T: T_fixed[i][j] = sum_k tmp[i][k] * R[j][k].
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      REAL accum = 0.0;
      for (int k = 0; k < 3; k++) {
        accum += tmp[i][k] * R[j][k];
      }
      T_fixed[i][j] = accum;
    }
  }
} // END FUNCTION so3_apply_R_to_tensorDD

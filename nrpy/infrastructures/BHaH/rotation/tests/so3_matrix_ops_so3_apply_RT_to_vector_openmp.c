#include "BHaH_defines.h"

/**
 * @brief Apply R^T to a vector using v_rot = R^T v_fixed.
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
 * @param[in] v_fixed Vector components in fixed basis.
 * @param[out] v_rot Vector components in rotating basis.
 */
void so3_apply_RT_to_vector(const REAL R[3][3], const REAL v_fixed[3], REAL v_rot[3]) {
  // Alias-safe copy: permits v_rot == v_fixed.
  const REAL v_in[3] = {v_fixed[0], v_fixed[1], v_fixed[2]};
  // v_rot[i] = sum_j (R^T)[i][j] v_in[j] = sum_j R[j][i] v_in[j].
  for (int i = 0; i < 3; i++) {
    REAL accum = 0.0;
    for (int j = 0; j < 3; j++) {
      accum += R[j][i] * v_in[j];
    }
    v_rot[i] = accum;
  }
} // END FUNCTION so3_apply_RT_to_vector

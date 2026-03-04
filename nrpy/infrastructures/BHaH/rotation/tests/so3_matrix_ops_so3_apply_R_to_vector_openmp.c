#include "BHaH_defines.h"

/**
 * @brief Apply R to a vector using v_fixed = R v_rot.
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
 * @param[in] v_rot Vector components in rotating basis.
 * @param[out] v_fixed Vector components in fixed basis.
 */
void so3_apply_R_to_vector(const REAL R[3][3], const REAL v_rot[3], REAL v_fixed[3]) {
  // Alias-safe copy: permits v_fixed == v_rot.
  const REAL v_in[3] = {v_rot[0], v_rot[1], v_rot[2]};
  // v_fixed[i] = sum_j R[i][j] v_in[j] with R[i][j] = row i, column j.
  for (int i = 0; i < 3; i++) {
    REAL accum = 0.0;
    for (int j = 0; j < 3; j++) {
      accum += R[i][j] * v_in[j];
    }
    v_fixed[i] = accum;
  }
} // END FUNCTION so3_apply_R_to_vector

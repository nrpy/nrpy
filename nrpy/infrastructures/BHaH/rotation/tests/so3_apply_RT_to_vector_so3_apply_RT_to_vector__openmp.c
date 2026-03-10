#include "BHaH_defines.h"

/**
 * @brief Apply the transpose of a rotation matrix to a vector in place.
 *
 * Convention:
 * - ``R`` maps rotating-frame components to fixed-frame components.
 * - ``R^T`` maps fixed-frame components to rotating-frame components.
 * - Vector update in Einstein notation: ``v^i_out = (R^T)^i{}_j v^j_in``.
 * - C layout statement: ``R[i][j]`` is row ``i``, column ``j``.
 *
 * This routine is alias-safe: it snapshots the input vector before overwriting
 * ``vU``.
 *
 * @param[in] R Rotation matrix whose transpose is applied.
 * @param[in,out] vU Vector updated in place with ``R^T * vU``.
 */
void so3_apply_RT_to_vector(const REAL R[3][3], REAL vU[3]) {
  // Snapshot the input vector so the in-place update is alias-safe.
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
  vU[0] = R[0][0] * v_in[0] + R[1][0] * v_in[1] + R[2][0] * v_in[2];
  vU[1] = R[0][1] * v_in[0] + R[1][1] * v_in[1] + R[2][1] * v_in[2];
  vU[2] = R[0][2] * v_in[0] + R[1][2] * v_in[1] + R[2][2] * v_in[2];
} // END FUNCTION so3_apply_RT_to_vector

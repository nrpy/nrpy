#include "BHaH_defines.h"

/**
 * Apply a rotation matrix to a rank-2 covariant tensor in place.
 *
 * Convention:
 * - ``R`` maps rotating-frame components to fixed-frame components.
 * - Tensor update in Einstein notation: ``T_out_ij = R_i^k R_j^l T_in_kl``.
 * - C layout statement: ``R[i][j]`` is row ``i``, column ``j``.
 *
 * This routine is alias-safe: it snapshots the input tensor before overwriting
 * ``tDD``. It writes all 9 components; callers may repack symmetric storage after
 * the call.
 *
 * @param[in] R Rotation matrix.
 * @param[in,out] tDD Tensor updated in place with ``R * tDD * R^T``.
 */
void so3_apply_R_to_tensorDD(const REAL R[3][3], REAL tDD[3][3]) {
  // Snapshot the input tensor so the in-place update is alias-safe.
  const REAL t_in[3][3] = {
      {tDD[0][0], tDD[0][1], tDD[0][2]},
      {tDD[1][0], tDD[1][1], tDD[1][2]},
      {tDD[2][0], tDD[2][1], tDD[2][2]},
  };
  const REAL tmp14 = R[1][1] * t_in[1][0];
  const REAL tmp15 = R[1][2] * t_in[2][0];
  const REAL tmp16 = R[1][0] * t_in[0][1];
  const REAL tmp17 = R[1][2] * t_in[2][1];
  const REAL tmp18 = R[1][0] * t_in[0][2];
  const REAL tmp19 = R[1][1] * t_in[1][2];
  const REAL tmp20 = R[1][0] * R[2][0] * t_in[0][0] + R[1][1] * R[2][1] * t_in[1][1] + R[1][2] * R[2][2] * t_in[2][2];
  const REAL tmp21 = R[2][1] * t_in[1][0];
  const REAL tmp22 = R[2][2] * t_in[2][0];
  const REAL tmp23 = R[2][0] * t_in[0][1];
  const REAL tmp24 = R[2][2] * t_in[2][1];
  const REAL tmp25 = R[2][0] * t_in[0][2];
  const REAL tmp26 = R[2][1] * t_in[1][2];
  const REAL tmp12 = R[0][0] * R[1][0] * t_in[0][0] + R[0][1] * R[1][1] * t_in[1][1] + R[0][2] * R[1][2] * t_in[2][2];
  const REAL tmp13 = R[0][0] * R[2][0] * t_in[0][0] + R[0][1] * R[2][1] * t_in[1][1] + R[0][2] * R[2][2] * t_in[2][2];
  tDD[0][0] = ((R[0][0]) * (R[0][0])) * t_in[0][0] + R[0][0] * R[0][1] * t_in[0][1] + R[0][0] * R[0][1] * t_in[1][0] +
              R[0][0] * R[0][2] * t_in[0][2] + R[0][0] * R[0][2] * t_in[2][0] + ((R[0][1]) * (R[0][1])) * t_in[1][1] +
              R[0][1] * R[0][2] * t_in[1][2] + R[0][1] * R[0][2] * t_in[2][1] + ((R[0][2]) * (R[0][2])) * t_in[2][2];
  tDD[0][1] = R[0][0] * R[1][1] * t_in[0][1] + R[0][0] * R[1][2] * t_in[0][2] + R[0][1] * R[1][0] * t_in[1][0] + R[0][1] * R[1][2] * t_in[1][2] +
              R[0][2] * R[1][0] * t_in[2][0] + R[0][2] * R[1][1] * t_in[2][1] + tmp12;
  tDD[0][2] = R[0][0] * R[2][1] * t_in[0][1] + R[0][0] * R[2][2] * t_in[0][2] + R[0][1] * R[2][0] * t_in[1][0] + R[0][1] * R[2][2] * t_in[1][2] +
              R[0][2] * R[2][0] * t_in[2][0] + R[0][2] * R[2][1] * t_in[2][1] + tmp13;
  tDD[1][0] = R[0][0] * tmp14 + R[0][0] * tmp15 + R[0][1] * tmp16 + R[0][1] * tmp17 + R[0][2] * tmp18 + R[0][2] * tmp19 + tmp12;
  tDD[1][1] = ((R[1][0]) * (R[1][0])) * t_in[0][0] + R[1][0] * tmp14 + R[1][0] * tmp15 + ((R[1][1]) * (R[1][1])) * t_in[1][1] + R[1][1] * tmp16 +
              R[1][1] * tmp17 + ((R[1][2]) * (R[1][2])) * t_in[2][2] + R[1][2] * tmp18 + R[1][2] * tmp19;
  tDD[1][2] = R[2][0] * tmp14 + R[2][0] * tmp15 + R[2][1] * tmp16 + R[2][1] * tmp17 + R[2][2] * tmp18 + R[2][2] * tmp19 + tmp20;
  tDD[2][0] = R[0][0] * tmp21 + R[0][0] * tmp22 + R[0][1] * tmp23 + R[0][1] * tmp24 + R[0][2] * tmp25 + R[0][2] * tmp26 + tmp13;
  tDD[2][1] = R[1][0] * tmp21 + R[1][0] * tmp22 + R[1][1] * tmp23 + R[1][1] * tmp24 + R[1][2] * tmp25 + R[1][2] * tmp26 + tmp20;
  tDD[2][2] = ((R[2][0]) * (R[2][0])) * t_in[0][0] + R[2][0] * tmp21 + R[2][0] * tmp22 + ((R[2][1]) * (R[2][1])) * t_in[1][1] + R[2][1] * tmp23 +
              R[2][1] * tmp24 + ((R[2][2]) * (R[2][2])) * t_in[2][2] + R[2][2] * tmp25 + R[2][2] * tmp26;
} // END FUNCTION: so3_apply_R_to_tensorDD

#include "BHaH_defines.h"

/**
 * @brief Update cumulative hat vectors by left-multiplying with a rotation matrix.
 *
 * Convention:
 * - ``R_old`` is the cumulative rotating-to-fixed matrix whose columns are the
 *   input hat vectors.
 * - ``R_new = DeltaR * R_old``.
 * - Each updated hat vector is one column of ``R_new``.
 * - C layout statement: ``DeltaR[i][j]`` is row ``i``, column ``j``.
 *
 * This routine is alias-safe: it snapshots all input hat vectors before
 * overwriting them in place.
 *
 * @param[in,out] xhatU Cumulative rotating-frame x-basis vector, updated in place.
 * @param[in,out] yhatU Cumulative rotating-frame y-basis vector, updated in place.
 * @param[in,out] zhatU Cumulative rotating-frame z-basis vector, updated in place.
 * @param[in] DeltaR Rotation matrix that left-multiplies the cumulative basis.
 */
void so3_left_multiply_hats_with_R(REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const REAL DeltaR[3][3]) {
  // Snapshot the input hat vectors so the in-place update is alias-safe.
  const REAL xhat_in[3] = {xhatU[0], xhatU[1], xhatU[2]};
  const REAL yhat_in[3] = {yhatU[0], yhatU[1], yhatU[2]};
  const REAL zhat_in[3] = {zhatU[0], zhatU[1], zhatU[2]};
  xhatU[0] = DeltaR[0][0] * xhat_in[0] + DeltaR[0][1] * xhat_in[1] + DeltaR[0][2] * xhat_in[2];
  xhatU[1] = DeltaR[1][0] * xhat_in[0] + DeltaR[1][1] * xhat_in[1] + DeltaR[1][2] * xhat_in[2];
  xhatU[2] = DeltaR[2][0] * xhat_in[0] + DeltaR[2][1] * xhat_in[1] + DeltaR[2][2] * xhat_in[2];
  yhatU[0] = DeltaR[0][0] * yhat_in[0] + DeltaR[0][1] * yhat_in[1] + DeltaR[0][2] * yhat_in[2];
  yhatU[1] = DeltaR[1][0] * yhat_in[0] + DeltaR[1][1] * yhat_in[1] + DeltaR[1][2] * yhat_in[2];
  yhatU[2] = DeltaR[2][0] * yhat_in[0] + DeltaR[2][1] * yhat_in[1] + DeltaR[2][2] * yhat_in[2];
  zhatU[0] = DeltaR[0][0] * zhat_in[0] + DeltaR[0][1] * zhat_in[1] + DeltaR[0][2] * zhat_in[2];
  zhatU[1] = DeltaR[1][0] * zhat_in[0] + DeltaR[1][1] * zhat_in[1] + DeltaR[1][2] * zhat_in[2];
  zhatU[2] = DeltaR[2][0] * zhat_in[0] + DeltaR[2][1] * zhat_in[1] + DeltaR[2][2] * zhat_in[2];
} // END FUNCTION so3_left_multiply_hats_with_R

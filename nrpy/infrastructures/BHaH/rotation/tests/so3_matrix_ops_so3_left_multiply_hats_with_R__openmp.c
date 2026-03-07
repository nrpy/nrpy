#include "BHaH_defines.h"

/**
 * @brief Update cumulative hats by left-multiplying with `DeltaR`.
 *
 * This applies `R_new = DeltaR * R_old`, where the hat vectors are the columns of
 * `R_old`.
 */
void so3_left_multiply_hats_with_R(REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const REAL DeltaR[3][3]) {
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

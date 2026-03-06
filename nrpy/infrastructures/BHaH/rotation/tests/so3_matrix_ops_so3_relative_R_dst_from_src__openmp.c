#include "BHaH_defines.h"

/**
 * @brief Compute `DeltaR_dst_from_src = R_dst^T R_src`.
 */
void so3_relative_R_dst_from_src(const REAL R_dst[3][3], const REAL R_src[3][3], REAL DeltaR[3][3]) {
  DeltaR[0][0] = R_dst[0][0] * R_src[0][0] + R_dst[1][0] * R_src[1][0] + R_dst[2][0] * R_src[2][0];
  DeltaR[0][1] = R_dst[0][0] * R_src[0][1] + R_dst[1][0] * R_src[1][1] + R_dst[2][0] * R_src[2][1];
  DeltaR[0][2] = R_dst[0][0] * R_src[0][2] + R_dst[1][0] * R_src[1][2] + R_dst[2][0] * R_src[2][2];
  DeltaR[1][0] = R_dst[0][1] * R_src[0][0] + R_dst[1][1] * R_src[1][0] + R_dst[2][1] * R_src[2][0];
  DeltaR[1][1] = R_dst[0][1] * R_src[0][1] + R_dst[1][1] * R_src[1][1] + R_dst[2][1] * R_src[2][1];
  DeltaR[1][2] = R_dst[0][1] * R_src[0][2] + R_dst[1][1] * R_src[1][2] + R_dst[2][1] * R_src[2][2];
  DeltaR[2][0] = R_dst[0][2] * R_src[0][0] + R_dst[1][2] * R_src[1][0] + R_dst[2][2] * R_src[2][0];
  DeltaR[2][1] = R_dst[0][2] * R_src[0][1] + R_dst[1][2] * R_src[1][1] + R_dst[2][2] * R_src[2][1];
  DeltaR[2][2] = R_dst[0][2] * R_src[0][2] + R_dst[1][2] * R_src[1][2] + R_dst[2][2] * R_src[2][2];
} // END FUNCTION so3_relative_R_dst_from_src

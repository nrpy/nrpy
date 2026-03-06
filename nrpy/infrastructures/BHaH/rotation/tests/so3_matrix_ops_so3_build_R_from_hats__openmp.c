#include "BHaH_defines.h"

/**
 * @brief Build `R` from cumulative hats.
 *
 * Convention:
 * - `R` maps rotating-frame components -> fixed-frame components.
 * - `R[:,0]=xhat`, `R[:,1]=yhat`, `R[:,2]=zhat`.
 */
void so3_build_R_from_hats(const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]) {
  R[0][0] = xhatU[0];
  R[0][1] = yhatU[0];
  R[0][2] = zhatU[0];
  R[1][0] = xhatU[1];
  R[1][1] = yhatU[1];
  R[1][2] = zhatU[1];
  R[2][0] = xhatU[2];
  R[2][1] = yhatU[2];
  R[2][2] = zhatU[2];
} // END FUNCTION so3_build_R_from_hats

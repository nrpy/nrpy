#include "BHaH_defines.h"

/**
 * @brief Apply `R^T` to a vector in place.
 */
void so3_apply_RT_to_vector(const REAL R[3][3], REAL vU[3]) {
  const REAL v_in[3] = {vU[0], vU[1], vU[2]};
  vU[0] = R[0][0] * v_in[0] + R[1][0] * v_in[1] + R[2][0] * v_in[2];
  vU[1] = R[0][1] * v_in[0] + R[1][1] * v_in[1] + R[2][1] * v_in[2];
  vU[2] = R[0][2] * v_in[0] + R[1][2] * v_in[1] + R[2][2] * v_in[2];
} // END FUNCTION so3_apply_RT_to_vector

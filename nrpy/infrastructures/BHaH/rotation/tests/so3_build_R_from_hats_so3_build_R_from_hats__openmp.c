#include "BHaH_defines.h"

/**
 * @brief Build a cumulative rotation matrix from cumulative hat vectors.
 *
 * Convention:
 * - ``R`` maps rotating-frame components to fixed-frame components.
 * - ``R[:,0]=xhat``, ``R[:,1]=yhat``, ``R[:,2]=zhat``.
 * - C layout statement: ``R[i][j]`` is row ``i``, column ``j``.
 *
 * @param[in] xhatU Cumulative rotating-frame x-basis vector in fixed components.
 * @param[in] yhatU Cumulative rotating-frame y-basis vector in fixed components.
 * @param[in] zhatU Cumulative rotating-frame z-basis vector in fixed components.
 * @param[out] R Rotation matrix assembled from the supplied hat vectors.
 */
void so3_build_R_from_hats(const REAL xhatU[3], const REAL yhatU[3], const REAL zhatU[3], REAL R[3][3]) {
  // Assemble the matrix column-by-column from cumulative hat vectors.
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

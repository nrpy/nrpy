#include "BHaH_defines.h"

/**
 * @brief Recover axis-angle data that maps one unit vector to another.
 *
 * Convention:
 * - ``aU`` and ``bU`` are interpreted as unit vectors in a common basis.
 * - ``nU`` and ``dphi`` satisfy ``bU = R(nU,dphi) aU``.
 * - Parallel input returns the deterministic identity pair
 *   ``nU=(1,0,0)``, ``dphi=0``.
 * - Antiparallel input returns ``dphi=pi`` with a deterministic axis
 *   orthogonal to ``aU``.
 *
 * @param[in] aU Source unit vector.
 * @param[in] bU Destination unit vector.
 * @param[out] nU Unit rotation axis.
 * @param[out] dphi Rotation angle in radians.
 */
void so3_find_nU_and_dphi_from_unit_vectors(const REAL aU[3], const REAL bU[3], REAL nU[3], REAL *restrict dphi) {
  REAL crossU[3];
  REAL dot;
  crossU[0] = aU[1] * bU[2] - aU[2] * bU[1];
  crossU[1] = -aU[0] * bU[2] + aU[2] * bU[0];
  crossU[2] = aU[0] * bU[1] - aU[1] * bU[0];
  dot = aU[0] * bU[0] + aU[1] * bU[1] + aU[2] * bU[2];

  const REAL cross_norm = sqrt(crossU[0] * crossU[0] + crossU[1] * crossU[1] + crossU[2] * crossU[2]);

  if (dot > 1.0)
    dot = 1.0;
  if (dot < -1.0)
    dot = -1.0;

  if (cross_norm < 1e-14) {
    if (dot > 0.0) {
      nU[0] = 1.0;
      nU[1] = 0.0;
      nU[2] = 0.0;
      *dphi = 0.0;
      return;
    } // END IF vectors are parallel.

    int basis_dir = 0;
    if (fabs(aU[1]) < fabs(aU[basis_dir]))
      basis_dir = 1;
    if (fabs(aU[2]) < fabs(aU[basis_dir]))
      basis_dir = 2;
    const REAL eU[3] = {basis_dir == 0 ? 1.0 : 0.0, basis_dir == 1 ? 1.0 : 0.0, basis_dir == 2 ? 1.0 : 0.0};
    const REAL axisU[3] = {
        aU[1] * eU[2] - aU[2] * eU[1],
        aU[2] * eU[0] - aU[0] * eU[2],
        aU[0] * eU[1] - aU[1] * eU[0],
    };
    const REAL axis_norm = sqrt(axisU[0] * axisU[0] + axisU[1] * axisU[1] + axisU[2] * axisU[2]);
    if (axis_norm < 1e-14) {
      fprintf(stderr, "ERROR: %s failed to build a deterministic antiparallel rotation axis.\n", __func__);
      exit(1);
    } // END IF deterministic antiparallel axis is degenerate.
    nU[0] = axisU[0] / axis_norm;
    nU[1] = axisU[1] / axis_norm;
    nU[2] = axisU[2] / axis_norm;
    *dphi = M_PI;
    return;
  } // END IF parallel/antiparallel branch.

  nU[0] = crossU[0] / cross_norm;
  nU[1] = crossU[1] / cross_norm;
  nU[2] = crossU[2] / cross_norm;
  *dphi = acos(dot);
} // END FUNCTION so3_find_nU_and_dphi_from_unit_vectors

#include "BHaH_defines.h"

/**
 * @brief Build a rotation matrix from axis-angle input.
 *
 * Convention:
 * - ``R`` maps rotating-frame components to fixed-frame components.
 * - ``v_fixed = R v_rot`` and ``v_rot = R^T v_fixed``.
 * - C layout statement: ``R[i][j]`` is row ``i``, column ``j``.
 *
 * If ``|dphi|`` is effectively zero or ``|nU|`` is negligible, this routine
 * returns the identity matrix as a deterministic no-op rotation.
 *
 * @param[in] nU Rotation axis, which need not be normalized.
 * @param[in] dphi Rotation angle in radians.
 * @param[out] R Rotation matrix computed from the axis-angle pair.
 */
void so3_axis_angle_to_R(const REAL nU[3], const REAL dphi, REAL R[3][3]) {
  const REAL n0 = nU[0];
  const REAL n1 = nU[1];
  const REAL n2 = nU[2];
  const REAL nnorm = sqrt(n0 * n0 + n1 * n1 + n2 * n2);
  const REAL safe_nnorm = (nnorm < 1e-300) ? 1.0 : nnorm;

  // Normalize the input axis before evaluating Rodrigues' formula.
  const REAL nx = n0 / safe_nnorm;
  const REAL ny = n1 / safe_nnorm;
  const REAL nz = n2 / safe_nnorm;
  const REAL tmp0 = cos(dphi);
  const REAL tmp2 = sin(dphi);
  const REAL tmp1 = 1 - tmp0;
  const REAL tmp5 = nx * ny * tmp1;
  const REAL tmp7 = nx * nz * tmp1;
  R[0][0] = ((nx) * (nx)) * tmp1 + tmp0;
  R[0][1] = -nz * tmp2 + tmp5;
  R[0][2] = ny * tmp2 + tmp7;
  R[1][0] = nz * tmp2 + tmp5;
  R[1][1] = ((ny) * (ny)) * tmp1 + tmp0;
  R[1][2] = -nx * tmp2 + ny * nz * tmp1;
  R[2][0] = -ny * tmp2 + tmp7;
  R[2][1] = nx * tmp2 + ny * nz * tmp1;
  R[2][2] = ((nz) * (nz)) * tmp1 + tmp0;

  if (fabs(dphi) < 1e-15 || nnorm < 1e-300) {
    // Deterministic identity fallback for degenerate axis-angle input.
    R[0][0] = 1.0;
    R[0][1] = 0.0;
    R[0][2] = 0.0;
    R[1][0] = 0.0;
    R[1][1] = 1.0;
    R[1][2] = 0.0;
    R[2][0] = 0.0;
    R[2][1] = 0.0;
    R[2][2] = 1.0;
  } // END IF: degenerate axis-angle input
} // END FUNCTION so3_axis_angle_to_R

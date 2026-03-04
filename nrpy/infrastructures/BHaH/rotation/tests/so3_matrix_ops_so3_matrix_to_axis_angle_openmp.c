#include "BHaH_defines.h"

/**
 * @brief Convert SO(3) matrix to axis-angle with robust pi-branch handling.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 *
 * Algorithm:
 * 1) c = clamp((trace(R)-1)/2, -1, 1), phi = acos(c)
 * 2) Small-angle branch: phi < 1e-12 -> n=(1,0,0), phi=0
 * 3) General branch: |pi-phi| >= 1e-8 -> n from skew(R)/(2 sin(phi))
 * 4) Pi branch: |pi-phi| < 1e-8 -> diagonal-based recovery plus deterministic
 *    sign resolution from symmetric off-diagonal terms
 * 5) Debug check: reconstruct R(n,phi) and abort if mismatch exceeds tolerance
 *
 * @param[in] R Input rotation matrix.
 * @param[out] nU Unit axis vector.
 * @param[out] dphi Rotation angle in radians.
 */
void so3_matrix_to_axis_angle(const REAL R[3][3], REAL nU[3], REAL *restrict dphi) {
  const REAL pi = acos(-1.0);
  const REAL traceR = R[0][0] + R[1][1] + R[2][2];
  REAL c = 0.5 * (traceR - 1.0);
  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;

  REAL phi = acos(c);
  *dphi = phi;

  if (phi < 1e-12) {
    nU[0] = 1.0;
    nU[1] = 0.0;
    nU[2] = 0.0;
    *dphi = 0.0;
  } else if (fabs(pi - phi) >= 1e-8) {
    const REAL sinphi = sin(phi);
    const REAL denom = 2.0 * sinphi;
    if (fabs(denom) < 1e-15) {
      fprintf(stderr, "ERROR in %s: unstable general branch (|2 sin(phi)| too small). phi=%.17e denom=%.17e\n", __func__, phi, denom);
      exit(1);
    }

    nU[0] = (R[2][1] - R[1][2]) / denom;
    nU[1] = (R[0][2] - R[2][0]) / denom;
    nU[2] = (R[1][0] - R[0][1]) / denom;
  } else {
    REAL n0 = 0.5 * (R[0][0] + 1.0);
    REAL n1 = 0.5 * (R[1][1] + 1.0);
    REAL n2 = 0.5 * (R[2][2] + 1.0);
    if (n0 < 0.0)
      n0 = 0.0;
    if (n1 < 0.0)
      n1 = 0.0;
    if (n2 < 0.0)
      n2 = 0.0;

    n0 = sqrt(n0);
    n1 = sqrt(n1);
    n2 = sqrt(n2);

    int imax = 0;
    if (n1 > n0 && n1 >= n2)
      imax = 1;
    else if (n2 > n0 && n2 > n1)
      imax = 2;

    if (imax == 0) {
      if (n1 > 0.0)
        n1 = copysign(n1, R[0][1] + R[1][0]);
      if (n2 > 0.0)
        n2 = copysign(n2, R[0][2] + R[2][0]);
    } else if (imax == 1) {
      if (n0 > 0.0)
        n0 = copysign(n0, R[0][1] + R[1][0]);
      if (n2 > 0.0)
        n2 = copysign(n2, R[1][2] + R[2][1]);
    } else {
      if (n0 > 0.0)
        n0 = copysign(n0, R[0][2] + R[2][0]);
      if (n1 > 0.0)
        n1 = copysign(n1, R[1][2] + R[2][1]);
    }

    nU[0] = n0;
    nU[1] = n1;
    nU[2] = n2;

    REAL nnorm = sqrt(nU[0] * nU[0] + nU[1] * nU[1] + nU[2] * nU[2]);
    if (nnorm < 1e-15) {
      // Deterministic fallback in fully degenerate pi-branch.
      if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) {
        nU[0] = 1.0;
        nU[1] = 0.0;
        nU[2] = 0.0;
      } else if (R[1][1] >= R[0][0] && R[1][1] >= R[2][2]) {
        nU[0] = 0.0;
        nU[1] = 1.0;
        nU[2] = 0.0;
      } else {
        nU[0] = 0.0;
        nU[1] = 0.0;
        nU[2] = 1.0;
      }
      nnorm = 1.0;
    }
    nU[0] /= nnorm;
    nU[1] /= nnorm;
    nU[2] /= nnorm;

    // Deterministic sign convention for pi-branch axis ambiguity.
    if (nU[0] < 0.0 || (fabs(nU[0]) <= 1e-16 && nU[1] < 0.0) || (fabs(nU[0]) <= 1e-16 && fabs(nU[1]) <= 1e-16 && nU[2] < 0.0)) {
      nU[0] = -nU[0];
      nU[1] = -nU[1];
      nU[2] = -nU[2];
    }
    *dphi = pi;
  }

  REAL nnorm = sqrt(nU[0] * nU[0] + nU[1] * nU[1] + nU[2] * nU[2]);
  if (nnorm < 1e-15) {
    fprintf(stderr, "ERROR in %s: could not recover nonzero rotation axis.\n", __func__);
    exit(1);
  }
  nU[0] /= nnorm;
  nU[1] /= nnorm;
  nU[2] /= nnorm;

  // Debug reconstruction check: R_check(n,phi) should match input R.
  const REAL cphi = cos(*dphi);
  const REAL sphi = sin(*dphi);
  const REAL one_minus_c = 1.0 - cphi;
  const REAL nx = nU[0];
  const REAL ny = nU[1];
  const REAL nz = nU[2];

  REAL Rcheck[3][3];
  Rcheck[0][0] = cphi + one_minus_c * nx * nx;
  Rcheck[0][1] = one_minus_c * nx * ny - sphi * nz;
  Rcheck[0][2] = one_minus_c * nx * nz + sphi * ny;
  Rcheck[1][0] = one_minus_c * ny * nx + sphi * nz;
  Rcheck[1][1] = cphi + one_minus_c * ny * ny;
  Rcheck[1][2] = one_minus_c * ny * nz - sphi * nx;
  Rcheck[2][0] = one_minus_c * nz * nx - sphi * ny;
  Rcheck[2][1] = one_minus_c * nz * ny + sphi * nx;
  Rcheck[2][2] = cphi + one_minus_c * nz * nz;

  REAL max_abs_err = 0.0;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      const REAL err = fabs(Rcheck[i][j] - R[i][j]);
      if (err > max_abs_err)
        max_abs_err = err;
    }
  }

  if (max_abs_err > 1e-10) {
    fprintf(stderr,
            "ERROR in %s: matrix->axis-angle reconstruction mismatch. max_abs_err=%.17e phi=%.17e "
            "n=(%.17e, %.17e, %.17e)\n",
            __func__, max_abs_err, *dphi, nU[0], nU[1], nU[2]);
    exit(1);
  }
} // END FUNCTION so3_matrix_to_axis_angle

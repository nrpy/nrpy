#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Compute axis-angle outputs from cumulative hats.
 *
 * This routine computes a single unrotation from the cumulative basis matrix and converts it to
 * one axis-angle pair in (@p nU_part1, @p dphi_part1). The second pair is set to
 * identity: (@p nU_part2=(1,0,0), @p dphi_part2=0).
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - Axis-angle pair is recovered directly from this R with robust near-pi handling.
 * - Einstein notation for the matrix map represented by part 1:
 *   v^i_fixed = R^i{}_j v^j_rot.
 *
 * @param[in] commondata Commondata structure containing cumulative regrid basis vectors.
 * @param[out] nU_part1 Rotation axis for first rotation.
 * @param[out] dphi_part1 Rotation angle for first rotation.
 * @param[out] nU_part2 Rotation axis for second rotation.
 * @param[out] dphi_part2 Rotation angle for second rotation.
 */
void unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame(const commondata_struct commondata, REAL nU_part1[3], REAL *restrict dphi_part1,
                                                              REAL nU_part2[3], REAL *restrict dphi_part2) {
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;

  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  REAL x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR_hats;
  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];
  x_norm = sqrt(((xhatU[0]) * (xhatU[0])) + ((xhatU[1]) * (xhatU[1])) + ((xhatU[2]) * (xhatU[2])));
  y_norm = sqrt(((yhatU[0]) * (yhatU[0])) + ((yhatU[1]) * (yhatU[1])) + ((yhatU[2]) * (yhatU[2])));
  z_norm = sqrt(((zhatU[0]) * (zhatU[0])) + ((zhatU[1]) * (zhatU[1])) + ((zhatU[2]) * (zhatU[2])));
  detR_hats = xhatU[0] * (yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1]) + xhatU[1] * (-yhatU[0] * zhatU[2] + yhatU[2] * zhatU[0]) +
              xhatU[2] * (yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0]);

  const int hats_invalid = (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol || fabs(x_norm - 1.0) > norm_tol ||
                            fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR_hats < det_tol);

  if (hats_invalid) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected at consuming boundary. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR_hats);
    exit(1);
  }

  REAL R[3][3];
  R[0][0] = xhatU[0];
  R[0][1] = yhatU[0];
  R[0][2] = zhatU[0];
  R[1][0] = xhatU[1];
  R[1][1] = yhatU[1];
  R[1][2] = zhatU[1];
  R[2][0] = xhatU[2];
  R[2][1] = yhatU[2];
  R[2][2] = zhatU[2];

  // Robust matrix-to-axis-angle conversion.
  const REAL pi = acos(-1.0);
  REAL traceR;
  traceR = R[0][0] + R[1][1] + R[2][2];

  REAL c = 0.5 * (traceR - 1.0);
  if (c > 1.0)
    c = 1.0;
  if (c < -1.0)
    c = -1.0;

  REAL phi = acos(c);
  *dphi_part1 = phi;

  if (phi < 1e-12) {
    nU_part1[0] = 1.0;
    nU_part1[1] = 0.0;
    nU_part1[2] = 0.0;
    *dphi_part1 = 0.0;
  } else if (fabs(pi - phi) >= 1e-8) {
    const REAL sinphi = sin(phi);
    const REAL denom = 2.0 * sinphi;
    if (fabs(denom) < 1e-15) {
      fprintf(stderr, "ERROR in %s: unstable general branch (|2 sin(phi)| too small). phi=%.17e denom=%.17e\n", __func__, phi, denom);
      exit(1);
    }

    const REAL tmp0 = (1.0 / (denom));
    nU_part1[0] = tmp0 * (-R[1][2] + R[2][1]);
    nU_part1[1] = tmp0 * (R[0][2] - R[2][0]);
    nU_part1[2] = tmp0 * (-R[0][1] + R[1][0]);

  } else {
    REAL n0, n1, n2;
    n0 = (1.0 / 2.0) * R[0][0] + 1.0 / 2.0;
    n1 = (1.0 / 2.0) * R[1][1] + 1.0 / 2.0;
    n2 = (1.0 / 2.0) * R[2][2] + 1.0 / 2.0;

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

    nU_part1[0] = n0;
    nU_part1[1] = n1;
    nU_part1[2] = n2;

    REAL nnorm = sqrt(nU_part1[0] * nU_part1[0] + nU_part1[1] * nU_part1[1] + nU_part1[2] * nU_part1[2]);
    if (nnorm < 1e-15) {
      // Deterministic fallback in fully degenerate pi-branch.
      if (R[0][0] >= R[1][1] && R[0][0] >= R[2][2]) {
        nU_part1[0] = 1.0;
        nU_part1[1] = 0.0;
        nU_part1[2] = 0.0;
      } else if (R[1][1] >= R[0][0] && R[1][1] >= R[2][2]) {
        nU_part1[0] = 0.0;
        nU_part1[1] = 1.0;
        nU_part1[2] = 0.0;
      } else {
        nU_part1[0] = 0.0;
        nU_part1[1] = 0.0;
        nU_part1[2] = 1.0;
      }
      nnorm = 1.0;
    }
    nU_part1[0] /= nnorm;
    nU_part1[1] /= nnorm;
    nU_part1[2] /= nnorm;

    // Deterministic sign convention for pi-branch axis ambiguity.
    if (nU_part1[0] < 0.0 || (fabs(nU_part1[0]) <= 1e-16 && nU_part1[1] < 0.0) ||
        (fabs(nU_part1[0]) <= 1e-16 && fabs(nU_part1[1]) <= 1e-16 && nU_part1[2] < 0.0)) {
      nU_part1[0] = -nU_part1[0];
      nU_part1[1] = -nU_part1[1];
      nU_part1[2] = -nU_part1[2];
    }
    *dphi_part1 = pi;
  }

  REAL nnorm = sqrt(nU_part1[0] * nU_part1[0] + nU_part1[1] * nU_part1[1] + nU_part1[2] * nU_part1[2]);
  if (nnorm < 1e-15) {
    fprintf(stderr, "ERROR in %s: could not recover nonzero rotation axis.\n", __func__);
    exit(1);
  }
  nU_part1[0] /= nnorm;
  nU_part1[1] /= nnorm;
  nU_part1[2] /= nnorm;

  // Debug reconstruction check: R_check(n,phi) should match input R.
  const REAL dphi_local = *dphi_part1;

  REAL Rcheck[3][3];
  const REAL tmp0 = cos(dphi_local);
  const REAL tmp2 = sin(dphi_local);
  const REAL tmp1 = 1 - tmp0;
  const REAL tmp5 = nU_part1[0] * nU_part1[1] * tmp1;
  const REAL tmp7 = nU_part1[0] * nU_part1[2] * tmp1;
  Rcheck[0][0] = ((nU_part1[0]) * (nU_part1[0])) * tmp1 + tmp0;
  Rcheck[0][1] = -nU_part1[2] * tmp2 + tmp5;
  Rcheck[0][2] = nU_part1[1] * tmp2 + tmp7;
  Rcheck[1][0] = nU_part1[2] * tmp2 + tmp5;
  Rcheck[1][1] = ((nU_part1[1]) * (nU_part1[1])) * tmp1 + tmp0;
  Rcheck[1][2] = -nU_part1[0] * tmp2 + nU_part1[1] * nU_part1[2] * tmp1;
  Rcheck[2][0] = -nU_part1[1] * tmp2 + tmp7;
  Rcheck[2][1] = nU_part1[0] * tmp2 + nU_part1[1] * nU_part1[2] * tmp1;
  Rcheck[2][2] = ((nU_part1[2]) * (nU_part1[2])) * tmp1 + tmp0;

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
            __func__, max_abs_err, *dphi_part1, nU_part1[0], nU_part1[1], nU_part1[2]);
    exit(1);
  }

  // Second output slot is deterministic identity.
  nU_part2[0] = 1.0;
  nU_part2[1] = 0.0;
  nU_part2[2] = 0.0;
  *dphi_part2 = 0.0;
} // END FUNCTION unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame

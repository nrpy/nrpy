#include "BHaH_defines.h"

/**
 * @brief Validate cumulative hats as a right-handed orthonormal basis; optionally repair.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 *
 * Validation triggers if any of:
 * - |x·y|, |x·z|, |y·z| > 1e-12
 * - |||x||-1|, |||y||-1|, |||z||-1| > 1e-12
 * - det(R) < 0.999999999999, where columns of R are (xhat,yhat,zhat)
 *
 * Deterministic repair policy when @p do_fix != 0:
 * 1) x <- normalize(x)
 * 2) y <- y - (x·y)x; y <- normalize(y)
 * 3) z <- x × y; z <- normalize(z)
 * 4) y <- z × x; y <- normalize(y)
 * 5) if det(R) < 0 then z <- -z and y <- z × x; y <- normalize(y)
 *
 * If post-fix invariants still fail, this routine aborts with diagnostics.
 *
 * @param[in,out] xhatU Candidate xhat vector (fixed basis).
 * @param[in,out] yhatU Candidate yhat vector (fixed basis).
 * @param[in,out] zhatU Candidate zhat vector (fixed basis).
 * @param[in] do_fix If nonzero, apply deterministic repair; otherwise report invalidity.
 * @return 0 if valid (possibly after repair), 1 if invalid and not repaired.
 */
int so3_validate_and_optionally_fix_hats(REAL xhatU[3], REAL yhatU[3], REAL zhatU[3], const int do_fix) {
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;
  const REAL tiny = 1e-300;

  REAL x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  REAL x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  REAL y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];

  REAL x_norm = sqrt(xhatU[0] * xhatU[0] + xhatU[1] * xhatU[1] + xhatU[2] * xhatU[2]);
  REAL y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  REAL z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);

  const REAL cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  REAL detR = xhatU[0] * cross_yz0 + xhatU[1] * cross_yz1 + xhatU[2] * cross_yz2;

  const int needs_fix = (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol || fabs(x_norm - 1.0) > norm_tol ||
                         fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (!needs_fix)
    return 0;

  if (!do_fix) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected and do_fix=0. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
    return 1;
  }

  // Deterministic Gram-Schmidt repair and right-handedness enforcement.
  if (x_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: xhat norm too small during fix (%.17e).\n", __func__, x_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    xhatU[i] /= x_norm;

  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  for (int i = 0; i < 3; i++)
    yhatU[i] -= x_dot_y * xhatU[i];

  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  if (y_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: yhat norm too small during fix (%.17e).\n", __func__, y_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    yhatU[i] /= y_norm;

  zhatU[0] = xhatU[1] * yhatU[2] - xhatU[2] * yhatU[1];
  zhatU[1] = xhatU[2] * yhatU[0] - xhatU[0] * yhatU[2];
  zhatU[2] = xhatU[0] * yhatU[1] - xhatU[1] * yhatU[0];

  z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);
  if (z_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: zhat norm too small during fix (%.17e).\n", __func__, z_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    zhatU[i] /= z_norm;

  yhatU[0] = zhatU[1] * xhatU[2] - zhatU[2] * xhatU[1];
  yhatU[1] = zhatU[2] * xhatU[0] - zhatU[0] * xhatU[2];
  yhatU[2] = zhatU[0] * xhatU[1] - zhatU[1] * xhatU[0];

  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  if (y_norm <= tiny) {
    fprintf(stderr, "ERROR in %s: yhat norm too small after z cross x (%.17e).\n", __func__, y_norm);
    exit(1);
  }
  for (int i = 0; i < 3; i++)
    yhatU[i] /= y_norm;

  const REAL post_cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL post_cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL post_cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  detR = xhatU[0] * post_cross_yz0 + xhatU[1] * post_cross_yz1 + xhatU[2] * post_cross_yz2;

  if (detR < 0.0) {
    for (int i = 0; i < 3; i++)
      zhatU[i] = -zhatU[i];

    yhatU[0] = zhatU[1] * xhatU[2] - zhatU[2] * xhatU[1];
    yhatU[1] = zhatU[2] * xhatU[0] - zhatU[0] * xhatU[2];
    yhatU[2] = zhatU[0] * xhatU[1] - zhatU[1] * xhatU[0];

    y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
    if (y_norm <= tiny) {
      fprintf(stderr, "ERROR in %s: yhat norm too small in right-handedness fix (%.17e).\n", __func__, y_norm);
      exit(1);
    }
    for (int i = 0; i < 3; i++)
      yhatU[i] /= y_norm;
  }

  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];
  x_norm = sqrt(xhatU[0] * xhatU[0] + xhatU[1] * xhatU[1] + xhatU[2] * xhatU[2]);
  y_norm = sqrt(yhatU[0] * yhatU[0] + yhatU[1] * yhatU[1] + yhatU[2] * yhatU[2]);
  z_norm = sqrt(zhatU[0] * zhatU[0] + zhatU[1] * zhatU[1] + zhatU[2] * zhatU[2]);

  const REAL final_cross_yz0 = yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1];
  const REAL final_cross_yz1 = yhatU[2] * zhatU[0] - yhatU[0] * zhatU[2];
  const REAL final_cross_yz2 = yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0];
  detR = xhatU[0] * final_cross_yz0 + xhatU[1] * final_cross_yz1 + xhatU[2] * final_cross_yz2;

  const int still_invalid = (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol || fabs(x_norm - 1.0) > norm_tol ||
                             fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (still_invalid) {
    fprintf(stderr,
            "ERROR in %s: post-fix hats remain invalid. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
    exit(1);
  }

  return 0;
} // END FUNCTION so3_validate_and_optionally_fix_hats

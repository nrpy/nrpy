#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

/**
 * @brief Map a Cartesian vector from rotating-frame coordinates to fixed-frame coordinates.
 *
 * Convention:
 * - R maps rotating-frame components -> fixed-frame components.
 * - R[:,0]=xhat, R[:,1]=yhat, R[:,2]=zhat (in fixed basis).
 * - v_fixed = R v_rot, v_rot = R^T v_fixed.
 * - T_fixed = R T_rot R^T, T_rot = R^T T_fixed R.
 * - DeltaR_dst_from_src = R_dst^T R_src.
 * - C layout statement for helper calls: R[i][j] is row i, column j.
 * - Einstein notation for this routine:
 *   x^i_fixed = R^i{}_j x^j_rot.
 * - This routine applies x_fixed = R * x_rot directly, i.e. it does not
 *   pass through an axis-angle decomposition.
 *
 * @param[in] commondata Commondata structure with cumulative regrid basis vectors.
 * @param[in,out] xCart Cartesian point/vector in rotating basis on input,
 *                     overwritten with fixed-basis components on output.
 */
void unrotate_xCart_to_fixed_frame(const commondata_struct commondata, REAL xCart[3]) {
  const REAL ortho_tol = 1e-12;
  const REAL norm_tol = 1e-12;
  const REAL det_tol = 0.999999999999;

  REAL xhatU[3] = {commondata.cumulative_regrid_xhatU[0], commondata.cumulative_regrid_xhatU[1], commondata.cumulative_regrid_xhatU[2]};
  REAL yhatU[3] = {commondata.cumulative_regrid_yhatU[0], commondata.cumulative_regrid_yhatU[1], commondata.cumulative_regrid_yhatU[2]};
  REAL zhatU[3] = {commondata.cumulative_regrid_zhatU[0], commondata.cumulative_regrid_zhatU[1], commondata.cumulative_regrid_zhatU[2]};
  REAL x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR;
  x_dot_y = xhatU[0] * yhatU[0] + xhatU[1] * yhatU[1] + xhatU[2] * yhatU[2];
  x_dot_z = xhatU[0] * zhatU[0] + xhatU[1] * zhatU[1] + xhatU[2] * zhatU[2];
  y_dot_z = yhatU[0] * zhatU[0] + yhatU[1] * zhatU[1] + yhatU[2] * zhatU[2];
  x_norm = sqrt(((xhatU[0]) * (xhatU[0])) + ((xhatU[1]) * (xhatU[1])) + ((xhatU[2]) * (xhatU[2])));
  y_norm = sqrt(((yhatU[0]) * (yhatU[0])) + ((yhatU[1]) * (yhatU[1])) + ((yhatU[2]) * (yhatU[2])));
  z_norm = sqrt(((zhatU[0]) * (zhatU[0])) + ((zhatU[1]) * (zhatU[1])) + ((zhatU[2]) * (zhatU[2])));
  detR = xhatU[0] * (yhatU[1] * zhatU[2] - yhatU[2] * zhatU[1]) + xhatU[1] * (-yhatU[0] * zhatU[2] + yhatU[2] * zhatU[0]) +
         xhatU[2] * (yhatU[0] * zhatU[1] - yhatU[1] * zhatU[0]);

  const int hats_invalid = (fabs(x_dot_y) > ortho_tol || fabs(x_dot_z) > ortho_tol || fabs(y_dot_z) > ortho_tol || fabs(x_norm - 1.0) > norm_tol ||
                            fabs(y_norm - 1.0) > norm_tol || fabs(z_norm - 1.0) > norm_tol || detR < det_tol);

  if (hats_invalid) {
    fprintf(stderr,
            "ERROR in %s: invalid hats detected at consuming boundary. "
            "dots=(%.17e, %.17e, %.17e) norms=(%.17e, %.17e, %.17e) det=%.17e\n",
            __func__, x_dot_y, x_dot_z, y_dot_z, x_norm, y_norm, z_norm, detR);
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

  // Alias-safe copy: permits xCart output to overwrite xCart input.
  const REAL x_in[3] = {xCart[0], xCart[1], xCart[2]};
  xCart[0] = R[0][0] * x_in[0] + R[0][1] * x_in[1] + R[0][2] * x_in[2];
  xCart[1] = R[1][0] * x_in[0] + R[1][1] * x_in[1] + R[1][2] * x_in[2];
  xCart[2] = R[2][0] * x_in[0] + R[2][1] * x_in[1] + R[2][2] * x_in[2];
} // END FUNCTION unrotate_xCart_to_fixed_frame

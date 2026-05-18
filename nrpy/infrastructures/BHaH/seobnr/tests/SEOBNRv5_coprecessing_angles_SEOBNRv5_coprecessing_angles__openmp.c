#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"

#include <gsl/gsl_spline.h>

typedef struct {
  REAL w;
  REAL x;
  REAL y;
  REAL z;
} quat_t;

static inline REAL seobnr_clamp_real(const REAL x, const REAL xmin, const REAL xmax) { return x < xmin ? xmin : (x > xmax ? xmax : x); }

static inline REAL seobnr_vec_dot3(const REAL ax, const REAL ay, const REAL az, const REAL bx, const REAL by, const REAL bz) {
  return ax * bx + ay * by + az * bz;
}

static inline REAL seobnr_vec_norm3(const REAL x, const REAL y, const REAL z) { return sqrt(x * x + y * y + z * z); }

static inline void seobnr_normalize3(REAL *restrict x, REAL *restrict y, REAL *restrict z) {
  const REAL norm = seobnr_vec_norm3(*x, *y, *z);
  if (norm > 0.0) {
    *x /= norm;
    *y /= norm;
    *z /= norm;
  }
}

static inline quat_t quat_normalize(quat_t q) {
  const REAL norm = sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
  if (norm > 0.0) {
    q.w /= norm;
    q.x /= norm;
    q.y /= norm;
    q.z /= norm;
  }
  return q;
}

static inline quat_t quat_conjugate(const quat_t q) { return (quat_t){q.w, -q.x, -q.y, -q.z}; }

static inline quat_t quat_negate(const quat_t q) { return (quat_t){-q.w, -q.x, -q.y, -q.z}; }

static inline quat_t quat_scale(const quat_t q, const REAL scale) { return (quat_t){q.w * scale, q.x * scale, q.y * scale, q.z * scale}; }

static inline REAL quat_dot(const quat_t a, const quat_t b) { return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z; }

static inline quat_t quat_multiply(const quat_t a, const quat_t b) {
  return (quat_t){a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z, a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
                  a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x, a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

static inline quat_t quat_exp_z(const REAL half_angle) { return (quat_t){cos(half_angle), 0.0, 0.0, sin(half_angle)}; }

static inline quat_t quat_sqrt_unit(quat_t q) {
  q = quat_normalize(q);
  const REAL w = seobnr_clamp_real(q.w, -1.0, 1.0);
  const REAL scalar = sqrt(0.5 * (1.0 + w));
  if (scalar > 1e-14) {
    const REAL scale = 0.5 / scalar;
    return quat_normalize((quat_t){scalar, q.x * scale, q.y * scale, q.z * scale});
  }
  const REAL vnorm = sqrt(q.x * q.x + q.y * q.y + q.z * q.z);
  if (vnorm > 0.0) {
    return (quat_t){0.0, q.x / vnorm, q.y / vnorm, q.z / vnorm};
  }
  return (quat_t){0.0, 1.0, 0.0, 0.0};
}

/**
 * Convert a rotation matrix to a unit quaternion.
 *
 * @param[in] rot_mat Rotation matrix.
 * @return Unit quaternion representing the same rotation.
 */
static inline quat_t quat_from_rotation_matrix(const REAL rot_mat[3][3]) {
  quat_t q;
  const REAL trace = rot_mat[0][0] + rot_mat[1][1] + rot_mat[2][2];
  if (trace > 0.0) {
    const REAL s = 2.0 * sqrt(trace + 1.0);
    q.w = 0.25 * s;
    q.x = (rot_mat[2][1] - rot_mat[1][2]) / s;
    q.y = (rot_mat[0][2] - rot_mat[2][0]) / s;
    q.z = (rot_mat[1][0] - rot_mat[0][1]) / s;
  } else if (rot_mat[0][0] > rot_mat[1][1] && rot_mat[0][0] > rot_mat[2][2]) {
    const REAL s = 2.0 * sqrt(1.0 + rot_mat[0][0] - rot_mat[1][1] - rot_mat[2][2]);
    q.w = (rot_mat[2][1] - rot_mat[1][2]) / s;
    q.x = 0.25 * s;
    q.y = (rot_mat[0][1] + rot_mat[1][0]) / s;
    q.z = (rot_mat[0][2] + rot_mat[2][0]) / s;
  } else if (rot_mat[1][1] > rot_mat[2][2]) {
    const REAL s = 2.0 * sqrt(1.0 + rot_mat[1][1] - rot_mat[0][0] - rot_mat[2][2]);
    q.w = (rot_mat[0][2] - rot_mat[2][0]) / s;
    q.x = (rot_mat[0][1] + rot_mat[1][0]) / s;
    q.y = 0.25 * s;
    q.z = (rot_mat[1][2] + rot_mat[2][1]) / s;
  } else {
    const REAL s = 2.0 * sqrt(1.0 + rot_mat[2][2] - rot_mat[0][0] - rot_mat[1][1]);
    q.w = (rot_mat[1][0] - rot_mat[0][1]) / s;
    q.x = (rot_mat[0][2] + rot_mat[2][0]) / s;
    q.y = (rot_mat[1][2] + rot_mat[2][1]) / s;
    q.z = 0.25 * s;
  }
  return quat_normalize(q);
} // END FUNCTION: quat_from_rotation_matrix

/**
 * Convert a unit quaternion to its rotation matrix.
 *
 * @param[in] q Input quaternion.
 * @param[out] rot_mat Rotation matrix.
 */
static inline void quat_to_rotation_matrix(quat_t q, REAL rot_mat[3][3]) {
  q = quat_normalize(q);
  const REAL ww = q.w * q.w;
  const REAL xx = q.x * q.x;
  const REAL yy = q.y * q.y;
  const REAL zz = q.z * q.z;
  const REAL wx = q.w * q.x;
  const REAL wy = q.w * q.y;
  const REAL wz = q.w * q.z;
  const REAL xy = q.x * q.y;
  const REAL xz = q.x * q.z;
  const REAL yz = q.y * q.z;

  rot_mat[0][0] = ww + xx - yy - zz;
  rot_mat[0][1] = 2.0 * (xy - wz);
  rot_mat[0][2] = 2.0 * (xz + wy);
  rot_mat[1][0] = 2.0 * (xy + wz);
  rot_mat[1][1] = ww - xx + yy - zz;
  rot_mat[1][2] = 2.0 * (yz - wx);
  rot_mat[2][0] = 2.0 * (xz - wy);
  rot_mat[2][1] = 2.0 * (yz + wx);
  rot_mat[2][2] = ww - xx - yy + zz;
} // END FUNCTION: quat_to_rotation_matrix

static inline void quat_image_of_z(const quat_t q_in, REAL *restrict zx, REAL *restrict zy, REAL *restrict zz) {
  REAL rot_mat[3][3];
  quat_to_rotation_matrix(q_in, rot_mat);
  *zx = rot_mat[0][2];
  *zy = rot_mat[1][2];
  *zz = rot_mat[2][2];
  seobnr_normalize3(zx, zy, zz);
}

static inline REAL angle_between_unit_vectors(const REAL ax, const REAL ay, const REAL az, const REAL bx, const REAL by, const REAL bz) {
  return acos(seobnr_clamp_real(seobnr_vec_dot3(ax, ay, az, bx, by, bz), -1.0, 1.0));
}

/**
 * Extract Z-Y-Z Euler angles from a unit quaternion.
 *
 * @param[in] q Input quaternion.
 * @param[out] alpha First Euler angle.
 * @param[out] beta Second Euler angle.
 * @param[out] gamma Third Euler angle.
 */
static inline void extract_zyz_euler_angles_from_quat(const quat_t q, REAL *restrict alpha, REAL *restrict beta, REAL *restrict gamma) {
  REAL rot_mat[3][3];
  quat_to_rotation_matrix(q, rot_mat);
  const REAL R33 = seobnr_clamp_real(rot_mat[2][2], -1.0, 1.0);
  *beta = acos(R33);
  if (R33 > 1.0 - 1e-14) {
    *alpha = atan2(rot_mat[1][0], rot_mat[0][0]);
    *gamma = 0.0;
  } else if (R33 < -1.0 + 1e-14) {
    *alpha = atan2(-rot_mat[1][0], -rot_mat[0][0]);
    *gamma = 0.0;
  } else {
    *alpha = atan2(rot_mat[1][2], rot_mat[0][2]);
    *gamma = atan2(rot_mat[2][1], -rot_mat[2][0]);
  }
} // END FUNCTION: extract_zyz_euler_angles_from_quat

static inline void sample_inspiral_time_and_omega(const commondata_struct *restrict commondata, const size_t i, REAL *restrict time_M,
                                                  REAL *restrict omega) {
  if (i < commondata->nsteps_low) {
    *time_M = commondata->dynamics_low[IDX(i, TIME)];
    *omega = commondata->dynamics_low[IDX(i, OMEGA)];
  } else {
    const size_t j = i - commondata->nsteps_low;
    *time_M = commondata->dynamics_fine[IDX(j, TIME)];
    *omega = commondata->dynamics_fine[IDX(j, OMEGA)];
  }
}

/**
 * Interpolate the fine-dynamics orbital frequency at a target time.
 *
 * @param[in] commondata Common data struct containing dynamics arrays.
 * @param[in] target_time Time at which to evaluate the frequency.
 * @return Interpolated orbital frequency.
 */
static inline REAL interpolate_fine_omega_at_time(const commondata_struct *restrict commondata, const REAL target_time) {
  if (commondata->nsteps_fine == 0) {
    return commondata->nsteps_low > 0 ? commondata->dynamics_low[IDX(commondata->nsteps_low - 1, OMEGA)] : commondata->initial_omega;
  }
  if (commondata->nsteps_fine < 2) {
    return commondata->dynamics_fine[IDX(0, OMEGA)];
  }
  for (size_t i = 0; i + 1 < commondata->nsteps_fine; i++) {
    const REAL t0 = commondata->dynamics_fine[IDX(i, TIME)];
    const REAL t1 = commondata->dynamics_fine[IDX(i + 1, TIME)];
    if (target_time <= t1 || i + 2 == commondata->nsteps_fine) {
      const REAL omega0 = commondata->dynamics_fine[IDX(i, OMEGA)];
      const REAL omega1 = commondata->dynamics_fine[IDX(i + 1, OMEGA)];
      const REAL denom = t1 - t0;
      if (fabs(denom) < 1e-15) {
        return omega0;
      }
      const REAL frac = (target_time - t0) / denom;
      return omega0 + frac * (omega1 - omega0);
    }
  }
  return commondata->dynamics_fine[IDX(commondata->nsteps_fine - 1, OMEGA)];
} // END FUNCTION: interpolate_fine_omega_at_time

/**
 * Compute not-a-knot cubic spline second derivatives at knots.
 *
 * @param[in] n Number of samples.
 * @param[in] x Strictly increasing sample coordinates.
 * @param[in] y Sample values.
 * @param[out] second_deriv Spline second derivatives.
 */
static void notaknot_spline_second_derivatives(const size_t n, const REAL *restrict x, const REAL *restrict y, REAL *restrict second_deriv) {
  for (size_t i = 0; i < n; i++) {
    second_deriv[i] = 0.0;
  }
  if (n < 4) {
    return;
  }

  const size_t m = n - 2;
  REAL *restrict lower = (REAL *)calloc(m, sizeof(REAL));
  REAL *restrict diag = (REAL *)calloc(m, sizeof(REAL));
  REAL *restrict upper = (REAL *)calloc(m, sizeof(REAL));
  REAL *restrict rhs = (REAL *)calloc(m, sizeof(REAL));
  if (lower == NULL || diag == NULL || upper == NULL || rhs == NULL) {
    fprintf(stderr, "Error: notaknot_spline_second_derivatives() malloc failed\n");
    exit(1);
  }

  const REAL h0 = x[1] - x[0];
  const REAL h1 = x[2] - x[1];
  const REAL d0 = (y[1] - y[0]) / h0;
  const REAL d1 = (y[2] - y[1]) / h1;
  diag[0] = h0 + 2.0 * h1;
  upper[0] = h1 - h0;
  rhs[0] = 6.0 * h1 * (d1 - d0) / (h0 + h1);

  for (size_t j = 1; j + 1 < m; j++) {
    const size_t i = j + 1;
    const REAL hm = x[i] - x[i - 1];
    const REAL hp = x[i + 1] - x[i];
    const REAL dm = (y[i] - y[i - 1]) / hm;
    const REAL dp = (y[i + 1] - y[i]) / hp;
    lower[j] = hm;
    diag[j] = 2.0 * (hm + hp);
    upper[j] = hp;
    rhs[j] = 6.0 * (dp - dm);
  }

  const REAL ha = x[n - 2] - x[n - 3];
  const REAL hb = x[n - 1] - x[n - 2];
  const REAL da = (y[n - 2] - y[n - 3]) / ha;
  const REAL db = (y[n - 1] - y[n - 2]) / hb;
  lower[m - 1] = ha - hb;
  diag[m - 1] = 2.0 * ha + hb;
  rhs[m - 1] = 6.0 * ha * (db - da) / (ha + hb);

  for (size_t i = 1; i < m; i++) {
    const REAL factor = lower[i] / diag[i - 1];
    diag[i] -= factor * upper[i - 1];
    rhs[i] -= factor * rhs[i - 1];
  }
  second_deriv[n - 2] = rhs[m - 1] / diag[m - 1];
  for (size_t j = m - 1; j-- > 0;) {
    second_deriv[j + 1] = (rhs[j] - upper[j] * second_deriv[j + 2]) / diag[j];
  }
  second_deriv[0] = ((h0 + h1) * second_deriv[1] - h0 * second_deriv[2]) / h1;
  second_deriv[n - 1] = ((ha + hb) * second_deriv[n - 2] - hb * second_deriv[n - 3]) / ha;

  free(lower);
  free(diag);
  free(upper);
  free(rhs);
} // END FUNCTION: notaknot_spline_second_derivatives

/**
 * Evaluate first derivatives of a cubic spline at its knots.
 *
 * @param[in] n Number of samples.
 * @param[in] x Strictly increasing sample coordinates.
 * @param[in] y Sample values.
 * @param[in] second_deriv Spline second derivatives.
 * @param[out] deriv Spline first derivatives at knots.
 */
static void notaknot_spline_derivatives_at_knots(const size_t n, const REAL *restrict x, const REAL *restrict y, const REAL *restrict second_deriv,
                                                 REAL *restrict deriv) {
  if (n < 2) {
    if (n == 1) {
      deriv[0] = 0.0;
    }
    return;
  }
  for (size_t i = 0; i + 1 < n; i++) {
    const REAL h = x[i + 1] - x[i];
    const REAL delta = (y[i + 1] - y[i]) / h;
    deriv[i] = delta - h * (2.0 * second_deriv[i] + second_deriv[i + 1]) / 6.0;
  }
  const REAL h = x[n - 1] - x[n - 2];
  const REAL delta = (y[n - 1] - y[n - 2]) / h;
  deriv[n - 1] = delta + h * (second_deriv[n - 2] + 2.0 * second_deriv[n - 1]) / 6.0;
} // END FUNCTION: notaknot_spline_derivatives_at_knots

/**
 * Integrate a cubic spline and return prefix integrals at knots.
 *
 * @param[in] n Number of samples.
 * @param[in] x Strictly increasing sample coordinates.
 * @param[in] y Sample values.
 * @param[in] second_deriv Spline second derivatives.
 * @param[out] integral Prefix integral values.
 */
static void notaknot_spline_integral_prefix(const size_t n, const REAL *restrict x, const REAL *restrict y, const REAL *restrict second_deriv,
                                            REAL *restrict integral) {
  if (n == 0)
    return;
  integral[0] = 0.0;
  for (size_t i = 0; i + 1 < n; i++) {
    const REAL h = x[i + 1] - x[i];
    const REAL area = 0.5 * h * (y[i] + y[i + 1]) - h * h * h * (second_deriv[i] + second_deriv[i + 1]) / 24.0;
    integral[i + 1] = integral[i] + area;
  }
} // END FUNCTION: notaknot_spline_integral_prefix

/**
 * Compute the pyseobnr-style minimal-rotation quaternion track.
 *
 * @param[in] n Number of samples.
 * @param[in] time_arr Time samples.
 * @param[in] omega_arr Orbital-frequency samples.
 * @param[in] LN_J_x_arr J-frame LN x samples.
 * @param[in] LN_J_y_arr J-frame LN y samples.
 * @param[in] LN_J_z_arr J-frame LN z samples.
 * @param[out] q_pre_arr Minimal-rotation quaternions.
 */
static void compute_py_style_coprecessing_rotator(const size_t n, const REAL *restrict time_arr, const REAL *restrict omega_arr,
                                                  const REAL *restrict LN_J_x_arr, const REAL *restrict LN_J_y_arr, const REAL *restrict LN_J_z_arr,
                                                  quat_t *restrict q_pre_arr) {
  if (n == 0)
    return;
  quat_t *restrict q_domega = (quat_t *)malloc(n * sizeof(quat_t));
  REAL *restrict q_component = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict spline_second = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict spline_deriv = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict omega_second = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict domega_dt = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict halfgammadot = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict halfgamma = (REAL *)malloc(n * sizeof(REAL));
  REAL *restrict halfgamma_second = (REAL *)malloc(n * sizeof(REAL));
  if (q_domega == NULL || q_component == NULL || spline_second == NULL || spline_deriv == NULL || omega_second == NULL || domega_dt == NULL ||
      halfgammadot == NULL || halfgamma == NULL || halfgamma_second == NULL) {
    fprintf(stderr, "Error: compute_py_style_coprecessing_rotator() malloc failed\n");
    exit(1);
  } // END IF: minimal-rotation workspace allocation failed

  REAL min_abs_z_alignment = 1.0;
  REAL min_abs_x_alignment = 1.0;
  for (size_t i = 0; i < n; i++) {
    const REAL z_alignment = fabs(fabs(LN_J_z_arr[i]) - 1.0);
    const REAL x_alignment = fabs(fabs(LN_J_x_arr[i]) - 1.0);
    if (z_alignment < min_abs_z_alignment)
      min_abs_z_alignment = z_alignment;
    if (x_alignment < min_abs_x_alignment)
      min_abs_x_alignment = x_alignment;
  } // END LOOP: for i over LN_J alignment checks

  quat_t aux_quat = {0.0, 0.0, 0.0, 0.0};
  const int use_aux_quat = min_abs_z_alignment < 1e-3;
  if (use_aux_quat) {
    aux_quat = (quat_t){0.0, 1.0, 0.0, 0.0};
    if (min_abs_x_alignment < 1e-3)
      aux_quat = (quat_t){0.0, 0.0, 1.0, 0.0};
  } // END IF: auxiliary quaternion needed near z alignment
  const quat_t z_quat = {0.0, 0.0, 0.0, 1.0};
  for (size_t i = 0; i < n; i++) {
    const quat_t LN_quat = {0.0, LN_J_x_arr[i], LN_J_y_arr[i], LN_J_z_arr[i]};
    if (use_aux_quat) {
      const quat_t step2 = quat_sqrt_unit(quat_negate(quat_multiply(LN_quat, aux_quat)));
      const quat_t step1 = quat_sqrt_unit(quat_negate(quat_multiply(aux_quat, z_quat)));
      q_pre_arr[i] = quat_multiply(step2, step1);
    } else {
      q_pre_arr[i] = quat_sqrt_unit(quat_negate(quat_multiply(LN_quat, z_quat)));
    } // END ELSE: direct z-to-LN construction is well conditioned
    q_pre_arr[i] = quat_normalize(q_pre_arr[i]);
    if (i > 0 && quat_dot(q_pre_arr[i], q_pre_arr[i - 1]) < 0.0)
      q_pre_arr[i] = quat_negate(q_pre_arr[i]);
  } // END LOOP: for i over initial pyseobnr-style quaternions

  notaknot_spline_second_derivatives(n, time_arr, omega_arr, omega_second);
  notaknot_spline_derivatives_at_knots(n, time_arr, omega_arr, omega_second, domega_dt);
  for (size_t iter = 0; iter < 2; iter++) {
    for (size_t i = 0; i < n; i++) {
      q_domega[i] = (quat_t){0.0, 0.0, 0.0, 0.0};
    } // END LOOP: for i over quaternion omega-derivative initialization
    for (size_t comp = 0; comp < 4; comp++) {
      for (size_t i = 0; i < n; i++) {
        q_component[i] = comp == 0 ? q_pre_arr[i].w : (comp == 1 ? q_pre_arr[i].x : (comp == 2 ? q_pre_arr[i].y : q_pre_arr[i].z));
      } // END LOOP: for i over one quaternion component
      notaknot_spline_second_derivatives(n, omega_arr, q_component, spline_second);
      notaknot_spline_derivatives_at_knots(n, omega_arr, q_component, spline_second, spline_deriv);
      for (size_t i = 0; i < n; i++) {
        if (comp == 0) {
          q_domega[i].w = spline_deriv[i];
        } else if (comp == 1) {
          q_domega[i].x = spline_deriv[i];
        } else if (comp == 2) {
          q_domega[i].y = spline_deriv[i];
        } else {
          q_domega[i].z = spline_deriv[i];
        } // END ELSE: storing z component omega derivative
      } // END LOOP: for i over quaternion omega-derivative samples
    } // END LOOP: for comp over quaternion components
    for (size_t i = 0; i < n; i++) {
      const quat_t qdot = quat_scale(q_domega[i], domega_dt[i]);
      const quat_t tmp = quat_multiply(quat_multiply(qdot, z_quat), quat_conjugate(q_pre_arr[i]));
      halfgammadot[i] = tmp.w;
    } // END LOOP: for i over minimal-rotation gamma derivatives
    notaknot_spline_second_derivatives(n, time_arr, halfgammadot, halfgamma_second);
    notaknot_spline_integral_prefix(n, time_arr, halfgammadot, halfgamma_second, halfgamma);
    for (size_t i = 0; i < n; i++) {
      q_pre_arr[i] = quat_multiply(q_pre_arr[i], quat_exp_z(halfgamma[i]));
      q_pre_arr[i] = quat_normalize(q_pre_arr[i]);
    } // END LOOP: for i over minimal-rotation gamma correction
  } // END LOOP: for iter over minimal-rotation correction passes

  free(q_domega);
  free(q_component);
  free(spline_second);
  free(spline_deriv);
  free(omega_second);
  free(domega_dt);
  free(halfgammadot);
  free(halfgamma);
  free(halfgamma_second);
} // END FUNCTION: compute_py_style_coprecessing_rotator

static inline REAL sanitize_omega_for_spin_splines(const REAL omega, const REAL omega_min, const REAL omega_max, const REAL abs_tol) {
  if (omega < omega_min && fabs(omega - omega_min) <= abs_tol)
    return omega_min;
  if (omega > omega_max && fabs(omega - omega_max) <= abs_tol)
    return omega_max;
  return omega;
} // END FUNCTION: sanitize_omega_for_spin_splines

static inline void check_omega_range(const char *restrict label, const REAL omega_val, const REAL omega_min, const REAL omega_max,
                                     const REAL time_val, const size_t idx_val) {
  if (omega_val < omega_min || omega_val > omega_max) {
    fprintf(stderr, "SEOBNRv5_coprecessing_angles: %s out of range at i=%zu, t=%.15e: %.15e not in [%.15e, %.15e]\n", label, idx_val,
            (double)time_val, (double)omega_val, (double)omega_min, (double)omega_max);
    exit(1);
  } // END IF: omega sample lies outside spin-spline range
} // END FUNCTION: check_omega_range

/**
 * Compute the inspiral J->P Euler-angle arrays for SEOBNRv5 coprecessing rotations.
 *
 * This routine follows the inspiral frame construction of Ramos-Buades et al.
 * (`SEOBNRv5PHM`, arXiv:2303.18046), especially Sec. III A, Eq. (15), and
 * Eqs. (16)-(17). The minimal-rotation condition is imposed using a quaternion
 * spline correction following Boyle, Owen, and Pfeiffer, Phys. Rev. D 84, 124011
 * (2011), arXiv:1110.2965, Eq. (18). The robust J-frame basis construction is
 * consistent with the earlier SEOBNR implementations of Babak, Taracchini,
 * and Buonanno, Phys. Rev. D 95, 024010 (2017), arXiv:1607.05661, and
 * Ossokine et al., Phys. Rev. D 102, 044055 (2020), arXiv:2004.09442.
 *
 * @param[in,out] commondata Common data struct updated with output Euler-angle arrays.
 */
void SEOBNRv5_coprecessing_angles(commondata_struct *restrict commondata) {
  const size_t n_low = commondata->nsteps_low;
  const size_t n_fine = commondata->nsteps_fine;
  const size_t n_insp = n_low + n_fine;
  FILE *fp_dbg = NULL;
  // Maintainer note: the orbital trajectory and spin-evolution splines are built by
  // separate numerical integrators. In practice their nominally identical start/end
  // frequencies can differ by O(1e-8) because of floating-point and adaptive-step
  // effects. We sanitize omega only within a tiny absolute tolerance so that we do
  // not fail spuriously at the spline boundaries. If larger excursions appear here,
  // the right long-term fix is to tighten the consistency of the two trajectories,
  // not to widen this tolerance.
  const REAL omega_boundary_tol = 1e-7;
  if (commondata->alpha_JP != NULL) {
    free(commondata->alpha_JP);
    commondata->alpha_JP = NULL;
  }
  if (commondata->beta_JP != NULL) {
    free(commondata->beta_JP);
    commondata->beta_JP = NULL;
  }
  if (commondata->gamma_JP != NULL) {
    free(commondata->gamma_JP);
    commondata->gamma_JP = NULL;
  }

  if (n_insp == 0) {
    if (fp_dbg != NULL)
      fclose(fp_dbg);
    return;
  } // END IF: no inspiral samples available for angle construction

  commondata->alpha_JP = (REAL *)malloc(n_insp * sizeof(REAL));
  commondata->beta_JP = (REAL *)malloc(n_insp * sizeof(REAL));
  commondata->gamma_JP = (REAL *)malloc(n_insp * sizeof(REAL));
  if (commondata->alpha_JP == NULL || commondata->beta_JP == NULL || commondata->gamma_JP == NULL) {
    fprintf(stderr, "Error: in SEOBNRv5_coprecessing_angles(), malloc() failed for Euler-angle arrays\n");
    exit(1);
  }

  // Evaluate J_f at the attachment time, matching the SEOBNRv5PHM construction
  // where the J-frame is tied to the attachment stage rather than the terminal
  // spin-dynamics sample.
  REAL omega_attach = interpolate_fine_omega_at_time(commondata, commondata->t_attach);
  omega_attach = sanitize_omega_for_spin_splines(omega_attach, commondata->omega_spin_min, commondata->omega_spin_max, omega_boundary_tol);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "attach %d %.15e %.15e %.15e %.15e\n", -1, commondata->t_attach, omega_attach, commondata->omega_spin_min,
            commondata->omega_spin_max);
    fprintf(fp_dbg, "ptr L_x %p %p\n", (void *)commondata->L_x.spline, (void *)commondata->L_x.acc);
    fprintf(fp_dbg, "ptr L_y %p %p\n", (void *)commondata->L_y.spline, (void *)commondata->L_y.acc);
    fprintf(fp_dbg, "ptr L_z %p %p\n", (void *)commondata->L_z.spline, (void *)commondata->L_z.acc);
    fprintf(fp_dbg, "ptr chi1_x %p %p\n", (void *)commondata->chi1_x_spline.spline, (void *)commondata->chi1_x_spline.acc);
    fprintf(fp_dbg, "ptr chi1_y %p %p\n", (void *)commondata->chi1_y_spline.spline, (void *)commondata->chi1_y_spline.acc);
    fprintf(fp_dbg, "ptr chi1_z %p %p\n", (void *)commondata->chi1_z_spline.spline, (void *)commondata->chi1_z_spline.acc);
    fprintf(fp_dbg, "ptr chi2_x %p %p\n", (void *)commondata->chi2_x_spline.spline, (void *)commondata->chi2_x_spline.acc);
    fprintf(fp_dbg, "ptr chi2_y %p %p\n", (void *)commondata->chi2_y_spline.spline, (void *)commondata->chi2_y_spline.acc);
    fprintf(fp_dbg, "ptr chi2_z %p %p\n", (void *)commondata->chi2_z_spline.spline, (void *)commondata->chi2_z_spline.acc);
    fflush(fp_dbg);
  }
  check_omega_range("omega_attach", omega_attach, commondata->omega_spin_min, commondata->omega_spin_max, commondata->t_attach, 0);

  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "about L_x %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL L_attach_x = gsl_spline_eval(commondata->L_x.spline, omega_attach, commondata->L_x.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done L_x %.15e\n", L_attach_x);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about L_y %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL L_attach_y = gsl_spline_eval(commondata->L_y.spline, omega_attach, commondata->L_y.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done L_y %.15e\n", L_attach_y);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about L_z %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL L_attach_z = gsl_spline_eval(commondata->L_z.spline, omega_attach, commondata->L_z.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done L_z %.15e\n", L_attach_z);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi1_x %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi1_attach_x = gsl_spline_eval(commondata->chi1_x_spline.spline, omega_attach, commondata->chi1_x_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi1_x %.15e\n", chi1_attach_x);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi1_y %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi1_attach_y = gsl_spline_eval(commondata->chi1_y_spline.spline, omega_attach, commondata->chi1_y_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi1_y %.15e\n", chi1_attach_y);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi1_z %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi1_attach_z = gsl_spline_eval(commondata->chi1_z_spline.spline, omega_attach, commondata->chi1_z_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi1_z %.15e\n", chi1_attach_z);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi2_x %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi2_attach_x = gsl_spline_eval(commondata->chi2_x_spline.spline, omega_attach, commondata->chi2_x_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi2_x %.15e\n", chi2_attach_x);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi2_y %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi2_attach_y = gsl_spline_eval(commondata->chi2_y_spline.spline, omega_attach, commondata->chi2_y_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi2_y %.15e\n", chi2_attach_y);
    fflush(fp_dbg);
    fprintf(fp_dbg, "about chi2_z %d %.15e\n", -1, omega_attach);
    fflush(fp_dbg);
  }
  const REAL chi2_attach_z = gsl_spline_eval(commondata->chi2_z_spline.spline, omega_attach, commondata->chi2_z_spline.acc);
  if (fp_dbg != NULL) {
    fprintf(fp_dbg, "done chi2_z %.15e\n", chi2_attach_z);
    fflush(fp_dbg);
  }
  // NRPy stores L in dimensionless units L/(M*mu) where M=1, so the stored L
  // is larger than pyseobnr's convention by 1/(M*mu) = 1/(m1*m2). Multiply by
  // mu = m1*m2 to restore the M*mu prefactor before forming J_f, matching the
  // pyseobnr convention: J_f = mu*L_vec + m1^2*chi1 + m2^2*chi2.
  const REAL mu_Jf = commondata->m1 * commondata->m2;
  commondata->J_f_x = mu_Jf * L_attach_x + commondata->m1 * commondata->m1 * chi1_attach_x + commondata->m2 * commondata->m2 * chi2_attach_x;
  commondata->J_f_y = mu_Jf * L_attach_y + commondata->m1 * commondata->m1 * chi1_attach_y + commondata->m2 * commondata->m2 * chi2_attach_y;
  commondata->J_f_z = mu_Jf * L_attach_z + commondata->m1 * commondata->m1 * chi1_attach_z + commondata->m2 * commondata->m2 * chi2_attach_z;
  const REAL J_f_norm = seobnr_vec_norm3(commondata->J_f_x, commondata->J_f_y, commondata->J_f_z);
  if (J_f_norm < 1e-15) {
    commondata->J_f_x = 0.0;
    commondata->J_f_y = 0.0;
    commondata->J_f_z = 1.0;
  } else {
    commondata->J_f_x /= J_f_norm;
    commondata->J_f_y /= J_f_norm;
    commondata->J_f_z /= J_f_norm;
  } // END ELSE: final angular momentum has a reliable direction

  // Build the J-frame basis from the final angular momentum direction, following
  // Eq. (15) of the SEOBNRv5PHM paper with the smooth x/y prescription used in
  // earlier SEOBNR precessing implementations.
  REAL e3J_x = commondata->J_f_x;
  REAL e3J_y = commondata->J_f_y;
  REAL e3J_z = commondata->J_f_z;
  seobnr_normalize3(&e3J_x, &e3J_y, &e3J_z);

  const REAL exdote3J = e3J_x;
  const REAL eydote3J = e3J_y;
  const REAL lambda_fac = 1.0 - fabs(exdote3J);
  REAL e1J_x, e1J_y, e1J_z;
  if (lambda_fac > 1e-4) {
    const REAL normfacx = sqrt(1.0 - exdote3J * exdote3J);
    e1J_x = (1.0 - exdote3J * e3J_x) / normfacx;
    e1J_y = (-exdote3J * e3J_y) / normfacx;
    e1J_z = (-exdote3J * e3J_z) / normfacx;
  } else if (lambda_fac < 1e-5) {
    const REAL normfacy = sqrt(1.0 - eydote3J * eydote3J);
    e1J_x = (-eydote3J * e3J_x) / normfacy;
    e1J_y = (1.0 - eydote3J * e3J_y) / normfacy;
    e1J_z = (-eydote3J * e3J_z) / normfacy;
  } else {
    const REAL weightx = (lambda_fac - 1e-5) / (1e-4 - 1e-5);
    const REAL weighty = 1.0 - weightx;
    const REAL normfacx = sqrt(1.0 - exdote3J * exdote3J);
    const REAL normfacy = sqrt(1.0 - eydote3J * eydote3J);
    e1J_x = weightx * (1.0 - exdote3J * e3J_x) / normfacx + weighty * (-eydote3J * e3J_x) / normfacy;
    e1J_y = weightx * (-exdote3J * e3J_y) / normfacx + weighty * (1.0 - eydote3J * e3J_y) / normfacy;
    e1J_z = weightx * (-exdote3J * e3J_z) / normfacx + weighty * (-eydote3J * e3J_z) / normfacy;
  } // END ELSE: blend x- and y-projection J-frame fallbacks
  seobnr_normalize3(&e1J_x, &e1J_y, &e1J_z);

  REAL e2J_x = e3J_y * e1J_z - e3J_z * e1J_y;
  REAL e2J_y = e3J_z * e1J_x - e3J_x * e1J_z;
  REAL e2J_z = e3J_x * e1J_y - e3J_y * e1J_x;
  seobnr_normalize3(&e2J_x, &e2J_y, &e2J_z);

  REAL *restrict time_arr = (REAL *)malloc(n_insp * sizeof(REAL));
  REAL *restrict omega_arr = (REAL *)malloc(n_insp * sizeof(REAL));
  REAL *restrict LN_J_x_arr = (REAL *)malloc(n_insp * sizeof(REAL));
  REAL *restrict LN_J_y_arr = (REAL *)malloc(n_insp * sizeof(REAL));
  REAL *restrict LN_J_z_arr = (REAL *)malloc(n_insp * sizeof(REAL));
  quat_t *restrict q_pre_arr = (quat_t *)malloc(n_insp * sizeof(quat_t));
  if (time_arr == NULL || omega_arr == NULL || LN_J_x_arr == NULL || LN_J_y_arr == NULL || LN_J_z_arr == NULL || q_pre_arr == NULL) {
    fprintf(stderr, "Error: in SEOBNRv5_coprecessing_angles(), malloc() failed for py-style rotation arrays\n");
    exit(1);
  }

  for (size_t i = 0; i < n_insp; i++) {
    REAL time_M = 0.0;
    REAL omega = 0.0;
    REAL lnhat_I_x = 0.0;
    REAL lnhat_I_y = 0.0;
    REAL lnhat_I_z = 1.0;
    sample_inspiral_time_and_omega(commondata, i, &time_M, &omega);
    omega = sanitize_omega_for_spin_splines(omega, commondata->omega_spin_min, commondata->omega_spin_max, omega_boundary_tol);
    check_omega_range("sample_i", omega, commondata->omega_spin_min, commondata->omega_spin_max, time_M, i);
    if (fp_dbg != NULL && (i < 5 || i + 5 >= n_insp)) {
      fprintf(fp_dbg, "sample %zu %.15e %.15e %.15e %.15e\n", i, time_M, omega, commondata->omega_spin_min, commondata->omega_spin_max);
      fflush(fp_dbg);
    }
    lnhat_I_x = gsl_spline_eval(commondata->lnhat_x.spline, omega, commondata->lnhat_x.acc);
    lnhat_I_y = gsl_spline_eval(commondata->lnhat_y.spline, omega, commondata->lnhat_y.acc);
    lnhat_I_z = gsl_spline_eval(commondata->lnhat_z.spline, omega, commondata->lnhat_z.acc);
    seobnr_normalize3(&lnhat_I_x, &lnhat_I_y, &lnhat_I_z);

    REAL LN_J_x = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e1J_x, e1J_y, e1J_z);
    REAL LN_J_y = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e2J_x, e2J_y, e2J_z);
    REAL LN_J_z = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e3J_x, e3J_y, e3J_z);
    seobnr_normalize3(&LN_J_x, &LN_J_y, &LN_J_z);
    time_arr[i] = time_M;
    omega_arr[i] = omega;
    LN_J_x_arr[i] = LN_J_x;
    LN_J_y_arr[i] = LN_J_y;
    LN_J_z_arr[i] = LN_J_z;
  } // END LOOP: for i over inspiral samples building J-frame inputs

  compute_py_style_coprecessing_rotator(n_insp, time_arr, omega_arr, LN_J_x_arr, LN_J_y_arr, LN_J_z_arr, q_pre_arr);

  for (size_t i = 0; i < n_insp; i++) {

    const quat_t q_total = quat_normalize(q_pre_arr[i]);
    extract_zyz_euler_angles_from_quat(q_total, &commondata->alpha_JP[i], &commondata->beta_JP[i], &commondata->gamma_JP[i]);

  } // END LOOP: for i over inspiral samples extracting Euler angles

  free(time_arr);
  free(omega_arr);
  free(LN_J_x_arr);
  free(LN_J_y_arr);
  free(LN_J_z_arr);
  free(q_pre_arr);
  if (fp_dbg != NULL) {
    fclose(fp_dbg);
  }
} // END FUNCTION: SEOBNRv5_coprecessing_angles

"""
Register CFunction for generating inspiral J->P Euler-angle arrays for SEOBNRv5.

This module constructs the inspiral co-precessing-frame Euler angles
``alpha_JP(t)``, ``beta_JP(t)``, and ``gamma_JP(t)`` from the precessing-spin
dynamics already stored in ``commondata``. The geometric setup follows
Sec. III A of Ramos-Buades et al., SEOBNRv5PHM (arXiv:2303.18046), especially
Eq. (15) for the J-frame and Eqs. (16)-(17) for the J->P rotation. The
minimal-rotation condition is imposed with a quaternion spline correction,
following Boyle, Owen, and Pfeiffer, Phys. Rev. D 84, 124011 (2011),
arXiv:1110.2965, Eq. (18).

The implementation is numerically consistent with the earlier SEOBNR
precessing-frame constructions described in Babak, Taracchini, and Buonanno,
Phys. Rev. D 95, 024010 (2017), arXiv:1607.05661, and Ossokine et al.,
Phys. Rev. D 102, 044055 (2020), arXiv:2004.09442.

Authors: Suchindram Dasgupta
         sd00113 at mix dot wvu dot edu
         Zachariah B. Etienne
         zachetie at gmail *dot com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_SEOBNRv5_coprecessing_angles(
    enable_forensic_diagnostics: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register CFunction for computing inspiral J->P Euler-angle arrays.

    The generated C function samples the precessing-spin dynamics on the
    inspiral waveform grid, builds the J-frame from the final angular momentum,
    enforces the minimal-rotation condition using a quaternion-based internal
    representation, and stores the resulting Euler angles in ``commondata``.

    :param enable_forensic_diagnostics: Whether to write development-only
        C-side validation dumps for pyseobnr comparisons.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> _ = register_CFunction_SEOBNRv5_coprecessing_angles()
    >>> generated_str = clang_format(cfc.CFunction_dict["SEOBNRv5_coprecessing_angles"].full_function)
    >>> _ = validate_strings(generated_str, "SEOBNRv5_coprecessing_angles__openmp", file_ext="c")
    >>> "notaknot_spline_second_derivatives" in generated_str
    True
    >>> "quat_sqrt_unit" in generated_str
    True
    >>> "quat_exp_z" in generated_str and "halfgamma" in generated_str
    True
    >>> "compute_py_style_coprecessing_rotator" in generated_str
    True
    >>> "rotation_max_dt_M" in generated_str
    False
    >>> "ceil(fabs(interval_dt) / rotation_max_dt_M)" in generated_str
    False
    >>> "q_total = quat_normalize(q_pre_arr[i]);" in generated_str
    True
    >>> "write_coprecessing_rotator_grid_variant_diagnostics" in generated_str
    False
    >>> "q_total = quat_multiply(q_pre, q_shift)" in generated_str
    False
    >>> "const quat_t q_shift =" in generated_str
    False
    >>> 'fopen("coprecessing_angles_debug.txt", "w")' in generated_str
    False
    >>> 'fopen("coprecessing_rotator_diagnostics.txt", "w")' in generated_str
    False
    >>> 'attach_inputs' in generated_str
    False
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    debug_file_init = (
        'fopen("coprecessing_angles_debug.txt", "w")'
        if enable_forensic_diagnostics
        else "NULL"
    )
    frame_invariant_file_init = (
        'fopen("frame_invariant_diagnostics.txt", "w")'
        if enable_forensic_diagnostics
        else "NULL"
    )
    spline_comparison_file_init = (
        'fopen("spline_comparison.txt", "w")' if enable_forensic_diagnostics else "NULL"
    )
    rotator_diagnostics_file_init = (
        'fopen("coprecessing_rotator_diagnostics.txt", "w")'
        if enable_forensic_diagnostics
        else "NULL"
    )
    forensic_setup = ""
    attach_input_forensics = ""
    rotator_forensics = ""
    loop_forensic_locals = ""
    q_pre_expression = "q_pre_arr[i]"
    frame_invariant_branch = ""
    forensic_close = ""
    if enable_forensic_diagnostics:
        attach_input_forensics = r"""
REAL omega_rISCO = interpolate_fine_omega_at_time(commondata, commondata->t_ISCO);
omega_rISCO = sanitize_omega_for_spin_splines(omega_rISCO, commondata->omega_spin_min, commondata->omega_spin_max, omega_boundary_tol);
check_omega_range("omega_rISCO", omega_rISCO, commondata->omega_spin_min, commondata->omega_spin_max, commondata->t_ISCO, 0);
if (fp_dbg != NULL) {
  REAL lnhat_rISCO_x = gsl_spline_eval(commondata->lnhat_x.spline, omega_rISCO, commondata->lnhat_x.acc);
  REAL lnhat_rISCO_y = gsl_spline_eval(commondata->lnhat_y.spline, omega_rISCO, commondata->lnhat_y.acc);
  REAL lnhat_rISCO_z = gsl_spline_eval(commondata->lnhat_z.spline, omega_rISCO, commondata->lnhat_z.acc);
  seobnr_normalize3(&lnhat_rISCO_x, &lnhat_rISCO_y, &lnhat_rISCO_z);
  const REAL chi1_rISCO_x = gsl_spline_eval(commondata->chi1_x_spline.spline, omega_rISCO, commondata->chi1_x_spline.acc);
  const REAL chi1_rISCO_y = gsl_spline_eval(commondata->chi1_y_spline.spline, omega_rISCO, commondata->chi1_y_spline.acc);
  const REAL chi1_rISCO_z = gsl_spline_eval(commondata->chi1_z_spline.spline, omega_rISCO, commondata->chi1_z_spline.acc);
  const REAL chi2_rISCO_x = gsl_spline_eval(commondata->chi2_x_spline.spline, omega_rISCO, commondata->chi2_x_spline.acc);
  const REAL chi2_rISCO_y = gsl_spline_eval(commondata->chi2_y_spline.spline, omega_rISCO, commondata->chi2_y_spline.acc);
  const REAL chi2_rISCO_z = gsl_spline_eval(commondata->chi2_z_spline.spline, omega_rISCO, commondata->chi2_z_spline.acc);
  const REAL chi1LN_rISCO = seobnr_vec_dot3(chi1_rISCO_x, chi1_rISCO_y, chi1_rISCO_z, lnhat_rISCO_x, lnhat_rISCO_y, lnhat_rISCO_z);
  const REAL chi2LN_rISCO = seobnr_vec_dot3(chi2_rISCO_x, chi2_rISCO_y, chi2_rISCO_z, lnhat_rISCO_x, lnhat_rISCO_y, lnhat_rISCO_z);
  const REAL ap_rISCO = commondata->m1 * chi1LN_rISCO + commondata->m2 * chi2LN_rISCO;
  const REAL am_rISCO = commondata->m1 * chi1LN_rISCO - commondata->m2 * chi2LN_rISCO;
  fprintf(fp_dbg, "attach_inputs %d %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n", -1, commondata->t_ISCO, commondata->Delta_t,
          commondata->t_attach, omega_rISCO, omega_attach, chi1LN_rISCO, chi2LN_rISCO, ap_rISCO, am_rISCO);
  fflush(fp_dbg);
}
"""
        rotator_forensics = r"""
FILE *fp_rotator = __ROTATOR_DIAGNOSTICS_FILE_INIT__;
if (fp_rotator != NULL) {
  fprintf(fp_rotator,
      "# label i time_M omega LN_J_x LN_J_y LN_J_z q_pre_w q_pre_x q_pre_y q_pre_z\n");
  dump_coprecessing_rotator_diagnostics_block(
      fp_rotator, "native_production_grid", n_insp, time_arr, omega_arr,
      LN_J_x_arr, LN_J_y_arr, LN_J_z_arr, q_pre_arr);
  write_coprecessing_rotator_grid_variant_diagnostics(
      fp_rotator, commondata, n_insp, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z,
      omega_boundary_tol);
  fclose(fp_rotator);
}
"""
        loop_forensic_locals = r"""
  const REAL time_M = time_arr[i];
  const REAL omega = omega_arr[i];
  const REAL LN_J_x = LN_J_x_arr[i];
  const REAL LN_J_y = LN_J_y_arr[i];
  const REAL LN_J_z = LN_J_z_arr[i];
  const quat_t q_pre = q_pre_arr[i];
"""
        q_pre_expression = "q_pre"
        forensic_setup = r"""
const REAL rotJ_matrix[3][3] = {
    {e1J_x, e2J_x, e3J_x},
    {e1J_y, e2J_y, e3J_y},
    {e1J_z, e2J_z, e3J_z}};
quat_t q_J2I = quat_from_rotation_matrix(rotJ_matrix);
const quat_t q_pre_ref = q_pre_arr[0];
if (quat_dot(q_pre_ref, q_J2I) < 0.0) {
  q_J2I = quat_negate(q_J2I);
}
const quat_t q_shift = quat_multiply(quat_conjugate(q_pre_ref), q_J2I);
if (fp_dbg != NULL) {
  fprintf(fp_dbg, "q_pre_ref %.15e %.15e %.15e %.15e\n", q_pre_ref.w, q_pre_ref.x, q_pre_ref.y, q_pre_ref.z);
  fprintf(fp_dbg, "q_J2I %.15e %.15e %.15e %.15e\n", q_J2I.w, q_J2I.x, q_J2I.y, q_J2I.z);
  fprintf(fp_dbg, "q_shift %.15e %.15e %.15e %.15e\n", q_shift.w, q_shift.x, q_shift.y, q_shift.z);
  fflush(fp_dbg);
}

FILE *fp_frame_inv = __FRAME_INVARIANT_FILE_INIT__;
REAL q_shift_z_x = 0.0;
REAL q_shift_z_y = 0.0;
REAL q_shift_z_z = 1.0;
quat_image_of_z(q_shift, &q_shift_z_x, &q_shift_z_y, &q_shift_z_z);
const REAL q_shift_z_angle = angle_between_unit_vectors(q_shift_z_x, q_shift_z_y, q_shift_z_z, 0.0, 0.0, 1.0);
if (fp_frame_inv != NULL) {
  fprintf(fp_frame_inv,
      "# q_shift_z_angle %.15e\n"
      "# i time_M omega LN_J_x LN_J_y LN_J_z"
      " q_pre_z_x q_pre_z_y q_pre_z_z q_total_z_x q_total_z_y q_total_z_z"
      " q_pre_z_angle q_total_z_angle\n",
      q_shift_z_angle);
}

// Development-only forensic dump; disabled in default generated C.
{
  const size_t n_spline_grid = 50000;
  FILE *fp_spline = __SPLINE_COMPARISON_FILE_INIT__;
  if (fp_spline != NULL) {
    fprintf(fp_spline,
        "# n_spline_grid %zu\n"
        "# omega lnhat_x lnhat_y lnhat_z"
        " chi1_x chi1_y chi1_z"
        " chi2_x chi2_y chi2_z"
        " L_x L_y L_z\n",
        n_spline_grid);
    const REAL omega_lo = commondata->omega_spin_min + omega_boundary_tol;
    const REAL omega_hi = commondata->omega_spin_max - omega_boundary_tol;
    const REAL log_lo   = log(omega_lo);
    const REAL log_hi   = log(omega_hi);
    for (size_t k = 0; k < n_spline_grid; k++) {
      const REAL frac    = (n_spline_grid > 1)
                           ? (REAL)k / (REAL)(n_spline_grid - 1) : 0.0;
      const REAL omega_k = exp(log_lo * (1.0 - frac) + log_hi * frac);
      const REAL ln_x = gsl_spline_eval(commondata->lnhat_x.spline, omega_k, commondata->lnhat_x.acc);
      const REAL ln_y = gsl_spline_eval(commondata->lnhat_y.spline, omega_k, commondata->lnhat_y.acc);
      const REAL ln_z = gsl_spline_eval(commondata->lnhat_z.spline, omega_k, commondata->lnhat_z.acc);
      const REAL c1x  = gsl_spline_eval(commondata->chi1_x_spline.spline, omega_k, commondata->chi1_x_spline.acc);
      const REAL c1y  = gsl_spline_eval(commondata->chi1_y_spline.spline, omega_k, commondata->chi1_y_spline.acc);
      const REAL c1z  = gsl_spline_eval(commondata->chi1_z_spline.spline, omega_k, commondata->chi1_z_spline.acc);
      const REAL c2x  = gsl_spline_eval(commondata->chi2_x_spline.spline, omega_k, commondata->chi2_x_spline.acc);
      const REAL c2y  = gsl_spline_eval(commondata->chi2_y_spline.spline, omega_k, commondata->chi2_y_spline.acc);
      const REAL c2z  = gsl_spline_eval(commondata->chi2_z_spline.spline, omega_k, commondata->chi2_z_spline.acc);
      const REAL Lx   = gsl_spline_eval(commondata->L_x.spline, omega_k, commondata->L_x.acc);
      const REAL Ly   = gsl_spline_eval(commondata->L_y.spline, omega_k, commondata->L_y.acc);
      const REAL Lz   = gsl_spline_eval(commondata->L_z.spline, omega_k, commondata->L_z.acc);
      fprintf(fp_spline,
          "%.15e %.15e %.15e %.15e"
          " %.15e %.15e %.15e"
          " %.15e %.15e %.15e"
          " %.15e %.15e %.15e\n",
          omega_k, ln_x, ln_y, ln_z,
          c1x, c1y, c1z,
          c2x, c2y, c2z,
          Lx, Ly, Lz);
    }
    fclose(fp_spline);
  }
}
"""
        frame_invariant_branch = r"""
  if (i == 0) {
    if (fp_frame_inv != NULL) {
      REAL q_pre_z_x = 0.0;
      REAL q_pre_z_y = 0.0;
      REAL q_pre_z_z = 1.0;
      REAL q_total_z_x = 0.0;
      REAL q_total_z_y = 0.0;
      REAL q_total_z_z = 1.0;
      quat_image_of_z(q_pre_ref, &q_pre_z_x, &q_pre_z_y, &q_pre_z_z);
      quat_image_of_z(q_total, &q_total_z_x, &q_total_z_y, &q_total_z_z);
      const REAL q_pre_z_angle = angle_between_unit_vectors(q_pre_z_x, q_pre_z_y, q_pre_z_z, LN_J_x, LN_J_y, LN_J_z);
      const REAL q_total_z_angle = angle_between_unit_vectors(q_total_z_x, q_total_z_y, q_total_z_z, LN_J_x, LN_J_y, LN_J_z);
      fprintf(fp_frame_inv,
          "%zu %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          (size_t)0, time_M, omega, LN_J_x, LN_J_y, LN_J_z,
          q_pre_z_x, q_pre_z_y, q_pre_z_z, q_total_z_x, q_total_z_y, q_total_z_z,
          q_pre_z_angle, q_total_z_angle);
    } // END IF: frame-invariant diagnostics file is open for first sample
  } else {
    if (fp_frame_inv != NULL) {
      REAL q_pre_z_x = 0.0;
      REAL q_pre_z_y = 0.0;
      REAL q_pre_z_z = 1.0;
      REAL q_total_z_x = 0.0;
      REAL q_total_z_y = 0.0;
      REAL q_total_z_z = 1.0;
      quat_image_of_z(q_pre, &q_pre_z_x, &q_pre_z_y, &q_pre_z_z);
      quat_image_of_z(q_total, &q_total_z_x, &q_total_z_y, &q_total_z_z);
      const REAL q_pre_z_angle = angle_between_unit_vectors(q_pre_z_x, q_pre_z_y, q_pre_z_z, LN_J_x, LN_J_y, LN_J_z);
      const REAL q_total_z_angle = angle_between_unit_vectors(q_total_z_x, q_total_z_y, q_total_z_z, LN_J_x, LN_J_y, LN_J_z);
      fprintf(fp_frame_inv,
          "%zu %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
          i, time_M, omega, LN_J_x, LN_J_y, LN_J_z,
          q_pre_z_x, q_pre_z_y, q_pre_z_z, q_total_z_x, q_total_z_y, q_total_z_z,
          q_pre_z_angle, q_total_z_angle);
    } // END IF: frame-invariant diagnostics file is open for later sample
  } // END ELSE: forensic frame-invariant sample is not first inspiral sample
"""
        forensic_close = r"""
if (fp_frame_inv != NULL) {
  fclose(fp_frame_inv);
}
"""
    prefunc = r"""
#include <gsl/gsl_spline.h>

typedef struct {
  REAL w;
  REAL x;
  REAL y;
  REAL z;
} quat_t;

static inline REAL seobnr_clamp_real(const REAL x, const REAL xmin, const REAL xmax) {
  return x < xmin ? xmin : (x > xmax ? xmax : x);
}

static inline REAL seobnr_vec_dot3(const REAL ax, const REAL ay, const REAL az,
                                   const REAL bx, const REAL by, const REAL bz) {
  return ax * bx + ay * by + az * bz;
}

static inline REAL seobnr_vec_norm3(const REAL x, const REAL y, const REAL z) {
  return sqrt(x * x + y * y + z * z);
}

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

static inline quat_t quat_conjugate(const quat_t q) {
  return (quat_t){q.w, -q.x, -q.y, -q.z};
}

static inline quat_t quat_negate(const quat_t q) {
  return (quat_t){-q.w, -q.x, -q.y, -q.z};
}

static inline quat_t quat_scale(const quat_t q, const REAL scale) {
  return (quat_t){q.w * scale, q.x * scale, q.y * scale, q.z * scale};
}

static inline REAL quat_dot(const quat_t a, const quat_t b) {
  return a.w * b.w + a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline quat_t quat_multiply(const quat_t a, const quat_t b) {
  return (quat_t){
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
      a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
      a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
      a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w};
}

static inline quat_t quat_exp_z(const REAL half_angle) {
  return (quat_t){cos(half_angle), 0.0, 0.0, sin(half_angle)};
}

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

static inline void quat_image_of_z(
    const quat_t q_in,
    REAL *restrict zx,
    REAL *restrict zy,
    REAL *restrict zz) {
  REAL rot_mat[3][3];
  quat_to_rotation_matrix(q_in, rot_mat);
  *zx = rot_mat[0][2];
  *zy = rot_mat[1][2];
  *zz = rot_mat[2][2];
  seobnr_normalize3(zx, zy, zz);
}

static inline REAL angle_between_unit_vectors(
    const REAL ax, const REAL ay, const REAL az,
    const REAL bx, const REAL by, const REAL bz) {
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
static inline void extract_zyz_euler_angles_from_quat(
    const quat_t q, REAL *restrict alpha, REAL *restrict beta, REAL *restrict gamma) {
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

static inline void sample_inspiral_time_and_omega(
    const commondata_struct *restrict commondata,
    const size_t i,
    REAL *restrict time_M,
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
static inline REAL interpolate_fine_omega_at_time(
    const commondata_struct *restrict commondata, const REAL target_time) {
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
static void notaknot_spline_second_derivatives(
    const size_t n,
    const REAL *restrict x,
    const REAL *restrict y,
    REAL *restrict second_deriv) {
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
static void notaknot_spline_derivatives_at_knots(
    const size_t n,
    const REAL *restrict x,
    const REAL *restrict y,
    const REAL *restrict second_deriv,
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
static void notaknot_spline_integral_prefix(
    const size_t n,
    const REAL *restrict x,
    const REAL *restrict y,
    const REAL *restrict second_deriv,
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
 * Linearly interpolate a monotonic tabulated series.
 *
 * @param[in] n Number of samples.
 * @param[in] x Monotonic sample coordinates.
 * @param[in] y Sample values.
 * @param[in] target Coordinate at which to interpolate.
 * @return Interpolated value.
 */
static REAL interpolate_monotonic_series(
    const size_t n,
    const REAL *restrict x,
    const REAL *restrict y,
    const REAL target) {
  if (n == 0) {
    return 0.0;
  }
  if (n == 1 || target <= x[0]) {
    return y[0];
  }
  for (size_t i = 0; i + 1 < n; i++) {
    if (target <= x[i + 1] || i + 2 == n) {
      const REAL denom = x[i + 1] - x[i];
      if (fabs(denom) < 1e-15) {
        return y[i];
      }
      const REAL frac = (target - x[i]) / denom;
      return y[i] + frac * (y[i + 1] - y[i]);
    }
  }
  return y[n - 1];
} // END FUNCTION: interpolate_monotonic_series

/**
 * Evaluate the orbital angular momentum direction in the J-frame.
 *
 * @param[in] commondata Common data struct containing spin-dynamics splines.
 * @param[in] omega Orbital frequency sample.
 * @param[in] e1J_x J-frame e1 x component.
 * @param[in] e1J_y J-frame e1 y component.
 * @param[in] e1J_z J-frame e1 z component.
 * @param[in] e2J_x J-frame e2 x component.
 * @param[in] e2J_y J-frame e2 y component.
 * @param[in] e2J_z J-frame e2 z component.
 * @param[in] e3J_x J-frame e3 x component.
 * @param[in] e3J_y J-frame e3 y component.
 * @param[in] e3J_z J-frame e3 z component.
 * @param[out] LN_J_x J-frame LN x component.
 * @param[out] LN_J_y J-frame LN y component.
 * @param[out] LN_J_z J-frame LN z component.
 */
static void evaluate_LN_J_from_omega(
    const commondata_struct *restrict commondata,
    const REAL omega,
    const REAL e1J_x, const REAL e1J_y, const REAL e1J_z,
    const REAL e2J_x, const REAL e2J_y, const REAL e2J_z,
    const REAL e3J_x, const REAL e3J_y, const REAL e3J_z,
    REAL *restrict LN_J_x,
    REAL *restrict LN_J_y,
    REAL *restrict LN_J_z) {
  REAL lnhat_I_x = gsl_spline_eval(commondata->lnhat_x.spline, omega, commondata->lnhat_x.acc);
  REAL lnhat_I_y = gsl_spline_eval(commondata->lnhat_y.spline, omega, commondata->lnhat_y.acc);
  REAL lnhat_I_z = gsl_spline_eval(commondata->lnhat_z.spline, omega, commondata->lnhat_z.acc);
  seobnr_normalize3(&lnhat_I_x, &lnhat_I_y, &lnhat_I_z);
  *LN_J_x = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e1J_x, e1J_y, e1J_z);
  *LN_J_y = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e2J_x, e2J_y, e2J_z);
  *LN_J_z = seobnr_vec_dot3(lnhat_I_x, lnhat_I_y, lnhat_I_z, e3J_x, e3J_y, e3J_z);
  seobnr_normalize3(LN_J_x, LN_J_y, LN_J_z);
} // END FUNCTION: evaluate_LN_J_from_omega

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
static void compute_py_style_coprecessing_rotator(
    const size_t n,
    const REAL *restrict time_arr,
    const REAL *restrict omega_arr,
    const REAL *restrict LN_J_x_arr,
    const REAL *restrict LN_J_y_arr,
    const REAL *restrict LN_J_z_arr,
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
  if (q_domega == NULL || q_component == NULL || spline_second == NULL || spline_deriv == NULL ||
      omega_second == NULL || domega_dt == NULL || halfgammadot == NULL || halfgamma == NULL || halfgamma_second == NULL) {
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

/**
 * Write a rotator diagnostic block for one sampled grid.
 *
 * @param[in,out] fp Open diagnostic file pointer.
 * @param[in] label Diagnostic block label.
 * @param[in] n Number of samples.
 * @param[in] time_arr Time samples.
 * @param[in] omega_arr Orbital-frequency samples.
 * @param[in] LN_J_x_arr J-frame LN x samples.
 * @param[in] LN_J_y_arr J-frame LN y samples.
 * @param[in] LN_J_z_arr J-frame LN z samples.
 * @param[in] q_pre_arr Minimal-rotation quaternions.
 */
static void dump_coprecessing_rotator_diagnostics_block(
    FILE *restrict fp,
    const char *restrict label,
    const size_t n,
    const REAL *restrict time_arr,
    const REAL *restrict omega_arr,
    const REAL *restrict LN_J_x_arr,
    const REAL *restrict LN_J_y_arr,
    const REAL *restrict LN_J_z_arr,
    const quat_t *restrict q_pre_arr) {
  if (fp == NULL) {
    return;
  }
  fprintf(fp, "# block %s n %zu\n", label, n);
  for (size_t i = 0; i < n; i++) {
    fprintf(fp,
        "%s %zu %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e %.15e\n",
        label, i, time_arr[i], omega_arr[i], LN_J_x_arr[i], LN_J_y_arr[i], LN_J_z_arr[i],
        q_pre_arr[i].w, q_pre_arr[i].x, q_pre_arr[i].y, q_pre_arr[i].z);
  }
  fflush(fp);
} // END FUNCTION: dump_coprecessing_rotator_diagnostics_block

static inline REAL sanitize_omega_for_spin_splines(
    const REAL omega,
    const REAL omega_min,
    const REAL omega_max,
    const REAL abs_tol);

/**
 * Build and write one forensic rotator grid-variant diagnostic.
 *
 * @param[in,out] fp Open diagnostic file pointer.
 * @param[in] label Diagnostic block label.
 * @param[in] commondata Common data struct containing splines.
 * @param[in] n Number of native samples.
 * @param[in] time_arr Native time samples.
 * @param[in] omega_arr Native orbital-frequency samples.
 * @param[in] e1J_x J-frame e1 x component.
 * @param[in] e1J_y J-frame e1 y component.
 * @param[in] e1J_z J-frame e1 z component.
 * @param[in] e2J_x J-frame e2 x component.
 * @param[in] e2J_y J-frame e2 y component.
 * @param[in] e2J_z J-frame e2 z component.
 * @param[in] e3J_x J-frame e3 x component.
 * @param[in] e3J_y J-frame e3 y component.
 * @param[in] e3J_z J-frame e3 z component.
 * @param[in] stride Native-grid stride.
 * @param[in] midpoint_refined Whether to insert midpoint samples.
 * @param[in] uniform_omega Whether to use a uniform omega grid.
 * @param[in] uniform_factor Uniform-grid multiplier.
 * @param[in] omega_boundary_tol Frequency boundary tolerance.
 */
static void write_one_coprecessing_rotator_grid_variant(
    FILE *restrict fp,
    const char *restrict label,
    const commondata_struct *restrict commondata,
    const size_t n,
    const REAL *restrict time_arr,
    const REAL *restrict omega_arr,
    const REAL e1J_x, const REAL e1J_y, const REAL e1J_z,
    const REAL e2J_x, const REAL e2J_y, const REAL e2J_z,
    const REAL e3J_x, const REAL e3J_y, const REAL e3J_z,
    const size_t stride,
    const int midpoint_refined,
    const int uniform_omega,
    const size_t uniform_factor,
    const REAL omega_boundary_tol) {
  size_t n_grid = 0;
  if (uniform_omega) {
    n_grid = uniform_factor * n;
  } else if (midpoint_refined) {
    n_grid = n > 0 ? 2 * n - 1 : 0;
  } else {
    n_grid = (n + stride - 1) / stride;
    if (n_grid > 0 && (n - 1) % stride != 0)
      n_grid++;
  } // END ELSE: use strided native grid
  if (n_grid < 4)
    return;

  REAL *restrict t_grid = (REAL *)malloc(n_grid * sizeof(REAL));
  REAL *restrict omega_grid = (REAL *)malloc(n_grid * sizeof(REAL));
  REAL *restrict LN_x = (REAL *)malloc(n_grid * sizeof(REAL));
  REAL *restrict LN_y = (REAL *)malloc(n_grid * sizeof(REAL));
  REAL *restrict LN_z = (REAL *)malloc(n_grid * sizeof(REAL));
  quat_t *restrict q_grid = (quat_t *)malloc(n_grid * sizeof(quat_t));
  if (t_grid == NULL || omega_grid == NULL || LN_x == NULL || LN_y == NULL || LN_z == NULL || q_grid == NULL) {
    fprintf(stderr, "Error: write_one_coprecessing_rotator_grid_variant() malloc failed\n");
    exit(1);
  } // END IF: rotator grid-variant allocation failed

  if (uniform_omega) {
    const REAL omega_min = omega_arr[0];
    const REAL omega_max = omega_arr[n - 1];
    for (size_t i = 0; i < n_grid; i++) {
      const REAL frac = n_grid > 1 ? (REAL)i / (REAL)(n_grid - 1) : 0.0;
      omega_grid[i] = omega_min + frac * (omega_max - omega_min);
      t_grid[i] = interpolate_monotonic_series(n, omega_arr, time_arr, omega_grid[i]);
    } // END LOOP: for i over uniform omega diagnostic samples
  } else if (midpoint_refined) {
    for (size_t i = 0; i < n; i++) {
      const size_t out = 2 * i;
      t_grid[out] = time_arr[i];
      omega_grid[out] = omega_arr[i];
      if (i + 1 < n) {
        t_grid[out + 1] = 0.5 * (time_arr[i] + time_arr[i + 1]);
        omega_grid[out + 1] = 0.5 * (omega_arr[i] + omega_arr[i + 1]);
      } // END IF: midpoint sample exists after native sample
    } // END LOOP: for i over midpoint-refined diagnostic samples
  } else {
    size_t out = 0;
    for (size_t i = 0; i < n; i += stride) {
      t_grid[out] = time_arr[i];
      omega_grid[out] = omega_arr[i];
      out++;
    } // END LOOP: for i over strided native diagnostic samples
    if (out == 0 || omega_grid[out - 1] != omega_arr[n - 1]) {
      t_grid[out] = time_arr[n - 1];
      omega_grid[out] = omega_arr[n - 1];
      out++;
    } // END IF: append final native diagnostic sample
    n_grid = out;
  } // END ELSE: use strided native diagnostic grid

  for (size_t i = 0; i < n_grid; i++) {
    omega_grid[i] = sanitize_omega_for_spin_splines(omega_grid[i], commondata->omega_spin_min, commondata->omega_spin_max, omega_boundary_tol);
    evaluate_LN_J_from_omega(commondata, omega_grid[i],
        e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z,
        &LN_x[i], &LN_y[i], &LN_z[i]);
  } // END LOOP: for i over diagnostic grid LN_J samples
  compute_py_style_coprecessing_rotator(n_grid, t_grid, omega_grid, LN_x, LN_y, LN_z, q_grid);
  dump_coprecessing_rotator_diagnostics_block(fp, label, n_grid, t_grid, omega_grid, LN_x, LN_y, LN_z, q_grid);

  free(t_grid);
  free(omega_grid);
  free(LN_x);
  free(LN_y);
  free(LN_z);
  free(q_grid);
} // END FUNCTION: write_one_coprecessing_rotator_grid_variant

/**
 * Write the suite of forensic rotator grid-variant diagnostics.
 *
 * @param[in,out] fp Open diagnostic file pointer.
 * @param[in] commondata Common data struct containing splines.
 * @param[in] n Number of native samples.
 * @param[in] time_arr Native time samples.
 * @param[in] omega_arr Native orbital-frequency samples.
 * @param[in] e1J_x J-frame e1 x component.
 * @param[in] e1J_y J-frame e1 y component.
 * @param[in] e1J_z J-frame e1 z component.
 * @param[in] e2J_x J-frame e2 x component.
 * @param[in] e2J_y J-frame e2 y component.
 * @param[in] e2J_z J-frame e2 z component.
 * @param[in] e3J_x J-frame e3 x component.
 * @param[in] e3J_y J-frame e3 y component.
 * @param[in] e3J_z J-frame e3 z component.
 * @param[in] omega_boundary_tol Frequency boundary tolerance.
 */
static void write_coprecessing_rotator_grid_variant_diagnostics(
    FILE *restrict fp,
    const commondata_struct *restrict commondata,
    const size_t n,
    const REAL *restrict time_arr,
    const REAL *restrict omega_arr,
    const REAL e1J_x, const REAL e1J_y, const REAL e1J_z,
    const REAL e2J_x, const REAL e2J_y, const REAL e2J_z,
    const REAL e3J_x, const REAL e3J_y, const REAL e3J_z,
    const REAL omega_boundary_tol) {
  write_one_coprecessing_rotator_grid_variant(fp, "native_drop_every_2", commondata, n, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z, 2, 0, 0, 1, omega_boundary_tol);
  write_one_coprecessing_rotator_grid_variant(fp, "native_drop_every_4", commondata, n, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z, 4, 0, 0, 1, omega_boundary_tol);
  write_one_coprecessing_rotator_grid_variant(fp, "native_midpoint_refined_2x", commondata, n, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z, 1, 1, 0, 1, omega_boundary_tol);
  write_one_coprecessing_rotator_grid_variant(fp, "uniform_omega_native_n", commondata, n, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z, 1, 0, 1, 1, omega_boundary_tol);
  write_one_coprecessing_rotator_grid_variant(fp, "uniform_omega_4x_native_n", commondata, n, time_arr, omega_arr,
      e1J_x, e1J_y, e1J_z, e2J_x, e2J_y, e2J_z, e3J_x, e3J_y, e3J_z, 1, 0, 1, 4, omega_boundary_tol);
} // END FUNCTION: write_coprecessing_rotator_grid_variant_diagnostics

static inline REAL sanitize_omega_for_spin_splines(
    const REAL omega,
    const REAL omega_min,
    const REAL omega_max,
    const REAL abs_tol) {
  if (omega < omega_min && fabs(omega - omega_min) <= abs_tol)
    return omega_min;
  if (omega > omega_max && fabs(omega - omega_max) <= abs_tol)
    return omega_max;
  return omega;
} // END FUNCTION: sanitize_omega_for_spin_splines

static inline void check_omega_range(
    const char *restrict label,
    const REAL omega_val,
    const REAL omega_min,
    const REAL omega_max,
    const REAL time_val,
    const size_t idx_val) {
  if (omega_val < omega_min || omega_val > omega_max) {
    fprintf(stderr,
            "SEOBNRv5_coprecessing_angles: %s out of range at i=%zu, t=%.15e: %.15e not in [%.15e, %.15e]\n",
            label, idx_val, (double)time_val, (double)omega_val,
            (double)omega_min, (double)omega_max);
    exit(1);
  } // END IF: omega sample lies outside spin-spline range
} // END FUNCTION: check_omega_range
"""
    if not enable_forensic_diagnostics:
        helper_start = prefunc.index(
            "/**\n * Linearly interpolate a monotonic tabulated series."
        )
        helper_end = prefunc.index(
            "/**\n * Compute the pyseobnr-style minimal-rotation quaternion track."
        )
        prefunc = prefunc[:helper_start] + prefunc[helper_end:]

        helper_start = prefunc.index("/**\n * Write a rotator diagnostic block")
        helper_end = prefunc.index(
            "static inline REAL sanitize_omega_for_spin_splines",
            prefunc.index(
                "} // END FUNCTION: write_coprecessing_rotator_grid_variant_diagnostics"
            ),
        )
        prefunc = prefunc[:helper_start] + prefunc[helper_end:]

    desc = r"""
Compute the inspiral J->P Euler-angle arrays for SEOBNRv5 coprecessing rotations.

This routine follows the inspiral frame construction of Ramos-Buades et al.
(`SEOBNRv5PHM`, arXiv:2303.18046), especially Sec. III A, Eq. (15), and
Eqs. (16)-(17). The minimal-rotation condition is imposed using a quaternion
spline correction following Boyle, Owen, and Pfeiffer, Phys. Rev. D 84, 124011
(2011), arXiv:1110.2965, Eq. (18). The robust J-frame basis construction is
consistent with the earlier SEOBNR implementations of Babak, Taracchini,
and Buonanno, Phys. Rev. D 95, 024010 (2017), arXiv:1607.05661, and
Ossokine et al., Phys. Rev. D 102, 044055 (2020), arXiv:2004.09442.

@param[in,out] commondata Common data struct updated with output Euler-angle arrays.
"""
    cfunc_type = "void"
    name = "SEOBNRv5_coprecessing_angles"
    params = "commondata_struct *restrict commondata"
    body = r"""
const size_t n_low = commondata->nsteps_low;
const size_t n_fine = commondata->nsteps_fine;
const size_t n_insp = n_low + n_fine;
FILE *fp_dbg = __DEBUG_FILE_INIT__;
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
  fprintf(fp_dbg, "attach %d %.15e %.15e %.15e %.15e\n", -1, commondata->t_attach, omega_attach, commondata->omega_spin_min, commondata->omega_spin_max);
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
__ATTACH_INPUT_FORENSICS__
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

__ROTATOR_FORENSICS__

__FORENSIC_SETUP__

for (size_t i = 0; i < n_insp; i++) {
__LOOP_FORENSIC_LOCALS__
  const quat_t q_total = quat_normalize(__Q_PRE_EXPRESSION__);
  extract_zyz_euler_angles_from_quat(q_total, &commondata->alpha_JP[i], &commondata->beta_JP[i], &commondata->gamma_JP[i]);
__FRAME_INVARIANT_BRANCH__
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
__FORENSIC_CLOSE__
"""
    body = body.replace("__DEBUG_FILE_INIT__", debug_file_init)
    body = body.replace("__ATTACH_INPUT_FORENSICS__", attach_input_forensics)
    body = body.replace("__ROTATOR_FORENSICS__", rotator_forensics)
    body = body.replace("__FORENSIC_SETUP__", forensic_setup)
    body = body.replace("__LOOP_FORENSIC_LOCALS__", loop_forensic_locals)
    body = body.replace("__Q_PRE_EXPRESSION__", q_pre_expression)
    body = body.replace("__FRAME_INVARIANT_BRANCH__", frame_invariant_branch)
    body = body.replace("__FORENSIC_CLOSE__", forensic_close)
    body = body.replace("__FRAME_INVARIANT_FILE_INIT__", frame_invariant_file_init)
    body = body.replace("__SPLINE_COMPARISON_FILE_INIT__", spline_comparison_file_init)
    body = body.replace(
        "__ROTATOR_DIAGNOSTICS_FILE_INIT__", rotator_diagnostics_file_init
    )
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

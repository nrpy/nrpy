"""
Register the generalized proper-circumference spin diagnostic for BHaHAHA.

The generated routine accepts the unit spin direction selected by the SpECTRE
surface-integral diagnostic. It constructs an equatorial great circle
orthogonal to that direction and one polar great circle containing it. On the
cell-centered angular grid it precomputes scalar line-element integrands from
the induced horizon metric and smooth analytic tangent fields, ghost-fills
those scalars, and interpolates them to midpoint quadrature samples. No angle
unwrapping or finite differencing of spherical angles is performed.

The polar meridian uses one deterministic transverse basis vector. All such
meridians agree for an axisymmetric Kerr horizon, but distorted horizons and
finite-resolution data can retain some dependence on that basis choice. The
Kerr circumference relation is inverted with a bracketed Newton/bisection
solve on the physical dimensionless-spin interval [0, 1].

Authors: Wesley Inselman
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.general_relativity.bhahaha.area as bhahaha_area
import nrpy.helpers.parallel_codegen as pcg


def register_CFunction_diagnostics_proper_circumferences_general(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C routine that computes equatorial/polar proper circumferences and a spin estimate from their ratio.

    The generated function precomputes scalar equatorial and polar
    line-element integrands from the induced metric, interpolates them along
    great circles defined by the supplied SpECTRE spin axis, integrates the
    proper circumferences, and inverts their Kerr ratio for a dimensionless
    spin magnitude and vector.

    :param enable_fd_functions: Whether to emit FD derivative helpers in the generated code (passed through to codegen).
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
/**
 * Complete elliptic integrals E(m) and K(m) via 8th-order midpoint quadrature.
 *
 * Computes the complete elliptic integrals in the **parameter** (Legendre) form
 *   E(m) = ∫₀^{π/2} sqrt(1 − m sin²θ) dθ,
 *   K(m) = ∫₀^{π/2} dθ / sqrt(1 − m sin²θ).
 *
 * The implementation uses a fixed 8th-order midpoint rule with periodic weights
 * over [0, π/2]. The sample count is fixed at 128 (power of two) for high accuracy
 * and good vectorization.
 *
 * @param m  Legendre **parameter** (a.k.a. k²). Real results require m ≤ 1.
 *           Negative m is supported (e.g. m ∈ [−1, 0] in this code path).
 *           @warning This is the parameter m, **not** the modulus k. If you have a
 *           modulus k, pass m = k*k.
 * @param[out] E  On return, E(m).
 * @param[out] K  On return, K(m) (diverges as m → 1⁻).
 *
 * @pre E and K are non-null.
 * @pre The internal weight generator expects the sample count to be compatible
 *      with the 8th-order periodic stencil (here fixed to 128).
 */
static void elliptic_E_and_K_integrals(const REAL m, REAL *restrict E, REAL *restrict K) {
  static const int N_sample_pts = 128; // Number of sample points for integration. Chosen for high precision.
  const REAL *restrict weights;        // Precomputed integration weights for accuracy in the midpoint method.
  int weight_stencil_size;             // Size of the weight stencil to correctly cycle through the weights.

  // Retrieve the integration weights based on the number of sample points (divisible by 8 -> 8th order).
  bah_diagnostics_integration_weights(N_sample_pts, N_sample_pts, &weights, &weight_stencil_size);

  const REAL a = 0.0;                            // Lower limit of integration (0 radians).
  const REAL b = M_PI / 2.0;                     // Upper limit of integration (pi/2 radians).
  const REAL h = (b - a) / ((REAL)N_sample_pts); // Step size for each subinterval based on sample points.

  REAL sum_E = 0.0; // Accumulator for the elliptic integral of the second kind E(m).
  REAL sum_K = 0.0; // Accumulator for the elliptic integral of the first kind K(m).

  // Parallelized loop to compute both integrals E(m) and K(m) using OpenMP.
  // #pragma omp parallel for reduction(+ : sum_E, sum_K) <- thread creation/destruction likely -> slower code here
  for (int i = 0; i < N_sample_pts; i++) {
    const REAL theta = a + ((REAL)i + 0.5) * h; // Compute the midpoint for the current subinterval.
    // Compute sin(theta). For optimization, we could use a lookup table if N_sample_pts is constant.
    const REAL sintheta = sin(theta);
    // Compute the integrands for the elliptic integrals of the second & first kinds (E(m) & K(m), respectively) at this midpoint.
    const REAL elliptic_E_integrand = sqrt(1.0 - m * sintheta * sintheta);
    const REAL elliptic_K_integrand = 1.0 / elliptic_E_integrand;
    // Update the running sums for E(m) & K(m), applying the corresponding weight.
    sum_E += weights[i % weight_stencil_size] * elliptic_E_integrand;
    sum_K += weights[i % weight_stencil_size] * elliptic_K_integrand;
  } // END LOOP: for i over sample points to compute both integrals

  // Multiply by the step size to complete the integration and store the results.
  *E = sum_E * h; // Elliptic integral of the second kind.
  *K = sum_K * h; // Elliptic integral of the first kind.
} // END FUNCTION: elliptic_E_and_K_integrals

/**
 * Estimates the spin parameter magnitude for equilibrium black holes based on the circumference ratio C_r.
 *
 * The physical Kerr range is C_r in [0.6080033617124899, 1]. Small endpoint
 * excursions are tolerance-clamped. Interior roots use an analytic initial
 * guess followed by Newton steps safeguarded by a shrinking bracket; unsafe
 * steps fall back to bisection.
 *
 * @param C_r_input The circumference ratio used to estimate the spin.
 * @return The estimated spin, or BHAHAHA_DIAGNOSTIC_UNAVAILABLE for invalid
 *         input or failure to converge.
 *
 */
static REAL compute_spin(const REAL C_r_input) {
  const REAL min_kerr_ratio = 0.6080033617124899;
  const REAL endpoint_tolerance = 1.0e-10;
  if (!isfinite(C_r_input) || C_r_input <= 0.0 ||
      C_r_input > 1.0 + endpoint_tolerance || C_r_input < min_kerr_ratio - endpoint_tolerance)
    return BHAHAHA_DIAGNOSTIC_UNAVAILABLE;

  const REAL C_r = fmin(1.0, fmax(min_kerr_ratio, C_r_input));
  if (1.0 - C_r <= endpoint_tolerance)
    return 0.0;
  if (C_r - min_kerr_ratio <= endpoint_tolerance)
    return 1.0;

  REAL lower = 0.0;
  REAL upper = 1.0;
  REAL spin = 0.5;
  const REAL spin_sq = 1.0 - (2.55 * C_r - 1.55) * (2.55 * C_r - 1.55);
  if (spin_sq > 0.0 && spin_sq < 1.0)
    spin = sqrt(spin_sq);
  spin = fmin(1.0 - 1.0e-8, fmax(1.0e-8, spin));

  const REAL abs_tolerance = 1.0e-12;
  const REAL rel_tolerance = 1.0e-10;
  const REAL residual_tolerance = 1.0e-13;
  const int max_its = 100;
  for (int it = 0; it < max_its; it++) {
    const REAL x = spin;
    const REAL one_minus_x_squared = 1.0 - x * x;
    if (!(one_minus_x_squared > 0.0) || !isfinite(one_minus_x_squared))
      return BHAHAHA_DIAGNOSTIC_UNAVAILABLE;

    const REAL root = sqrt(one_minus_x_squared);
    const REAL elliptic_parameter = -(x * x) / ((1.0 + root) * (1.0 + root));
    REAL E, K;
    elliptic_E_and_K_integrals(elliptic_parameter, &E, &K);
    if (!isfinite(E) || !isfinite(K))
      return BHAHAHA_DIAGNOSTIC_UNAVAILABLE;

    const REAL residual = (1.0 + root) * E / M_PI - C_r;
    if (!isfinite(residual))
      return BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
    if (fabs(residual) <= residual_tolerance)
      return x;
    if (residual > 0.0)
      lower = x;
    else
      upper = x;

    const REAL bracket_width = upper - lower;
    if (bracket_width <= abs_tolerance + rel_tolerance * fmax(fabs(lower), fabs(upper)))
      return 0.5 * (lower + upper);

    // Complete a Newton-Raphson iteration, then safeguard it with the known bracket.
"""
    prefunc += ccg.c_codegen(
        bhahaha_area.spin_NewtonRaphson(),
        "const REAL x_np1",
        include_braces=False,
    )
    prefunc += r"""
    const REAL safeguarded_lower = lower + 0.25 * bracket_width;
    const REAL safeguarded_upper = upper - 0.25 * bracket_width;
    if (isfinite(x_np1) && x_np1 > safeguarded_lower && x_np1 < safeguarded_upper)
      spin = x_np1;
    else
      spin = 0.5 * (lower + upper);
  } // END LOOP: bracketed Newton/bisection inversion

  return BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
} // END FUNCTION: compute_spin

// Apply inner BCs for a selection of gridfunctions
// Note: Nxx_plus_2NGHOSTS2 is needed for IDX4pt()
static void apply_inner_bc_for_selected_gfs(bc_struct *restrict bcstruct, REAL *restrict metric_data_gfs, const int Nxx_plus_2NGHOSTS0,
                                            const int Nxx_plus_2NGHOSTS1, const int Nxx_plus_2NGHOSTS2, const int *which_gfs, const int num_gfs) {

  const bc_info_struct *bc_info = &bcstruct->bc_info;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2

  // Apply boundary conditions at inner boundary points for the selected gridfunctions.
#pragma omp parallel for collapse(2)
  for (int gf_idx = 0; gf_idx < num_gfs; gf_idx++) {
    for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
      const int which_gf = which_gfs[gf_idx];
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt; // Destination point index.
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt; // Source point index for copying.

      // Extract the i0, i1, and i2 indices from dstpt and srcpt.
      //  -> idx3 = i + Nx0*(j + Nx1*k)
      //  -> i = mod(idx3, Nx0)
      //   tmp = (idx3-i)/Nx0
      //  -> j = mod(tmp, Nx1)
      //  -> k = (tmp-j)/Nx1
      const int dst_i0 = dstpt % Nxx_plus_2NGHOSTS0;
      const int dsttmp = (dstpt - dst_i0) / Nxx_plus_2NGHOSTS0;
      const int dst_i1 = dsttmp % Nxx_plus_2NGHOSTS1;
      const int dst_i2 = (dsttmp - dst_i1) / Nxx_plus_2NGHOSTS1;

      const int src_i0 = srcpt % Nxx_plus_2NGHOSTS0;
      const int srctmp = (srcpt - src_i0) / Nxx_plus_2NGHOSTS0;
      const int src_i1 = srctmp % Nxx_plus_2NGHOSTS1;
      const int src_i2 = (srctmp - src_i1) / Nxx_plus_2NGHOSTS1;

      if (dst_i0 == NGHOSTS) {
        metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(dst_i1, dst_i2)] = metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(src_i1, src_i2)];
      }
    } // END LOOP: for pt over inner boundary points
  } // END LOOP: for which_gf over gridfunctions
} // END FUNCTION: apply_inner_bc_for_selected_gfs

// ---------- small vector & math helpers (pure C) ----------
// These are intentionally tiny & inlined: they live in tight OpenMP loops,
// and we want predictable codegen and good autovectorization.
static inline REAL dot3(const REAL a[3], const REAL b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline void cross3(const REAL a[3], const REAL b[3], REAL c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
} // END FUNCTION: cross3
static inline REAL norm3(const REAL a[3]) { return sqrt(dot3(a,a)); }
static inline void normalize3(REAL a[3]) { const REAL n = norm3(a); if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; } }
// Build an orthonormal basis {e1, e2} orthogonal to s; avoid near-collinearity for numerical stability.
static inline void build_basis_from_s(const REAL s[3], REAL e1[3], REAL e2[3]) {
  REAL a[3] = {1.0, 0.0, 0.0};
  if (fabs(dot3(a,s)) > 0.9) { a[0]=0.0; a[1]=1.0; a[2]=0.0; }
  cross3(s, a, e1); normalize3(e1);
  cross3(s, e1, e2); normalize3(e2);
} // END FUNCTION: build_basis_from_s
// Cartesian unit vector -> spherical angles (theta, phi); inputs here are unit by construction.
static inline void cart_to_sph(const REAL r[3], REAL *theta, REAL *phi) {
  const REAL x=r[0], y=r[1], z=r[2];
  REAL zz = z; if (zz < -1.0) zz = -1.0; if (zz > 1.0) zz = 1.0;
  *theta = acos(zz);
  *phi   = atan2(y, x);
} // END FUNCTION: cart_to_sph

// Midpoint integrate over alpha with precomputed 8th-order weights (N_angle == Nxx2 is typically divisible by 8 -> 8th order).
static inline REAL integrate_over_alpha(const REAL *vals, int N_angle, REAL d_alpha) {
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);
  REAL sum = 0.0;
  // #pragma omp parallel for reduction(+ : sum) // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++)
    sum += vals[i] * weights[i % weight_stencil_size];
  return sum * d_alpha;
} // END FUNCTION: integrate_over_alpha
/* ---------- end helpers ---------- */
"""
    desc = """
Compute proper circumferences and Kerr spin magnitude about the supplied
SpECTRE spin direction using the induced metric on r = h(theta, phi).

Line element integrated (alpha parametrizes the great circle):
    ds = sqrt( q_tt * (dtheta/dalpha)^2 + 2*q_tp*(dtheta/dalpha)*(dphi/dalpha) + q_pp * (dphi/dalpha)^2 ) dalpha,
with q_tt = (sqrt_qtt)^2, q_pp = (sqrt_qpp)^2, and q_tp as-is.

@param[in,out] commondata Pointer to common data structure containing shared parameters and settings.
@param[in,out] griddata Pointer to grid data structures for each grid, containing parameters and gridfunctions.
@return BHAHAHA_SUCCESS, an interpolation error, a circumference allocation
error, or DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT.
@note Uses OpenMP for parallel loops; performs interpolation and midpoint integration over angular samples.
@note One deterministic polar meridian is used. Distorted or finite-resolution
horizons can therefore retain orientation sensitivity about the spin axis.

Precomputation strategy on the (theta,phi) grid at fixed i0=NGHOSTS:
  q_{theta theta}, q_{phi phi}, and q_{theta phi} are generated through
  area.py::circumference_metric_roots. Scalar f_eq(theta,phi) and
  f_pol(theta,phi) integrands are then formed with analytic tangent fields.
  Only those scalars are ghost-filled and interpolated, avoiding angle branch
  cuts and finite differencing of spherical angles.
"""
    cfunc_type = "int"
    name = "diagnostics_proper_circumferences_general"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;
  bhahaha_diags->BHAHAHA_CIRC_GENERAL_POLAR = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
  bhahaha_diags->BHAHAHA_CIRC_GENERAL_EQUATOR = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
  bhahaha_diags->BHAHAHA_CIRC_GENERAL_SPIN = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
  bhahaha_diags->spin_chi_x_circumference = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
  bhahaha_diags->spin_chi_y_circumference = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;
  bhahaha_diags->spin_chi_z_circumference = BHAHAHA_DIAGNOSTIC_UNAVAILABLE;

  // Require the direction selected by the SpECTRE diagnostic. Never invent a
  // fallback direction, because that would change the meaning of this result.
  REAL s[3] = {
      bhahaha_diags->BHAHAHA_SPIN_AXIS_X,
      bhahaha_diags->BHAHAHA_SPIN_AXIS_Y,
      bhahaha_diags->BHAHAHA_SPIN_AXIS_Z};
  if (s[0] == BHAHAHA_DIAGNOSTIC_UNAVAILABLE || s[1] == BHAHAHA_DIAGNOSTIC_UNAVAILABLE ||
      s[2] == BHAHAHA_DIAGNOSTIC_UNAVAILABLE || !isfinite(s[0]) || !isfinite(s[1]) || !isfinite(s[2]))
    return DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT;
  const REAL s_norm = norm3(s);
  if (!(s_norm > 1.0e-14) || !isfinite(s_norm))
    return DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT;
  s[0] /= s_norm;
  s[1] /= s_norm;
  s[2] /= s_norm;

  // Build an orthonormal basis {e1, e2} orthogonal to s. The polar circle is
  // the single deterministic meridian in the plane spanned by s and e1.
  REAL e1[3], e2[3];
  build_basis_from_s(s, e1, e2);
  REAL nvec[3];
  cross3(s, e1, nvec);

  // which_gf in {0: sqrt(q_tt), 1: sqrt(q_pp), 2: q_tp, 3: f_eq(θ,φ), 4: f_pol(θ,φ)};
  // diagonals as sqrt for stability; f_eq/f_pol are scalar line-element integrands.
  const int NUM_DIAG_GFS = 5;
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.
  const int N_angle = griddata[grid].params.Nxx2;
  if (Nxx_plus_2NGHOSTS0 <= 0 || Nxx_plus_2NGHOSTS1 <= 0 || Nxx_plus_2NGHOSTS2 <= 0 || N_angle <= 0)
    return DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT;

  REAL *restrict metric_data_gfs = NULL;
  BHAH_MALLOC(metric_data_gfs,
              (size_t)Nxx_plus_2NGHOSTS0 * (size_t)Nxx_plus_2NGHOSTS1 *
                  (size_t)Nxx_plus_2NGHOSTS2 * (size_t)NUM_DIAG_GFS * sizeof(REAL));
  if (metric_data_gfs == NULL)
    return DIAG_PROPER_CIRCUM_MALLOC_ERROR;

  // Compute sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} across the entire (theta,phi) grid at i0=NGHOSTS.
  {
    // Extract pointers to auxiliary and evolved gridfunctions, and coordinate arrays.
    const REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
    const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
    const REAL *restrict xx[3] = {griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]};
    // Inverse grid spacings for theta and phi directions.
    const REAL invdxx1 = griddata[grid].params.invdxx1;
    const REAL invdxx2 = griddata[grid].params.invdxx2;
    const int i0 = NGHOSTS; // Fixed index for radial coordinate (r) where r=h(theta,phi) lives.

    // Loop over angular grid points (theta and phi) to compute:
    // 1. sqrt(q_{theta theta}), stored to metric_data_gfs[IDX4(0,...)], and
    // 2. sqrt(q_{phi phi}),    stored to metric_data_gfs[IDX4(1,...)] at each point (theta, phi).
    // 3. q_{theta phi},       stored to metric_data_gfs[IDX4(2,...)]. Keep sign; no sqrt.
    // Clever IDX math ensures this 2D computation stays within the memory bounds of the 3D allocation.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED REAL xx2 = xx[2][i2]; // Phi coordinate at index i2.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED REAL xx1 = xx[1][i1]; // Theta coordinate at index i1.
"""
    # Generate the three outputs from the induced 2-metric roots:
    body += ccg.c_codegen(
        bhahaha_area.circumference_metric_roots(),
        [
            "metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)]",  # sqrt(q_tt)
            "metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)]",  # sqrt(q_pp)
            "metric_data_gfs[IDX4pt(2, 0) + IDX2(i1, i2)]",  # q_tp
        ],
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )
    body += r"""
      } // END LOOP: for i1 over theta points on the horizon surface
    } // END LOOP: for i2 over phi points on the horizon surface

  } // END BLOCK: compute q-metric root gridfunctions on the horizon surface

  // Uniform alpha step for midpoint samples in [-pi, pi).
  const REAL d_alpha = (M_PI - (-M_PI)) / ((REAL)N_angle);

  REAL dst_pts[N_angle][2];
  REAL theta[N_angle];
  REAL phi[N_angle];
  REAL integrand[N_angle];

  // ================================================================
  // Precompute integrand scalars f_eq(θ,φ) and f_pol(θ,φ) on the i0=NGHOSTS slab.
  // Each uses a smooth tangent field on the unit sphere:
  //   equator: t = normalize(s × r),   polar: t = normalize((s × e1) × r)
  // and the spherical basis (θ̂, φ̂) to build dθ/dalpha = t·θ̂, dφ/dalpha = (t·φ̂)/sinθ.
  // This avoids any finite differencing of angles and eliminates branch cuts.
  // ================================================================
#pragma omp parallel for
  for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
    const REAL phi_c = griddata[grid].xx[2][i2];
    const REAL sinph = sin(phi_c), cosph = cos(phi_c);
    for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
      const REAL theta_c = griddata[grid].xx[1][i1];
      const REAL sinth = sin(theta_c), costh = cos(theta_c);
      // unit position and spherical basis
      const REAL rx = sinth * cosph, ry = sinth * sinph, rz = costh;
      const REAL thx = costh * cosph, thy = costh * sinph, thz = -sinth;
      const REAL phx = -sinph, phy = cosph, phz = 0.0;

      // load q-metric pieces (recall diagonals are stored as square-roots)
      const REAL sqtt = metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)];
      const REAL sqpp = metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)];
      const REAL qtp = metric_data_gfs[IDX4pt(2, 0) + IDX2(i1, i2)];
      const REAL qtt = sqtt * sqtt;
      const REAL qpp = sqpp * sqpp;

      // ---------- equator tangent: t = normalize(s × r) ----------
      REAL tex = s[1] * rz - s[2] * ry;
      REAL tey = s[2] * rx - s[0] * rz;
      REAL tez = s[0] * ry - s[1] * rx;
      REAL tnorm = sqrt(tex * tex + tey * tey + tez * tez);
      if (tnorm > 1e-14) {
        tex /= tnorm;
        tey /= tnorm;
        tez /= tnorm;
      } else {
        tex = tey = tez = 0.0;
      }
      const REAL dth_eq = tex * thx + tey * thy + tez * thz;
      const REAL dph_eq = (tex * phx + tey * phy + tez * phz) / fmax(1e-14, sinth);
      const REAL feq = sqrt(qtt * dth_eq * dth_eq + 2.0 * qtp * dth_eq * dph_eq + qpp * dph_eq * dph_eq);
      metric_data_gfs[IDX4pt(3, 0) + IDX2(i1, i2)] = feq;

      // ---------- polar tangent: plane normal n = s × e1; t = normalize(n × r) ----------
      REAL tpx = nvec[1] * rz - nvec[2] * ry;
      REAL tpy = nvec[2] * rx - nvec[0] * rz;
      REAL tpz = nvec[0] * ry - nvec[1] * rx;
      tnorm = sqrt(tpx * tpx + tpy * tpy + tpz * tpz);
      if (tnorm > 1e-14) {
        tpx /= tnorm;
        tpy /= tnorm;
        tpz /= tnorm;
      } else {
        tpx = tpy = tpz = 0.0;
      }
      const REAL dth_pol = tpx * thx + tpy * thy + tpz * thz;
      const REAL dph_pol = (tpx * phx + tpy * phy + tpz * phz) / fmax(1e-14, sinth);
      const REAL fpol = sqrt(qtt * dth_pol * dth_pol + 2.0 * qtp * dth_pol * dph_pol + qpp * dph_pol * dph_pol);
      metric_data_gfs[IDX4pt(4, 0) + IDX2(i1, i2)] = fpol;
    } // END LOOP: for i1 over theta
  } // END LOOP: for i2 over phi

  // Apply inner boundary conditions to the newly computed scalar integrands f_eq (which_gf=3)
  // and f_pol (which_gf=4), so their ghost zones are valid prior to interpolation.
  {
    const int which_gfs_2[2] = {3, 4};
    apply_inner_bc_for_selected_gfs(&griddata[grid].bcstruct, metric_data_gfs, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                    which_gfs_2, 2);
  } // END BLOCK: apply inner boundary conditions to equatorial and polar integrands

  // Convenience base pointers for scalar integrands:
  const REAL *restrict src_feq = &metric_data_gfs[IDX4pt(3, 0)];
  const REAL *restrict src_fpol = &metric_data_gfs[IDX4pt(4, 0)];

  // Some compilers warn on constness; build an explicit coords pointer matching prototype.
  REAL *restrict(*coords)[3] = (REAL *restrict(*)[3]) & griddata[grid].xx;

  // ================================================================
  // Equatorial great circle: orthogonal to s.
  //   Parameterization: r(alpha) = cos(alpha) e1 + sin(alpha) e2, alpha in [-pi, pi).
  //   Convert to (theta, phi), interpolate precomputed scalar integrand, and integrate.
  // ================================================================
  // #pragma omp parallel for // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++) {
    const REAL alpha = -M_PI + ((REAL)i + 0.5) * d_alpha;
    REAL rvec[3] = {cos(alpha) * e1[0] + sin(alpha) * e2[0], cos(alpha) * e1[1] + sin(alpha) * e2[1], cos(alpha) * e1[2] + sin(alpha) * e2[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  } // END LOOP: for i over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed equator integrand f_eq onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (REAL *restrict *)(*coords), src_feq, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END IF: equatorial-integrand interpolation failed
  } // END BLOCK: interpolate equatorial integrand onto great-circle samples

  const REAL C_equator = integrate_over_alpha(integrand, N_angle, d_alpha);

  // ================================================================
  // Polar great circle: contains s.
  //   Parameterization: r(alpha) = cos(alpha) s + sin(alpha) e1, alpha in [-pi, pi).
  //   Proceed as above.
  // ================================================================
  // #pragma omp parallel for reduction(+ : sum) // <- N_angle ~64 => thread creation/destruction makes OMP SLOWER
  for (int i = 0; i < N_angle; i++) {
    const REAL alpha = -M_PI + ((REAL)i + 0.5) * d_alpha;
    REAL rvec[3] = {cos(alpha) * s[0] + sin(alpha) * e1[0], cos(alpha) * s[1] + sin(alpha) * e1[1], cos(alpha) * s[2] + sin(alpha) * e1[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  } // END LOOP: for i over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed polar integrand f_pol onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (REAL *restrict *)(*coords), src_fpol, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END IF: polar-integrand interpolation failed
  } // END BLOCK: interpolate polar integrand onto great-circle samples

  const REAL C_polar = integrate_over_alpha(integrand, N_angle, d_alpha);

  if (!(C_equator > 0.0) || !(C_polar > 0.0) || !isfinite(C_equator) || !isfinite(C_polar)) {
    free(metric_data_gfs);
    return DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT;
  } // END IF: integrated circumferences are invalid
  const REAL ratio = C_polar / C_equator;
  const REAL chi_circumference = compute_spin(ratio);
  if (chi_circumference == BHAHAHA_DIAGNOSTIC_UNAVAILABLE || !isfinite(chi_circumference) ||
      chi_circumference < 0.0 || chi_circumference > 1.0) {
    free(metric_data_gfs);
    return DIAG_PROPER_CIRCUM_GENERAL_INVALID_INPUT;
  } // END IF: circumference ratio has no physical Kerr inversion

  bhahaha_diags->BHAHAHA_CIRC_GENERAL_POLAR = C_polar;
  bhahaha_diags->BHAHAHA_CIRC_GENERAL_EQUATOR = C_equator;
  bhahaha_diags->BHAHAHA_CIRC_GENERAL_SPIN = chi_circumference;
  bhahaha_diags->spin_chi_x_circumference = chi_circumference * s[0];
  bhahaha_diags->spin_chi_y_circumference = chi_circumference * s[1];
  bhahaha_diags->spin_chi_z_circumference = chi_circumference * s[2];

  free(metric_data_gfs);

  return BHAHAHA_SUCCESS; // Return success status code.
"""
    cfc.register_CFunction(
        subdirectory="",
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

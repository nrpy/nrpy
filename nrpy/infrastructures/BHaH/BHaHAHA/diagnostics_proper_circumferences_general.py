# BHaHAHA/diagnostics_proper_circumferences.py
"""
Register and generate a C routine that computes equatorial and polar proper circumferences of a horizon surface.

We integrate the arclength along great circles on the 2D surface r = h(theta, phi) using the induced 2-metric q_{AB}.
Let alpha parameterize a great circle on the unit sphere (so theta = theta(alpha), phi = phi(alpha)) and let overdots
denote derivatives w.r.t. alpha. The physical line element we integrate is:

    ds = sqrt( q_{theta theta} * (dtheta/dalpha)^2
               + 2 * q_{theta phi} * (dtheta/dalpha) * (dphi/dalpha)
               +     q_{phi phi}   * (dphi/dalpha)^2 ) dalpha.

Design choices that might look odd at first glance, and why they are intentional:

* Store/interpolate sqrt(q_{theta theta}) and sqrt(q_{phi phi}) but not sqrt(q_{theta phi}).
  The diagonals are nonnegative “scale factors” along coordinate lines, so interpolating their square roots reduces
  dynamic range and suppresses negative overshoots from interpolation; we then square them at use sites to recover q_tt, q_pp.
  The off-diagonal q_{theta phi} can be positive or negative and enters linearly in the cross term; taking a square root would
  both be undefined for negative values and discard the sign that the cross term needs.

* Work on the i0 = NGHOSTS radial slice.
  In this infrastructure, r = h(theta, phi) lives at a fixed interior radial index (the horizon layer). Precomputing
  the induced 2-metric ingredients on that slice keeps memory traffic low and avoids 3D interpolation for a 2D integral.

* Unwrap phi before differencing.
  Great circles cross the branch cut at +/- pi; naive finite differences across that cut generate O(2*pi) jumps that
  corrupt the ds integrand. We explicitly unwrap phi(alpha) prior to forming central differences to maintain continuity.

* Use 8th-order midpoint integration with precomputed weights.
  The ds integrand is smooth along alpha; a fixed stencil midpoint rule with tuned weights offers excellent accuracy per
  evaluation, vectorizes well, and avoids the overhead/complexity of adaptive quadrature in tight inner loops.

* Use stack-allocated variable-length arrays (VLAs) for per-angle work buffers.
  N_angle is known at runtime, and VLAs keep allocation overhead minimal and cache locality high inside OpenMP regions.

* Estimate spin from C_polar/C_equator with a safeguarded Newton–Raphson.
  We start from an analytic approximation (Alcubierre et al., Eq. 5.3) and refine with a bounded/robust NR step computed
  from elliptic integrals; we clamp out-of-range updates and provide a clear failure code when convergence is not achieved.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
from nrpy.equations.general_relativity.bhahaha.ExpansionFunctionTheta import (
    ExpansionFunctionTheta,
)
from nrpy.infrastructures import BHaH


# circumferential_arclength: Calculates the differential arclength element along a specified direction (theta or phi).
def circumference_metric_roots() -> Tuple[sp.Expr, sp.Expr, sp.Expr]:
    """
    Compute the induced-metric ingredients for the proper-circumference integrand.

    We construct q_{AB} on the embedded surface r = h(theta, phi) in terms of the ambient spatial metric gamma_{ij}
    and angular derivatives of h. Explicitly,
        q_{theta theta} = gamma_{rr} h_{,theta}^2 + 2 gamma_{r theta} h_{,theta} + gamma_{theta theta}
        q_{phi phi}     = gamma_{rr} h_{,phi}^2  + 2 gamma_{r phi}  h_{,phi}  + gamma_{phi phi}
        q_{theta phi}   = gamma_{rr} h_{,theta} h_{,phi} + gamma_{r theta} h_{,phi}
                          + gamma_{r phi} h_{,theta} + gamma_{theta phi}.

    We return sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} (not sqrt) for the following reasons:
      (1) The diagonals are nonnegative scale factors; interpolating their square roots improves conditioning and guards
          against tiny negative overshoots that would otherwise cause NaNs in the outer sqrt(ds^2).
      (2) The cross-term q_{theta phi} changes sign and enters linearly in the ds integrand; preserving its sign is essential.

    :return: A tuple of SymPy expressions (sqrt_qtt, sqrt_qpp, qtp) evaluated with xx0 -> h substitutions applied.
    """
    Th = ExpansionFunctionTheta["Spherical"]
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD")

    # Induced 2-metric components on r = h(theta,phi), expressed on the spherical reference metric.
    qtt = (
        Th.gammaDD[0][0] * h_dD[1] ** 2
        + 2 * Th.gammaDD[0][1] * h_dD[1]
        + Th.gammaDD[1][1]
    )
    qpp = (
        Th.gammaDD[0][0] * h_dD[2] ** 2
        + 2 * Th.gammaDD[0][2] * h_dD[2]
        + Th.gammaDD[2][2]
    )
    qtp = (
        Th.gammaDD[0][0] * h_dD[1] * h_dD[2]
        + Th.gammaDD[0][1] * h_dD[2]
        + Th.gammaDD[0][2] * h_dD[1]
        + Th.gammaDD[1][2]
    )

    # Replace the symbolic radius placeholder xx0 with h in all expressions, then
    # return sqrt(diagonals) and raw off-diagonal as discussed above.
    sqrt_qtt = sp.sqrt(qtt.replace(sp.sympify("xx0"), h))
    sqrt_qpp = sp.sqrt(qpp.replace(sp.sympify("xx0"), h))
    qtp = qtp.replace(sp.sympify("xx0"), h)
    return sqrt_qtt, sqrt_qpp, qtp


def register_CFunction_diagnostics_proper_circumferences(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C routine that computes equatorial/polar proper circumferences and a spin estimate from their ratio.

    This generates a C function that:
      * precomputes sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} on the (theta,phi) grid at i0=NGHOSTS,
      * samples great circles defined by the provided spin axis s^i,
      * builds the line element
            ds = sqrt( q_tt*(dtheta/dalpha)^2 + 2*q_tp*(dtheta/dalpha)*(dphi/dalpha) + q_pp*(dphi/dalpha)^2 ) dalpha,
        with q_tt = (sqrt_qtt)^2 and q_pp = (sqrt_qpp)^2,
      * integrates ds over alpha in [-pi, pi) to obtain C_equator and C_polar, and
      * estimates a spin magnitude a/M from C_polar/C_equator.

    :param enable_fd_functions: Whether to emit FD derivative helpers in the generated code (passed through to codegen).
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    DocTests:
    >>> import nrpy.grid as gri
    >>> _ = gri.register_gridfunctions("hh")[0]
    >>> env = register_CFunction_diagnostics_proper_circumferences()
    Setting up ExpansionFunctionThetaClass[Spherical]...
    Setting up reference_metric[Spherical]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
/**
 * Elliptic integrals and helpers for circumference diagnostics.
 *
 * The ds integrand used later is:
 *   ds = sqrt( q_tt * (dtheta/dalpha)^2 + 2*q_tp*(dtheta/dalpha)*(dphi/dalpha) + q_pp*(dphi/dalpha)^2 ) dalpha,
 * with q_tt = (sqrt_qtt)^2, q_pp = (sqrt_qpp)^2, and q_tp as-is. We store sqrt(diagonal) for stability/positivity.
 */

/**
 * Computes the complete elliptic integrals of the second kind E(k) and the first kind K(k)
 * using 8th-order midpoint integration.
 *
 * Rationale: The integrands are smooth and bounded on [0, pi/2]; a fixed high-order midpoint stencil with precomputed
 *            weights delivers excellent accuracy per function evaluation, vectorizes well, and avoids adaptive overhead
 *            in tight loops. 128 sample points is a conservative default providing ample headroom for convergence tests.
 *
 * @param k - The elliptic modulus parameter (0 <= k <= 1).
 * @param E - Pointer to store the result for the elliptic integral of the second kind E(k).
 * @param K - Pointer to store the result for the elliptic integral of the first kind K(k).
 *
 * @note The function is parallelized with OpenMP to improve performance on large arrays.
 */
static void elliptic_E_and_K_integrals(const REAL k, REAL *restrict E, REAL *restrict K) {
  static const int N_sample_pts = 128; // High-accuracy default; power-of-two to simplify weight cycling.
  const REAL *restrict weights;        // 8th-order midpoint weights (periodic stencil).
  int weight_stencil_size;             // Size of the periodic weight stencil.

  // Retrieve the integration weights based on the number of sample points (divisible by 8 -> 8th order).
  bah_diagnostics_integration_weights(N_sample_pts, N_sample_pts, &weights, &weight_stencil_size);

  const REAL a = 0.0;                            // Lower limit (0).
  const REAL b = M_PI / 2.0;                     // Upper limit (pi/2).
  const REAL h = (b - a) / ((REAL)N_sample_pts); // Midpoint step.

  REAL sum_E = 0.0; // Accumulator for E(k).
  REAL sum_K = 0.0; // Accumulator for K(k).

  // Parallelized midpoint rule with periodic weights.
#pragma omp parallel for reduction(+ : sum_E, sum_K)
  for (int i = 0; i < N_sample_pts; i++) {
    const REAL theta = a + ((REAL)i + 0.5) * h; // midpoint
    const REAL sintheta = sin(theta);
    const REAL elliptic_E_integrand = sqrt(1.0 - k * sintheta * sintheta);
    const REAL elliptic_K_integrand = 1.0 / elliptic_E_integrand;
    sum_E += weights[i % weight_stencil_size] * elliptic_E_integrand;
    sum_K += weights[i % weight_stencil_size] * elliptic_K_integrand;
  }

  *E = sum_E * h;
  *K = sum_K * h;
}

/**
 * Estimates the spin parameter magnitude for equilibrium black holes based on the circumference ratio C_r.
 *
 * Rationale & safeguards:
 *   * Start from an analytic approximation (Alcubierre et al., Eq. 5.3) to land near the root basin.
 *   * Refine with Newton–Raphson using elliptic integrals (Eq. 5.2), but clamp any out-of-range iterates
 *     into [0,1] and terminate if we step negative; this avoids excursions where the modulus or square roots
 *     would be undefined or numerically fragile.
 *
 * @param C_r The circumference ratio C_polar/C_equator used to estimate the spin.
 * @return    The estimated spin parameter. Returns -10.0 if C_r is out of valid bounds or if convergence fails.
 */
REAL compute_spin(const REAL C_r) {
  // Reject obviously invalid input.
  if (C_r > 1)
    return -10.0;

  // Conservative initial guess (updated below by analytic approximation).
  REAL spin = 0.9;

  // Analytic approximation (Alcubierre et al., Eq. 5.3).
  const REAL spin_sq = 1 - (2.55 * C_r - 1.55) * (2.55 * C_r - 1.55);
  if (spin_sq >= 0 && spin_sq < 1)
    spin = sqrt(spin_sq);

  const REAL rel_tolerance = 1e-7; // Tight tolerance is cheap here.
  REAL rel_diff = 1e10;            // Start far away.
  const int max_its = 20;          // Hard cap to avoid rare pathologies.
  int it = 0;

  // Safeguarded Newton–Raphson iteration.
  while (rel_diff > rel_tolerance && it < max_its) {
    const REAL x = spin;
    REAL E, K;

    // Compute elliptic integrals at the current iterate; modulus choice matches the analytic mapping.
    elliptic_E_and_K_integrals(-((x * x) / pow(1 + sqrt(1 - (x * x)), 2)), &E, &K);

    // One Newton–Raphson update (generated by NRPy from a symbolic expression).
"""
    # Newton-Raphson is a bit more robust; also our initial guess is pretty good, so typically we need only a few iterations.
    prefunc += ccg.c_codegen(
        BHaH.BHaHAHA.area.spin_NewtonRaphson(), "const REAL x_np1", include_braces=False
    )
    prefunc += r"""
    if (x_np1 > 1.0) {
      // Reflect back into [0,1] if we overshoot; keeps modulus well-defined.
      spin = 2.0 - x_np1;
    } else if (x_np1 < 0) {
      // Bail out if we step negative; report failure below.
      it = max_its;
      break;
    } else {
      // Accept step and measure relative change.
      rel_diff = fabs(x_np1 - x) / x;
      spin = x_np1;
    }
    it++;
  }

  // Report a clear failure code if we couldn't converge.
  if (it >= max_its) {
    spin = -10.0;
  }
  return spin;
}

#ifndef BHAHAHA_FAILURE
#define BHAHAHA_FAILURE (-1)
#endif

// ---------- small vector & math helpers (pure C) ----------
// These are intentionally tiny & inlined: they live in tight OpenMP loops,
// and we want predictable codegen and good autovectorization.
static inline REAL dot3(const REAL a[3], const REAL b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline void cross3(const REAL a[3], const REAL b[3], REAL c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
static inline REAL norm3(const REAL a[3]) { return sqrt(dot3(a,a)); }
static inline void normalize3(REAL a[3]) { const REAL n = norm3(a); if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; } }
// Wrap dphi into (-pi, pi]; helps robustly unwrap cumulative angles later.
static inline REAL wrap_dphi(REAL dphi) { while (dphi > M_PI) dphi -= 2.0*M_PI; while (dphi <= -M_PI) dphi += 2.0*M_PI; return dphi; }
// Build an orthonormal basis {e1, e2} orthogonal to s; avoid near-collinearity for numerical stability.
static inline void build_basis_from_s(const REAL s[3], REAL e1[3], REAL e2[3]) {
  REAL a[3] = {1.0, 0.0, 0.0};
  if (fabs(dot3(a,s)) > 0.9) { a[0]=0.0; a[1]=1.0; a[2]=0.0; }
  cross3(s, a, e1); normalize3(e1);
  cross3(s, e1, e2); normalize3(e2);
}
// Cartesian unit vector -> spherical angles (theta, phi); inputs here are unit by construction.
static inline void cart_to_sph(const REAL r[3], REAL *theta, REAL *phi) {
  const REAL x=r[0], y=r[1], z=r[2];
  REAL zz = z; if (zz < -1.0) zz = -1.0; if (zz > 1.0) zz = 1.0;
  *theta = acos(zz);
  *phi   = atan2(y, x);
}

// Compute dtheta/dalpha and dphi/dalpha with central differences.
// IMPORTANT: We unwrap phi(alpha) first to avoid the +/- pi branch cut. Without unwrapping,
// differences across the cut produce O(2*pi) jumps that would pollute the ds integrand.
static inline void derivs_wrt_alpha(const REAL *theta, const REAL *phi,
                                    REAL *dtheta_dalpha, REAL *dphi_dalpha,
                                    int N, REAL d_alpha) {
  REAL *phi_unwrap = (REAL *)malloc((size_t)N * sizeof(REAL));
  if (phi_unwrap == NULL) {
    // Fallback: difference across wrapped phi using wrap_dphi on the fly.
    for (int i = 0; i < N; i++) {
      const int im = (i - 1 + N) % N, ip = (i + 1) % N;
      dtheta_dalpha[i] = (theta[ip] - theta[im]) / (2.0*d_alpha);
      dphi_dalpha[i]   = wrap_dphi(phi[ip] - phi[im]) / (2.0*d_alpha);
    }
    return;
  }
  // Cumulative unwrap: ensures continuity over the full loop.
  phi_unwrap[0] = phi[0];
  for (int i = 1; i < N; i++) {
    const REAL d = wrap_dphi(phi[i] - phi[i-1]);
    phi_unwrap[i] = phi_unwrap[i-1] + d;
  }
  for (int i = 0; i < N; i++) {
    const int im = (i - 1 + N) % N, ip = (i + 1) % N;
    dtheta_dalpha[i] = (theta[ip] - theta[im]) / (2.0*d_alpha);
    dphi_dalpha[i]   = (phi_unwrap[ip] - phi_unwrap[im]) / (2.0*d_alpha);
  }
  free(phi_unwrap);
}

// Midpoint integrate over alpha with precomputed 8th-order weights.
//  We keep this minimal and branchless inside the OpenMP loop for best vectorization.
static inline REAL integrate_over_alpha(const REAL *vals, int N_angle, REAL d_alpha) {
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);
  REAL sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < N_angle; i++) sum += vals[i] * weights[i % weight_stencil_size];
  return sum * d_alpha;
}

// interpolate sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} to dst_pts[(theta,phi)].
//  Note the diagonals are stored as their square roots to enforce positivity & improve conditioning.
static inline int interp_qroots(const griddata_struct *grid,
                                REAL *restrict (*coords)[3],
                                const REAL *restrict src_qtt,
                                const REAL *restrict src_qpp,
                                const REAL *restrict src_qtp,
                                const REAL (*dst_pts)[2], int N_angle,
                                REAL *sqrt_qtt, REAL *sqrt_qpp, REAL *qtp) {
  const int Nxx_plus_2NGHOSTS1 = grid->params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = grid->params.Nxx_plus_2NGHOSTS2;
  const REAL dxx1 = grid->params.dxx1;
  const REAL dxx2 = grid->params.dxx2;

  int err = bah_interpolation_2d_general__uniform_src_grid(
      NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
      (REAL *restrict *)(*coords), src_qtt, N_angle, dst_pts, sqrt_qtt);
  if (err != BHAHAHA_SUCCESS) return err;

  err = bah_interpolation_2d_general__uniform_src_grid(
      NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
      (REAL *restrict *)(*coords), src_qpp, N_angle, dst_pts, sqrt_qpp);
  if (err != BHAHAHA_SUCCESS) return err;

  err = bah_interpolation_2d_general__uniform_src_grid(
      NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
      (REAL *restrict *)(*coords), src_qtp, N_angle, dst_pts, qtp);
  return err;
}
/* ---------- end helpers ---------- */
"""
    desc = """
Computes proper circumferences along the equator and polar directions with respect to the provided spin axis,
using the induced 2-metric q_{AB} on the surface r = h(theta,phi).

Line element integrated (alpha parametrizes the great circle):
    ds = sqrt( q_tt * (dtheta/dalpha)^2 + 2*q_tp*(dtheta/dalpha)*(dphi/dalpha) + q_pp * (dphi/dalpha)^2 ) dalpha,
with q_tt = (sqrt_qtt)^2, q_pp = (sqrt_qpp)^2, and q_tp as-is. Storing sqrt(diagonals) improves interpolation stability
and preserves positivity; keeping q_tp raw preserves the sign required by the cross term.

@param commondata - Pointer to common data structure containing shared parameters and settings.
@param griddata - Pointer to grid data structures for each grid, containing parameters and gridfunctions.
@return - Status code indicating success or type of error (e.g., BHAHAHA_SUCCESS or INITIAL_DATA_MALLOC_ERROR).
@note - Uses OpenMP for parallel loops; performs interpolation and midpoint integration over angular samples.

Precomputation strategy on the (theta,phi) grid at fixed i0=NGHOSTS:
  sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} are generated via NRPy's SymPy expressions
  (BHaHAHA/area.py::circumference_metric_roots), then interpolated along great circles defined by the unit spin axis s^i.

We finally store C_polar, C_equator, their ratio (C_polar/C_equator), and an estimated spin a/M based on the ratio.
"""
    cfunc_type = "int"
    name = "diagnostics_proper_circumferences"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  const int NUM_DIAG_GFS = 3; // which_gf in {0: sqrt(q_tt), 1: sqrt(q_pp), 2: q_tp}; diagonals as sqrt for stability.
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  REAL *restrict metric_data_gfs;
  // Single heap allocation sized for a 2D (theta,phi) slab at i0 = NGHOSTS.
  BHAH_MALLOC(metric_data_gfs, Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS);

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
        list(circumference_metric_roots()),
        [
            "metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)]",  # sqrt(q_tt)
            "metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)]",  # sqrt(q_pp)
            "metric_data_gfs[IDX4pt(2, 0) + IDX2(i1, i2)]",  # q_tp
        ],
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )
    body += r"""
      } // END LOOP over i1 (theta)
    } // END LOOP over i2 (phi)

    // Apply inner boundary conditions to the computed gridfunctions.
    {
      bc_struct *restrict bcstruct = &griddata[grid].bcstruct; // Retrieve boundary condition structure for the grid.

      // Unpack bc_info from bcstruct.
      const bc_info_struct *bc_info = &bcstruct->bc_info;

      // Apply boundary conditions at inner boundary points for the selected gridfunctions.
#pragma omp parallel for collapse(2)
      for (int which_gf = 0; which_gf < 3; which_gf++) {
        for (int pt = 0; pt < bc_info->num_inner_boundary_points; pt++) {
          const int dstpt = bcstruct->inner_bc_array[pt].dstpt; // Destination point index.
          const int srcpt = bcstruct->inner_bc_array[pt].srcpt; // Source point index for copying.

          // Decode (i0,i1,i2) from the flattened idx3; only copy if we are on the i0 == NGHOSTS slab.
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
        } // END LOOP over inner boundary points
      } // END LOOP over gridfunctions
    } // END application of inner boundary conditions
  } // END computation of q-metric root gridfunctions

  // Number of angular points to sample over 2 pi radians (controls resolution of great-circle integrals).
  const int N_angle = griddata[grid].params.Nxx2;
  // Uniform alpha step for midpoint samples in [-pi, pi).
  const REAL d_alpha = (M_PI - (-M_PI)) / ((REAL)N_angle);

  // C99 variable-length arrays on the stack (per-thread in OpenMP); avoids heap overhead in tight loops.
  REAL dst_pts[N_angle][2];
  REAL theta[N_angle];
  REAL phi[N_angle];
  REAL dth[N_angle];
  REAL dph[N_angle];
  REAL sqrt_qtt[N_angle];
  REAL sqrt_qpp[N_angle];
  REAL qtp[N_angle];
  REAL integrand[N_angle];

  // Normalize spin axis; if zero-length is provided, fall back to z-axis for determinism.
  REAL s[3] = {
    commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_X, // unit spin direction, x-component
    commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Y, // unit spin direction, y-component
    commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Z  // unit spin direction, z-component
  };
  REAL s_norm = sqrt(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
  if (s_norm > 0.0) {
    s[0] /= s_norm;
    s[1] /= s_norm;
    s[2] /= s_norm;
  } else {
    s[0] = 0.0;
    s[1] = 0.0;
    s[2] = 1.0;
  } // fallback: z-axis
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_X = s[0];
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Y = s[1];
  commondata->bhahaha_diagnostics->BHAHAHA_SPIN_AXIS_Z = s[2];

  // Build an orthonormal basis {e1, e2} orthogonal to s, to define great circles.
  REAL e1[3], e2[3];
  build_basis_from_s(s, e1, e2);

  // Convenience base pointers for interpolation sources; diagonals are stored as square roots by design.
  const REAL *restrict src_qtt = &metric_data_gfs[IDX4pt(0, 0)]; // sqrt(q_{theta theta})
  const REAL *restrict src_qpp = &metric_data_gfs[IDX4pt(1, 0)]; // sqrt(q_{phi phi})
  const REAL *restrict src_qtp = &metric_data_gfs[IDX4pt(2, 0)]; // q_{theta phi}

  // Some compilers warn on constness; build an explicit coords pointer matching prototype.
  REAL *restrict(*coords)[3] = (REAL *restrict(*)[3]) & griddata[grid].xx;

  // ================================================================
  // Equatorial great circle: orthogonal to s.
  //   Parameterization: r(alpha) = cos(alpha) e1 + sin(alpha) e2, alpha in [-pi, pi).
  //   Convert to (theta, phi), interpolate q-roots, compute derivatives wrt alpha,
  //   build the full line element including q_{theta phi}, and integrate.
  // ================================================================
#pragma omp parallel for
  for (int i = 0; i < N_angle; i++) {
    const REAL alpha = -M_PI + ((REAL)i + 0.5) * d_alpha;
    REAL rvec[3] = {cos(alpha) * e1[0] + sin(alpha) * e2[0],
                    cos(alpha) * e1[1] + sin(alpha) * e2[1],
                    cos(alpha) * e1[2] + sin(alpha) * e2[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  }

  {
    // Interpolate sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} onto the equator points.
    const int err = interp_qroots(&griddata[grid], coords, src_qtt, src_qpp, src_qtp,
                                  dst_pts, N_angle, sqrt_qtt, sqrt_qpp, qtp);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    }
  }

  // Compute dtheta/dalpha and dphi/dalpha along the equator; phi is unwrapped inside this helper.
  derivs_wrt_alpha(theta, phi, dth, dph, N_angle, d_alpha);

  // Build the proper line-element integrand and integrate to obtain C_equator.
  // NOTE: q_tt and q_pp are reconstructed by squaring the stored square roots; q_tp is used as-is to preserve sign.
#pragma omp parallel for
  for (int i = 0; i < N_angle; i++) {
    const REAL qtt = sqrt_qtt[i] * sqrt_qtt[i];
    const REAL qpp = sqrt_qpp[i] * sqrt_qpp[i];
    // ds = sqrt( qtt * dth^2 + 2*qtp * dth*dph + qpp * dph^2 )
    integrand[i] = sqrt(qtt * (dth[i] * dth[i]) + 2.0 * qtp[i] * (dth[i] * dph[i]) + qpp * (dph[i] * dph[i]));
  }
  const REAL C_equator = integrate_over_alpha(integrand, N_angle, d_alpha);

  // ================================================================
  // Polar great circle: contains s.
  //   Parameterization: r(alpha) = cos(alpha) s + sin(alpha) e1, alpha in [-pi, pi).
  //   Proceed as above.
  // ================================================================
#pragma omp parallel for
  for (int i = 0; i < N_angle; i++) {
    const REAL alpha = -M_PI + ((REAL)i + 0.5) * d_alpha;
    REAL rvec[3] = {cos(alpha) * s[0] + sin(alpha) * e1[0],
                    cos(alpha) * s[1] + sin(alpha) * e1[1],
                    cos(alpha) * s[2] + sin(alpha) * e1[2]};
    cart_to_sph(rvec, &theta[i], &phi[i]);
    dst_pts[i][0] = theta[i];
    dst_pts[i][1] = phi[i];
  }

  {
    // Interpolate sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} onto the polar points.
    const int err = interp_qroots(&griddata[grid], coords, src_qtt, src_qpp, src_qtp,
                                  dst_pts, N_angle, sqrt_qtt, sqrt_qpp, qtp);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    }
  }

  // Compute dtheta/dalpha and dphi/dalpha along the polar great circle.
  derivs_wrt_alpha(theta, phi, dth, dph, N_angle, d_alpha);

  // Build the proper line-element integrand and integrate to obtain C_polar.
#pragma omp parallel for
  for (int i = 0; i < N_angle; i++) {
    const REAL qtt = sqrt_qtt[i] * sqrt_qtt[i];
    const REAL qpp = sqrt_qpp[i] * sqrt_qpp[i];
    // ds = sqrt( qtt * dth^2 + 2*qtp * dth*dph + qpp * dph^2 )
    integrand[i] = sqrt(qtt * (dth[i] * dth[i]) + 2.0 * qtp[i] * (dth[i] * dph[i]) + qpp * (dph[i] * dph[i]));
  }
  const REAL C_polar = integrate_over_alpha(integrand, N_angle, d_alpha);

  // Guard against invalid integrals and store outputs.
  if (!(C_equator > 0.0) || !isfinite(C_equator) || !isfinite(C_polar)) {
    commondata->error_flag = BHAHAHA_FAILURE;
    free(metric_data_gfs);
    return commondata->error_flag;
  }

  // Compute ratio and spin estimate; write results to macro-backed fields.
  const REAL ratio = C_polar / C_equator; // C_polar / C_equator wrt spin axis
  const REAL a_est = compute_spin(ratio); // compute_spin(circumf_ratio_polar_over_equator)

  commondata->bhahaha_diagnostics->BHAHAHA_OUT_CIRC_POLAR_FIELD   = C_polar;     // polar_circ
  commondata->bhahaha_diagnostics->BHAHAHA_OUT_CIRC_EQUATOR_FIELD = C_equator;   // equator_circ
  commondata->bhahaha_diagnostics->BHAHAHA_OUT_CIRC_RATIO_FIELD   = ratio;       // circumf_ratio_polar_over_equator
  commondata->bhahaha_diagnostics->BHAHAHA_OUT_SPIN_FIELD         = a_est;       // spin_a_along_spin_axis_from_circumf

  // Free the only heap allocation from this scope.
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

    results = doctest.testmod()

    if results.failed > 0:
        raise RuntimeError(
            f"Doctest failed: {results.failed} of {results.attempted} test(s)"
        )
    print(f"Doctest passed: All {results.attempted} test(s) passed")

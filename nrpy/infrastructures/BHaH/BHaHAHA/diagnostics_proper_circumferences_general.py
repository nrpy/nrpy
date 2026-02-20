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
  The diagonals are nonnegative "scale factors" along coordinate lines, so interpolating their square roots reduces
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

**Implementation note / justification (update):**
  The spherical grid is cell-centered in (theta,phi), so no poles are sampled. We take advantage of this smoothness by
  precomputing the *scalar* line-element integrands f_eq(theta,phi) and f_pol(theta,phi) on the i0=NGHOSTS slab,
  each built from the local 2-metric and a smooth tangent field for the corresponding great-circle family.
  During integration we then interpolate just this scalar along alpha midpoints and apply the same midpoint weights.

  Why this is preferable here:
   - Consistent accuracy: eliminates 2nd-order angle differencing while retaining high-order midpoint integration.
   - Fewer interpolations: 1 scalar interpolant instead of 3 fields (sqrt_qtt, sqrt_qpp, q_tp) per sample.
   - Robustness: no branch cuts or unwrapping are needed since we avoid finite differencing of angles.
   - Vectorization-friendly: the precompute step is a tight 2D loop; integration becomes a single scalar gather.

Authors: Wesley Inselman
         Zachariah B. Etienne
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


def register_CFunction_diagnostics_proper_circumferences_general(
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
    >>> env = register_CFunction_diagnostics_proper_circumferences_general()
    Setting up ExpansionFunctionThetaClass[Spherical]...
    Setting up reference_metric[Spherical]...
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
static void elliptic_E_and_K_integrals(const REAL k, REAL *restrict E, REAL *restrict K) {
  static const int N_sample_pts = 128; // Number of sample points for integration. Chosen for high precision.
  const REAL *restrict weights;        // Precomputed integration weights for accuracy in the midpoint method.
  int weight_stencil_size;             // Size of the weight stencil to correctly cycle through the weights.

  // Retrieve the integration weights based on the number of sample points (divisible by 8 -> 8th order).
  bah_diagnostics_integration_weights(N_sample_pts, N_sample_pts, &weights, &weight_stencil_size);

  const REAL a = 0.0;                            // Lower limit of integration (0 radians).
  const REAL b = M_PI / 2.0;                     // Upper limit of integration (pi/2 radians).
  const REAL h = (b - a) / ((REAL)N_sample_pts); // Step size for each subinterval based on sample points.

  REAL sum_E = 0.0; // Accumulator for the elliptic integral of the second kind E(k).
  REAL sum_K = 0.0; // Accumulator for the elliptic integral of the first kind K(k).

  // Parallelized loop to compute both integrals E(k) and K(k) using OpenMP.
  // #pragma omp parallel for reduction(+ : sum_E, sum_K) <- thread creation/destruction likely -> slower code here
  for (int i = 0; i < N_sample_pts; i++) {
    const REAL theta = a + ((REAL)i + 0.5) * h; // Compute the midpoint for the current subinterval.
    // Compute sin(theta). For optimization, we could use a lookup table if N_sample_pts is constant.
    const REAL sintheta = sin(theta);
    // Compute the integrands for the elliptic integrals of the second & first kinds (E(k) & K(k), respectively) at this midpoint.
    const REAL elliptic_E_integrand = sqrt(1.0 - k * sintheta * sintheta);
    const REAL elliptic_K_integrand = 1.0 / elliptic_E_integrand;
    // Update the running sums for E(k) & K(k), applying the corresponding weight.
    sum_E += weights[i % weight_stencil_size] * elliptic_E_integrand;
    sum_K += weights[i % weight_stencil_size] * elliptic_K_integrand;
  } // END LOOP over sample points to compute both integrals

  // Multiply by the step size to complete the integration and store the results.
  *E = sum_E * h; // Elliptic integral of the second kind.
  *K = sum_K * h; // Elliptic integral of the first kind.
} // END FUNCTION: elliptic_E_and_K_integrals

/**
 * Estimates the spin parameter magnitude for equilibrium black holes based on the circumference ratio C_r.
 *
 * Rationale & safeguards:
 *   * Start from an analytic approximation (Alcubierre et al., Eq. 5.3) to land near the root basin.
 *   * Refine with Newton–Raphson using elliptic integrals (Eq. 5.2), but clamp any out-of-range iterates
 *     into [0,1] and terminate if we step negative; this avoids excursions where the modulus or square roots
 *     would be undefined or numerically fragile.
 *
 * @param C_r The circumference ratio parameter used to estimate the spin.
 * @return    The estimated spin parameter. Returns -10.0 if C_r is out of valid bounds or if convergence fails.
 *
 */
static REAL compute_spin(const REAL C_r) {
  // Validate the input parameter. Return an error code if C_r exceeds the valid range.
  if (C_r > 1)
    return -10.0;

  // Turns out, this is a conservative spin estimate.
  REAL spin = 0.9;

  // Refine the initial guess using an analytical approximation based on Eq. 5.3 of Alcubierre et al arXiv:gr-qc/0411149.
  const REAL spin_sq = 1 - (2.55 * C_r - 1.55) * (2.55 * C_r - 1.55);
  if (spin_sq >= 0 && spin_sq < 1)
    spin = sqrt(spin_sq);

  const REAL rel_tolerance = 1e-7; // Desired relative tolerance for convergence.
  REAL rel_diff = 1e10;            // Initialize relative difference to a large value.
  const int max_its = 20;          // Maximum number of iterations to prevent infinite loops.
  int it = 0;                      // Iteration counter.

  // Iteratively refine the spin estimate until the relative difference is within tolerance or max iterations are reached.
  while (rel_diff > rel_tolerance && it < max_its) {
    const REAL x = spin;
    REAL E, K;

    // Compute the elliptic integrals E and K based on the current spin estimate.
    elliptic_E_and_K_integrals(-((x * x) / pow(1 + sqrt(1 - (x * x)), 2)), &E, &K);

    // Next complete a Newton-Raphson iteration to improve the spin estimate.
"""
    # Newton-Raphson is a bit more robust; also our initial guess is pretty good, so typically we need only a few iterations.
    prefunc += ccg.c_codegen(
        BHaH.BHaHAHA.area.spin_NewtonRaphson(), "const REAL x_np1", include_braces=False
    )
    prefunc += r"""
    if (x_np1 > 1.0) {
      // Adjust the spin estimate to remain within valid bounds.
      spin = 2.0 - x_np1;
    } else if (x_np1 < 0) {
      // Terminate iteration if the spin magnitude estimate becomes negative.
      it = max_its;
      break;
    } else {
      // Calculate the relative difference and update the spin estimate.
      rel_diff = fabs(x_np1 - x) / x;
      spin = x_np1;
    } // END spin adjustment to go back in-bounds
    it++;
  } // END WHILE: Refining spin estimate until convergence or maximum iterations

  // Assign spin=-10 if the Newton-Raphson did not converge within the allowed iterations.
  if (it >= max_its)
    spin = -10.0;

  return spin;
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
    } // END LOOP over inner boundary points
  } // END LOOP over gridfunctions
} // END apply_inner_bc_for_selected_gfs()

// ---------- small vector & math helpers (pure C) ----------
// These are intentionally tiny & inlined: they live in tight OpenMP loops,
// and we want predictable codegen and good autovectorization.
static inline REAL dot3(const REAL a[3], const REAL b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline void cross3(const REAL a[3], const REAL b[3], REAL c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
} // END FUNCTION cross3()
static inline REAL norm3(const REAL a[3]) { return sqrt(dot3(a,a)); }
static inline void normalize3(REAL a[3]) { const REAL n = norm3(a); if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; } }
// Build an orthonormal basis {e1, e2} orthogonal to s; avoid near-collinearity for numerical stability.
static inline void build_basis_from_s(const REAL s[3], REAL e1[3], REAL e2[3]) {
  REAL a[3] = {1.0, 0.0, 0.0};
  if (fabs(dot3(a,s)) > 0.9) { a[0]=0.0; a[1]=1.0; a[2]=0.0; }
  cross3(s, a, e1); normalize3(e1);
  cross3(s, e1, e2); normalize3(e2);
} // END FUNCTION build_basis_from_s()
// Cartesian unit vector -> spherical angles (theta, phi); inputs here are unit by construction.
static inline void cart_to_sph(const REAL r[3], REAL *theta, REAL *phi) {
  const REAL x=r[0], y=r[1], z=r[2];
  REAL zz = z; if (zz < -1.0) zz = -1.0; if (zz > 1.0) zz = 1.0;
  *theta = acos(zz);
  *phi   = atan2(y, x);
} // END FUNCTION cart_to_sph()

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
} // END FUNCTION integrate_over_alpha()
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
  (BHaHAHA/area.py::circumference_metric_roots). In addition, we precompute scalar integrands
  f_eq(theta,phi) and f_pol(theta,phi) using smooth great-circle tangent fields, so the integration phase
  only interpolates a single scalar per sample. This avoids finite differencing of angles entirely
  and aligns the derivative accuracy with the high-order midpoint quadrature.
"""
    cfunc_type = "int"
    name = "diagnostics_proper_circumferences_general"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  // which_gf in {0: sqrt(q_tt), 1: sqrt(q_pp), 2: q_tp, 3: f_eq(θ,φ), 4: f_pol(θ,φ)};
  // diagonals as sqrt for stability; f_eq/f_pol are scalar line-element integrands.
  const int NUM_DIAG_GFS = 5;
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  REAL *restrict metric_data_gfs;
  // Single heap allocation sized for a 2D (theta,phi) slab at i0 = NGHOSTS.
  BHAH_MALLOC(metric_data_gfs, Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(REAL));

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

    // Apply inner boundary conditions to q-metric gridfunctions sqrt(qtt), sqrt(qpp), and qtp:
    {
      const int which_gfs_1[3] = {0, 1, 2};
      apply_inner_bc_for_selected_gfs(&griddata[grid].bcstruct, metric_data_gfs, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                      which_gfs_1, 3);
    } // END application of inner boundary conditions
  } // END computation of q-metric root gridfunctions

  // Number of angular points to sample over 2 pi radians (controls resolution of great-circle integrals).
  const int N_angle = griddata[grid].params.Nxx2;
  // Uniform alpha step for midpoint samples in [-pi, pi).
  const REAL d_alpha = (M_PI - (-M_PI)) / ((REAL)N_angle);

  REAL dst_pts[N_angle][2];
  REAL theta[N_angle];
  REAL phi[N_angle];
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
  // Normal of the polar great-circle plane (spanned by s and e1):
  REAL nvec[3];
  cross3(s, e1, nvec);

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
    } // END LOOP over theta
  } // END LOOP over phi

  // Apply inner boundary conditions to the newly computed scalar integrands f_eq (which_gf=3)
  // and f_pol (which_gf=4), so their ghost zones are valid prior to interpolation.
  {
    const int which_gfs_2[2] = {3, 4};
    apply_inner_bc_for_selected_gfs(&griddata[grid].bcstruct, metric_data_gfs, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2,
                                    which_gfs_2, 2);
  } // END application of inner boundary conditions for f_eq/f_pol

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
  } // END LOOP over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed equator integrand f_eq onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (REAL *restrict *)(*coords), src_feq, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END error check
  } // END interpolation step across equatorial circumference

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
  } // END LOOP over alpha: cell-centered sampling 2 pi across N_angle points

  {
    // Interpolate precomputed polar integrand f_pol onto the alpha-midpoints.
    int err =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, griddata[grid].params.dxx1, griddata[grid].params.dxx2, Nxx_plus_2NGHOSTS1,
                                                       Nxx_plus_2NGHOSTS2, (REAL *restrict *)(*coords), src_fpol, N_angle, dst_pts, integrand);
    if (err != BHAHAHA_SUCCESS) {
      free(metric_data_gfs);
      return err;
    } // END error check
  } // END interpolation step across polar circumference

  const REAL C_polar = integrate_over_alpha(integrand, N_angle, d_alpha);

  // Compute ratio and spin estimate; write results to macro-backed fields.
  const REAL ratio = C_polar / C_equator; // C_polar / C_equator wrt spin axis
  const REAL a_est = compute_spin(ratio); // compute_spin(circumf_ratio_polar_over_equator)

  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_POLAR = C_polar;     // polar_circ
  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_EQUATOR = C_equator; // equator_circ
  commondata->bhahaha_diagnostics->BHAHAHA_CIRC_GENERAL_SPIN = a_est;        // spin_a_along_spin_axis_from_circumf

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

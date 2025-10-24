# BHaHAHA/diagnostics_proper_circumferences.py
"""
Register the C function for computing polar & equatorial circumferences of the horizon.
Needs: coordinate centroid of the horizon.

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
    Build the three quantities needed to evaluate the proper-circumference arclength integrand.

    I.e., these quantities:
    sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi},

    where q_{ab} is the induced 2-metric on r = h(theta,phi).

    These are pure functions of the ambient metric components gamma_{ij}, the radius h(theta,phi),
    and its angular derivatives h_{,theta}, h_{,phi}:

    q_{theta theta} = gamma_{rr} h_{,theta}^2 + 2 gamma_{r theta} h_{,theta} + gamma_{theta theta}
    q_{phi phi}   = gamma_{rr} h_{,phi}^2  + 2 gamma_{r phi}  h_{,phi}  + gamma_{phi phi}
    q_{theta phi}  = gamma_{rr} h_{,theta} h_{,phi} + gamma_{r theta} h_{,phi} + gamma_{r phi} h_{,theta} + gamma_{theta phi}

    The square-roots are applied only to the diagonal elements, because the line element along a curve
    requires q_{theta theta}, q_{phi phi}, and q_{theta phi} separately.
    """
    Th = ExpansionFunctionTheta["Spherical"]
    h = sp.Symbol("hh", real=True)
    h_dD = ixp.declarerank1("hh_dD")

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

    sqrt_qtt = sp.sqrt(qtt.replace(sp.sympify("xx0"), h))
    sqrt_qpp = sp.sqrt(qpp.replace(sp.sympify("xx0"), h))
    qtp = qtp.replace(sp.sympify("xx0"), h)
    return sqrt_qtt, sqrt_qpp, qtp


def register_CFunction_diagnostics_proper_circumferences(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the C function for computing minimum, maximum, and mean radius *from the coordinate centroid* of the horizon.
    Needs: coordinate centroid of the horizon.

    :param enable_fd_functions: Whether to enable finite difference functions, defaults to True.
    :return: An NRPyEnv_type object if registration is successful, otherwise None.

    DocTests:
    >>> import nrpy.grid as gri
    >>> _ = gri.register_gridfunctions("hh")[0]
    >>> env = register_CFunction_diagnostics_proper_circumferences()
    Setting up reference_metric[Spherical]...
    Setting up ExpansionFunctionThetaClass[Spherical]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = """
/**
 * Computes the complete elliptic integrals of the second kind E(k) and the first kind K(k)
 * using 8th-order midpoint integration.
 *
 * @param k - The elliptic modulus parameter (0 <= k <= 1).
 * @param E - Pointer to store the result for the elliptic integral of the second kind E(k).
 * @param K - Pointer to store the result for the elliptic integral of the first kind K(k).
 *
 * This function uses a midpoint integration method with a fixed number of sample points (128)
 * to ensure high accuracy, leveraging precomputed weights for efficiency.
 *
 * @note The function is parallelized with OpenMP to improve performance on large arrays.
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
#pragma omp parallel for reduction(+ : sum_E, sum_K)
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
  } // END PARALLEL FOR: Loop through sample points to compute both integrals

  // Multiply by the step size to complete the integration and store the results.
  *E = sum_E * h; // Elliptic integral of the second kind.
  *K = sum_K * h; // Elliptic integral of the first kind.
} // END FUNCTION: elliptic_E_and_K_integrals

/**
 * Estimates the spin parameter magnitude for equilibrium black holes based on the circumference ratio C_r.
 *
 * @param C_r The circumference ratio parameter used to estimate the spin.
 * @return    The estimated spin parameter. Returns -10.0 if C_r is out of valid bounds or if convergence fails.
 *
 * This function implements an iterative method to estimate the spin parameter using
 * elliptic integrals. It starts with an initial guess and refines it to achieve a desired
 * relative tolerance. The method is based on Eq. 5.2 of Alcubierre et al. (arXiv:gr-qc/0411149).
 */
REAL compute_spin(const REAL C_r) {
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
    }
    it++;
  } // END WHILE: Refining spin estimate until convergence or maximum iterations

  // Assign spin=-10 if the Newton-Raphson did not converge within the allowed iterations.
  if (it >= max_its) {
    spin = -10.0;
  }
  return spin;
}

#ifndef BHAHAHA_FAILURE
#define BHAHAHA_FAILURE (-1)
#endif

/* ---------- small vector & math helpers (pure C) ---------- */
static inline REAL dot3(const REAL a[3], const REAL b[3]) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]; }
static inline void cross3(const REAL a[3], const REAL b[3], REAL c[3]) {
  c[0] = a[1]*b[2] - a[2]*b[1];
  c[1] = a[2]*b[0] - a[0]*b[2];
  c[2] = a[0]*b[1] - a[1]*b[0];
}
static inline REAL norm3(const REAL a[3]) { return sqrt(dot3(a,a)); }
static inline void normalize3(REAL a[3]) { const REAL n = norm3(a); if (n > 0.0) { a[0]/=n; a[1]/=n; a[2]/=n; } }
static inline REAL wrap_dphi(REAL dphi) { while (dphi > M_PI) dphi -= 2.0*M_PI; while (dphi <= -M_PI) dphi += 2.0*M_PI; return dphi; }
static inline void build_basis_from_s(const REAL s[3], REAL e1[3], REAL e2[3]) {
  REAL a[3] = {1.0, 0.0, 0.0};
  if (fabs(dot3(a,s)) > 0.9) { a[0]=0.0; a[1]=1.0; a[2]=0.0; }
  cross3(s, a, e1); normalize3(e1);
  cross3(s, e1, e2); normalize3(e2);
}
static inline void cart_to_sph(const REAL r[3], REAL *theta, REAL *phi) {
  const REAL x=r[0], y=r[1], z=r[2];
  REAL zz = z; if (zz < -1.0) zz = -1.0; if (zz > 1.0) zz = 1.0;
  *theta = acos(zz);
  *phi   = atan2(y, x);
}

/* Compute dtheta/dalpha and dphi/dalpha with central differences.
 * phi(alpha) is unwrapped first to avoid the branch cut at +/- pi: without unwrapping,
 * differences across the cut would produce O(2*pi) jumps that corrupt the line element. */
static inline void derivs_wrt_alpha(const REAL *theta, const REAL *phi,
                                    REAL *dtheta_dalpha, REAL *dphi_dalpha,
                                    int N, REAL d_alpha) {
  REAL *phi_unwrap = (REAL *)malloc((size_t)N * sizeof(REAL));
  if (phi_unwrap == NULL) {
    for (int i = 0; i < N; i++) {
      const int im = (i - 1 + N) % N, ip = (i + 1) % N;
      dtheta_dalpha[i] = (theta[ip] - theta[im]) / (2.0*d_alpha);
      dphi_dalpha[i]   = wrap_dphi(phi[ip] - phi[im]) / (2.0*d_alpha);
    }
    return;
  }
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

static inline REAL integrate_over_alpha(const REAL *vals, int N_angle, REAL d_alpha) {
  const REAL *restrict weights;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);
  REAL sum = 0.0;
#pragma omp parallel for reduction(+ : sum)
  for (int i = 0; i < N_angle; i++) sum += vals[i] * weights[i % weight_stencil_size];
  return sum * d_alpha;
}

/* DRY helper: interpolate sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi} to dst_pts[ (theta,phi) ] */
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

@param commondata - Pointer to common data structure containing shared parameters and settings.
@param griddata - Pointer to grid data structures for each grid, containing parameters and gridfunctions.
@return - Status code indicating success or type of error (e.g., BHAHAHA_SUCCESS or INITIAL_DATA_MALLOC_ERROR).
@note - This function uses OpenMP for parallel loops and performs interpolation and integration over grid data.

We precompute on the (theta,phi) grid at fixed i0=NGHOSTS the following three quantities:
  sqrt(q_{theta theta}), sqrt(q_{phi phi}), and q_{theta phi},
generated via NRPy's SymPy expressions (BHaHAHA/area.py::circumference_metric_roots).

Then, along great circles defined by the unit spin axis s^i, we interpolate these quantities,
build the proper line element sqrt(q_tt (dtheta/dalpha)^2 + 2 q_tp (dtheta/dalpha)(dphi/dalpha) + q_pp (dphi/dalpha)^2),
and integrate over alpha in [-pi,pi) to obtain:
  - equator_circ: great circle orthogonal to s^i,
  - polar_circ:  great circle containing s^i.

We finally store the ratio C_polar/C_equator and an estimated spin a/M.
"""
    cfunc_type = "int"
    name = "diagnostics_proper_circumferences"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  const int NUM_DIAG_GFS = 3; // which_gf in {0: sqrt(q_tt), 1: sqrt(q_pp), 2: q_tp}
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  REAL *restrict metric_data_gfs;
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
    const int i0 = NGHOSTS; // Fixed index for radial coordinate (r).

    // Loop over angular grid points (theta and phi) to compute:
    // 1. sqrt(q_{theta theta}), stored to metric_data_gfs[IDX4(0,...)], and
    // 2. sqrt(q_{phi phi}), stored to metric_data_gfs[IDX4(1,...)] at each point (theta, phi).
    // 3. q_{theta phi}, stored to metric_data_gfs[IDX4(2,...)] at each point (theta, phi).
    // Notice we do clever indexing to ensure this 2D computation stays within the memory bounds of (3D) metric_data_gfs.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED REAL xx2 = xx[2][i2]; // Phi coordinate at index i2.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED REAL xx1 = xx[1][i1]; // Theta coordinate at index i1.
"""
    # Generate the three outputs from the induced 2-metric roots:
    body += ccg.c_codegen(
        list(BHaH.BHaHAHA.area.circumference_metric_roots()),
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

          // Apply boundary condition if at the radial interior point (i0 == NGHOSTS).
          if (dst_i0 == NGHOSTS) {
            metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(dst_i1, dst_i2)] = metric_data_gfs[IDX4pt(which_gf, 0) + IDX2(src_i1, src_i2)];
          }
        } // END LOOP over inner boundary points
      } // END LOOP over gridfunctions
    } // END application of inner boundary conditions
  } // END computation of q-metric root gridfunctions

  // Number of angular points to sample over 2 pi radians.
  const int N_angle = griddata[grid].params.Nxx2;
  // Compute the angular increment for sampling over 2 pi radians, in alpha parameter.
  const REAL d_alpha = (M_PI - (-M_PI)) / ((REAL)N_angle);

  // C99 variable-length arrays on the stack (no heap allocations here).
  REAL dst_pts[N_angle][2];
  REAL theta[N_angle];
  REAL phi[N_angle];
  REAL dth[N_angle];
  REAL dph[N_angle];
  REAL sqrt_qtt[N_angle];
  REAL sqrt_qpp[N_angle];
  REAL qtp[N_angle];
  REAL integrand[N_angle];

  // Next: normalize spin axis and write it back through the macro-backed fields.
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

  // Convenience base pointers for interpolation sources:
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

  // Compute dtheta/dalpha and dphi/dalpha along the equator, with phi unwrapped to avoid branch-cut artifacts.
  derivs_wrt_alpha(theta, phi, dth, dph, N_angle, d_alpha);

  // Build the proper line-element integrand and integrate to obtain C_equator.
#pragma omp parallel for
  for (int i = 0; i < N_angle; i++) {
    const REAL qtt = sqrt_qtt[i] * sqrt_qtt[i];
    const REAL qpp = sqrt_qpp[i] * sqrt_qpp[i];
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

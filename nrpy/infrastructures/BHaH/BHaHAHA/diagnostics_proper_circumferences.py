# BHaHAHA/diagnostics_proper_circumferences.py
"""
Register the C function for computing polar & equatorial circumferences of the horizon.
Needs: coordinate centroid of the horizon.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.infrastructures import BHaH


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
"""
    desc = """
Computes proper circumferences along the equator and polar directions for apparent horizon diagnostics.

@param commondata - Pointer to common data structure containing shared parameters and settings.
@param griddata - Pointer to grid data structures for each grid, containing parameters and gridfunctions.
@return - Status code indicating success or type of error (e.g., BHAHAHA_SUCCESS or INITIAL_DATA_MALLOC_ERROR).
@note - This function uses OpenMP for parallel loops and performs interpolation and integration over grid data.
"""
    cfunc_type = "int"
    name = "diagnostics_proper_circumferences"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )
    body = r"""
  const int NUM_DIAG_GFS = 2;
  const int grid = 0;
  // Extract grid dimensions, including ghost zones, for each coordinate direction. Needed for IDX4() macro.
  const int Nxx_plus_2NGHOSTS0 = griddata[grid].params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata[grid].params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata[grid].params.Nxx_plus_2NGHOSTS2;
  const int NUM_THETA = Nxx_plus_2NGHOSTS1; // Needed for IDX2() macro.

  REAL *restrict metric_data_gfs;
  BHAH_MALLOC(metric_data_gfs, Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2 * NUM_DIAG_GFS * sizeof(REAL));

  // Compute the line element in the phi direction (sqrt(q_{phi phi})) across the entire grid.
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
    // Notice we do clever indexing to ensure this 2D computation stays within the memory bounds of (3D) metric_data_gfs.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      const MAYBE_UNUSED REAL xx2 = xx[2][i2]; // Phi coordinate at index i2.
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        const MAYBE_UNUSED REAL xx1 = xx[1][i1]; // Theta coordinate at index i1.
"""
    body += ccg.c_codegen(
        [
            BHaH.BHaHAHA.area.circumferential_arclength(direction="theta"),
            BHaH.BHaHAHA.area.circumferential_arclength(direction="phi"),
        ],
        [
            "metric_data_gfs[IDX4pt(0, 0) + IDX2(i1, i2)]",
            "metric_data_gfs[IDX4pt(1, 0) + IDX2(i1, i2)]",
        ],
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )
    body += r"""
      } // END LOOP over i1 (theta)
    } // END LOOP over i2 (phi)

    // Apply inner boundary conditions to the computed sqrt(q_{phi phi}) gridfunction.
    {
      bc_struct *restrict bcstruct = &griddata[grid].bcstruct; // Retrieve boundary condition structure for the grid.

      // Unpack bc_info from bcstruct.
      const bc_info_struct *bc_info = &bcstruct->bc_info;

      // Apply boundary conditions at inner boundary points for the selected gridfunctions.
#pragma omp parallel for collapse(2)
      for (int which_gf = 0; which_gf < 2; which_gf++) {
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
  } // END computation of line element gridfunction

  // Grid spacings in theta and phi directions.
  const REAL dxx1 = griddata[grid].params.dxx1;
  const REAL dxx2 = griddata[grid].params.dxx2;
  // Number of angular points to sample over 2 pi radians.
  const int N_angle = griddata[grid].params.Nxx2;
  // Allocate arrays for destination points (theta, phi) and circumference values.
  REAL(*dst_pts)[2] = malloc(N_angle * sizeof(*dst_pts));
  REAL *restrict circumference = malloc(N_angle * sizeof(REAL));
  // Check for successful memory allocation.
  if (dst_pts == NULL || circumference == NULL) {
    if (dst_pts != NULL)
      free(dst_pts);
    if (circumference != NULL)
      free(circumference);
    commondata->error_flag = DIAG_PROPER_CIRCUM_MALLOC_ERROR; // Return error code if allocation fails.
    return commondata->error_flag;
  }

  // Compute the angular increment for sampling over 2 pi radians, whether it be in theta (xz & yz planes) or phi (xy-plane).
  const REAL d_angle = (M_PI - (-M_PI)) / ((REAL)N_angle);

  // Equatorial (xy-plane) circumference first
  {
    // Initialize destination points along the equator (theta = pi/2) for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      dst_pts[i2][0] = M_PI / 2;                                   // Equator: theta = pi/2.
      dst_pts[i2][1] = -M_PI + ((REAL)i2 + (1.0 / 2.0)) * d_angle; // Equator: phi = [-pi, pi].
    } // END LOOP over phi angles

    // Interpolate sqrt(q_{phi phi}) values onto the equator points to compute the circumference;
    //   note that sqrt(q_{phi phi}) is stored in metric_data_gfs[IDX4(1,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(1, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total xy-plane circumference by integrating over the sampled points.
    REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP over ic
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->xy_plane_circumference = sum_circumference * d_angle;
  } // END xy-plane circumference

  // Polar (xz-plane) circumference next
  {
    // Initialize destination points along the xz-plane for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      if (i2 < N_angle / 2) {
        // First half: Theta from 0 to pi, phi = 0
        dst_pts[i2][0] = ((REAL)i2 + 0.5) * (M_PI / ((REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = 0.0;
      } else {
        // Second half: Theta from pi back to 0, phi = -pi
        dst_pts[i2][0] = ((REAL)(N_angle - i2) - 0.5) * (M_PI / ((REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = -M_PI; // phi spans from [-pi, pi), so instead of interpolating at phi=pi, must interpolate at phi=-pi.
      } // END IF theta is going from 0 to pi or vice-versa.
    } // END LOOP over angle

    // Interpolate sqrt(q_{theta theta}) values onto the polar (xz-plane) points to compute the circumference;
    //   note that sqrt(q_{theta theta}) is stored in metric_data_gfs[IDX4(0,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(0, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total xz-plane circumference by integrating over the sampled points.
    REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP over ic
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->xz_plane_circumference = sum_circumference * d_angle;
  } // END xz-plane circumference

  // Polar (yz-plane) circumference next
  {
    // Initialize destination points along the yz-plane for interpolation.
#pragma omp parallel for
    for (int i2 = 0; i2 < N_angle; i2++) {
      if (i2 < N_angle / 2) {
        // First half: Theta from 0 to pi, phi = pi/2
        dst_pts[i2][0] = ((REAL)i2 + 0.5) * (M_PI / ((REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = M_PI / 2.0;
      } else {
        // Second half: Theta from pi back to 0, phi = -pi/2
        dst_pts[i2][0] = ((REAL)(N_angle - i2) - 0.5) * (M_PI / ((REAL)(N_angle) / 2.0));
        dst_pts[i2][1] = -M_PI / 2.0;
      } // END IF theta is going from 0 to pi or vice-versa.
    } // END LOOP over angle

    // Interpolate sqrt(q_{theta theta}) values onto the polar (yz-plane) points to compute the circumference;
    //   note that sqrt(q_{theta theta}) is stored in metric_data_gfs[IDX4(0,...)]
    const int error =
        bah_interpolation_2d_general__uniform_src_grid(NinterpGHOSTS, dxx1, dxx2, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, griddata[grid].xx,
                                                       &metric_data_gfs[IDX4pt(0, 0)], N_angle, dst_pts, circumference);
    if (error != BHAHAHA_SUCCESS)
      return error;

    // Retrieve integration weights for numerical integration over the sampled points.
    const REAL *restrict weights;
    int weight_stencil_size;
    bah_diagnostics_integration_weights(N_angle, N_angle, &weights, &weight_stencil_size);

    // Compute the total yz-plane circumference by integrating over the sampled points.
    REAL sum_circumference = 0.0;
#pragma omp parallel for reduction(+ : sum_circumference)
    for (int ic = 0; ic < N_angle; ic++) {
      const REAL weight = weights[ic % weight_stencil_size]; // Integration weight for this point.
      sum_circumference += circumference[ic] * weight;
    } // END LOOP over ic
    // Multiply the sum by d[angle]
    commondata->bhahaha_diagnostics->yz_plane_circumference = sum_circumference * d_angle;
  } // END yz-plane circumference

  // Next estimate spin parameter magnitudes, valid for equilibrium BHs only.
  //   Based on Eq 5.2 of Alcubierre et al arXiv:gr-qc/0411149.
  {
    REAL C_xz_yz = commondata->bhahaha_diagnostics->xz_plane_circumference / commondata->bhahaha_diagnostics->yz_plane_circumference;
    REAL C_xy_yz = commondata->bhahaha_diagnostics->xy_plane_circumference / commondata->bhahaha_diagnostics->yz_plane_circumference;
    REAL C_yz_xz = commondata->bhahaha_diagnostics->yz_plane_circumference / commondata->bhahaha_diagnostics->xz_plane_circumference;
    REAL C_xy_xz = commondata->bhahaha_diagnostics->xy_plane_circumference / commondata->bhahaha_diagnostics->xz_plane_circumference;
    REAL C_xz_xy = commondata->bhahaha_diagnostics->xz_plane_circumference / commondata->bhahaha_diagnostics->xy_plane_circumference;
    REAL C_yz_xy = commondata->bhahaha_diagnostics->yz_plane_circumference / commondata->bhahaha_diagnostics->xy_plane_circumference;

    // Compute spin estimates for each direction based on different circumference ratios
    commondata->bhahaha_diagnostics->spin_a_x_from_xz_over_yz_prop_circumfs = compute_spin(C_xz_yz);
    commondata->bhahaha_diagnostics->spin_a_x_from_xy_over_yz_prop_circumfs = compute_spin(C_xy_yz);
    commondata->bhahaha_diagnostics->spin_a_y_from_yz_over_xz_prop_circumfs = compute_spin(C_yz_xz);
    commondata->bhahaha_diagnostics->spin_a_y_from_xy_over_xz_prop_circumfs = compute_spin(C_xy_xz);
    commondata->bhahaha_diagnostics->spin_a_z_from_xz_over_xy_prop_circumfs = compute_spin(C_xz_xy);
    commondata->bhahaha_diagnostics->spin_a_z_from_yz_over_xy_prop_circumfs = compute_spin(C_yz_xy);
  }
  // Free allocated memory for destination points, circumference values, and metric_data_gfs.
  free(dst_pts);
  free(circumference);
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

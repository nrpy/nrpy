#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation_lagrange_uniform.h"
#include <stdint.h>

/**
 * Map one raw extended-grid stencil node back into the stored payload.
 *
 * The interpolation coefficients are always built on the raw uniform stencil.
 * This helper only recovers the in-domain payload index for one raw node by
 * applying the dataset-specific spherical reflection rules directly. The key
 * idea is that the polynomial nodes stay on the mathematically uniform
 * extended stencil, while only the payload reads are reflected back into the
 * stored two-plane data.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param[in] r_ext Raw stencil radial coordinate, possibly outside the interior.
 * @param[in] theta_ext Raw stencil polar coordinate, possibly outside the interior.
 * @param[in] i2_base Logical phi index for the selected stored reference plane.
 * @param[out] i0i1i2_map In-domain logical payload index recovered by reflection.
 * @return Status code indicating whether the mapped node lands inside the payload.
 */
static int azimuthal_symmetry_spatial_lagrange_map_extended_node(const params_struct *restrict params, const REAL r_ext, const REAL theta_ext,
                                                                 const int i2_base, int i0i1i2_map[3]) {
  const int payload_i0_start = NGHOSTS;
  const int payload_i0_end = NGHOSTS + params->Nxx0;
  const int payload_i1_start = NGHOSTS;
  const int payload_i1_end = NGHOSTS + params->Nxx1;
  const int payload_i2_start = NGHOSTS;
  const int payload_i2_end = NGHOSTS + params->Nxx2;
  const REAL theta_tol = (REAL)1.0e-12;
  const REAL r_endpoint_tol = (REAL)(1.0e-12 * fmax(1.0, fabs((double)params->dxx0)));
  const REAL theta_endpoint_tol = (REAL)(1.0e-12 * fmax(1.0, fabs((double)params->dxx1)));
  const REAL theta_max_supported = (REAL)(params->xxmin1 + (((payload_i1_end - 1 - NGHOSTS) + 0.5) * params->dxx1));
  const REAL r_max_supported = (REAL)(params->xxmin0 + (((payload_i0_end - 1 - NGHOSTS) + 0.5) * params->dxx0));
  REAL r_map = r_ext;
  REAL theta_map = theta_ext;
  REAL theta_index_map = theta_ext;
  int phi_toggle = 0;
  int i0_map;
  int i1_map;
  int i2_map;

  // Step 1: Reflect negative-radius raw nodes through the origin.
  if (r_map < (REAL)0.0) {
    r_map = -r_map;
    theta_map = (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - theta_map;
    phi_toggle ^= 1;
  } // END IF: raw stencil node crossed the inner radial boundary

  // Step 2: Fold theta back into the physical [0, pi] range.
  if (theta_map < (REAL)0.0) {
    theta_map = -theta_map;
    phi_toggle ^= 1;
  } else if (theta_map > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI) {
    theta_map = (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI) - theta_map;
    phi_toggle ^= 1;
  } // END IF: raw stencil node crossed a polar boundary

  // Step 3: Only accept theta values that now land inside the physical range.
  if (theta_map < -theta_tol || theta_map > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + theta_tol) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: reflected theta remained materially out of range
  if (theta_map < (REAL)0.0) {
    theta_map = (REAL)0.0;
  } else if (theta_map > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI) {
    theta_map = (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI;
  } // END IF: theta exceeded the physical range only by roundoff

  // Step 4: Repair only small endpoint roundoff before converting to logical indices.
  if (r_map > r_max_supported && fabs((double)(r_map - r_max_supported)) <= r_endpoint_tol) {
    r_map = r_max_supported;
  } // END IF: radial coordinate landed on the last stored cell center to roundoff
  theta_index_map = theta_map;
  if (theta_index_map > theta_max_supported && theta_index_map <= theta_max_supported + (REAL)(0.5 * params->dxx1) + theta_endpoint_tol) {
    theta_index_map = theta_max_supported;
  } // END IF: theta landed in the last physical half-cell and should read the endpoint payload cell

  // Step 5: Convert the reflected coordinates and phi toggle to payload indices.
  i0_map = (int)((r_map - params->xxmin0) / params->dxx0 + (REAL)NGHOSTS);
  i1_map = (int)((theta_index_map - params->xxmin1) / params->dxx1 + (REAL)NGHOSTS);
  if (phi_toggle == 0) {
    i2_map = i2_base;
  } else {
    i2_map = payload_i2_end - 1 - (i2_base - payload_i2_start);
  } // END IF: apply the pi shift through the stored two-plane layout

  // Step 6: Confirm that the remapped logical index lies inside the stored payload.
  if (i0_map < payload_i0_start || i0_map >= payload_i0_end || i1_map < payload_i1_start || i1_map >= payload_i1_end || i2_map < payload_i2_start ||
      i2_map >= payload_i2_end) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: reflected stencil node could not be supported by the payload
  i0i1i2_map[0] = i0_map;
  i0i1i2_map[1] = i1_map;
  i0i1i2_map[2] = i2_map;
  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_map_extended_node

/**
 * Accumulate one serialized tensor record onto the common stored phi plane.
 *
 * The remap logic guarantees that the payload contributes either from the same
 * stored phi plane or from the opposite plane separated by pi. The same-plane
 * case is an identity, while the opposite-plane case is diagonal in Cartesian
 * components and therefore reduces to serialized sign flips.
 *
 * @param weight Weighted 2D Lagrange coefficient for this stencil node.
 * @param apply_pi_shift Whether the node came from the opposite stored plane.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[out] tensor_ref Accumulator on the selected common stored phi plane.
 */
static void azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_to_common_plane(
    const REAL weight, const int apply_pi_shift, const double *restrict tensor_record,
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT]) {
  const double *restrict gamma_record = tensor_record + AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;
  static const REAL g4dd_pi_sign[10] = {1.0, -1.0, -1.0, 1.0, 1.0, 1.0, -1.0, 1.0, -1.0, 1.0};
  static const REAL gamma4udd_pi_sign[40] = {1.0,  -1.0, -1.0, 1.0,  1.0,  1.0,  -1.0, 1.0, -1.0, 1.0,  -1.0, 1.0,  1.0, -1.0,
                                             -1.0, -1.0, 1.0,  -1.0, 1.0,  -1.0, -1.0, 1.0, 1.0,  -1.0, -1.0, -1.0, 1.0, -1.0,
                                             1.0,  -1.0, 1.0,  -1.0, -1.0, 1.0,  1.0,  1.0, -1.0, 1.0,  -1.0, 1.0};

  // Step 1: Same-plane nodes need only weighted accumulation.
  if (apply_pi_shift == 0) {
    for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++) {
      tensor_ref[comp] += weight * (REAL)tensor_record[comp];
    } // END LOOP: for comp over serialized metric payload components
    for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT; comp++) {
      tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] += weight * (REAL)gamma_record[comp];
    } // END LOOP: for comp over serialized Christoffel payload components
  } else {
    // Step 2: Opposite-plane nodes differ only by the pi-rotation sign pattern.
    for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++) {
      tensor_ref[comp] += weight * g4dd_pi_sign[comp] * (REAL)tensor_record[comp];
    } // END LOOP: for comp over pi-rotated serialized metric components
    for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT; comp++) {
      tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] += weight * gamma4udd_pi_sign[comp] * (REAL)gamma_record[comp];
    } // END LOOP: for comp over pi-rotated serialized Christoffel components
  } // END IF: choose same-plane or opposite-plane accumulation path
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_to_common_plane

/**
 * Rotate Cartesian-basis metric and Christoffel components about the z axis.
 *
 * The payload stores tensors on one of exactly two native reference azimuthal
 * planes. This helper does not interpolate in `phi`; instead, axisymmetry
 * recovers the target azimuth by an active spatial rotation through
 * `delta_phi` after the `(r, theta)` interpolation has completed. This keeps
 * the interpolation strictly two-dimensional in the stored data while still
 * returning tensors at the requested azimuth, embedded in 4D so that:
 *
 *   x'^i = Lambda^i_j x^j
 *   g'_{mu nu} = (Lambda^T)^alpha_mu (Lambda^T)^beta_nu g_{alpha beta}
 *   Gamma'^alpha_{mu nu} =
 *       Lambda^alpha_beta (Lambda^T)^gamma_mu (Lambda^T)^delta_nu
 *       Gamma^beta_{gamma delta}
 *
 * The time and z directions are unchanged, while the x-y block is the usual
 * planar rotation matrix.
 *
 * @param delta_phi Active rotation angle from stored plane to target azimuth.
 * @param[in] g4dd_ref Upper-triangular metric components on the stored plane.
 * @param[in] gamma_ref Serialized Christoffel components on the stored plane.
 * @param[out] g4dd_rot Rotated upper-triangular metric components.
 * @param[out] gamma_rot Rotated serialized Christoffel components.
 */
static void azimuthal_symmetry_spatial_lagrange_rotate_about_z(const REAL delta_phi,
                                                               const REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
                                                               const REAL gamma_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT],
                                                               REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
                                                               REAL gamma_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT]) {
  const REAL cos_delta = cos(delta_phi);
  const REAL sin_delta = sin(delta_phi);

  // Step 1: Rotate only the serialized outputs required by the caller.
  g4dd_rot[0] = g4dd_ref[0];
  g4dd_rot[1] = cos_delta * g4dd_ref[1] - sin_delta * g4dd_ref[2];
  g4dd_rot[2] = sin_delta * g4dd_ref[1] + cos_delta * g4dd_ref[2];
  g4dd_rot[3] = g4dd_ref[3];
  g4dd_rot[4] = cos_delta * cos_delta * g4dd_ref[4] - 2.0 * cos_delta * sin_delta * g4dd_ref[5] + sin_delta * sin_delta * g4dd_ref[7];
  g4dd_rot[5] = cos_delta * sin_delta * g4dd_ref[4] - sin_delta * sin_delta * g4dd_ref[5] + cos_delta * cos_delta * g4dd_ref[5] -
                cos_delta * sin_delta * g4dd_ref[7];
  g4dd_rot[6] = cos_delta * g4dd_ref[6] - sin_delta * g4dd_ref[8];
  g4dd_rot[7] = sin_delta * sin_delta * g4dd_ref[4] + 2.0 * cos_delta * sin_delta * g4dd_ref[5] + cos_delta * cos_delta * g4dd_ref[7];
  g4dd_rot[8] = sin_delta * g4dd_ref[6] + cos_delta * g4dd_ref[8];
  g4dd_rot[9] = g4dd_ref[9];
  gamma_rot[0] = gamma_ref[0];
  gamma_rot[1] = cos_delta * gamma_ref[1] - sin_delta * gamma_ref[2];
  gamma_rot[2] = sin_delta * gamma_ref[1] + cos_delta * gamma_ref[2];
  gamma_rot[3] = gamma_ref[3];
  gamma_rot[4] = cos_delta * cos_delta * gamma_ref[4] - 2.0 * cos_delta * sin_delta * gamma_ref[5] + sin_delta * sin_delta * gamma_ref[7];
  gamma_rot[5] = cos_delta * sin_delta * gamma_ref[4] - sin_delta * sin_delta * gamma_ref[5] + cos_delta * cos_delta * gamma_ref[5] -
                 cos_delta * sin_delta * gamma_ref[7];
  gamma_rot[6] = cos_delta * gamma_ref[6] - sin_delta * gamma_ref[8];
  gamma_rot[7] = sin_delta * sin_delta * gamma_ref[4] + 2.0 * cos_delta * sin_delta * gamma_ref[5] + cos_delta * cos_delta * gamma_ref[7];
  gamma_rot[8] = sin_delta * gamma_ref[6] + cos_delta * gamma_ref[8];
  gamma_rot[9] = gamma_ref[9];
  gamma_rot[10] = cos_delta * gamma_ref[10] - sin_delta * gamma_ref[20];
  gamma_rot[11] = cos_delta * cos_delta * gamma_ref[11] - cos_delta * sin_delta * gamma_ref[12] - cos_delta * sin_delta * gamma_ref[21] +
                  sin_delta * sin_delta * gamma_ref[22];
  gamma_rot[12] = cos_delta * sin_delta * gamma_ref[11] + cos_delta * cos_delta * gamma_ref[12] - sin_delta * sin_delta * gamma_ref[21] -
                  cos_delta * sin_delta * gamma_ref[22];
  gamma_rot[13] = cos_delta * gamma_ref[13] - sin_delta * gamma_ref[23];
  gamma_rot[14] = cos_delta * cos_delta * cos_delta * gamma_ref[14] - 2.0 * cos_delta * cos_delta * sin_delta * gamma_ref[15] +
                  cos_delta * sin_delta * sin_delta * gamma_ref[17] - cos_delta * cos_delta * sin_delta * gamma_ref[24] +
                  2.0 * cos_delta * sin_delta * sin_delta * gamma_ref[25] - sin_delta * sin_delta * sin_delta * gamma_ref[27];
  gamma_rot[15] = cos_delta * cos_delta * sin_delta * gamma_ref[14] - cos_delta * sin_delta * sin_delta * gamma_ref[15] +
                  cos_delta * cos_delta * cos_delta * gamma_ref[15] - cos_delta * cos_delta * sin_delta * gamma_ref[17] -
                  cos_delta * sin_delta * sin_delta * gamma_ref[24] + sin_delta * sin_delta * sin_delta * gamma_ref[25] -
                  cos_delta * cos_delta * sin_delta * gamma_ref[25] + cos_delta * sin_delta * sin_delta * gamma_ref[27];
  gamma_rot[16] = cos_delta * cos_delta * gamma_ref[16] - cos_delta * sin_delta * gamma_ref[18] - cos_delta * sin_delta * gamma_ref[26] +
                  sin_delta * sin_delta * gamma_ref[28];
  gamma_rot[17] = cos_delta * sin_delta * sin_delta * gamma_ref[14] + 2.0 * cos_delta * cos_delta * sin_delta * gamma_ref[15] +
                  cos_delta * cos_delta * cos_delta * gamma_ref[17] - sin_delta * sin_delta * sin_delta * gamma_ref[24] -
                  2.0 * cos_delta * sin_delta * sin_delta * gamma_ref[25] - cos_delta * cos_delta * sin_delta * gamma_ref[27];
  gamma_rot[18] = cos_delta * sin_delta * gamma_ref[16] + cos_delta * cos_delta * gamma_ref[18] - sin_delta * sin_delta * gamma_ref[26] -
                  cos_delta * sin_delta * gamma_ref[28];
  gamma_rot[19] = cos_delta * gamma_ref[19] - sin_delta * gamma_ref[29];
  gamma_rot[20] = sin_delta * gamma_ref[10] + cos_delta * gamma_ref[20];
  gamma_rot[21] = cos_delta * sin_delta * gamma_ref[11] - sin_delta * sin_delta * gamma_ref[12] + cos_delta * cos_delta * gamma_ref[21] -
                  cos_delta * sin_delta * gamma_ref[22];
  gamma_rot[22] = sin_delta * sin_delta * gamma_ref[11] + cos_delta * sin_delta * gamma_ref[12] + cos_delta * sin_delta * gamma_ref[21] +
                  cos_delta * cos_delta * gamma_ref[22];
  gamma_rot[23] = sin_delta * gamma_ref[13] + cos_delta * gamma_ref[23];
  gamma_rot[24] = cos_delta * cos_delta * sin_delta * gamma_ref[14] - 2.0 * cos_delta * sin_delta * sin_delta * gamma_ref[15] +
                  sin_delta * sin_delta * sin_delta * gamma_ref[17] + cos_delta * cos_delta * cos_delta * gamma_ref[24] -
                  2.0 * cos_delta * cos_delta * sin_delta * gamma_ref[25] + cos_delta * sin_delta * sin_delta * gamma_ref[27];
  gamma_rot[25] = cos_delta * sin_delta * sin_delta * gamma_ref[14] - sin_delta * sin_delta * sin_delta * gamma_ref[15] +
                  cos_delta * cos_delta * sin_delta * gamma_ref[15] - cos_delta * sin_delta * sin_delta * gamma_ref[17] +
                  cos_delta * cos_delta * sin_delta * gamma_ref[24] - cos_delta * sin_delta * sin_delta * gamma_ref[25] +
                  cos_delta * cos_delta * cos_delta * gamma_ref[25] - cos_delta * cos_delta * sin_delta * gamma_ref[27];
  gamma_rot[26] = cos_delta * sin_delta * gamma_ref[16] - sin_delta * sin_delta * gamma_ref[18] + cos_delta * cos_delta * gamma_ref[26] -
                  cos_delta * sin_delta * gamma_ref[28];
  gamma_rot[27] = sin_delta * sin_delta * sin_delta * gamma_ref[14] + 2.0 * cos_delta * sin_delta * sin_delta * gamma_ref[15] +
                  cos_delta * cos_delta * sin_delta * gamma_ref[17] + cos_delta * sin_delta * sin_delta * gamma_ref[24] +
                  2.0 * cos_delta * cos_delta * sin_delta * gamma_ref[25] + cos_delta * cos_delta * cos_delta * gamma_ref[27];
  gamma_rot[28] = sin_delta * sin_delta * gamma_ref[16] + cos_delta * sin_delta * gamma_ref[18] + cos_delta * sin_delta * gamma_ref[26] +
                  cos_delta * cos_delta * gamma_ref[28];
  gamma_rot[29] = sin_delta * gamma_ref[19] + cos_delta * gamma_ref[29];
  gamma_rot[30] = gamma_ref[30];
  gamma_rot[31] = cos_delta * gamma_ref[31] - sin_delta * gamma_ref[32];
  gamma_rot[32] = sin_delta * gamma_ref[31] + cos_delta * gamma_ref[32];
  gamma_rot[33] = gamma_ref[33];
  gamma_rot[34] = cos_delta * cos_delta * gamma_ref[34] - 2.0 * cos_delta * sin_delta * gamma_ref[35] + sin_delta * sin_delta * gamma_ref[37];
  gamma_rot[35] = cos_delta * sin_delta * gamma_ref[34] - sin_delta * sin_delta * gamma_ref[35] + cos_delta * cos_delta * gamma_ref[35] -
                  cos_delta * sin_delta * gamma_ref[37];
  gamma_rot[36] = cos_delta * gamma_ref[36] - sin_delta * gamma_ref[38];
  gamma_rot[37] = sin_delta * sin_delta * gamma_ref[34] + 2.0 * cos_delta * sin_delta * gamma_ref[35] + cos_delta * cos_delta * gamma_ref[37];
  gamma_rot[38] = sin_delta * gamma_ref[36] + cos_delta * gamma_ref[38];
  gamma_rot[39] = gamma_ref[39];
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_rotate_about_z

/**
 * Interpolate Cartesian geodesic tensors at one spatial position.
 *
 * The caller supplies a trusted spatial context and already-mapped time-slice
 * payload pointers. The helper builds one native `(r, theta)` stencil, remaps the
 * payload-read indices once with explicit spherical reflections, reads tensor
 * components directly from mapped payload memory for each requested slice,
 * rotates each remapped node back to one common stored phi plane, interpolates
 * there, rotates to the target azimuth, and writes flat per-slice outputs.
 *
 * This helper assumes the stored payload has constant native grid spacing in
 * `r` and `theta`, and that axisymmetry is represented by exactly two stored phi
 * planes. It therefore performs no interpolation in `phi`: it interpolates only
 * on the uniform `(r, theta)` stencil. Reflected nodes that land on the opposite
 * stored phi plane are first rotated back to the selected stored reference plane
 * before weighted accumulation, and the final interpolated Cartesian-basis
 * tensors are then rotated from that reference plane to the target azimuth.
 *
 * @param[in] context Trusted spatial context.
 * @param[in] commondata Common runtime parameters.
 * @param[in] params Generated BHaH grid parameters.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param num_target_slices Number of mapped slice payload pointers.
 * @param[in] slice_payloads Mapped slice payload pointers.
 * @param[out] g4dd_out Flat metric output.
 * @param[out] gamma4udd_out Flat Christoffel output.
 * @return Status code indicating success or the interpolation failure reason.
 *
 * @note Each `slice_payloads` entry must point to the beginning of one mapped
 * 3D-grid payload and remain valid for the duration of this call. The spatial
 * stencil half-width is read from
 * `commondata->numerical_spacetime_spatial_interp_order`; the actual number of
 * radial and polar Lagrange nodes is `2*n+1`.
 */
int azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(const azimuthal_symmetry_spatial_lagrange_context_struct *restrict context,
                                                                      const commondata_struct *restrict commondata,
                                                                      const params_struct *restrict params, const REAL x, const REAL y, const REAL z,
                                                                      const int num_target_slices, const double *const *restrict slice_payloads,
                                                                      REAL *restrict g4dd_out, REAL *restrict gamma4udd_out) {
  // Step 1: Bind the trusted parsed context.
  const REAL xCart[3] = {x, y, z};
  REAL xx_target[3];
  int center_idx[3];

  // Step 2: Convert the target Cartesian point to native spherical coordinates.
  Cart_to_xx_and_nearest_i0i1i2__rfm__Spherical(params, xCart, xx_target, center_idx);

  const REAL target_r = xx_target[0];
  const REAL target_theta = xx_target[1];
  REAL target_phi = xx_target[2];
  const REAL rho_sq = x * x + y * y;
  const REAL origin_epsilon = 1.0e-14;
  const REAL axis_rho_epsilon = 1.0e-14;
  const int n_interp_ghosts = commondata->numerical_spacetime_spatial_interp_order;
  const long int interp_order_long = 2L * (long int)n_interp_ghosts + 1L;
  const REAL phi0 = (REAL)context->stored_phi_samples[0];
  const REAL phi1 = (REAL)context->stored_phi_samples[1];
  const REAL r_max_supported = (REAL)(params->xxmin0 + (((params->Nxx0 - 1) + 0.5) * params->dxx0));
  const REAL r_support_tol = (REAL)(1.0e-12 * fmax(1.0, fabs((double)r_max_supported)));
  int phi_plane;
  int i2_base;
  REAL phi_ref = 0.0;
  REAL phi_delta = phi1 - phi0;

  if (!isfinite((double)target_r) || !isfinite((double)target_theta) || !isfinite((double)target_phi)) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  }
  if (target_phi >= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - (REAL)1.0e-12 &&
      target_phi <= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + (REAL)1.0e-12) {
    target_phi -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  }
  if (target_phi < (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - 1.0e-12) ||
      target_phi > (REAL)(AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + 1.0e-12)) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  }

  // Reject the origin/axis because the azimuth is degenerate and this helper
  // recovers phi through rotation, not through a separate interpolation.
  if (target_r <= origin_epsilon || rho_sq <= axis_rho_epsilon * axis_rho_epsilon) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  }
  while (phi_delta <= (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI))
    phi_delta += (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  while (phi_delta > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)
    phi_delta -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (params->Nxx2 != 2 || !isfinite((double)phi0) || !isfinite((double)phi1) ||
      fabs((double)(fabs((double)phi_delta) - AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)) > 1.0e-12) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  }
  if (n_interp_ghosts < 0 || n_interp_ghosts > AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      interp_order_long > (long int)params->Nxx0 || interp_order_long > (long int)params->Nxx1) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  }

  // Step 3: Select one of the two stored phi planes and build 2D interpolation
  // coefficients. No interpolation in phi is performed anywhere in this helper;
  // phi is handled later by rotating tensors from the stored plane to the
  // target azimuth.
  const int interp_order = (int)interp_order_long;
  phi_plane = target_phi < (REAL)0.0 ? 0 : 1;
  i2_base = NGHOSTS + phi_plane;
  phi_ref = (REAL)context->stored_phi_samples[phi_plane];

  REAL inv_denom[interp_order];
  REAL src_r_stencil[interp_order];
  REAL src_theta_stencil[interp_order];
  REAL diffs_r[interp_order];
  REAL diffs_theta[interp_order];
  REAL coeff_r[interp_order];
  REAL coeff_theta[interp_order];
  uint64_t mapped_point_index[interp_order][interp_order];
  int mapped_phi_plane[interp_order][interp_order];
  const REAL normalization_2d = pow((REAL)(params->dxx0 * params->dxx1), -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int i0_raw = center_idx[0] + (u - n_interp_ghosts);
    // Keep the polynomial nodes on the raw extended grid, even if the later
    // payload reads are recovered back into the interior.
    src_r_stencil[u] = (REAL)(params->xxmin0 + (((i0_raw - NGHOSTS) + 0.5) * params->dxx0));
  } // END LOOP: for u over radial stencil nodes
  for (int v = 0; v < interp_order; v++) {
    const int i1_raw = center_idx[1] + (v - n_interp_ghosts);
    src_theta_stencil[v] = (REAL)(params->xxmin1 + (((i1_raw - NGHOSTS) + 0.5) * params->dxx1));
  } // END LOOP: for v over theta stencil nodes

  if (src_r_stencil[interp_order - 1] > r_max_supported + r_support_tol) {
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  } // END IF: largest raw radial stencil node exceeded the supported payload radius

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(interp_order, target_r, src_r_stencil, diffs_r);
  compute_diffs_xi(interp_order, target_theta, src_theta_stencil, diffs_theta);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_r, coeff_r);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_theta, coeff_theta);

  // Step 4: Map each raw stencil node once for this photon target.
  for (int v = 0; v < interp_order; v++) {
    const REAL theta_ext = src_theta_stencil[v];

    for (int u = 0; u < interp_order; u++) {
      const REAL r_ext = src_r_stencil[u];
      int i0i1i2_map[3] = {-1, -1, -1};

      const int map_status = azimuthal_symmetry_spatial_lagrange_map_extended_node(params, r_ext, theta_ext, i2_base, i0i1i2_map);
      if (map_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS) {
        return map_status;
      }
      const uint64_t j0 = (uint64_t)(i0i1i2_map[0] - NGHOSTS);
      const uint64_t j1 = (uint64_t)(i0i1i2_map[1] - NGHOSTS);
      const uint64_t j2 = (uint64_t)(i0i1i2_map[2] - NGHOSTS);
      const int mapped_phi = i0i1i2_map[2] - NGHOSTS;
      mapped_point_index[v][u] = j0 + (uint64_t)params->Nxx0 * (j1 + (uint64_t)params->Nxx1 * j2);
      if (mapped_phi < 0 || mapped_phi >= 2)
        return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
      mapped_phi_plane[v][u] = mapped_phi;
    } // END LOOP: for u over radial stencil nodes during remap precompute
  } // END LOOP: for v over theta stencil nodes during remap precompute

  // Step 5: Interpolate each requested time slice independently.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const double *restrict slice_payload = slice_payloads[which_slice];
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT] = {0};

    // Step 5.a: Read, basis-align, and accumulate the weighted stencil nodes for this slice.
    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload + mapped_point_index[v][u] * (uint64_t)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT + 3ULL;
        const REAL weight = normalization_2d * coeff_theta[v] * coeff_r[u];
        const int apply_pi_shift = mapped_phi_plane[v][u] != phi_plane;

        // Same-plane nodes accumulate directly. Opposite-plane nodes differ
        // only by the diagonal pi-rotation sign pattern in Cartesian basis.
        azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_to_common_plane(weight, apply_pi_shift, tensor_record, tensor_ref);
      } // END LOOP: for u over radial stencil nodes
    } // END LOOP: for v over theta stencil nodes

    // Step 5.b: Rotate the interpolated Cartesian-basis tensors to the target azimuth.
    azimuthal_symmetry_spatial_lagrange_rotate_about_z(target_phi - phi_ref, &tensor_ref[0],
                                                       &tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
                                                       &g4dd_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
                                                       &gamma4udd_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT]);
  } // END LOOP: for which_slice over requested slice indices

  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical

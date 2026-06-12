#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation_lagrange_uniform.h"
#include <stdint.h>

/**
 * Map one full-payload storage-grid index triplet to the serialized point order.
 *
 * The combined numerical payload stores one point record for every storage-grid
 * location, including ghost zones, in i2-major, i0-fast order. For the current
 * writer contract the raw stencil indices are already storage indices, so this
 * helper only validates bounds and converts them to the serialized payload
 * point-record index.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param i0_storage Storage-grid x0 index, including ghost zones.
 * @param i1_storage Storage-grid x1 index, including ghost zones.
 * @param i2_storage Storage-grid x2 index, including ghost zones.
 * @param[out] point_index Serialized zero-based point-record index.
 * @return `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS` if the storage
 * indices are inside the full payload, otherwise
 * `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL`.
 */
static int azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(const params_struct *restrict params, const int i0_storage,
                                                                                     const int i1_storage, const int i2_storage,
                                                                                     uint64_t *restrict point_index) {
  const int payload_i0_count = params->Nxx_plus_2NGHOSTS0;
  const int payload_i1_count = params->Nxx_plus_2NGHOSTS1;
  const int payload_i2_count = params->Nxx_plus_2NGHOSTS2;

  if (i0_storage < 0 || i0_storage >= payload_i0_count || i1_storage < 0 || i1_storage >= payload_i1_count || i2_storage < 0 ||
      i2_storage >= payload_i2_count)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

  *point_index = (uint64_t)i0_storage + (uint64_t)payload_i0_count * ((uint64_t)i1_storage + (uint64_t)payload_i1_count * (uint64_t)i2_storage);
  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices

/**
 * Accumulate one serialized tensor record onto the selected stored phi plane.
 *
 * The ghost-zone-aware payload already stores Cartesian-basis tensor records at
 * the storage-grid location that the raw stencil node requested. After the
 * direct payload read no additional per-node axisymmetry transform is needed
 * before weighted accumulation on the selected stored phi plane.
 *
 * @param weight Weighted 2D Lagrange coefficient for this stencil node.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[in,out] tensor_ref Spatial interpolation accumulator on the
 * selected common stored phi plane.
 */
static void
azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct(const REAL weight, const double *restrict tensor_record,
                                                                    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT]) {
  const double *restrict gamma_record = tensor_record + AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT;

  for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++)
    tensor_ref[comp] += weight * (REAL)tensor_record[comp];
  for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT; comp++)
    tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] += weight * (REAL)gamma_record[comp];
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct

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
 * payload pointers. The helper builds one native `(r, theta)` stencil, converts
 * each raw storage-grid stencil node to the corresponding full-payload point
 * record index once, reads tensor components directly from mapped payload memory
 * for each requested slice, interpolates on the selected stored phi plane,
 * rotates to the target azimuth, and writes flat per-slice outputs.
 *
 * This helper assumes the stored payload has constant native grid spacing in
 * `r` and `theta`, and that axisymmetry is represented by exactly two stored phi
 * planes. It therefore performs no interpolation in `phi`: it interpolates only
 * on the uniform `(r, theta)` stencil. The payload is expected to include
 * ghost-zone point records, and those records are treated as authoritative
 * Cartesian-basis tensor values at their storage-grid locations. The final
 * interpolated Cartesian-basis tensors are then rotated from the selected stored
 * reference plane to the target azimuth.
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
 * @return `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS` on success,
 * `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET` for invalid target
 * coordinates, or
 * `AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL` for
 * unsupported stencil or payload geometry.
 *
 * @note Each `slice_payloads` entry must point to the beginning of one mapped
 * ghost-zone-inclusive 3D-grid payload and remain valid for the duration of this
 * call. The spatial stencil half-width is read from
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
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS0 || interp_order_long > (long int)params->Nxx_plus_2NGHOSTS1) {
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
  int i0_storage_stencil[interp_order];
  int i1_storage_stencil[interp_order];
  const REAL normalization_2d = pow((REAL)(params->dxx0 * params->dxx1), -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int i0_raw = center_idx[0] + (u - n_interp_ghosts);
    i0_storage_stencil[u] = i0_raw;
    // Keep the polynomial nodes on raw ghost-zone-inclusive storage-grid
    // indices.
    src_r_stencil[u] = (REAL)(params->xxmin0 + (((i0_raw - NGHOSTS) + 0.5) * params->dxx0));
  } // END LOOP: for u over radial stencil nodes
  for (int v = 0; v < interp_order; v++) {
    const int i1_raw = center_idx[1] + (v - n_interp_ghosts);
    i1_storage_stencil[v] = i1_raw;
    src_theta_stencil[v] = (REAL)(params->xxmin1 + (((i1_raw - NGHOSTS) + 0.5) * params->dxx1));
  } // END LOOP: for v over theta stencil nodes

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(interp_order, target_r, src_r_stencil, diffs_r);
  compute_diffs_xi(interp_order, target_theta, src_theta_stencil, diffs_theta);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_r, coeff_r);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_theta, coeff_theta);

  // Step 4: Convert each raw storage-grid stencil node to one payload point index.
  for (int v = 0; v < interp_order; v++) {
    for (int u = 0; u < interp_order; u++) {
      const int index_status = azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(
          params, i0_storage_stencil[u], i1_storage_stencil[v], i2_base, &mapped_point_index[v][u]);
      if (index_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
        return index_status;
    } // END LOOP: for u over radial stencil nodes during payload-index precompute
  } // END LOOP: for v over theta stencil nodes during payload-index precompute

  // Step 5: Interpolate each requested time slice independently.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const double *restrict slice_payload = slice_payloads[which_slice];
    REAL tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_TENSOR_COMPONENT_COUNT] = {0};

    // Step 5.a: Read ghost-zone-aware payload records and directly accumulate
    // the weighted stencil nodes for this slice.
    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload + mapped_point_index[v][u] * (uint64_t)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT + 3ULL;
        const REAL weight = normalization_2d * coeff_theta[v] * coeff_r[u];

        azimuthal_symmetry_spatial_lagrange_accumulate_tensor_record_direct(weight, tensor_record, tensor_ref);
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

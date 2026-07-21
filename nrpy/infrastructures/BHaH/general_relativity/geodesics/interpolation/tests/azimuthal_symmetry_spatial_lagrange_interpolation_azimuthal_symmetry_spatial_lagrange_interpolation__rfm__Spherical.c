#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "interpolation/interpolation_lagrange_uniform.h"
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
  REAL g_ref[4][4] = {{0}};
  REAL g_dst[4][4] = {{0}};
  REAL gamma_ref_full[4][4][4] = {{{0}}};
  REAL gamma_dst_full[4][4][4] = {{{0}}};
  REAL rotation_matrix[4][4] = {{0}};
  REAL rotation_matrix_transpose[4][4] = {{0}};
  const REAL cos_delta = cos(delta_phi);
  const REAL sin_delta = sin(delta_phi);
  int gamma_index = 0;

  // Step 1: Unpack the stored upper-triangular metric components.
  g_ref[0][0] = g4dd_ref[0];
  g_ref[0][1] = g_ref[1][0] = g4dd_ref[1];
  g_ref[0][2] = g_ref[2][0] = g4dd_ref[2];
  g_ref[0][3] = g_ref[3][0] = g4dd_ref[3];
  g_ref[1][1] = g4dd_ref[4];
  g_ref[1][2] = g_ref[2][1] = g4dd_ref[5];
  g_ref[1][3] = g_ref[3][1] = g4dd_ref[6];
  g_ref[2][2] = g4dd_ref[7];
  g_ref[2][3] = g_ref[3][2] = g4dd_ref[8];
  g_ref[3][3] = g4dd_ref[9];

  // Step 2: Unpack the serialized Christoffels, mirroring lower-index symmetry.
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
        gamma_ref_full[alpha][mu][nu] = gamma_ref[gamma_index];
        gamma_ref_full[alpha][nu][mu] = gamma_ref[gamma_index];
        gamma_index++;
      } // END LOOP: for nu over upper-triangular lower-index Christoffel entries
    } // END LOOP: for mu over lower Christoffel index rows
  } // END LOOP: for alpha over upper Christoffel index

  // Step 3: Build the active z-axis rotation matrix in Cartesian components.
  rotation_matrix[0][0] = 1.0;
  rotation_matrix[1][1] = cos_delta;
  rotation_matrix[1][2] = -sin_delta;
  rotation_matrix[2][1] = sin_delta;
  rotation_matrix[2][2] = cos_delta;
  rotation_matrix[3][3] = 1.0;

  // Step 4: For pure rotations, the inverse acting on lower indices is Lambda^T.
  for (int a = 0; a < 4; a++) {
    for (int b = 0; b < 4; b++) {
      rotation_matrix_transpose[a][b] = rotation_matrix[b][a];
    } // END LOOP: for b over rotation-matrix columns
  } // END LOOP: for a over rotation-matrix rows

  // Step 5: Rotate the covariant metric with the lower-index transform.
  for (int mu = 0; mu < 4; mu++) {
    for (int nu = 0; nu < 4; nu++) {
      REAL sum = 0.0;
      for (int a = 0; a < 4; a++) {
        for (int b = 0; b < 4; b++) {
          sum += rotation_matrix_transpose[a][mu] * rotation_matrix_transpose[b][nu] * g_ref[a][b];
        } // END LOOP: for b over source metric columns
      } // END LOOP: for a over source metric rows
      g_dst[mu][nu] = sum;
    } // END LOOP: for nu over destination metric columns
  } // END LOOP: for mu over destination metric rows

  // Step 6: Rotate the connection with one upper and two lower tensor factors.
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = 0; nu < 4; nu++) {
        REAL sum = 0.0;
        for (int b = 0; b < 4; b++) {
          for (int g = 0; g < 4; g++) {
            for (int d = 0; d < 4; d++) {
              sum += rotation_matrix[alpha][b] * rotation_matrix_transpose[g][mu] * rotation_matrix_transpose[d][nu] * gamma_ref_full[b][g][d];
            } // END LOOP: for d over source lower Christoffel column index
          } // END LOOP: for g over source lower Christoffel row index
        } // END LOOP: for b over source upper Christoffel index
        gamma_dst_full[alpha][mu][nu] = sum;
      } // END LOOP: for nu over destination lower Christoffel column index
    } // END LOOP: for mu over destination lower Christoffel row index
  } // END LOOP: for alpha over destination upper Christoffel index

  // Step 7: Repack into the writer's serialized ordering.
  g4dd_rot[0] = g_dst[0][0];
  g4dd_rot[1] = g_dst[0][1];
  g4dd_rot[2] = g_dst[0][2];
  g4dd_rot[3] = g_dst[0][3];
  g4dd_rot[4] = g_dst[1][1];
  g4dd_rot[5] = g_dst[1][2];
  g4dd_rot[6] = g_dst[1][3];
  g4dd_rot[7] = g_dst[2][2];
  g4dd_rot[8] = g_dst[2][3];
  g4dd_rot[9] = g_dst[3][3];

  gamma_index = 0;
  for (int alpha = 0; alpha < 4; alpha++) {
    for (int mu = 0; mu < 4; mu++) {
      for (int nu = mu; nu < 4; nu++) {
        gamma_rot[gamma_index] = gamma_dst_full[alpha][mu][nu];
        gamma_index++;
      } // END LOOP: for nu over serialized upper-triangular Christoffel entries
    } // END LOOP: for mu over serialized Christoffel row index
  } // END LOOP: for alpha over serialized Christoffel upper index
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
  Cart_to_xx_and_nearest_i0i1i2_assume_valid(params, xCart, xx_target, center_idx);

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
    REAL g4dd_node_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL gamma4udd_node_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT];
    REAL g4dd_node_common[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL gamma4udd_node_common[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT];

    // Step 5.a: Read, basis-align, and accumulate the weighted stencil nodes for this slice.
    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload + mapped_point_index[v][u] * (uint64_t)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT + 3ULL;
        const REAL weight = normalization_2d * coeff_theta[v] * coeff_r[u];
        const REAL node_phi_ref = (REAL)context->stored_phi_samples[mapped_phi_plane[v][u]];

        for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++) {
          g4dd_node_ref[comp] = (REAL)tensor_record[comp];
        } // END LOOP: for comp over serialized metric payload components
        for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT; comp++) {
          gamma4udd_node_ref[comp] = (REAL)tensor_record[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp];
        } // END LOOP: for comp over serialized Christoffel payload components

        // Rotate each remapped node back to the selected stored reference
        // plane so weighted accumulation never mixes Cartesian bases from
        // different azimuthal orientations.
        azimuthal_symmetry_spatial_lagrange_rotate_about_z(phi_ref - node_phi_ref, g4dd_node_ref, gamma4udd_node_ref, g4dd_node_common,
                                                           gamma4udd_node_common);

        for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++) {
          tensor_ref[comp] += weight * g4dd_node_common[comp];
        } // END LOOP: for comp over metric components on the common reference plane
        for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_GAMMA_COMPONENT_COUNT; comp++) {
          tensor_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT + comp] += weight * gamma4udd_node_common[comp];
        } // END LOOP: for comp over Christoffel components on the common reference plane
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

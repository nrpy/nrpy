#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include "differentiate_interpolation_lagrange_uniform.h"
#include "interpolation_lagrange_uniform.h"
#include <math.h>
#include <stdint.h>

/**
 * Map one full-payload storage-grid index triplet to the serialized point order.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param i0_storage Storage-grid x0 index, including ghost zones.
 * @param i1_storage Storage-grid x1 index, including ghost zones.
 * @param i2_storage Storage-grid x2 index, including ghost zones.
 * @param[out] point_index Serialized zero-based point-record index.
 * @return Interpolation status code.
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
 * Evaluate the inverse reference-metric coordinate Jacobian at one native point.
 *
 * The generated entries satisfy
 * `inverse_jacobian[native_direction][cartesian_direction] =
 * d xx^native_direction / d xCart^cartesian_direction`.
 *
 * @param[in] params Generated BHaH parameter struct.
 * @param[in] xx Native target coordinates.
 * @param[out] inverse_jacobian Native-with-respect-to-Cartesian Jacobian.
 */
static void azimuthal_symmetry_spatial_lagrange_inverse_jacobian(const params_struct *restrict params, const REAL xx[3],
                                                                 REAL inverse_jacobian[3][3]) {
  const REAL xx0 = xx[0];
  const REAL xx1 = xx[1];
  const REAL xx2 = xx[2];
  const REAL tmp0 = cos(xx1);
  const REAL tmp1 = (1.0 / (params->SINHWRHO));
  const REAL tmp7 = params->AMPLRHO - params->rho_slope;
  const REAL tmp10 = (1.0 / (params->SINHWZ));
  const REAL tmp16 = sin(xx1);
  const REAL tmp6 = (1.0 / (exp(tmp1) - exp(-tmp1)));
  const REAL tmp14 = (params->AMPLZ - params->z_slope) / (exp(tmp10) - exp(-tmp10));
  const REAL tmp3 = exp(tmp1 * xx0);
  const REAL tmp4 = exp(-tmp1 * xx0);
  const REAL tmp8 = tmp6 * tmp7 * ((xx0) * (xx0));
  const REAL tmp12 = exp(tmp10 * xx2);
  const REAL tmp13 = exp(-tmp10 * xx2);
  const REAL tmp5 = tmp3 - tmp4;
  const REAL tmp15 = params->z_slope + tmp14 * ((xx2) * (xx2)) * (tmp10 * tmp12 + tmp10 * tmp13) + 2 * tmp14 * xx2 * (tmp12 - tmp13);
  const REAL tmp9 = params->rho_slope * xx0 + tmp5 * tmp8;
  const REAL tmp17 = params->rho_slope + 2 * tmp5 * tmp6 * tmp7 * xx0 + tmp8 * (tmp1 * tmp3 + tmp1 * tmp4);
  const REAL tmp19 = ((tmp16) * (tmp16)) * tmp17 * tmp9;
  const REAL tmp20 = ((tmp0) * (tmp0)) * tmp17 * tmp9;
  const REAL tmp21 = (1.0 / (tmp15 * tmp19 + tmp15 * tmp20));
  const REAL tmp23 = tmp15 * tmp21 * tmp9;
  const REAL tmp24 = tmp15 * tmp17 * tmp21;
  inverse_jacobian[0][0] = tmp0 * tmp23;
  inverse_jacobian[0][1] = tmp16 * tmp23;
  inverse_jacobian[0][2] = 0;
  inverse_jacobian[1][0] = -tmp16 * tmp24;
  inverse_jacobian[1][1] = tmp0 * tmp24;
  inverse_jacobian[1][2] = 0;
  inverse_jacobian[2][0] = 0;
  inverse_jacobian[2][1] = 0;
  inverse_jacobian[2][2] = tmp21 * (tmp19 + tmp20);
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_inverse_jacobian

/**
 * Accumulate metric values and two native derivatives from one payload record.
 *
 * The record begins with the ten Cartesian-basis metric components after its
 * three Cartesian coordinate values have been skipped by the caller.
 *
 * @param value_weight Two-dimensional Lagrange value weight.
 * @param interp_0_derivative_weight Derivative weight in the first native interpolation direction.
 * @param interp_1_derivative_weight Derivative weight in the second native interpolation direction.
 * @param[in] tensor_record Serialized payload record beginning at the metric fields.
 * @param[in,out] g4dd_ref Interpolated metric on the selected stored phi plane.
 * @param[in,out] g4dd_interp_0_ref First native interpolation-direction derivative.
 * @param[in,out] g4dd_interp_1_ref Second native interpolation-direction derivative.
 */
static void azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct(
    const REAL value_weight, const REAL interp_0_derivative_weight, const REAL interp_1_derivative_weight, const double *restrict tensor_record,
    REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],

    REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT]) {
  for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++) {
    const REAL metric_value = (REAL)tensor_record[comp];
    g4dd_ref[comp] += value_weight * metric_value;

    g4dd_interp_0_ref[comp] += interp_0_derivative_weight * metric_value;
    g4dd_interp_1_ref[comp] += interp_1_derivative_weight * metric_value;
  } // END LOOP: for comp over serialized metric components
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct

/**
 * Rotate the metric and native metric derivatives about the z axis.
 *
 * The two differentiated Lagrange reconstructions rotate as ordinary metric
 * tensors. The native phi derivative is generated analytically by
 * differentiating the same z-axis rotation used for the metric values. No
 * interpolation in phi is performed.
 *
 * @param delta_phi Active rotation angle from stored plane to target azimuth.
 * @param[in] g4dd_ref Metric components on the selected stored phi plane.
 * @param[in] g4dd_interp_0_ref First interpolated native derivative on that plane.
 * @param[in] g4dd_interp_1_ref Second interpolated native derivative on that plane.
 * @param[out] g4dd_rot Metric components at the target azimuth.
 * @param[out] g4dd_native_derivatives_rot Native derivatives at the target azimuth.
 */
static void azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z(
    const REAL delta_phi, const REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    const REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    const REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT],
    REAL g4dd_native_derivatives_rot[3][AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT]) {
  const REAL cos_delta = cos(delta_phi);
  const REAL sin_delta = sin(delta_phi);

  for (int native_direction = 0; native_direction < 3; native_direction++)
    for (int comp = 0; comp < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; comp++)
      g4dd_native_derivatives_rot[native_direction][comp] = 0.0;

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
  g4dd_native_derivatives_rot[0][0] = g4dd_interp_0_ref[0];
  g4dd_native_derivatives_rot[0][1] = cos_delta * g4dd_interp_0_ref[1] - sin_delta * g4dd_interp_0_ref[2];
  g4dd_native_derivatives_rot[0][2] = sin_delta * g4dd_interp_0_ref[1] + cos_delta * g4dd_interp_0_ref[2];
  g4dd_native_derivatives_rot[0][3] = g4dd_interp_0_ref[3];
  g4dd_native_derivatives_rot[0][4] = cos_delta * cos_delta * g4dd_interp_0_ref[4] - 2.0 * cos_delta * sin_delta * g4dd_interp_0_ref[5] +
                                      sin_delta * sin_delta * g4dd_interp_0_ref[7];
  g4dd_native_derivatives_rot[0][5] = cos_delta * sin_delta * g4dd_interp_0_ref[4] - sin_delta * sin_delta * g4dd_interp_0_ref[5] +
                                      cos_delta * cos_delta * g4dd_interp_0_ref[5] - cos_delta * sin_delta * g4dd_interp_0_ref[7];
  g4dd_native_derivatives_rot[0][6] = cos_delta * g4dd_interp_0_ref[6] - sin_delta * g4dd_interp_0_ref[8];
  g4dd_native_derivatives_rot[0][7] = sin_delta * sin_delta * g4dd_interp_0_ref[4] + 2.0 * cos_delta * sin_delta * g4dd_interp_0_ref[5] +
                                      cos_delta * cos_delta * g4dd_interp_0_ref[7];
  g4dd_native_derivatives_rot[0][8] = sin_delta * g4dd_interp_0_ref[6] + cos_delta * g4dd_interp_0_ref[8];
  g4dd_native_derivatives_rot[0][9] = g4dd_interp_0_ref[9];
  g4dd_native_derivatives_rot[1][0] = 0.0;
  g4dd_native_derivatives_rot[1][1] = -sin_delta * g4dd_ref[1] - cos_delta * g4dd_ref[2];
  g4dd_native_derivatives_rot[1][2] = cos_delta * g4dd_ref[1] - sin_delta * g4dd_ref[2];
  g4dd_native_derivatives_rot[1][3] = 0.0;
  g4dd_native_derivatives_rot[1][4] = -2.0 * cos_delta * sin_delta * g4dd_ref[4] + 2.0 * sin_delta * sin_delta * g4dd_ref[5] -
                                      2.0 * cos_delta * cos_delta * g4dd_ref[5] + 2.0 * cos_delta * sin_delta * g4dd_ref[7];
  g4dd_native_derivatives_rot[1][5] = -sin_delta * sin_delta * g4dd_ref[4] + cos_delta * cos_delta * g4dd_ref[4] -
                                      4.0 * cos_delta * sin_delta * g4dd_ref[5] + sin_delta * sin_delta * g4dd_ref[7] -
                                      cos_delta * cos_delta * g4dd_ref[7];
  g4dd_native_derivatives_rot[1][6] = -sin_delta * g4dd_ref[6] - cos_delta * g4dd_ref[8];
  g4dd_native_derivatives_rot[1][7] = 2.0 * cos_delta * sin_delta * g4dd_ref[4] - 2.0 * sin_delta * sin_delta * g4dd_ref[5] +
                                      2.0 * cos_delta * cos_delta * g4dd_ref[5] - 2.0 * cos_delta * sin_delta * g4dd_ref[7];
  g4dd_native_derivatives_rot[1][8] = cos_delta * g4dd_ref[6] - sin_delta * g4dd_ref[8];
  g4dd_native_derivatives_rot[1][9] = 0.0;
  g4dd_native_derivatives_rot[2][0] = g4dd_interp_1_ref[0];
  g4dd_native_derivatives_rot[2][1] = cos_delta * g4dd_interp_1_ref[1] - sin_delta * g4dd_interp_1_ref[2];
  g4dd_native_derivatives_rot[2][2] = sin_delta * g4dd_interp_1_ref[1] + cos_delta * g4dd_interp_1_ref[2];
  g4dd_native_derivatives_rot[2][3] = g4dd_interp_1_ref[3];
  g4dd_native_derivatives_rot[2][4] = cos_delta * cos_delta * g4dd_interp_1_ref[4] - 2.0 * cos_delta * sin_delta * g4dd_interp_1_ref[5] +
                                      sin_delta * sin_delta * g4dd_interp_1_ref[7];
  g4dd_native_derivatives_rot[2][5] = cos_delta * sin_delta * g4dd_interp_1_ref[4] - sin_delta * sin_delta * g4dd_interp_1_ref[5] +
                                      cos_delta * cos_delta * g4dd_interp_1_ref[5] - cos_delta * sin_delta * g4dd_interp_1_ref[7];
  g4dd_native_derivatives_rot[2][6] = cos_delta * g4dd_interp_1_ref[6] - sin_delta * g4dd_interp_1_ref[8];
  g4dd_native_derivatives_rot[2][7] = sin_delta * sin_delta * g4dd_interp_1_ref[4] + 2.0 * cos_delta * sin_delta * g4dd_interp_1_ref[5] +
                                      cos_delta * cos_delta * g4dd_interp_1_ref[7];
  g4dd_native_derivatives_rot[2][8] = sin_delta * g4dd_interp_1_ref[6] + cos_delta * g4dd_interp_1_ref[8];
  g4dd_native_derivatives_rot[2][9] = g4dd_interp_1_ref[9];
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z

/**
 * Interpolate the selected Cartesian geometry bundle.
 *
 * The helper reads the ten metric values and the selected secondary payload from
 * the full ghost-zone-inclusive records. ``g4DD`` reconstructs all spatial metric
 * derivatives and sets the ten temporal slots to zero. ``g4DD_d0`` reads and
 * rotates the ten stored temporal metric derivatives while reconstructing the
 * spatial derivatives. ``GammaUDD`` interpolates and rotates the forty stored
 * Christoffel components.
 *
 * Axisymmetry is handled without phi interpolation. The interpolated Cartesian
 * quantities are rotated from one stored reference-phi plane to the target
 * azimuth. The metric derivative index is transformed to Cartesian with the
 * inverse Jacobian generated from
 * `reference_metric[CoordSystem].Jac_dUrfm_dDCartUD` for the two metric methods.
 * The Christoffel rotation uses the constant Cartesian z-axis transformation.
 *
 * The neutral `rhs_geometry_out` buffer contains forty values per slice. Its
 * meaning is fixed by the selected code-generation method: metric derivatives
 * for the two metric methods or Christoffel components for ``GammaUDD``.
 *
 * @param[in] context Trusted spatial context.
 * @param[in] commondata Common interpolation parameters.
 * @param[in] params Generated BHaH grid and reference-metric parameters.
 * @param x Cartesian x coordinate.
 * @param y Cartesian y coordinate.
 * @param z Cartesian z coordinate.
 * @param num_target_slices Number of mapped slice payload pointers.
 * @param[in] slice_payloads Mapped ghost-zone-inclusive slice payload pointers.
 * @param[out] g4dd_out Flat metric output, ten values per slice.
 * @param[out] g4dd_dD_out Metric-derivative output, forty values per slice.
 * @return Interpolation status code.
 */
int azimuthal_symmetry_spatial_lagrange_interpolation__rfm__SinhCylindricalv2n2(
    const azimuthal_symmetry_spatial_lagrange_context_struct *restrict context, const commondata_struct *restrict commondata,
    const params_struct *restrict params, const REAL x, const REAL y, const REAL z, const int num_target_slices,
    const double *const *restrict slice_payloads, REAL *restrict g4dd_out, REAL *restrict rhs_geometry_out) {
  // Step 1: Validate pointers, then convert the target Cartesian point to
  // native coordinates.
  if (context == NULL || commondata == NULL || params == NULL || slice_payloads == NULL || g4dd_out == NULL || rhs_geometry_out == NULL ||
      num_target_slices <= 0)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  const REAL xCart[3] = {x, y, z};
  REAL xx_target[3];
  int center_idx[3];
  Cart_to_xx_and_nearest_i0i1i2__rfm__SinhCylindricalv2n2(params, xCart, xx_target, center_idx);

  const REAL target_interp_0 = xx_target[0];
  const REAL target_interp_1 = xx_target[2];
  REAL target_phi = xx_target[1];
  const REAL rho_sq = x * x + y * y;
  const REAL origin_epsilon = 1.0e-14;
  const REAL axis_rho_epsilon = 1.0e-14;
  const int n_interp_ghosts = commondata->numerical_spacetime_spatial_interp_order;
  const long int interp_order_long = 2L * (long int)n_interp_ghosts + 1L;
  const REAL phi0 = (REAL)context->stored_phi_samples[0];
  const REAL phi1 = (REAL)context->stored_phi_samples[1];
  REAL phi_delta = phi1 - phi0;

  if (!isfinite((double)target_interp_0) || !isfinite((double)target_interp_1) || !isfinite((double)target_phi))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  if (target_phi >= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - (REAL)1.0e-12 &&
      target_phi <= (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + (REAL)1.0e-12)
    target_phi -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (target_phi < (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI - 1.0e-12) ||
      target_phi > (REAL)(AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI + 1.0e-12))
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  // The azimuthal rotation and inverse cylindrical/spherical Jacobians are
  // singular on the symmetry axis, so retain the existing axis rejection.
  if (target_interp_0 <= origin_epsilon || rho_sq <= axis_rho_epsilon * axis_rho_epsilon)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;
  while (phi_delta <= (REAL)(-AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI))
    phi_delta += (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  while (phi_delta > (REAL)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)
    phi_delta -= (REAL)(2.0 * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI);
  if (params->Nxx1 != 2 || !isfinite((double)phi0) || !isfinite((double)phi1) ||
      fabs((double)(fabs((double)phi_delta) - AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_PI)) > 1.0e-12)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;
  if (n_interp_ghosts < 0 || n_interp_ghosts > AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      interp_order_long > (long int)params->Nxx_plus_2NGHOSTS0 || interp_order_long > (long int)params->Nxx_plus_2NGHOSTS2)
    return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

  // Step 2: Select one stored phi plane and construct value and derivative
  // coefficients on the two-dimensional native stencil.
  const int interp_order = (int)interp_order_long;
  const int phi_plane = target_phi < (REAL)0.0 ? 0 : 1;
  const int phi_plane_storage_index = NGHOSTS + phi_plane;
  const REAL phi_ref = (REAL)context->stored_phi_samples[phi_plane];

  REAL inv_denom[interp_order];
  REAL src_interp_0_stencil[interp_order];
  REAL src_interp_1_stencil[interp_order];
  REAL diffs_interp_0[interp_order];
  REAL diffs_interp_1[interp_order];
  REAL coeff_interp_0[interp_order];
  REAL coeff_interp_1[interp_order];
  REAL derivative_coeff_interp_0[interp_order];
  REAL derivative_coeff_interp_1[interp_order];
  uint64_t mapped_point_index[interp_order][interp_order];
  int interp_0_storage_stencil[interp_order];
  int interp_1_storage_stencil[interp_order];
  const REAL normalization_2d = pow((REAL)(params->dxx0 * params->dxx2), -(interp_order - 1));

  for (int u = 0; u < interp_order; u++) {
    const int interp_0_raw = center_idx[0] + (u - n_interp_ghosts);
    interp_0_storage_stencil[u] = interp_0_raw;
    src_interp_0_stencil[u] = (REAL)(params->xxmin0 + (((interp_0_raw - NGHOSTS) + 0.5) * params->dxx0));
  } // END LOOP: for u over first interpolation-dimension stencil nodes
  for (int v = 0; v < interp_order; v++) {
    const int interp_1_raw = center_idx[2] + (v - n_interp_ghosts);
    interp_1_storage_stencil[v] = interp_1_raw;
    src_interp_1_stencil[v] = (REAL)(params->xxmin2 + (((interp_1_raw - NGHOSTS) + 0.5) * params->dxx2));
  } // END LOOP: for v over second interpolation-dimension stencil nodes

  compute_inv_denom(interp_order, inv_denom);
  compute_diffs_xi(interp_order, target_interp_0, src_interp_0_stencil, diffs_interp_0);
  compute_diffs_xi(interp_order, target_interp_1, src_interp_1_stencil, diffs_interp_1);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_interp_0, coeff_interp_0);
  compute_lagrange_basis_coeffs_xi(interp_order, inv_denom, diffs_interp_1, coeff_interp_1);
  compute_lagrange_basis_derivative_coeffs_xi(interp_order, inv_denom, diffs_interp_0, derivative_coeff_interp_0);
  compute_lagrange_basis_derivative_coeffs_xi(interp_order, inv_denom, diffs_interp_1, derivative_coeff_interp_1);

  // Step 3: Convert each raw storage-grid stencil node to one payload index.
  for (int v = 0; v < interp_order; v++) {
    for (int u = 0; u < interp_order; u++) {
      const int i0_storage = interp_0_storage_stencil[u];
      const int i1_storage = phi_plane_storage_index;
      const int i2_storage = interp_1_storage_stencil[v];
      const int index_status = azimuthal_symmetry_spatial_lagrange_point_index_from_full_payload_indices(params, i0_storage, i1_storage, i2_storage,
                                                                                                         &mapped_point_index[v][u]);
      if (index_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
        return index_status;
    } // END LOOP: for u over first interpolation-dimension payload indices
  } // END LOOP: for v over second interpolation-dimension payload indices

  // Step 4: Evaluate the native-with-respect-to-Cartesian Jacobian once for
  // this target position. It is common to every requested time slice.
  REAL inverse_jacobian[3][3];
  azimuthal_symmetry_spatial_lagrange_inverse_jacobian(params, xx_target, inverse_jacobian);
  for (int native_direction = 0; native_direction < 3; native_direction++)
    for (int cartesian_direction = 0; cartesian_direction < 3; cartesian_direction++)
      if (!isfinite((double)inverse_jacobian[native_direction][cartesian_direction]))
        return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_INVALID_TARGET;

  // Step 5: Interpolate and differentiate each requested time slice.
  for (int which_slice = 0; which_slice < num_target_slices; which_slice++) {
    const double *restrict slice_payload = slice_payloads[which_slice];
    REAL g4dd_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_interp_0_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_interp_1_ref[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_rot[AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL g4dd_native_derivatives_rot[3][AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];

    if (slice_payload == NULL)
      return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_UNSUPPORTED_STENCIL;

    for (int v = 0; v < interp_order; v++) {
      for (int u = 0; u < interp_order; u++) {
        const double *restrict tensor_record =
            slice_payload + mapped_point_index[v][u] * (uint64_t)AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_RECORD_COMPONENT_COUNT + 3ULL;
        const REAL value_weight = normalization_2d * coeff_interp_1[v] * coeff_interp_0[u];
        const REAL interp_0_derivative_weight = normalization_2d * coeff_interp_1[v] * derivative_coeff_interp_0[u];
        const REAL interp_1_derivative_weight = normalization_2d * derivative_coeff_interp_1[v] * coeff_interp_0[u];

        azimuthal_symmetry_spatial_lagrange_accumulate_metric_record_direct(value_weight, interp_0_derivative_weight, interp_1_derivative_weight,
                                                                            tensor_record, g4dd_ref, g4dd_interp_0_ref, g4dd_interp_1_ref);
      } // END LOOP: for u over first interpolation-dimension values
    } // END LOOP: for v over second interpolation-dimension values

    azimuthal_symmetry_spatial_lagrange_rotate_metric_and_derivatives_about_z(target_phi - phi_ref, g4dd_ref, g4dd_interp_0_ref, g4dd_interp_1_ref,
                                                                              g4dd_rot, g4dd_native_derivatives_rot);

    REAL *restrict g4dd_slice_out = &g4dd_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT];
    REAL *restrict metric_derivative_slice_out =
        &rhs_geometry_out[which_slice * AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_METRIC_DERIVATIVE_COMPONENT_COUNT];

    // Step 5.a: Transform only the derivative index from native coordinates
    // to Cartesian coordinates, with direction order (dt, dx, dy, dz).
    for (int metric_component = 0; metric_component < AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_RT_G4_COMPONENT_COUNT; metric_component++) {
      g4dd_slice_out[metric_component] = g4dd_rot[metric_component];
      metric_derivative_slice_out[4 * metric_component] = 0.0;
      for (int cartesian_direction = 0; cartesian_direction < 3; cartesian_direction++) {
        REAL cartesian_derivative = 0.0;
        for (int native_direction = 0; native_direction < 3; native_direction++)
          cartesian_derivative +=
              inverse_jacobian[native_direction][cartesian_direction] * g4dd_native_derivatives_rot[native_direction][metric_component];
        metric_derivative_slice_out[4 * metric_component + cartesian_direction + 1] = cartesian_derivative;
      } // END LOOP: for cartesian_direction over x, y, z
    } // END LOOP: for metric_component over symmetric g4DD components
  } // END LOOP: for which_slice over requested slice indices

  return AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS;
} // END FUNCTION: azimuthal_symmetry_spatial_lagrange_interpolation__rfm__SinhCylindricalv2n2

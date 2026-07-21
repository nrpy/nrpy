#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "BHaH_defines.h"

/**
 *  Computes the 10 unique components of the KerrSchild_Cartesian metric $g_{mu nu}$ for a photon particle.
 *     @param commondata Struct containing global spacetime parameters.
 *     @param f_local Thread-local array containing the 1D flattened state vector.
 *     @param metric_local Thread-local array where the symmetric metric components are stored.
 */
BHAH_HD_INLINE void g4DD_metric_KerrSchild_Cartesian(const commondata_struct *restrict commondata, const double *restrict f_local,
                                                     double *restrict metric_local) {
#include "set_CodeParameters.h"
  //==========================================
  // METRIC EVALUATION & THREAD-LOCAL UNPACKING
  //==========================================
  // Extract spatial coordinates $x^i$ and compute $g_{mu nu}$.
  // Unpack position coordinates $x^i$ from the thread-local state vector.
  // Evaluated at compile time for state vector size: 9
  const double x = f_local[1];
  const double y = f_local[2];
  const double z = f_local[3];

  const REAL tmp0 = ((a_spin) * (a_spin));
  const REAL tmp4 = ((z) * (z));
  const REAL tmp6 = (1.0 / 2.0) * tmp4 + (1.0 / 2.0) * ((x) * (x)) + (1.0 / 2.0) * ((y) * (y)) +
                    (1.0 / 2.0) * sqrt(4 * tmp0 * tmp4 + ((-tmp0 + tmp4 + ((x) * (x)) + ((y) * (y))) * (-tmp0 + tmp4 + ((x) * (x)) + ((y) * (y)))));
  const REAL tmp7 = -1.0 / 2.0 * tmp0 + tmp6;
  const REAL tmp12 = (1.0 / 2.0) * tmp0 + tmp6;
  const REAL tmp8 = 2 * M_scale / (tmp0 * tmp4 + ((tmp7) * (tmp7)));
  const REAL tmp10 = sqrt(tmp7);
  const REAL tmp13 = (1.0 / (tmp12));
  const REAL tmp9 = pow(tmp7, 3.0 / 2.0) * tmp8;
  const REAL tmp11 = a_spin * y + tmp10 * x;
  const REAL tmp15 = -a_spin * x + tmp10 * y;
  const REAL tmp16 = tmp7 * tmp8 * z;
  const REAL tmp17 = tmp9 / ((tmp12) * (tmp12));
  metric_local[0] = tmp9 - 1;
  metric_local[1] = tmp11 * tmp13 * tmp9;
  metric_local[2] = tmp13 * tmp15 * tmp9;
  metric_local[3] = tmp16;
  metric_local[4] = ((tmp11) * (tmp11)) * tmp17 + 1;
  metric_local[5] = tmp11 * tmp15 * tmp17;
  metric_local[6] = tmp11 * tmp13 * tmp16;
  metric_local[7] = ((tmp15) * (tmp15)) * tmp17 + 1;
  metric_local[8] = tmp13 * tmp15 * tmp16;
  metric_local[9] = tmp10 * tmp4 * tmp8 + 1;

} // END FUNCTION: g4DD_metric_KerrSchild_Cartesian

#include "BHaH_defines.h"

/**
 *  Computes the 40 unique Christoffel symbols $Gamma^alpha_{mu nu}$ for the KerrSchild_Cartesian metric.
 *     @param commondata Struct containing global spacetime parameters.
 *     @param f_local Thread-local array containing the 1D flattened state vector.
 *     @param Gamma_local Thread-local array where the connection components are stored.
 */
BHAH_HD_INLINE void connections_KerrSchild_Cartesian(const commondata_struct *restrict commondata, const double *restrict f_local,
                                                     double *restrict Gamma_local) {
#include "set_CodeParameters.h"
  //==========================================
  // CONNECTION EVALUATION & THREAD-LOCAL UNPACKING
  //==========================================
  // Extract spatial coordinates $x^i$ and compute $Gamma^alpha_{mu nu}$.
  // Unpack position coordinates $x^i$ from the thread-local state vector.
  // Evaluated at compile time for state vector size: 9
  const double x = f_local[1];
  const double y = f_local[2];
  const double z = f_local[3];

  const REAL tmp1 = ((a_spin) * (a_spin));
  const REAL tmp2 = ((z) * (z));
  const REAL tmp21 = 2 * x;
  const REAL tmp30 = 2 * y;
  const REAL tmp34 = 2 * z;
  const REAL tmp6 = -tmp1 + tmp2 + ((x) * (x)) + ((y) * (y));
  const REAL tmp7 = sqrt(4 * tmp1 * tmp2 + ((tmp6) * (tmp6)));
  const REAL tmp8 = (1.0 / (tmp7));
  const REAL tmp11 = (1.0 / 2.0) * tmp2 + (1.0 / 2.0) * tmp7 + (1.0 / 2.0) * ((x) * (x)) + (1.0 / 2.0) * ((y) * (y));
  const REAL tmp9 = tmp6 * tmp8;
  const REAL tmp12 = -1.0 / 2.0 * tmp1 + tmp11;
  const REAL tmp35 = tmp8 * (4 * tmp1 * z + tmp34 * tmp6);
  const REAL tmp44 = (1.0 / 2.0) * tmp1 + tmp11;
  const REAL tmp13 = sqrt(tmp12);
  const REAL tmp22 = tmp21 * tmp9 + tmp21;
  const REAL tmp31 = tmp30 * tmp9 + tmp30;
  const REAL tmp37 = pow(tmp12, 3.0 / 2.0);
  const REAL tmp41 = tmp9 * x + x;
  const REAL tmp45 = (1.0 / ((tmp44) * (tmp44)));
  const REAL tmp50 = (1.0 / (tmp44));
  const REAL tmp53 = (1.0 / 2.0) * tmp9 * x + (1.0 / 2.0) * x;
  const REAL tmp60 = (1.0 / 2.0) * tmp9 * y + (1.0 / 2.0) * y;
  const REAL tmp63 = tmp9 * y + y;
  const REAL tmp67 = (1.0 / 4.0) * tmp35 + (1.0 / 2.0) * z;
  const REAL tmp70 = (1.0 / 2.0) * tmp35 + z;
  const REAL tmp14 = 2 * tmp13;
  const REAL tmp16 = tmp1 * tmp2 + ((tmp12) * (tmp12));
  const REAL tmp39 = -tmp1 * tmp34 - tmp12 * (tmp34 + tmp35);
  const REAL tmp43 = a_spin * y + tmp13 * x;
  const REAL tmp54 = (1.0 / (tmp13));
  const REAL tmp73 = -a_spin * x + tmp13 * y;
  const REAL tmp95 = -tmp34 - tmp35;
  const REAL tmp17 = (1.0 / (tmp16));
  const REAL tmp23 = (1.0 / ((tmp16) * (tmp16)));
  const REAL tmp51 = tmp43 * tmp50;
  const REAL tmp75 = tmp50 * tmp73;
  const REAL tmp76 = tmp13 + tmp54 * tmp60 * y;
  const REAL tmp92 = tmp43 * tmp45;
  const REAL tmp99 = tmp45 * tmp73;
  const REAL tmp108 = -a_spin + tmp53 * tmp54 * y;
  const REAL tmp18 = M_scale * tmp17;
  const REAL tmp25 = 2 * M_scale * tmp23;
  const REAL tmp46 = 2 * M_scale * tmp17;
  const REAL tmp56 = tmp13 + tmp53 * tmp54 * x;
  const REAL tmp61 = a_spin + tmp54 * tmp60 * x;
  const REAL tmp82 = M_scale * tmp23 * tmp34;
  const REAL tmp89 = ((tmp43) * (tmp43)) * tmp45;
  const REAL tmp97 = tmp73 * tmp92;
  const REAL tmp111 = tmp45 * ((tmp73) * (tmp73));
  const REAL tmp19 = tmp14 * tmp18;
  const REAL tmp26 = pow(tmp12, 5.0 / 2.0) * tmp25;
  const REAL tmp40 = tmp25 * tmp37 * tmp39;
  const REAL tmp47 = tmp37 * tmp46;
  const REAL tmp80 = tmp18 * tmp34;
  const REAL tmp83 = tmp12 * tmp39 * tmp82;
  const REAL tmp114 = ((tmp12) * (tmp12)) * tmp82;
  const REAL tmp20 = tmp19 * ((3.0 / 2.0) * tmp9 * x + (3.0 / 2.0) * x);
  const REAL tmp27 = tmp22 * tmp26;
  const REAL tmp29 = tmp19 * ((3.0 / 2.0) * tmp9 * y + (3.0 / 2.0) * y);
  const REAL tmp32 = tmp26 * tmp31;
  const REAL tmp36 = tmp19 * ((3.0 / 4.0) * tmp35 + (3.0 / 2.0) * z);
  const REAL tmp62 = tmp47 * tmp50;
  const REAL tmp69 = tmp12 * tmp18 * tmp50 * tmp67;
  const REAL tmp81 = tmp70 * tmp80;
  const REAL tmp87 = tmp47 / ((tmp44) * (tmp44) * (tmp44));
  const REAL tmp93 = tmp12 * tmp18 * tmp67;
  const REAL tmp102 = tmp12 * tmp46 * tmp50;
  const REAL tmp104 = tmp13 * tmp18 * tmp50 * tmp67 * z;
  const REAL tmp105 = tmp12 * tmp80;
  const REAL tmp116 = tmp63 * tmp80;
  const REAL tmp117 = tmp114 * tmp31;
  const REAL tmp119 = tmp12 * tmp50 * tmp80;
  const REAL tmp49 = tmp43 * tmp45 * tmp47;
  const REAL tmp74 = tmp45 * tmp47 * tmp73;
  const REAL tmp84 = tmp12 * tmp46 + tmp81 + tmp83;
  const REAL tmp88 = ((tmp43) * (tmp43)) * tmp87;
  const REAL tmp96 = tmp43 * tmp73 * tmp87;
  const REAL tmp106 = -tmp105 * tmp70;
  const REAL tmp110 = ((tmp73) * (tmp73)) * tmp87;
  const REAL tmp58 = tmp20 * tmp51 - tmp27 * tmp51 - tmp41 * tmp49 + tmp47 * tmp50 * tmp56;
  const REAL tmp65 = tmp29 * tmp51 - tmp32 * tmp51 - tmp49 * tmp63 + tmp61 * tmp62;
  const REAL tmp72 = tmp21 * tmp69 + tmp36 * tmp51 + tmp40 * tmp51 - tmp49 * tmp70;
  const REAL tmp77 = tmp29 * tmp75 - tmp32 * tmp75 + tmp62 * tmp76 - tmp63 * tmp74;
  const REAL tmp78 = tmp30 * tmp69 + tmp36 * tmp75 + tmp40 * tmp75 - tmp70 * tmp74;
  const REAL tmp98 = tmp29 * tmp97 - tmp31 * tmp96 - tmp32 * tmp97 + tmp49 * tmp76 + tmp61 * tmp74;
  const REAL tmp101 = tmp21 * tmp93 * tmp99 + tmp30 * tmp92 * tmp93 + tmp36 * tmp97 + tmp40 * tmp97 + tmp95 * tmp96;
  const REAL tmp107 = tmp102 * tmp43 + tmp104 * tmp21 + tmp106 * tmp92 + tmp51 * tmp81 + tmp51 * tmp83;
  const REAL tmp112 = tmp102 * tmp73 + tmp104 * tmp30 + tmp106 * tmp99 + tmp75 * tmp81 + tmp75 * tmp83;
  Gamma_local[0] = 0;
  Gamma_local[1] = tmp20 - tmp27;
  Gamma_local[2] = tmp29 - tmp32;
  Gamma_local[3] = tmp36 + tmp40;
  Gamma_local[4] = tmp58;
  Gamma_local[5] = tmp65;
  Gamma_local[6] = tmp72;
  Gamma_local[7] = tmp77;
  Gamma_local[8] = tmp78;
  Gamma_local[9] = tmp84;
  Gamma_local[10] = 0;
  Gamma_local[11] = tmp58;
  Gamma_local[12] = tmp65;
  Gamma_local[13] = tmp72;
  Gamma_local[14] = tmp20 * tmp89 - tmp22 * tmp88 - tmp27 * tmp89 + tmp49 * (tmp14 + tmp21 * tmp53 * tmp54);
  Gamma_local[15] = tmp29 * tmp89 - tmp31 * tmp88 - tmp32 * tmp89 + tmp49 * (2 * a_spin + tmp21 * tmp54 * tmp60);
  Gamma_local[16] = tmp36 * tmp89 + tmp40 * tmp89 + tmp88 * tmp95 + 4 * tmp92 * tmp93 * x;
  Gamma_local[17] = tmp98;
  Gamma_local[18] = tmp101;
  Gamma_local[19] = tmp107;
  Gamma_local[20] = 0;
  Gamma_local[21] = tmp108 * tmp62 + tmp20 * tmp75 - tmp27 * tmp75 - tmp41 * tmp74;
  Gamma_local[22] = tmp77;
  Gamma_local[23] = tmp78;
  Gamma_local[24] = tmp108 * tmp49 + tmp20 * tmp97 - tmp22 * tmp96 - tmp27 * tmp97 + tmp56 * tmp74;
  Gamma_local[25] = tmp98;
  Gamma_local[26] = tmp101;
  Gamma_local[27] = -tmp110 * tmp31 + tmp111 * tmp29 - tmp111 * tmp32 + tmp74 * (tmp14 + tmp30 * tmp54 * tmp60);
  Gamma_local[28] = tmp110 * tmp95 + tmp111 * tmp36 + tmp111 * tmp40 + 4 * tmp93 * tmp99 * y;
  Gamma_local[29] = tmp112;
  Gamma_local[30] = 0;
  Gamma_local[31] = -tmp114 * tmp22 + tmp41 * tmp80;
  Gamma_local[32] = tmp116 - tmp117;
  Gamma_local[33] = tmp84;
  Gamma_local[34] = -tmp105 * tmp41 * tmp92 + tmp105 * tmp50 * tmp56 - tmp114 * tmp22 * tmp51 + tmp41 * tmp51 * tmp80;
  Gamma_local[35] = -tmp105 * tmp63 * tmp92 + tmp116 * tmp51 - tmp117 * tmp51 + tmp119 * tmp61;
  Gamma_local[36] = tmp107;
  Gamma_local[37] = -tmp105 * tmp63 * tmp99 + tmp116 * tmp75 - tmp117 * tmp75 + tmp119 * tmp76;
  Gamma_local[38] = tmp112;
  Gamma_local[39] = M_scale * tmp14 * tmp2 * tmp23 * tmp39 + 4 * tmp13 * tmp18 * z + tmp2 * tmp46 * tmp54 * tmp67;

} // END FUNCTION: connections_KerrSchild_Cartesian

/**
 * Interpolate piecewise static/numerical spacetime tensors for one photon chunk.
 *
 * The caller supplies an active numerical time window, a spatial interpolation
 * context, and one chunk of photon states in the same Structure-of-Arrays bundle
 * layout used by the analytic geodesic interpolation kernel. This CPU wrapper
 * parallelizes over rays, spatially interpolates the first selected numerical slice
 * for times at or below `t_numerical_initial`, freezes the numerical spacetime to the
 * final selected slice for times at or above that final slice time, and otherwise
 * reconstructs one mixed temporal stencil. Lower missing stencil nodes use a
 * fresh spatial interpolation of the first selected numerical slice, whereas upper
 * missing nodes reuse the final selected numerical slice. Every reconstructed
 * stencil uses ordinary temporal interpolation.
 *
 * The design goal is to let all photons in the chunk reuse the same mapped
 * numerical-spacetime payload window rather than loading numerical grids
 * independently ray-by-ray.
 *
 * @param[in] commondata Common runtime parameters.
 * @param[in] params Generated BHaH grid parameters for the mapped numerical data.
 * @param[in] spatial_context Trusted azimuthal-symmetry spatial interpolation context.
 * @param[in] numerical_window Active mapped numerical time-window manager.
 * @param[in] d_f_bundle Photon state bundle.
 * @param[out] d_metric_bundle Destination metric bundle.
 * @param[out] d_metric_derivative_bundle Destination metric-derivative bundle, or NULL.
 * @param chunk_size Number of active rays in the chunk.
 * @param stream_idx Analytic-kernel compatibility argument; ignored on CPU.
 *
 * @pre `numerical_window` must already map all slices needed by the active chunk.
 */
void numerical_interpolation(const commondata_struct *restrict commondata, const params_struct *restrict params,
                             const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context,
                             const NumericalTimeWindowManager *restrict numerical_window, const double *restrict d_f_bundle,

                             double *restrict d_metric_bundle, double *restrict d_metric_derivative_bundle, const long int chunk_size,
                             const int stream_idx) {
  (void)stream_idx;

#define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_METRIC_DERIVATIVE(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define G4_SLICE(base_ptr, slice_idx, comp_idx) ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT + (comp_idx)])
#define METRIC_DERIVATIVE_SLICE(base_ptr, slice_idx, comp_idx)                                                                                       \
  ((base_ptr)[(slice_idx) * TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT + (comp_idx)])

  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  const REAL t_numerical_initial = (REAL)commondata->t_numerical_initial;
  const REAL t_numerical_end = (REAL)commondata->t_numerical_end;
  const REAL dt_numerical_spacetime_data = (REAL)commondata->dt_numerical_spacetime_data;
  // Step 0: Validate the mapped numerical window pointers before reading any
  // of their fields.
  if (numerical_window == NULL || numerical_window->slice_times == NULL || numerical_window->num_time_slices < 1ULL ||
      numerical_window->time_slice_stride == 0ULL) {
#pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_metric_derivative_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++)
          d_metric_derivative_bundle[IDX_METRIC_DERIVATIVE(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid mapped numerical window metadata
    return;
  } // END IF: mapped numerical window metadata was unavailable
  const uint64_t first_slice_index = 0ULL;
  const uint64_t final_slice_index = time_window_manager_numerical_final_selected_slice(numerical_window);
  const REAL t_final_numerical_slice = (REAL)numerical_window->slice_times[final_slice_index];
  const REAL selected_slice_dt = dt_numerical_spacetime_data * (REAL)numerical_window->time_slice_stride;
  // The user-facing t_numerical_end must identify the final selected numerical
  // slice to within two selected time steps.
  const REAL t_numerical_end_tolerance = 2.0 * selected_slice_dt;
  // The mapped numerical window and the temporal helper must agree on the
  // centered temporal stencil width before any variable-length arrays are
  // sized from runtime data.
  if (temporal_half_width < 0 || temporal_half_width > TEMPORAL_LAGRANGE_INTERP_MAX_HALF_WIDTH ||
      temporal_half_width != numerical_window->temporal_interp_half_width || !isfinite((double)t_numerical_initial) ||
      !isfinite((double)t_numerical_end) || !isfinite((double)dt_numerical_spacetime_data) || dt_numerical_spacetime_data <= 0.0 ||
      !isfinite((double)selected_slice_dt) || selected_slice_dt <= 0.0 || !isfinite((double)t_final_numerical_slice) ||
      t_numerical_initial >= t_numerical_end) {
#pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_metric_derivative_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++)
          d_metric_derivative_bundle[IDX_METRIC_DERIVATIVE(comp, i)] = NAN;
    } // END LOOP: for i over rays after invalid temporal bounds or stencil metadata
    return;
  } // END IF: temporal bounds or mapped stencil metadata were invalid
  if (fabs((double)(t_numerical_end - t_final_numerical_slice)) > (double)t_numerical_end_tolerance) {
    fprintf(stderr,
            "ERROR: commondata->t_numerical_end=%e differs from the final selected numerical slice time=%e "
            "by more than 2*selected_slice_dt=%e.\n",
            (double)t_numerical_end, (double)t_final_numerical_slice, (double)t_numerical_end_tolerance);
    exit(1);
  } // END IF: t_numerical_end was inconsistent with the final selected numerical slice time
  const int temporal_num_points = 2 * temporal_half_width + 1;
  if (temporal_num_points != numerical_window->temporal_interp_num_points) {
#pragma omp parallel for
    for (long int i = 0; i < chunk_size; i++) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_metric_derivative_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++)
          d_metric_derivative_bundle[IDX_METRIC_DERIVATIVE(comp, i)] = NAN;
    } // END LOOP: for i over rays after inconsistent temporal stencil size
    return;
  } // END IF: runtime temporal stencil size did not match the mapped numerical window

#pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    double f_local[9];
    for (int comp = 0; comp < 9; comp++)
      f_local[comp] = d_f_bundle[IDX_F(comp, i)];
    const REAL t = (REAL)f_local[0];
    const REAL x = (REAL)f_local[1];
    const REAL y = (REAL)f_local[2];
    const REAL z = (REAL)f_local[3];
    int ray_failed = 0;
    uint64_t available_slice_indices[temporal_num_points];
    REAL available_slice_times[temporal_num_points];
    const double *available_slice_payloads[temporal_num_points];
    int num_available_slices = 0;
    REAL missing_slice_times[temporal_num_points];
    int num_missing_slices = 0;
    REAL g4dd_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL g4dd_dD_available[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT];
    REAL full_slice_times[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL g4dd_dD_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL g4dd_dD_local[TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT];
    REAL g4dd_missing_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT] = {0};
    REAL g4dd_dD_missing_local[TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT] = {0};

    // Step 1: Dispatch directly to a static numerical endpoint when this ray
    // is outside the mixed temporal-interpolation region. At and below
    // t_numerical_initial, interpolate the first stored slice in space only, so its
    // metric time derivatives remain zero. At or above the final selected
    // slice time, likewise use the final slice in space only.
    if (t <= t_numerical_initial) {
      const double *first_slice_payloads[1];
      first_slice_payloads[0] = time_window_manager_numerical_grid_ptr(numerical_window, first_slice_index);
      if (first_slice_payloads[0] == NULL) {
        ray_failed = 1;
      } else {
        const int spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(spatial_context, commondata, params, x, y, z, 1,
                                                                                                     first_slice_payloads, g4dd_local, g4dd_dD_local);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
      } // END ELSE: first selected numerical slice was available for direct spatial interpolation
    } else if (t >= t_final_numerical_slice) {
      const double *final_slice_payloads[1];
      final_slice_payloads[0] = time_window_manager_numerical_grid_ptr(numerical_window, final_slice_index);
      if (final_slice_payloads[0] == NULL) {
        ray_failed = 1;
      } else {
        const int spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(spatial_context, commondata, params, x, y, z, 1,
                                                                                                     final_slice_payloads, g4dd_local, g4dd_dD_local);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
      } // END ELSE: final selected numerical slice was available for direct spatial interpolation
    } else {
      // Step 2: Recover one adaptive numerical stencil from the slot-level
      // time window shared by the whole chunk.
      const int window_status = time_window_manager_numerical_stencil_for_time(
          numerical_window, (double)t, temporal_half_width, available_slice_indices, available_slice_times, available_slice_payloads,
          &num_available_slices, missing_slice_times, &num_missing_slices);
      if (window_status != TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS || num_available_slices <= 0 ||
          num_available_slices + num_missing_slices != temporal_num_points) {
        ray_failed = 1;
      } else {
        // Step 3: Interpolate only the available mapped numerical slices in
        // space at the photon position.
        const int spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
            spatial_context, commondata, params, x, y, z, num_available_slices, available_slice_payloads, g4dd_available, g4dd_dD_available);
        if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
          ray_failed = 1;
        else {
          int missing_is_upper = 0;
          int missing_is_lower = 0;

          // Step 4: Centered lower-boundary stencils pad their missing negative
          // time nodes from the static first numerical slice. Freeze the final
          // selected numerical slice for ordinary upper missing nodes.
          if (num_missing_slices > 0) {
            missing_is_upper = missing_slice_times[0] > available_slice_times[num_available_slices - 1];
            missing_is_lower = missing_slice_times[num_missing_slices - 1] < available_slice_times[0];

            if (missing_is_upper == missing_is_lower) {
              ray_failed = 1;
            } else {
              if (missing_is_upper) {
                const int last_available_slot = num_available_slices - 1;
                if (available_slice_indices[last_available_slot] != final_slice_index) {
                  ray_failed = 1;
                } else {
                  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                    g4dd_missing_local[comp] = G4_SLICE(g4dd_available, last_available_slot, comp);
                  } // END LOOP: for comp over upper frozen metric components
                  for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++) {
                    g4dd_dD_missing_local[comp] = METRIC_DERIVATIVE_SLICE(g4dd_dD_available, last_available_slot, comp);
                  } // END LOOP: for comp over upper frozen metric-derivative components
                } // END ELSE: upper missing stencil edge reached the final selected numerical slice
              } else {
                const double *first_slice_payloads[1];
                first_slice_payloads[0] = time_window_manager_numerical_grid_ptr(numerical_window, first_slice_index);
                if (first_slice_payloads[0] == NULL) {
                  ray_failed = 1;
                } else {
                  const int first_slice_spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
                      spatial_context, commondata, params, x, y, z, 1, first_slice_payloads, g4dd_missing_local, g4dd_dD_missing_local);
                  if (first_slice_spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS)
                    ray_failed = 1;
                } // END ELSE: first selected numerical slice was available for lower stencil padding
              } // END ELSE: lower missing stencil edge uses the first selected numerical slice
            } // END ELSE: missing stencil edge classification was usable
          } // END IF: one stencil edge must be synthesized

          if (!ray_failed) {
            // Step 5: Reconstruct the full ordered temporal stencil expected
            // by the temporal interpolation helper.
            if (num_missing_slices > 0 && missing_is_lower) {
              for (int s = 0; s < num_missing_slices; s++) {
                full_slice_times[s] = missing_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over missing first-slice metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++) {
                  METRIC_DERIVATIVE_SLICE(g4dd_dD_slices, s, comp) = g4dd_dD_missing_local[comp];
                } // END LOOP: for comp over missing first-slice metric derivatives
              } // END LOOP: for s over lower missing static first-slice nodes
              for (int s = 0; s < num_available_slices; s++) {
                const int full_slot = num_missing_slices + s;
                full_slice_times[full_slot] = available_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) = G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components after lower first-slice padding
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++) {
                  METRIC_DERIVATIVE_SLICE(g4dd_dD_slices, full_slot, comp) = METRIC_DERIVATIVE_SLICE(g4dd_dD_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric derivatives after lower first-slice padding
              } // END LOOP: for s over available numerical stencil nodes after lower first-slice padding
            } else {
              for (int s = 0; s < num_available_slices; s++) {
                full_slice_times[s] = available_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, s, comp) = G4_SLICE(g4dd_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++) {
                  METRIC_DERIVATIVE_SLICE(g4dd_dD_slices, s, comp) = METRIC_DERIVATIVE_SLICE(g4dd_dD_available, s, comp);
                } // END LOOP: for comp over mapped numerical metric-derivative components
              } // END LOOP: for s over available numerical stencil nodes
              for (int s = 0; s < num_missing_slices; s++) {
                const int full_slot = num_available_slices + s;
                full_slice_times[full_slot] = missing_slice_times[s];
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
                  G4_SLICE(g4dd_slices, full_slot, comp) = g4dd_missing_local[comp];
                } // END LOOP: for comp over upper missing frozen metric components
                for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++) {
                  METRIC_DERIVATIVE_SLICE(g4dd_dD_slices, full_slot, comp) = g4dd_dD_missing_local[comp];
                } // END LOOP: for comp over upper missing frozen metric-derivative components
              } // END LOOP: for s over upper missing stencil nodes
            } // END ELSE: upper missing edge or fully numerical centered stencil
          } // END IF: full ordered stencil reconstruction remained valid

          if (!ray_failed) {
            // Step 6: Interpolate the reconstructed stencil in physical time.
            const int temporal_status =
                temporal_lagrange_interpolation(commondata, full_slice_times, g4dd_slices, g4dd_dD_slices, t, g4dd_local, g4dd_dD_local);
            if (temporal_status != TEMPORAL_LAGRANGE_INTERP_SUCCESS)
              ray_failed = 1;
          } // END IF: reconstructed stencil was ready for temporal interpolation
        } // END ELSE: spatial interpolation succeeded for the mapped numerical stencil subset
      } // END ELSE: adaptive stencil query succeeded for this mixed ray
    } // END ELSE: photon required numerical or mixed interpolation

    if (ray_failed) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      if (d_metric_derivative_bundle != NULL)
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++)
          d_metric_derivative_bundle[IDX_METRIC_DERIVATIVE(comp, i)] = NAN;
      continue;
    } // END IF: at least one interpolation stage failed for this ray

    for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++)
      d_metric_bundle[IDX_METRIC(comp, i)] = (double)g4dd_local[comp];
    if (d_metric_derivative_bundle != NULL)
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_METRIC_DERIVATIVE_COMPONENT_COUNT; comp++)
        d_metric_derivative_bundle[IDX_METRIC_DERIVATIVE(comp, i)] = (double)g4dd_dD_local[comp];
  } // END LOOP: for i over rays in chunk

#undef IDX_F
#undef IDX_METRIC
#undef IDX_METRIC_DERIVATIVE
#undef G4_SLICE
#undef METRIC_DERIVATIVE_SLICE
} // END FUNCTION: numerical_interpolation

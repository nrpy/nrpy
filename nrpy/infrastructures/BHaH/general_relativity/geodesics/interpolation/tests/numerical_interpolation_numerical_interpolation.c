#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#include <math.h>
#include <stdint.h>

/**
 * Interpolate numerical-spacetime tensors for one photon chunk.
 *
 * The caller supplies an active numerical time window, a spatial interpolation
 * context, and one chunk of photon states in the same Structure-of-Arrays bundle
 * layout used by the analytic geodesic interpolation kernel. This CPU wrapper
 * parallelizes over rays, selects each photon's mapped temporal stencil, performs
 * spatial interpolation on every stencil slice, performs temporal interpolation
 * at the photon coordinate time, and writes the final metric and optional
 * Christoffel bundles.
 *
 * @param[in] commondata Common runtime parameters.
 * @param[in] params Generated BHaH grid parameters for the mapped numerical data.
 * @param[in] spatial_context Trusted azimuthal-symmetry spatial interpolation context.
 * @param[in] numerical_window Active mapped numerical time-window manager.
 * @param[in] d_f_bundle Photon state bundle.
 * @param[out] d_metric_bundle Destination metric bundle.
 * @param[out] d_connection_bundle Destination Christoffel bundle, or NULL.
 * @param chunk_size Number of active rays in the chunk.
 * @param stream_idx Analytic-kernel compatibility argument; ignored on CPU.
 *
 * @pre `numerical_window` must already map all slices needed by the active chunk.
 */
void numerical_interpolation(const commondata_struct *restrict commondata, const params_struct *restrict params,
                             const azimuthal_symmetry_spatial_lagrange_context_struct *restrict spatial_context,
                             const NumericalTimeWindowManager *restrict numerical_window, const double *restrict d_f_bundle,
                             double *restrict d_metric_bundle, double *restrict d_connection_bundle, const long int chunk_size,
                             const int stream_idx) {
  (void)stream_idx;

#define IDX_F(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_METRIC(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))
#define IDX_CONN(c, ray_id) ((c) * BUNDLE_CAPACITY + (ray_id))

  const int temporal_half_width = commondata->numerical_spacetime_temporal_interp_order;
  const int temporal_num_points = 2 * temporal_half_width + 1;

#pragma omp parallel for
  for (long int i = 0; i < chunk_size; i++) {
    const REAL t = (REAL)d_f_bundle[IDX_F(0, i)];
    const REAL x = (REAL)d_f_bundle[IDX_F(1, i)];
    const REAL y = (REAL)d_f_bundle[IDX_F(2, i)];
    const REAL z = (REAL)d_f_bundle[IDX_F(3, i)];
    uint64_t slice_indices[temporal_num_points];
    REAL slice_times[temporal_num_points];
    const double *slice_payloads[temporal_num_points];
    REAL g4dd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_slices[temporal_num_points * TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];
    REAL g4dd_local[TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT];
    REAL gamma4udd_local[TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT];

    const int window_status =
        numerical_time_window_manager_stencil_for_time(numerical_window, (double)t, temporal_half_width, slice_indices, slice_times, slice_payloads);
    if (window_status != NUMERICAL_TIME_WINDOW_MANAGER_SUCCESS) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric failure outputs
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection failure outputs
      } // END IF: connection output bundle was requested
      continue;
    } // END IF: temporal stencil was not inside the mapped numerical window

    const int spatial_status = azimuthal_symmetry_spatial_lagrange_interpolation__rfm__Spherical(
        spatial_context, commondata, params, x, y, z, temporal_num_points, slice_payloads, g4dd_slices, gamma4udd_slices);
    if (spatial_status != AZIMUTHAL_SYMMETRY_SPATIAL_LAGRANGE_INTERP_SUCCESS) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
        d_metric_bundle[IDX_METRIC(comp, i)] = NAN;
      } // END LOOP: for comp over metric failure outputs
      if (d_connection_bundle != NULL) {
        for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
          d_connection_bundle[IDX_CONN(comp, i)] = NAN;
        } // END LOOP: for comp over connection failure outputs
      } // END IF: connection output bundle was requested
      continue;
    } // END IF: spatial interpolation failed for this ray

    temporal_lagrange_interpolation(commondata, slice_times, g4dd_slices, gamma4udd_slices, t, g4dd_local, gamma4udd_local);

    for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_G4_COMPONENT_COUNT; comp++) {
      d_metric_bundle[IDX_METRIC(comp, i)] = (double)g4dd_local[comp];
    } // END LOOP: for comp over final metric components
    if (d_connection_bundle != NULL) {
      for (int comp = 0; comp < TEMPORAL_LAGRANGE_INTERP_GAMMA_COMPONENT_COUNT; comp++) {
        d_connection_bundle[IDX_CONN(comp, i)] = (double)gamma4udd_local[comp];
      } // END LOOP: for comp over final Christoffel components
    } // END IF: connection output bundle was requested
  } // END LOOP: for i over rays in chunk

#undef IDX_F
#undef IDX_METRIC
#undef IDX_CONN
} // END FUNCTION: numerical_interpolation

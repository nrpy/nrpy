"""
Module for generating Cartesian initial conditions for the photon integration pipeline.

This module defines the C function responsible for initializing the spatial
coordinates and spatial momentum of photons on a camera window.
"""

import nrpy.c_function as cfc
import nrpy.params as par


def set_initial_conditions_cartesian(spacetime_name: str) -> None:
    """
    Register the C function to set up Cartesian initial conditions for photons.

    :param spacetime_name: The specific metric or spacetime identifier (e.g., 'Kerr').

    Examples
        >>> set_initial_conditions_cartesian("Kerr")
    """
    par.register_CodeParameter(
        "int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True
    )
    par.register_CodeParameter(
        "REAL", __name__, "t_start", 100, commondata=True, add_to_parfile=True
    )

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "math.h",
        "stdio.h",
        "stdlib.h",
        "omp.h",
    ]

    desc = f"""@brief Initializes Cartesian starting conditions for photons in {spacetime_name}.
    Detailed algorithm: Computes the geometric basis vectors of the camera window,
    maps individual photon rays to pixel coordinates, and populates the initial
    flattened Structure of Arrays (SoA) state vector."""

    name = f"set_initial_conditions_cartesian_{spacetime_name}"

    params = (
        "const commondata_struct *restrict commondata, "
        "long int num_rays, "
        "PhotonStateSoA *restrict all_photons, "
        "double window_center_out[3], "
        "double n_x_out[3], "
        "double n_y_out[3], "
        "double n_z_out[3]"
    )

    body = f"""
    // Spatial coordinates of the observer's camera in the global Cartesian frame.
    const double camera_pos[3] = {{commondata->camera_pos_x, commondata->camera_pos_y, commondata->camera_pos_z}};
    // Geometric center of the projection window in the global Cartesian frame.
    const double window_center[3] = {{commondata->window_center_x, commondata->window_center_y, commondata->window_center_z}};

    // Vector normal to the camera window (line of sight direction).
    double n_z[3] = {{window_center[0] - camera_pos[0], window_center[1] - camera_pos[1], window_center[2] - camera_pos[2]}};
    // Magnitude of the line of sight vector used for normalization.
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]);
    for(int i=0; i<3; i++) n_z[i] /= mag_n_z;

    // User-defined upward orientation vector for the camera frame.
    const double guide_up[3] = {{commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}};

    // Orthonormal basis vector describing the horizontal axis of the camera frame.
    double n_x[3] = {{n_z[1]*guide_up[2] - n_z[2]*guide_up[1], n_z[2]*guide_up[0] - n_z[0]*guide_up[2], n_z[0]*guide_up[1] - n_z[1]*guide_up[0]}};
    // Magnitude of the horizontal basis vector used for normalization.
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);

    // Fallback mechanism to prevent cross-product singularities at the poles.
    if (mag_n_x < 1e-9) {{
        // Alternative upward direction if the primary guide aligns with the line of sight.
        double alternative_up[3] = {{0.0, 1.0, 0.0}};
        if (fabs(n_z[1]) > 0.999) {{ alternative_up[1] = 0.0; alternative_up[2] = 1.0; }}
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    }}

    for(int i=0; i<3; i++) n_x[i] /= mag_n_x;

    // Orthonormal basis vector describing the vertical axis of the camera frame.
    double n_y[3] = {{n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]}};

    for(int i=0; i<3; i++) {{
        window_center_out[i] = window_center[i];
        n_x_out[i] = n_x[i]; n_y_out[i] = n_y[i]; n_z_out[i] = n_z[i];
    }}

    #ifdef USE_GPU
        // Device-allocated compacted list of photon indices mapping back to global memory.
        int *req_photon_ids = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        // Device-allocated temporary positional data array for interpolation mapping.
        double *req_pos = (double *)omp_target_alloc(sizeof(double) * 4 * BUNDLE_CAPACITY, omp_get_default_device());
        // Device-allocated storage for the 10 independent components of the symmetric 4D metric tensor.
        double *metric_g4DD = (double *)omp_target_alloc(sizeof(double) * 10 * BUNDLE_CAPACITY, omp_get_default_device());
    #else
        // Host-allocated compacted list of photon indices mapping back to global memory.
        int *req_photon_ids = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
        // Host-allocated temporary positional data array for interpolation mapping.
        double *req_pos = (double *)malloc(sizeof(double) * 4 * BUNDLE_CAPACITY);
        // Host-allocated storage for the 10 independent components of the symmetric 4D metric tensor.
        double *metric_g4DD = (double *)malloc(sizeof(double) * 10 * BUNDLE_CAPACITY);
    #endif

    if (!req_photon_ids || !req_pos || !metric_g4DD) {{
        // Ensures critical memory allocation succeeds before proceeding.
        fprintf(stderr, "FATAL: Failed to allocate flattened initialization batch arrays.\\n");
        exit(1);
    }}

    #ifdef USE_GPU
    #pragma omp target data map(to: commondata[0:1], camera_pos[0:3], window_center[0:3], n_x[0:3], n_y[0:3], n_z[0:3])
    #endif
    {{
        for (long int start_idx = 0; start_idx < num_rays; start_idx += BUNDLE_CAPACITY) {{
            // Determine the valid number of threads for the current batch.
            long int current_chunk_size = (num_rays - start_idx < BUNDLE_CAPACITY) ? (num_rays - start_idx) : BUNDLE_CAPACITY;

            #ifdef USE_GPU
                #pragma omp target teams distribute parallel for map(to: all_photons) has_device_addr(all_photons->f) is_device_ptr(req_photon_ids, req_pos)
            #else
                #pragma omp parallel for
            #endif
            for (long int c = 0; c < current_chunk_size; c++) {{
                // The global index of the photon currently being initialized.
                long int i = start_idx + c;

                // The vertical pixel coordinate index derived from the 1D global ray identifier.
                const int row = i / commondata->scan_density;
                // The horizontal pixel coordinate index derived from the 1D global ray identifier.
                const int col = i % commondata->scan_density;

                // Local physical distance along the camera frame's horizontal axis.
                const double x_pix = -commondata->window_width/2.0 + (col + 0.5) * (commondata->window_width / commondata->scan_density);
                // Local physical distance along the camera frame's vertical axis.
                const double y_pix = -commondata->window_height/2.0 + (row + 0.5) * (commondata->window_height / commondata->scan_density);

                // Global Cartesian intersection point on the projection window.
                const double target_pos[3] = {{
                    window_center[0] + x_pix*n_x[0] + y_pix*n_y[0],
                    window_center[1] + x_pix*n_x[1] + y_pix*n_y[1],
                    window_center[2] + x_pix*n_x[2] + y_pix*n_y[2]
                }};

                all_photons->f[IDX_GLOBAL(0, i, num_rays)] = commondata->t_start;
                all_photons->f[IDX_GLOBAL(1, i, num_rays)] = camera_pos[0];
                all_photons->f[IDX_GLOBAL(2, i, num_rays)] = camera_pos[1];
                all_photons->f[IDX_GLOBAL(3, i, num_rays)] = camera_pos[2];

                // Unnormalized geometric trajectory vector connecting the camera to the target pixel.
                const double V_x = target_pos[0] - camera_pos[0];
                const double V_y = target_pos[1] - camera_pos[1];
                const double V_z = target_pos[2] - camera_pos[2];

                // Normalization scalar to generate a unit direction vector.
                const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y + V_z*V_z);

                all_photons->f[IDX_GLOBAL(5, i, num_rays)] = V_x * inv_mag_V;
                all_photons->f[IDX_GLOBAL(6, i, num_rays)] = V_y * inv_mag_V;
                all_photons->f[IDX_GLOBAL(7, i, num_rays)] = V_z * inv_mag_V;
                all_photons->f[IDX_GLOBAL(8, i, num_rays)] = 0.0; // Affine parameter

                req_photon_ids[c] = (int)i;
                for(int m=0; m<4; m++) {{
                    req_pos[IDX_LOCAL(m, c, BUNDLE_CAPACITY)] = all_photons->f[IDX_GLOBAL(m, i, num_rays)];
                }}
            }}

            placeholder_interpolation_engine_{spacetime_name}(commondata, (int)current_chunk_size, req_photon_ids, req_pos, metric_g4DD, NULL);

            #ifdef USE_GPU
                #pragma omp target teams distribute parallel for map(to: all_photons) has_device_addr(all_photons->f) is_device_ptr(metric_g4DD)
            #else
                #pragma omp parallel for
            #endif
            for (long int c = 0; c < current_chunk_size; c++) {{
                // The global index of the photon for p0 constraint solving.
                long int i = start_idx + c;
                p0_reverse(metric_g4DD, all_photons->f, num_rays, BUNDLE_CAPACITY, (int)i, (int)c, &all_photons->f[IDX_GLOBAL(4, i, num_rays)]);
            }}
        }}
    }}

    #ifdef USE_GPU
        omp_target_free(req_photon_ids, omp_get_default_device());
        omp_target_free(req_pos, omp_get_default_device());
        omp_target_free(metric_g4DD, omp_get_default_device());
    #else
        free(req_photon_ids); free(req_pos); free(metric_g4DD);
    #endif
    """

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
    )

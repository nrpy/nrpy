# set_initial_conditions_cartesian.py
import nrpy.c_function as cfc
import nrpy.params as par

def set_initial_conditions_cartesian(spacetime_name: str) -> None:
    par.register_CodeParameter("int", __name__, "scan_density", 500, commondata=True, add_to_parfile=True)
    par.register_CodeParameter("REAL", __name__, "t_start", 100, commondata=True, add_to_parfile=True)

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "math.h", "stdio.h", "stdlib.h", "omp.h"]
    desc = f"Sets up Cartesian initial conditions for photons in the {spacetime_name} spacetime."
    name = f"set_initial_conditions_cartesian_{spacetime_name}"
    params = "const commondata_struct *restrict commondata, long int num_rays, PhotonStateSoA *restrict all_photons, double window_center_out[3], double n_x_out[3], double n_y_out[3], double n_z_out[3]"

    body = f"""
    const double camera_pos[3] = {{commondata->camera_pos_x, commondata->camera_pos_y, commondata->camera_pos_z}};
    const double window_center[3] = {{commondata->window_center_x, commondata->window_center_y, commondata->window_center_z}};

    // Construct basis vectors on host CPU
    double n_z[3] = {{window_center[0] - camera_pos[0], window_center[1] - camera_pos[1], window_center[2] - camera_pos[2]}};
    double mag_n_z = sqrt(n_z[0]*n_z[0] + n_z[1]*n_z[1] + n_z[2]*n_z[2]);
    for(int i=0; i<3; i++) n_z[i] /= mag_n_z;

    const double guide_up[3] = {{commondata->window_up_vec_x, commondata->window_up_vec_y, commondata->window_up_vec_z}};
    double n_x[3] = {{n_z[1]*guide_up[2] - n_z[2]*guide_up[1], n_z[2]*guide_up[0] - n_z[0]*guide_up[2], n_z[0]*guide_up[1] - n_z[1]*guide_up[0]}};
    double mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    
    // Pole Singularity Protection
    if (mag_n_x < 1e-9) {{
        double alternative_up[3] = {{0.0, 1.0, 0.0}};
        if (fabs(n_z[1]) > 0.999) {{ alternative_up[1] = 0.0; alternative_up[2] = 1.0; }}
        n_x[0] = alternative_up[1]*n_z[2] - alternative_up[2]*n_z[1];
        n_x[1] = alternative_up[2]*n_z[0] - alternative_up[0]*n_z[2];
        n_x[2] = alternative_up[0]*n_z[1] - alternative_up[1]*n_z[0];
        mag_n_x = sqrt(n_x[0]*n_x[0] + n_x[1]*n_x[1] + n_x[2]*n_x[2]);
    }}

    for(int i=0; i<3; i++) n_x[i] /= mag_n_x;
    double n_y[3] = {{n_z[1]*n_x[2] - n_z[2]*n_x[1], n_z[2]*n_x[0] - n_z[0]*n_x[2], n_z[0]*n_x[1] - n_z[1]*n_x[0]}};

    for(int i=0; i<3; i++) {{
        window_center_out[i] = window_center[i];
        n_x_out[i] = n_x[i]; n_y_out[i] = n_y[i]; n_z_out[i] = n_z[i];
    }}

    #ifdef USE_GPU
        int *req_photon_ids = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        double *req_pos = (double *)omp_target_alloc(sizeof(double) * 4 * BUNDLE_CAPACITY, omp_get_default_device());
        double *metric_g4DD = (double *)omp_target_alloc(sizeof(double) * 10 * BUNDLE_CAPACITY, omp_get_default_device()); 
    #else
        int *req_photon_ids = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
        double *req_pos = (double *)malloc(sizeof(double) * 4 * BUNDLE_CAPACITY);
        double *metric_g4DD = (double *)malloc(sizeof(double) * 10 * BUNDLE_CAPACITY); 
    #endif

    // Allocation Safety Check
    if (!req_photon_ids || !req_pos || !metric_g4DD) {{ 
        fprintf(stderr, "Error: Failed to allocate flattened initialization batch arrays.\\n"); 
        exit(1); 
    }}

    for (long int start_idx = 0; start_idx < num_rays; start_idx += BUNDLE_CAPACITY) {{
        long int current_chunk_size = (num_rays - start_idx < BUNDLE_CAPACITY) ? (num_rays - start_idx) : BUNDLE_CAPACITY;

        #ifdef USE_GPU
            #pragma omp target teams distribute parallel for \
                        map(to: all_photons, commondata[0:1]) \
                        map(to: camera_pos[0:3], window_center[0:3], n_x[0:3], n_y[0:3], n_z[0:3]) \
                        has_device_addr(all_photons.f) is_device_ptr(req_photon_ids, req_pos)
        #else
            #pragma omp parallel for
        #endif
        for (long int c = 0; c < current_chunk_size; c++) {{
            long int i = start_idx + c;
            const int row = i / commondata->scan_density;
            const int col = i % commondata->scan_density;
            const double x_pix = -commondata->window_width/2.0 + (col + 0.5) * (commondata->window_width / commondata->scan_density);
            const double y_pix = -commondata->window_height/2.0 + (row + 0.5) * (commondata->window_height / commondata->scan_density);
            const double target_pos[3] = {{
                window_center[0] + x_pix*n_x[0] + y_pix*n_y[0],
                window_center[1] + x_pix*n_x[1] + y_pix*n_y[1],
                window_center[2] + x_pix*n_x[2] + y_pix*n_y[2]
            }};
            all_photons->f[IDX_GLOBAL(0, i, num_rays)] = commondata->t_start;       
            all_photons->f[IDX_GLOBAL(1, i, num_rays)] = camera_pos[0]; 
            all_photons->f[IDX_GLOBAL(2, i, num_rays)] = camera_pos[1]; 
            all_photons->f[IDX_GLOBAL(3, i, num_rays)] = camera_pos[2]; 
            const double V_x = target_pos[0] - camera_pos[0], V_y = target_pos[1] - camera_pos[1], V_z = target_pos[2] - camera_pos[2];
            const double inv_mag_V = 1.0 / sqrt(V_x*V_x + V_y*V_y+ V_z*V_z);
            all_photons->f[IDX_GLOBAL(5, i, num_rays)] = V_x * inv_mag_V; 
            all_photons->f[IDX_GLOBAL(6, i, num_rays)] = V_y * inv_mag_V; 
            all_photons->f[IDX_GLOBAL(7, i, num_rays)] = V_z * inv_mag_V; 
            all_photons->f[IDX_GLOBAL(8, i, num_rays)] = 0.0;
            req_photon_ids[c] = (int)i;
            for(int m=0; m<4; m++) req_pos[IDX_LOCAL(m, c, BUNDLE_CAPACITY)] = all_photons->f[IDX_GLOBAL(m, i, num_rays)];
        }}

        placeholder_interpolation_engine_{spacetime_name}(commondata, (int)current_chunk_size, req_photon_ids, req_pos, metric_g4DD, NULL);

        #ifdef USE_GPU
            #pragma omp target teams distribute parallel for map(to: all_photons) has_device_addr(all_photons.f) is_device_ptr(metric_g4DD)
        #else
            #pragma omp parallel for
        #endif
        for (long int c = 0; c < current_chunk_size; c++) {{
            long int i = start_idx + c;
            p0_reverse(metric_g4DD, all_photons->f, num_rays, BUNDLE_CAPACITY, (int)i, (int)c, &all_photons->f[IDX_GLOBAL(4, i, num_rays)]);
        }}
    }}

    #ifdef USE_GPU
        omp_target_free(req_photon_ids, omp_get_default_device()); omp_target_free(req_pos, omp_get_default_device()); omp_target_free(metric_g4DD, omp_get_default_device());
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
        include_CodeParameters_h=False
    )
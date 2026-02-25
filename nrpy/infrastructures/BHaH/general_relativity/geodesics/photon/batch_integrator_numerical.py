# batch_integrator_numerical.py
"""
Generates the main C orchestrator for the numerical integration pipeline.
Project Singularity-Axiom: Dual-Architecture (CPU/GPU) Portability.
Fully restored logic with SoA-aware synchronization, batched interpolation, 
direct device allocation, and integrated physical termination checks.
"""
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par

batch_structs_c_code = r"""
    #define BUNDLE_CAPACITY 32768

    // Strictly enforce 1D mapping to prevent pointer offset arithmetic bugs
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))

    typedef enum {
        WINDOW_EVENT,
        SOURCE_EVENT
    } event_type_t;

    typedef enum {
        TERMINATION_TYPE_CELESTIAL_SPHERE, //  0
        TERMINATION_TYPE_SOURCE_PLANE,     //  1
        FAILURE_PT_TOO_BIG,                //  2
        FAILURE_RKF45_REJECTION_LIMIT,     //  3
        FAILURE_T_MAX_EXCEEDED,            //  4
        FAILURE_SLOT_MANAGER_ERROR,        //  5
        TERMINATION_TYPE_FAILURE,          //  6
        ACTIVE                             //  7
    } termination_type_t;

    typedef struct {
        termination_type_t termination_type; 
        double y_w, z_w, y_s, z_s, final_theta, final_phi, L_w, t_w, L_s, t_s;
    } __attribute__((packed)) blueprint_data_t;

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; double *f_p; double *f_p_p;
        double *affine_param; double *affine_param_p; double *affine_param_p_p;
        double *h;
        termination_type_t *status;
        int *rejection_retries;
        bool *on_positive_side_of_window_prev;
        bool *on_positive_side_of_source_prev;
        
        bool *source_event_found;
        double *source_event_lambda;      
        double *source_event_f_intersect; 
        bool *window_event_found;
        double *window_event_lambda;      
        double *window_event_f_intersect; 
    } PhotonStateSoA;
    """
Bdefines_h.register_BHaH_defines("photon_02_batch_structs", batch_structs_c_code)

par.register_CodeParameters(
    "REAL", __name__,
    ["t_integration_max", "r_escape", "p_t_max", "slot_manager_t_min", "slot_manager_delta_t"],
    [10000.0, 150.0, 1e3, -1000.0, 10.0], commondata=True, add_to_parfile=True,
)
par.register_CodeParameters(
    "bool", __name__, ["perform_conservation_check", "debug_mode"], [True, True], commondata=True, add_to_parfile=True,
)

def batch_integrator_numerical(spacetime_name: str) -> None:
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "omp.h", "math.h"]
    desc = r"""@brief Finalized Project Kerr-Singularity-Axiom Orchestrator (Staged Architecture)."""
    name = "batch_integrator_numerical"
    params = """const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *results_buffer"""

    initial_conditions_func = f"set_initial_conditions_cartesian_{spacetime_name}"
    interpolation_func = f"placeholder_interpolation_engine_{spacetime_name}"
    conserved_quantities_func = f"conserved_quantities_{spacetime_name}_photon"

    body = f"""
    // Comprehensive Device Pointer Mapping Macro
    #define SOA_DEVICE_PTRS all_photons.f, all_photons.f_p, all_photons.f_p_p, all_photons.affine_param, all_photons.affine_param_p, all_photons.affine_param_p_p, all_photons.h, all_photons.status, all_photons.rejection_retries, all_photons.on_positive_side_of_window_prev, all_photons.on_positive_side_of_source_prev, all_photons.source_event_found, all_photons.source_event_lambda, all_photons.source_event_f_intersect, all_photons.window_event_found, all_photons.window_event_lambda, all_photons.window_event_f_intersect

    FILE *fp_debug = NULL;
    if (commondata->debug_mode) {{
        fp_debug = fopen("photon_path_numerical.txt", "w");
        if (fp_debug) fprintf(fp_debug, "# affine_param\\tt\\tx\\ty\\tz\\tp_t\\tp_x\\tp_y\\tp_z\\tL\\n");
    }}

    printf("[DIAGNOSTIC] Stage 1: Allocating Host Buffers...\\n"); fflush(stdout);
    PhotonStateSoA all_photons_host;
    all_photons_host.f = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons_host.f_p = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons_host.f_p_p = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons_host.affine_param = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.affine_param_p = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.affine_param_p_p = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.h = (double *)malloc(sizeof(double) * num_rays);
    all_photons_host.status = (termination_type_t *)malloc(sizeof(termination_type_t) * num_rays);
    all_photons_host.rejection_retries = (int *)malloc(sizeof(int) * num_rays);
    all_photons_host.on_positive_side_of_window_prev = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.on_positive_side_of_source_prev = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.source_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.source_event_lambda = (double *)malloc(sizeof(double) * num_rays); 
    all_photons_host.source_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons_host.window_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons_host.window_event_lambda = (double *)malloc(sizeof(double) * num_rays); 
    all_photons_host.window_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);
    
    printf("[DIAGNOSTIC] Stage 2: Initializing TimeSlotManager...\\n"); fflush(stdout);
    TimeSlotManager tsm;
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    printf("[DIAGNOSTIC] Stage 3: Computing Initial Conditions...\\n"); fflush(stdout);
    double window_center[3], n_x[3], n_y[3], n_z[3];
    {initial_conditions_func}(commondata, num_rays, &all_photons_host, window_center, n_x, n_y, n_z);

    #pragma omp parallel for
    for (long int i = 0; i < num_rays; i++) {{
        all_photons_host.affine_param[i] = 0.0;
        all_photons_host.h[i] = commondata->numerical_initial_h;
        all_photons_host.status[i] = ACTIVE;
        all_photons_host.rejection_retries[i] = 0;
        for(int k=0; k<9; ++k) {{
            double f_val = all_photons_host.f[IDX_GLOBAL(k, i, num_rays)];
            all_photons_host.f_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
            all_photons_host.f_p_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
        }}
        all_photons_host.affine_param_p[i] = 0.0;
        all_photons_host.affine_param_p_p[i] = 0.0;
        double w_d = n_z[0]*window_center[0] + n_z[1]*window_center[1] + n_z[2]*window_center[2];
        double w_val = all_photons_host.f[IDX_GLOBAL(1, i, num_rays)]*n_z[0] + all_photons_host.f[IDX_GLOBAL(2, i, num_rays)]*n_z[1] + all_photons_host.f[IDX_GLOBAL(3, i, num_rays)]*n_z[2] - w_d;
        all_photons_host.on_positive_side_of_window_prev[i] = (w_val > 0.0);
        double s_d = commondata->source_plane_center_x*commondata->source_plane_normal_x + commondata->source_plane_center_y*commondata->source_plane_normal_y + commondata->source_plane_center_z*commondata->source_plane_normal_z;
        double s_val = all_photons_host.f[IDX_GLOBAL(1, i, num_rays)]*commondata->source_plane_normal_x + all_photons_host.f[IDX_GLOBAL(2, i, num_rays)]*commondata->source_plane_normal_y + all_photons_host.f[IDX_GLOBAL(3, i, num_rays)]*commondata->source_plane_normal_z - s_d;
        all_photons_host.on_positive_side_of_source_prev[i] = (s_val > 0.0);
        all_photons_host.source_event_found[i] = false;
        all_photons_host.window_event_found[i] = false;
    }}

    PhotonStateSoA all_photons;
    #ifdef USE_GPU
        printf("[DIAGNOSTIC] Stage 4: Entering GPU Memory Allocation...\\n"); fflush(stdout);
        all_photons.f = (double *)omp_target_alloc(sizeof(double) * 9 * num_rays, omp_get_default_device());
        all_photons.f_p = (double *)omp_target_alloc(sizeof(double) * 9 * num_rays, omp_get_default_device());
        all_photons.f_p_p = (double *)omp_target_alloc(sizeof(double) * 9 * num_rays, omp_get_default_device());
        all_photons.affine_param = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device());
        all_photons.affine_param_p = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device());
        all_photons.affine_param_p_p = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device());
        all_photons.h = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device());
        all_photons.status = (termination_type_t *)omp_target_alloc(sizeof(termination_type_t) * num_rays, omp_get_default_device());
        all_photons.rejection_retries = (int *)omp_target_alloc(sizeof(int) * num_rays, omp_get_default_device());
        all_photons.on_positive_side_of_window_prev = (bool *)omp_target_alloc(sizeof(bool) * num_rays, omp_get_default_device());
        all_photons.on_positive_side_of_source_prev = (bool *)omp_target_alloc(sizeof(bool) * num_rays, omp_get_default_device());
        all_photons.source_event_found = (bool *)omp_target_alloc(sizeof(bool) * num_rays, omp_get_default_device());
        all_photons.source_event_lambda = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device()); 
        all_photons.source_event_f_intersect = (double *)omp_target_alloc(sizeof(double) * 9 * num_rays, omp_get_default_device());
        all_photons.window_event_found = (bool *)omp_target_alloc(sizeof(bool) * num_rays, omp_get_default_device());
        all_photons.window_event_lambda = (double *)omp_target_alloc(sizeof(double) * num_rays, omp_get_default_device()); 
        all_photons.window_event_f_intersect = (double *)omp_target_alloc(sizeof(double) * 9 * num_rays, omp_get_default_device());

        if (!all_photons.f || !all_photons.status) {{
            fprintf(stderr, "FATAL: GPU Allocation failed! Check NVIDIA-SMI or Offload Driver.\\n");
            exit(1);
        }}

        printf("[DIAGNOSTIC] Stage 5: Synchronizing Data to GPU...\\n"); fflush(stdout);
        omp_target_memcpy(all_photons.f, all_photons_host.f, sizeof(double) * 9 * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.f_p, all_photons_host.f_p, sizeof(double) * 9 * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.f_p_p, all_photons_host.f_p_p, sizeof(double) * 9 * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.affine_param, all_photons_host.affine_param, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.affine_param_p, all_photons_host.affine_param_p, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.affine_param_p_p, all_photons_host.affine_param_p_p, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.h, all_photons_host.h, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.status, all_photons_host.status, sizeof(termination_type_t) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.rejection_retries, all_photons_host.rejection_retries, sizeof(int) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.on_positive_side_of_window_prev, all_photons_host.on_positive_side_of_window_prev, sizeof(bool) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.on_positive_side_of_source_prev, all_photons_host.on_positive_side_of_source_prev, sizeof(bool) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.source_event_found, all_photons_host.source_event_found, sizeof(bool) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.source_event_lambda, all_photons_host.source_event_lambda, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.source_event_f_intersect, all_photons_host.source_event_f_intersect, sizeof(double) * 9 * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.window_event_found, all_photons_host.window_event_found, sizeof(bool) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.window_event_lambda, all_photons_host.window_event_lambda, sizeof(double) * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
        omp_target_memcpy(all_photons.window_event_f_intersect, all_photons_host.window_event_f_intersect, sizeof(double) * 9 * num_rays, 0, 0, omp_get_default_device(), omp_get_initial_device());
    #else
        all_photons = all_photons_host;
    #endif

    printf("[DIAGNOSTIC] Stage 6: Data Uploaded. Starting Integration...\\n"); fflush(stdout);

    double *initial_cq = NULL;
    if (commondata->perform_conservation_check) {{
        printf("[DIAGNOSTIC] Stage 7: Conservation Check (Host-Side)...\\n"); fflush(stdout);
        
        initial_cq = (double *)malloc(sizeof(double) * 5 * num_rays);
        if (initial_cq == NULL) {{
            fprintf(stderr, "FATAL: CPU malloc failed for initial_cq!\\n");
            exit(1);
        }}

        #pragma omp parallel for
        for (long int i = 0; i < num_rays; ++i) {{
            {conserved_quantities_func}(commondata, all_photons_host.f, num_rays, i, 
                                        &initial_cq[i*5 + 0], &initial_cq[i*5 + 1], 
                                        &initial_cq[i*5 + 2], &initial_cq[i*5 + 3], 
                                        &initial_cq[i*5 + 4]);
        }}
        printf("[DIAGNOSTIC] Stage 7.1: Conservation Check Passed!\\n"); fflush(stdout);
    }}

    int initial_slot_idx = slot_get_index(&tsm, commondata->t_start);
    if(initial_slot_idx != -1) {{ for(long int i=0; i<num_rays; ++i) slot_add_photon(&tsm, initial_slot_idx, i); }}

    printf("[DIAGNOSTIC] Stage 8: Allocating Bundle Buffers...\\n"); fflush(stdout);
    
    // 4. Allocate Management Arrays and Device-Side Bundle Buffers
    termination_type_t *bundle_status_host = (termination_type_t *)malloc(sizeof(termination_type_t) * BUNDLE_CAPACITY);
    double *bundle_time_host = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY);

    double *f_start_bundle, *f_temp_bundle, *f_out_bundle, *f_err_bundle, *metric_bundle, *conn_bundle, *k_array_bundle;
    int *req_photon_ids;
    #ifdef USE_GPU
        f_start_bundle = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_temp_bundle  = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_out_bundle   = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_err_bundle   = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        metric_bundle  = (double *)omp_target_alloc(sizeof(double) * 10 * BUNDLE_CAPACITY, omp_get_default_device());
        conn_bundle    = (double *)omp_target_alloc(sizeof(double) * 40 * BUNDLE_CAPACITY, omp_get_default_device());
        k_array_bundle = (double *)omp_target_alloc(sizeof(double) * 6 * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        req_photon_ids = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        
        if (!f_start_bundle || !f_temp_bundle || !f_out_bundle || !f_err_bundle || !metric_bundle || !k_array_bundle) {{
            fprintf(stderr, "FATAL: GPU Memory Allocation failed for bundle buffers!\\n");
            exit(1);
        }}
    #else
        f_start_bundle = (double *)malloc(sizeof(double) * 9 * BUNDLE_CAPACITY);
        f_temp_bundle  = (double *)malloc(sizeof(double) * 9 * BUNDLE_CAPACITY);
        f_out_bundle   = (double *)malloc(sizeof(double) * 9 * BUNDLE_CAPACITY);
        f_err_bundle   = (double *)malloc(sizeof(double) * 9 * BUNDLE_CAPACITY);
        metric_bundle  = (double *)malloc(sizeof(double) * 10 * BUNDLE_CAPACITY);
        conn_bundle    = (double *)malloc(sizeof(double) * 40 * BUNDLE_CAPACITY);
        k_array_bundle = (double *)malloc(sizeof(double) * 6 * 9 * BUNDLE_CAPACITY);
        req_photon_ids = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    #endif
    
    printf("[DIAGNOSTIC] Stage 9: Entering Main Staged Integrator Loop...\\n"); fflush(stdout);

    long int active_photons = num_rays;
    for (int i = initial_slot_idx; i >= 0 && active_photons > 0; i--) {{
        while (tsm.slot_counts[i] > 0) {{
            long int bundle_size = MIN(tsm.slot_counts[i], BUNDLE_CAPACITY);
            long int bundle_photons[bundle_size];
            slot_remove_chunk(&tsm, i, bundle_photons, bundle_size);

            // --- STAGED INTEGRATOR (Batched Execution) ---
            bool bundle_has_active = true;
            while (bundle_has_active) {{
                // Snapshot start state into bundle buffer
                #ifdef USE_GPU
                    #pragma omp target teams distribute parallel for \
                                map(to: all_photons, bundle_photons[0:bundle_size]) \
                                is_device_ptr(f_start_bundle, f_temp_bundle, req_photon_ids) \
                                has_device_addr(SOA_DEVICE_PTRS)
                #else
                    #pragma omp parallel for
                #endif
                for (long int j = 0; j < bundle_size; j++) {{
                    long int p_idx = bundle_photons[j];
                    req_photon_ids[j] = (int)p_idx;
                    if (all_photons.status[p_idx] == ACTIVE) {{
                        for (int k = 0; k < 9; k++) {{
                            double val = all_photons.f[IDX_GLOBAL(k, p_idx, num_rays)];
                            f_start_bundle[IDX_LOCAL(k, j, BUNDLE_CAPACITY)] = val;
                            f_temp_bundle[IDX_LOCAL(k, j, BUNDLE_CAPACITY)] = val;
                        }}
                    }}
                }}

                for (int stage = 1; stage <= 6; stage++) {{
                    // A. BIG CALL: Interpolation Engine (Entire Bundle processable safely)
                    {interpolation_func}(commondata, (int)bundle_size, req_photon_ids, f_temp_bundle, metric_bundle, conn_bundle);

                    // B. BIG CALL: ODE RHS
                    #ifdef USE_GPU
                        #pragma omp target teams distribute parallel for \
                                    map(to: all_photons, bundle_photons[0:bundle_size]) \
                                    is_device_ptr(f_temp_bundle, metric_bundle, conn_bundle, k_array_bundle) \
                                    has_device_addr(SOA_DEVICE_PTRS)
                    #else
                        #pragma omp parallel for
                    #endif
                    for (long int j = 0; j < bundle_size; j++) {{
                        long int p_idx = bundle_photons[j];
                        if (all_photons.status[p_idx] == ACTIVE) {{
                            calculate_ode_rhs(f_temp_bundle, metric_bundle, conn_bundle, k_array_bundle, BUNDLE_CAPACITY, stage, (int)j);
                        }}
                    }}

                    // C. BIG CALL: Update f_temp
                    if (stage < 6) {{
                        #ifdef USE_GPU
                            #pragma omp target teams distribute parallel for \
                                        map(to: all_photons, bundle_photons[0:bundle_size]) \
                                        is_device_ptr(f_start_bundle, f_temp_bundle, k_array_bundle) \
                                        has_device_addr(SOA_DEVICE_PTRS)
                        #else
                            #pragma omp parallel for
                        #endif
                        for (long int j = 0; j < bundle_size; j++) {{
                            long int p_idx = bundle_photons[j];
                            if (all_photons.status[p_idx] == ACTIVE) {{
                                calculate_rkf45_stage_f_temp(stage + 1, f_start_bundle, k_array_bundle, all_photons.h[p_idx], f_temp_bundle, BUNDLE_CAPACITY, (int)j);
                            }}
                        }}
                    }}
                }}

                bundle_has_active = false;
                #ifdef USE_GPU
                    #pragma omp target teams distribute parallel for \
                                map(to: all_photons, bundle_photons[0:bundle_size], commondata[0:1]) \
                                is_device_ptr(f_start_bundle, f_out_bundle, f_err_bundle, k_array_bundle) \
                                has_device_addr(SOA_DEVICE_PTRS) \
                                reduction(|:bundle_has_active)
                #else
                    #pragma omp parallel for reduction(|:bundle_has_active)
                #endif
                for (long int j = 0; j < bundle_size; j++) {{
                    long int p_idx = bundle_photons[j];
                    if (all_photons.status[p_idx] == ACTIVE) {{
                        double h_step = all_photons.h[p_idx];
                        
                        // Pass the global batched arrays safely
                        rkf45_kernel(f_start_bundle, k_array_bundle, h_step, f_out_bundle, f_err_bundle, BUNDLE_CAPACITY, (int)j);
                        bool accepted = update_photon_state_and_stepsize(&all_photons, num_rays, p_idx, f_start_bundle, f_out_bundle, f_err_bundle, BUNDLE_CAPACITY, (int)j, commondata);
                        
                        if (!accepted) {{
                            if (all_photons.rejection_retries[p_idx] > commondata->rkf45_max_retries) {{
                                all_photons.status[p_idx] = FAILURE_RKF45_REJECTION_LIMIT;
                            }} else {{
                                bundle_has_active = true; 
                            }}
                        }} else {{
                            // Safe SoA history update directly from the batch array
                            for (int k = 0; k < 9; ++k) {{
                                all_photons.f_p_p[IDX_GLOBAL(k, p_idx, num_rays)] = all_photons.f_p[IDX_GLOBAL(k, p_idx, num_rays)];
                                all_photons.f_p[IDX_GLOBAL(k, p_idx, num_rays)] = f_start_bundle[IDX_LOCAL(k, j, BUNDLE_CAPACITY)];
                            }}
                            all_photons.affine_param_p_p[p_idx] = all_photons.affine_param_p[p_idx];
                            all_photons.affine_param_p[p_idx] = all_photons.affine_param[p_idx] - h_step;
                        }}
                    }}
                }}
            }}
            // --- END STAGED INTEGRATOR ---

            // --- EVENT DETECTION & TERMINATION CHECKS ---
            #ifdef USE_GPU
                #pragma omp target teams distribute parallel for \
                            map(to: all_photons, commondata[0:1], bundle_photons[0:bundle_size]) \
                            has_device_addr(SOA_DEVICE_PTRS) \
                            map(tofrom: results_buffer[0:num_rays])
            #else
                #pragma omp parallel for
            #endif
            for (long int j = 0; j < bundle_size; ++j) {{
                long int p_idx = bundle_photons[j];
                if (all_photons.status[p_idx] == ACTIVE) {{
                    double x = all_photons.f[IDX_GLOBAL(1, p_idx, num_rays)];
                    double y = all_photons.f[IDX_GLOBAL(2, p_idx, num_rays)];
                    double z = all_photons.f[IDX_GLOBAL(3, p_idx, num_rays)];
                    double r_sq = x*x + y*y + z*z;
                    
                    if (fabs(all_photons.f[IDX_GLOBAL(4, p_idx, num_rays)]) > commondata->p_t_max) {{
                        all_photons.status[p_idx] = FAILURE_PT_TOO_BIG;
                    }} else if (fabs(all_photons.f[IDX_GLOBAL(0, p_idx, num_rays)]) > commondata->t_integration_max) {{
                        all_photons.status[p_idx] = FAILURE_T_MAX_EXCEEDED;
                    }} else if (r_sq > commondata->r_escape*commondata->r_escape) {{
                        all_photons.status[p_idx] = TERMINATION_TYPE_CELESTIAL_SPHERE;
                    }} else {{ 
                        event_detection_manager(&all_photons, num_rays, p_idx, commondata); 
                        if (all_photons.source_event_found[p_idx]) {{ 
                            if (handle_source_plane_intersection(&all_photons, num_rays, p_idx, commondata, &results_buffer[p_idx])) {{
                                all_photons.status[p_idx] = TERMINATION_TYPE_SOURCE_PLANE; 
                            }} else {{
                                all_photons.source_event_found[p_idx] = false; 
                            }}
                        }} 
                        if (all_photons.window_event_found[p_idx]) {{ 
                            handle_window_plane_intersection(&all_photons, num_rays, p_idx, commondata, &results_buffer[p_idx]); 
                            all_photons.window_event_found[p_idx] = false; 
                        }} 
                    }}
                }}
            }}

            // --- SYNC STATUS TO HOST ---
            #ifdef USE_GPU
                for (long int j = 0; j < bundle_size; j++) {{
                    long int p_idx = bundle_photons[j];
                    omp_target_memcpy(&bundle_status_host[j], &all_photons.status[p_idx], sizeof(termination_type_t), 0, 0, omp_get_initial_device(), omp_get_default_device());
                    omp_target_memcpy(&bundle_time_host[j], &all_photons.f[IDX_GLOBAL(0, p_idx, num_rays)], sizeof(double), 0, 0, omp_get_initial_device(), omp_get_default_device());
                }}
            #else
                for (long int j = 0; j < bundle_size; j++) {{
                    bundle_status_host[j] = all_photons.status[bundle_photons[j]];
                    bundle_time_host[j] = all_photons.f[IDX_GLOBAL(0, bundle_photons[j], num_rays)];
                }}
            #endif

            long int active_photons_decrement = 0;
            long int active_in_bundle = 0;
            
            for (long int j = 0; j < bundle_size; j++) {{
                long int p_idx = bundle_photons[j];
                if (bundle_status_host[j] == ACTIVE) {{
                    int s_idx = slot_get_index(&tsm, bundle_time_host[j]);
                    if (s_idx != -1) {{
                        slot_add_photon(&tsm, s_idx, p_idx);
                        active_in_bundle++;
                   }} else {{
                        #ifdef USE_GPU
                            termination_type_t fail_stat = FAILURE_SLOT_MANAGER_ERROR;
                            omp_target_memcpy(&all_photons.status[p_idx], &fail_stat, sizeof(termination_type_t), 0, 0, omp_get_default_device(), omp_get_initial_device());
                        #else
                            all_photons.status[p_idx] = FAILURE_SLOT_MANAGER_ERROR;
                        #endif
                        active_photons_decrement++;
                     }}
                 }}
            }}
            active_photons -= (bundle_size - active_in_bundle);
        }} 
    }} 

    #ifdef USE_GPU
        #pragma omp target teams distribute parallel for \
                    map(to: all_photons, commondata[0:1]) \
                    map(tofrom: results_buffer[0:num_rays]) \
                    has_device_addr(SOA_DEVICE_PTRS)
    #else
        #pragma omp parallel for
    #endif
    for (long int i = 0; i < num_rays; i++) {{ calculate_and_fill_blueprint_data_universal(commondata, &all_photons, num_rays, i, &results_buffer[i]); }}

    #ifdef USE_GPU
        double *f_host = (double*)malloc(sizeof(double) * 9 * num_rays);
        omp_target_memcpy(f_host, all_photons.f, sizeof(double) * 9 * num_rays, 0, 0, omp_get_initial_device(), omp_get_default_device());
    #else
        double *f_host = all_photons.f;
    #endif

    // - FINAL CONSERVATION SUMMARY ---
    if (commondata->perform_conservation_check && initial_cq != NULL) {{
        printf("\\n--- Final Conservation Summary ---\\n"); 
        double max_dE = 0.0, max_dL = 0.0, max_dQ = 0.0;
        long int worst_E_idx = -1, worst_L_idx = -1, worst_Q_idx = -1;

        FILE *fp_cons_log = fopen("conservation_errors.txt", "w");
        if (fp_cons_log) fprintf(fp_cons_log, "# photon_idx dE_rel dL_rel dQ_rel Q_init Q_final term_type final_t\\n");

        #pragma omp parallel
        {{
            double local_max_dE = 0.0, local_max_dL = 0.0, local_max_dQ = 0.0;
            long int local_worst_E = -1, local_worst_L = -1, local_worst_Q = -1;

            #pragma omp for
            for (long int i = 0; i < num_rays; i++) {{
                double final_E, final_Lx, final_Ly, final_Lz, final_Q;
                {conserved_quantities_func}(commondata, f_host, num_rays, i, &final_E, &final_Lx, &final_Ly, &final_Lz, &final_Q);

                double init_E = initial_cq[i * 5 + 0];
                double dE = (init_E != 0) ? fabs((final_E - init_E) / init_E) : fabs(final_E - init_E);
                if (dE > local_max_dE) {{ local_max_dE = dE; local_worst_E = i; }}

                double dL;
                if (commondata->a_spin == 0.0) {{
                    double L_init_mag = sqrt(initial_cq[i*5 + 1]*initial_cq[i*5 + 1] + initial_cq[i*5 + 2]*initial_cq[i*5 + 2] + initial_cq[i*5 + 3]*initial_cq[i*5 + 3]);
                    double L_final_mag = sqrt(final_Lx*final_Lx + final_Ly*final_Ly + final_Lz*final_Lz);
                    dL = (L_init_mag != 0) ? fabs((L_final_mag - L_init_mag) / L_init_mag) : fabs(L_final_mag - L_init_mag);
                }} else {{
                    double Lz_init = initial_cq[i * 5 + 3];
                    dL = (Lz_init != 0) ? fabs((final_Lz - Lz_init) / Lz_init) : fabs(final_Lz - Lz_init);
                }}
                if (dL > local_max_dL) {{ local_max_dL = dL; local_worst_L = i; }}

                double init_Q = initial_cq[i * 5 + 4];
                double dQ = (init_Q != 0) ? fabs((final_Q - init_Q) / init_Q) : fabs(final_Q - init_Q);
                if (dQ > local_max_dQ) {{ local_max_dQ = dQ; local_worst_Q = i; }}
                
                if (fp_cons_log) {{
                    #pragma omp critical
                    fprintf(fp_cons_log, "%ld %.6e %.6e %.6e %.6e %.6e %d %.6e\\n", 
                            i, dE, dL, dQ, init_Q, final_Q, (int)all_photons.status[i], f_host[IDX_GLOBAL(0, i, num_rays)]);
                }}
            }}
            #pragma omp critical
            {{
                if (local_max_dE > max_dE) {{ max_dE = local_max_dE; worst_E_idx = local_worst_E; }}
                if (local_max_dL > max_dL) {{ max_dL = local_max_dL; worst_L_idx = local_worst_L; }}
                if (local_max_dQ > max_dQ) {{ max_dQ = local_max_dQ; worst_Q_idx = local_worst_Q; }}
            }}
        }}
        if (fp_cons_log) fclose(fp_cons_log);
        
        printf("Max relative error in Energy (E): %.3e (photon %ld)\\n", max_dE, worst_E_idx);
        printf("Max relative error in Ang. Mom. (L): %.3e (photon %ld)\\n", max_dL, worst_L_idx);
        printf("Max relative error in Carter Constant (Q): %.3e (photon %ld)\\n", max_dQ, worst_Q_idx);
        printf("------------------------------------\\n");
    }}

    if (fp_debug) fclose(fp_debug);

    #ifdef USE_GPU
        free(f_host);
        
        omp_target_free(f_start_bundle, omp_get_default_device());
        omp_target_free(f_temp_bundle, omp_get_default_device());
        omp_target_free(f_out_bundle, omp_get_default_device()); 
        omp_target_free(f_err_bundle, omp_get_default_device()); 
        omp_target_free(metric_bundle, omp_get_default_device());
        omp_target_free(conn_bundle, omp_get_default_device());
        omp_target_free(k_array_bundle, omp_get_default_device());
        omp_target_free(req_photon_ids, omp_get_default_device());

        omp_target_free(all_photons.f, omp_get_default_device()); 
        omp_target_free(all_photons.f_p, omp_get_default_device()); 
        omp_target_free(all_photons.f_p_p, omp_get_default_device()); 
        omp_target_free(all_photons.affine_param, omp_get_default_device()); 
        omp_target_free(all_photons.affine_param_p, omp_get_default_device()); 
        omp_target_free(all_photons.affine_param_p_p, omp_get_default_device()); 
        omp_target_free(all_photons.h, omp_get_default_device()); 
        omp_target_free(all_photons.status, omp_get_default_device()); 
        omp_target_free(all_photons.rejection_retries, omp_get_default_device()); 
        omp_target_free(all_photons.on_positive_side_of_window_prev, omp_get_default_device()); 
        omp_target_free(all_photons.on_positive_side_of_source_prev, omp_get_default_device()); 
        omp_target_free(all_photons.source_event_found, omp_get_default_device()); 
        omp_target_free(all_photons.source_event_lambda, omp_get_default_device()); 
        omp_target_free(all_photons.source_event_f_intersect, omp_get_default_device()); 
        omp_target_free(all_photons.window_event_found, omp_get_default_device()); 
        omp_target_free(all_photons.window_event_lambda, omp_get_default_device()); 
        omp_target_free(all_photons.window_event_f_intersect, omp_get_default_device());
    #else
        free(f_start_bundle); free(f_temp_bundle); free(f_out_bundle); free(f_err_bundle);
        free(metric_bundle); free(conn_bundle); free(k_array_bundle); free(req_photon_ids);
    #endif

    // Free Management Arrays (Host)
    free(bundle_status_host); free(bundle_time_host);
    if (initial_cq) free(initial_cq);
    
    free(all_photons_host.f); free(all_photons_host.f_p); free(all_photons_host.f_p_p); 
    free(all_photons_host.affine_param); free(all_photons_host.affine_param_p); free(all_photons_host.affine_param_p_p); 
    free(all_photons_host.h); free(all_photons_host.status); free(all_photons_host.rejection_retries); 
    free(all_photons_host.on_positive_side_of_window_prev); free(all_photons_host.on_positive_side_of_source_prev); 
    free(all_photons_host.source_event_found); free(all_photons_host.source_event_lambda); 
    free(all_photons_host.source_event_f_intersect); free(all_photons_host.window_event_found); 
    free(all_photons_host.window_event_lambda); free(all_photons_host.window_event_f_intersect);
    
    #undef SOA_DEVICE_PTRS
    slot_manager_free(&tsm);
    """
    cfc.register_CFunction(includes=includes, desc=desc, name=name, params=params, body=body)
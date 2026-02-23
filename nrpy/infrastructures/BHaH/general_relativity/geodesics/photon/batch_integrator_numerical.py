# batch_integrator_numerical.py
"""
Generates the main C orchestrator for the numerical integration pipeline.
Corrected to prevent Window Plane termination and improve debug logging.
Updated to restore conservation checks using the clean SoA architecture.
Includes final blueprint data extraction fix.

Author: Dalton J. Moone
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
        double *f;                 
        double *f_p;               
        double *f_p_p;             
        double *affine_param;
        double *affine_param_p;
        double *affine_param_p_p;
        double *h;
        termination_type_t *status;
        int *rejection_retries;
        bool *on_positive_side_of_window_prev;
        bool *on_positive_side_of_source_prev;
        
        // Flattened Event Data
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
    desc = r"""@brief Main orchestrator with strict SIMT-compatible SoA architecture and simple conservation checks."""
    name = "batch_integrator_numerical"
    params = """const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *results_buffer"""

    initial_conditions_func = f"set_initial_conditions_cartesian_{spacetime_name}"
    interpolation_func = f"placeholder_interpolation_engine_{spacetime_name}"
    conserved_quantities_func = f"conserved_quantities_{spacetime_name}_photon"

    body = f"""
    // --- DEBUG MODE SETUP ---
    FILE *fp_debug = NULL;
    if (commondata->debug_mode) {{
        fp_debug = fopen("photon_path_numerical.txt", "w");
        if (fp_debug) {{
            fprintf(fp_debug, "# affine_param\\tt\\tx\\ty\\tz\\tp_t\\tp_x\\tp_y\\tp_z\\tL\\n");
        }}
    }}

    printf("Initializing SoA architecture for %ld photon states...\\n", num_rays);
    
    PhotonStateSoA all_photons;
    all_photons.f = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons.f_p = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons.f_p_p = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons.affine_param = (double *)malloc(sizeof(double) * num_rays);
    all_photons.affine_param_p = (double *)malloc(sizeof(double) * num_rays);
    all_photons.affine_param_p_p = (double *)malloc(sizeof(double) * num_rays);
    all_photons.h = (double *)malloc(sizeof(double) * num_rays);
    all_photons.status = (termination_type_t *)malloc(sizeof(termination_type_t) * num_rays);
    all_photons.rejection_retries = (int *)malloc(sizeof(int) * num_rays);
    all_photons.on_positive_side_of_window_prev = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons.on_positive_side_of_source_prev = (bool *)malloc(sizeof(bool) * num_rays);
    
    all_photons.source_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons.source_event_lambda = (double *)malloc(sizeof(double) * num_rays); 
    all_photons.source_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);
    all_photons.window_event_found = (bool *)malloc(sizeof(bool) * num_rays);
    all_photons.window_event_lambda = (double *)malloc(sizeof(double) * num_rays); 
    all_photons.window_event_f_intersect = (double *)malloc(sizeof(double) * 9 * num_rays);
    
    long int active_photons = num_rays;

    TimeSlotManager tsm;
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    double window_center[3], n_x[3], n_y[3], n_z[3];
    {initial_conditions_func}(commondata, num_rays, &all_photons, window_center, n_x, n_y, n_z);

    #pragma omp parallel for
    for (long int i = 0; i < num_rays; i++) {{
        all_photons.affine_param[i] = 0.0;
        all_photons.h[i] = commondata->numerical_initial_h;
        all_photons.status[i] = ACTIVE;
        all_photons.rejection_retries[i] = 0;

        for(int k=0; k<9; ++k) {{
            double f_val = all_photons.f[IDX_GLOBAL(k, i, num_rays)];
            all_photons.f_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
            all_photons.f_p_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
        }}
        all_photons.affine_param_p[i] = all_photons.affine_param[i];
        all_photons.affine_param_p_p[i] = all_photons.affine_param[i];

        // Inline window plane check
        double w_d = n_z[0]*window_center[0] + n_z[1]*window_center[1] + n_z[2]*window_center[2];
        double w_val = all_photons.f[IDX_GLOBAL(1, i, num_rays)] * n_z[0] +
                       all_photons.f[IDX_GLOBAL(2, i, num_rays)] * n_z[1] +
                       all_photons.f[IDX_GLOBAL(3, i, num_rays)] * n_z[2] - w_d;
        all_photons.on_positive_side_of_window_prev[i] = (w_val > 0.0);

        // Inline source plane check
        double s_d = commondata->source_plane_center_x*commondata->source_plane_normal_x + 
                     commondata->source_plane_center_y*commondata->source_plane_normal_y + 
                     commondata->source_plane_center_z*commondata->source_plane_normal_z;
        double s_val = all_photons.f[IDX_GLOBAL(1, i, num_rays)] * commondata->source_plane_normal_x +
                       all_photons.f[IDX_GLOBAL(2, i, num_rays)] * commondata->source_plane_normal_y +
                       all_photons.f[IDX_GLOBAL(3, i, num_rays)] * commondata->source_plane_normal_z - s_d;
        all_photons.on_positive_side_of_source_prev[i] = (s_val > 0.0);

        all_photons.source_event_found[i] = false;
        all_photons.window_event_found[i] = false;
    }}

    // --- LOG INITIAL STATE FOR PHOTON 0 ---
    if (commondata->debug_mode && fp_debug) {{
        fprintf(fp_debug, "%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\n",
                0.0, 
                all_photons.f[IDX_GLOBAL(0, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(1, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(2, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(3, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(4, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(5, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(6, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(7, 0, num_rays)],
                all_photons.f[IDX_GLOBAL(8, 0, num_rays)]);
        fflush(fp_debug);
    }}

    // --- INITIAL CONSERVATION CHECK ---
    double *initial_cq = NULL;
    if (commondata->perform_conservation_check) {{
        printf("Performing initial conservation checks for all %ld rays...\\n", num_rays);
        initial_cq = (double *)malloc(sizeof(double) * 5 * num_rays);
        
        #pragma omp parallel for
        for (long int i = 0; i < num_rays; ++i) {{
            {conserved_quantities_func}(commondata, all_photons.f, num_rays, i,
                                        &initial_cq[i * 5 + 0], &initial_cq[i * 5 + 1],
                                        &initial_cq[i * 5 + 2], &initial_cq[i * 5 + 3],
                                        &initial_cq[i * 5 + 4]);
        }}
    }}

    int initial_slot_idx = slot_get_index(&tsm, commondata->t_start);
    if(initial_slot_idx != -1) {{
        for(long int i=0; i<num_rays; ++i) {{ slot_add_photon(&tsm, initial_slot_idx, i); }}
    }}

    int *req_photon_ids = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    double *req_pos = (double *)malloc(sizeof(double) * 4 * BUNDLE_CAPACITY);
    double *metric_g4DD = (double *)malloc(sizeof(double) * 10 * BUNDLE_CAPACITY);
    double *conn_GammaUDD = (double *)malloc(sizeof(double) * 40 * BUNDLE_CAPACITY);
    
    double *k_array = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY * 6 * 9);
    double *f_start = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY * 9);
    double *f_temp  = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY * 9);
    double *f_out   = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY * 9);
    double *f_err   = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY * 9);

    int *retry_keep_mask = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    int *retry_scan_idx  = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    int *event_keep_mask = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    int *event_scan_idx  = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
    long int *compacted_active = (long int *)malloc(sizeof(long int) * BUNDLE_CAPACITY);

    // === MAIN EPOCH LOOP ===
    for (int i = initial_slot_idx; i >= 0 && active_photons > 0; i--) {{
        while (tsm.slot_counts[i] > 0) {{
            long int bundle_size = MIN(tsm.slot_counts[i], BUNDLE_CAPACITY);
            long int bundle_photons[bundle_size];
            slot_remove_chunk(&tsm, i, bundle_photons, bundle_size);

            long int needs_retry_indices[bundle_size];
            long int next_retry_indices[bundle_size]; 
            long int needs_retry_count = bundle_size;
            for(long int j=0; j<bundle_size; ++j) needs_retry_indices[j] = bundle_photons[j];

            while (needs_retry_count > 0) {{
                #pragma omp parallel for
                for(long int j=0; j<needs_retry_count; ++j) {{
                    long int photon_idx = needs_retry_indices[j];
                    for(int k=0; k<9; ++k) f_start[IDX_LOCAL(k, j, BUNDLE_CAPACITY)] = all_photons.f[IDX_GLOBAL(k, photon_idx, num_rays)];
                }}

                for (int stage = 1; stage <= 6; ++stage) {{
                    #pragma omp parallel for
                    for (long int j = 0; j < needs_retry_count; ++j) {{
                        calculate_rkf45_stage_f_temp(stage, f_start, k_array, all_photons.h[needs_retry_indices[j]], f_temp, BUNDLE_CAPACITY, j);
                        req_photon_ids[j] = needs_retry_indices[j];
                        for(int m=0; m<4; m++) req_pos[IDX_LOCAL(m, j, BUNDLE_CAPACITY)] = f_temp[IDX_LOCAL(m, j, BUNDLE_CAPACITY)];
                    }}
                    {interpolation_func}(commondata, needs_retry_count, req_photon_ids, req_pos, metric_g4DD, conn_GammaUDD);
                    #pragma omp parallel for
                    for (long int j = 0; j < needs_retry_count; ++j) calculate_ode_rhs(f_temp, metric_g4DD, conn_GammaUDD, k_array, BUNDLE_CAPACITY, stage, j);
                }}

                long int current_retry_count = needs_retry_count;
                needs_retry_count = 0;

                #pragma omp parallel for
                for (long int j = 0; j < current_retry_count; ++j) {{
                    long int photon_idx = needs_retry_indices[j];
                    rkf45_kernel(f_start, k_array, all_photons.h[photon_idx], f_out, f_err, BUNDLE_CAPACITY, j);

                    double h_taken = all_photons.h[photon_idx]; 
                    bool step_accepted = update_photon_state_and_stepsize(&all_photons, num_rays, photon_idx, f_start, f_out, f_err, BUNDLE_CAPACITY, j, commondata);

                    if (step_accepted) {{
                        for(int k=0; k<9; ++k) {{
                            all_photons.f_p_p[IDX_GLOBAL(k, photon_idx, num_rays)] = all_photons.f_p[IDX_GLOBAL(k, photon_idx, num_rays)];
                            all_photons.f_p[IDX_GLOBAL(k, photon_idx, num_rays)] = f_start[IDX_LOCAL(k, j, BUNDLE_CAPACITY)];
                        }}
                        all_photons.affine_param_p_p[photon_idx] = all_photons.affine_param_p[photon_idx];
                        all_photons.affine_param_p[photon_idx] = all_photons.affine_param[photon_idx] - h_taken;
                        retry_keep_mask[j] = 0;

                        // --- LOG PATH FOR PHOTON 0 ---
                        if (commondata->debug_mode && photon_idx == 0 && fp_debug) {{
                            #pragma omp critical
                            {{
                                fprintf(fp_debug, "%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\t%.15e\\n",
                                        all_photons.affine_param[0],
                                        all_photons.f[IDX_GLOBAL(0, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(1, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(2, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(3, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(4, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(5, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(6, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(7, 0, num_rays)],
                                        all_photons.f[IDX_GLOBAL(8, 0, num_rays)]);
                                fflush(fp_debug);
                            }}
                        }}
                    }} else {{
                        retry_keep_mask[j] = (all_photons.rejection_retries[photon_idx] > commondata->rkf45_max_retries) ? (all_photons.status[photon_idx] = FAILURE_RKF45_REJECTION_LIMIT, 0) : 1;
                    }}
                }}
                
                if (current_retry_count > 0) {{
                    retry_scan_idx[0] = 0;
                    for (long int j = 1; j < current_retry_count; ++j) retry_scan_idx[j] = retry_scan_idx[j-1] + retry_keep_mask[j-1];
                    needs_retry_count = retry_scan_idx[current_retry_count-1] + retry_keep_mask[current_retry_count-1];
                    for (long int j = 0; j < current_retry_count; ++j) if (retry_keep_mask[j]) next_retry_indices[ retry_scan_idx[j] ] = needs_retry_indices[j];
                    for (long int j = 0; j < needs_retry_count; ++j) needs_retry_indices[j] = next_retry_indices[j];
                }}
            }} 

            #pragma omp parallel for
            for (long int j = 0; j < bundle_size; ++j) {{
                long int photon_idx = bundle_photons[j];
                if (all_photons.status[photon_idx] == ACTIVE) {{
                    const double p_t = all_photons.f[IDX_GLOBAL(4, photon_idx, num_rays)];
                    const double x = all_photons.f[IDX_GLOBAL(1, photon_idx, num_rays)];
                    const double y = all_photons.f[IDX_GLOBAL(2, photon_idx, num_rays)];
                    const double z = all_photons.f[IDX_GLOBAL(3, photon_idx, num_rays)];
                    const double r_sq = (x*x) + (y*y) + (z*z);

                    if (fabs(p_t) > commondata->p_t_max) {{
                        all_photons.status[photon_idx] = FAILURE_PT_TOO_BIG;
                    }} else if (fabs(all_photons.f[IDX_GLOBAL(0, photon_idx, num_rays)]) > commondata->t_integration_max) {{
                        all_photons.status[photon_idx] = FAILURE_T_MAX_EXCEEDED;
                    }} else if (r_sq > (commondata->r_escape * commondata->r_escape)) {{
                        all_photons.status[photon_idx] = TERMINATION_TYPE_CELESTIAL_SPHERE;
                        if (r_sq > 1e-9) {{
                            results_buffer[photon_idx].final_theta = acos(z / sqrt(r_sq));
                            results_buffer[photon_idx].final_phi = atan2(y, x);
                        }}
                    }} else {{
                        event_detection_manager(&all_photons, num_rays, photon_idx, commondata);
                        
                        if (all_photons.source_event_found[photon_idx]) {{
                            if (handle_source_plane_intersection(&all_photons, num_rays, photon_idx, commondata, &results_buffer[photon_idx]))
                                all_photons.status[photon_idx] = TERMINATION_TYPE_SOURCE_PLANE;
                            else all_photons.source_event_found[photon_idx] = false; 
                        }}
                        
                        if (all_photons.window_event_found[photon_idx]) {{
                            handle_window_plane_intersection(&all_photons, num_rays, photon_idx, commondata, &results_buffer[photon_idx]);
                            all_photons.window_event_found[photon_idx] = false; 
                        }}
                    }}
                }}
                event_keep_mask[j] = (all_photons.status[photon_idx] == ACTIVE) ? 1 : 0;
            }}
            
            event_scan_idx[0] = 0;
            for (long int j = 1; j < bundle_size; ++j) event_scan_idx[j] = event_scan_idx[j-1] + event_keep_mask[j-1];
            long int active_in_bundle = event_scan_idx[bundle_size-1] + event_keep_mask[bundle_size-1];
            active_photons -= (bundle_size - active_in_bundle);
            for (long int j = 0; j < bundle_size; ++j) if (event_keep_mask[j]) compacted_active[ event_scan_idx[j] ] = bundle_photons[j];
            for (long int j = 0; j < active_in_bundle; ++j) {{
                long int p_idx = compacted_active[j];
                int s_idx = slot_get_index(&tsm, all_photons.f[IDX_GLOBAL(0, p_idx, num_rays)]);
                if (s_idx != -1) slot_add_photon(&tsm, s_idx, p_idx);
                else {{ all_photons.status[p_idx] = FAILURE_SLOT_MANAGER_ERROR; active_photons--; }}
            }}
        }} 
    }} 

    // --- POPULATE BLUEPRINT BUFFER ---
    #pragma omp parallel for
    for (long int i = 0; i < num_rays; i++) {{
        calculate_and_fill_blueprint_data_universal(
            commondata, &all_photons, num_rays, i, &results_buffer[i]
        );
    }}

    // --- FINAL CONSERVATION CHECK ---
    if (commondata->perform_conservation_check && initial_cq != NULL) {{
        printf("\\n--- Conservation Check Summary ---\\n");
        double max_dE = 0.0, max_dL = 0.0, max_dQ = 0.0;
        long int worst_E_idx = -1, worst_L_idx = -1, worst_Q_idx = -1;
        
        FILE *fp_cons_log = fopen("conservation_errors.txt", "w");
        
        struct ConsErrData {{
            double dE, dL, dQ, init_Q, final_Q, final_t;
            int status;
        }};
        struct ConsErrData *err_data = NULL;
        if (fp_cons_log) {{
            fprintf(fp_cons_log, "# photon_idx dE_rel dL_rel dQ_rel Q_init Q_final term_type final_t\\n");
            err_data = (struct ConsErrData *)malloc(sizeof(struct ConsErrData) * num_rays);
        }}

        #pragma omp parallel
        {{
            double local_max_dE = 0.0, local_max_dL = 0.0, local_max_dQ = 0.0;
            long int local_worst_E = -1, local_worst_L = -1, local_worst_Q = -1;
            
            #pragma omp for
            for (long int i = 0; i < num_rays; i++) {{
                double final_E, final_Lx, final_Ly, final_Lz, final_Q;
                
                {conserved_quantities_func}(commondata, all_photons.f, num_rays, i,
                                            &final_E, &final_Lx, &final_Ly, &final_Lz, &final_Q);

                double init_E = initial_cq[i * 5 + 0];
                double dE = (init_E != 0) ? fabs((final_E - init_E) / init_E) : fabs(final_E - init_E);
                if (dE > local_max_dE) {{ local_max_dE = dE; local_worst_E = i; }}

                double dL;
                if (commondata->a_spin == 0.0) {{
                    double L_init_mag = sqrt(initial_cq[i*5 + 1]*initial_cq[i*5 + 1] + 
                                             initial_cq[i*5 + 2]*initial_cq[i*5 + 2] + 
                                             initial_cq[i*5 + 3]*initial_cq[i*5 + 3]);
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

                if (err_data) {{
                    err_data[i].dE = dE;
                    err_data[i].dL = dL;
                    err_data[i].dQ = dQ;
                    err_data[i].init_Q = init_Q;
                    err_data[i].final_Q = final_Q;
                    err_data[i].status = all_photons.status[i];
                    err_data[i].final_t = all_photons.f[IDX_GLOBAL(0, i, num_rays)];
                }}
            }}
            
            #pragma omp critical
            {{
                if (local_max_dE > max_dE) {{ max_dE = local_max_dE; worst_E_idx = local_worst_E; }}
                if (local_max_dL > max_dL) {{ max_dL = local_max_dL; worst_L_idx = local_worst_L; }}
                if (local_max_dQ > max_dQ) {{ max_dQ = local_max_dQ; worst_Q_idx = local_worst_Q; }}
            }}
        }}

        if (fp_cons_log) {{
            if (err_data) {{
                for (long int i = 0; i < num_rays; i++) {{
                    fprintf(fp_cons_log, "%ld %.6e %.6e %.6e %.6e %.6e %d %.6e\\n", 
                            i, err_data[i].dE, err_data[i].dL, err_data[i].dQ, 
                            err_data[i].init_Q, err_data[i].final_Q, err_data[i].status, err_data[i].final_t);
                }}
                free(err_data);
            }}
            fclose(fp_cons_log);
            printf("Full conservation error data saved to 'conservation_errors.txt'.\\n");
        }}

        printf("Max relative error in Energy (E): %.3e (photon %ld)\\n", max_dE, worst_E_idx);
        printf("Max relative error in Ang. Mom. (L): %.3e (photon %ld)\\n", max_dL, worst_L_idx);
        printf("Max relative error in Carter Constant (Q): %.3e (photon %ld)\\n", max_dQ, worst_Q_idx);
        printf("------------------------------------\\n");
        free(initial_cq);
    }}

    if (fp_debug) fclose(fp_debug);
    free(req_photon_ids); free(req_pos); free(metric_g4DD); free(conn_GammaUDD);
    free(k_array); free(f_start); free(f_temp); free(f_out); free(f_err);
    free(retry_keep_mask); free(retry_scan_idx); free(event_keep_mask); free(event_scan_idx); free(compacted_active);
    free(all_photons.f); free(all_photons.f_p); free(all_photons.f_p_p);
    free(all_photons.affine_param); free(all_photons.affine_param_p); free(all_photons.affine_param_p_p);
    free(all_photons.h); free(all_photons.status); free(all_photons.rejection_retries);
    free(all_photons.on_positive_side_of_window_prev); free(all_photons.on_positive_side_of_source_prev);
    free(all_photons.source_event_found); free(all_photons.source_event_lambda); free(all_photons.source_event_f_intersect);
    free(all_photons.window_event_found); free(all_photons.window_event_lambda); free(all_photons.window_event_f_intersect);
    slot_manager_free(&tsm);
    """
    cfc.register_CFunction(includes=includes, desc=desc, name=name, params=params, body=body)
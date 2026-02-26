"""
Orchestration module for the Project Singularity-Axiom numerical integration pipeline.

This module generates the high-level C orchestrator responsible for managing the 
life cycle of photon trajectories in curved spacetimes. It implements a 
dual-architecture (CPU/GPU) compatible integration loop that strictly adheres 
to the following architectural constraints:

1. Memory Architecture: Utilizes a flattened Structure of Arrays (SoA) via 
   the `PhotonStateSoA` and `blueprint_data_t` structures. Access is governed 
   by the `IDX_GLOBAL` and `IDX_LOCAL` macros to ensure cache-efficient 
   strided access and prevent pointer arithmetic errors.
2. The 9-Component Mapping: The state vector `f` is strictly immutable: 
   f[0]:t; f[1,2,3]:x,y,z; f[4]:p_t; f[5,6,7]:p_x,p_y,p_z; f[8]:lambda.
3. Temporal Binning: Employs a lock-free `TimeSlotManager` to bin active 
   rays by physical coordinate time (t), enabling synchronized batch processing.
4. Stream Compaction: Implements a staged integration pattern where photons 
   are processed in dense bundles (BUNDLE_CAPACITY = 32768) to maintain high GPU 
   occupancy, utilizing a local reverse-map to circumvent OpenMP race conditions.
"""

### STEP 1: MODULE IMPORTS & SETUP ###
import nrpy.c_function as cfc
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par

### STEP 2: C-STRUCT & PARAMETER REGISTRATION ###
batch_structs_c_code = r"""
    #define BUNDLE_CAPACITY 32768 // Maximum number of photons processed per batch to fit within L1/L2 cache.

    // Strictly enforce 1D mapping to prevent pointer offset arithmetic bugs
    #define IDX_GLOBAL(component, ray_id, num_rays) ((component) * (num_rays) + (ray_id))
    #define IDX_LOCAL(component, batch_id, batch_size) ((component) * (batch_size) + (batch_id))

    // Defines the physical planes where a photon trajectory might terminate.
    typedef enum {
        WINDOW_EVENT, // Intersection with the observer's camera window.
        SOURCE_EVENT  // Intersection with the emission source plane.
    } event_type_t;

    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        TERMINATION_TYPE_CELESTIAL_SPHERE, // 0: Photon escaped to infinity (exceeded r_escape).
        TERMINATION_TYPE_SOURCE_PLANE,     // 1: Photon successfully hit the source emission plane.
        FAILURE_PT_TOO_BIG,                // 2: Integration failed due to unbounded momentum $p_t$.
        FAILURE_RKF45_REJECTION_LIMIT,     // 3: Adaptive step-size rejected too many consecutive times.
        FAILURE_T_MAX_EXCEEDED,            // 4: Integration exceeded maximum allowable physical time.
        FAILURE_SLOT_MANAGER_ERROR,        // 5: TimeSlotManager failed to allocate or retrieve the photon.
        TERMINATION_TYPE_FAILURE,          // 6: Generic unclassified numerical failure.
        ACTIVE                             // 7: Photon is currently undergoing integration.
    } termination_type_t;

    // Stores the final physical properties of a photon upon integration termination.
    typedef struct {
        termination_type_t termination_type; // The exit condition of the photon.
        double y_w; // Local y-coordinate intersection on the observer window.
        double z_w; // Local z-coordinate intersection on the observer window.
        double y_s; // Local y-coordinate intersection on the source plane.
        double z_s; // Local z-coordinate intersection on the source plane.
        double final_theta; // Final polar angle $\theta$ at termination.
        double final_phi;   // Final azimuthal angle $\phi$ at termination.
        double L_w; // Affine parameter $\lambda$ at the window intersection.
        double t_w; // Physical coordinate time $t$ at the window intersection.
        double L_s; // Affine parameter $\lambda$ at the source intersection.
        double t_s; // Physical coordinate time $t$ at the source intersection.
    } __attribute__((packed)) blueprint_data_t;

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; // Flattened state vector: 9 components (t, x, y, z, p_t, p_x, p_y, p_z, aux).
        double *f_p; // State vector at the previous integration step.
        double *f_p_p; // State vector at two integration steps prior.
        double *affine_param; // Current affine parameter $\lambda$ for the trajectory.
        double *affine_param_p; // Affine parameter $\lambda$ at the previous step.
        double *affine_param_p_p; // Affine parameter $\lambda$ at two steps prior.
        double *h; // Current adaptive step size for the RKF45 integrator.
        termination_type_t *status; // Current physical/numerical status of the photon.
        int *rejection_retries; // Counter for consecutive RKF45 error tolerance rejections.
        
        // Event Detection State Flags
        bool *on_positive_side_of_window_prev; // True if photon was previously 'above' the window plane.
        bool *on_positive_side_of_source_prev; // True if photon was previously 'above' the source plane.
        
        bool *source_event_found; // Flag indicating a source plane intersection was detected.
        double *source_event_lambda; // Exact affine parameter $\lambda$ of the source intersection.
        double *source_event_f_intersect; // Interpolated 9-component state vector at the source intersection.
        
        bool *window_event_found; // Flag indicating an observer window intersection was detected.
        double *window_event_lambda; // Exact affine parameter $\lambda$ of the window intersection.
        double *window_event_f_intersect; // Interpolated 9-component state vector at the window intersection.
    } PhotonStateSoA;
    """
Bdefines_h.register_BHaH_defines("photon_02_batch_structs", batch_structs_c_code)

par.register_CodeParameters(
    "REAL", __name__,
    ["t_integration_max", "r_escape", "p_t_max", "slot_manager_t_min", "slot_manager_delta_t", "numerical_initial_h"],
    [10000.0, 150.0, 1e3, -1000.0, 10.0, 0.1], commondata=True, add_to_parfile=True,
)
par.register_CodeParameters(
    "bool", __name__, ["perform_conservation_check", "debug_mode"], [True, True], commondata=True, add_to_parfile=True,
)


### STEP 3: C-CODE GENERATION ###
def batch_integrator_numerical(spacetime_name: str) -> None:
    """
    Generate the C orchestrator for the staged batched integration pipeline.

    This function metaprograms a hardware-agnostic C orchestrator that manages 
    the lifecycle of photon batches. It handles device memory allocation (using 
    OpenMP target memory APIs), initializes trajectories based on the specified 
    metric, and executes a synchronized stream-compaction loop for the RKF45 
    integrator stages to maximize hardware occupancy.

    Args:
        spacetime_name (str): The identifier for the spacetime metric (e.g., 'Kerr'), 
                              used for linkage with generated physics kernels.

    Returns:
        None. The generated function is registered via `cfc.register_CFunction`.

    Raises:
        ValueError: If the spacetime_name is not supported by the BHaH environment.

    Example:
        >>> # Typical invocation for a Kerr black hole simulation:
        >>> batch_integrator_numerical("Kerr")
    """
    # 1. Define C-Function metadata in order of appearance
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "omp.h", "math.h"]
    
    desc = r"""@brief Finalized Project Singularity-Axiom Orchestrator (Staged Architecture).
    
    This function serves as the master entry point for the ray-tracing pipeline. 
    It performs the following high-level operations:
    1. Memory Management: Allocates host and device-side SoA memory for 
       num_rays trajectories, strictly following the 9-component f-vector mapping.
    2. Initialization: Calls the spacetime-specific initial conditions.
    3. Staged Integration: Executes a temporal binning loop via the TimeSlotManager 
       and utilizes stream compaction (with a dense-to-sparse reverse map) to safely 
       execute RKF45 stages on dense GPU bundles without OpenMP race conditions.
    4. Event Processing: Computes intersections with bounding geometries.
    5. Integrity Validation: Computes relative errors in conserved quantities."""
    
    cfunc_type = "void"
    name = "batch_integrator_numerical"
    params = "const commondata_struct *restrict commondata, long int num_rays, blueprint_data_t *restrict results_buffer"

    # Dependency function names
    initial_conditions_func = f"set_initial_conditions_cartesian_{spacetime_name}"
    interpolation_func = f"placeholder_interpolation_engine_{spacetime_name}"
    conserved_quantities_func = f"conserved_quantities_{spacetime_name}_photon"


    # 2. Build the C body with internal descriptive comments and the Preamble pattern
    body = f"""
    // Comprehensive Device Pointer Mapping Macro for OpenMP Offloading
    #define SOA_DEVICE_PTRS all_photons.f, all_photons.f_p, all_photons.f_p_p, all_photons.affine_param, all_photons.affine_param_p, all_photons.affine_param_p_p, all_photons.h, all_photons.status, all_photons.rejection_retries, all_photons.on_positive_side_of_window_prev, all_photons.on_positive_side_of_source_prev, all_photons.source_event_found, all_photons.source_event_lambda, all_photons.source_event_f_intersect, all_photons.window_event_found, all_photons.window_event_lambda, all_photons.window_event_f_intersect

    FILE *fp_debug = NULL; // File pointer for outputting sequential trajectory data for post-run analysis.
    if (commondata->debug_mode) {{
        fp_debug = fopen("photon_path_numerical.txt", "w");
        if (fp_debug) fprintf(fp_debug, "# affine_param\\tt\\tx\\ty\\tz\\tp_t\\tp_x\\tp_y\\tp_z\\tL\\n");
    }}

    PhotonStateSoA all_photons_host; // Host-side replica of the flattened master storage representing the 9-component map.
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
    
    TimeSlotManager tsm; // Lock-free temporal arena for binning active rays by physical coordinate time (t).
    slot_manager_init(&tsm, commondata->slot_manager_t_min, commondata->t_start + 1.0, commondata->slot_manager_delta_t, num_rays);

    double window_center[3]; // Cartesian coordinates mapping the geometric center of the camera window.
    double n_x[3]; // Orthonormal basis vector describing the horizontal axis of the camera frame.
    double n_y[3]; // Orthonormal basis vector describing the vertical axis of the camera frame.
    double n_z[3]; // Orthonormal basis vector normal to the camera window (line of sight).
    {initial_conditions_func}(commondata, num_rays, &all_photons_host, window_center, n_x, n_y, n_z);

    // Compute plane scalars for event detection initialization
    const double window_plane_distance = n_z[0]*window_center[0] + n_z[1]*window_center[1] + n_z[2]*window_center[2]; // Scalar distance from origin to the window plane.
    const double source_plane_distance = commondata->source_plane_center_x*commondata->source_plane_normal_x + 
                                  commondata->source_plane_center_y*commondata->source_plane_normal_y + 
                                  commondata->source_plane_center_z*commondata->source_plane_normal_z; // Scalar distance to the physical source plane.

    // TARGET: CPU INITIALIZATION
    // Preparing initial values for all trajectory permutations on host before offloading.
    #pragma omp parallel for
    for (long int i = 0; i < num_rays; i++) {{
        // Preamble: Unpack initial Cartesian coordinates for this ray
        const double initial_x = all_photons_host.f[IDX_GLOBAL(1, i, num_rays)]; // Cartesian x-coordinate
        const double initial_y = all_photons_host.f[IDX_GLOBAL(2, i, num_rays)]; // Cartesian y-coordinate
        const double initial_z = all_photons_host.f[IDX_GLOBAL(3, i, num_rays)]; // Cartesian z-coordinate

        all_photons_host.affine_param[i] = 0.0;
        all_photons_host.h[i] = commondata->numerical_initial_h;
        all_photons_host.status[i] = ACTIVE;
        all_photons_host.rejection_retries[i] = 0;
        
        for(int k=0; k<9; ++k) {{
            const double f_val = all_photons_host.f[IDX_GLOBAL(k, i, num_rays)]; // Base component value representing the unperturbed state.
            all_photons_host.f_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
            all_photons_host.f_p_p[IDX_GLOBAL(k, i, num_rays)] = f_val;
        }}
        
        all_photons_host.affine_param_p[i] = 0.0;
        all_photons_host.affine_param_p_p[i] = 0.0;

        // Evaluate initial plane orientations (Event crossing detection flags)
        const double window_eval = initial_x*n_z[0] + initial_y*n_z[1] + initial_z*n_z[2] - window_plane_distance; // Geometric evaluation against the window normal.
        all_photons_host.on_positive_side_of_window_prev[i] = (window_eval > 0.0);
        
        const double source_eval = initial_x*commondata->source_plane_normal_x + initial_y*commondata->source_plane_normal_y + initial_z*commondata->source_plane_normal_z - source_plane_distance; // Geometric evaluation against the source normal.
        all_photons_host.on_positive_side_of_source_prev[i] = (source_eval > 0.0);
        
        all_photons_host.source_event_found[i] = false;
        all_photons_host.window_event_found[i] = false;
    }}

    /*
     * TARGET: HYBRID EXECUTION DIRECTIVE
     * What: Initializing the results buffer for all individual photons.
     * Why: Ensures clean baseline structs before the device offloading map process begins.
     * How: GPU leverages teams distribute parallel for; CPU falls back to standard OpenMP thread distribution.
     */
    #ifdef USE_GPU
        #pragma omp target teams distribute parallel for map(tofrom: results_buffer[0:num_rays])
    #else
        #pragma omp parallel for
    #endif
    for (long int i = 0; i < num_rays; i++) {{
        results_buffer[i].termination_type = ACTIVE;
        results_buffer[i].y_w = NAN;
        results_buffer[i].z_w = NAN;
        results_buffer[i].y_s = NAN;
        results_buffer[i].z_s = NAN;
        results_buffer[i].final_theta = NAN;
        results_buffer[i].final_phi = NAN;
        results_buffer[i].L_w = NAN;
        results_buffer[i].t_w = NAN;
        results_buffer[i].L_s = NAN;
        results_buffer[i].t_s = NAN;
    }}

    PhotonStateSoA all_photons; // Device-mapped master structure for OpenMP offloading.
    /*
     * TARGET: GPU DEVICE ALLOCATION
     * What: Mapping the unified memory layout from Host CPU to the GPU Device.
     * Why: Crucial for the Structure of Arrays (SoA) layout to fit contiguously inside VRAM for maximum global memory bandwidth.
     * How: Uses omp_target_alloc to create explicit allocations, followed by memcpy commands.
     */
    #ifdef USE_GPU
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

    double *initial_cq = NULL; // Pointer array for storing baseline conserved quantities (E, L, Q) for tracking drift.
    if (commondata->perform_conservation_check) {{
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
    }}

    int initial_slot_idx = slot_get_index(&tsm, commondata->t_start); // Calculates the starting temporal bin index.
    if(initial_slot_idx != -1) {{ 
        for(long int i=0; i<num_rays; ++i) slot_add_photon(&tsm, initial_slot_idx, i); 
    }}

    // Management Arrays and Device-Side Bundle Buffers
    termination_type_t *bundle_status_host = (termination_type_t *)malloc(sizeof(termination_type_t) * BUNDLE_CAPACITY); // Host buffer replicating status array chunk.
    double *bundle_time_host = (double *)malloc(sizeof(double) * BUNDLE_CAPACITY); // Host buffer mapping physical times for binning updates.

    double *f_start_bundle; // Snapshot of the state vector prior to the RKF45 step attempt.
    double *f_temp_bundle;  // Intermediate state vector modified during RKF45 intermediate stages.
    double *f_out_bundle;   // Resulting candidate state vector if the RKF45 step is accepted.
    double *f_err_bundle;   // Calculated local truncation error from the RKF45 embedded method.
    double *metric_bundle;  // Local batch storage for the spacetime metric components.
    double *conn_bundle;    // Local batch storage for the Christoffel symbols (connection coefficients).
    double *k_array_bundle; // 6-stage Runge-Kutta evaluation constants.
    int *req_photon_ids;    // Compacted list of photon indices mapping back to global memory.
    int *compacted_ids;     // Dense list of active threads globally requiring an integration step.
    int *compacted_local_ids; // Reverse map linking dense array indices back to the sparse bundle index.
    bool *needs_step_bundle; // Boolean mask indicating which rays in the local bundle still require processing.
    
    /*
     * TARGET: GPU CACHE-ALIGNED BUNDLE BUFFERS
     * What: Allocating dense execution arrays limited to BUNDLE_CAPACITY.
     * Why: Prevents L1/L2 cache spillage on the GPU, guaranteeing optimal thread execution without divergence.
     */
    #ifdef USE_GPU
        f_start_bundle = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_temp_bundle  = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_out_bundle   = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        f_err_bundle   = (double *)omp_target_alloc(sizeof(double) * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        metric_bundle  = (double *)omp_target_alloc(sizeof(double) * 10 * BUNDLE_CAPACITY, omp_get_default_device());
        conn_bundle    = (double *)omp_target_alloc(sizeof(double) * 40 * BUNDLE_CAPACITY, omp_get_default_device());
        k_array_bundle = (double *)omp_target_alloc(sizeof(double) * 6 * 9 * BUNDLE_CAPACITY, omp_get_default_device());
        req_photon_ids = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        compacted_ids  = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        compacted_local_ids = (int *)omp_target_alloc(sizeof(int) * BUNDLE_CAPACITY, omp_get_default_device());
        needs_step_bundle = (bool *)omp_target_alloc(sizeof(bool) * BUNDLE_CAPACITY, omp_get_default_device());
        
        if (!f_start_bundle || !f_temp_bundle || !compacted_ids || !compacted_local_ids || !needs_step_bundle) {{
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
        compacted_ids  = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
        compacted_local_ids = (int *)malloc(sizeof(int) * BUNDLE_CAPACITY);
        needs_step_bundle = (bool *)malloc(sizeof(bool) * BUNDLE_CAPACITY);
    #endif

    long int active_photons = num_rays; // Counter tracking the global number of photons actively undergoing integration.
    
    for (int i = initial_slot_idx; i >= 0 && active_photons > 0; i--) {{
        while (tsm.slot_counts[i] > 0) {{
            long int bundle_size = MIN(tsm.slot_counts[i], BUNDLE_CAPACITY); // The dynamic size of the current photon batch, capped at hardware capacity limits.
            long int bundle_photons[bundle_size]; // Local array storing the global indices of the sparse photons in the current batch.
            slot_remove_chunk(&tsm, i, bundle_photons, bundle_size);

            /*
             * TARGET: HYBRID EXECUTION DIRECTIVE
             * What: Initializing the local needs_step mask.
             * Why: We assume all photons in a newly pulled bundle require an initial RKF45 step.
             */
            #ifdef USE_GPU
                #pragma omp target teams distribute parallel for is_device_ptr(needs_step_bundle)
            #else
                #pragma omp parallel for
            #endif
            for (long int j = 0; j < bundle_size; j++) {{ needs_step_bundle[j] = true; }}

            // --- STAGED INTEGRATOR (Synchronized Stream Compaction) ---
            bool bundle_has_active = true; // Boolean flag governing continuous iteration until all threads resolve numerical steps.
            
            while (bundle_has_active) {{
                int active_count = 0; // Local counter tracking the dense list of active threads currently requiring an RKF45 integration step.

                // 1. Stream Compaction: Build dense list of work for the GPU
                /*
                 * TARGET: HYBRID EXECUTION DIRECTIVE
                 * What: Stream compaction using OpenMP atomic captures.
                 * Why: Thread divergence heavily penalizes GPU execution. Compacting active threads onto a dense map ensures Warp saturation.
                 */
                #ifdef USE_GPU
                    #pragma omp target teams distribute parallel for map(tofrom: active_count) \\
                                map(to: all_photons, bundle_photons[0:bundle_size]) \\
                                is_device_ptr(needs_step_bundle, compacted_ids, compacted_local_ids) \\
                                has_device_addr(SOA_DEVICE_PTRS)
                #else
                    #pragma omp parallel for
                #endif
                for (long int j = 0; j < bundle_size; j++) {{
                    long int p_idx = bundle_photons[j]; // Global identifier for the evaluated photon.
                    if (needs_step_bundle[j] && all_photons.status[p_idx] == ACTIVE) {{
                        int pos; // Compacted dense array index assigned via hardware atomic increment.
                        #pragma omp atomic capture
                        pos = active_count++;
                        compacted_ids[pos] = (int)p_idx;
                        compacted_local_ids[pos] = (int)j;
                    }}
                }}

                if (active_count == 0) break;

                // 2. Snapshot start state
                #ifdef USE_GPU
                    #pragma omp target teams distribute parallel for \\
                                map(to: all_photons) \\
                                is_device_ptr(f_start_bundle, f_temp_bundle, compacted_ids) \\
                                has_device_addr(SOA_DEVICE_PTRS)
                #else
                    #pragma omp parallel for
                #endif
                for (int i = 0; i < active_count; i++) {{
                    long int p_idx = compacted_ids[i]; // Unpack the mapped global photon index.
                    for (int k = 0; k < 9; k++) {{
                        double val = all_photons.f[IDX_GLOBAL(k, p_idx, num_rays)]; // Base coordinate extraction following the 9-Component Map.
                        f_start_bundle[IDX_LOCAL(k, i, BUNDLE_CAPACITY)] = val;
                        f_temp_bundle[IDX_LOCAL(k, i, BUNDLE_CAPACITY)] = val;
                    }}
                }}

                // 3. RKF45 Stages
                for (int stage = 1; stage <= 6; stage++) {{ // Loop over the 6 intermediate stages of the Runge-Kutta-Fehlberg 4(5) method.
                    {interpolation_func}(commondata, active_count, compacted_ids, f_temp_bundle, metric_bundle, conn_bundle);

                    /*
                     * TARGET: HYBRID EXECUTION DIRECTIVE
                     * What: Invoking the physics core RHS kernels over dense arrays.
                     * Why: We ensure the math ops run over tightly packed IDs so that unmasked threads don't execute useless FLOPS.
                     */
                    #ifdef USE_GPU
                        #pragma omp target teams distribute parallel for \\
                                    map(to: all_photons) \\
                                    is_device_ptr(f_start_bundle, f_temp_bundle, metric_bundle, conn_bundle, k_array_bundle, compacted_ids) \\
                                    has_device_addr(SOA_DEVICE_PTRS)
                    #else
                        #pragma omp parallel for
                    #endif
                    for (int i = 0; i < active_count; i++) {{
                        long int p_idx = compacted_ids[i]; // Unpack mapped ID.
                        calculate_ode_rhs(f_temp_bundle, metric_bundle, conn_bundle, k_array_bundle, BUNDLE_CAPACITY, stage, i);
                        if (stage < 6) {{
                            calculate_rkf45_stage_f_temp(stage + 1, f_start_bundle, k_array_bundle, all_photons.h[p_idx], f_temp_bundle, BUNDLE_CAPACITY, i);
                        }}
                    }}
                }}

                // 4. Masked Acceptance Logic
                bundle_has_active = false;
                /*
                 * TARGET: HYBRID EXECUTION DIRECTIVE
                 * What: Computing acceptance of intermediate RKF45 steps using a logical OR reduction.
                 * Why: Allows adaptive step sizing to bound relative numerical errors, rejecting steps across threads seamlessly.
                 */
                #ifdef USE_GPU
                    #pragma omp target teams distribute parallel for reduction(|:bundle_has_active) \\
                                map(to: all_photons, commondata[0:1]) \\
                                is_device_ptr(f_start_bundle, f_out_bundle, f_err_bundle, k_array_bundle, needs_step_bundle, compacted_ids, compacted_local_ids) \\
                                has_device_addr(SOA_DEVICE_PTRS)
                #else
                    #pragma omp parallel for reduction(|:bundle_has_active)
                #endif
                for (int i = 0; i < active_count; i++) {{
                    long int p_idx = compacted_ids[i]; // Unpack global ID.
                    int j = compacted_local_ids[i]; // Reverse map referencing the sparse bundle index.
                    
                    double h_step = all_photons.h[p_idx]; // Adaptive step size assigned per-trajectory.
                    rkf45_kernel(f_start_bundle, k_array_bundle, h_step, f_out_bundle, f_err_bundle, BUNDLE_CAPACITY, i);
                    bool accepted = update_photon_state_and_stepsize(&all_photons, num_rays, p_idx, f_start_bundle, f_out_bundle, f_err_bundle, BUNDLE_CAPACITY, i, commondata); // Validates L-infinity limits.

                    if (!accepted) {{
                        if (all_photons.rejection_retries[p_idx] > commondata->rkf45_max_retries) {{
                            all_photons.status[p_idx] = FAILURE_RKF45_REJECTION_LIMIT;
                            needs_step_bundle[j] = false;
                        }} else {{
                            bundle_has_active = true; 
                        }}
                    }} else {{
                        needs_step_bundle[j] = false; // Step bounded correctly within structural limits.
                        for (int k = 0; k < 9; ++k) {{
                            all_photons.f_p_p[IDX_GLOBAL(k, p_idx, num_rays)] = all_photons.f_p[IDX_GLOBAL(k, p_idx, num_rays)];
                            all_photons.f_p[IDX_GLOBAL(k, p_idx, num_rays)] = f_start_bundle[IDX_LOCAL(k, i, BUNDLE_CAPACITY)];
                        }}
                        all_photons.affine_param_p_p[p_idx] = all_photons.affine_param_p[p_idx];
                        all_photons.affine_param_p[p_idx] = all_photons.affine_param[p_idx] - h_step;
                    }}
                }}
            }} // end while (bundle_has_active)
            // --- END STAGED INTEGRATOR ---
            
            // --- EVENT DETECTION & TERMINATION CHECKS ---
            /*
             * TARGET: HYBRID EXECUTION DIRECTIVE
             * What: Event evaluation based on the final updated geometry properties.
             * Why: We decouple root-finding events from stream-compaction physics to prevent divergent branching.
             */
            #ifdef USE_GPU
                #pragma omp target teams distribute parallel for \\
                            map(to: all_photons, commondata[0:1], bundle_photons[0:bundle_size]) \\
                            has_device_addr(SOA_DEVICE_PTRS) \\
                            map(tofrom: results_buffer[0:num_rays])
            #else
                #pragma omp parallel for
            #endif
            for (long int j = 0; j < bundle_size; ++j) {{
                long int p_idx = bundle_photons[j]; // The global index of the evaluated photon.
                if (all_photons.status[p_idx] == ACTIVE) {{
                    
                    // Preamble: Unpack structural arrays to highly descriptive local variables.
                    const double pos_x = all_photons.f[IDX_GLOBAL(1, p_idx, num_rays)]; // The current Cartesian x-coordinate of the photon.
                    const double pos_y = all_photons.f[IDX_GLOBAL(2, p_idx, num_rays)]; // The current Cartesian y-coordinate of the photon.
                    const double pos_z = all_photons.f[IDX_GLOBAL(3, p_idx, num_rays)]; // The current Cartesian z-coordinate of the photon.
                    const double p_t_val = all_photons.f[IDX_GLOBAL(4, p_idx, num_rays)]; // Temporal momentum component $p_t$, dictating conserved energy metric.
                    const double t_val = all_photons.f[IDX_GLOBAL(0, p_idx, num_rays)]; // The physical coordinate time $t$ of the trajectory.
                    
                    const double r_sq = pos_x*pos_x + pos_y*pos_y + pos_z*pos_z; // Squared radial distance defining celestial escape bounds.
                    
                    if (fabs(p_t_val) > commondata->p_t_max) {{
                        all_photons.status[p_idx] = FAILURE_PT_TOO_BIG;
                    }} else if (fabs(t_val) > commondata->t_integration_max) {{
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
                            if (isnan(results_buffer[p_idx].t_w)) {{
                                handle_window_plane_intersection(&all_photons, num_rays, p_idx, commondata, &results_buffer[p_idx]); 
                            }}
                            all_photons.window_event_found[p_idx] = false; 
                        }}
                    }}
                }}
            }}

            // --- SYNC STATUS TO HOST ---
            #ifdef USE_GPU
                for (long int j = 0; j < bundle_size; j++) {{
                    long int p_idx = bundle_photons[j]; // Extracting mapping index.
                    omp_target_memcpy(&bundle_status_host[j], &all_photons.status[p_idx], sizeof(termination_type_t), 0, 0, omp_get_initial_device(), omp_get_default_device());
                    omp_target_memcpy(&bundle_time_host[j], &all_photons.f[IDX_GLOBAL(0, p_idx, num_rays)], sizeof(double), 0, 0, omp_get_initial_device(), omp_get_default_device());
                }}
            #else
                for (long int j = 0; j < bundle_size; j++) {{
                    bundle_status_host[j] = all_photons.status[bundle_photons[j]];
                    bundle_time_host[j] = all_photons.f[IDX_GLOBAL(0, bundle_photons[j], num_rays)];
                }}
            #endif

            long int active_photons_decrement = 0; // Cumulative reduction counter.
            long int active_in_bundle = 0; // Retracking photons looping in current slot.
            
            for (long int j = 0; j < bundle_size; j++) {{
                long int p_idx = bundle_photons[j]; // Map check ID.
                if (bundle_status_host[j] == ACTIVE) {{
                    int s_idx = slot_get_index(&tsm, bundle_time_host[j]); // Calculated index logic binding coordinate time.
                    if (s_idx != -1) {{
                        slot_add_photon(&tsm, s_idx, p_idx);
                        active_in_bundle++;
                   }} else {{
                        #ifdef USE_GPU
                            termination_type_t fail_stat = FAILURE_SLOT_MANAGER_ERROR; // Explicit bounds break fault logic.
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

    /*
     * TARGET: HYBRID EXECUTION DIRECTIVE
     * What: Finalization routing where raw components are evaluated into abstract representations.
     * Why: Pushes blueprint derivations (final polar angles, intercepts) on device arrays right before concluding.
     */
    #ifdef USE_GPU
        #pragma omp target teams distribute parallel for \\
                    map(to: all_photons, commondata[0:1]) \\
                    map(tofrom: results_buffer[0:num_rays]) \\
                    has_device_addr(SOA_DEVICE_PTRS)
    #else
        #pragma omp parallel for
    #endif
    for (long int i = 0; i < num_rays; i++) {{ calculate_and_fill_blueprint_data_universal(&all_photons, num_rays, i, &results_buffer[i]); }}

    #ifdef USE_GPU
        double *f_host = (double*)malloc(sizeof(double) * 9 * num_rays); // Host-side extraction buffer for validation outputs.
        omp_target_memcpy(f_host, all_photons.f, sizeof(double) * 9 * num_rays, 0, 0, omp_get_initial_device(), omp_get_default_device());
    #else
        double *f_host = all_photons.f;
    #endif

    // --- FINAL CONSERVATION INTEGRITY REPORT ---
    // If enabled, we perform a final pass over the terminated rays to quantify 
    // numerical drift in the constants of motion (E, L, Q).
    if (commondata->perform_conservation_check && initial_cq != NULL) {{
        printf("\\n--- Final Conservation Summary ---\\n"); 
        
        // These variables track the single worst-performing photon for each conserved quantity.
        double max_dE = 0.0, max_dL = 0.0, max_dQ = 0.0;
        long int worst_E_idx = -1, worst_L_idx = -1, worst_Q_idx = -1;

        FILE *fp_cons_log = fopen("conservation_errors.txt", "w");
        if (fp_cons_log) fprintf(fp_cons_log, "# photon_idx dE_rel dL_rel dQ_rel Q_init Q_final term_type final_t\\n");

        #pragma omp parallel
        {{
            // Thread-local variables to avoid frequent atomic updates/contention.
            double local_max_dE = 0.0, local_max_dL = 0.0, local_max_dQ = 0.0;
            long int local_worst_E = -1, local_worst_L = -1, local_worst_Q = -1;

            #pragma omp for
            for (long int i = 0; i < num_rays; i++) {{
                double final_E, final_Lx, final_Ly, final_Lz, final_Q;
                {conserved_quantities_func}(commondata, f_host, num_rays, i, &final_E, &final_Lx, &final_Ly, &final_Lz, &final_Q);

                // 1. Energy Conservation (E = -p_t)
                double init_E = initial_cq[i * 5 + 0];
                double dE = (init_E != 0) ? fabs((final_E - init_E) / init_E) : fabs(final_E - init_E);
                if (dE > local_max_dE) {{ local_max_dE = dE; local_worst_E = i; }}

                // 2. Angular Momentum Conservation (L)
                double dL;
                if (commondata->a_spin == 0.0) {{
                    // For Schwarzschild, the total magnitude |L| is conserved.
                    double L_init_mag = sqrt(initial_cq[i*5 + 1]*initial_cq[i*5 + 1] + initial_cq[i*5 + 2]*initial_cq[i*5 + 2] + initial_cq[i*5 + 3]*initial_cq[i*5 + 3]);
                    double L_final_mag = sqrt(final_Lx*final_Lx + final_Ly*final_Ly + final_Lz*final_Lz);
                    dL = (L_init_mag != 0) ? fabs((L_final_mag - L_init_mag) / L_init_mag) : fabs(L_final_mag - L_init_mag);
                }} else {{
                    // For Kerr, only the axial component L_z is strictly conserved.
                    double L_z_init = initial_cq[i * 5 + 3];
                    dL = (L_z_init != 0) ? fabs((final_Lz - L_z_init) / L_z_init) : fabs(final_Lz - L_z_init);
                }}
                if (dL > local_max_dL) {{ local_max_dL = dL; local_worst_L = i; }}

                // 3. Carter Constant Conservation (Q)
                double init_Q = initial_cq[i * 5 + 4];
                double dQ = (init_Q != 0) ? fabs((final_Q - init_Q) / init_Q) : fabs(final_Q - init_Q);
                if (dQ > local_max_dQ) {{ local_max_dQ = dQ; local_worst_Q = i; }}
                
                // Optional: Write every photon's error to a log for post-processing/histograms.
                if (fp_cons_log) {{
                    #pragma omp critical
                    fprintf(fp_cons_log, "%ld %.6e %.6e %.6e %.6e %.6e %d %.6e\\n", 
                            i, dE, dL, dQ, init_Q, final_Q, (int)all_photons.status[i], f_host[IDX_GLOBAL(0, i, num_rays)]);
                }}
            }}
            
            // Critical section to merge thread-local 'worst' results into the global maximums.
            #pragma omp critical
            {{
                if (local_max_dE > max_dE) {{ max_dE = local_max_dE; worst_E_idx = local_worst_E; }}
                if (local_max_dL > max_dL) {{ max_dL = local_max_dL; worst_L_idx = local_worst_L; }}
                if (local_max_dQ > max_dQ) {{ max_dQ = local_max_dQ; worst_Q_idx = local_worst_Q; }}
            }}
        }}
        if (fp_cons_log) fclose(fp_cons_log);
        
        // Output the peak errors found across the entire ray batch to stdout.
        printf("Max relative error in Energy (E): %.3e (photon %ld)\\n", max_dE, worst_E_idx);
        printf("Max relative error in Ang. Mom. (L): %.3e (photon %ld)\\n", max_dL, worst_L_idx);
        printf("Max relative error in Carter Constant (Q): %.3e (photon %ld)\\n", max_dQ, worst_Q_idx);
        printf("------------------------------------\\n");
    }}

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
        omp_target_free(compacted_ids, omp_get_default_device());
        omp_target_free(compacted_local_ids, omp_get_default_device()); // Prompt 2 Patch Application
        omp_target_free(needs_step_bundle, omp_get_default_device());

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
        free(compacted_ids); free(compacted_local_ids); free(needs_step_bundle); // Prompt 2 Patch Application
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
    
    # 3. Register the function with the NRPy environment following strict formatting template
    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type, 
        name=name,
        params=params,
        body=body
    )
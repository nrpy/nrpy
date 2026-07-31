# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/set_initial_conditions_kernel.py
"""
Defines the shared metric-driven C initializer for photon batches.

Analytical and numerical spacetime paths both provide one covariant observer
metric evaluated at the common observer event.  The initializer builds one
metric-orthonormal observer tetrad per tile-batch call, initializes every ray
in that tile with a complete past-directed null four-momentum, and uses each
ray's normalized image-sample coordinates for ray construction. Batch struct
definitions are registered separately so analytical single-ray code can retain
its smaller state layout.

Single coalesced memory writes prevent thread serialization and ensure aligned cache
access. An explicit hardware error synchronization trap prevents silent link-time
symbol failures caused by compiling with -rdc=true. Hydrating pinned memory via data
bus seeds the Time Slot Manager. Evaluating the initial side of the independent
nonterminal and terminal planes natively prevents redundant device memory allocation
and data transfers.
Thread identification boundaries prevent out-of-bounds access for threads exceeding
the active chunk. Parallelized batch processing distributes execution across threads.
Processing memory in static bundles protects hardware limits. A synchronization
transfer updates the master Structure of Arrays state.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.general_relativity.geodesics.geodesics as geo
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.infrastructures.BHaH.BHaH_defines_h as Bdefines_h
import nrpy.params as par
from nrpy.helpers.loop import loop

# Define the C-structs required for the simulation pipeline.
# These are registered here to ensure they appear in BHaH_defines.h before
# the initialization kernel is compiled.
batch_structs_c_code = r"""
    #include <stddef.h>

    // Maximum number of photons processed per batch to fit within L1/L2 cache.
    #define BUNDLE_CAPACITY 524288

    // Defines the physical event surfaces tracked during integration.
    typedef enum {
        NON_TERMINAL_PLANE_EVENT, // Intersection with the nonterminal plane.
        TERMINAL_PLANE_EVENT  // Intersection with the terminal plane.
    } event_type_t; // END ENUM: event_type_t

    // Defines the specific exit condition for a photon's integration loop.
    typedef enum {
        STOP_CONDITION_COORD_RADIUS_EXCEEDED, // 0: Coordinate-radius stop condition was reached.
        STOP_CONDITION_TERMINAL_PLANE, // 1: Photon hit the terminal plane.
        STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED, // 2: Evolution-measure stop condition was reached.
        FAILURE_RKF45_REJECTION_LIMIT, // 3: Adaptive step size was rejected too many times.
        STOP_CONDITION_T_MAX_EXCEEDED, // 4: Maximum physical coordinate time was exceeded.
        FAILURE_SLOT_MANAGER_ERROR, // 5: TimeSlotManager failed to handle the photon.
        FAILURE_GENERIC, // 6: Generic unclassified numerical failure.
        ACTIVE, // 7: Photon is currently undergoing integration.
        REJECTED, // 8: Photon RKF45 step was rejected.
    } termination_type_t; // END ENUM: termination_type_t

    // Native same-build metadata for one serialized blueprint tile.
    #define BLUEPRINT_MAGIC "NRPYBP01"
    // Schema v6 stores fields of view once in the tile header and stores each
    // ray's normalized image sample directly in the record.
    #define BLUEPRINT_SCHEMA_VERSION 6U
    typedef struct {
        char magic[8];
        uint32_t native_schema_version;
        uint32_t header_size;
        uint32_t record_size;
        uint32_t tx;
        uint32_t ty;
        uint32_t tiles_w;
        uint32_t tiles_h;
        uint64_t record_count;
        double alpha_w; // Width-direction field of view in radians.
        double alpha_h; // Height-direction field of view in radians.
    } __attribute__((packed)) blueprint_header_t; // END STRUCT: blueprint_header_t

    // Stores the final physical properties of a photon upon integration termination.
    typedef struct {
        int32_t termination_type; // Fixed-width serialized exit condition.
        double y_nt; // Local horizontal coordinate on nonterminal plane.
        double z_nt; // Local vertical coordinate on nonterminal plane.
        double y_t; // Local horizontal coordinate on terminal plane.
        double z_t; // Local vertical coordinate on terminal plane.
        double final_theta; // Final polar angle $\theta$ at termination.
        double final_phi;   // Final azimuthal angle $\phi$ at termination.
        double non_terminal_plane_lambda; // Affine parameter at nonterminal-plane intersection.
        double non_terminal_plane_t; // Coordinate time at nonterminal-plane intersection.
        double L_f; // Affine parameter $\lambda$ when the photon terminated.
        double t_f; // Physical coordinate time $t$ when the photon terminated.
        double image_width_fraction; // Normalized width coordinate in [0,1].
        double image_height_fraction; // Normalized height coordinate in [0,1].
    } __attribute__((packed)) blueprint_data_t; // END STRUCT: blueprint_data_t

#if defined(__cplusplus)
    static_assert(sizeof(blueprint_header_t) == 60, "blueprint_header_t size changed");
    static_assert(offsetof(blueprint_header_t, magic) == 0, "blueprint header magic offset changed");
    static_assert(offsetof(blueprint_header_t, native_schema_version) == 8, "blueprint header version offset changed");
    static_assert(offsetof(blueprint_header_t, header_size) == 12, "blueprint header size offset changed");
    static_assert(offsetof(blueprint_header_t, record_size) == 16, "blueprint header record size offset changed");
    static_assert(offsetof(blueprint_header_t, tx) == 20, "blueprint header tx offset changed");
    static_assert(offsetof(blueprint_header_t, ty) == 24, "blueprint header ty offset changed");
    static_assert(offsetof(blueprint_header_t, tiles_w) == 28, "blueprint header tiles_w offset changed");
    static_assert(offsetof(blueprint_header_t, tiles_h) == 32, "blueprint header tiles_h offset changed");
    static_assert(offsetof(blueprint_header_t, record_count) == 36, "blueprint header count offset changed");
    static_assert(offsetof(blueprint_header_t, alpha_w) == 44, "blueprint alpha_w offset changed");
    static_assert(offsetof(blueprint_header_t, alpha_h) == 52, "blueprint alpha_h offset changed");
    static_assert(sizeof(blueprint_data_t) == 100, "blueprint_data_t size changed");
    static_assert(offsetof(blueprint_data_t, termination_type) == 0, "blueprint termination offset changed");
    static_assert(offsetof(blueprint_data_t, y_nt) == 4, "blueprint y_nt offset changed");
    static_assert(offsetof(blueprint_data_t, z_nt) == 12, "blueprint z_nt offset changed");
    static_assert(offsetof(blueprint_data_t, y_t) == 20, "blueprint y_t offset changed");
    static_assert(offsetof(blueprint_data_t, z_t) == 28, "blueprint z_t offset changed");
    static_assert(offsetof(blueprint_data_t, final_theta) == 36, "blueprint final_theta offset changed");
    static_assert(offsetof(blueprint_data_t, final_phi) == 44, "blueprint final_phi offset changed");
    static_assert(offsetof(blueprint_data_t, non_terminal_plane_lambda) == 52, "blueprint non_terminal_plane_lambda offset changed");
    static_assert(offsetof(blueprint_data_t, non_terminal_plane_t) == 60, "blueprint non_terminal_plane_t offset changed");
    static_assert(offsetof(blueprint_data_t, L_f) == 68, "blueprint L_f offset changed");
    static_assert(offsetof(blueprint_data_t, t_f) == 76, "blueprint t_f offset changed");
    static_assert(offsetof(blueprint_data_t, image_width_fraction) == 84, "blueprint width fraction offset changed");
    static_assert(offsetof(blueprint_data_t, image_height_fraction) == 92, "blueprint height fraction offset changed");
#else
    _Static_assert(sizeof(blueprint_header_t) == 60, "blueprint_header_t size changed");
    _Static_assert(offsetof(blueprint_header_t, magic) == 0, "blueprint header magic offset changed");
    _Static_assert(offsetof(blueprint_header_t, native_schema_version) == 8, "blueprint header version offset changed");
    _Static_assert(offsetof(blueprint_header_t, header_size) == 12, "blueprint header size offset changed");
    _Static_assert(offsetof(blueprint_header_t, record_size) == 16, "blueprint header record size offset changed");
    _Static_assert(offsetof(blueprint_header_t, tx) == 20, "blueprint header tx offset changed");
    _Static_assert(offsetof(blueprint_header_t, ty) == 24, "blueprint header ty offset changed");
    _Static_assert(offsetof(blueprint_header_t, tiles_w) == 28, "blueprint header tiles_w offset changed");
    _Static_assert(offsetof(blueprint_header_t, tiles_h) == 32, "blueprint header tiles_h offset changed");
    _Static_assert(offsetof(blueprint_header_t, record_count) == 36, "blueprint header count offset changed");
    _Static_assert(offsetof(blueprint_header_t, alpha_w) == 44, "blueprint alpha_w offset changed");
    _Static_assert(offsetof(blueprint_header_t, alpha_h) == 52, "blueprint alpha_h offset changed");
    _Static_assert(sizeof(blueprint_data_t) == 100, "blueprint_data_t size changed");
    _Static_assert(offsetof(blueprint_data_t, termination_type) == 0, "blueprint termination offset changed");
    _Static_assert(offsetof(blueprint_data_t, y_nt) == 4, "blueprint y_nt offset changed");
    _Static_assert(offsetof(blueprint_data_t, z_nt) == 12, "blueprint z_nt offset changed");
    _Static_assert(offsetof(blueprint_data_t, y_t) == 20, "blueprint y_t offset changed");
    _Static_assert(offsetof(blueprint_data_t, z_t) == 28, "blueprint z_t offset changed");
    _Static_assert(offsetof(blueprint_data_t, final_theta) == 36, "blueprint final_theta offset changed");
    _Static_assert(offsetof(blueprint_data_t, final_phi) == 44, "blueprint final_phi offset changed");
    _Static_assert(offsetof(blueprint_data_t, non_terminal_plane_lambda) == 52, "blueprint non_terminal_plane_lambda offset changed");
    _Static_assert(offsetof(blueprint_data_t, non_terminal_plane_t) == 60, "blueprint non_terminal_plane_t offset changed");
    _Static_assert(offsetof(blueprint_data_t, L_f) == 68, "blueprint L_f offset changed");
    _Static_assert(offsetof(blueprint_data_t, t_f) == 76, "blueprint t_f offset changed");
    _Static_assert(offsetof(blueprint_data_t, image_width_fraction) == 84, "blueprint width fraction offset changed");
    _Static_assert(offsetof(blueprint_data_t, image_height_fraction) == 92, "blueprint height fraction offset changed");
#endif

    // ==========================================
    // Flattened SoA Struct (Master Storage)
    // ==========================================
    typedef struct {
        double *f; // Flattened state vector: t, x, y, z, mode-specific energy, p^x, p^y, p^z, aux.
        double *f_p; // State vector at the previous integration step.
        double *f_p_p; // State vector at two integration steps prior.
        double *integration_param; // Current mode-dependent integration parameter.
        double *integration_param_p; // Integration parameter at the previous step.
        double *integration_param_p_p; // Integration parameter at two steps prior.
        double *h; // Current adaptive step size $h$ for the RKF45 integrator.
        termination_type_t *status; // Current physical/numerical status of the photon.
        int *rejection_retries; // Counter for consecutive RKF45 error tolerance rejections.

        // Event Detection State Flags (Persistence Layer for Batch C)
        bool *on_positive_side_of_non_terminal_plane_prev; // Previous side of nonterminal plane.
        bool *on_positive_side_of_terminal_plane_prev; // True if photon was previously 'above' the terminal plane.

        bool *terminal_plane_event_found; // Flag indicating a terminal plane intersection was detected.
        double *terminal_plane_event_lambda; // Exact affine parameter $\lambda$ at terminal crossing.
        double *terminal_plane_event_f_intersect; // State vector at terminal crossing.

        bool *non_terminal_plane_event_found; // Nonterminal-plane intersection lock.
        double *non_terminal_plane_event_lambda; // Affine parameter $\lambda$ at nonterminal plane.
        double *non_terminal_plane_event_f_intersect; // State at nonterminal-plane intersection.
    } PhotonStateSoA; // END STRUCT: PhotonStateSoA
"""


def register_photon_batch_structs() -> None:
    """
    Register the shared photon batch structs and blueprint schema.

    Batch and single-ray generators call this before registering code that
    consumes the shared PhotonStateSoA. The shared definition includes the
    event-side and interpolated-intersection fields used by the initializer and
    event-detection kernels.
    """
    if "photon_batch_structs" not in par.glb_extras_dict.get("BHaH_defines", {}):
        Bdefines_h.register_BHaH_defines("photon_batch_structs", batch_structs_c_code)


def set_initial_conditions_kernel(normalized_eom: bool = False) -> None:
    """
    Register shared photon initialization from one observer-event metric.

    The caller must provide the ten independent covariant metric components
    ``(g00, g01, g02, g03, g11, g12, g13, g22, g23, g33)`` evaluated at the
    observer position and ``t_start``.  This function constructs one observer
    tetrad on the host, validates it, and passes its sixteen contravariant
    components to every ray in the current tile.

    Direct momentum state convention before normalized conversion:

    ``f[4] = p^0`` and ``f[5:8] = (p^1, p^2, p^3)``.

    For normalized evolution, the existing downstream conversion replaces
    these four slots with ``u = log(abs(alpha*p^0))`` and ``Pi_i``.  The
    initializer therefore always constructs the complete direct momentum
    first.  No algebraic temporal-momentum recovery is needed.

    :param normalized_eom: Whether to initialize the normalized photon state
        layout.
    """
    # Step 1: Register tile-sampling state.  Tile indices are the only mutable
    # tile state.  Each tile origin is derived below from the active indices;
    # no origin or tile pixel dimension is stored in commondata.
    par.register_CodeParameters(
        "int",
        __name__,
        [
            "tile_index_width",
            "tile_index_height",
        ],
        [0, 0],
        commondata=True,
        add_to_parfile=False,
    )
    par.register_CodeParameter(
        "int",
        __name__,
        "scan_density",
        1,
        commondata=True,
        add_to_parfile=True,
        description="Width-side ray-sample count per tile.",
    )
    par.register_CodeParameters(
        "REAL",
        __name__,
        [
            "observer_x",
            "observer_y",
            "observer_z",
            "observer_look_forward_x",
            "observer_look_forward_y",
            "observer_look_forward_z",
            "observer_up_x",
            "observer_up_y",
            "observer_up_z",
            "alpha_w",
            "alpha_h",
            "initial_h",
            "t_start",
        ],
        [
            51.0,
            0.0,
            0.0,
            -1.0,
            0.0,
            0.0,
            0.0,
            0.0,
            1.0,
            1.0,
            1.0,
            0.1,
            100.0,
        ],
        commondata=True,
        add_to_parfile=True,
    )
    # Plane geometry parameters are owned by the corresponding intersection
    # handlers. The initializer only consumes their commondata values when it
    # initializes event-side history for batch paths.

    # Step 2: Pull universal C and q expressions from the equation layer.  The
    # generated ray loop supplies screen offsets through C variables named h
    # and v, matching the symbols declared by geodesics.py.
    observer_C_expr, observer_q_expr = (
        geo.GeodesicEquations.photon_observer_ray_C_and_q()
    )
    observer_ray_math = ccg.c_codegen(
        [observer_C_expr, observer_q_expr],
        ["ray_C", "ray_q"],
        enable_cse=True,
        enable_simd=False,
        verbose=False,
        include_braces=False,
    )

    parallelization = par.parval_from_str("parallelization")
    cd_access = parallel_utils.get_commondata_access(parallelization)
    initial_state_parameter = "0.0" if normalized_eom else f"{cd_access}t_start"
    normalized_tracker_initialization = (
        """
    // Normalized evolution uses coordinate time as its external integration
    // parameter.  The affine parameter remains in f[0].
    for (long int ray = 0; ray < num_rays; ++ray) {
        all_photons->integration_param[ray] = commondata->t_start;
        all_photons->integration_param_p[ray] = commondata->t_start;
        all_photons->integration_param_p_p[ray] = commondata->t_start;
    } // END LOOP: initialize normalized coordinate-time trackers
"""
        if normalized_eom
        else ""
    )

    # Step 3: Describe kernel arguments.  Tetrad components are passed as
    # scalars so CUDA kernels do not need a device pointer to host stack data.
    tetrad_argument_names = [
        f"observer_tetrad_{basis}_{mu}" for basis in range(4) for mu in range(4)
    ]
    arg_dict = {
        "num_rays": "const long int",
        "d_f_bundle": "double *restrict",
        "d_h_bundle": "double *restrict",
        "scan_density_height": "const int",
        "start_idx": "const long int",
        "chunk_size": "const long int",
    }
    arg_dict.update({name: "const double" for name in tetrad_argument_names})
    if parallelization != "cuda":
        arg_dict["commondata"] = "const commondata_struct *restrict"

    # Step 4: Select architecture-specific thread setup.
    if parallelization == "cuda":
        loop_preamble = r"""
    //==========================================
    // CUDA THREAD IDENTIFICATION
    //==========================================
    const long int c = blockIdx.x * blockDim.x + threadIdx.x;
    if (c >= chunk_size) return;
"""
        loop_postamble = ""
    else:
        loop_preamble = r"""
    //==========================================
    // OPENMP LOOP ARCHITECTURE
    //==========================================
    #pragma omp parallel for
    for (long int c = 0; c < chunk_size; ++c) {
"""
        loop_postamble = "    } // END LOOP: initialize current ray bundle"

    # Step 5: Build per-ray initialization.  Only pixel coordinates, C/q, and
    # the linear combination of the precomputed tetrad occur in this loop.
    tetrad_loads = "\n".join(
        f"    const double e{basis}[4] = "
        f"{{observer_tetrad_{basis}_0, observer_tetrad_{basis}_1, "
        f"observer_tetrad_{basis}_2, observer_tetrad_{basis}_3}};"
        for basis in range(4)
    )
    core_math = r"""
    //==========================================
    // MACRO DEFINITIONS
    //==========================================
    #define IDX_F(component, ray_id) ((component) * BUNDLE_CAPACITY + (ray_id))
    #define IDX_H(ray_id) (ray_id)

    //==========================================
    // OBSERVER TETRAD LOAD
    //==========================================
    // Each e_a array stores contravariant components e_a^mu in ordering
    // (t, x, y, z).  These values were constructed once from g_mu_nu at the
    // observer event and are identical for every ray in this bundle.
__TETRAD_LOADS__

    //==========================================
    // PIXEL CENTER AND ANGULAR OFFSETS
    //==========================================
    // Tile indices partition the common normalized image domain [0,1] in each
    // direction.  The host validates the derived height-side sample count and
    // passes it as a scalar so CUDA and CPU use identical sample ordering.
    const int scan_density_width = __CD_ACCESS__scan_density;
    const int local_col = (int)(c % scan_density_width);
    const int local_row = (int)(c / scan_density_width);
    const double a =
        ((double)__CD_ACCESS__tile_index_width +
         ((double)local_col + 0.5) / (double)scan_density_width) /
        (double)__CD_ACCESS__tiles_width;
    const double b =
        ((double)__CD_ACCESS__tile_index_height +
         ((double)local_row + 0.5) / (double)scan_density_height) /
        (double)__CD_ACCESS__tiles_height;

    // Pixel-center sampling uses half-cell offsets. Increasing the column
    // increases the width offset and therefore moves the ray toward e_3;
    // increasing the row increases the height offset and moves it toward e_2.
    const double h = (2.0 * a - 1.0) * tan(0.5 * __CD_ACCESS__alpha_w);
    const double v = (2.0 * b - 1.0) * tan(0.5 * __CD_ACCESS__alpha_h);

    // Generate C and q from geodesics.py.  ray_C enforces nullness; ray_q
    // sets the observer-frame energy magnitude to one.
    double ray_C = 0.0;
    double ray_q = 0.0;
__OBSERVER_RAY_MATH__

    //==========================================
    // INITIAL POSITION AND COMPLETE MOMENTUM
    //==========================================
    d_f_bundle[IDX_F(0, c)] = __INITIAL_STATE_PARAMETER__; // lambda or t.
    d_f_bundle[IDX_F(1, c)] = __CD_ACCESS__observer_x;
    d_f_bundle[IDX_F(2, c)] = __CD_ACCESS__observer_y;
    d_f_bundle[IDX_F(3, c)] = __CD_ACCESS__observer_z;
    // Reverse ray tracing evolves the past-directed vector
    // k^mu = -e_0^mu + (e_1^mu + v e_2^mu + h e_3^mu)/C.
    // Writing -ray_q*ray_C exposes the paper's q chi construction and avoids
    // solving a separate quadratic for p^0.
    const double past_e0_coefficient = -ray_q * ray_C;
    d_f_bundle[IDX_F(4, c)] =
        past_e0_coefficient * e0[0] +
        ray_q * (e1[0] + v * e2[0] + h * e3[0]);
    d_f_bundle[IDX_F(5, c)] =
        past_e0_coefficient * e0[1] +
        ray_q * (e1[1] + v * e2[1] + h * e3[1]);
    d_f_bundle[IDX_F(6, c)] =
        past_e0_coefficient * e0[2] +
        ray_q * (e1[2] + v * e2[2] + h * e3[2]);
    d_f_bundle[IDX_F(7, c)] =
        past_e0_coefficient * e0[3] +
        ray_q * (e1[3] + v * e2[3] + h * e3[3]);

    // Affine-distance diagnostic starts at zero in both direct and normalized
    // state layouts.
    d_f_bundle[IDX_F(8, c)] = 0.0;
    d_h_bundle[IDX_H(c)] = __CD_ACCESS__initial_h;

    #undef IDX_F
    #undef IDX_H
""".replace("__TETRAD_LOADS__", tetrad_loads)
    core_math = core_math.replace("__CD_ACCESS__", cd_access)
    core_math = core_math.replace(
        "__INITIAL_STATE_PARAMETER__", initial_state_parameter
    )
    core_math = core_math.replace("__OBSERVER_RAY_MATH__", observer_ray_math)
    kernel_body = f"{loop_preamble}\n{core_math}\n{loop_postamble}"

    # Step 6: Register the generated per-ray kernel and launch metadata.
    launch_dict = (
        {
            "threads_per_block": ["BHAH_THREADS_IN_X_DIR_DEFAULT", "1", "1"],
            "blocks_per_grid": [
                "(chunk_size + BHAH_THREADS_IN_X_DIR_DEFAULT - 1) / "
                "BHAH_THREADS_IN_X_DIR_DEFAULT",
                "1",
                "1",
            ],
            "stream": "stream_idx",
        }
        if parallelization == "cuda"
        else None
    )
    kernel_prefunc, launch_code = parallel_utils.generate_kernel_and_launch_code(
        kernel_name="set_initial_conditions_kernel",
        kernel_body=kernel_body,
        arg_dict_cuda=arg_dict,
        arg_dict_host=arg_dict,
        parallelization=parallelization,
        launch_dict=launch_dict,
        cfunc_decorators="__global__" if parallelization == "cuda" else "",
    )

    # Step 7: Generate host/device staging transfers.  Host validation after
    # each transfer fails immediately if any initialized state is non-finite.
    finite_check = r"""
        for (long int check_c = 0; check_c < chunk_size; ++check_c) {
            const long int check_i = start_idx + check_c;
            for (int check_component = 0; check_component < 9; ++check_component) {
                const double check_value =
                    all_photons->f[check_component * num_rays + check_i];
                if (!isfinite(check_value)) {
                    fprintf(
                        stderr,
                        "ERROR: numerical observer initialization produced a non-finite "
                        "state at ray %ld, component %d.\n",
                        check_i,
                        check_component);
                    exit(EXIT_FAILURE);
                }
            }
            if (!isfinite(all_photons->h[check_i])) {
                fprintf(
                    stderr,
                    "ERROR: numerical observer initialization produced a non-finite "
                    "step size at ray %ld.\n",
                    check_i);
                exit(EXIT_FAILURE);
            }
        }
"""
    if parallelization == "cuda":
        sync_and_transfer_code = r"""
        #ifdef DEBUG
        cudaDeviceSynchronize();
        cudaError_t err = cudaGetLastError();
        if (err != cudaSuccess) {
            fprintf(
                stderr,
                    "ERROR: numerical observer initialization kernel failed at start %ld: %s\n",
                start_idx,
                cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }
        #endif

        for (int component = 0; component < 9; ++component) {
            cudaMemcpy(
                all_photons->f + component * num_rays + start_idx,
                d_f_bundle + component * BUNDLE_CAPACITY,
                sizeof(double) * chunk_size,
                cudaMemcpyDeviceToHost);
        }
        cudaMemcpy(
            all_photons->h + start_idx,
            d_h_bundle,
            sizeof(double) * chunk_size,
            cudaMemcpyDeviceToHost);
"""
    else:
        sync_and_transfer_code = r"""
        for (int component = 0; component < 9; ++component) {
            memcpy(
                all_photons->f + component * num_rays + start_idx,
                d_f_bundle + component * BUNDLE_CAPACITY,
                sizeof(double) * chunk_size);
        }
        memcpy(
            all_photons->h + start_idx,
            d_h_bundle,
            sizeof(double) * chunk_size);
"""
    sync_and_transfer_code += finite_check

    host_loop_code = loop(
        idx_var="start_idx",
        lower_bound="0",
        upper_bound="num_rays",
        increment="BUNDLE_CAPACITY",
        pragma="",
        idx_type="long int",
        loop_body=f"""
        const long int chunk_size = NRPYMIN(num_rays - start_idx, BUNDLE_CAPACITY);
        {launch_code}
        {sync_and_transfer_code}
""",
    )

    # Step 8: Embed metric helpers.  The e_3 construction uses metric
    # Gram-Schmidt from a right-direction seed, so no inverse metric is needed
    # for this initializer.  The determinant still checks singular/non-Lorentzian
    # input before any tetrad normalization occurs.
    tetrad_helpers = r"""
static double nrpy_photon_metric_inner(
    const double metric[4][4], const double left[4], const double right[4])
{
    double result = 0.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            result += metric[mu][nu] * left[mu] * right[nu];
        }
    }
    return result;
}

static double nrpy_photon_metric_scale(const double metric[4][4])
{
    double scale = 1.0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            scale = fmax(scale, fabs(metric[mu][nu]));
        }
    }
    return scale;
}

static double nrpy_photon_metric_determinant(
    const double metric[4][4], const double pivot_tolerance)
{
    double work[4][4];
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = 0; nu < 4; ++nu) {
            work[mu][nu] = metric[mu][nu];
        }
    }

    double determinant = 1.0;
    for (int pivot_col = 0; pivot_col < 4; ++pivot_col) {
        int pivot_row = pivot_col;
        double pivot_abs = fabs(work[pivot_row][pivot_col]);
        for (int row = pivot_col + 1; row < 4; ++row) {
            const double candidate_abs = fabs(work[row][pivot_col]);
            if (candidate_abs > pivot_abs) {
                pivot_row = row;
                pivot_abs = candidate_abs;
            }
        }
        if (!isfinite(pivot_abs) || pivot_abs <= pivot_tolerance) {
            return NAN;
        }
        if (pivot_row != pivot_col) {
            for (int nu = 0; nu < 4; ++nu) {
                const double temporary = work[pivot_col][nu];
                work[pivot_col][nu] = work[pivot_row][nu];
                work[pivot_row][nu] = temporary;
            }
            determinant = -determinant;
        }
        determinant *= work[pivot_col][pivot_col];
        for (int row = pivot_col + 1; row < 4; ++row) {
            const double factor = work[row][pivot_col] / work[pivot_col][pivot_col];
            for (int nu = pivot_col + 1; nu < 4; ++nu) {
                work[row][nu] -= factor * work[pivot_col][nu];
            }
        }
    }
    return determinant;
}

static bool nrpy_photon_metric_orthogonalize_spacelike(
    const double metric[4][4],
    const double e0[4],
    const double e1[4],
    const double e2[4],
    const double seed[4],
    const double norm_tolerance,
    double output[4])
{
    // Project seed away from timelike e_0.  Since e_0.e_0 = -1, the projection
    // uses a plus sign: w = seed + (seed.e_0)e_0.
    const double seed_e0 = nrpy_photon_metric_inner(metric, seed, e0);
    for (int mu = 0; mu < 4; ++mu) {
        output[mu] = seed[mu] + seed_e0 * e0[mu];
    }

    // Project away from already-normalized spacelike e_1 and e_2.  Their norms
    // are +1, so both projections use minus signs.
    const double output_e1 = nrpy_photon_metric_inner(metric, output, e1);
    const double output_e2 = nrpy_photon_metric_inner(metric, output, e2);
    for (int mu = 0; mu < 4; ++mu) {
        output[mu] -= output_e1 * e1[mu] + output_e2 * e2[mu];
    }

    // A small or non-finite norm means seed is degenerate after metric
    // projection.  Caller may try a fallback seed.
    const double output_norm = nrpy_photon_metric_inner(metric, output, output);
    if (!isfinite(output_norm) || output_norm <= norm_tolerance) {
        return false;
    }
    const double inverse_norm = 1.0 / sqrt(output_norm);
    if (!isfinite(inverse_norm)) {
        return false;
    }
    for (int mu = 0; mu < 4; ++mu) {
        output[mu] *= inverse_norm;
    }
    return true;
}

static void nrpy_photon_cross_spatial(
    const double left[4], const double right[4], double output[4])
{
    output[0] = 0.0;
    output[1] = left[2] * right[3] - left[3] * right[2];
    output[2] = left[3] * right[1] - left[1] * right[3];
    output[3] = left[1] * right[2] - left[2] * right[1];
}

static void nrpy_photon_construct_observer_tetrad(
    const commondata_struct *commondata,
    const double metric_components[10],
    double observer_tetrad[4][4])
{
    // Step 1: Expand the ten independent symmetric g_mu_nu components into a
    // full four-dimensional covariant metric matrix.
    double metric[4][4] = {{0.0}};
    int component = 0;
    for (int mu = 0; mu < 4; ++mu) {
        for (int nu = mu; nu < 4; ++nu) {
            metric[mu][nu] = metric_components[component];
            metric[nu][mu] = metric_components[component];
            ++component;
        }
    }

    // Step 2: Reject non-finite or singular metrics before tetrad algebra.
    double metric_scale = nrpy_photon_metric_scale(metric);
    if (!isfinite(metric_scale) || metric_scale <= 0.0) {
        fprintf(stderr, "ERROR: observer metric has invalid scale.\n");
        exit(EXIT_FAILURE);
    }
    const double metric_entry_tolerance = 1.0e-14 * metric_scale;
    const double determinant_tolerance =
        1.0e-14 * metric_scale * metric_scale * metric_scale * metric_scale;
    const double determinant =
        nrpy_photon_metric_determinant(metric, metric_entry_tolerance);
    if (!isfinite(determinant) || determinant >= -determinant_tolerance) {
        fprintf(
            stderr,
            "ERROR: observer metric is singular or not Lorentzian (-+++); "
            "det(g)=% .17e.\n",
            determinant);
        exit(EXIT_FAILURE);
    }

    // Step 3: Build e_0 from stationary numerical-coordinate seed v_0=(1,0,0,0).
    const double v0[4] = {1.0, 0.0, 0.0, 0.0};
    const double v0_norm = nrpy_photon_metric_inner(metric, v0, v0);
    if (!isfinite(v0_norm) || v0_norm >= -metric_entry_tolerance) {
        fprintf(
            stderr,
            "ERROR: stationary observer seed is not timelike; "
            "g(v0,v0)=% .17e.\n",
            v0_norm);
        exit(EXIT_FAILURE);
    }
    const double inverse_e0_norm = 1.0 / sqrt(-v0_norm);
    if (!isfinite(inverse_e0_norm)) {
        fprintf(stderr, "ERROR: stationary observer normalization is non-finite.\n");
        exit(EXIT_FAILURE);
    }
    double e0[4];
    for (int mu = 0; mu < 4; ++mu) {
        e0[mu] = v0[mu] * inverse_e0_norm;
    }

    // Step 4: Build e_1 from the observer's supplied look-forward direction.
    // This is a direction vector, not a point and not a plane normal. It is
    // shared by every tile; no tile center enters the tetrad.
    const double v1[4] = {
        0.0,
        commondata->observer_look_forward_x,
        commondata->observer_look_forward_y,
        commondata->observer_look_forward_z,
    };
    const double forward_euclidean_norm =
        v1[1] * v1[1] + v1[2] * v1[2] + v1[3] * v1[3];
    if (!isfinite(forward_euclidean_norm) || forward_euclidean_norm <= 1.0e-28) {
        fprintf(
            stderr,
            "ERROR: observer look-forward direction is zero; cannot construct e_1.\n");
        exit(EXIT_FAILURE);
    }
    const double v1_e0 = nrpy_photon_metric_inner(metric, v1, e0);
    double e1_seed[4];
    for (int mu = 0; mu < 4; ++mu) {
        e1_seed[mu] = v1[mu] + v1_e0 * e0[mu];
    }
    const double e1_norm = nrpy_photon_metric_inner(metric, e1_seed, e1_seed);
    const double vector_norm_tolerance = 1.0e-14 * metric_scale;
    if (!isfinite(e1_norm) || e1_norm <= vector_norm_tolerance) {
        fprintf(
            stderr,
            "ERROR: forward observer seed has non-positive or tiny metric norm; "
            "g(w1,w1)=% .17e.\n",
            e1_norm);
        exit(EXIT_FAILURE);
    }
    const double inverse_e1_norm = 1.0 / sqrt(e1_norm);
    double e1[4];
    for (int mu = 0; mu < 4; ++mu) {
        e1[mu] = e1_seed[mu] * inverse_e1_norm;
    }

    // Step 5: Build e_2 from the supplied up seed using metric Gram-Schmidt.
    // If up is nearly parallel to e_1, try coordinate-axis seeds, still using
    // the same spacetime-metric projection and normalization.
    const double up_seed_candidates[4][4] = {
        {0.0, commondata->observer_up_x, commondata->observer_up_y, commondata->observer_up_z},
        {0.0, 1.0, 0.0, 0.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
    };
    double selected_up_seed[4] = {0.0, 0.0, 0.0, 0.0};
    double e2[4] = {0.0, 0.0, 0.0, 0.0};
    // Passing zero_spacelike_seed as the third basis vector disables the e_2
    // projection while reusing the same helper for the e_0/e_1 step.
    const double zero_spacelike_seed[4] = {0.0, 0.0, 0.0, 0.0};
    bool found_e2 = false;
    for (int candidate = 0; candidate < 4 && !found_e2; ++candidate) {
        found_e2 = nrpy_photon_metric_orthogonalize_spacelike(
            metric,
            e0,
            e1,
            zero_spacelike_seed,
            up_seed_candidates[candidate],
            vector_norm_tolerance,
            e2);
        if (found_e2) {
            for (int mu = 0; mu < 4; ++mu) {
                selected_up_seed[mu] = up_seed_candidates[candidate][mu];
            }
        }
    }
    if (!found_e2) {
        fprintf(stderr, "ERROR: no usable metric-orthogonal up seed for e_2.\n");
        exit(EXIT_FAILURE);
    }

    // Step 6: Construct e_3 by metric Gram-Schmidt from the Euclidean right
    // seed up x forward.  The cross product supplies orientation only; metric
    // orthogonalization supplies the actual spacetime normalization.
    double right_seed[4];
    nrpy_photon_cross_spatial(selected_up_seed, v1, right_seed);
    double e3[4];
    bool found_e3 = nrpy_photon_metric_orthogonalize_spacelike(
        metric,
        e0,
        e1,
        e2,
        right_seed,
        vector_norm_tolerance,
        e3);
    if (!found_e3) {
        fprintf(stderr, "ERROR: no usable metric-orthogonal right seed for e_3.\n");
        exit(EXIT_FAILURE);
    }

    // Preserve expected image handedness.  This dot product does not construct
    // e_3; it only detects whether metric Gram-Schmidt selected the opposite
    // sign from the supplied right-direction seed.
    const double right_orientation =
        e3[1] * right_seed[1] + e3[2] * right_seed[2] + e3[3] * right_seed[3];
    if (!isfinite(right_orientation)) {
        fprintf(stderr, "ERROR: right-direction orientation check is non-finite.\n");
        exit(EXIT_FAILURE);
    }
    if (right_orientation < 0.0) {
        for (int mu = 0; mu < 4; ++mu) {
            e3[mu] = -e3[mu];
        }
    }

    // Step 7: Validate all sixteen metric inner products against eta_ab.
    const double expected_signature[4] = {-1.0, 1.0, 1.0, 1.0};
    const double tetrad_tolerance = 1.0e-9;
    const double *tetrad[4] = {e0, e1, e2, e3};
    for (int a = 0; a < 4; ++a) {
        for (int b = 0; b < 4; ++b) {
            const double inner = nrpy_photon_metric_inner(metric, tetrad[a], tetrad[b]);
            const double expected = (a == b) ? expected_signature[a] : 0.0;
            if (!isfinite(inner) || fabs(inner - expected) > tetrad_tolerance) {
                fprintf(
                    stderr,
                    "ERROR: observer tetrad validation failed for (%d,%d): "
                    "inner=% .17e expected=% .17e.\n",
                    a,
                    b,
                    inner,
                    expected);
                exit(EXIT_FAILURE);
            }
        }
    }

    // Step 8: Return observer tetrad with documented indexing observer_tetrad[a][mu].
    for (int a = 0; a < 4; ++a) {
        for (int mu = 0; mu < 4; ++mu) {
            observer_tetrad[a][mu] = tetrad[a][mu];
        }
    }
}
"""

    # Step 9: Build the host orchestration around one observer metric and one
    # tetrad.  Plane-side flags remain initialized for current event detection;
    # they do not define ray direction.
    tetrad_scalar_declarations = "\n".join(
        f"    const double observer_tetrad_{basis}_{mu} = "
        f"observer_tetrad_out[{basis}][{mu}];"
        for basis in range(4)
        for mu in range(4)
    )
    body = r"""
    //==========================================
    // OBSERVER METRIC AND TETRAD SETUP
    //==========================================
    // observer_metric uses symmetric covariant ordering:
    // (g00,g01,g02,g03,g11,g12,g13,g22,g23,g33).
    if (observer_metric == NULL || observer_tetrad_out == NULL) {
        fprintf(stderr, "ERROR: numerical initializer received a null observer input.\n");
        exit(EXIT_FAILURE);
    }
    nrpy_photon_construct_observer_tetrad(
        commondata,
        observer_metric,
        observer_tetrad_out);

    // These scalar copies let each architecture pass the immutable per-tile
    // tetrad into its ray kernel without recomputing metric algebra per ray.
__TETRAD_SCALAR_DECLARATIONS__

    //==========================================
    // INPUT VALIDATION AND TILE SAMPLE GEOMETRY
    //==========================================
    const int tiles_width = commondata->tiles_width;
    const int tiles_height = commondata->tiles_height;
    const int tile_index_width = commondata->tile_index_width;
    const int tile_index_height = commondata->tile_index_height;
    const int scan_density_width = commondata->scan_density;
    if (tiles_width <= 0 || tiles_height <= 0 ||
        tile_index_width < 0 || tile_index_width >= tiles_width ||
        tile_index_height < 0 || tile_index_height >= tiles_height ||
        scan_density_width <= 0) {
        fprintf(
            stderr,
            "ERROR: invalid tile indices, tile counts, or scan_density.\n");
        exit(EXIT_FAILURE);
    }
    const double pi = acos(-1.0);
    if (!isfinite(commondata->alpha_w) || !isfinite(commondata->alpha_h) ||
        commondata->alpha_w <= 0.0 || commondata->alpha_h <= 0.0 ||
        commondata->alpha_w >= pi || commondata->alpha_h >= pi) {
        fprintf(stderr, "ERROR: alpha_w and alpha_h must lie strictly between 0 and pi.\n");
        exit(EXIT_FAILURE);
    }
    const double projected_tile_height_to_width =
        tan(0.5 * commondata->alpha_h) * (double)tiles_width /
        (tan(0.5 * commondata->alpha_w) * (double)tiles_height);
    const double scan_density_height_real =
        (double)scan_density_width * projected_tile_height_to_width;
    if (!isfinite(scan_density_height_real) ||
        scan_density_height_real < 1.0 ||
        scan_density_height_real > (double)INT_MAX) {
        fprintf(stderr, "ERROR: derived scan-density height is invalid.\n");
        exit(EXIT_FAILURE);
    }
    const int scan_density_height = (int)llround(scan_density_height_real);
    if (scan_density_height <= 0) {
        fprintf(stderr, "ERROR: derived scan-density height rounded to zero.\n");
        exit(EXIT_FAILURE);
    }
    if ((long int)scan_density_height >
        LONG_MAX / (long int)scan_density_width) {
        fprintf(stderr, "ERROR: tile ray count overflowed long int.\n");
        exit(EXIT_FAILURE);
    }
    const long int expected_ray_count =
        (long int)scan_density_width * (long int)scan_density_height;
    if (num_rays != expected_ray_count) {
        fprintf(
            stderr,
            "ERROR: initializer received %ld rays; expected %ld from scan density.\n",
            num_rays,
            expected_ray_count);
        exit(EXIT_FAILURE);
    }

    //==========================================
    // INITIAL EVENT-SIDE FLAGS
    //==========================================
    // Event detection uses each plane's own normal and center. The observer
    // direction and observer up seed are deliberately not used here.
    const double non_terminal_plane_normal_norm = sqrt(
        commondata->non_terminal_plane_normal_x * commondata->non_terminal_plane_normal_x +
        commondata->non_terminal_plane_normal_y * commondata->non_terminal_plane_normal_y +
        commondata->non_terminal_plane_normal_z * commondata->non_terminal_plane_normal_z);
    const double terminal_plane_normal_norm = sqrt(
        commondata->terminal_plane_normal_x * commondata->terminal_plane_normal_x +
        commondata->terminal_plane_normal_y * commondata->terminal_plane_normal_y +
        commondata->terminal_plane_normal_z * commondata->terminal_plane_normal_z);
    if (!isfinite(non_terminal_plane_normal_norm) ||
        !isfinite(terminal_plane_normal_norm) ||
        non_terminal_plane_normal_norm <= 1.0e-14 ||
        terminal_plane_normal_norm <= 1.0e-14) {
        fprintf(stderr, "ERROR: an event-plane normal is degenerate.\n");
        exit(EXIT_FAILURE);
    }
    if (!isfinite(commondata->terminal_plane_min_coord_radius) ||
        !isfinite(commondata->terminal_plane_max_coord_radius) ||
        commondata->terminal_plane_min_coord_radius < 0.0 ||
        commondata->terminal_plane_max_coord_radius <
            commondata->terminal_plane_min_coord_radius) {
        fprintf(
            stderr,
            "ERROR: terminal-plane coordinate-radius bounds are invalid.\n");
        exit(EXIT_FAILURE);
    }
    const double non_terminal_plane_side_value =
        commondata->non_terminal_plane_normal_x *
            (commondata->observer_x - commondata->non_terminal_plane_center_x) +
        commondata->non_terminal_plane_normal_y *
            (commondata->observer_y - commondata->non_terminal_plane_center_y) +
        commondata->non_terminal_plane_normal_z *
            (commondata->observer_z - commondata->non_terminal_plane_center_z);
    const bool init_non_terminal_plane_side = (non_terminal_plane_side_value > 0.0);
    const double terminal_plane_side_value =
        commondata->terminal_plane_normal_x *
            (commondata->observer_x - commondata->terminal_plane_center_x) +
        commondata->terminal_plane_normal_y *
            (commondata->observer_y - commondata->terminal_plane_center_y) +
        commondata->terminal_plane_normal_z *
            (commondata->observer_z - commondata->terminal_plane_center_z);
    const bool init_terminal_plane_side = (terminal_plane_side_value >= 0.0);
    if (all_photons->on_positive_side_of_non_terminal_plane_prev != NULL &&
        all_photons->on_positive_side_of_terminal_plane_prev != NULL) {
        for (long int plane_i = 0; plane_i < num_rays; ++plane_i) {
            all_photons->on_positive_side_of_non_terminal_plane_prev[plane_i] = init_non_terminal_plane_side;
            all_photons->on_positive_side_of_terminal_plane_prev[plane_i] = init_terminal_plane_side;
        }
    }

    //==========================================
    // STAGING ALLOCATION AND RAY INITIALIZATION
    //==========================================
    double *d_f_bundle;
    double *d_h_bundle;
    BHAH_MALLOC_DEVICE(d_f_bundle, sizeof(double) * 9 * BUNDLE_CAPACITY);
    BHAH_MALLOC_DEVICE(d_h_bundle, sizeof(double) * BUNDLE_CAPACITY);

    __HOST_LOOP_CODE__

    BHAH_FREE_DEVICE(d_f_bundle);
    BHAH_FREE_DEVICE(d_h_bundle);
__NORMALIZED_TRACKER_INITIALIZATION__
"""
    body = body.replace("__TETRAD_SCALAR_DECLARATIONS__", tetrad_scalar_declarations)
    body = body.replace("__HOST_LOOP_CODE__", str(host_loop_code))
    body = body.replace(
        "__NORMALIZED_TRACKER_INITIALIZATION__", normalized_tracker_initialization
    )
    if parallelization != "cuda":
        # The CPU backend has no device-allocation macro in generated headers;
        # use the ordinary host allocator for the temporary initializer bundle.
        body = body.replace("BHAH_MALLOC_DEVICE", "BHAH_MALLOC")
        body = body.replace("BHAH_FREE_DEVICE", "BHAH_FREE")
    prefunc = f"{tetrad_helpers}\n{kernel_prefunc}"

    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<limits.h>",
        "<stdbool.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]
    if parallelization == "cuda":
        includes.append("cuda_intrinsics.h")

    desc = r""" Initializes photons from a metric observer tetrad.

    The caller supplies one interpolated covariant observer metric at the common
    observer event.  This function constructs and validates one tetrad
    ``observer_tetrad[a][mu] = e_a^mu`` per tile-batch call.  The ray kernel
    performs only normalized sample-coordinate mapping, C/q evaluation, and
    the linear tetrad combination for the past-directed momentum.

    @param[in] commondata Observer, image, plane, and integration parameters.
    @param num_rays Number of rays in the current tile.
    @param[in,out] all_photons Host Structure of Arrays state and event flags.
    @param[in] observer_metric Ten independent covariant metric components at the
        observer event, ordered ``(g00,g01,g02,g03,g11,g12,g13,g22,g23,g33)``.
    @param[out] observer_tetrad_out Output tetrad with indexing ``[a][mu]``.
    @param stream_idx CUDA stream identifier when CUDA is enabled.
    """
    cfunc_type = "void"
    name = "set_initial_conditions_kernel"
    stream_param = "const int stream_idx" if parallelization == "cuda" else ""
    params = (
        "const commondata_struct *restrict commondata, "
        "long int num_rays, "
        "PhotonStateSoA *restrict all_photons, "
        "const double observer_metric[10], "
        "double observer_tetrad_out[4][4]"
    )
    if stream_param:
        params += f", {stream_param}"

    cfc.register_CFunction(
        prefunc=prefunc,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

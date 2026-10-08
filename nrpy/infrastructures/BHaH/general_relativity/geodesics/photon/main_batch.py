# nrpy/infrastructures/BHaH/general_relativity/geodesics/photon/main_batch.py
"""
Defines the main() C function for the geodesic integration pipeline.

This module acts as the master orchestrator, handling the explicit population of
RKF45 integrator configurations, simulation defaults, memory allocation for the
trajectory results, and the final serialization of the light blueprint.

The architecture implements pure angular-sample tiling. Observer parameters define
one common tetrad and field of view; each tile changes only its two tile indices.
The independent nonterminal and terminal planes are never shifted by tiling.
Normalized ray-sample coordinates are written during blueprint serialization;
plane coordinates remain event diagnostics.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import logging

import nrpy.c_function as cfc
import nrpy.params as par

# Configure logging to adhere to static analysis and project standards
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def main(
    spacetime_name: str,
    integrator_mode: str = "Analytical",
    normalized_eom: bool = False,
    axisymmetric_about_z: bool = False,
) -> None:
    """
    Register the master orchestrator C function for the geodesic integrator.

    This version implements Spatial Domain Decomposition (Tiling).

    :param spacetime_name: Metric name for numerical integration.
    :param integrator_mode: Batch-integrator implementation to emit.
    :param normalized_eom: Whether numerical evolution uses coordinate time as
        its integration parameter.
    :param axisymmetric_about_z: Whether numerical evolution includes the
        optional axial-angular-momentum diagnostic.
    :raises ValueError: If `integrator_mode` is not supported.
    :raises ValueError: If normalized evolution is requested for analytical mode.
    :raises ValueError: If the axial-angular-momentum diagnostic is requested
        for analytical mode.
    """
    if integrator_mode not in ("Analytical", "Numerical"):
        raise ValueError(
            "integrator_mode must be either 'Analytical' or 'Numerical'; "
            f"found '{integrator_mode}'."
        )
    if normalized_eom and integrator_mode != "Numerical":
        raise ValueError("normalized_eom is supported only for numerical mode.")
    if axisymmetric_about_z and integrator_mode != "Numerical":
        raise ValueError("axisymmetric_about_z is supported only for numerical mode.")

    diagnostic_flag_name = (
        "perform_normalization_check"
        if integrator_mode == "Numerical"
        else "perform_conservation_check"
    )
    diagnostic_flag_label = (
        "Normalization Check"
        if integrator_mode == "Numerical"
        else "Conservation Check"
    )
    if integrator_mode == "Numerical":
        Lz_path_argument = ", Lz_bin_name" if axisymmetric_about_z else ""
        batch_integrator_call = (
            "batch_integrator_numerical(&commondata, num_rays, results_buffer, "
            "norm_abs_bin_name, norm_abs_non_terminal_bin_name, "
            "synthetic_slice_usage_bin_name, "
            "non_terminal_crossings_path, terminal_crossings_path"
            f"{Lz_path_argument});"
        )
    else:
        batch_integrator_call = (
            "batch_integrator_analytical(&commondata, num_rays, results_buffer, "
            "non_terminal_crossings_path, terminal_crossings_path);"
        )
    analytical_plane_Lz_absence = (
        ""
        if integrator_mode == "Numerical"
        else """
            // Analytical batches do not evaluate axial angular momentum at
            // the first accepted state after the nonterminal plane.
            for (long int i = 0; i < num_rays; ++i) {
                results_buffer[i].non_terminal_post_step_Lz = NAN;
                results_buffer[i].non_terminal_post_step_t = NAN;
                results_buffer[i].non_terminal_post_step_distance = NAN;
            }
"""
    )
    numerical_startup_validation = (
        """
    if (time_window_manager_numerical_validate_batch_startup_from_bin(&commondata) !=
        TIME_WINDOW_MANAGER_NUMERICAL_SUCCESS) {
        return 1;
    } // END IF: numerical startup validation failed
"""
        if integrator_mode == "Numerical"
        else ""
    )
    serialization_desc = (
        " When normalization checks are enabled in numerical mode, each tile also "
        "writes a matching raw 'light_blueprint_norm_abs_XX_YY.bin' sidecar file "
        "and a sparse 'light_blueprint_norm_abs_non_terminal_XX_YY.bin' sidecar "
        "containing photon-index/norm-error records for nonterminal crossings. "
        "When synthetic slice-use checks are enabled, each tile writes "
        "'light_blueprint_synthetic_slice_usage_XX_YY.bin' with 26-byte records: "
        "a uint64_t photon index, two uint8_t use flags, and two f64 request times; "
        "false flags have NaN times."
        if integrator_mode == "Numerical"
        else ""
    )
    serialization_desc += (
        " Each tile also writes sparse 108-byte records to "
        "'light_blueprint_non_terminal_crossings_XX_YY.bin' and "
        "'light_blueprint_terminal_crossings_XX_YY.bin'."
    )
    if axisymmetric_about_z:
        serialization_desc += (
            " When the L_z check is enabled, each tile writes "
            "'light_blueprint_Lz_XX_YY.bin' records containing photon index, "
            "initial L_z, and final L_z."
        )
    numerical_sidecar_path_setup = (
        """
            char norm_abs_bin_name[256];
            if (commondata.perform_normalization_check) {
                snprintf(
                    norm_abs_bin_name,
                    sizeof(norm_abs_bin_name),
                    "light_blueprint_norm_abs_%02d_%02d.bin",
                    tx,
                    ty);
            } else {
                norm_abs_bin_name[0] = '\\0';
            } // END ELSE: disable normalization sidecar
            char norm_abs_non_terminal_bin_name[256];
            if (commondata.perform_normalization_check) {
                snprintf(
                    norm_abs_non_terminal_bin_name,
                    sizeof(norm_abs_non_terminal_bin_name),
                    "light_blueprint_norm_abs_non_terminal_%02d_%02d.bin",
                    tx,
                    ty);
            } else {
                norm_abs_non_terminal_bin_name[0] = '\\0';
            } // END ELSE: disable nonterminal normalization sidecar
            char synthetic_slice_usage_bin_name[256];
            if (commondata.perform_synthetic_slice_check) {
                snprintf(
                    synthetic_slice_usage_bin_name,
                    sizeof(synthetic_slice_usage_bin_name),
                    "light_blueprint_synthetic_slice_usage_%02d_%02d.bin",
                    tx,
                    ty);
            } else {
                synthetic_slice_usage_bin_name[0] = '\\0';
            } // END ELSE: disable synthetic slice-use sidecar
"""
        if integrator_mode == "Numerical"
        else ""
    )
    crossing_path_setup = r"""
            char non_terminal_crossings_path[256];
            char terminal_crossings_path[256];
            const int non_terminal_name_status = snprintf(
                non_terminal_crossings_path,
                sizeof(non_terminal_crossings_path),
                "light_blueprint_non_terminal_crossings_%02d_%02d.bin", tx, ty);
            const int terminal_name_status = snprintf(
                terminal_crossings_path,
                sizeof(terminal_crossings_path),
                "light_blueprint_terminal_crossings_%02d_%02d.bin", tx, ty);
            if (non_terminal_name_status < 0 ||
                (size_t)non_terminal_name_status >= sizeof(non_terminal_crossings_path) ||
                terminal_name_status < 0 ||
                (size_t)terminal_name_status >= sizeof(terminal_crossings_path)) {
                fprintf(stderr, "ERROR: plane crossing filename is too long.\n");
                free(results_buffer);
                return 1;
            }
"""
    if axisymmetric_about_z:
        numerical_sidecar_path_setup += r"""
            char Lz_bin_name[256];
            if (commondata.perform_Lz_check) {
                snprintf(
                    Lz_bin_name,
                    sizeof(Lz_bin_name),
                    "light_blueprint_Lz_%02d_%02d.bin",
                    tx,
                    ty);
            } else {
                Lz_bin_name[0] = '\0';
            } // END ELSE: disable axial-angular-momentum output
"""
    analytic_spacetime_telemetry = (
        """
    printf("--- Analytic Spacetime Physics ---\\n");
    printf("Mass Scale (M): %.2f\\n", commondata.M_scale);
    printf("Spin Parameter (a): %.2f\\n", commondata.a_spin);
"""
        if integrator_mode == "Analytical"
        else ""
    )
    numerical_spacetime_telemetry = (
        """
    printf("--- Numerical Spacetime ---\\n");
    printf("Data File: %s\\n", commondata.numerical_spacetime_bin_path);
    printf("Configured final numerical time: %.2f\\n", commondata.t_numerical_end);
    printf("Slice Spacing / Stride: %.6f / %d\\n", commondata.dt_numerical_spacetime_data, commondata.numerical_spacetime_time_slice_stride);
    printf("RKF45 Time-Window Lookahead: %.2f\\n", commondata.rkf45_max_delta_t);
"""
        if integrator_mode == "Numerical"
        else ""
    )
    eom_telemetry = (
        """
    printf("Equations of Motion: Normalized coordinate-time evolution\\n");
    printf("Log-Energy Tolerance: %e\\n", commondata.rkf45_log_energy_tolerance);
"""
        if normalized_eom
        else """
    printf("Equations of Motion: Affine-parameter geodesic evolution\\n");
"""
    )
    angular_momentum_telemetry = (
        '    printf("L_z Check: %d\\n", commondata.perform_Lz_check);\n'
        if axisymmetric_about_z
        else ""
    )
    synthetic_slice_usage_telemetry = (
        '    printf("Synthetic slice-use check: %d\\n", commondata.perform_synthetic_slice_check);\n'
        if integrator_mode == "Numerical"
        else ""
    )

    # Step 1: Register pure angular-grid tiling parameters. Physical plane
    # parameters belong to the initializer and event handlers, not to tiling.
    par.register_CodeParameter(
        "int", __name__, "tiles_width", 1, commondata=True, add_to_parfile=True
    )
    par.register_CodeParameter(
        "int",
        __name__,
        "tiles_height",
        1,
        commondata=True,
        add_to_parfile=True,
    )
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "limits.h",
        "stdio.h",
        "stdlib.h",
        "math.h",
    ]

    desc = f""" Master entry point for the integrator.

    Algorithm:
    1. Initializes core data structures and sets default physical constants.
    2. Parses command-line arguments and parameter files to override defaults.
    3. Performs numerical startup validation once when numerical mode is selected.
    4. Validates and normalizes the independent event-plane normals.
    5. Loops through a grid of tiles (tx, ty), changing only active tile indices.
    6. Dispatches the selected {integrator_mode.lower()} batch integrator for the {spacetime_name} metric per tile.
    7. Serializes each tile to a versioned native 'light_blueprint_XX_YY.bin'.{serialization_desc}"""

    cfunc_type = "int"
    name = "main"
    params = "int argc, const char *argv[]"

    # Step 2: Build the C body
    body = f"""
    //==========================================
    // CORE STRUCTURE INITIALIZATION
    //==========================================
    commondata_struct commondata;

    //==========================================
    // PARAMETER PARSING
    //==========================================
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);

{numerical_startup_validation}

    // Validate each enabled plane after reading the parameter file.
    const bool plane_enabled[2] = {{
        commondata.non_terminal_plane_enabled,
        commondata.terminal_plane_enabled
    }};
    const char *plane_name[2] = {{"nonterminal", "terminal"}};
    const double plane_normal[2][3] = {{
        {{commondata.non_terminal_plane_normal_x,
          commondata.non_terminal_plane_normal_y,
          commondata.non_terminal_plane_normal_z}},
        {{commondata.terminal_plane_normal_x,
          commondata.terminal_plane_normal_y,
          commondata.terminal_plane_normal_z}}
    }};
    const double plane_up[2][3] = {{
        {{commondata.non_terminal_plane_up_x,
          commondata.non_terminal_plane_up_y,
          commondata.non_terminal_plane_up_z}},
        {{commondata.terminal_plane_up_x,
          commondata.terminal_plane_up_y,
          commondata.terminal_plane_up_z}}
    }};
    for (int plane = 0; plane < 2; ++plane) {{
        if (!plane_enabled[plane]) continue;
        const double nx = plane_normal[plane][0];
        const double ny = plane_normal[plane][1];
        const double nz = plane_normal[plane][2];
        const double ux = plane_up[plane][0];
        const double uy = plane_up[plane][1];
        const double uz = plane_up[plane][2];
        const double normal_mag = sqrt(nx * nx + ny * ny + nz * nz);
        const double up_mag = sqrt(ux * ux + uy * uy + uz * uz);
        if (!isfinite(normal_mag) || !isfinite(up_mag) ||
            normal_mag < 1.0e-2 || up_mag < 1.0e-2) {{
            fprintf(stderr,
                "FATAL: %s plane requires finite normal and up magnitudes >= 1e-2; normal=(%.17g, %.17g, %.17g), up=(%.17g, %.17g, %.17g).\\n",
                plane_name[plane], nx, ny, nz, ux, uy, uz);
            return 1;
        }} // END IF: invalid plane vector magnitude
        const double alignment = fabs((nx * ux + ny * uy + nz * uz) /
                                      (normal_mag * up_mag));
        if (!isfinite(alignment) || alignment > 0.25) {{
            fprintf(stderr,
                "FATAL: %s plane requires |normal dot up|/(|normal||up|) <= 0.25; value=%.17g, normal=(%.17g, %.17g, %.17g), up=(%.17g, %.17g, %.17g).\\n",
                plane_name[plane], alignment, nx, ny, nz, ux, uy, uz);
            return 1;
        }} // END IF: plane vectors insufficiently perpendicular
        if (plane == 0) {{
            commondata.non_terminal_plane_normal_x = nx / normal_mag;
            commondata.non_terminal_plane_normal_y = ny / normal_mag;
            commondata.non_terminal_plane_normal_z = nz / normal_mag;
        }} // END IF: normalize nonterminal plane normal
        else {{
            commondata.terminal_plane_normal_x = nx / normal_mag;
            commondata.terminal_plane_normal_y = ny / normal_mag;
            commondata.terminal_plane_normal_z = nz / normal_mag;
        }} // END ELSE: normalize terminal plane normal
    }} // END LOOP: enabled plane validation

    //==========================================
    // TELEMETRY
    //==========================================
    printf("=============================================\\n");
    printf("  Geodesic Engine\\n");
    printf("=============================================\\n");
    printf("Spacetime Metric: {spacetime_name}\\n");
{analytic_spacetime_telemetry}
{numerical_spacetime_telemetry}

    printf("--- Observer ---\\n");
    printf("Position (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.observer_x, commondata.observer_y, commondata.observer_z);
    printf("Look-forward (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.observer_look_forward_x,
           commondata.observer_look_forward_y,
           commondata.observer_look_forward_z);
    printf("Up seed (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.observer_up_x, commondata.observer_up_y, commondata.observer_up_z);
    printf("FOV (width x height): %.6f x %.6f radians\\n",
           commondata.alpha_w, commondata.alpha_h);
    printf("Tile grid: %d x %d tiles, scan density width=%d\\n",
           commondata.tiles_width, commondata.tiles_height,
           commondata.scan_density);

    printf("--- Nonterminal Plane ---\\n");
    printf("Center (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.non_terminal_plane_center_x,
           commondata.non_terminal_plane_center_y,
           commondata.non_terminal_plane_center_z);
    printf("Normal (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.non_terminal_plane_normal_x,
           commondata.non_terminal_plane_normal_y,
           commondata.non_terminal_plane_normal_z);
    printf("Up seed (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.non_terminal_plane_up_x,
           commondata.non_terminal_plane_up_y,
           commondata.non_terminal_plane_up_z);

    printf("--- Terminal Plane ---\\n");
    printf("Center (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.terminal_plane_center_x,
           commondata.terminal_plane_center_y,
           commondata.terminal_plane_center_z);
    printf("Normal (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.terminal_plane_normal_x,
           commondata.terminal_plane_normal_y,
           commondata.terminal_plane_normal_z);
    printf("Up seed (x, y, z): %.2f, %.2f, %.2f\\n",
           commondata.terminal_plane_up_x,
           commondata.terminal_plane_up_y,
           commondata.terminal_plane_up_z);
    printf("Coordinate-radius bounds: %.2f -> %.2f\\n",
           commondata.terminal_plane_min_coord_radius,
           commondata.terminal_plane_max_coord_radius);

    printf("--- Temporal & Boundary Conditions ---\\n");
    printf("Start Time (t_start): %.2f\\n", commondata.t_start);
    printf("Escape Radius (r_escape): %.2f\\n", commondata.r_escape);
    printf("Evolution Measure Limit: %.2f\\n", commondata.evolution_measure_max);

    printf("--- Batch Integrator & RKF45 Settings ---\\n");
    printf("Initial Step Size: %.2f\\n", commondata.initial_h);
    printf("Min / Max Step: %e / %.2f\\n",
           commondata.rkf45_h_min, commondata.rkf45_h_max);
    printf("Abs / Rel Tolerance: %e / %e\\n",
           commondata.rkf45_absolute_error_tolerance,
           commondata.rkf45_error_tolerance);
{eom_telemetry}
    printf("Adaptive Step Damping: 0.90\\n");
    printf("Max Retries: %d\\n", commondata.rkf45_max_retries);
    printf("{diagnostic_flag_label}: %d\\n", commondata.{diagnostic_flag_name});
{synthetic_slice_usage_telemetry}
{angular_momentum_telemetry}
    printf("Slot Manager (Delta t / Min t): %.2f / %.2f\\n",
           commondata.slot_manager_delta_t, commondata.slot_manager_t_min);

    //==========================================
    // ANGULAR SAMPLE-GRID VALIDATION
    //==========================================
    const int scan_density_width = commondata.scan_density;
    const double pi = acos(-1.0);
    if (commondata.tiles_width <= 0 || commondata.tiles_height <= 0 ||
        scan_density_width <= 0 || !isfinite(commondata.alpha_w) ||
        !isfinite(commondata.alpha_h) || commondata.alpha_w <= 0.0 ||
        commondata.alpha_h <= 0.0 || commondata.alpha_w >= pi ||
        commondata.alpha_h >= pi) {{
        fprintf(stderr,
                "FATAL: tile counts, scan_density, and fields of view must be valid.\\n");
        return 1;
    }} // END IF: invalid angular sampling
    const double projected_tile_height_to_width =
        tan(0.5 * commondata.alpha_h) * (double)commondata.tiles_width /
        (tan(0.5 * commondata.alpha_w) * (double)commondata.tiles_height);
    const double scan_density_height_real =
        (double)scan_density_width * projected_tile_height_to_width;
    if (!isfinite(scan_density_height_real) ||
        scan_density_height_real < 1.0 ||
        scan_density_height_real > (double)INT_MAX) {{
        fprintf(stderr, "FATAL: derived scan-density height is invalid.\\n");
        return 1;
    }} // END IF: invalid height sample count
    const int scan_density_height = (int)llround(scan_density_height_real);
    if (scan_density_height <= 0 ||
        (long int)scan_density_height >
            LONG_MAX / (long int)scan_density_width) {{
        fprintf(stderr, "FATAL: derived tile ray count is invalid.\\n");
        return 1;
    }} // END IF: invalid tile ray count
    const uint64_t record_count =
        (uint64_t)scan_density_width * (uint64_t)scan_density_height;
    if (record_count > (uint64_t)LONG_MAX) {{
        fprintf(stderr, "FATAL: ray count does not fit in a host long int.\\n");
        return 1;
    }} // END IF: ray count overflow
    const long int num_rays = (long int)record_count;
    blueprint_data_t *restrict results_buffer =
        (blueprint_data_t *)malloc(sizeof(blueprint_data_t) * num_rays);
    if (results_buffer == NULL) {{
        fprintf(stderr, "FATAL: failed to allocate %ld rays.\\n", num_rays);
        return 1;
    }} // END IF: result allocation failed

    //==========================================
    // ANGULAR SAMPLE TILING LOOP
    //==========================================
    for (int ty = 0; ty < commondata.tiles_height; ++ty) {{
        for (int tx = 0; tx < commondata.tiles_width; ++tx) {{
            // Only active tile indices change. Observer setup, fields of view,
            // tetrad inputs, and both physical event planes stay fixed.
            commondata.tile_index_width = tx;
            commondata.tile_index_height = ty;

            printf("[Tile %02d,%02d] samples=%ld\\n", tx, ty, num_rays);

{numerical_sidecar_path_setup}
{crossing_path_setup}
            for (long int i = 0; i < num_rays; ++i) {{
                results_buffer[i] = (blueprint_data_t){{0}};
                results_buffer[i].termination_type = FAILURE_GENERIC;
                results_buffer[i].non_terminal_post_step_Lz = NAN;
                results_buffer[i].non_terminal_post_step_t = NAN;
                results_buffer[i].non_terminal_post_step_distance = NAN;
            }} // END LOOP: initialize tile results

            {batch_integrator_call}
{analytical_plane_Lz_absence}

            // Store normalized ray-sample coordinates before serialization.
            // They are derived from stable ray ordinal and tile indices, so
            // integration does not carry image metadata in PhotonStateSoA.
            for (long int i = 0; i < num_rays; ++i) {{
                const int local_col = (int)(i % scan_density_width);
                const int local_row = (int)(i / scan_density_width);
                results_buffer[i].image_width_fraction =
                    ((double)tx +
                     ((double)local_col + 0.5) / (double)scan_density_width) /
                    (double)commondata.tiles_width;
                results_buffer[i].image_height_fraction =
                    ((double)ty +
                     ((double)local_row + 0.5) / (double)scan_density_height) /
                    (double)commondata.tiles_height;
            }} // END LOOP: serialize ray coordinates

            char bin_name[256];
            const int name_status = snprintf(
                bin_name, sizeof(bin_name),
                "light_blueprint_%02d_%02d.bin", tx, ty);
            if (name_status < 0 || (size_t)name_status >= sizeof(bin_name)) {{
                fprintf(stderr, "ERROR: blueprint filename is too long.\\n");
                free(results_buffer);
                return 1;
            }} // END IF: blueprint filename overflow

            blueprint_header_t header = {{0}};
            memcpy(header.magic, BLUEPRINT_MAGIC, sizeof(header.magic));
            header.native_schema_version = BLUEPRINT_SCHEMA_VERSION;
            header.header_size = (uint32_t)sizeof(blueprint_header_t);
            header.record_size = (uint32_t)sizeof(blueprint_data_t);
            header.tx = (uint32_t)tx;
            header.ty = (uint32_t)ty;
            header.tiles_w = (uint32_t)commondata.tiles_width;
            header.tiles_h = (uint32_t)commondata.tiles_height;
            header.record_count = record_count;
            header.alpha_w = commondata.alpha_w;
            header.alpha_h = commondata.alpha_h;

            if (write_blueprint_file(bin_name, &header, results_buffer, record_count) != 0) {{
                free(results_buffer);
                return 1;
            }} // END IF: blueprint write failed
        }} // END LOOP: tile columns
    }} // END LOOP: tile rows

    free(results_buffer);
    printf("Simulation successful. Data stored in native blueprint tiles.\\n");
    return 0;
    """

    prefunc = r"""
/**
 * Write one tile's blueprint header and photon records.
 *
 * @param[in] filename Output binary-file path.
 * @param[in] header Header with tile dimensions and ray count.
 * @param[in] records Photon records in tile sample order.
 * @param record_count Number of photon records to write.
 * @return 0 on success, 1 on file failure.
 */
static int write_blueprint_file(
    const char *filename,
    const blueprint_header_t *header,
    const blueprint_data_t *records,
    uint64_t record_count) {
    FILE *blueprint_file = fopen(filename, "wb");
    if (blueprint_file == NULL) {
        fprintf(stderr, "ERROR: Could not open '%s' for writing.\\n", filename);
        return 1;
    } // END IF: blueprint open failed
    if (fwrite(header, sizeof(*header), 1, blueprint_file) != 1) {
        fprintf(stderr, "ERROR: Could not write header to '%s'.\\n", filename);
        fclose(blueprint_file);
        return 1;
    } // END IF: blueprint header write failed
    if (fwrite(records, sizeof(*records), (size_t)record_count, blueprint_file) !=
        (size_t)record_count) {
        fprintf(stderr, "ERROR: Could not write payload to '%s'.\\n", filename);
        fclose(blueprint_file);
        return 1;
    } // END IF: photon record write failed
    if (fclose(blueprint_file) != 0) {
        fprintf(stderr, "ERROR: Could not close '%s'.\\n", filename);
        return 1;
    } // END IF: blueprint close failed
    return 0;
} // END FUNCTION: write_blueprint_file
"""

    # Step 3: Register the function with the NRPy environment
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
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

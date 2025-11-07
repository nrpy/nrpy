"""
Generate the main() C function, the entry point for the executable.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def main_c_function() -> None:
    """
    Generate the main() C function, the entry point for the executable.

    This function acts as the master orchestrator. It initializes parameters,
    determines the run mode (debug or production), and dispatches to the
    appropriate top-level orchestrator.
    """
    name = "main"
    cfunc_type = "int"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<string.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]
    desc = "Main entry point for the massive particle geodesic integrator."
    params = "int argc, const char *argv[]"

    debug_body = r"""
        #include "set_CodeParameters-nopointer.h"

        particle_initial_state_t initial_state;
        initial_state.id = 0;
        const char *filename = "particle_debug_initial_conditions.txt";
        FILE *fp_in = fopen(filename, "r");
        if (fp_in == NULL) {
            printf("File '%s' not found. Creating it with default values.n", filename);
            fp_in = fopen(filename, "w");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not create '%s'.n", filename); return 1; }
            fprintf(fp_in, "# Format: t_initial pos_x pos_y pos_z u_x u_y u_z\n");
            fprintf(fp_in, "0.0 10.0  0.0  0.0   0.0  0.363232  0.0\n");
            fclose(fp_in);
            fp_in = fopen(filename, "r");
            if (fp_in == NULL) { fprintf(stderr, "Error: Could not re-open '%s'.n", filename); return 1; }
        }
        char line[256];
        while (fgets(line, sizeof(line), fp_in)) {
            if (line[0] != '#') {
                sscanf(line, "%lf %lf %lf %lf %lf %lf %lf",
                       &initial_state.pos[0], &initial_state.pos[1], &initial_state.pos[2], &initial_state.pos[3],
                       &initial_state.u_spatial[0], &initial_state.u_spatial[1], &initial_state.u_spatial[2]);
                break;
            }
        }
        fclose(fp_in);
        printf("--- Single Particle Debug Run ---\n");
        printf("  pos = (t=%.4f, x=%.4f, y=%.4f, z=%.4f)\n", initial_state.pos[0], initial_state.pos[1], initial_state.pos[2], initial_state.pos[3]);
        double y_start[8], y_final[8];
        set_initial_conditions_massive(&initial_state, &commondata, &params, y_start);
        printf("\nInitial State Vector (y_start):\n");
        printf("  t=%.2f, x=%.2f, y=%.2f, z=%.2f\n", y_start[0], y_start[1], y_start[2], y_start[3]);
        printf("  u^t=%.4f, u^x=%.4f, u^y=%.4f, u^z=%.4f\n\n", y_start[4], y_start[5], y_start[6], y_start[7]);
        if(commondata.perform_conservation_check) {
            double E_i, Lx_i, Ly_i, Lz_i, Q_i, E_f, Lx_f, Ly_f, Lz_f, Q_f;
            check_conservation_massive(&commondata, &params, &metric, y_start, &E_i, &Lx_i, &Ly_i, &Lz_i, &Q_i);
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
            check_conservation_massive(&commondata, &params, &metric, y_final, &E_f, &Lx_f, &Ly_f, &Lz_f, &Q_f);
        } else {
            integrate_single_particle_DEBUG(&commondata, &params, &metric, y_start, y_final);
        }
        printf("\nDebug run finished. Trajectory saved to 'massive_particle_path.txt'.\n");
    """

    body = f"""
    commondata_struct commondata;
    params_struct params;
    metric_params metric;
    commondata_struct_set_to_default(&commondata);
    cmdline_input_and_parfile_parser(&commondata, argc, argv);
    metric.type = (commondata.a_spin == 0.0) ? Schwarzschild : Kerr;
    if (commondata.run_in_debug_mode) {{
        {debug_body}
    }} else {{
        printf("----------------------------------------\\n");
        printf("Massive Particle Integrator\\n");
        printf("----------------------------------------\\n");
        printf("Metric Settings: %s (a=%.2f, M=%.2f)\\n", (metric.type == Kerr) ? "Kerr" : "Schwarzschild", commondata.a_spin, commondata.M_scale);
        printf("Integration Settings: Max Time = %.1f M, Snapshot Every = %.1f M\\n", commondata.t_final, commondata.snapshot_every_t);
        printf("Initial Conditions: Type = %s, Particles = %d x %d, Radius = %.2f to %.2f M\\n", commondata.initial_conditions_type, commondata.disk_num_r, commondata.disk_num_phi, commondata.disk_r_min, commondata.disk_r_max);
        printf("----------------------------------------\\n\\n");
        run_mass_integrator_production(&commondata, &params, &metric);
        printf("\\nProduction run finished successfully.\\n");
        // NOTE: The post-processing system() call has been REMOVED from the C code.
        // The Python script that calls this executable is now responsible for post-processing.
    }}
    return 0;
    """

    cfc.register_CFunction(
        name=name,
        cfunc_type=cfunc_type,
        includes=includes,
        desc=desc,
        params=params,
        body=body
    )
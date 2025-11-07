"""
Generate the C function for a full production run of the mass integrator.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def run_mass_integrator_production() -> None:
    """
    Generate the top-level orchestrator for a full production run.

    This C function sets up the particle ensemble, then enters the main time
    loop. At each step, it saves a snapshot and then calls
    `integrate_single_particle` for each active particle in a parallel loop.
    """
    name = "run_mass_integrator_production"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<sys/stat.h>",
        "<string.h>",
        "<stdio.h>",
        "<stdlib.h>",
    ]
    desc = "Orchestrates the full production run for the mass integrator."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  const metric_params *restrict metric"""
    body = r"""
    const int max_particles = commondata->disk_num_r * commondata->disk_num_phi;
    double *y_initial_states = (double *)malloc(sizeof(double) * 8 * max_particles);
    if (!y_initial_states) { fprintf(stderr, "Error: Failed to allocate memory for initial states.n"); exit(1); }
    int num_particles = 0;
    printf("Generating initial conditions of type: %s\n", commondata->initial_conditions_type);
    if (strcmp(commondata->initial_conditions_type, "SpiralGalaxy") == 0) {
        num_particles = generate_spiral_galaxy_initial_conditions(commondata, params, y_initial_states);
    } else if (strcmp(commondata->initial_conditions_type, "BarredFlocculentSpiral") == 0) {
        num_particles = generate_barred_flocculent_spiral_ic(commondata, params, y_initial_states);
    } else {
        if (strcmp(commondata->initial_conditions_type, "KeplerianDisk") != 0) {
            printf("Warning: Unrecognized type '%s'. Defaulting to KeplerianDisk.n", commondata->initial_conditions_type);
        }
        num_particles = generate_disk_initial_conditions(commondata, params, y_initial_states);
    }
    printf("Generated %d particles.\n", num_particles);
    mass_particle_state_t *particle_states = (mass_particle_state_t *)malloc(sizeof(mass_particle_state_t) * num_particles);
    if (!particle_states) { fprintf(stderr, "Error: Failed to allocate memory for particle states.n"); exit(1); }
    for (int i=0; i<num_particles; i++) {
        double *y_start = &y_initial_states[i*8];
        particle_states[i].id = i;
        particle_states[i].pos[0] = y_start[1]; particle_states[i].pos[1] = y_start[2]; particle_states[i].pos[2] = y_start[3];
        particle_states[i].u[0] = y_start[4]; particle_states[i].u[1] = y_start[5]; particle_states[i].u[2] = y_start[6]; particle_states[i].u[3] = y_start[7];
        const double r = sqrt(y_start[1]*y_start[1] + y_start[2]*y_start[2]);
        particle_states[i].lambda_rest = commondata->disk_lambda_rest_at_r_min * pow(r / commondata->disk_r_min, 0.75);
        particle_states[i].j_intrinsic = (float)pow(r / commondata->disk_r_min, -3.0);
    }
    free(y_initial_states);
    mkdir(commondata->output_folder, 0755);
    int snapshot_count = 0;
    for (double current_t = 0; current_t <= commondata->t_final; current_t += commondata->snapshot_every_t) {
        char filename[200];
        snprintf(filename, 200, "%s/mass_blueprint_t_%04d.bin", commondata->output_folder, snapshot_count);
        printf("Saving snapshot: %s (t=%.2f)\n", filename, current_t);
        FILE *fp_out = fopen(filename, "wb");
        if (!fp_out) { exit(1); }
        int active_particles = 0;
        for(int i=0; i<num_particles; i++) if (!isnan(particle_states[i].pos[0])) active_particles++;
        fwrite(&active_particles, sizeof(int), 1, fp_out);
        for (int i=0; i<num_particles; i++) if (!isnan(particle_states[i].pos[0])) fwrite(&particle_states[i], sizeof(mass_particle_state_t), 1, fp_out);
        fclose(fp_out);
        snapshot_count++;
        if (current_t >= commondata->t_final) break;
        const double t_next_snapshot = current_t + commondata->snapshot_every_t;
        #pragma omp parallel for
        for (int i = 0; i < num_particles; i++) {
            if (isnan(particle_states[i].pos[0])) continue;
            double y_particle[8];
            y_particle[0] = current_t;
            y_particle[1] = particle_states[i].pos[0]; y_particle[2] = particle_states[i].pos[1]; y_particle[3] = particle_states[i].pos[2];
            y_particle[4] = particle_states[i].u[0]; y_particle[5] = particle_states[i].u[1]; y_particle[6] = particle_states[i].u[2]; y_particle[7] = particle_states[i].u[3];
            int status = integrate_single_particle(commondata, params, metric, y_particle[0], t_next_snapshot, y_particle);
            if (status != 0) {
                particle_states[i].pos[0] = NAN;
            } else {
                particle_states[i].pos[0] = y_particle[1]; particle_states[i].pos[1] = y_particle[2]; particle_states[i].pos[2] = y_particle[3];
                particle_states[i].u[0] = y_particle[4]; particle_states[i].u[1] = y_particle[5]; particle_states[i].u[2] = y_particle[6]; particle_states[i].u[3] = y_particle[7];
            }
        }
    }
    free(particle_states);
    """

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True
    )
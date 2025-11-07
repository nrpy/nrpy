"""
Generate the C function for a barred, clumpy ("flocculent") spiral galaxy.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def generate_barred_flocculent_spiral_ic() -> None:
    """
    Generate the C orchestrator for a barred, clumpy ("flocculent") spiral.

    This function uses rejection sampling to create a complex geometry with a
    dense central bulge, a rectangular bar, and clumpy spiral arms.
    """
    name = "generate_barred_flocculent_spiral_ic"
    cfunc_type = "int"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdlib.h>",
        "<time.h>",
    ]
    desc = "Generates the initial state for all particles in a barred flocculent spiral galaxy."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    unsigned int seed = 42;
    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;
    const double r_min = commondata->disk_r_min, r_max = commondata->disk_r_max;
    const double bar_len = commondata->bar_length, bar_width = commondata->bar_length * commondata->bar_aspect_ratio;
    const double bulge_rad = commondata->bulge_radius;
    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;
    for (int i = 0; i < num_particles_total; i++) {
        double r, phi, x, y;
        seed += i;
        while (1) {
            double random_val = (double)rand_r(&seed) / RAND_MAX;
            r = sqrt(random_val * (r_max*r_max - r_min*r_min) + r_min*r_min);
            phi = 2.0 * M_PI * ((double)rand_r(&seed) / RAND_MAX);
            x = r * cos(phi); y = r * sin(phi);
            double acceptance_prob = 0.0;
            if (r < bulge_rad) {
                acceptance_prob = commondata->bulge_density_factor * commondata->arm_particle_density;
            } else if (fabs(x) < bar_len / 2.0 && fabs(y) < bar_width / 2.0) {
                acceptance_prob = commondata->bar_density_factor * commondata->arm_particle_density;
            } else {
                if (((double)rand_r(&seed) / RAND_MAX) < commondata->arm_particle_density) {
                    double theta_base = (1.0 / arm_tightness) * log(r / r_min);
                    double delta_phi_raw = phi - theta_base;
                    double angle_between_arms = 2.0 * M_PI / num_arms;
                    double delta_phi_folded = fmod(delta_phi_raw, angle_between_arms);
                    if (delta_phi_folded > 0.5 * angle_between_arms) delta_phi_folded -= angle_between_arms;
                    else if (delta_phi_folded < -0.5 * angle_between_arms) delta_phi_folded += angle_between_arms;
                    double arm_profile = exp(-0.5 * SQR(delta_phi_folded / (arm_tightness * 0.5)));
                    double clump_modulation = cos(commondata->arm_clumpiness_factor * (theta_base - phi));
                    clump_modulation *= cos(num_arms * phi * commondata->arm_clump_size);
                    acceptance_prob = arm_profile + 0.5 * clump_modulation;
                }
            }
            if (((double)rand_r(&seed) / RAND_MAX) < acceptance_prob) break;
        }
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        double *y_state = &y_initial_states[particle_count * 8];
        y_state[0] = 0.0; y_state[1] = x; y_state[2] = y; y_state[3] = 0.0;
        y_state[4] = ut_at_r; y_state[5] = -y * uphi_at_r; y_state[6] = x * uphi_at_r; y_state[7] = 0.0;
        particle_count++;
    }
    return particle_count;
    """

    cfc.register_CFunction(
        name=name,
        cfunc_type=cfunc_type,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
        include_CodeParameters_h=True
    )
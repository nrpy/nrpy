"""
Generate the C function to create a "grand design" spiral galaxy.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def generate_spiral_galaxy_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a "grand design" spiral galaxy.

    This function uses a random placement strategy based on a logarithmic
    spiral formula to create a spiral galaxy pattern.
    """
    name = "generate_spiral_galaxy_initial_conditions"
    cfunc_type = "int"
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "<math.h>",
        "<stdlib.h>",
    ]
    desc = "Generates the complete y[8] initial state for all particles in a spiral galaxy disk."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    int particle_count = 0;
    const int num_particles_total = commondata->disk_num_r * commondata->disk_num_phi;
    const int num_arms = commondata->spiral_galaxy_num_arms;
    const double arm_tightness = commondata->spiral_galaxy_arm_tightness;
    srand(42);
    for (int i = 0; i < num_particles_total; i++) {
        const double r = commondata->disk_r_min + (commondata->disk_r_max - commondata->disk_r_min) * sqrt((double)rand() / RAND_MAX);
        const double theta_base = (1.0 / arm_tightness) * log(r / commondata->disk_r_min);
        const int arm_index = rand() % num_arms;
        const double arm_offset = (2.0 * M_PI / num_arms) * arm_index;
        const double phi_spread = (M_PI / num_arms) * 0.2 * ((double)rand() / RAND_MAX - 0.5);
        const double phi = theta_base + arm_offset + phi_spread;
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        double *y = &y_initial_states[particle_count * 8];
        y[0] = 0.0; y[1] = r * cos(phi); y[2] = r * sin(phi); y[3] = 0.0;
        y[4] = ut_at_r; y[5] = -y[2] * uphi_at_r; y[6] = y[1] * uphi_at_r; y[7] = 0.0;
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
"""
Generate the C function to create a simple, axisymmetric Keplerian disk.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def generate_disk_initial_conditions() -> None:
    """
    Generate a C orchestrator to create a simple, axisymmetric Keplerian disk.

    This C function populates a large array of initial state vectors by looping
    through a grid of radii `r` and azimuthal angles `phi`.
    """
    name = "generate_disk_initial_conditions"
    cfunc_type = "int"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Generates the complete y[8] initial state for all particles in a Keplerian disk."
    params = """const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double *restrict y_initial_states"""
    body = r"""
    int particle_count = 0;
    const double dr = (commondata->disk_num_r > 1) ? (commondata->disk_r_max - commondata->disk_r_min) / (commondata->disk_num_r - 1) : 0;
    for (int i = 0; i < commondata->disk_num_r; i++) {
        const double r = commondata->disk_r_min + i * dr;
        const int num_phi_at_r = (commondata->disk_num_phi > 1) ? (int)(commondata->disk_num_phi * (r / commondata->disk_r_max)) : 1;
        if (num_phi_at_r == 0) continue;
        const double dphi = 2.0 * M_PI / num_phi_at_r;
        double ut_at_r, uphi_at_r;
        calculate_ut_uphi_from_r(r, commondata, params, &ut_at_r, &uphi_at_r);
        for (int j = 0; j < num_phi_at_r; j++) {
            const double phi = j * dphi;
            double *y = &y_initial_states[particle_count * 8];
            y[0] = 0.0; y[1] = r * cos(phi); y[2] = r * sin(phi); y[3] = 0.0;
            y[4] = ut_at_r; y[5] = -y[2] * uphi_at_r; y[6] = y[1] * uphi_at_r; y[7] = 0.0;
            particle_count++;
        }
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
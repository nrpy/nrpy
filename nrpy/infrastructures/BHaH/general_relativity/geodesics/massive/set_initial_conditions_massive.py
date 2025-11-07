"""
Generate the C function to set the full initial 8-component state vector.

Author: Dalton J. Moone
"""
import nrpy.c_function as cfc


def set_initial_conditions_massive() -> None:
    """
    Generate the C engine to set the full initial 8-component state vector.

    This orchestrator takes a particle's initial position, calls the
    `calculate_ut_uphi_from_r` helper, and then transforms the resulting u^phi
    into Cartesian 4-velocity components u^x and u^y.
    """
    name = "set_initial_conditions_massive"
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>"]
    desc = "Sets the initial 8-component state vector for a massive particle."
    params = """const particle_initial_state_t *restrict initial_state,
                  const commondata_struct *restrict commondata,
                  const params_struct *restrict params,
                  double y_out[8]"""
    body = r"""
    const double x = initial_state->pos[1];
    const double y = initial_state->pos[2];
    const double z = initial_state->pos[3];
    const double r = sqrt(x*x + y*y + z*z);
    double ut, uphi;
    calculate_ut_uphi_from_r(r, commondata, params, &ut, &uphi);
    if (isnan(ut)) {
        for(int i=0; i<8; i++) y_out[i] = NAN;
        return;
    }
    y_out[0] = initial_state->pos[0];
    y_out[1] = x;
    y_out[2] = y;
    y_out[3] = z;
    y_out[4] = ut;
    y_out[5] = -y * uphi;
    y_out[6] =  x * uphi;
    y_out[7] = 0.0;
    """

    cfc.register_CFunction(
        name=name,
        includes=includes,
        desc=desc,
        params=params,
        body=body,
    )
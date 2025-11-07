"""
Generates the C engine for handling a physical accretion disk intersection.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def handle_disk_intersection() -> None:
    """
    Generate and register the C engine for processing disk intersections.

    This function generates the C code that is called by the finalizer when a
    photon has terminated on the physical accretion disk. It orchestrates the
    full radiative transfer calculation by getting the metric, lowering indices,
    and calling the core radiative transfer engine.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "<math.h>"]
    desc = r"""@brief Performs the final radiative transfer calculation for a disk intersection.
    @details This engine is called by the finalizer. It takes the photon's final
             state and the stored nearest-neighbor data, computes the observed
             intensity and wavelength, and populates the final blueprint record.
    @param[in]  final_y              The photon's final 9-component state vector.
    @param[in]  nearest_neighbor     Pointer to the data of the disk particle that was hit.
    @param[in]  commondata           Pointer to commondata struct with runtime parameters.
    @param[in]  params               Pointer to params struct (unused, for signature compatibility).
    @param[in]  metric               Pointer to the metric_params struct specifying the metric type.
    @param[out] final_blueprint_data Pointer to the final output struct to be populated.
    """
    name = "handle_disk_intersection"
    params = """
    const double final_y[9],
    const MassiveParticle *restrict nearest_neighbor,
    const commondata_struct *restrict commondata, const params_struct *restrict params,
    const metric_params *restrict metric,
    blueprint_data_t *restrict final_blueprint_data
    """

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    body = r"""
    // --- Step 1: Get metric at the photon's final position (the intersection point) ---
    metric_struct g4DD;
    g4DD_metric(commondata, params, metric, final_y, &g4DD);

    // --- Step 2: Lower the indices of the photon's 4-momentum and the neighbor's 4-velocity ---
    // This is required for the dot products in the radiative transfer equation.
    const double g_munu[4][4] = {
        {g4DD.g00, g4DD.g01, g4DD.g02, g4DD.g03},
        {g4DD.g01, g4DD.g11, g4DD.g12, g4DD.g13},
        {g4DD.g02, g4DD.g12, g4DD.g22, g4DD.g23},
        {g4DD.g03, g4DD.g13, g4DD.g23, g4DD.g33}
    };

    double photon_p_mu[4] = {0,0,0,0};
    double disk_u_mu[4] = {0,0,0,0};
    for(int mu=0; mu<4; mu++) {
        for(int nu=0; nu<4; nu++) {
            photon_p_mu[mu] += g_munu[mu][nu] * final_y[nu+4]; // p_μ = g_μν p^ν
            disk_u_mu[mu] += g_munu[mu][nu] * nearest_neighbor->u[nu]; // u_μ = g_μν u^ν
        }
    }

    // --- Step 3: Call the core radiative transfer physics engine ---
    double temp_stokes_I;
    double temp_lambda_observed;
    calculate_radiative_transfer(photon_p_mu, disk_u_mu,
                                 nearest_neighbor->j_intrinsic, nearest_neighbor->lambda_rest,
                                 &temp_stokes_I, &temp_lambda_observed);

    // --- Step 4: Populate the final blueprint with the physical results ---
    final_blueprint_data->stokes_I = temp_stokes_I;
    final_blueprint_data->lambda_observed = temp_lambda_observed;

    // --- Step 5: Populate diagnostic information from the intersection ---
    final_blueprint_data->y_s = nearest_neighbor->pos[0]; // x-pos of neighbor
    final_blueprint_data->z_s = nearest_neighbor->pos[1]; // y-pos of neighbor
    final_blueprint_data->t_s = final_y[0]; // time of intersection
    final_blueprint_data->L_s = final_y[8]; // path length at intersection
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )
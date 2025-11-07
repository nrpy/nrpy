"""
Generates the C engine for relativistic radiative transfer.

Author: Dalton J. Moone
"""

import nrpy.c_function as cfc


def radiative_transfer_engine() -> None:
    """
    Generate and register the C engine for radiative transfer physics.

    This function generates the C code that implements the relativistic
    radiative transfer equation for an optically thin emitter. It calculates
    the observed intensity and wavelength based on the redshift factor, which
    accounts for both gravitational redshift and the relativistic Doppler effect.
    """
    # Per project standards, define local variables for all register_CFunction args.
    includes = ["BHaH_defines.h", "<math.h>"]
    desc = r"""@brief Calculates the observed intensity and wavelength from the disk and photon state.

    This engine implements the relativistic radiative transfer equation for an
    optically thin source. It computes the observed physical quantities based on
    the redshift factor g = E_obs / E_emit.

    @param[in]  photon_p_mu       The photon's covariant 4-momentum at the emission point.
    @param[in]  disk_u_mu         The fluid's covariant 4-velocity at the emission point.
    @param[in]  disk_j_intrinsic  The fluid's intrinsic emissivity.
    @param[in]  disk_lambda_rest  The rest-frame emission wavelength.
    @param[out] stokes_I          Pointer to store the calculated observed intensity.
    @param[out] lambda_observed   Pointer to store the calculated observed wavelength.
    """
    name = "calculate_radiative_transfer"
    params = """
    const double photon_p_mu[4], const double disk_u_mu[4],
    const float disk_j_intrinsic, const double disk_lambda_rest,
    double *stokes_I, double *lambda_observed
    """

    # The body is algorithmic, not symbolic, so it is defined as a raw C string.
    # Step 1: Implement the radiative transfer equations.
    # Redshift factor: g = (-p_μ u^μ)_obs / (-p_μ u^μ)_disk
    # Observed Intensity: I_obs = j_int * g^3
    # Observed Wavelength: λ_obs = λ_rest / g
    body = r"""
    // The observer is assumed to be at rest in the coordinate frame far away,
    // so their 4-velocity is u_obs^μ = (1, 0, 0, 0).
    // The metric is Minkowski far away, so u_obs_μ = (-1, 0, 0, 0).
    // The photon momentum is p_μ.
    // Therefore, (-p_μ u^μ)_obs = - (p_0 * -1) = p_0.
    const double p_mu_u_mu_obs = photon_p_mu[0];

    // Calculate (-p_μ u^μ)_disk
    const double p_mu_u_mu_disk = - (photon_p_mu[0] * disk_u_mu[0] +
                                     photon_p_mu[1] * disk_u_mu[1] +
                                     photon_p_mu[2] * disk_u_mu[2] +
                                     photon_p_mu[3] * disk_u_mu[3]);

    // Redshift (Doppler) factor g = E_obs / E_disk = (-p_μ u^μ)_obs / (-p_μ u^μ)_disk
    const double doppler_factor = p_mu_u_mu_obs / p_mu_u_mu_disk;

    // Observed intensity I_obs = j_intrinsic * g^3
    *stokes_I = disk_j_intrinsic * doppler_factor * doppler_factor * doppler_factor;

    // Observed wavelength λ_obs = λ_rest / g
    *lambda_observed = disk_lambda_rest / doppler_factor;
    """

    # Register the C function.
    cfc.register_CFunction(
        includes=includes, desc=desc, name=name, params=params, body=body
    )
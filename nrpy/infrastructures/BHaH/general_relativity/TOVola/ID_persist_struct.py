"""
Initialize and set BH@H's persistent initial data structure ID_struct for TOVola initial data.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import nrpy.params as par

par.register_CodeParameters(
    "REAL",
    "TOVola",
    [
        "initial_central_density",
        "poly_eos_K",
        "poly_eos_Gamma",
        "initial_ode_step_size",
        "min_step_size",
        "max_step_size",
        "ode_error_limit",
    ],
    [0.125, 1.0, 2.0, 1e-20, 1e-20, 1.0, 1e-6],
    commondata=True,
)
par.register_CodeParameters(
    "int",
    "TOVola",
    [
        "ode_max_steps",
        "interpolation_stencil_size",
        "max_interpolation_stencil_size",
    ],
    [5000000, 11, 13],
    commondata=True,
)


def ID_persist_str() -> str:
    """
    Return contents of ID_persist_struct for TOVola initial data.

    :return: ID_persist_struct contents.
    """
    return r"""
  // The following arrays store stellar information at all numpoints_arr radii:
  REAL *restrict r_Schw_arr; // Stellar radial coordinate in units of Schwarzschild radius
  REAL *restrict rho_baryon_arr; // Baryonic mass density
  REAL *restrict rho_energy_arr;    // Mass-energy density
  REAL *restrict P_arr;  // Pressure
  REAL *restrict M_arr;  // Integrated rest mass
  REAL *restrict expnu_arr;    // Metric quantity
  REAL *restrict exp4phi_arr;  // Metric quantity
  REAL *restrict r_iso_arr;  // Isotropic radial coordinate, in literature, sometimes called rbar
  int numpoints_arr; // Number of radii stored in the arrays
"""

"""
Python module for

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
import nrpy.c_function as cfc
import nrpy.params as par

par.register_CodeParameter(
    "char[100]",
    "TwoPunctures",
    "TP_BBH_description",
    "No name",
    commondata=True,
    add_to_parfile=True,
    add_to_set_CodeParameters_h=False,
)
par.register_CodeParameters(
    "int",
    "TwoPunctures",
    ["TP_npoints_A", "TP_npoints_B", "TP_npoints_phi"],
    [-1, -1, -1],
    commondata=True,
)
par.register_CodeParameters(
    "REAL",
    "TwoPunctures",
    ["TP_bare_mass_m", "TP_bare_mass_M"],
    [-1, -1],
    commondata=True,
)
par.register_CodeParameters(
    "REAL",
    __name__,
    [
        "bbhzx_BH_M_chix",
        "bbhzx_BH_M_chiy",
        "bbhzx_BH_M_chiz",
        "bbhzx_BH_m_chix",
        "bbhzx_BH_m_chiy",
        "bbhzx_BH_m_chiz",
    ],
    [0, 0, 0, 0, 0, 0],
    commondata=True,
)


def register_CFunction_initialize_ID_persist_struct() -> None:
    """
    Register C function initialize_ID_persist_struct(), which populates ID_persist_struct with defaults, and overrides defaults with commondata.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Initialize ID_persist_struct: populate some with defaults; set others with inputs from commondata;
    set up initial_p_t and initial_p_r if not set in parfile; and finally rotate inputs from xy to zx plane."""
    c_type = "void"
    name = "initialize_ID_persist_struct"
    params = "commondata_struct *restrict commondata, ID_persist_struct *restrict par"
    body = r"""  // Step 1: Set default values

  // psi^-2 initial lapse performs better than twopunctures-averaged.
  //   HOWEVER, "twopunctures-averaged" is default in the TwoPunctures
  //   Einstein Toolkit thorn.
  snprintf(par->initial_lapse, 100, "psi^n"); // Initial lapse.
  par->initial_lapse_psi_exponent = -2.0;     // Exponent n for psi^-n initial lapse profile

  // psi^-2 initial lapse is actually better,
  // The behavior of the lapse at the outer boundary is far worse with this choice.
  //      In fact you may witness code crashes.
  // snprintf(par->initial_lapse, 100, "twopunctures-averaged"); // Initial lapse.

  par->Newton_maxit = 5;       // Maximum number of Newton iterations
  par->Newton_tol = 1.0e-10;   // Tolerance for Newton solver
  par->TP_epsilon = 1e-6;      // A small number to smooth out singularities at the puncture locations
  par->adm_tol = 1.0e-10;      // Tolerance of ADM masses when give_bare_mass=no
  par->TP_Tiny = 0.0;          // Tiny number to avoid nans near or at the pucture locations
  par->TP_Extend_Radius = 0.0; // Radius of an extended spacetime instead of the puncture

  par->verbose = true;         // Print screen output while solving
  par->keep_u_around = false;  // Keep the variable u around after solving
  par->swap_xz = true;         // Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
  par->use_sources = false;    // Use sources?
  par->rescale_sources = true; // If sources are used - rescale them after solving?
  // par->use_external_initial_guess; // Set initial guess by external function?
  par->do_residuum_debug_output = false;               // Output debug information about the residuum
  par->do_initial_debug_output = false;                // Output debug information about initial guess
  par->multiply_old_lapse = false;                     // Multiply the old lapse with the new one
  par->solve_momentum_constraint = false;              // Solve for momentum constraint?
  snprintf(par->grid_setup_method, 100, "evaluation"); // How to fill the 3D grid from the spectral grid. evaluation is the only accurate one; don't use the other option!

  par->give_bare_mass = false; // User provides bare masses rather than target ADM masses

  // Any of the above parameters will be overwritten depending on the ID type.

  //fprintf(stderr, "Setting up %s TwoPunctures initial data...\n", commondata->TP_BBH_description);

  // Step 2: Set values from commondata
  // Step 2.a: Set TwoPunctures grid
  if (commondata->TP_npoints_A < 0 || commondata->TP_npoints_B < 0 || commondata->TP_npoints_phi < 0) {
    // Default TwoPunctures grid setup
    if (commondata->mass_ratio >= 2.0) {
      // higher mass ratios need higher spectral resolution on TwoPunctures grid.
      par->npoints_A = 66;   // Number of coefficients in the compactified radial direction
      par->npoints_B = 66;   // Number of coefficients in the angular direction
      par->npoints_phi = 28; // Number of coefficients in the phi direction
    } else {
      par->npoints_A = 48;   // Number of coefficients in the compactified radial direction
      par->npoints_B = 48;   // Number of coefficients in the angular direction
      par->npoints_phi = 20; // Number of coefficients in the phi direction
    }
  } else {
    par->npoints_A = commondata->TP_npoints_A;
    par->npoints_B = commondata->TP_npoints_B;
    par->npoints_phi = commondata->TP_npoints_phi;
  }

  // Step 2.b: Set initial tangential and radial momenta if initialized to -1.0
  //           (-1.0 is *exactly* representable in single/double precision)
  if (commondata->initial_p_t == -1.0 || commondata->initial_p_r == -1.0) {
    NRPyPN_quasicircular_momenta(commondata);
  }

  // Step 2.c: Prepare inputs for TwoPunctures, which will set up the binary
  //           to be orbiting on the *xz* plane. Here we convert all inputs
  //           to be consistent with this assumption.

  // IMPORTANT: The inputs for commondata->bbhxy_BH_{m,M}_chi{x,y,z}
  //   assume the BHs are initially (instantaneously) orbiting on
  //   the xy plane. So does NRPyPN:
  // However, TwoPunctures below assumes the BHs are orbiting
  //   in the xz plane (par->swap_xz).
  // IMPORTANT: Flipping x<->z will change the sign of the
  //   yhat direction. This sign flip is done at the bottom of
  //   this function.
  // Inputs: xy-plane. Outputs: zx-plane:
  //  z = x, x = y, y = z
  commondata->bbhzx_BH_m_chiz = commondata->bbhxy_BH_m_chix;
  commondata->bbhzx_BH_m_chix = commondata->bbhxy_BH_m_chiy;
  commondata->bbhzx_BH_m_chiy = commondata->bbhxy_BH_m_chiz;

  commondata->bbhzx_BH_M_chiz = commondata->bbhxy_BH_M_chix;
  commondata->bbhzx_BH_M_chix = commondata->bbhxy_BH_M_chiy;
  commondata->bbhzx_BH_M_chiy = commondata->bbhxy_BH_M_chiz;

  fprintf(stderr, "NRPyPN: Found p_t, p_r = %.8f %.8f\n", commondata->initial_p_t, commondata->initial_p_r);
  // For q=1, spins=0, diameter_of_separation=4.0, p_r is negative.
  //    This is likely due to the separation being too small.
  //    We assume below that p_r is positive, so let's fix it:
  if (commondata->initial_p_r < 0.0)
    commondata->initial_p_r *= -1.0;

  // UNIVERSAL PARAMETERS:
  {
    const REAL q = commondata->mass_ratio;
    const REAL p_t = commondata->initial_p_t;
    const REAL p_r = commondata->initial_p_r;
    commondata->mass_M = q / (1.0 + q);
    commondata->mass_m = 1.0 / (1.0 + q);
    par->target_M_plus = commondata->mass_M;  // MORE MASSIVE: target ADM mass for m+
    par->target_M_minus = commondata->mass_m; // LESS MASSIVE: target ADM mass for m-

    if (commondata->TP_bare_mass_M < 0 || commondata->TP_bare_mass_m < 0) {
      par->give_bare_mass = false;
    } else {
      par->give_bare_mass = true;
      par->par_m_plus = commondata->TP_bare_mass_M;
      par->par_m_minus = commondata->TP_bare_mass_m;
    }

    if (par->give_bare_mass == false) {
      // Set initial guesses for bare masses. Typically
      //   these are roughly 0.9*target_M_plus/minus, but
      //   the actual values may vary. The Newton-Raphson
      //   solver will sort it out, no worries :)
      par->par_m_plus = 0.9 * par->target_M_plus;
      par->par_m_minus = 0.9 * par->target_M_minus;
    }

    par->par_b = 0.5 * commondata->initial_sep; // z coordinate of the m+ puncture *on TwoPunctures grid*
    const REAL grid_dist_from_origin_BH_m = commondata->initial_sep * par->target_M_plus;
    const REAL grid_dist_from_origin_BH_M = commondata->initial_sep * par->target_M_minus;
    //   Initialize center_offset to zero before setting z component
    for (int ii = 0; ii < 3; ii++)
      par->center_offset[ii] = 0.0;
    par->center_offset[2] = -(grid_dist_from_origin_BH_m - grid_dist_from_origin_BH_M) * 0.5;

    // PUNCTURES ORBIT IN THE X-Z PLANE:
    //   Initialize linear momenta to zero before setting initial radial and tangential momenta
    for (int ii = 0; ii < 3; ii++)
      par->par_P_plus[ii] = par->par_P_minus[ii] = 0.0;
    par->par_P_plus[0] = -p_r;  // momentum of the m+ puncture
    par->par_P_minus[0] = +p_r; // momentum of the m- puncture

    par->par_P_plus[2] = +p_t;  // momentum of the m+ puncture
    par->par_P_minus[2] = -p_t; // momentum of the m- puncture
    {
      // Dimensionless spin parameter chi = J/M^2 --> J = chi * M^2
      par->par_S_minus[0] = commondata->bbhzx_BH_m_chix * par->target_M_minus * par->target_M_minus;
      par->par_S_minus[1] = commondata->bbhzx_BH_m_chiy * par->target_M_minus * par->target_M_minus;
      par->par_S_minus[2] = commondata->bbhzx_BH_m_chiz * par->target_M_minus * par->target_M_minus;

      par->par_S_plus[0] = commondata->bbhzx_BH_M_chix * par->target_M_plus * par->target_M_plus;
      par->par_S_plus[1] = commondata->bbhzx_BH_M_chiy * par->target_M_plus * par->target_M_plus;
      par->par_S_plus[2] = commondata->bbhzx_BH_M_chiz * par->target_M_plus * par->target_M_plus;
    }
    // Since we flip x<->z, the sign of the y-component of spin must
    //   flip in order to keep a right-handed coordinate system,
    //   consistent with assumptions in NRPyPN.
    par->par_S_minus[1] *= -1.0;
    par->par_S_plus[1] *= -1.0;
  }

  fprintf(stderr, "#################################\n");
  fprintf(stderr, "-={ INITIAL BINARY PARAMETERS }=-\n");
  fprintf(stderr, "M=1 (sum of individual ADM masses as defined in TwoPunctures)\n");
  fprintf(stderr, "d_initial/M = %.15f, q = %.15f\n", commondata->initial_sep, commondata->mass_ratio);
  fprintf(stderr, "bbhxy_BH_m_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_m_chix, commondata->bbhxy_BH_m_chiy, commondata->bbhxy_BH_m_chiz);
  fprintf(stderr, "bbhxy_BH_M_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_M_chix, commondata->bbhxy_BH_M_chiy, commondata->bbhxy_BH_M_chiz);
  fprintf(stderr, "p_t = %.15f, p_r = %.15f\n", commondata->initial_p_t, commondata->initial_p_r);
  fprintf(stderr, "TP resolution: %d  %d  %d\n", par->npoints_A, par->npoints_B, par->npoints_phi);
  fprintf(stderr, "#################################\n");
"""
    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        body=body,
    )

"""
Initialize and set BH@H's persistent initial data structure ID_struct for TwoPunctures initial data.

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


def ID_persist_str() -> str:
    """Return contents of ID_persist_struct for TwoPunctures initial data."""
    return r"""
  derivs v;    // stores coefficients
  derivs cf_v; // stores coefficients

  REAL mp,mm, mp_adm,mm_adm, E, J1,J2,J3;

  int npoints_A; //Number of coefficients in the compactified radial direction
  int npoints_B; //Number of coefficients in the angular direction
  int npoints_phi; //Number of coefficients in the phi direction
  int Newton_maxit; //Maximum number of Newton iterations

  REAL adm_tol; //Tolerance of ADM masses when give_bare_mass=no
  REAL Newton_tol; //Tolerance for Newton solver
  REAL TP_epsilon; //A small number to smooth out singularities at the puncture locations
  REAL TP_Tiny; //Tiny number to avoid nans near or at the pucture locations
  REAL TP_Extend_Radius; //Radius of an extended spacetime instead of the puncture
  REAL par_b; //x coordinate of the m+ puncture STEERABLE=always
  REAL par_m_plus; //mass of the m+ puncture STEERABLE = ALWAYS
  REAL par_m_minus; //mass of the m- puncture STEERABLE = ALWAYS
  REAL target_M_plus; //target ADM mass for m+
  REAL target_M_minus; //target ADM mass for m-
  REAL par_P_plus[3]; //momentum of the m+ puncture
  REAL par_P_minus[3]; //momentum of the m- puncture
  REAL par_S_plus[3]; //spin of the m+ puncture
  REAL par_S_minus[3]; //spin of the m- puncture
  REAL center_offset[3]; //offset b=0 to position (x,y,z)
  REAL initial_lapse_psi_exponent; //Exponent n for psi^-n initial lapse profile

  bool verbose; //Print screen output while solving
  bool keep_u_around; //Keep the variable u around after solving
  bool give_bare_mass; //User provides bare masses rather than target ADM masses
  bool swap_xz; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
  bool use_sources; //Use sources?
  bool rescale_sources; //If sources are used - rescale them after solving?
  //bool use_external_initial_guess; //Set initial guess by external function?
  bool do_residuum_debug_output; //Output debug information about the residuum
  bool do_initial_debug_output; //Output debug information about initial guess
  bool multiply_old_lapse; //Multiply the old lapse with the new one
  bool solve_momentum_constraint; //Solve for momentum constraint?

  char grid_setup_method[100]; //How to fill the 3D grid from the spectral grid
  char initial_lapse[100]; //initial lapse
"""


def register_CFunction_initialize_ID_persist_struct() -> None:
    """
    Register C function initialize_ID_persist_struct().

    This function populates ID_persist_struct with defaults, and overrides defaults with commondata.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """
IMPORTANT: We set up initial data in TwoPunctures assuming the BBH is initially in the xy-plane,
           as that is what TwoPunctures was designed to do in the Toolkit, AND this is what NRPyPN
           assumes. Also the "swap_xz" option in TwoPunctures results in a non-right-handed
           coordinate system. So we set up initial data in xy, then ***at the end*** we simply
           rotate xy to zx via: x->z , y->x, z->y.
Initialize ID_persist_struct: populate some with defaults; set others with inputs from commondata;
set up initial_p_t and initial_p_r if not set in parfile.
"""
    c_type = "void"
    name = "initialize_ID_persist_struct"
    params = "commondata_struct *restrict commondata, ID_persist_struct *restrict par"
    body = r"""  // Step 1: Set default TwoPunctures values

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
  par->swap_xz = false;         // Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
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

  // fprintf(stderr, "Setting up %s TwoPunctures initial data...\n", commondata->TP_BBH_description);

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

  if (commondata->initial_p_t == -1.0 || commondata->initial_p_r == -1.0) {
    // Step 2.b: Set initial tangential and radial momenta if initialized to -1.0
    //           (-1.0 is *exactly* representable in single/double precision)
    // The inputs for commondata->bbh_physical_params.chi_BH_{m,M}
    //   assume the BHs are initially (instantaneously) orbiting on
    //   the xy plane. So does NRPyPN:
    NRPyPN_quasicircular_momenta(commondata);

    fprintf(stderr, "NRPyPN: Found |p_t|, |p_r| = %.8f %.8f\n", fabs(commondata->initial_p_t), fabs(commondata->initial_p_r));
  }

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
    par->center_offset[0] = -(grid_dist_from_origin_BH_m - grid_dist_from_origin_BH_M) * 0.5;

    //   Initialize linear momenta to zero before setting initial radial and tangential momenta
    for (int ii = 0; ii < 3; ii++)
      par->par_P_plus[ii] = par->par_P_minus[ii] = 0.0;
    par->par_P_plus[0] = -fabs(p_r);  // momentum of the m+ puncture
    par->par_P_minus[0] = +fabs(p_r); // momentum of the m- puncture

    par->par_P_plus[1] = +fabs(p_t);  // momentum of the m+ puncture
    par->par_P_minus[1] = -fabs(p_t); // momentum of the m- puncture
    for (int ii = 0; ii < 3; ii++) {
      const REAL bbhxy_BH_m_chi[3] = {commondata->bbhxy_BH_m_chix, commondata->bbhxy_BH_m_chiy, commondata->bbhxy_BH_m_chiz };
      const REAL bbhxy_BH_M_chi[3] = {commondata->bbhxy_BH_M_chix, commondata->bbhxy_BH_M_chiy, commondata->bbhxy_BH_M_chiz };

      // Dimensionless spin parameter chi = J/M^2 --> J = chi * M^2
      par->par_S_minus[ii] = bbhxy_BH_m_chi[ii] * par->target_M_minus * par->target_M_minus;
      par->par_S_plus[ii] = bbhxy_BH_M_chi[ii] * par->target_M_plus * par->target_M_plus;
    }
  }

  fprintf(stderr, "#################################\n");
  fprintf(stderr, "-={ INITIAL BINARY PARAMETERS }=-\n");
  fprintf(stderr, "M=1 (sum of individual ADM masses as defined in TwoPunctures)\n");
  fprintf(stderr, "d_initial/M = %.15f, q = %.15f\n", commondata->initial_sep, commondata->mass_ratio);
  fprintf(stderr, "bbhxy_BH_m_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_m_chix, commondata->bbhxy_BH_m_chiy, commondata->bbhxy_BH_m_chiz);
  fprintf(stderr, "bbhxy_BH_M_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_M_chix, commondata->bbhxy_BH_M_chiy, commondata->bbhxy_BH_M_chiz);
  fprintf(stderr, "|p_t| = %.15f, |p_r| = %.15f\n", fabs(commondata->initial_p_t), fabs(commondata->initial_p_r));
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

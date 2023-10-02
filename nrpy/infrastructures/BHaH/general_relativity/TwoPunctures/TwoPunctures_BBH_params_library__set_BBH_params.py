"""
Python module for

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import namedtuple  # Standard Python: Enable namedtuple data type
from typing import Dict
import nrpy.c_function as cfc

bbh_params = namedtuple(
    "bbh_params",
    "q d M_chix M_chiy M_chiz m_chix m_chiy m_chiz p_t p_r bare_mass_m_plus bare_mass_m_minus NA NB Nphi",
)
bbh_params_dict: Dict[str, bbh_params] = {}


def add_to_bbh_params_dict(
    name,
    q=1.0,
    d=2.0,
    M_chix=0.0,
    M_chiy=0.0,
    M_chiz=0.0,
    m_chix=0.0,
    m_chiy=0.0,
    m_chiz=0.0,
    p_t=0.0,
    p_r=0.0,
    bare_mass_m_plus=-1.0,
    bare_mass_m_minus=-1.0,
    NA=-1,
    NB=-1,
    Nphi=-1,
):
    bbh_params_dict[name] = bbh_params(
        q,
        d,
        M_chix,
        M_chiy,
        M_chiz,
        m_chix,
        m_chiy,
        m_chiz,
        p_t,
        p_r,
        bare_mass_m_plus,
        bare_mass_m_minus,
        NA,
        NB,
        Nphi,
    )


add_to_bbh_params_dict(
    "GW150914ET",
    q=36.0 / 29.0,
    d=10.0,
    M_chiy=+0.31,
    m_chiy=-0.46,
    p_t=0.09530152296974252,
    p_r=0.00084541526517121,
    NA=30,
    NB=30,
    Nphi=16,
)
add_to_bbh_params_dict(
    "BHB_q1_chi_m_0.3__chi_M_0_sep_11p8",
    q=1.0,
    d=5.9 * 2.0,
    m_chiy=0.075 / (0.5 * 0.5),  # J=0.075 is given, and chi = J/m^2, with m=0.5
    p_t=0.0852865,
    p_r=0.000515022,
    bare_mass_m_plus=0.48811120218500131,
    bare_mass_m_minus=0.4698442439908046,
    NA=64,
    NB=64,
    Nphi=44,
)
add_to_bbh_params_dict(
    "q1sep10",
    q=1.0,
    d=10.0,
    p_t=0.0962578089026658,
    p_r=0.00100787185295814,  # from NRPyPN
    NA=44,
    NB=44,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q1sep8",
    q=1.0,
    d=8.0,
    p_t=0.112845235097096,
    p_r=0.00228434381143799,  # from NRPyPN
    NA=44,
    NB=44,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q4sep11RITBBH0088",
    q=4.0,
    d=11.0,
    p_t=0.0578913282748551,
    p_r=0.000302008048634052,  # from NRPyPN
    NA=48,
    NB=48,
    Nphi=20,
)
add_to_bbh_params_dict(
    "q6sep13_SXSBBH0166",
    q=6.0,
    d=13.0,
    p_t=0.0396619891606381,
    p_r=0.000101374427664241,  # from NRPyPN
    NA=66,
    NB=66,
    Nphi=28,
)
add_to_bbh_params_dict(
    "q4sep12p5_SXSBBH0167",
    q=4.0,
    d=12.5,
    p_t=0.0531333569625095,
    p_r=0.000195987020545860,  # from NRPyPN
    NA=48,
    NB=48,
    Nphi=20,
)
add_to_bbh_params_dict(
    "q1sep12-1810.00036",
    # FROM TABLE V OF https://arxiv.org/pdf/1810.00036.pdf
    q=1.0,
    d=12.0,
    p_t=0.850686e-1,
    p_r=0.468113e-3,
    NA=44,
    NB=44,
    Nphi=24,
)
add_to_bbh_params_dict(
    "QC0_mclachlan",
    # From einsteinanalysis/Extract/test/qc0-mclachlan.par:
    q=1.0,
    d=1.168642873 * 2,
    p_t=0.3331917498,
    p_r=0.0,
    bare_mass_m_plus=0.453,
    bare_mass_m_minus=0.453,
    NA=24,
    NB=24,
    Nphi=18,
)
add_to_bbh_params_dict(
    "QC0",
    q=1.0,
    d=1.168642873 * 2,
    p_t=0.3331917498,
    p_r=0.0,
    bare_mass_m_plus=0.453,
    bare_mass_m_minus=0.453,
)
add_to_bbh_params_dict(
    "TP_BL",
    # TwoPunctures version of Brill-Lindquist
    # (trivial solve, but beware orbital separation=0 won't work)
    q=3.0,
    d=14.5,
    p_t=0.0,
    p_r=0.0,
    NA=10,
    NB=10,
    Nphi=6,
)


def add_to_Cfunction_dict_TwoPunctures_BBH_params_library__set_BBH_params():
    includes = ["NRPy_basic_defines.h", "NRPy_function_prototypes.h"]
    desc = """This function provides both
  * a library of input BBH parameters for TwoPunctures, as well as
  * a means to generate quasicircular BBH input parameters from arbitrary physical inputs.
"""
    c_type = "void"
    name = "TwoPunctures_BBH_params_library__set_BBH_params"
    params = "commondata_struct *restrict commondata, char *restrict TP_ID_type, ID_persist_struct *restrict par"
    body = r"""  // First set default parameters:

  // psi^-2 initial lapse performs better than twopunctures-averaged.
  //   HOWEVER, "twopunctures-averaged" is default in the TwoPunctures
  //   Einstein Toolkit thorn.
  snprintf(par->initial_lapse, 100, "psi^n"); //Initial lapse.
  par->initial_lapse_psi_exponent = -2.0; //Exponent n for psi^-n initial lapse profile

  // psi^-2 initial lapse is actually better,
  // The behavior of the lapse at the outer boundary is far worse with this choice.
  //      In fact you may witness code crashes.
  //snprintf(par->initial_lapse, 100, "twopunctures-averaged"); //Initial lapse.

  par->Newton_maxit = 5; //Maximum number of Newton iterations
  par->Newton_tol = 1.0e-10; //Tolerance for Newton solver
  par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
  par->adm_tol = 1.0e-10; //Tolerance of ADM masses when give_bare_mass=no
  par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
  par->TP_Extend_Radius = 0.0; //Radius of an extended spacetime instead of the puncture

  par->verbose = true; //Print screen output while solving
  par->keep_u_around = false; //Keep the variable u around after solving
  par->swap_xz = true; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
  par->use_sources = false; //Use sources?
  par->rescale_sources = true; //If sources are used - rescale them after solving?
  //par->use_external_initial_guess; //Set initial guess by external function?
  par->do_residuum_debug_output = false; //Output debug information about the residuum
  par->do_initial_debug_output = false; //Output debug information about initial guess
  par->multiply_old_lapse = false; //Multiply the old lapse with the new one
  par->solve_momentum_constraint = false; //Solve for momentum constraint?
  snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid. evaluation is the only accurate one; don't use the other option!

  par->give_bare_mass = false; //User provides bare masses rather than target ADM masses

  // Any of the above parameters will be overwritten depending on the ID type.

  fprintf(stderr, "Setting up %s TwoPunctures initial data...\n", TP_ID_type);

  {
    // Default TwoPunctures grid setup
    if(commondata->bbh_physical_params.mass_ratio >= 2.0) {
      // higher mass ratios need higher spectral resolution on TwoPunctures grid.
      par->npoints_A   = 66; //Number of coefficients in the compactified radial direction
      par->npoints_B   = 66; //Number of coefficients in the angular direction
      par->npoints_phi = 28; //Number of coefficients in the phi direction
    } else {
      par->npoints_A   = 48; //Number of coefficients in the compactified radial direction
      par->npoints_B   = 48; //Number of coefficients in the angular direction
      par->npoints_phi = 20; //Number of coefficients in the phi direction
    }
    if(strcmp(TP_ID_type, "NRPyPN_cheapTP")==0) { // NRPyPN, but with "cheap" TwoPunctures
      par->npoints_A   = 30; //Number of coefficients in the compactified radial direction
      par->npoints_B   = 30; //Number of coefficients in the angular direction
      par->npoints_phi = 16; //Number of coefficients in the phi direction
    }
  }

  if(strcmp(TP_ID_type, "NRPyPN")==0) {
    // The inputs for commondata->bbh_physical_params.chi_BH_{m,M}
    //   assume the BHs are initially (instantaneously) orbiting on
    //   the xy plane. So does NRPyPN:
    set_physical_params_NRPyPN(commondata->bbh_physical_params.mass_ratio,
                               commondata->bbh_physical_params.chi_BH_M,
                               commondata->bbh_physical_params.chi_BH_m,
                               commondata->bbh_physical_params.initial_orbital_separation,
                               &commondata->bbh_physical_params.initial_p_t,
                               &commondata->bbh_physical_params.initial_p_r);

    // However, TwoPunctures below assumes the BHs are orbiting
    //   in the xz plane (par->swap_xz).
    // IMPORTANT: Flipping x<->z will change the sign of the
    //   yhat direction. This sign flip is done at the bottom of
    //   this function.
    // Inputs: xy-plane. Outputs: zx-plane:
    //  z = x, x = y, y = z
    REAL zx_plane_chi_BH_m[3], zx_plane_chi_BH_M[3];
    zx_plane_chi_BH_m[2] = commondata->bbh_physical_params.chi_BH_m[0];
    zx_plane_chi_BH_m[0] = commondata->bbh_physical_params.chi_BH_m[1];
    zx_plane_chi_BH_m[1] = commondata->bbh_physical_params.chi_BH_m[2];

    zx_plane_chi_BH_M[2] = commondata->bbh_physical_params.chi_BH_M[0];
    zx_plane_chi_BH_M[0] = commondata->bbh_physical_params.chi_BH_M[1];
    zx_plane_chi_BH_M[1] = commondata->bbh_physical_params.chi_BH_M[2];

    for(int ii=0;ii<3;ii++) {
      commondata->bbh_physical_params.chi_BH_m[ii] = zx_plane_chi_BH_m[ii];
      commondata->bbh_physical_params.chi_BH_M[ii] = zx_plane_chi_BH_M[ii];
    }

    fprintf(stderr, "NRPyPN: Found p_t, p_r = %.8f %.8f\n",
            commondata->bbh_physical_params.initial_p_t,
            commondata->bbh_physical_params.initial_p_r);
    // For q=1, spins=0, diameter_of_separation=4.0, p_r is negative.
    //    This is likely due to the separation being too small.
    //    We assume below that p_r is positive, so let's fix it:
    if(commondata->bbh_physical_params.initial_p_r < 0.0) commondata->bbh_physical_params.initial_p_r *= -1.0;
  }"""
    # Loop through the bbh_params_dict items to populate the body string
    for key, item in bbh_params_dict.items():
        body += f""" else if(strcmp(TP_ID_type, "{key}")==0) {{
        commondata->bbh_physical_params.mass_ratio = {item.q};
        commondata->bbh_physical_params.initial_orbital_separation = {item.d};
    """
        if item.M_chix != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_M[0] = {item.M_chix};\n"
            )
        if item.M_chiy != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_M[1] = {item.M_chiy};\n"
            )
        if item.M_chiz != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_M[2] = {item.M_chiz};\n"
            )

        if item.m_chix != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_m[0] = {item.m_chix};\n"
            )
        if item.m_chiy != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_m[1] = {item.m_chiy};\n"
            )
        if item.m_chiz != 0.0:
            body += (
                f"    commondata->bbh_physical_params.chi_BH_m[2] = {item.m_chiz};\n"
            )

        body += f"    commondata->bbh_physical_params.initial_p_t = {item.p_t};\n"
        body += f"    commondata->bbh_physical_params.initial_p_r = {item.p_r};\n"

        if item.bare_mass_m_plus != -1.0 and item.bare_mass_m_minus != -1.0:
            body += f"    par->give_bare_mass = true; //User provides bare masses rather than target ADM masses\n"
            body += f"    par->par_m_plus = {item.bare_mass_m_plus};\n"
            body += f"    par->par_m_minus = {item.bare_mass_m_minus};\n"

        if item.NA > 0:
            body += f"    par->npoints_A   = {item.NA};\n"
        if item.NB > 0:
            body += f"    par->npoints_B   = {item.NB};\n"
        if item.Nphi > 0:
            body += f"    par->npoints_phi = {item.Nphi};\n"

        body += "  }"

    body += r""" else {
        fprintf(stderr, "Error: did not recognize TwoPunctures ID type = %s\n", TP_ID_type);
        fprintf(stderr, "       Please pick one of the following library cases, or choose your own parameters:\n");
        fprintf(stderr,"""

    keys_str = ",".join([f"{key}\\n" for key in bbh_params_dict.keys()])
    body += f'"{keys_str}\\n"'
    body += r""");
        exit(1);
      }"""
    body += r"""

  // UNIVERSAL PARAMETERS:
  const REAL q   = commondata->bbh_physical_params.mass_ratio;
  const REAL p_t = commondata->bbh_physical_params.initial_p_t;
  const REAL p_r = commondata->bbh_physical_params.initial_p_r;
  commondata->bbh_physical_params.mass_M =   q / (1.0 + q);
  commondata->bbh_physical_params.mass_m = 1.0 / (1.0 + q);
  par->target_M_plus  = commondata->bbh_physical_params.mass_M; //MORE MASSIVE: target ADM mass for m+
  par->target_M_minus = commondata->bbh_physical_params.mass_m; //LESS MASSIVE: target ADM mass for m-
  if(par->give_bare_mass == false) {
    // Set initial guesses for bare masses. Typically
    //   these are roughly 0.9*target_M_plus/minus, but
    //   the actual values may vary. The Newton-Raphson
    //   solver will sort it out, no worries :)
    par->par_m_plus  = 0.9*par->target_M_plus;
    par->par_m_minus = 0.9*par->target_M_minus;
  }
  par->par_b = 0.5 * commondata->bbh_physical_params.initial_orbital_separation; //x coordinate of the m+ puncture STEERABLE=always
  const REAL grid_dist_from_origin_BH_m = commondata->bbh_physical_params.initial_orbital_separation*par->target_M_plus;
  const REAL grid_dist_from_origin_BH_M = commondata->bbh_physical_params.initial_orbital_separation*par->target_M_minus;
  //   Initialize center_offset to zero before setting z component
  for(int ii=0;ii<3;ii++) par->center_offset[ii] = 0.0;
  par->center_offset[2] = -(grid_dist_from_origin_BH_m - grid_dist_from_origin_BH_M) * 0.5;

  // PUNCTURES ORBIT IN THE X-Z PLANE:
  //   Initialize linear momenta to zero before setting initial radial and tangential momenta
  for(int ii=0;ii<3;ii++) par->par_P_plus[ii] = par->par_P_minus[ii] = 0.0;
  par->par_P_plus[0]    = -p_r; //momentum of the m+ puncture
  par->par_P_minus[0]   = +p_r; //momentum of the m- puncture

  par->par_P_plus[2]    = +p_t; //momentum of the m+ puncture
  par->par_P_minus[2]   = -p_t; //momentum of the m- puncture
  for(int ii=0;ii<3;ii++) {
    // Dimensionless spin parameter chi = J/M^2 --> J = chi * M^2
    par->par_S_minus[ii] = commondata->bbh_physical_params.chi_BH_m[ii] * par->target_M_minus * par->target_M_minus;
    par->par_S_plus[ii]  = commondata->bbh_physical_params.chi_BH_M[ii] * par->target_M_plus  * par->target_M_plus;
  }
  // Since we flip x<->z, the sign of the y-component of spin must
  //   flip in order to keep a right-handed coordinate system,
  //   consistent with assumptions in NRPyPN.
  par->par_S_minus[1] *= -1.0;
  par->par_S_plus[1]  *= -1.0;

  // Copy mass ratio information to mass_ratio_for_grid_centroid_offset, which is used for setting up the grids.
  commondata->mass_ratio_for_grid_centroid_offset = commondata->bbh_physical_params.mass_ratio;
"""
    cfc.register_CFunction(
        includes=includes, desc=desc, c_type=c_type, name=name, params=params, body=body
    )

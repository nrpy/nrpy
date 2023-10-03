"""
Python module for registering all TwoPunctures functions within NRPy+'s
    CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from pathlib import Path
import shutil

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import ID_persist_struct
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import CoordTransf
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Equations
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import FuncAndJacobian
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Newton
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_interp
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_solve
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_utilities


def copy_TwoPunctures_h(project_Path: Path) -> None:
    """
    Copies TwoPunctures.h into project directory.

    :param project_Path: The path of the project directory where the file will be copied.
    """
    try:
        # only Python 3.7+ has importlib.resources
        import importlib.resources  # pylint: disable=E1101,C0415

        source_path = (
            # pylint: disable=E1101
            importlib.resources.files(
                "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures"
            )
            / "TwoPunctures.h"
        )
        shutil.copy(str(source_path), str(project_Path))
    except ImportError:  # Fallback to resource_filename for older Python versions
        # pylint: disable=E1101,C0415
        from pkg_resources import resource_filename  # type: ignore

        source_path = resource_filename(
            "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures",
            "TwoPunctures.h",
        )
        shutil.copy(source_path, str(project_Path))  # type: ignore


def ID_persist_str():
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


def register_CFunction_bbh_initial_params():
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = """Output initial data parameters for binary black hole system."""
    c_type = "void"
    name = "bbh_initial_params"
    params = (
        "ID_persist_struct *restrict ID_persist, commondata_struct *restrict commondata"
    )
    body = r"""  fprintf(stderr, "#################################\n");
  fprintf(stderr, "-={ INITIAL BINARY PARAMETERS }=-\n");
  fprintf(stderr, "M=1 (sum of individual ADM masses as defined in TwoPunctures)\n");
  fprintf(stderr, "d_initial/M = %.15f, q = %.15f\n", commondata->initial_orbital_separation,commondata->mass_ratio);
  fprintf(stderr, "bbhxy_BH_m_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_m_chi[0],commondata->bbhxy_BH_m_chi[1],commondata->bbhxy_BH_m_chi[2]);
  fprintf(stderr, "bbhxy_BH_M_chi = %.15f %.15f %.15f\n", commondata->bbhxy_BH_M_chi[0],commondata->bbhxy_BH_M_chi[1],commondata->bbhxy_BH_M_chi[2]);
  fprintf(stderr, "p_t = %.15f, p_r = %.15f\n", commondata->initial_p_t,commondata->initial_p_r);
  fprintf(stderr, "TP resolution: %d  %d  %d\n", ID_persist->npoints_A,ID_persist->npoints_B,ID_persist->npoints_phi);
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


def register_CFunction_TwoPunctures_initialize_parameters():
    includes = [
        "../NRPy_basic_defines.h",
        "../NRPy_function_prototypes.h",
        "TwoPunctures.h",
    ]
    desc = (
        "Initialize either prescribed or custom (NRPyPN-based) TwoPunctures parameters"
    )
    name = "TwoPunctures_initialize_parameters"
    params = """char *ID_type, ID_persist_struct *par,
                                        const REAL q, const REAL chi_BH_m[3], const REAL chi_BH_M[3], const REAL diameter_of_separation,
                                        const REAL p_t_override, const REAL p_r_override"""
    body = r"""
  // SET DEFAULT INPUT PARAMETERS
  {
    par->npoints_A = 30; //Number of coefficients in the compactified radial direction
    par->npoints_B = 30; //Number of coefficients in the angular direction
    par->npoints_phi = 16; //Number of coefficients in the phi direction
    par->Newton_maxit = 5; //Maximum number of Newton iterations

    par->adm_tol = 1.0e-10; //Tolerance of ADM masses when give_bare_mass=no
    par->Newton_tol = 1.0e-10; //Tolerance for Newton solver
    par->TP_epsilon = 0.0; //A small number to smooth out singularities at the puncture locations
    par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
    par->TP_Extend_Radius = 0.0; //Radius of an extended spacetime instead of the puncture
    par->par_b = 1.0; //x coordinate of the m+ puncture STEERABLE=always
    par->par_m_plus    =   q / (1.0 + q); //mass of the m+ puncture STEERABLE = ALWAYS
    par->par_m_minus   = 1.0 / (1.0 + q); //mass of the m- puncture STEERABLE = ALWAYS
    par->target_M_plus =   q / (1.0 + q); //target ADM mass for m+
    par->target_M_minus = 1.0 / (1.0 + q); //target ADM mass for m-
    for(int i=0;i<3;i++) par->par_P_plus[i]    = 0.0; //momentum of the m+ puncture
    for(int i=0;i<3;i++) par->par_P_minus[i]   = 0.0; //momentum of the m- puncture
    for(int i=0;i<3;i++) par->par_S_plus[i]    = 0.0; //spin of the m+ puncture
    for(int i=0;i<3;i++) par->par_S_minus[i]   = 0.0; //spin of the m- puncture
    for(int i=0;i<3;i++) par->center_offset[i] = 0.0; //offset b=0 to position (x,y,z)
    // FIXME: Probably want to experiment with psi^{-2}; might
    //        reduce eccentricity and/or constraint violations.
    //snprintf(par->initial_lapse,100,"psi^n"); //Initial lapse.
    par->initial_lapse_psi_exponent = -2.0; //Exponent n for psi^-n initial lapse profile

    // To be more consistent with standard assumptions
    //   (e.g., made in Einstein Toolkit examples):
    snprintf(par->initial_lapse,100,"twopunctures-averaged"); //Initial lapse.

    par->verbose = true; //Print screen output while solving
    par->keep_u_around = false; //Keep the variable u around after solving
    par->give_bare_mass = false; //User provides bare masses rather than target ADM masses
    par->swap_xz = true; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
    par->use_sources = false; //Use sources?
    par->rescale_sources = true; //If sources are used - rescale them after solving?
    //par->use_external_initial_guess; //Set initial guess by external function?
    par->do_residuum_debug_output = false; //Output debug information about the residuum
    par->do_initial_debug_output = false; //Output debug information about initial guess
    par->multiply_old_lapse = false; //Multiply the old lapse with the new one
    par->solve_momentum_constraint = false; //Solve for momentum constraint?
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
  }

  if(strcmp(ID_type, "NRPyPN")==0) {
    REAL p_t, p_r;
    set_physical_params_NRPyPN(q, chi_BH_m, chi_BH_M, diameter_of_separation, &p_t, &p_r);

    par->Newton_tol = 1.0e-10; //Tolerance for Newton solver
    par->verbose = true; //Print screen output while solving
    par->npoints_A   = 48; //Number of coefficients in the compactified radial direction
    par->npoints_B   = 48; //Number of coefficients in the angular direction
    par->npoints_phi = 16; //Number of coefficients in the phi direction
    par->target_M_plus  =   q / (1.0 + q); //MORE MASSIVE: target ADM mass for m+
    par->target_M_minus = 1.0 / (1.0 + q); //LESS MASSIVE: target ADM mass for m-
    par->par_b = diameter_of_separation * 0.5 / q; //x coordinate of the m+ puncture STEERABLE=always

    // PUNCTURES ORBIT IN THE X-Z PLANE:
    par->par_P_plus[0]    = -p_r; //momentum of the m+ puncture
    par->par_P_minus[0]   = +p_r; //momentum of the m- puncture

    par->par_P_plus[2]    = +p_t; //momentum of the m+ puncture
    par->par_P_minus[2]   = -p_t; //momentum of the m- puncture
    fprintf(stderr, "NRPyPN: Found p_r, p_t = %e %e\n", p_r, p_t);
    for(int ii=0;ii<3;ii++) {
      // Dimensionless spin parameter chi = J/M^2 --> J = chi * M^2
      par->par_S_minus[ii] = chi_BH_m[ii] * par->target_M_minus * par->target_M_minus;
      par->par_S_plus[ii]  = chi_BH_M[ii] * par->target_M_plus  * par->target_M_plus;
    }
    // Since we flip x<->z, the sign of the y-component of spin must
    //   flip in order to keep a right-handed coordinate system,
    //   consistent with assumptions in NRPyPN.
    par->par_S_minus[1] *= -1.0;
    par->par_S_plus[1]  *= -1.0;

  } else if(strcmp(ID_type, "BrillLindquist")==0) {
    par->target_M_plus  =   q / (1.0 + q); //MORE MASSIVE: target ADM mass for m+
    par->target_M_minus = 1.0 / (1.0 + q); //LESS MASSIVE: target ADM mass for m-
    par->par_b = diameter_of_separation * 0.5 / q; //x coordinate of the m+ puncture STEERABLE=always
    // Momenta are zero.
  } else if(strcmp(ID_type, "QC0")==0) {
    //  equal-mass, one BH has dimensionless spin param = +0.3
    fprintf(stderr, "SETTING UP QC0 INITIAL DATA\n");

    // FROM ETK's einsteinanalysis/Extract/test/qc0-mclachlan.par :
    snprintf(par->initial_lapse,100,"twopunctures-averaged"); //Initial lapse.
    // INCREASED FROM DEFAULTS IN THE .par file:
    //par->npoints_A = 44; //Number of coefficients in the compactified radial direction
    //par->npoints_B = 44; //Number of coefficients in the angular direction
    //par->npoints_phi = 24; //Number of coefficients in the phi direction

    par->npoints_A = 24; //Number of coefficients in the compactified radial direction
    par->npoints_B = 24; //Number of coefficients in the angular direction
    par->npoints_phi = 18; //Number of coefficients in the phi direction

    // END MODIFICATIONS
    par->par_b = 1.168642873; //x coordinate of the m+ puncture STEERABLE=always
    par->give_bare_mass = true; //User provides bare masses rather than target ADM masses
    par->par_m_plus = 0.453; //mass of the m+ puncture STEERABLE = ALWAYS
    par->par_m_minus = 0.453; //mass of the m- puncture STEERABLE = ALWAYS
    // MODIFIED TO MOVE UP THE Z-AXIS INSTEAD OF UP THE Y-AXIS:
    // PUNCTURES ORBIT IN X-Z PLANE:
    par->par_P_plus[2]    =  0.3331917498; //momentum of the m+ puncture
    par->par_P_minus[2]   = -0.3331917498; //momentum of the m- puncture
    ////////par->swap_xz = true; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
    // ORIGINAL:
    //par->par_P_plus[1]    =  0.3331917498; //momentum of the m+ puncture
    //par->par_P_minus[1]   = -0.3331917498; //momentum of the m- puncture

    // ADJUSTED FROM DEFAULTS IN THE .par file:
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
    // END MODIFICATIONS

    // DECREASED FROM DEFAULTS IN THE .par file:
    par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
    par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
    // END MODIFICATIONS
    //  equal-mass, one BH has dimensionless spin param = +0.3
  } else if(strcmp(ID_type, "BHB with one 0.3 spin BH")==0) {
    // OVERRIDE CERTAIN DEFAULTS SO WE HAVE A BHB:
    //  equal-mass, one BH has dimensionless spin param = +0.3
    fprintf(stderr, "SETTING UP BHB WITH ONE-0.3 SPIN BH INITIAL DATA\n");

    par->Newton_tol = 1.0e-10; //Tolerance for Newton solver
    par->verbose = true; //Print screen output while solving
    par->npoints_A = 64; //Number of coefficients in the compactified radial direction
    par->npoints_B = 64; //Number of coefficients in the angular direction
    par->npoints_phi = 44; //Number of coefficients in the phi direction
    par->target_M_plus = 0.5; //target ADM mass for m+
    par->target_M_minus = 0.5; //target ADM mass for m-
    par->par_m_plus  = 0.48811120218500131; //mass of the m+ puncture STEERABLE = ALWAYS
    par->par_m_minus = 0.4698442439908046; //mass of the m- puncture STEERABLE = ALWAYS
    par->par_b = 5.9; //x coordinate of the m+ puncture STEERABLE=always
    par->par_P_plus[0]    = -0.000515022; //momentum of the m+ puncture
    par->par_P_minus[0]   =  0.000515022; //momentum of the m- puncture
    par->par_P_plus[1]    =  0.0852865; //momentum of the m+ puncture
    par->par_P_minus[1]   = -0.0852865; //momentum of the m- puncture
    par->par_S_minus[2]   = 0.075; //spin of the m- puncture: dimensionless spin parameter J/M^2 = 0.075/(0.5^2) = 0.3
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
    par->give_bare_mass = true; //User provides bare masses rather than target ADM masses
    par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
  } else if(strcmp(ID_type, "q1sep10xz")==0) {
    fprintf(stderr, "SETTING UP %s INITIAL DATA\n", ID_type);

    // FROM ETK's einsteinanalysis/Extract/test/qc0-mclachlan.par :
    snprintf(par->initial_lapse,100,"twopunctures-averaged"); //Initial lapse.
    // INCREASED FROM DEFAULTS IN THE .par file:
    par->npoints_A = 44; //Number of coefficients in the compactified radial direction
    par->npoints_B = 44; //Number of coefficients in the angular direction
    par->npoints_phi = 28; //Number of coefficients in the phi direction
    // ORIGINAL, res a bit too low
    //par->npoints_A = 24; //Number of coefficients in the compactified radial direction
    //par->npoints_B = 24; //Number of coefficients in the angular direction
    //par->npoints_phi = 18; //Number of coefficients in the phi direction

    // END MODIFICATIONS
    par->par_b = 5.0; //x coordinate of the m+ puncture STEERABLE=always
    par->give_bare_mass = false; //User provides bare masses rather than target ADM masses
    par->target_M_plus = 0.5; //target ADM mass for m+
    par->target_M_minus = 0.5; //target ADM mass for m-

    // MODIFIED TO MOVE UP THE Z-AXIS INSTEAD OF UP THE Y-AXIS:
    // PUNCTURES ORBIT IN X-Z PLANE:
    // Be sure that P_r is inward directed, and that the center of mass is set at r=0.
    // |P_t| =  0.0962578089026658
    // |P_r| =  0.00100787185295814
    par->swap_xz = true; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
    par->par_P_plus[0]    = -0.00100787185295814; // x (actually z)-momentum of the m+ puncture
    par->par_P_minus[0]   = +0.00100787185295814; // x (actually z)-momentum of the m- puncture
    par->par_P_plus[2]    = +0.0962578089026658;  // z (actually x)-momentum of the m+ puncture
    par->par_P_minus[2]   = -0.0962578089026658;  // z (actually x)-momentum of the m- puncture

    // ORIGINAL:
    ///////par->par_P_plus[1]    =  0.3331917498; //momentum of the m+ puncture
    ///////par->par_P_minus[1]   = -0.3331917498; //momentum of the m- puncture

    // ADJUSTED FROM DEFAULTS IN THE .par file:
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
    // END MODIFICATIONS

    // DECREASED FROM DEFAULTS IN THE .par file:
    par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
    par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
    // END MODIFICATIONS
  } else if(strcmp(ID_type, "q1sep8xz")==0) {
    fprintf(stderr, "SETTING UP %s INITIAL DATA\n", ID_type);

    // FROM ETK's einsteinanalysis/Extract/test/qc0-mclachlan.par :
    snprintf(par->initial_lapse,100,"twopunctures-averaged"); //Initial lapse.
    // INCREASED FROM DEFAULTS IN THE .par file:
    par->npoints_A = 44; //Number of coefficients in the compactified radial direction
    par->npoints_B = 44; //Number of coefficients in the angular direction
    par->npoints_phi = 28; //Number of coefficients in the phi direction

    // END MODIFICATIONS
    par->par_b = 4.0; //x coordinate of the m+ puncture STEERABLE=always
    par->give_bare_mass = false; //User provides bare masses rather than target ADM masses
    par->target_M_plus = 0.5; //target ADM mass for m+
    par->target_M_minus = 0.5; //target ADM mass for m-

    // |P_t| =  0.112845235097096
    // |P_r| =  0.00228434381143799
    par->swap_xz = true; //Swap x and z coordinates when interpolating, so that the black holes are separated in the z direction
    par->par_P_plus[0]    = -0.00228434381143799; // x (actually z)-momentum of the m+ puncture
    par->par_P_minus[0]   = +0.00228434381143799; // x (actually z)-momentum of the m- puncture
    par->par_P_plus[2]    = +0.112845235097096;  // z (actually x)-momentum of the m+ puncture
    par->par_P_minus[2]   = -0.112845235097096;  // z (actually x)-momentum of the m- puncture

    // ORIGINAL:
    ///////par->par_P_plus[1]    =  0.3331917498; //momentum of the m+ puncture
    ///////par->par_P_minus[1]   = -0.3331917498; //momentum of the m- puncture

    // ADJUSTED FROM DEFAULTS IN THE .par file:
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
    // END MODIFICATIONS

    // DECREASED FROM DEFAULTS IN THE .par file:
    par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
    par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
    // END MODIFICATIONS
  }  else if(strcmp(ID_type, "QC0separation12-1810.00036")==0) {
    fprintf(stderr,"SETTING UP %s INITIAL DATA\n",ID_type);

    // FROM ETK's einsteinanalysis/Extract/test/qc0-mclachlan.par :
    snprintf(par->initial_lapse,100,"twopunctures-averaged"); //Initial lapse.
    // INCREASED FROM DEFAULTS IN THE .par file:

    par->npoints_A = 44; //Number of coefficients in the compactified radial direction
    par->npoints_B = 44; //Number of coefficients in the angular direction
    par->npoints_phi = 24; //Number of coefficients in the phi direction
    //par->npoints_A = 24; //Number of coefficients in the compactified radial direction
    //par->npoints_B = 24; //Number of coefficients in the angular direction
    //par->npoints_phi = 18; //Number of coefficients in the phi direction

    // END MODIFICATIONS
    par->par_b = 6.0; //x coordinate of the m+ puncture STEERABLE=always
    par->give_bare_mass = false; //User provides bare masses rather than target ADM masses
    par->target_M_plus = 0.5; //target ADM mass for m+
    par->target_M_minus = 0.5; //target ADM mass for m-
    // MODIFIED TO MOVE UP THE Z-AXIS INSTEAD OF UP THE Y-AXIS:
    // PUNCTURES ORBIT IN X-Z PLANE:
    // FROM TABLE V OF https://arxiv.org/pdf/1810.00036.pdf
    par->par_P_plus[2]    = +0.850686e-1; //tangential momentum of the m+ puncture
    par->par_P_minus[2]   = -0.850686e-1; //tangential momentum of the m- puncture

    par->par_P_plus[0]    = -0.468113e-3; //radial momentum of the m+ puncture
    par->par_P_minus[0]   = +0.468113e-3; //radial momentum of the m- puncture

    // ADJUSTED FROM DEFAULTS IN THE .par file:
    snprintf(par->grid_setup_method,100,"evaluation"); //How to fill the 3D grid from the spectral grid
    // END MODIFICATIONS

    // DECREASED FROM DEFAULTS IN THE .par file:
    par->TP_epsilon = 1e-6; //A small number to smooth out singularities at the puncture locations
    par->TP_Tiny = 0.0; //Tiny number to avoid nans near or at the pucture locations
    // END MODIFICATIONS
  } else {
    fprintf(stderr, "SETTING UP DEFAULT INITIAL DATA\n");
  }
"""
    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        includes=includes,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )


def register_C_functions():
    register_CFunction_bbh_initial_params()
    register_CFunction_TwoPunctures_initialize_parameters()

    ID_persist_struct.register_CFunction_initialize_ID_persist_struct()

    # Original TwoPunctures functions:
    CoordTransf.register_CFunction_TP_CoordTransf()
    Equations.register_CFunction_TP_Equations()
    FuncAndJacobian.register_CFunction_FuncAndJacobian()
    Newton.register_CFunction_TP_Newton()
    TP_interp.register_CFunction_TP_Interp()
    TP_solve.register_CFunction_TP_solve()
    TP_utilities.register_CFunction_TP_utilities()

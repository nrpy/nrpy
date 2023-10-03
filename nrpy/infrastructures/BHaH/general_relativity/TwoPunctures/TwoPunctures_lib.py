"""
Python module for registering all TwoPunctures functions within NRPy+'s
    CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from pathlib import Path
import shutil

from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import ID_persist_struct
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import CoordTransf
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Equations
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import FuncAndJacobian
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Newton
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_interp
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_solve
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_utilities


def copy_TwoPunctures_header_files(TwoPunctures_Path: Path) -> None:
    """
    Copies TwoPunctures.h and TP_utilities.h into project directory.

    :param project_Path: The path of the project directory where the file will be copied.
    """
    TwoPunctures_Path.mkdir(parents=True, exist_ok=True)

    try:
        # only Python 3.7+ has importlib.resources
        from importlib import resources  # pylint: disable=E1101,C0415

        for header_file in ["TwoPunctures.h", "TP_utilities.h"]:
            source_path = (
                # pylint: disable=E1101
                resources.files(
                    "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures"
                )
                / header_file
            )
            shutil.copy(str(source_path), str(TwoPunctures_Path / header_file))
    except ImportError:  # Fallback to resource_filename for older Python versions
        # pylint: disable=E1101,C0415
        from pkg_resources import resource_filename  # type: ignore

        for header_file in ["TwoPunctures.h", "TP_utilities.h"]:
            source_path = resource_filename(
                "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures",
                header_file,
            )
            shutil.copy(source_path, str(TwoPunctures_Path))  # type: ignore


def ID_persist_str() -> str:
    """
    Return contents of ID_persist_struct for TwoPunctures initial data.
    """
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


def register_C_functions() -> None:
    """Register all C functions needed for TwoPunctures solve."""
    ID_persist_struct.register_CFunction_initialize_ID_persist_struct()

    # Original TwoPunctures functions:
    CoordTransf.register_CFunction_TP_CoordTransf()
    Equations.register_CFunction_TP_Equations()
    FuncAndJacobian.register_CFunction_TP_FuncAndJacobian()
    Newton.register_CFunction_TP_Newton()
    TP_interp.register_CFunction_TP_Interp()
    TP_solve.register_CFunction_TP_solve()
    TP_utilities.register_CFunction_TP_utilities()

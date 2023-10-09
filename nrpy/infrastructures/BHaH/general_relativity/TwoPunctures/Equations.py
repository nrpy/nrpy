"""
Register TwoPunctures code Equations.c.

Note that this C file is just a collection of functions for
basic equations. The registered function TP_Equations() is unused.

TwoPunctures creates initial for two puncture black holes using a
single domain spectral method.  This method is described in
Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy, "A single-domain
spectral method for black hole puncture data", PRD 70, 064011 (2004),
arXiv:gr-qc/0404056.

Code originally from Einstein Toolkit's EinsteinInitialData
arrangement:
https://bitbucket.org/einsteintoolkit/einsteininitialdata

License: Lesser GNU Public License, version 2.0+

Authors: Marcus Ansorg
         Erik Schnetter
         Frank Löffler
         Zachariah B. Etienne (pasted original code into Python strings)
         zachetie **at** gmail **dot* com
"""
import nrpy.c_function as cfc


def register_CFunction_TP_Equations() -> None:
    """Register C function TP_Equations(), provides basic equations functions for TwoPunctures solver."""
    desc = "Equations.c from TwoPunctures. Note that this C file is just a collection of functions for basic equations. TP_Equations() is unused."
    # prefunc contains most of the source code
    prefunc = f"// {desc}\n\n"
    prefunc += r"""/* TwoPunctures:  File  "Equations.c"*/

#include "TP_utilities.h"
#include "TwoPunctures.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) */
/* U.d1[ivar]   = U[ivar]_x;  */
/* U.d2[ivar]   = U[ivar]_y;  */
/* U.d3[ivar]   = U[ivar]_z;  */
/* U.d11[ivar]  = U[ivar]_xx; */
/* U.d12[ivar]  = U[ivar]_xy; */
/* U.d13[ivar]  = U[ivar]_xz;*/
/* U.d22[ivar]  = U[ivar]_yy;*/
/* U.d23[ivar]  = U[ivar]_yz;*/
/* U.d33[ivar]  = U[ivar]_zz;*/

REAL BY_KKofxyz(ID_persist_struct par, REAL x, REAL y, REAL z) {
  int i, j;
  REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm, Aij, AijAij, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par.par_b) * (x - par.par_b) + y * y + z * z;
  r2_minus = (x + par.par_b) * (x + par.par_b) + y * y + z * z;
  r_plus = sqrt(r2_plus);
  r_minus = sqrt(r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par.par_b) / r_plus;
  n_minus[0] = (x + par.par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++) {
    np_Pp += n_plus[i] * par.par_P_plus[i];
    nm_Pm += n_minus[i] * par.par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par.par_S_plus[2] - n_plus[2] * par.par_S_plus[1];
  np_Sp[1] = n_plus[2] * par.par_S_plus[0] - n_plus[0] * par.par_S_plus[2];
  np_Sp[2] = n_plus[0] * par.par_S_plus[1] - n_plus[1] * par.par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par.par_S_minus[2] - n_minus[2] * par.par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par.par_S_minus[0] - n_minus[0] * par.par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par.par_S_minus[1] - n_minus[1] * par.par_S_minus[0];
  AijAij = 0;
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) { /* Bowen-York-Curvature :*/
      Aij = +1.5 * (par.par_P_plus[i] * n_plus[j] + par.par_P_plus[j] * n_plus[i] + np_Pp * n_plus[i] * n_plus[j]) / r2_plus +
            1.5 * (par.par_P_minus[i] * n_minus[j] + par.par_P_minus[j] * n_minus[i] + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus - 3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus -
            3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
        Aij -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
      AijAij += Aij * Aij;
    }
  }

  return AijAij;
}

void BY_Aijofxyz(ID_persist_struct par, REAL x, REAL y, REAL z, REAL Aij[3][3]) {
  int i, j;
  REAL r_plus, r2_plus, r3_plus, r_minus, r2_minus, r3_minus, np_Pp, nm_Pm, n_plus[3], n_minus[3], np_Sp[3], nm_Sm[3];

  r2_plus = (x - par.par_b) * (x - par.par_b) + y * y + z * z;
  r2_minus = (x + par.par_b) * (x + par.par_b) + y * y + z * z;
  r2_plus = sqrt(pow(r2_plus, 2) + pow(par.TP_epsilon, 4));
  r2_minus = sqrt(pow(r2_minus, 2) + pow(par.TP_epsilon, 4));
  if (r2_plus < pow(par.TP_Tiny, 2))
    r2_plus = pow(par.TP_Tiny, 2);
  if (r2_minus < pow(par.TP_Tiny, 2))
    r2_minus = pow(par.TP_Tiny, 2);
  r_plus = sqrt(r2_plus);
  r_minus = sqrt(r2_minus);
  r3_plus = r_plus * r2_plus;
  r3_minus = r_minus * r2_minus;

  n_plus[0] = (x - par.par_b) / r_plus;
  n_minus[0] = (x + par.par_b) / r_minus;
  n_plus[1] = y / r_plus;
  n_minus[1] = y / r_minus;
  n_plus[2] = z / r_plus;
  n_minus[2] = z / r_minus;

  /* dot product: np_Pp = (n_+).(P_+); nm_Pm = (n_-).(P_-) */
  np_Pp = 0;
  nm_Pm = 0;
  for (i = 0; i < 3; i++) {
    np_Pp += n_plus[i] * par.par_P_plus[i];
    nm_Pm += n_minus[i] * par.par_P_minus[i];
  }
  /* cross product: np_Sp[i] = [(n_+) x (S_+)]_i; nm_Sm[i] = [(n_-) x (S_-)]_i*/
  np_Sp[0] = n_plus[1] * par.par_S_plus[2] - n_plus[2] * par.par_S_plus[1];
  np_Sp[1] = n_plus[2] * par.par_S_plus[0] - n_plus[0] * par.par_S_plus[2];
  np_Sp[2] = n_plus[0] * par.par_S_plus[1] - n_plus[1] * par.par_S_plus[0];
  nm_Sm[0] = n_minus[1] * par.par_S_minus[2] - n_minus[2] * par.par_S_minus[1];
  nm_Sm[1] = n_minus[2] * par.par_S_minus[0] - n_minus[0] * par.par_S_minus[2];
  nm_Sm[2] = n_minus[0] * par.par_S_minus[1] - n_minus[1] * par.par_S_minus[0];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) { /* Bowen-York-Curvature :*/
      Aij[i][j] = +1.5 * (par.par_P_plus[i] * n_plus[j] + par.par_P_plus[j] * n_plus[i] + np_Pp * n_plus[i] * n_plus[j]) / r2_plus +
                  1.5 * (par.par_P_minus[i] * n_minus[j] + par.par_P_minus[j] * n_minus[i] + nm_Pm * n_minus[i] * n_minus[j]) / r2_minus -
                  3.0 * (np_Sp[i] * n_plus[j] + np_Sp[j] * n_plus[i]) / r3_plus - 3.0 * (nm_Sm[i] * n_minus[j] + nm_Sm[j] * n_minus[i]) / r3_minus;
      if (i == j)
        Aij[i][j] -= +1.5 * (np_Pp / r2_plus + nm_Pm / r2_minus);
    }
  }
}

/*-----------------------------------------------------------*/
/********           Nonlinear Equations            ***********/
/*-----------------------------------------------------------*/
void NonLinEquations(ID_persist_struct par, REAL rho_adm, REAL A, REAL B, REAL X, REAL R, REAL x, REAL r, REAL phi, REAL y, REAL z, derivs U, REAL *values) {
  REAL r_plus, r_minus, psi, psi2, psi4, psi7;

  r_plus = sqrt((x - par.par_b) * (x - par.par_b) + y * y + z * z);
  r_minus = sqrt((x + par.par_b) * (x + par.par_b) + y * y + z * z);

  psi = 1. + 0.5 * par.par_m_plus / r_plus + 0.5 * par.par_m_minus / r_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi7 = psi * psi2 * psi4;

  values[0] = U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * BY_KKofxyz(par, x, y, z) / psi7 + 2.0 * Pi / psi2 / psi * rho_adm;
}

/*-----------------------------------------------------------*/
/********               Linear Equations           ***********/
/*-----------------------------------------------------------*/
void LinEquations(ID_persist_struct par, REAL A, REAL B, REAL X, REAL R, REAL x, REAL r, REAL phi, REAL y, REAL z, derivs dU, derivs U, REAL *values) {
  REAL r_plus, r_minus, psi, psi2, psi4, psi8;

  r_plus = sqrt((x - par.par_b) * (x - par.par_b) + y * y + z * z);
  r_minus = sqrt((x + par.par_b) * (x + par.par_b) + y * y + z * z);

  psi = 1. + 0.5 * par.par_m_plus / r_plus + 0.5 * par.par_m_minus / r_minus + U.d0[0];
  psi2 = psi * psi;
  psi4 = psi2 * psi2;
  psi8 = psi4 * psi4;

  values[0] = dU.d11[0] + dU.d22[0] + dU.d33[0] - 0.875 * BY_KKofxyz(par, x, y, z) / psi8 * dU.d0[0];
}
"""

    name = "TP_Equations"
    params = """ """
    body = r"""
  // This space intentionally left blank:
  // NRPy+ requires C files share the same name as the primary function within the C file.
  // As TwoPunctures was developed outside the NRPy+ infrastructure, it does not abide by this rule.
  // Hence the need for a function with the same name as the output file.
"""

    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )

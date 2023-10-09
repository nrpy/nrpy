"""
Register TwoPunctures code Newton.c.

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


def register_CFunction_TP_Newton() -> None:
    """
    Register C function TP_Newton().

    This function provides Newton-Raphson algorithm needed for adjusting bare masses to achieve target BH masses.
    """
    desc = "Newton.c from TwoPunctures. Note that this C file is a collection of functions for Newton-Raphson method, with TP_Newton() being housed here."
    # prefunc contains most of the source code
    prefunc = f"// {desc}\n\n"
    prefunc += r"""/* TwoPunctures:  File  "Newton.c"*/

#include "TP_utilities.h"
#include "TwoPunctures.h"
#include <ctype.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int bicgstab(ID_persist_struct par, int const nvar, int const n1, int const n2, int const n3, derivs v, derivs dv, int const output, int const itmax, REAL const tol,
                    REAL *restrict const normres);
static REAL norm_inf(REAL const *restrict const F, int const ntotal);
static void relax(REAL *restrict const dv, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols, int **cols, REAL **JFD);
/*
  static void
  resid (REAL * restrict const res,
  int const ntotal,
  REAL const * restrict const dv,
  REAL const * restrict const rhs,
  int const * restrict const ncols,
  int **cols,
  REAL **JFD);
*/
static void LineRelax_al(REAL *restrict const dv, int const j, int const k, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols,
                         int **cols, REAL **JFD);
static void LineRelax_be(REAL *restrict const dv, int const i, int const k, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols,
                         int **cols, REAL **JFD);
/* --------------------------------------------------------------------------*/
static REAL norm_inf(REAL const *restrict const F, int const ntotal) {
  REAL dmax = -1;
#pragma omp parallel
  {
    REAL dmax1 = -1;
#pragma omp for
    for (int j = 0; j < ntotal; j++)
      if (fabs(F[j]) > dmax1)
        dmax1 = fabs(F[j]);
#pragma omp critical
    if (dmax1 > dmax)
      dmax = dmax1;
  }
  return dmax;
}
/* --------------------------------------------------------------------------*/
/* UNUSED
   static void
   resid (REAL * restrict const res,
   int const ntotal,
   REAL const * restrict const dv,
   REAL const * restrict const rhs,
   int const * restrict const ncols,
   int **cols,
   REAL **JFD)
   {
   #pragma omp parallel for
   for (int i = 0; i < ntotal; i++)
   {
   REAL JFDdv_i = 0;
   for (int m = 0; m < ncols[i]; m++)
   JFDdv_i += JFD[i][m] * dv[cols[i][m]];
   res[i] = rhs[i] - JFDdv_i;
   }
   }
*/
/* -------------------------------------------------------------------------*/
static void LineRelax_al(REAL *restrict const dv, int const j, int const k, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols,
                         int **cols, REAL **JFD) {
  int i, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n1);
  gsl_vector *e = gsl_vector_alloc(n1 - 1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n1 - 1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n1);     /* rhs */
  gsl_vector *x = gsl_vector_alloc(n1);     /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++) {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (i = 0; i < n1; i++) {
      Ip = Index(ivar, i + 1, j, k, nvar, n1, n2, n3);
      Ic = Index(ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index(ivar, i - 1, j, k, nvar, n1, n2, n3);
      gsl_vector_set(b, i, rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++) {
        col = cols[Ic][m];
        if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, i) -= JFD[Ic][m] * dv[col];
        else {
          if (col == Im && i > 0)
            gsl_vector_set(f, i - 1, JFD[Ic][m]);
          if (col == Ic)
            gsl_vector_set(diag, i, JFD[Ic][m]);
          if (col == Ip && i < n1 - 1)
            gsl_vector_set(e, i, JFD[Ic][m]);
        }
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (i = 0; i < n1; i++) {
      Ic = Index(ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, i);
    }
  }

  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void LineRelax_be(REAL *restrict const dv, int const i, int const k, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols,
                         int **cols, REAL **JFD) {
  int j, m, Ic, Ip, Im, col, ivar;

  gsl_vector *diag = gsl_vector_alloc(n2);
  gsl_vector *e = gsl_vector_alloc(n2 - 1); /* above diagonal */
  gsl_vector *f = gsl_vector_alloc(n2 - 1); /* below diagonal */
  gsl_vector *b = gsl_vector_alloc(n2);     /* rhs */
  gsl_vector *x = gsl_vector_alloc(n2);     /* solution vector */

  for (ivar = 0; ivar < nvar; ivar++) {
    gsl_vector_set_zero(diag);
    gsl_vector_set_zero(e);
    gsl_vector_set_zero(f);
    for (j = 0; j < n2; j++) {
      Ip = Index(ivar, i, j + 1, k, nvar, n1, n2, n3);
      Ic = Index(ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index(ivar, i, j - 1, k, nvar, n1, n2, n3);
      gsl_vector_set(b, j, rhs[Ic]);
      for (m = 0; m < ncols[Ic]; m++) {
        col = cols[Ic][m];
        if (col != Ip && col != Ic && col != Im)
          *gsl_vector_ptr(b, j) -= JFD[Ic][m] * dv[col];
        else {
          if (col == Im && j > 0)
            gsl_vector_set(f, j - 1, JFD[Ic][m]);
          if (col == Ic)
            gsl_vector_set(diag, j, JFD[Ic][m]);
          if (col == Ip && j < n2 - 1)
            gsl_vector_set(e, j, JFD[Ic][m]);
        }
      }
    }
    gsl_linalg_solve_tridiag(diag, e, f, b, x);
    for (j = 0; j < n2; j++) {
      Ic = Index(ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = gsl_vector_get(x, j);
    }
  }
  gsl_vector_free(diag);
  gsl_vector_free(e);
  gsl_vector_free(f);
  gsl_vector_free(b);
  gsl_vector_free(x);
}

/* --------------------------------------------------------------------------*/
static void relax(REAL *restrict const dv, int const nvar, int const n1, int const n2, int const n3, REAL const *restrict const rhs, int const *restrict const ncols, int **cols, REAL **JFD) {
  int i, j, k, n;

  for (k = 0; k < n3; k = k + 2) {
    for (n = 0; n < N_PlaneRelax; n++) {
#pragma omp parallel for schedule(dynamic)
      for (i = 2; i < n1; i = i + 2)
        LineRelax_be(dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (i = 1; i < n1; i = i + 2)
        LineRelax_be(dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 1; j < n2; j = j + 2)
        LineRelax_al(dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 0; j < n2; j = j + 2)
        LineRelax_al(dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
  for (k = 1; k < n3; k = k + 2) {
    for (n = 0; n < N_PlaneRelax; n++) {
#pragma omp parallel for schedule(dynamic)
      for (i = 0; i < n1; i = i + 2)
        LineRelax_be(dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (i = 1; i < n1; i = i + 2)
        LineRelax_be(dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 1; j < n2; j = j + 2)
        LineRelax_al(dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (j = 0; j < n2; j = j + 2)
        LineRelax_al(dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
}

/* --------------------------------------------------------------------------*/
/* ZACH SAYS: THIS FUNCTION IS COMMENTED OUT WITHIN Newton(), SO WE COMMENT IT OUT HERE AS WELL.
   void
   TestRelax (ID_persist_struct par,
   int nvar, int n1, int n2, int n3, derivs v,
   REAL *dv)
   {
   int ntotal = n1 * n2 * n3 * nvar, **cols, *ncols, maxcol =
   StencilSize * nvar, j;
   REAL *F, *res, **JFD;
   derivs u;

   F = dvector (0, ntotal - 1);
   res = dvector (0, ntotal - 1);
   allocate_derivs (&u, ntotal);

   JFD = dmatrix (0, ntotal - 1, 0, maxcol - 1);
   cols = imatrix (0, ntotal - 1, 0, maxcol - 1);
   ncols = ivector (0, ntotal - 1);

   F_of_v (par, nvar, n1, n2, n3, v, F, u);

   SetMatrix_JFD (par, nvar, n1, n2, n3, u, ncols, cols, JFD);

   for (j = 0; j < ntotal; j++)
   dv[j] = 0;
   resid (res, ntotal, dv, F, ncols, cols, JFD);
   fprintf(stderr,"Before: |F|=%20.15e\n", (double) norm1 (res, ntotal));
   fflush(stdout);
   for (j = 0; j < NRELAX; j++)
   {
   relax (dv, nvar, n1, n2, n3, F, ncols, cols, JFD);	// solves JFD*sh = s
   if (j % Step_Relax == 0)
   {
   resid (res, ntotal, dv, F, ncols, cols, JFD);
   fprintf(stderr,"j=%d\t |F|=%20.15e\n", j, (double) norm1 (res, ntotal));
   fflush(stdout);
   }
   }

   resid (res, ntotal, dv, F, ncols, cols, JFD);
   fprintf(stderr,"After: |F|=%20.15e\n", (double) norm1 (res, ntotal));
   fflush(stdout);

   free_dvector (F, 0, ntotal - 1);
   free_dvector (res, 0, ntotal - 1);
   free_derivs (&u, ntotal);

   free_dmatrix (JFD, 0, ntotal - 1, 0, maxcol - 1);
   free_imatrix (cols, 0, ntotal - 1, 0, maxcol - 1);
   free_ivector (ncols, 0, ntotal - 1);
   }
*/

/* --------------------------------------------------------------------------*/
static int bicgstab(ID_persist_struct par, int const nvar, int const n1, int const n2, int const n3, derivs v, derivs dv, int const output, int const itmax, REAL const tol,
                    REAL *restrict const normres) {
  int ntotal = n1 * n2 * n3 * nvar, ii;
  REAL alpha = 0, beta = 0;
  REAL rho = 0, rho1 = 1, rhotol = 1e-50;
  REAL omega = 0, omegatol = 1e-50;
  REAL *p, *rt, *s, *t, *r, *vv;
  REAL **JFD;
  int **cols, *ncols, maxcol = StencilSize * nvar;
  REAL *F;
  derivs u, ph, sh;

  F = dvector(0, ntotal - 1);
  allocate_derivs(&u, ntotal);

  JFD = dmatrix(0, ntotal - 1, 0, maxcol - 1);
  cols = imatrix(0, ntotal - 1, 0, maxcol - 1);
  ncols = ivector(0, ntotal - 1);

  F_of_v(par, nvar, n1, n2, n3, v, F, u);
  SetMatrix_JFD(par, nvar, n1, n2, n3, u, ncols, cols, JFD);

  /* temporary storage */
  r = dvector(0, ntotal - 1);
  p = dvector(0, ntotal - 1);
  allocate_derivs(&ph, ntotal);
  /*      ph  = dvector(0, ntotal-1);*/
  rt = dvector(0, ntotal - 1);
  s = dvector(0, ntotal - 1);
  allocate_derivs(&sh, ntotal);
  /*      sh  = dvector(0, ntotal-1);*/
  t = dvector(0, ntotal - 1);
  vv = dvector(0, ntotal - 1);

  /* check */
  if (output == 1) {
    fprintf(stderr, "bicgstab:  itmax %d, tol %e\n", itmax, (double)tol);
    fflush(stdout);
  }

  /* compute initial residual rt = r = F - J*dv */
  J_times_dv(par, nvar, n1, n2, n3, dv, r, u);
#pragma omp parallel for
  for (int j = 0; j < ntotal; j++)
    rt[j] = r[j] = F[j] - r[j];

  *normres = norm2(r, ntotal);
  if (output == 1) {
    fprintf(stderr, "bicgstab: %5d  %10.3e\n", 0, (double)*normres);
    fflush(stdout);
  }

  if (*normres > tol) {
    /* cgs iteration */
    for (ii = 0; ii < itmax; ii++) {
      rho = scalarproduct(rt, r, ntotal);
      if (fabs(rho) < rhotol)
        break;

      /* compute direction vector p */
      if (ii == 0) {
#pragma omp parallel for
        for (int j = 0; j < ntotal; j++)
          p[j] = r[j];
      } else {
        beta = (rho / rho1) * (alpha / omega);
#pragma omp parallel for
        for (int j = 0; j < ntotal; j++)
          p[j] = r[j] + beta * (p[j] - omega * vv[j]);
      }

      /* compute direction adjusting vector ph and scalar alpha */
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
        ph.d0[j] = 0;
      for (int j = 0; j < NRELAX; j++) /* solves JFD*ph = p by relaxation*/
        relax(ph.d0, nvar, n1, n2, n3, p, ncols, cols, JFD);

      J_times_dv(par, nvar, n1, n2, n3, ph, vv, u); /* vv=J*ph*/
      alpha = rho / scalarproduct(rt, vv, ntotal);
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
        s[j] = r[j] - alpha * vv[j];

      /* early check of tolerance */
      *normres = norm2(s, ntotal);
      if (*normres <= tol) {
#pragma omp parallel for
        for (int j = 0; j < ntotal; j++)
          dv.d0[j] += alpha * ph.d0[j];
        if (output == 1) {
          fprintf(stderr, "bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", ii + 1, (double)*normres, (double)alpha, (double)beta, (double)omega);
          fflush(stdout);
        }
        break;
      }

      /* compute stabilizer vector sh and scalar omega */
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++)
        sh.d0[j] = 0;
      for (int j = 0; j < NRELAX; j++) /* solves JFD*sh = s by relaxation*/
        relax(sh.d0, nvar, n1, n2, n3, s, ncols, cols, JFD);

      J_times_dv(par, nvar, n1, n2, n3, sh, t, u); /* t=J*sh*/
      omega = scalarproduct(t, s, ntotal) / scalarproduct(t, t, ntotal);

      /* compute new solution approximation */
#pragma omp parallel for
      for (int j = 0; j < ntotal; j++) {
        dv.d0[j] += alpha * ph.d0[j] + omega * sh.d0[j];
        r[j] = s[j] - omega * t[j];
      }
      /* are we done? */
      *normres = norm2(r, ntotal);
      if (output == 1) {
        fprintf(stderr, "bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n", ii + 1, (double)*normres, (double)alpha, (double)beta, (double)omega);
        fflush(stdout);
      }
      if (*normres <= tol)
        break;
      rho1 = rho;
      if (fabs(omega) < omegatol)
        break;
    }
  }
  /* free temporary storage */
  free_dvector(r, 0, ntotal - 1);
  free_dvector(p, 0, ntotal - 1);
  /*      free_dvector(ph,  0, ntotal-1);*/
  free_derivs(&ph, ntotal);
  free_dvector(rt, 0, ntotal - 1);
  free_dvector(s, 0, ntotal - 1);
  /*      free_dvector(sh,  0, ntotal-1);*/
  free_derivs(&sh, ntotal);
  free_dvector(t, 0, ntotal - 1);
  free_dvector(vv, 0, ntotal - 1);

  free_dvector(F, 0, ntotal - 1);
  free_derivs(&u, ntotal);

  free_dmatrix(JFD, 0, ntotal - 1, 0, maxcol - 1);
  free_imatrix(cols, 0, ntotal - 1, 0, maxcol - 1);
  free_ivector(ncols, 0, ntotal - 1);

  if (*normres <= tol)
    return 0;

  /* iteration failed */
  if (ii > itmax)
    return -1;

  /* breakdown */
  if (fabs(rho) < rhotol)
    return -10;
  if (fabs(omega) < omegatol)
    return -11;

  /* success! */
  return ii + 1;
}

/* -------------------------------------------------------------------*/
"""

    name = "TP_Newton"
    params = """ID_persist_struct par,
            int const nvar, int const n1, int const n2, int const n3,
            derivs v,
            REAL const tol, int const itmax"""

    body = r"""  int ntotal = n1 * n2 * n3 * nvar, it; // , ii
  REAL *F, dmax, normres;
  derivs u, dv;

  F = dvector(0, ntotal - 1);
  allocate_derivs(&dv, ntotal);
  allocate_derivs(&u, ntotal);

  /*         TestRelax(nvar, n1, n2, n3, v, dv.d0); */
  it = 0;
  dmax = 1;
  while (dmax > tol && it < itmax) {
    if (it == 0) {
      F_of_v(par, nvar, n1, n2, n3, v, F, u);
      dmax = norm_inf(F, ntotal);
    }
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      dv.d0[j] = 0;

    if (par.verbose == 1) {
      fprintf(stderr, "Newton: it=%d \t |F|=%e\n", it, (double)dmax);
      fprintf(stderr, "bare mass: mp=%g \t mm=%g\n", (double)par.par_m_plus, (double)par.par_m_minus);
    }

    fflush(stdout);
    // ii =
    bicgstab(par, nvar, n1, n2, n3, v, dv, par.verbose, 100, dmax * 1.e-3, &normres);
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      v.d0[j] -= dv.d0[j];
    F_of_v(par, nvar, n1, n2, n3, v, F, u);
    dmax = norm_inf(F, ntotal);
    it += 1;
  }
  if (itmax == 0) {
    F_of_v(par, nvar, n1, n2, n3, v, F, u);
    dmax = norm_inf(F, ntotal);
  }

  if (par.verbose == 1)
    fprintf(stderr, "Newton: it=%d \t |F|=%e \n", it, (double)dmax);

  fflush(stdout);

  free_dvector(F, 0, ntotal - 1);
  free_derivs(&dv, ntotal);
  free_derivs(&u, ntotal);
"""

    cfc.register_CFunction(
        subdirectory="TwoPunctures",
        prefunc=prefunc,
        desc=desc,
        name=name,
        params=params,
        body=body,
    )

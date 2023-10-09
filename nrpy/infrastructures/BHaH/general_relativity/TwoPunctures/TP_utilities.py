"""
Register TwoPunctures code TP_utilities.c.

Note that this C file is just a collection of functions for
memory allocation/deallocation, vector operations, etc.
The registered function TP_utilities() is unused.

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


def register_CFunction_TP_utilities() -> None:
    """Register C function TP_utilities(), a collection of functions for memory allocation/deallocation, vector operations, etc."""
    desc = "TP_utilities.c from TwoPunctures. Note that this C file is just a collection of functions for memory allocation/deallocation, vector operations, etc. TP_utilities() is unused."
    # prefunc contains most of the source code
    prefunc = "// " + desc + "\n\n"
    prefunc += r"""/* TwoPunctures:  File  "TP_utilities.c"*/

#include "TP_utilities.h"
#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

/*---------------------------------------------------------------------------*/
int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *retval;

  retval = malloc(sizeof(int) * (nh - nl + 1));
  if (retval == NULL) {
    fprintf(stderr, "allocation failure in ivector()");
    exit(1);
  }

  return retval - nl;
}

/*---------------------------------------------------------------------------*/
REAL *dvector(long nl, long nh)
/* allocate a REAL vector with subscript range v[nl..nh] */
{
  REAL *retval;

  retval = malloc(sizeof(REAL) * (nh - nl + 1));
  if (retval == NULL) {
    fprintf(stderr, "allocation failure in dvector()");
    exit(1);
  }

  return retval - nl;
}

/*---------------------------------------------------------------------------*/
int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  int **retval;

  retval = malloc(sizeof(int *) * (nrh - nrl + 1));
  if (retval == NULL) {
    fprintf(stderr, "allocation failure (1) in imatrix()");
    exit(1);
  }

  /* get all memory for the matrix in on chunk */
  retval[0] = malloc(sizeof(int) * (nrh - nrl + 1) * (nch - ncl + 1));
  if (retval[0] == NULL) {
    fprintf(stderr, "allocation failure (2) in imatrix()");
    exit(1);
  }

  /* apply column and row offsets */
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows */
  long width = (nch - ncl + 1);
  for (long i = nrl + 1; i <= nrh; i++)
    retval[i] = retval[i - 1] + width;
  assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);

  return retval;
}

/*---------------------------------------------------------------------------*/
REAL **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  REAL **retval;

  retval = malloc(sizeof(REAL *) * (nrh - nrl + 1));
  if (retval == NULL) {
    fprintf(stderr, "allocation failure (1) in dmatrix()");
    exit(1);
  }

  /* get all memory for the matrix in on chunk */
  retval[0] = malloc(sizeof(REAL) * (nrh - nrl + 1) * (nch - ncl + 1));
  if (retval[0] == NULL) {
    fprintf(stderr, "allocation failure (2) in dmatrix()");
    exit(1);
  }

  /* apply column and row offsets */
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows */
  long width = (nch - ncl + 1);
  for (long i = nrl + 1; i <= nrh; i++)
    retval[i] = retval[i - 1] + width;
  assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);

  return retval;
}

/*---------------------------------------------------------------------------*/
REAL ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a REAL 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  REAL ***retval;

  /* get memory for index structures */
  retval = malloc(sizeof(REAL **) * (nrh - nrl + 1));
  if (retval == NULL) {
    fprintf(stderr, "allocation failure (1) in d3tensor()");
    exit(1);
  }

  retval[0] = malloc(sizeof(REAL *) * (nrh - nrl + 1) * (nch - ncl + 1));
  if (retval[0] == NULL) {
    fprintf(stderr, "allocation failure (2) in d3tensor()");
    exit(1);
  }

  /* get all memory for the tensor in on chunk */
  retval[0][0] = malloc(sizeof(REAL) * (nrh - nrl + 1) * (nch - ncl + 1) * (ndh - ndl + 1));
  if (retval[0][0] == NULL) {
    fprintf(stderr, "allocation failure (3) in d3tensor()");
    exit(1);
  }

  /* apply all offsets */
  retval[0][0] -= ndl;
  retval[0] -= ncl;
  retval -= nrl;

  /* slice chunk into rows and columns */
  long width = (nch - ncl + 1);
  long depth = (ndh - ndl + 1);
  for (long j = ncl + 1; j <= nch; j++) { /* first row of columns */
    retval[nrl][j] = retval[nrl][j - 1] + depth;
  }
  assert(retval[nrl][nch] - retval[nrl][ncl] == (nch - ncl) * depth);
  for (long i = nrl + 1; i <= nrh; i++) {
    retval[i] = retval[i - 1] + width;
    retval[i][ncl] = retval[i - 1][ncl] + width * depth; /* first cell in column */
    for (long j = ncl + 1; j <= nch; j++) {
      retval[i][j] = retval[i][j - 1] + depth;
    }
    assert(retval[i][nch] - retval[i][ncl] == (nch - ncl) * depth);
  }
  assert(retval[nrh] - retval[nrl] == (nrh - nrl) * width);
  assert(&retval[nrh][nch][ndh] - &retval[nrl][ncl][ndl] == (nrh - nrl + 1) * (nch - ncl + 1) * (ndh - ndl + 1) - 1);

  return retval;
}

/*--------------------------------------------------------------------------*/
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free(v + nl);
}

/*--------------------------------------------------------------------------*/
void free_dvector(REAL *v, long nl, long nh)
/* free an double vector allocated with dvector() */
{
  free(v + nl);
}

/*--------------------------------------------------------------------------*/
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free(m[nrl] + ncl);
  free(m + nrl);
}

/*--------------------------------------------------------------------------*/
void free_dmatrix(REAL **m, long nrl, long nrh, long ncl, long nch)
/* free a REAL matrix allocated by dmatrix() */
{
  free(m[nrl] + ncl);
  free(m + nrl);
}

/*--------------------------------------------------------------------------*/
void free_d3tensor(REAL ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* free a REAL d3tensor allocated by d3tensor() */
{
  free(t[nrl][ncl] + ndl);
  free(t[nrl] + ncl);
  free(t + nrl);
}

/*--------------------------------------------------------------------------*/
int minimum2(int i, int j) {
  int result = i;
  if (j < result)
    result = j;
  return result;
}

/*-------------------------------------------------------------------------*/
int minimum3(int i, int j, int k) {
  int result = i;
  if (j < result)
    result = j;
  if (k < result)
    result = k;
  return result;
}

/*--------------------------------------------------------------------------*/
int maximum2(int i, int j) {
  int result = i;
  if (j > result)
    result = j;
  return result;
}

/*--------------------------------------------------------------------------*/
int maximum3(int i, int j, int k) {
  int result = i;
  if (j > result)
    result = j;
  if (k > result)
    result = k;
  return result;
}

/*--------------------------------------------------------------------------*/
int pow_int(int mantisse, int exponent) {
  int i, result = 1;

  for (i = 1; i <= exponent; i++)
    result *= mantisse;

  return result;
}

/*--------------------------------------------------------------------------*/
void chebft_Zeros(REAL u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.4) of 2nd edition C++ NR */
{
  int k, j, isignum;
  REAL fac, sum, Pion, *c;

  c = dvector(0, n);
  Pion = Pi / n;
  if (inv == 0) {
    fac = 2.0 / n;
    isignum = 1;
    for (j = 0; j < n; j++) {
      sum = 0.0;
      for (k = 0; k < n; k++)
        sum += u[k] * cos(Pion * j * (k + 0.5));
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
  } else {
    for (j = 0; j < n; j++) {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++) {
        sum += u[k] * cos(Pion * (j + 0.5) * k) * isignum;
        isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
#if 0
    if (fabs(c[j]) < 5.e-16)
      u[j] = 0.0;
    else
#endif
    u[j] = c[j];
  free_dvector(c, 0, n);
}

/* --------------------------------------------------------------------------*/

void chebft_Extremes(REAL u[], int n, int inv)
/* eq. 5.8.7 and 5.8.8 at x = (5.8.5) of 2nd edition C++ NR */
{
  int k, j, isignum, N = n - 1;
  REAL fac, sum, PioN, *c;

  c = dvector(0, N);
  PioN = Pi / N;
  if (inv == 0) {
    fac = 2.0 / N;
    isignum = 1;
    for (j = 0; j < n; j++) {
      sum = 0.5 * (u[0] + u[N] * isignum);
      for (k = 1; k < N; k++)
        sum += u[k] * cos(PioN * j * k);
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
    c[N] = 0.5 * c[N];
  } else {
    for (j = 0; j < n; j++) {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++) {
        sum += u[k] * cos(PioN * j * k) * isignum;
        isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
    u[j] = c[j];
  free_dvector(c, 0, N);
}

/* --------------------------------------------------------------------------*/

void chder(REAL *c, REAL *cder, int n) {
  int j;

  cder[n] = 0.0;
  cder[n - 1] = 0.0;
  for (j = n - 2; j >= 0; j--)
    cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
}

/* --------------------------------------------------------------------------*/
REAL chebev(REAL a, REAL b, REAL c[], int m, REAL x)
/* eq. 5.8.11 of C++ NR (2nd ed) */
{
  int j;
  REAL djp2, djp1, dj; /* d_{j+2}, d_{j+1} and d_j */
  REAL y;

  /* rescale input to lie within [-1,1] */
  y = 2 * (x - 0.5 * (b + a)) / (b - a);

  dj = djp1 = 0;
  for (j = m - 1; j >= 1; j--) {
    /* advance the coefficients */
    djp2 = djp1;
    djp1 = dj;
    dj = 2 * y * djp1 - djp2 + c[j];
  }

  return y * dj - djp1 + 0.5 * c[0];
}

/* --------------------------------------------------------------------------*/
void fourft(REAL *u, int N, int inv)
/* a (slow) Fourier transform, seems to be just eq. 12.1.6 and 12.1.9 of C++ NR (2nd ed) */
{
  int l, k, iy, M;
  REAL x, x1, fac, Pi_fac, *a, *b;

  M = N / 2;
  a = dvector(0, M);
  b = dvector(1, M);   /* Actually: b=vector(1,M-1) but this is problematic if M=1*/
  a[0] = a[M] = 1e100; // <- Initialize to sidestep (bogus) compiler warning.
  fac = 1. / M;
  Pi_fac = Pi * fac;
  if (inv == 0) {
    for (l = 0; l <= M; l++) {
      a[l] = 0;
      if (l > 0 && l < M)
        b[l] = 0;
      x1 = Pi_fac * l;
      for (k = 0; k < N; k++) {
        x = x1 * k;
        a[l] += fac * u[k] * cos(x);
        if (l > 0 && l < M)
          b[l] += fac * u[k] * sin(x);
      }
    }
    u[0] = a[0];
    u[M] = a[M];
    for (l = 1; l < M; l++) {
      u[l] = a[l];
      u[l + M] = b[l];
    }
  } else {
    a[0] = u[0];
    a[M] = u[M];
    for (l = 1; l < M; l++) {
      a[l] = u[l];
      b[l] = u[M + l];
    }
    iy = 1;
    for (k = 0; k < N; k++) {
      u[k] = 0.5 * (a[0] + a[M] * iy);
      x1 = Pi_fac * k;
      for (l = 1; l < M; l++) {
        x = x1 * l;
        u[k] += a[l] * cos(x) + b[l] * sin(x);
      }
      iy = -iy;
    }
  }
  free_dvector(a, 0, M);
  free_dvector(b, 1, M);
}

/* -----------------------------------------*/
void fourder(REAL u[], REAL du[], int N) {
  int l, M, lpM;

  M = N / 2;
  du[0] = 0.;
  du[M] = 0.;
  for (l = 1; l < M; l++) {
    lpM = l + M;
    du[l] = u[lpM] * l;
    du[lpM] = -u[l] * l;
  }
}

/* -----------------------------------------*/
void fourder2(REAL u[], REAL d2u[], int N) {
  int l, l2, M, lpM;

  d2u[0] = 0.;
  M = N / 2;
  for (l = 1; l <= M; l++) {
    l2 = l * l;
    lpM = l + M;
    d2u[l] = -u[l] * l2;
    if (l < M)
      d2u[lpM] = -u[lpM] * l2;
  }
}

/* ----------------------------------------- */
REAL fourev(REAL *u, int N, REAL x) {
  int l, M = N / 2;
  REAL xl, result;

  result = 0.5 * (u[0] + u[M] * cos(x * M));
  for (l = 1; l < M; l++) {
    xl = x * l;
    result += u[l] * cos(xl) + u[M + l] * sin(xl);
  }
  return result;
}

/* ------------------------------------------------------------------------*/
REAL norm1(REAL *v, int n) {
  int i;
  REAL result = -1;

  for (i = 0; i < n; i++)
    if (fabs(v[i]) > result)
      result = fabs(v[i]);

  return result;
}

/* -------------------------------------------------------------------------*/
REAL norm2(REAL *v, int n) {
  int i;
  REAL result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * v[i];

  return sqrt(result);
}

/* -------------------------------------------------------------------------*/
REAL scalarproduct(REAL *v, REAL *w, int n) {
  int i;
  REAL result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * w[i];

  return result;
}

/* -------------------------------------------------------------------------*/
"""

    name = "TP_utilities"
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

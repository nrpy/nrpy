/*
TwoPunctures creates initial for two puncture black holes using a
single domain spectral method.  This method is described in
Marcus Ansorg, Bernd Brügmann, Wolfgang Tichy, "A single-domain
spectral method for black hole puncture data", PRD 70, 064011 (2004),
arXiv:gr-qc/0404056.

Code originally from Einstein Toolkit's EinsteinInitialData
arrangement:
https://bitbucket.org/einsteintoolkit/einsteininitialdata

License: Lesser GNU Public License, version 2.0+
*/

/* TwoPunctures:  File  "utilities.h"*/

#include <math.h>
#include <stdbool.h>
#define REAL double

#include "../BHaH_defines.h"
#include "TwoPunctures.h"

// #define REAL double
#define CCTK_EQUALS(A, B) (strncmp(A, B, 100) == 0)

#define Pi 3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164 /* Pi/2*/
#define Piq 0.78539816339744830961566084582 /* Pi/4*/

#define TINY 1.0e-20
#define SWAP(a, b)                                                                                                                                                                                     \
  {                                                                                                                                                                                                    \
    temp = (a);                                                                                                                                                                                        \
    (a) = (b);                                                                                                                                                                                         \
    (b) = temp;                                                                                                                                                                                        \
  }

#define nrerror TP_nrerror
#define ivector TP_ivector
#define dvector TP_dvector
#define imatrix TP_imatrix
#define dmatrix TP_dmatrix
#define d3tensor TP_d3tensor
#define free_ivector TP_free_ivector
#define free_dvector TP_free_dvector
#define free_imatrix TP_free_imatrix
#define free_dmatrix TP_free_dmatrix
#define free_d3tensor TP_free_d3tensor

void nrerror(char error_text[]);
int *ivector(long nl, long nh);
REAL *dvector(long nl, long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
REAL **dmatrix(long nrl, long nrh, long ncl, long nch);
REAL ***d3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
void free_ivector(int *v, long nl, long nh);
void free_dvector(REAL *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix(REAL **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor(REAL ***t, long nrl, long nrh, long ncl, long nch, long ndl, long ndh);

int minimum2(int i, int j);
int minimum3(int i, int j, int k);
int maximum2(int i, int j);
int maximum3(int i, int j, int k);
int pow_int(int mantisse, int exponent);

void chebft_Zeros(REAL u[], int n, int inv);
void chebft_Extremes(REAL u[], int n, int inv);
void chder(REAL *c, REAL *cder, int n);
REAL chebev(REAL a, REAL b, REAL c[], int m, REAL x);
void fourft(REAL *u, int N, int inv);
void fourder(REAL u[], REAL du[], int N);
void fourder2(REAL u[], REAL d2u[], int N);
REAL fourev(REAL *u, int N, REAL x);

REAL norm1(REAL *v, int n);
REAL norm2(REAL *v, int n);
REAL scalarproduct(REAL *v, REAL *w, int n);

/* Routines in  "TwoPunctures.c"*/
REAL TestSolution(REAL A, REAL B, REAL X, REAL R, REAL phi);
void TestVector_w(REAL *par, int nvar, int n1, int n2, int n3, REAL *w);

/* Routines in  "FuncAndJacobian.c"*/
int Index(int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs(derivs *v, int n);
void free_derivs(derivs *v, int n);
void Derivatives_AB3(int nvar, int n1, int n2, int n3, derivs v);
void F_of_v(ID_persist_struct par, int nvar, int n1, int n2, int n3, derivs v, REAL *F, derivs u);
void J_times_dv(ID_persist_struct par, int nvar, int n1, int n2, int n3, derivs dv, REAL *Jdv, derivs u);
void JFD_times_dv(ID_persist_struct par, int i, int j, int k, int nvar, int n1, int n2, int n3, derivs dv, derivs u, REAL *values);
void SetMatrix_JFD(ID_persist_struct par, int nvar, int n1, int n2, int n3, derivs u, int *ncols, int **cols, REAL **Matrix);
REAL PunctEvalAtArbitPosition(REAL *v, int ivar, REAL A, REAL B, REAL phi, int nvar, int n1, int n2, int n3);
void calculate_derivs(int i, int j, int k, int ivar, int nvar, int n1, int n2, int n3, derivs v, derivs vv);
REAL interpol(REAL a, REAL b, REAL c, derivs v);
REAL PunctTaylorExpandAtArbitPosition(ID_persist_struct par, int ivar, int nvar, int n1, int n2, int n3, derivs v, REAL x, REAL y, REAL z);
REAL PunctIntPolAtArbitPosition(ID_persist_struct par, int ivar, int nvar, int n1, int n2, int n3, derivs v, REAL x, REAL y, REAL z);
void SpecCoef(int n1, int n2, int n3, int ivar, REAL *v, REAL *cf);
REAL PunctEvalAtArbitPositionFast(REAL *v, int ivar, REAL A, REAL B, REAL phi, int nvar, int n1, int n2, int n3);
REAL PunctIntPolAtArbitPositionFast(ID_persist_struct par, int ivar, int nvar, int n1, int n2, int n3, derivs v, REAL x, REAL y, REAL z);

/* Routines in  "CoordTransf.c"*/
void AB_To_XR(int nvar, REAL A, REAL B, REAL *X, REAL *R, derivs U);
void C_To_c(ID_persist_struct par, int nvar, REAL X, REAL R, REAL *x, REAL *r, derivs U);
void rx3_To_xyz(int nvar, REAL x, REAL r, REAL phi, REAL *y, REAL *z, derivs U);

/* Routines in  "Equations.c"*/
REAL BY_KKofxyz(ID_persist_struct par, REAL x, REAL y, REAL z);
void BY_Aijofxyz(ID_persist_struct par, REAL x, REAL y, REAL z, REAL Aij[3][3]);
void NonLinEquations(ID_persist_struct par, REAL rho_adm, REAL A, REAL B, REAL X, REAL R, REAL x, REAL r, REAL phi, REAL y, REAL z, derivs U, REAL *values);
void LinEquations(ID_persist_struct par, REAL A, REAL B, REAL X, REAL R, REAL x, REAL r, REAL phi, REAL y, REAL z, derivs dU, derivs U, REAL *values);

/* Routines in  "Newton.c"*/
void TestRelax(ID_persist_struct par, int nvar, int n1, int n2, int n3, derivs v, REAL *dv);
void TP_Newton(ID_persist_struct par, int nvar, int n1, int n2, int n3, derivs v, REAL tol, int itmax);

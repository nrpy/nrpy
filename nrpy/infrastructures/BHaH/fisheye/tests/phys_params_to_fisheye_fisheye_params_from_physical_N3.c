#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifndef STANDALONE
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#else
typedef double REAL;
typedef struct __commondata_struct__ {
  REAL fisheye_phys_L;
  REAL fisheye_phys_r_trans1;
  REAL fisheye_phys_r_trans2;
  REAL fisheye_phys_r_trans3;
  REAL fisheye_phys_w_trans1;
  REAL fisheye_phys_w_trans2;
  REAL fisheye_phys_w_trans3;
} commondata_struct;
typedef struct __params_struct__ {
  REAL fisheye_a0;
  REAL fisheye_a1;
  REAL fisheye_a2;
  REAL fisheye_a3;
  REAL fisheye_R1;
  REAL fisheye_R2;
  REAL fisheye_R3;
  REAL fisheye_s1;
  REAL fisheye_s2;
  REAL fisheye_s3;
  REAL fisheye_c;
} params_struct;
#endif

#ifndef NRPY_PHYSICAL_FISHEYE_SOLVER_HELPERS_H
#define NRPY_PHYSICAL_FISHEYE_SOLVER_HELPERS_H
typedef double fisheye_solver_real;
#include <float.h>

/**
 * Compute log(cosh(x)) stably for large |x|.
 *
 * Uses log(cosh(x)) = |x| + log(1 + exp(-2|x|)) - log(2).
 *
 * @param[in] x Input value.
 * @return Stable log(cosh(x)).
 */
static inline fisheye_solver_real logcosh_stable(const fisheye_solver_real x) {
  const fisheye_solver_real ax = fabs(x);
  return ax + log1p(exp((fisheye_solver_real)-2.0 * ax)) - log((fisheye_solver_real)2.0);
} // END FUNCTION: logcosh_stable

/**
 * Evaluate the unscaled fisheye radial map at radius r.
 *
 * @param[in] r Radius at which to evaluate the map.
 * @param[in] a Fisheye plateau factors.
 * @param[in] R Fisheye transition centers.
 * @param[in] s Fisheye transition widths.
 * @param[in] N Number of fisheye transitions.
 * @return Unscaled fisheye radius.
 */
static inline fisheye_solver_real rbar_unscaled(const fisheye_solver_real r, const fisheye_solver_real a[], const fisheye_solver_real R[],
                                                const fisheye_solver_real s[], const int N) {
  fisheye_solver_real rb = a[N] * r;
  for (int i = 0; i < N; i++) {
    const fisheye_solver_real delta_a = a[i] - a[i + 1];
    const fisheye_solver_real denom = (fisheye_solver_real)2.0 * tanh(R[i] / s[i]);
    const fisheye_solver_real u = (r + R[i]) / s[i];
    const fisheye_solver_real v = (r - R[i]) / s[i];
    const fisheye_solver_real term = (delta_a * s[i]) / denom * (logcosh_stable(u) - logcosh_stable(v));
    rb += term;
  } // END LOOP: for i over fisheye transitions
  return rb;
} // END FUNCTION: rbar_unscaled

/**
 * Evaluate physical fisheye constraints and the current outer-boundary scale factor.
 *
 * @param[in] L Outer physical boundary radius.
 * @param[in] r_trans Physical transition centers.
 * @param[in] w_trans Physical transition widths.
 * @param[in] a Fisheye plateau factors.
 * @param[in] R Trial fisheye transition centers.
 * @param[in] s Trial fisheye transition widths.
 * @param[in] N Number of fisheye transitions.
 * @param[out] F Constraint residuals.
 * @param[out] c_out Optional scale-factor output.
 * @return 0 on success, nonzero on invalid or non-finite scaling.
 */
static inline int evaluate_constraints(const fisheye_solver_real L, const fisheye_solver_real r_trans[], const fisheye_solver_real w_trans[],
                                       const fisheye_solver_real a[], const fisheye_solver_real R[], const fisheye_solver_real s[], const int N,
                                       fisheye_solver_real F[], fisheye_solver_real *c_out) {
  const fisheye_solver_real rbar_L = rbar_unscaled(L, a, R, s, N);
  if (!(isfinite(rbar_L)) || !(rbar_L > (fisheye_solver_real)0.0))
    return 1;
  const fisheye_solver_real c = L / rbar_L;
  if (!(isfinite(c)) || !(c > (fisheye_solver_real)0.0))
    return 1;

  for (int i = 0; i < N; i++) {
    const fisheye_solver_real Rm = R[i] - s[i];
    const fisheye_solver_real Rp = R[i] + s[i];
    const fisheye_solver_real Rphys_R = c * rbar_unscaled(R[i], a, R, s, N);
    const fisheye_solver_real Rphys_Rp = c * rbar_unscaled(Rp, a, R, s, N);
    const fisheye_solver_real Rphys_Rm = c * rbar_unscaled(Rm, a, R, s, N);
    F[2 * i + 0] = Rphys_R - r_trans[i];
    F[2 * i + 1] = (Rphys_Rp - Rphys_Rm) - w_trans[i];
  } // END LOOP: for i over fisheye constraints

  if (c_out)
    *c_out = c;
  return 0;
} // END FUNCTION: evaluate_constraints

/**
 * Solve a dense linear system with Gaussian elimination and partial pivoting.
 *
 * @param[in] n Linear-system dimension.
 * @param[in] A_in Row-major n x n matrix.
 * @param[in] b_in Right-hand side vector.
 * @param[out] x_out Solution vector.
 * @return 0 on success, nonzero on invalid dimensions or singular matrix.
 */
static inline int solve_linear_system(const int n, const fisheye_solver_real *A_in, const fisheye_solver_real *b_in, fisheye_solver_real *x_out) {
  enum { MAX_LINEAR_N = 64 };
  if (n <= 0 || n > MAX_LINEAR_N) {
    fprintf(stderr, "ERROR: fisheye solve_linear_system requires 0 < n <= %d (got n=%d)\\n", MAX_LINEAR_N, n);
    return 1;
  } // END IF: invalid linear-system dimension
  fisheye_solver_real A[MAX_LINEAR_N][MAX_LINEAR_N + 1];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      A[i][j] = A_in[i * n + j];
    A[i][n] = b_in[i];
  } // END LOOP: for i over matrix rows

  for (int k = 0; k < n; k++) {
    int pivot = k;
    fisheye_solver_real maxval = fabs(A[k][k]);
    for (int i = k + 1; i < n; i++) {
      const fisheye_solver_real val = fabs(A[i][k]);
      if (val > maxval) {
        maxval = val;
        pivot = i;
      }
    }
    if (!(maxval > (fisheye_solver_real)0.0))
      return 1;
    if (pivot != k) {
      for (int j = k; j < n + 1; j++) {
        const fisheye_solver_real tmp = A[k][j];
        A[k][j] = A[pivot][j];
        A[pivot][j] = tmp;
      } // END LOOP: for j over pivot-row entries
    } // END IF: pivot row must be swapped
    const fisheye_solver_real inv_pivot = (fisheye_solver_real)1.0 / A[k][k];
    for (int j = k; j < n + 1; j++)
      A[k][j] *= inv_pivot;
    for (int i = 0; i < n; i++) {
      if (i == k)
        continue;
      const fisheye_solver_real factor = A[i][k];
      for (int j = k; j < n + 1; j++)
        A[i][j] -= factor * A[k][j];
    } // END LOOP: for i over elimination rows
  } // END LOOP: for k over Gaussian-elimination pivots
  for (int i = 0; i < n; i++)
    x_out[i] = A[i][n];
  return 0;
} // END FUNCTION: solve_linear_system

/**
 * Solve physical fisheye inputs for internal transition centers and widths.
 *
 * @param[in] L Outer physical boundary radius.
 * @param[in] a Fisheye plateau factors.
 * @param[in] r_trans Physical transition centers.
 * @param[in] w_trans Physical transition widths.
 * @param[in] NTRANS Number of fisheye transitions.
 * @param[in] max_iter Maximum Newton iterations.
 * @param[out] R_out Internal transition centers.
 * @param[out] s_out Internal transition widths.
 * @param[out] c_out Outer-boundary scale factor.
 * @return 0 on success, nonzero on invalid input, Newton failure, or invalid REAL outputs.
 */
static inline int solve_physical_fisheye_params(const fisheye_solver_real L, const fisheye_solver_real a[], const fisheye_solver_real r_trans[],
                                                const fisheye_solver_real w_trans[], const int NTRANS, const int max_iter,
                                                fisheye_solver_real R_out[], fisheye_solver_real s_out[], fisheye_solver_real *c_out) {
  const int NUNK = 2 * NTRANS;

  for (int i = 0; i < NTRANS + 1; i++) {
    if (!(a[i] > (fisheye_solver_real)0.0))
      return 1;
  } // END LOOP: for i over plateau factors
  if (!(L > (fisheye_solver_real)0.0))
    return 1;
  for (int i = 0; i < NTRANS; i++) {
    if (!(r_trans[i] > (fisheye_solver_real)0.0) || !(w_trans[i] > (fisheye_solver_real)0.0))
      return 1;
    if (i > 0 && !(r_trans[i] > r_trans[i - 1]))
      return 1;
  } // END LOOP: for i over transition inputs

  fisheye_solver_real approx_R[NTRANS];
  fisheye_solver_real approx_s[NTRANS];
  for (int i = 0; i < NTRANS; i++) {
    approx_R[i] = r_trans[i];
    approx_s[i] = (fisheye_solver_real)0.5 * w_trans[i];
  } // END LOOP: for i over initial approximations
  const fisheye_solver_real approx_rbar_L = rbar_unscaled(L, a, approx_R, approx_s, NTRANS);
  if (!(isfinite(approx_rbar_L)) || !(approx_rbar_L > (fisheye_solver_real)0.0))
    return 1;
  const fisheye_solver_real approx_c = L / approx_rbar_L;
  if (!(isfinite(approx_c)) || !(approx_c > (fisheye_solver_real)0.0))
    return 1;

  fisheye_solver_real x[NUNK];
  for (int i = 0; i < NTRANS; i++) {
    x[2 * i + 0] = r_trans[i] / approx_c;
    x[2 * i + 1] = ((fisheye_solver_real)0.5 * w_trans[i]) / approx_c;
  } // END LOOP: for i over initial guess

  const fisheye_solver_real relative_tol = (fisheye_solver_real)1e-10;
  const fisheye_solver_real eps = (fisheye_solver_real)1e-6;

  int converged = 0;
  fisheye_solver_real c = (fisheye_solver_real)1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    // Step 1: Build R and s arrays from x.
    fisheye_solver_real R[NTRANS];
    fisheye_solver_real s[NTRANS];
    for (int i = 0; i < NTRANS; i++) {
      R[i] = x[2 * i + 0];
      s[i] = x[2 * i + 1];
      if (!(R[i] > (fisheye_solver_real)0.0) || !(s[i] > (fisheye_solver_real)0.0))
        return 1;
      if (!(R[i] > s[i]))
        return 1;
    } // END LOOP: for i over Newton state

    fisheye_solver_real F[NUNK];
    if (evaluate_constraints(L, r_trans, w_trans, a, R, s, NTRANS, F, &c))
      return 1;

    fisheye_solver_real relative_Fnorm = (fisheye_solver_real)0.0;
    for (int i = 0; i < NTRANS; i++) {
      relative_Fnorm += fabs(F[2 * i + 0]) / r_trans[i];
      relative_Fnorm += fabs(F[2 * i + 1]) / w_trans[i];
    } // END LOOP: for i over normalized residuals
    if (relative_Fnorm < relative_tol) {
      converged = 1;
      break;
    }

    // Step 2: Compute the numerical Jacobian.
    fisheye_solver_real J[NUNK][NUNK];
    for (int j = 0; j < NUNK; j++) {
      const fisheye_solver_real xj = x[j];
      const fisheye_solver_real dx = eps * (fabs(xj) + (fisheye_solver_real)1.0);
      x[j] = xj + dx;

      fisheye_solver_real Rp[NTRANS];
      fisheye_solver_real sp[NTRANS];
      for (int i = 0; i < NTRANS; i++) {
        Rp[i] = x[2 * i + 0];
        sp[i] = x[2 * i + 1];
      } // END LOOP: for i over perturbed state

      fisheye_solver_real Fp[NUNK];
      if (evaluate_constraints(L, r_trans, w_trans, a, Rp, sp, NTRANS, Fp, NULL))
        return 1;
      for (int i = 0; i < NUNK; i++) {
        J[i][j] = (Fp[i] - F[i]) / dx;
      } // END LOOP: for i over Jacobian rows
      x[j] = xj;
    } // END LOOP: for j over Jacobian columns

    fisheye_solver_real delta[NUNK];
    fisheye_solver_real minusF[NUNK];
    for (int i = 0; i < NUNK; i++)
      minusF[i] = -F[i];

    if (solve_linear_system(NUNK, &J[0][0], minusF, delta))
      return 1;

    // Step 3: Apply a damped step to enforce positivity and improve robustness.
    fisheye_solver_real alpha = (fisheye_solver_real)1.0;
    int accepted = 0;
    for (int attempt = 0; attempt < 12; attempt++) {
      int ok = 1;
      for (int i = 0; i < NUNK; i++) {
        const fisheye_solver_real trial = x[i] + alpha * delta[i];
        if (!(trial > (fisheye_solver_real)0.0)) {
          ok = 0;
          break;
        }
      } // END LOOP: for i over trial values
      if (ok) {
        for (int i = 0; i < NTRANS; i++) {
          const fisheye_solver_real Rtrial = x[2 * i + 0] + alpha * delta[2 * i + 0];
          const fisheye_solver_real strial = x[2 * i + 1] + alpha * delta[2 * i + 1];
          if (!(Rtrial > strial)) {
            ok = 0;
            break;
          }
        } // END LOOP: for i over width checks
      } // END IF: trial values are positive
      if (ok) {
        fisheye_solver_real x_trial[NUNK];
        for (int i = 0; i < NUNK; i++)
          x_trial[i] = x[i] + alpha * delta[i];
        fisheye_solver_real Rt[NTRANS];
        fisheye_solver_real st[NTRANS];
        for (int i = 0; i < NTRANS; i++) {
          Rt[i] = x_trial[2 * i + 0];
          st[i] = x_trial[2 * i + 1];
        } // END LOOP: for i over trial state
        fisheye_solver_real F_trial[NUNK];
        if (evaluate_constraints(L, r_trans, w_trans, a, Rt, st, NTRANS, F_trial, NULL))
          return 1;
        fisheye_solver_real relative_Fnorm_trial = (fisheye_solver_real)0.0;
        for (int i = 0; i < NTRANS; i++) {
          relative_Fnorm_trial += fabs(F_trial[2 * i + 0]) / r_trans[i];
          relative_Fnorm_trial += fabs(F_trial[2 * i + 1]) / w_trans[i];
        } // END LOOP: for i over trial residuals
        if (relative_Fnorm_trial < relative_Fnorm) {
          for (int i = 0; i < NUNK; i++)
            x[i] = x_trial[i];
          accepted = 1;
          break;
        } // END IF: trial reduces residual norm
      } // END IF: trial state passes validation
      alpha *= (fisheye_solver_real)0.5;
    } // END LOOP: for attempt over damping attempts
    if (!accepted)
      return 1;
  } // END LOOP: for iter over Newton iterations

  if (!converged)
    return 1;

  // Validate every REAL conversion before writing any output.
  const long double output_max = sizeof(REAL) == sizeof(float) ? FLT_MAX : sizeof(REAL) == sizeof(double) ? DBL_MAX : LDBL_MAX;
  if (!isfinite(c) || !(c > 0) || c > output_max)
    return 1;
  if (!((REAL)c > 0))
    return 1;
  for (int i = 0; i < NTRANS; i++) {
    const fisheye_solver_real Ri = x[2 * i + 0];
    const fisheye_solver_real si = x[2 * i + 1];
    if (!isfinite(Ri) || !isfinite(si) || !(Ri > 0) || !(si > 0) || Ri > output_max || si > output_max)
      return 1;
    if (!((REAL)si > 0) || !((REAL)Ri > (REAL)si))
      return 1;
  } // END LOOP: validate stored transition parameters

  for (int i = 0; i < NTRANS; i++) {
    R_out[i] = x[2 * i + 0];
    s_out[i] = x[2 * i + 1];
  } // END LOOP: for i over solved transitions
  *c_out = c;
  return 0;
} // END FUNCTION: solve_physical_fisheye_params

#endif // NRPY_PHYSICAL_FISHEYE_SOLVER_HELPERS_H

/**
 * Compute fisheye internal parameters (R_i, s_i) and compute c from physical fisheye inputs
 * (r_trans_i, w_trans_i, L) for an N-transition fisheye.
 *
 * Physical/input parameter meanings:
 * - fisheye_phys_a{i}: user-specified fisheye relative plateau/stretch factor a_i.
 * - fisheye_phys_L: outer physical boundary radius.
 * - fisheye_phys_r_trans{i}: physical radius where transition i is centered.
 * - fisheye_phys_w_trans{i}: physical width of transition i.
 */
int fisheye_params_from_physical_N3(const commondata_struct *restrict commondata, params_struct *restrict params) {
  const int NTRANS = 3;

  const fisheye_solver_real L = commondata->fisheye_phys_L;
  const fisheye_solver_real a[4] = {params->fisheye_a0, params->fisheye_a1, params->fisheye_a2, params->fisheye_a3};
  const fisheye_solver_real r_trans[3] = {commondata->fisheye_phys_r_trans1, commondata->fisheye_phys_r_trans2, commondata->fisheye_phys_r_trans3};
  const fisheye_solver_real w_trans[3] = {commondata->fisheye_phys_w_trans1, commondata->fisheye_phys_w_trans2, commondata->fisheye_phys_w_trans3};

  const int max_iter = 80;
  fisheye_solver_real R[NTRANS];
  fisheye_solver_real s[NTRANS];
  fisheye_solver_real c = (fisheye_solver_real)1.0;
  if (solve_physical_fisheye_params(L, a, r_trans, w_trans, NTRANS, max_iter, R, s, &c))
    return 1;

  // Step 4: Commit results into params.
  for (int i = 0; i < NTRANS; i++) {
    const REAL Ri = (REAL)R[i];
    const REAL si = (REAL)s[i];
    switch (i) {

    case 0:
      params->fisheye_R1 = Ri;
      params->fisheye_s1 = si;
      break;

    case 1:
      params->fisheye_R2 = Ri;
      params->fisheye_s2 = si;
      break;

    case 2:
      params->fisheye_R3 = Ri;
      params->fisheye_s3 = si;
      break;

    } // END SWITCH: assign fisheye transition parameters
  } // END LOOP: for i over transition outputs
  params->fisheye_c = (REAL)c;
  return 0;
} // END FUNCTION: fisheye_params_from_physical_N3

#ifdef STANDALONE

#define FISHEYE_GRID_N 70

/**
 * Write a standalone two-dimensional fisheye grid visualization dataset.
 *
 * @param[in] fname Output text filename.
 * @return 0 on success, nonzero on invalid solver state or file I/O failure.
 */
static int write_fisheye_grid_txt(const char *fname) {
  const int N = FISHEYE_GRID_N;
  const REAL L = 1.400000000000000e+01;
  // Physical parameters:
  //   r_trans[i]: physical transition centers; w_trans[i]: physical widths.
  const fisheye_solver_real r_trans[] = {4.500000000000000e+00, 8.000000000000000e+00, 1.200000000000000e+01};
  const fisheye_solver_real w_trans[] = {1.000000000000000e+00, 1.500000000000000e+00, 1.500000000000000e+00};
  const fisheye_solver_real a[] = {1.000000000000000e+00, 2.000000000000000e+00, 4.000000000000000e+00, 8.000000000000000e+00};
  const int NTRANS = (int)(sizeof(r_trans) / sizeof(r_trans[0]));
  const int max_iter = 60;
  fisheye_solver_real R[NTRANS];
  fisheye_solver_real s[NTRANS];
  fisheye_solver_real c = (fisheye_solver_real)1.0;
  if (solve_physical_fisheye_params(L, a, r_trans, w_trans, NTRANS, max_iter, R, s, &c))
    return 1;

  const REAL dx = (REAL)2.0 * L / (REAL)(N - 1);

  FILE *fp = fopen(fname, "w");
  if (!fp)
    return 1;
  for (int i = 0; i < N; i++) {
    const REAL X = -L + dx * (REAL)i;
    for (int j = 0; j < N; j++) {
      const REAL Y = -L + dx * (REAL)j;
      const REAL r = sqrt(X * X + Y * Y);
      REAL Xb = 0.0, Yb = 0.0;
      if (r > (REAL)0.0) {
        const REAL Rphys = c * rbar_unscaled(r, a, R, s, NTRANS);
        const REAL scale = Rphys / r;
        Xb = X * scale;
        Yb = Y * scale;
      } // END IF: nonzero radius mapping
      fprintf(fp, "%d %d %.15e %.15e\n", i, j, Xb, Yb);
    } // END LOOP: for j over grid y
  } // END LOOP: for i over grid x
  fclose(fp);
  return 0;
} // END FUNCTION: write_fisheye_grid_txt

/**
 * Run the standalone fisheye-grid writer and print a plotting command.
 *
 * @return EXIT_SUCCESS on success, EXIT_FAILURE otherwise.
 */
int main(void) {
  if (write_fisheye_grid_txt("fisheye_grid.txt") != 0) {
    fprintf(stderr, "ERROR: write_fisheye_grid_txt failed.\n");
    return EXIT_FAILURE;
  } // END IF: standalone grid write failed
  printf("Wrote fisheye_grid.txt\n");
  printf("Plot command:\n");
  printf("python3 -c \"import numpy as np, matplotlib.pyplot as plt; d=np.loadtxt('fisheye_grid.txt'); "
         "N=int(d[:,0].max()+1); Xb=d[:,2].reshape(N,N); Yb=d[:,3].reshape(N,N); "
         "[plt.plot(Xb[i,:],Yb[i,:],c='k',lw=0.45) for i in range(N)]; "
         "[plt.plot(Xb[:,j],Yb[:,j],c='k',lw=0.45) for j in range(N)]; "
         "r_trans=[4.500000000000000e+00, 8.000000000000000e+00, 1.200000000000000e+01]; "
         "theta=np.linspace(0,2*np.pi,400); "
         "[plt.plot(r*np.cos(theta), r*np.sin(theta), c='k', lw=2.5, ls=(0,(2,4))) for r in r_trans]; "
         "plt.gca().set_aspect('equal','box'); plt.axis('off'); plt.show()\"\n");
  return EXIT_SUCCESS;
} // END FUNCTION: main

#endif // STANDALONE

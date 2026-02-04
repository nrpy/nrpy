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
  REAL fisheye_phys_w_trans1;
  REAL fisheye_phys_w_trans2;
  REAL fisheye_a0;
  REAL fisheye_a1;
  REAL fisheye_a2;
  REAL fisheye_R1;
  REAL fisheye_R2;
  REAL fisheye_s1;
  REAL fisheye_s2;
  REAL fisheye_c;
} commondata_struct;
#endif

static inline REAL rbar_unscaled(const REAL r, const REAL a[], const REAL R[], const REAL s[], const int N) {
  REAL rb = a[N] * r;
  for (int i = 0; i < N; i++) {
    const REAL delta_a = a[i] - a[i + 1];
    const REAL denom = (REAL)2.0 * tanh(R[i] / s[i]);
    const REAL term = (delta_a * s[i]) / denom * log(cosh((r + R[i]) / s[i]) / cosh((r - R[i]) / s[i]));
    rb += term;
  }
  return rb;
}

static inline void evaluate_constraints(const REAL L, const REAL r_trans[], const REAL w_trans[], const REAL a[], const REAL R[], const REAL s[],
                                        const int N, REAL F[], REAL *c_out) {
  const REAL rbar_L = rbar_unscaled(L, a, R, s, N);
  const REAL c = L / rbar_L;

  for (int i = 0; i < N; i++) {
    const REAL Rm = R[i] - s[i];
    const REAL Rp = R[i] + s[i];
    const REAL Rphys_R = c * rbar_unscaled(R[i], a, R, s, N);
    const REAL Rphys_Rp = c * rbar_unscaled(Rp, a, R, s, N);
    const REAL Rphys_Rm = c * rbar_unscaled(Rm, a, R, s, N);
    F[2 * i + 0] = Rphys_R - r_trans[i];
    F[2 * i + 1] = (Rphys_Rp - Rphys_Rm) - w_trans[i];
  }

  if (c_out)
    *c_out = c;
}

static inline int solve_linear_system(const int n, const REAL *A_in, const REAL *b_in, REAL *x_out) {
  // Gaussian elimination with partial pivoting; A_in is row-major n x n.
  REAL A[n][n + 1];
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++)
      A[i][j] = A_in[i * n + j];
    A[i][n] = b_in[i];
  }

  for (int k = 0; k < n; k++) {
    int pivot = k;
    REAL maxval = fabs(A[k][k]);
    for (int i = k + 1; i < n; i++) {
      const REAL val = fabs(A[i][k]);
      if (val > maxval) {
        maxval = val;
        pivot = i;
      }
    }
    if (!(maxval > (REAL)0.0))
      return 1;
    if (pivot != k) {
      for (int j = k; j < n + 1; j++) {
        const REAL tmp = A[k][j];
        A[k][j] = A[pivot][j];
        A[pivot][j] = tmp;
      }
    }
    const REAL inv_pivot = (REAL)1.0 / A[k][k];
    for (int j = k; j < n + 1; j++)
      A[k][j] *= inv_pivot;
    for (int i = 0; i < n; i++) {
      if (i == k)
        continue;
      const REAL factor = A[i][k];
      for (int j = k; j < n + 1; j++)
        A[i][j] -= factor * A[k][j];
    }
  }
  for (int i = 0; i < n; i++)
    x_out[i] = A[i][n];
  return 0;
}

/**
 * Compute fisheye internal parameters (R_i, s_i, c) from physical fisheye inputs
 * (r_trans_i, w_trans_i, L) for an N-transition fisheye.
 *
 * Physical parameter meanings (physical radius):
 * - fisheye_phys_L: outer physical boundary radius.
 * - fisheye_phys_r_trans{i}: physical radius where transition i is centered.
 * - fisheye_phys_w_trans{i}: physical width of transition i.
 */
int fisheye_params_from_physical_N2(commondata_struct *restrict commondata) {
  const int NTRANS = 2;
  const int NUNK = 2 * NTRANS;

  const REAL L = commondata->fisheye_phys_L;
  const REAL a[3] = {commondata->fisheye_a0, commondata->fisheye_a1, commondata->fisheye_a2};
  const REAL r_trans[2] = {commondata->fisheye_phys_r_trans1, commondata->fisheye_phys_r_trans2};
  const REAL w_trans[2] = {commondata->fisheye_phys_w_trans1, commondata->fisheye_phys_w_trans2};

  if (!(L > (REAL)0.0))
    return 1;
  for (int i = 0; i < NTRANS; i++) {
    if (!(r_trans[i] > (REAL)0.0) || !(w_trans[i] > (REAL)0.0))
      return 1;
    if (i > 0 && !(r_trans[i] > r_trans[i - 1]))
      return 1;
  }

  REAL x[NUNK];
  for (int i = 0; i < NTRANS; i++) {
    x[2 * i + 0] = r_trans[i];
    x[2 * i + 1] = (REAL)0.5 * w_trans[i];
  }

  const REAL tol = (REAL)1e-12;
  const int max_iter = 80;
  const REAL eps = (REAL)1e-6;

  int converged = 0;
  REAL c = (REAL)1.0;

  for (int iter = 0; iter < max_iter; iter++) {
    // Build R and s arrays from x
    REAL R[NTRANS];
    REAL s[NTRANS];
    for (int i = 0; i < NTRANS; i++) {
      R[i] = x[2 * i + 0];
      s[i] = x[2 * i + 1];
      if (!(R[i] > (REAL)0.0) || !(s[i] > (REAL)0.0))
        return 1;
    }

    REAL F[NUNK];
    evaluate_constraints(L, r_trans, w_trans, a, R, s, NTRANS, F, &c);

    REAL Fnorm = (REAL)0.0;
    for (int i = 0; i < NUNK; i++)
      Fnorm += fabs(F[i]);
    if (Fnorm < tol) {
      converged = 1;
      break;
    }

    // Numerical Jacobian
    REAL J[NUNK][NUNK];
    for (int j = 0; j < NUNK; j++) {
      const REAL xj = x[j];
      const REAL dx = eps * (fabs(xj) + (REAL)1.0);
      x[j] = xj + dx;

      REAL Rp[NTRANS];
      REAL sp[NTRANS];
      for (int i = 0; i < NTRANS; i++) {
        Rp[i] = x[2 * i + 0];
        sp[i] = x[2 * i + 1];
      }

      REAL Fp[NUNK];
      evaluate_constraints(L, r_trans, w_trans, a, Rp, sp, NTRANS, Fp, NULL);
      for (int i = 0; i < NUNK; i++) {
        J[i][j] = (Fp[i] - F[i]) / dx;
      }
      x[j] = xj;
    }

    REAL delta[NUNK];
    REAL minusF[NUNK];
    for (int i = 0; i < NUNK; i++)
      minusF[i] = -F[i];

    if (solve_linear_system(NUNK, &J[0][0], minusF, delta))
      return 1;

    // Damped step to enforce positivity and improve robustness
    REAL alpha = (REAL)1.0;
    int accepted = 0;
    for (int attempt = 0; attempt < 12; attempt++) {
      int ok = 1;
      for (int i = 0; i < NUNK; i++) {
        const REAL trial = x[i] + alpha * delta[i];
        if (!(trial > (REAL)0.0)) {
          ok = 0;
          break;
        }
      }
      if (ok) {
        REAL x_trial[NUNK];
        for (int i = 0; i < NUNK; i++)
          x_trial[i] = x[i] + alpha * delta[i];
        REAL Rt[NTRANS];
        REAL st[NTRANS];
        for (int i = 0; i < NTRANS; i++) {
          Rt[i] = x_trial[2 * i + 0];
          st[i] = x_trial[2 * i + 1];
        }
        REAL F_trial[NUNK];
        evaluate_constraints(L, r_trans, w_trans, a, Rt, st, NTRANS, F_trial, NULL);
        REAL Fnorm_trial = (REAL)0.0;
        for (int i = 0; i < NUNK; i++)
          Fnorm_trial += fabs(F_trial[i]);
        if (Fnorm_trial < Fnorm) {
          for (int i = 0; i < NUNK; i++)
            x[i] = x_trial[i];
          accepted = 1;
          break;
        }
      }
      alpha *= (REAL)0.5;
    }
    if (!accepted)
      return 1;
  }

  if (!converged)
    return 1;

  // Commit results into commondata
  for (int i = 0; i < NTRANS; i++) {
    const REAL Ri = x[2 * i + 0];
    const REAL si = x[2 * i + 1];
    switch (i) {

    case 0:
      commondata->fisheye_R1 = Ri;
      commondata->fisheye_s1 = si;
      break;

    case 1:
      commondata->fisheye_R2 = Ri;
      commondata->fisheye_s2 = si;
      break;
    }
  }
  commondata->fisheye_c = c;
  return 0;
} // END FUNCTION fisheye_params_from_physical_N2

#ifdef STANDALONE

#define FISHEYE_GRID_N 70

static int write_fisheye_grid_txt(const char *fname) {
  const int N = FISHEYE_GRID_N;
  const REAL L = 1.000000000000000e+01;
  // Physical parameters:
  //   r_trans[i]: physical transition centers; w_trans[i]: physical widths.
  const REAL r_trans[] = {4.500000000000000e+00, 8.000000000000000e+00};
  const REAL w_trans[] = {1.000000000000000e+00, 1.500000000000000e+00};
  const REAL a[] = {1.000000000000000e+00, 2.000000000000000e+00, 4.000000000000000e+00};
  const int NTRANS = (int)(sizeof(r_trans) / sizeof(r_trans[0]));
  const int NUNK = 2 * NTRANS;

  REAL x[NUNK];
  for (int i = 0; i < NTRANS; i++) {
    x[2 * i + 0] = r_trans[i];
    x[2 * i + 1] = (REAL)0.5 * w_trans[i];
  }
  REAL c = (REAL)1.0;
  const REAL tol = (REAL)1e-12;
  const int max_iter = 60;

  int converged = 0;
  for (int iter = 0; iter < max_iter; iter++) {
    for (int i = 0; i < NUNK; i++) {
      if (!(x[i] > (REAL)0.0))
        return 1;
    }

    REAL R[NTRANS];
    REAL s[NTRANS];
    for (int i = 0; i < NTRANS; i++) {
      R[i] = x[2 * i + 0];
      s[i] = x[2 * i + 1];
    }
    REAL F[NUNK];
    evaluate_constraints(L, r_trans, w_trans, a, R, s, NTRANS, F, &c);

    REAL Fnorm = (REAL)0.0;
    for (int i = 0; i < NUNK; i++)
      Fnorm += fabs(F[i]);
    if (Fnorm < tol) {
      converged = 1;
      break;
    }

    REAL J[NUNK][NUNK];
    const REAL eps = (REAL)1e-6;
    for (int j = 0; j < NUNK; j++) {
      REAL xj = x[j];
      const REAL dx = eps * (fabs(xj) + (REAL)1.0);
      x[j] = xj + dx;

      REAL Rp[NTRANS];
      REAL sp[NTRANS];
      for (int i = 0; i < NTRANS; i++) {
        Rp[i] = x[2 * i + 0];
        sp[i] = x[2 * i + 1];
      }
      REAL Fp[NUNK];
      evaluate_constraints(L, r_trans, w_trans, a, Rp, sp, NTRANS, Fp, NULL);
      for (int i = 0; i < NUNK; i++) {
        J[i][j] = (Fp[i] - F[i]) / dx;
      }
      x[j] = xj;
    }

    REAL delta[NUNK];
    REAL minusF[NUNK];
    for (int i = 0; i < NUNK; i++)
      minusF[i] = -F[i];
    if (solve_linear_system(NUNK, &J[0][0], minusF, delta))
      return 1;

    REAL alpha = (REAL)1.0;
    int accepted = 0;
    for (int attempt = 0; attempt < 12; attempt++) {
      int ok = 1;
      REAL x_trial[NUNK];
      for (int i = 0; i < NUNK; i++) {
        x_trial[i] = x[i] + alpha * delta[i];
        if (!(x_trial[i] > (REAL)0.0))
          ok = 0;
      }
      if (ok) {
        REAL Rt[NTRANS];
        REAL st[NTRANS];
        for (int i = 0; i < NTRANS; i++) {
          Rt[i] = x_trial[2 * i + 0];
          st[i] = x_trial[2 * i + 1];
        }
        REAL F_trial[NUNK];
        evaluate_constraints(L, r_trans, w_trans, a, Rt, st, NTRANS, F_trial, NULL);
        REAL Fnorm_trial = (REAL)0.0;
        for (int i = 0; i < NUNK; i++)
          Fnorm_trial += fabs(F_trial[i]);
        if (Fnorm_trial < Fnorm) {
          for (int i = 0; i < NUNK; i++)
            x[i] = x_trial[i];
          accepted = 1;
          break;
        }
      }
      alpha *= (REAL)0.5;
    }
    if (!accepted)
      return 1;
  }
  if (!converged)
    return 1;

  REAL R[NTRANS];
  REAL s[NTRANS];
  for (int i = 0; i < NTRANS; i++) {
    R[i] = x[2 * i + 0];
    s[i] = x[2 * i + 1];
  }
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
      }
      fprintf(fp, "%d %d %.15e %.15e\n", i, j, Xb, Yb);
    }
  }
  fclose(fp);
  return 0;
}

int main(void) {
  if (write_fisheye_grid_txt("fisheye_grid.txt") != 0) {
    fprintf(stderr, "ERROR: write_fisheye_grid_txt failed.\n");
    return EXIT_FAILURE;
  }
  printf("Wrote fisheye_grid.txt\n");
  printf("Plot command:\n");
  printf("python3 -c \"import numpy as np, matplotlib.pyplot as plt; d=np.loadtxt('fisheye_grid.txt'); "
         "N=int(d[:,0].max()+1); Xb=d[:,2].reshape(N,N); Yb=d[:,3].reshape(N,N); "
         "[plt.plot(Xb[i,:],Yb[i,:],c='k',lw=0.45) for i in range(N)]; "
         "[plt.plot(Xb[:,j],Yb[:,j],c='k',lw=0.45) for j in range(N)]; "
         "r_trans=[4.500000000000000e+00, 8.000000000000000e+00]; "
         "theta=np.linspace(0,2*np.pi,400); "
         "[plt.plot(r*np.cos(theta), r*np.sin(theta), c='k', lw=2.5, ls=(0,(2,4))) for r in r_trans]; "
         "plt.gca().set_aspect('equal','box'); plt.axis('off'); plt.show()\"\n");
  return EXIT_SUCCESS;
}

#endif // STANDALONE

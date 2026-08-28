#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef int32_t akv_blas_int;

#define IDX(i, j, ld) ((i) + (j) * (ld))

static int akv_lsame(const char *c, char want) {
  char got = c == NULL ? '\0' : c[0];
  if (got >= 'a' && got <= 'z') {
    got = (char)(got - 'a' + 'A');
  }
  if (want >= 'a' && want <= 'z') {
    want = (char)(want - 'a' + 'A');
  }
  return got == want;
}

static double *akv_alloc(size_t n) {
  return (double *)calloc(n == 0 ? 1 : n, sizeof(double));
}

static double akv_sym_get(const double *a, int lda, int i, int j) {
  return i <= j ? a[IDX(i, j, lda)] : a[IDX(j, i, lda)];
}

static void akv_swap_double(double *a, double *b) {
  const double t = *a;
  *a = *b;
  *b = t;
}

void bah_akv_internal_dcopy(akv_blas_int *n, const double *x, akv_blas_int *incx,
                            double *y, akv_blas_int *incy) {
  for (akv_blas_int i = 0; i < *n; i++) {
    y[i * *incy] = x[i * *incx];
  }
}

void bah_akv_internal_daxpy(akv_blas_int *n, const double *alpha,
                            const double *x, akv_blas_int *incx, double *y,
                            akv_blas_int *incy) {
  for (akv_blas_int i = 0; i < *n; i++) {
    y[i * *incy] += *alpha * x[i * *incx];
  }
}

void bah_akv_internal_dscal(akv_blas_int *n, const double *alpha, double *x,
                            akv_blas_int *incx) {
  for (akv_blas_int i = 0; i < *n; i++) {
    x[i * *incx] *= *alpha;
  }
}

double bah_akv_internal_ddot(akv_blas_int *n, const double *x,
                             akv_blas_int *incx, double *y,
                             akv_blas_int *incy) {
  double sum = 0.0;
  for (akv_blas_int i = 0; i < *n; i++) {
    sum += x[i * *incx] * y[i * *incy];
  }
  return sum;
}

void bah_akv_internal_dgemv(const char *transa, akv_blas_int *m,
                            akv_blas_int *n, const double *alpha,
                            const double *a, akv_blas_int *lda,
                            const double *x, akv_blas_int *incx,
                            const double *beta, double *y,
                            akv_blas_int *incy) {
  const int trans = !akv_lsame(transa, 'N');
  const int rows = trans ? *n : *m;
  const int inner = trans ? *m : *n;
  for (int i = 0; i < rows; i++) {
    double sum = 0.0;
    for (int k = 0; k < inner; k++) {
      const double aval = trans ? a[IDX(k, i, *lda)] : a[IDX(i, k, *lda)];
      sum += aval * x[k * *incx];
    }
    y[i * *incy] = *alpha * sum + *beta * y[i * *incy];
  }
}

void bah_akv_internal_dgemm(const char *transa, const char *transb,
                            akv_blas_int *m, akv_blas_int *n,
                            akv_blas_int *k, const double *alpha,
                            const double *a, akv_blas_int *lda,
                            const double *b, akv_blas_int *ldb,
                            const double *beta, double *c,
                            akv_blas_int *ldc) {
  const int ta = !akv_lsame(transa, 'N');
  const int tb = !akv_lsame(transb, 'N');
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *m; i++) {
      double sum = 0.0;
      for (int l = 0; l < *k; l++) {
        const double av = ta ? a[IDX(l, i, *lda)] : a[IDX(i, l, *lda)];
        const double bv = tb ? b[IDX(j, l, *ldb)] : b[IDX(l, j, *ldb)];
        sum += av * bv;
      }
      c[IDX(i, j, *ldc)] = *alpha * sum + *beta * c[IDX(i, j, *ldc)];
    }
  }
}

void bah_akv_internal_dsymm(const char *side, const char *uplo,
                            akv_blas_int *m, akv_blas_int *n,
                            const double *alpha, const double *a,
                            akv_blas_int *lda, const double *b,
                            akv_blas_int *ldb, const double *beta, double *c,
                            akv_blas_int *ldc) {
  (void)uplo;
  const int left = akv_lsame(side, 'L');
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *m; i++) {
      double sum = 0.0;
      const int inner = left ? *m : *n;
      for (int k = 0; k < inner; k++) {
        const double av = left ? akv_sym_get(a, *lda, i, k)
                               : akv_sym_get(a, *lda, k, j);
        const double bv = left ? b[IDX(k, j, *ldb)] : b[IDX(i, k, *ldb)];
        sum += av * bv;
      }
      c[IDX(i, j, *ldc)] = *alpha * sum + *beta * c[IDX(i, j, *ldc)];
    }
  }
}

static double akv_tri_get(const double *a, int lda, int i, int j,
                          const char *uplo, const char *transa,
                          const char *diag) {
  int r = i;
  int c = j;
  if (!akv_lsame(transa, 'N')) {
    r = j;
    c = i;
  }
  if (r == c && akv_lsame(diag, 'U')) {
    return 1.0;
  }
  if (akv_lsame(uplo, 'U') && r > c) {
    return 0.0;
  }
  if (akv_lsame(uplo, 'L') && r < c) {
    return 0.0;
  }
  return a[IDX(r, c, lda)];
}

void bah_akv_internal_dtrmm(const char *side, const char *uplo,
                            const char *transa, const char *diag,
                            akv_blas_int *m, akv_blas_int *n,
                            const double *alpha, const double *a,
                            akv_blas_int *lda, double *b,
                            akv_blas_int *ldb) {
  double *out = akv_alloc((size_t)*m * (size_t)*n);
  const int left = akv_lsame(side, 'L');
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *m; i++) {
      double sum = 0.0;
      const int inner = left ? *m : *n;
      for (int k = 0; k < inner; k++) {
        sum += left ? akv_tri_get(a, *lda, i, k, uplo, transa, diag) *
                          b[IDX(k, j, *ldb)]
                    : b[IDX(i, k, *ldb)] *
                          akv_tri_get(a, *lda, k, j, uplo, transa, diag);
      }
      out[IDX(i, j, *m)] = *alpha * sum;
    }
  }
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *m; i++) {
      b[IDX(i, j, *ldb)] = out[IDX(i, j, *m)];
    }
  }
  free(out);
}

static int akv_solve_dense(int n, double *a, double *b) {
  for (int k = 0; k < n; k++) {
    int piv = k;
    double pivabs = fabs(a[IDX(k, k, n)]);
    for (int i = k + 1; i < n; i++) {
      const double v = fabs(a[IDX(i, k, n)]);
      if (v > pivabs) {
        piv = i;
        pivabs = v;
      }
    }
    if (pivabs == 0.0) {
      return k + 1;
    }
    if (piv != k) {
      for (int j = k; j < n; j++) {
        akv_swap_double(&a[IDX(k, j, n)], &a[IDX(piv, j, n)]);
      }
      akv_swap_double(&b[k], &b[piv]);
    }
    for (int i = k + 1; i < n; i++) {
      const double f = a[IDX(i, k, n)] / a[IDX(k, k, n)];
      a[IDX(i, k, n)] = 0.0;
      for (int j = k + 1; j < n; j++) {
        a[IDX(i, j, n)] -= f * a[IDX(k, j, n)];
      }
      b[i] -= f * b[k];
    }
  }
  for (int i = n - 1; i >= 0; i--) {
    double sum = b[i];
    for (int j = i + 1; j < n; j++) {
      sum -= a[IDX(i, j, n)] * b[j];
    }
    b[i] = sum / a[IDX(i, i, n)];
  }
  return 0;
}

void bah_akv_internal_dtrsm(const char *side, const char *uplo,
                            const char *transa, const char *diag,
                            akv_blas_int *m, akv_blas_int *n, double *alpha,
                            double *a, akv_blas_int *lda, double *b,
                            akv_blas_int *ldb) {
  const int left = akv_lsame(side, 'L');
  const int dim = left ? *m : *n;
  double *amat = akv_alloc((size_t)dim * (size_t)dim);
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < dim; i++) {
      amat[IDX(i, j, dim)] =
          left ? akv_tri_get(a, *lda, i, j, uplo, transa, diag)
               : akv_tri_get(a, *lda, j, i, uplo, transa, diag);
    }
  }
  const int nrhs = left ? *n : *m;
  double *rhs = akv_alloc((size_t)dim);
  for (int r = 0; r < nrhs; r++) {
    for (int i = 0; i < dim; i++) {
      rhs[i] = *alpha * (left ? b[IDX(i, r, *ldb)] : b[IDX(r, i, *ldb)]);
    }
    double *acopy = akv_alloc((size_t)dim * (size_t)dim);
    memcpy(acopy, amat, sizeof(double) * (size_t)dim * (size_t)dim);
    (void)akv_solve_dense(dim, acopy, rhs);
    free(acopy);
    for (int i = 0; i < dim; i++) {
      if (left) {
        b[IDX(i, r, *ldb)] = rhs[i];
      } else {
        b[IDX(r, i, *ldb)] = rhs[i];
      }
    }
  }
  free(rhs);
  free(amat);
}

void bah_akv_internal_dpotrf(const char *uplo, akv_blas_int *n, double *a,
                             akv_blas_int *lda, akv_blas_int *info) {
  double *full = akv_alloc((size_t)*n * (size_t)*n);
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      full[IDX(i, j, *n)] = akv_sym_get(a, *lda, i, j);
    }
  }
  *info = 0;
  for (int j = 0; j < *n; j++) {
    double diag = full[IDX(j, j, *n)];
    for (int k = 0; k < j; k++) {
      diag -= full[IDX(j, k, *n)] * full[IDX(j, k, *n)];
    }
    if (diag <= 0.0) {
      *info = j + 1;
      free(full);
      return;
    }
    full[IDX(j, j, *n)] = sqrt(diag);
    for (int i = j + 1; i < *n; i++) {
      double sum = full[IDX(i, j, *n)];
      for (int k = 0; k < j; k++) {
        sum -= full[IDX(i, k, *n)] * full[IDX(j, k, *n)];
      }
      full[IDX(i, j, *n)] = sum / full[IDX(j, j, *n)];
    }
  }
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      if (akv_lsame(uplo, 'L')) {
        a[IDX(i, j, *lda)] = i >= j ? full[IDX(i, j, *n)] : 0.0;
      } else {
        a[IDX(i, j, *lda)] = i <= j ? full[IDX(j, i, *n)] : 0.0;
      }
    }
  }
  free(full);
}

void bah_akv_internal_dgetrf(akv_blas_int *m, akv_blas_int *n, double *a,
                             akv_blas_int *lda, akv_blas_int *ipiv,
                             akv_blas_int *info) {
  *info = 0;
  const int mn = *m < *n ? *m : *n;
  for (int k = 0; k < mn; k++) {
    int piv = k;
    double pivabs = fabs(a[IDX(k, k, *lda)]);
    for (int i = k + 1; i < *m; i++) {
      const double v = fabs(a[IDX(i, k, *lda)]);
      if (v > pivabs) {
        piv = i;
        pivabs = v;
      }
    }
    ipiv[k] = piv + 1;
    if (pivabs == 0.0 && *info == 0) {
      *info = k + 1;
    }
    if (piv != k) {
      for (int j = 0; j < *n; j++) {
        akv_swap_double(&a[IDX(k, j, *lda)], &a[IDX(piv, j, *lda)]);
      }
    }
    if (k < *m && a[IDX(k, k, *lda)] != 0.0) {
      for (int i = k + 1; i < *m; i++) {
        a[IDX(i, k, *lda)] /= a[IDX(k, k, *lda)];
        for (int j = k + 1; j < *n; j++) {
          a[IDX(i, j, *lda)] -= a[IDX(i, k, *lda)] * a[IDX(k, j, *lda)];
        }
      }
    }
  }
}

void bah_akv_internal_dgetrs(const char *trans, akv_blas_int *n,
                             akv_blas_int *nrhs, double *a, akv_blas_int *lda,
                             akv_blas_int *ipiv, double *b,
                             akv_blas_int *ldb, akv_blas_int *info) {
  *info = 0;
  if (!akv_lsame(trans, 'N')) {
    *info = -1;
    return;
  }
  for (int r = 0; r < *nrhs; r++) {
    for (int k = 0; k < *n; k++) {
      const int piv = ipiv[k] - 1;
      if (piv != k) {
        akv_swap_double(&b[IDX(k, r, *ldb)], &b[IDX(piv, r, *ldb)]);
      }
    }
    for (int i = 1; i < *n; i++) {
      double sum = b[IDX(i, r, *ldb)];
      for (int k = 0; k < i; k++) {
        sum -= a[IDX(i, k, *lda)] * b[IDX(k, r, *ldb)];
      }
      b[IDX(i, r, *ldb)] = sum;
    }
    for (int i = *n - 1; i >= 0; i--) {
      double sum = b[IDX(i, r, *ldb)];
      for (int k = i + 1; k < *n; k++) {
        sum -= a[IDX(i, k, *lda)] * b[IDX(k, r, *ldb)];
      }
      b[IDX(i, r, *ldb)] = sum / a[IDX(i, i, *lda)];
    }
  }
}

void bah_akv_internal_dsytrf(const char *uplo, akv_blas_int *n, double *a,
                             akv_blas_int *lda, akv_blas_int *ipiv,
                             double *work, akv_blas_int *lwork,
                             akv_blas_int *info) {
  (void)uplo;
  if (*lwork == -1) {
    work[0] = (double)(*n > 1 ? *n : 1);
    *info = 0;
    return;
  }
  akv_blas_int m = *n;
  bah_akv_internal_dgetrf(&m, n, a, lda, ipiv, info);
}

void bah_akv_internal_dsytrs(const char *uplo, akv_blas_int *n,
                             akv_blas_int *nrhs, double *a, akv_blas_int *lda,
                             akv_blas_int *ipiv, double *b,
                             akv_blas_int *ldb, akv_blas_int *info) {
  (void)uplo;
  bah_akv_internal_dgetrs("N", n, nrhs, a, lda, ipiv, b, ldb, info);
}

static uint32_t akv_seed_to_uint(const akv_blas_int *iseed) {
  uint32_t s = 0x6d2b79f5u;
  for (int i = 0; i < 4; i++) {
    s = s * 1664525u + (uint32_t)iseed[i] + 1013904223u;
  }
  return s == 0 ? 1u : s;
}

void bah_akv_internal_dlarnv(akv_blas_int *idist, akv_blas_int *iseed,
                             akv_blas_int *n, double *x) {
  uint32_t s = akv_seed_to_uint(iseed);
  for (int i = 0; i < *n; i++) {
    s = s * 1664525u + 1013904223u;
    double u = ((double)(s >> 8) + 1.0) / 16777217.0;
    if (*idist == 2) {
      u = 2.0 * u - 1.0;
    } else if (*idist == 3) {
      s = s * 1664525u + 1013904223u;
      const double v = ((double)(s >> 8) + 1.0) / 16777217.0;
      u = sqrt(-2.0 * log(u)) * cos(6.2831853071795864769 * v);
    }
    x[i] = u;
  }
  for (int i = 3; i >= 0; i--) {
    iseed[i] = (akv_blas_int)(s & 4095u);
    s >>= 3;
  }
}

static void akv_jacobi_eigen(int n, double *a, double *w, double *v) {
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < n; i++) {
      v[IDX(i, j, n)] = i == j ? 1.0 : 0.0;
    }
  }
  const int max_iter = n > 0 ? 80 * n * n : 1;
  for (int iter = 0; iter < max_iter; iter++) {
    int p = 0;
    int q = n > 1 ? 1 : 0;
    double max_off = 0.0;
    for (int j = 1; j < n; j++) {
      for (int i = 0; i < j; i++) {
        const double off = fabs(a[IDX(i, j, n)]);
        if (off > max_off) {
          max_off = off;
          p = i;
          q = j;
        }
      }
    }
    if (max_off < 1.0e-13) {
      break;
    }
    const double app = a[IDX(p, p, n)];
    const double aqq = a[IDX(q, q, n)];
    const double apq = a[IDX(p, q, n)];
    const double tau = (aqq - app) / (2.0 * apq);
    const double t =
        (tau >= 0.0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau * tau));
    const double c = 1.0 / sqrt(1.0 + t * t);
    const double s = t * c;
    for (int k = 0; k < n; k++) {
      if (k != p && k != q) {
        const double akp = a[IDX(k, p, n)];
        const double akq = a[IDX(k, q, n)];
        a[IDX(k, p, n)] = c * akp - s * akq;
        a[IDX(p, k, n)] = a[IDX(k, p, n)];
        a[IDX(k, q, n)] = s * akp + c * akq;
        a[IDX(q, k, n)] = a[IDX(k, q, n)];
      }
    }
    a[IDX(p, p, n)] = app - t * apq;
    a[IDX(q, q, n)] = aqq + t * apq;
    a[IDX(p, q, n)] = 0.0;
    a[IDX(q, p, n)] = 0.0;
    for (int k = 0; k < n; k++) {
      const double vkp = v[IDX(k, p, n)];
      const double vkq = v[IDX(k, q, n)];
      v[IDX(k, p, n)] = c * vkp - s * vkq;
      v[IDX(k, q, n)] = s * vkp + c * vkq;
    }
  }
  for (int i = 0; i < n; i++) {
    w[i] = a[IDX(i, i, n)];
  }
  for (int i = 0; i < n - 1; i++) {
    int best = i;
    for (int j = i + 1; j < n; j++) {
      if (w[j] < w[best]) {
        best = j;
      }
    }
    if (best != i) {
      akv_swap_double(&w[i], &w[best]);
      for (int r = 0; r < n; r++) {
        akv_swap_double(&v[IDX(r, i, n)], &v[IDX(r, best, n)]);
      }
    }
  }
}

static int akv_selected(const char *range, int idx0, double eig, double vl,
                        double vu, int il, int iu) {
  if (akv_lsame(range, 'A')) {
    return 1;
  }
  if (akv_lsame(range, 'I')) {
    return idx0 + 1 >= il && idx0 + 1 <= iu;
  }
  return eig >= vl && eig <= vu;
}

void bah_akv_internal_dsyevx(const char *jobz, const char *range,
                             const char *uplo, akv_blas_int *n, double *a,
                             akv_blas_int *lda, double *vl, double *vu,
                             akv_blas_int *il, akv_blas_int *iu,
                             double *abstol, akv_blas_int *m, double *w,
                             double *z, akv_blas_int *ldz, double *work,
                             akv_blas_int *lwork, akv_blas_int *iwork,
                             akv_blas_int *ifail, akv_blas_int *info) {
  (void)uplo;
  (void)abstol;
  (void)iwork;
  if (*lwork == -1) {
    work[0] = (double)(8 * (*n > 1 ? *n : 1));
    *info = 0;
    return;
  }
  double *full = akv_alloc((size_t)*n * (size_t)*n);
  double *vec = akv_alloc((size_t)*n * (size_t)*n);
  double *eval = akv_alloc((size_t)*n);
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      full[IDX(i, j, *n)] = akv_sym_get(a, *lda, i, j);
    }
  }
  akv_jacobi_eigen(*n, full, eval, vec);
  *m = 0;
  for (int j = 0; j < *n; j++) {
    if (akv_selected(range, j, eval[j], *vl, *vu, *il, *iu)) {
      w[*m] = eval[j];
      if (akv_lsame(jobz, 'V')) {
        for (int i = 0; i < *n; i++) {
          z[IDX(i, *m, *ldz)] = vec[IDX(i, j, *n)];
        }
        ifail[*m] = 0;
      }
      (*m)++;
    }
  }
  *info = 0;
  free(eval);
  free(vec);
  free(full);
}

static int akv_cholesky_lower(int n, double *b) {
  for (int j = 0; j < n; j++) {
    double diag = b[IDX(j, j, n)];
    for (int k = 0; k < j; k++) {
      diag -= b[IDX(j, k, n)] * b[IDX(j, k, n)];
    }
    if (diag <= 0.0) {
      return j + 1;
    }
    b[IDX(j, j, n)] = sqrt(diag);
    for (int i = j + 1; i < n; i++) {
      double sum = b[IDX(i, j, n)];
      for (int k = 0; k < j; k++) {
        sum -= b[IDX(i, k, n)] * b[IDX(j, k, n)];
      }
      b[IDX(i, j, n)] = sum / b[IDX(j, j, n)];
    }
  }
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < j; i++) {
      b[IDX(i, j, n)] = 0.0;
    }
  }
  return 0;
}

void bah_akv_internal_dsygvx(akv_blas_int *itype, const char *jobz,
                             const char *range, const char *uplo,
                             akv_blas_int *n, double *a, akv_blas_int *lda,
                             double *b, akv_blas_int *ldb, double *vl,
                             double *vu, akv_blas_int *il, akv_blas_int *iu,
                             double *abstol, akv_blas_int *m, double *w,
                             double *z, akv_blas_int *ldz, double *work,
                             akv_blas_int *lwork, akv_blas_int *iwork,
                             akv_blas_int *ifail, akv_blas_int *info) {
  (void)itype;
  (void)uplo;
  (void)abstol;
  (void)iwork;
  if (*lwork == -1) {
    work[0] = (double)(8 * (*n > 1 ? *n : 1));
    *info = 0;
    return;
  }
  double *amat = akv_alloc((size_t)*n * (size_t)*n);
  double *bmat = akv_alloc((size_t)*n * (size_t)*n);
  double *linv = akv_alloc((size_t)*n * (size_t)*n);
  double *tmp = akv_alloc((size_t)*n * (size_t)*n);
  double *c = akv_alloc((size_t)*n * (size_t)*n);
  double *vec = akv_alloc((size_t)*n * (size_t)*n);
  double *eval = akv_alloc((size_t)*n);
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      amat[IDX(i, j, *n)] = akv_sym_get(a, *lda, i, j);
      bmat[IDX(i, j, *n)] = akv_sym_get(b, *ldb, i, j);
    }
  }
  const int chol_info = akv_cholesky_lower(*n, bmat);
  if (chol_info != 0) {
    *info = *n + chol_info;
    goto cleanup;
  }
  for (int col = 0; col < *n; col++) {
    for (int i = 0; i < *n; i++) {
      double sum = i == col ? 1.0 : 0.0;
      for (int k = 0; k < i; k++) {
        sum -= bmat[IDX(i, k, *n)] * linv[IDX(k, col, *n)];
      }
      linv[IDX(i, col, *n)] = sum / bmat[IDX(i, i, *n)];
    }
  }
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      double sum = 0.0;
      for (int k = 0; k < *n; k++) {
        sum += linv[IDX(i, k, *n)] * amat[IDX(k, j, *n)];
      }
      tmp[IDX(i, j, *n)] = sum;
    }
  }
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      double sum = 0.0;
      for (int k = 0; k < *n; k++) {
        sum += tmp[IDX(i, k, *n)] * linv[IDX(j, k, *n)];
      }
      c[IDX(i, j, *n)] = 0.5 * (sum + sum);
    }
  }
  akv_jacobi_eigen(*n, c, eval, vec);
  *m = 0;
  for (int j = 0; j < *n; j++) {
    if (akv_selected(range, j, eval[j], *vl, *vu, *il, *iu)) {
      w[*m] = eval[j];
      if (akv_lsame(jobz, 'V')) {
        for (int i = 0; i < *n; i++) {
          double sum = 0.0;
          for (int k = 0; k < *n; k++) {
            sum += linv[IDX(k, i, *n)] * vec[IDX(k, j, *n)];
          }
          z[IDX(i, *m, *ldz)] = sum;
        }
        double norm = 0.0;
        for (int p = 0; p < *n; p++) {
          for (int q = 0; q < *n; q++) {
            norm += z[IDX(p, *m, *ldz)] * akv_sym_get(b, *ldb, p, q) *
                    z[IDX(q, *m, *ldz)];
          }
        }
        norm = sqrt(fabs(norm));
        if (norm > 0.0) {
          for (int i = 0; i < *n; i++) {
            z[IDX(i, *m, *ldz)] /= norm;
          }
        }
        ifail[*m] = 0;
      }
      (*m)++;
    }
  }
  *info = 0;
cleanup:
  free(eval);
  free(vec);
  free(c);
  free(tmp);
  free(linv);
  free(bmat);
  free(amat);
}

void bah_akv_internal_dgesvd(const char *jobu, const char *jobvt,
                             akv_blas_int *m, akv_blas_int *n, double *a,
                             akv_blas_int *lda, double *s, double *u,
                             akv_blas_int *ldu, double *vt,
                             akv_blas_int *ldvt, double *work,
                             akv_blas_int *lwork, akv_blas_int *info) {
  if (*lwork == -1) {
    work[0] = (double)(8 * ((*m > *n ? *m : *n) + 1));
    *info = 0;
    return;
  }
  const int minmn = *m < *n ? *m : *n;
  double *ata = akv_alloc((size_t)*n * (size_t)*n);
  double *vec = akv_alloc((size_t)*n * (size_t)*n);
  double *eval = akv_alloc((size_t)*n);
  for (int j = 0; j < *n; j++) {
    for (int i = 0; i < *n; i++) {
      double sum = 0.0;
      for (int k = 0; k < *m; k++) {
        sum += a[IDX(k, i, *lda)] * a[IDX(k, j, *lda)];
      }
      ata[IDX(i, j, *n)] = sum;
    }
  }
  akv_jacobi_eigen(*n, ata, eval, vec);
  for (int i = 0; i < *n / 2; i++) {
    akv_swap_double(&eval[i], &eval[*n - 1 - i]);
    for (int r = 0; r < *n; r++) {
      akv_swap_double(&vec[IDX(r, i, *n)], &vec[IDX(r, *n - 1 - i, *n)]);
    }
  }
  for (int j = 0; j < minmn; j++) {
    s[j] = sqrt(eval[j] > 0.0 ? eval[j] : 0.0);
  }
  if (!akv_lsame(jobvt, 'N')) {
    for (int j = 0; j < *n; j++) {
      for (int i = 0; i < *n; i++) {
        vt[IDX(i, j, *ldvt)] = vec[IDX(j, i, *n)];
      }
    }
  }
  if (!akv_lsame(jobu, 'N')) {
    for (int j = 0; j < minmn; j++) {
      for (int i = 0; i < *m; i++) {
        double sum = 0.0;
        for (int k = 0; k < *n; k++) {
          sum += a[IDX(i, k, *lda)] * vec[IDX(k, j, *n)];
        }
        u[IDX(i, j, *ldu)] = s[j] > 0.0 ? sum / s[j] : 0.0;
      }
    }
  }
  *info = 0;
  free(eval);
  free(vec);
  free(ata);
}

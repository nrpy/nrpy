"""
Scaffold registration and implementation module.

This code registers the AKV C scaffold implements the equations
through SymPy to generate efficient C code,
referencing the callbacks defined in the scaffold.

Author: Wesley Inselman
"""

from __future__ import annotations

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity.bhahaha.approx_killing_vector_spin import (
    ApproxKillingSpinClass,
)
from nrpy.infrastructures import BHaH


def register_CFunction_diagnostics_approx_killing_vector_spin(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BHaHAHA approximate Killing vector (AKV) spin diagnostic C function.

    This generates C code for the AKV spin diagnostic, including code for the
    horizon surface element and the L1-reduced AKV integrands, then registers
    the resulting 'diagnostics_approx_killing_vector_spin' C function with NRPy. During the
    parallel-codegen registration phase, the function records the call and
    returns 'None' without generating code.

    :param enable_fd_functions: Whether to have NRPy emit finite-difference helper
        functions when generating C code for expressions that require finite
        differencing.

    :return: 'None' during the parallel-codegen registration phase. Otherwise,
        returns the populated NRPy environment after registering the generated
        C function.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Robustly generate C for sqrt(det q) from area3() across NRPy versions.
    # Some versions accept a single outvar string; others require a list-of-lvalues.
    area_codegen = ccg.c_codegen(
        BHaH.BHaHAHA.area.area3(),
        "const REAL sqrt_det_q",
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
        enable_cse=True,
        include_braces=True,
        verbose=False,
    )

    # Build AKV L1-reduced per-point integrands.
    # NOTE: ApproxKillingSpinClass defines Hmn/Nmn/Jm as *densities* (they include sqrt(q)).
    # The AKV scaffold contract expects the callback to return *integrands with no quadrature weight*.
    # In this BHaHAHA pipeline, w2d is built including sqrt_det_q, so by default we return
    # density-free integrands by dividing by sqrt(q) at the callback boundary (C-side macro control).
    akv = ApproxKillingSpinClass("Spherical_rfm_precompute")
    Hraw = akv.Hmn_integrand
    Nraw = akv.Nmn_integrand
    Jraw = akv.Jm_integrand

    akv_exprs = (
        [Hraw[m][n] for m in range(3) for n in range(3)]
        + [Nraw[m][n] for m in range(3) for n in range(3)]
        + [Jraw[m] for m in range(3)]
    )

    akv_outvars = (
        [f"const REAL akv_Hraw_{m}_{n}" for m in range(3) for n in range(3)]
        + [f"const REAL akv_Nraw_{m}_{n}" for m in range(3) for n in range(3)]
        + [f"const REAL akv_Jraw_{m}" for m in range(3)]
    )

    akv_codegen = ccg.c_codegen(
        akv_exprs,
        akv_outvars,
        enable_cse=True,
        include_braces=True,
        verbose=False,
        enable_fd_codegen=True,
        enable_fd_functions=enable_fd_functions,
    )

    prefunc = (
        r"""
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdbool.h>

#ifndef AKV_REAL
#define AKV_REAL REAL
#endif

#ifndef BHAHAHA_AKV_WEIGHTS_INCLUDE_DX
/* In this codebase, bah_diagnostics_integration_weights returns *dimensionless* repeating
 * stencil coefficients (see bah_diagnostics_integration_weights.c). Therefore the correct
 * default is to include physical spacing factors (dx1*dx2) in w2d.
 * If a downstream fork changes the weights routine to include spacing, set to 0 to avoid
 * double-counting.
 */
#define BHAHAHA_AKV_WEIGHTS_INCLUDE_DX 1
#endif

#ifndef BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION
#define BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION 1
#endif
// Quality-flag extension bits used by this wrapper (must not conflict with scaffold bits).
// If conflicts arise in your build, set BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION=0.
#ifndef BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT
#define BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT 4
#endif
#ifndef BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT
#define BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT 5
#endif

#if BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION
// Sanity checks for extension-bit configuration (avoid undefined shifts).
#if (BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT < 0) || (BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT > 30)
#error "BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT must be in [0,30]"
#endif
#if (BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT < 0) || (BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT > 24)
#error "BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT must be in [0,24]"
#endif
#endif  // BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION
#ifndef BHAHAHA_AKV_CALLBACK_DENSITY_FREE
#define BHAHAHA_AKV_CALLBACK_DENSITY_FREE 1
#endif
#ifndef BHAHAHA_AKV_EXPR_ARE_DENSITIES
#define BHAHAHA_AKV_EXPR_ARE_DENSITIES 1
#endif
#ifndef BHAHAHA_AKV_AH_SLAB_I0
#define BHAHAHA_AKV_AH_SLAB_I0 NGHOSTS
#endif

#ifndef BHAHAHA_AKV_ENABLE_SPIN_VECTOR
#define BHAHAHA_AKV_ENABLE_SPIN_VECTOR 0
#endif




/* w2d is constructed as:
 *   w2d = sqrt_det_q * wtheta * wphi * (optional dx1*dx2).
 * If w2d includes sqrt_det_q (as below), then density-form expressions must be converted
 * to density-free integrands in the callback by dividing by sqrt_det_q_cb.
 * This pairing is controlled by:
 *   - BHAHAHA_AKV_CALLBACK_DENSITY_FREE (default 1)
 *   - BHAHAHA_AKV_EXPR_ARE_DENSITIES   (default 1)
 */

#ifndef BHAHAHA_AKV_ENABLE_GRIDPOINT
#define BHAHAHA_AKV_ENABLE_GRIDPOINT 0
#endif


/* TLS portability.
 * We use TLS only for the per-thread w2d cache. There is deliberately no TLS/global fallback
 * for stub-usage flags (those are accumulated per-call via ctx->stub_flags_ptr).
 */
#if defined(__STDC_VERSION__) && (__STDC_VERSION__ >= 201112L)
  #define BHAHAHA_TLS _Thread_local
  #define BHAHAHA_HAVE_TLS 1
#elif defined(__GNUC__) || defined(__clang__)
  #define BHAHAHA_TLS __thread
  #define BHAHAHA_HAVE_TLS 1
#else
  #define BHAHAHA_TLS /* none */
  #define BHAHAHA_HAVE_TLS 0
#endif

#ifndef BHAHAHA_BHAHAHA_NRPY_AKV_SCAFFOLD_INCLUDED
#define BHAHAHA_BHAHAHA_NRPY_AKV_SCAFFOLD_INCLUDED 1

/*
 * Start of scaffold
 *
 * Compilable AKV diagnostic scaffold with two workflows:
 *   (A) L1-reduced AKV: assemble 3x3 generalized eigenproblem (H x = λ N x) and solve.
 *   (B) Gridpoint-basis AKV: assemble dense K,B from bilinear-form/stencil callbacks, enforce mean-zero
 *       via an explicit (N-1) basis, solve generalized eigenproblem in that subspace, and postprocess.
 *
 * Geometric/equation content is intentionally excluded and must be supplied via callbacks.
 *
 * Contracts:
 *
 * 1) L1-reduced integrands:
 *    - The per-point Jm integrand MUST NOT include the global 1/(8π) factor.
 *    - The callback returns integrands (no quadrature weight applied inside callback).
 *
 * 2) Gridpoint stencil/bilinear assembly callback:
 *    - i_row is an interior DOF index in [0, N_theta*N_phi).
 *    - DOF ordering is row-major over interior points:
 *        it = i_row / N_phi, ip = i_row % N_phi,
 *      where it=0..N_theta-1 and ip=0..N_phi-1 are interior indices.
 *    - A canonical mapping to full ghosted-grid indices is:
 *        it_full = it + NG_theta, ip_full = ip + NG_phi, p = it_full*N_phi_tot + ip_full.
 *    - Neighbor indices in j_cols must use the same interior DOF indexing as i_row.
 *    - For a given i_row, (i_row, j) pairs should be unique. If duplicates appear, they are summed.
 *
 * 3) Dense storage:
 *    - Dense matrices are stored column-major: M[i + N*j].
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef USE_LAPACKE
#include <lapacke.h>
#endif

// -----------------------------
// Error codes
// -----------------------------
typedef enum {
  AKV_SUCCESS = 0,
  AKV_ERR_NULLPTR = 1,
  AKV_ERR_BADPARAM = 2,
  AKV_ERR_ALLOC = 3,
  AKV_ERR_EIGEN_FAIL = 4,
  AKV_ERR_NOT_IMPLEMENTED = 5
} akv_error_t;

// -----------------------------
// Method selection
// -----------------------------
typedef enum {
  AKV_METHOD_L1_REDUCED = 0,
  AKV_METHOD_GRIDPOINT  = 1
} akv_method_t;

// -----------------------------
// Horizon grid description
// -----------------------------
typedef struct {
  int N_theta, N_phi;     // interior
  int NG_theta, NG_phi;   // ghost zones
  AKV_REAL dtheta, dphi;

  int N_theta_tot;
  int N_phi_tot;
  int N_tot;

  // Quadrature weights on the full (ghosted) grid (ghosts unused in integration).
  const AKV_REAL *w; // length N_tot

  // Optional per-point geometry pointer required by user callbacks.
  const void *geom;
} akv_horizon_grid_t;

// -----------------------------
// Diagnostics output
// -----------------------------
typedef struct {
  AKV_REAL akv_lambda[3];
  AKV_REAL akv_J[3];
  AKV_REAL akv_a[3];
  AKV_REAL akv_spin_vec[3];
  AKV_REAL akv_eig_gap_43;
  AKV_REAL akv_eig_resid[3];

  // Quality flag bitfield (0 = OK):
  //   bit0: residual tolerance exceeded for any reported mode
  //   bit1: small eigenvalue gap (gridpoint method)
  //   bit2: solver/library unavailable
  //   bit3: B not SPD / regularization retry used / near-singular handling invoked
  int akv_quality_flag;

  akv_method_t method_used;

  // L1-reduced: generalized eigenvectors (columns) in rigid-rotation coefficient basis.
  AKV_REAL akv_l1_evecs[3][3];

  // Gridpoint: store the lowest three full-space mean-zero eigenvectors z(p) (interior DOF space).
  // These are allocated by the solver; call akv_diagnostics_free() to release.
  int gp_N;                 // = N_theta*N_phi if available
  AKV_REAL *gp_z[3];        // each length gp_N, or NULL if not computed/available
} akv_diagnostics_t;

void akv_diagnostics_free(akv_diagnostics_t *d) {
  if (!d) return;
  for (int a = 0; a < 3; a++) {
    free(d->gp_z[a]);
    d->gp_z[a] = NULL;
  }
  d->gp_N = 0;
}

// -----------------------------
// Runtime parameters
// -----------------------------
typedef struct {
  akv_method_t method;

  int l1_choose_index;

  int full_num_eigs;     // >=4 recommended for gap diagnostic
  AKV_REAL eig_tol;      // e.g., 1e-10 to 1e-8
  AKV_REAL reg_eps;      // initial diagonal regularization for B_red
  AKV_REAL reg_eps_max;  // maximum regularization allowed
  int reg_max_tries;     // retry count for SPD failures
  AKV_REAL gap_ratio_thresh; // e.g., 1.5
  AKV_REAL horizon_area; // the area of the black hole horizon

  bool build_spin_vector;

  // If true and stencil callback is absent, allow debug-only O(N^3) assembly via applyK/applyB.
  bool allow_debug_assembly;
} akv_params_t;

// -----------------------------
// Callbacks (equations intentionally omitted)
// -----------------------------

// L1-reduced per-point integrands (no quadrature weight applied inside callback).
// Jm integrand must NOT include the global 1/(8π) factor.
typedef void (*akv_eval_l1_integrands_f)(
    int p,
    const akv_horizon_grid_t *grid,
    AKV_REAL Hmn[3][3],
    AKV_REAL Nmn[3][3],
    AKV_REAL Jm[3]);

// Optional mapping: given 3 coefficients, compute background-frame spin vector S_out[3].
typedef void (*akv_map_to_spinvec_f)(
    const AKV_REAL coeffs[3],
    const akv_horizon_grid_t *grid,
    AKV_REAL S_out[3]);

// Optional sign-fix evaluator for L1-reduced modes.
// Return a reference scalar that should be positive after sign convention is applied.
typedef AKV_REAL (*akv_eval_sign_l1_f)(
    const AKV_REAL coeffs[3],
    const AKV_REAL Jm_assembled[3],
    const akv_horizon_grid_t *grid);

// Gridpoint operator hooks (optional debug-only assembly and/or residual checks without dense matrices).
typedef void (*akv_apply_K_f)(const akv_horizon_grid_t *grid, const AKV_REAL *z_in, AKV_REAL *Kz_out);
typedef void (*akv_apply_B_f)(const akv_horizon_grid_t *grid, const AKV_REAL *z_in, AKV_REAL *Bz_out);

// Gridpoint-basis stencil/bilinear-form row evaluation for dense matrix assembly.
typedef int (*akv_eval_row_stencil_f)(
    int i_row,
    const akv_horizon_grid_t *grid,
    int max_entries,
    int *j_cols,
    AKV_REAL *K_vals,
    AKV_REAL *B_vals);

// Optional sign-fix evaluator for gridpoint modes, given the full-space mean-zero eigenvector z (interior DOF space).
typedef AKV_REAL (*akv_eval_sign_full_f)(
    const AKV_REAL *z_full,
    int N_full,
    const akv_horizon_grid_t *grid);

// Optional angular momentum integrand callback for gridpoint eigenmodes (loop scaffold).
// Return the per-point integrand (no weight, no 1/(8π)); caller multiplies by w[p].
typedef AKV_REAL (*akv_eval_J_integrand_full_f)(
    int p,
    const akv_horizon_grid_t *grid,
    const AKV_REAL *z_full,
    int N_full);

// -----------------------------
// Indexing helpers
// -----------------------------
static inline int akv_idx2_full(const akv_horizon_grid_t *g, int it_full, int ip_full) {
  return it_full * g->N_phi_tot + ip_full;
}

// Interior DOF index -> interior (it,ip)
static inline void akv_dof_to_it_ip(const akv_horizon_grid_t *g, int dof, int *it, int *ip) {
  (void)g;
  *it = dof / g->N_phi;
  *ip = dof - (*it) * g->N_phi;
}

// Interior (it,ip) -> interior DOF index
static inline int akv_it_ip_to_dof(const akv_horizon_grid_t *g, int it, int ip) {
  return it * g->N_phi + ip;
}

// Interior DOF -> full ghosted p
static inline int akv_dof_to_full_p(const akv_horizon_grid_t *g, int dof) {
  int it, ip;
  akv_dof_to_it_ip(g, dof, &it, &ip);
  int it_full = it + g->NG_theta;
  int ip_full = ip + g->NG_phi;
  return akv_idx2_full(g, it_full, ip_full);
}

// -----------------------------
// Small linear algebra helpers
// -----------------------------
static inline void mat3_zero(AKV_REAL A[3][3]) { memset(A, 0, 9 * sizeof(AKV_REAL)); }
static inline void vec3_zero(AKV_REAL v[3]) { memset(v, 0, 3 * sizeof(AKV_REAL)); }

static inline void mat3_add_inplace(AKV_REAL A[3][3], const AKV_REAL B[3][3], AKV_REAL alpha) {
  for (int i = 0; i < 3; i++) for (int j = 0; j < 3; j++) A[i][j] += alpha * B[i][j];
}
static inline void vec3_add_inplace(AKV_REAL a[3], const AKV_REAL b[3], AKV_REAL alpha) {
  for (int i = 0; i < 3; i++) a[i] += alpha * b[i];
}
static inline AKV_REAL vec3_dot(const AKV_REAL a[3], const AKV_REAL b[3]) {
  return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
static inline AKV_REAL vecn_dot(const AKV_REAL *a, const AKV_REAL *b, int n) {
  AKV_REAL s = 0;
  for (int i = 0; i < n; i++) s += a[i]*b[i];
  return s;
}
static inline AKV_REAL vecn_norm2(const AKV_REAL *a, int n) {
  return sqrt(vecn_dot(a, a, n));
}
static inline void vecn_scale(AKV_REAL *a, int n, AKV_REAL alpha) {
  for (int i = 0; i < n; i++) a[i] *= alpha;
}
static inline void vecn_axpy(AKV_REAL *y, const AKV_REAL *x, int n, AKV_REAL alpha) {
  for (int i = 0; i < n; i++) y[i] += alpha * x[i];
}

static void akv_compute_dimensionless_spins(const AKV_REAL J[3], AKV_REAL horizon_area, AKV_REAL a_out[3]) {
  a_out[0] = 0.0;
  a_out[1] = 0.0;
  a_out[2] = 0.0;

  // Whitepaper normalization:
  // M_irr^2 = A/(16*pi), M_H^2 = M_irr^2 + J_ref^2/(4*M_irr^2), a_i = J_i/M_H^2.
  if (!(horizon_area > 0.0)) return;
  const AKV_REAL M_irr2 = horizon_area / (16.0 * (AKV_REAL)M_PI);
  if (!(M_irr2 > 0.0)) return;

  const AKV_REAL Jref2 = J[0]*J[0] + J[1]*J[1] + J[2]*J[2];
  const AKV_REAL M_H2 = M_irr2 + Jref2 / (4.0 * M_irr2);
  if (!(M_H2 > 0.0)) return;

  const AKV_REAL inv_M_H2 = 1.0 / M_H2;
  a_out[0] = J[0] * inv_M_H2;
  a_out[1] = J[1] * inv_M_H2;
  a_out[2] = J[2] * inv_M_H2;
}

static void mat3_mul_vec(const AKV_REAL A[3][3], const AKV_REAL x[3], AKV_REAL y[3]) {
  for (int i = 0; i < 3; i++) {
    y[i] = 0;
    for (int j = 0; j < 3; j++) y[i] += A[i][j] * x[j];
  }
}

static bool chol3_lower(const AKV_REAL N[3][3], AKV_REAL L[3][3]) {
  mat3_zero(L);
  AKV_REAL a00 = N[0][0];
  if (a00 <= 0) return false;
  L[0][0] = sqrt(a00);

  L[1][0] = N[1][0] / L[0][0];
  L[2][0] = N[2][0] / L[0][0];

  AKV_REAL a11 = N[1][1] - L[1][0]*L[1][0];
  if (a11 <= 0) return false;
  L[1][1] = sqrt(a11);

  L[2][1] = (N[2][1] - L[2][0]*L[1][0]) / L[1][1];

  AKV_REAL a22 = N[2][2] - L[2][0]*L[2][0] - L[2][1]*L[2][1];
  if (a22 <= 0) return false;
  L[2][2] = sqrt(a22);

  return true;
}

static void tri3_solve_lower(const AKV_REAL L[3][3], const AKV_REAL b[3], AKV_REAL y[3]) {
  y[0] = b[0] / L[0][0];
  y[1] = (b[1] - L[1][0]*y[0]) / L[1][1];
  y[2] = (b[2] - L[2][0]*y[0] - L[2][1]*y[1]) / L[2][2];
}
static void tri3_solve_upper_from_lowerT(const AKV_REAL L[3][3], const AKV_REAL y[3], AKV_REAL x[3]) {
  x[2] = y[2] / L[2][2];
  x[1] = (y[1] - L[2][1]*x[2]) / L[1][1];
  x[0] = (y[0] - L[1][0]*x[1] - L[2][0]*x[2]) / L[0][0];
}

static void generalized_to_standard3(const AKV_REAL H[3][3], const AKV_REAL L[3][3], AKV_REAL C[3][3]) {
  AKV_REAL A[3][3];

  for (int i = 0; i < 3; i++) {
    AKV_REAL h[3] = { H[0][i], H[1][i], H[2][i] };
    AKV_REAL a[3];
    tri3_solve_lower(L, h, a); // a = L^-1 h
    A[0][i] = a[0]; A[1][i] = a[1]; A[2][i] = a[2]; // A = L^-1 H
  }

  for (int i = 0; i < 3; i++) {
    AKV_REAL c[3];
    tri3_solve_lower(L, A[i], c); // c = L^-1 a^T
                                  //   = L^-1 h L^-T (H is symmetric)
    C[0][i] = c[0]; C[1][i] = c[1]; C[2][i] = c[2];
  }

  // Fixes asymmetry that may be caused by floating point precision.
  for (int i = 0; i < 3; i++) for (int j = i+1; j < 3; j++) {
    AKV_REAL s = (C[i][j] + C[j][i])/2;
    C[i][j] = s; C[j][i] = s;
  }
}

static void jacobi_eig_sym3(const AKV_REAL A_in[3][3], AKV_REAL eval[3], AKV_REAL V[3][3]) {
  AKV_REAL A[3][3];
  memcpy(A, A_in, 9*sizeof(AKV_REAL));

  mat3_zero(V);
  V[0][0]=1; V[1][1]=1; V[2][2]=1;

  const int max_sweeps = 50;
  const AKV_REAL tol = 1e-14;

  for (int sweep = 0; sweep < max_sweeps; sweep++) {
    int p=0,q=1;
    AKV_REAL maxv = fabs(A[0][1]);
    if (fabs(A[0][2]) > maxv) { maxv = fabs(A[0][2]); p=0; q=2; }
    if (fabs(A[1][2]) > maxv) { maxv = fabs(A[1][2]); p=1; q=2; }
    if (maxv < tol) break;

    AKV_REAL app = A[p][p];
    AKV_REAL aqq = A[q][q];
    AKV_REAL apq = A[p][q];

    AKV_REAL tau = (aqq - app) / (2.0 * apq);
    AKV_REAL t = (tau >= 0 ? 1.0 : -1.0) / (fabs(tau) + sqrt(1.0 + tau*tau));
    AKV_REAL c = 1.0 / sqrt(1.0 + t*t);
    AKV_REAL s = t * c;

    A[p][p] = app - t*apq;
    A[q][q] = aqq + t*apq;
    A[p][q] = 0.0; A[q][p] = 0.0;

    for (int r = 0; r < 3; r++) {
      if (r == p || r == q) continue;
      AKV_REAL arp = A[r][p];
      AKV_REAL arq = A[r][q];
      A[r][p] = c*arp - s*arq; A[p][r] = A[r][p];
      A[r][q] = c*arq + s*arp; A[q][r] = A[r][q];
    }

    for (int r = 0; r < 3; r++) {
      AKV_REAL vrp = V[r][p];
      AKV_REAL vrq = V[r][q];
      V[r][p] = c*vrp - s*vrq;
      V[r][q] = c*vrq + s*vrp;
    }
  }

  eval[0] = A[0][0];
  eval[1] = A[1][1];
  eval[2] = A[2][2];
}

static void sort_eigs3(AKV_REAL eval[3], AKV_REAL V[3][3]) {
  for (int i = 0; i < 3; i++) {
    int k = i;
    for (int j = i+1; j < 3; j++) if (eval[j] < eval[k]) k = j;
    if (k != i) {
      AKV_REAL te = eval[i]; eval[i] = eval[k]; eval[k] = te;
      for (int r = 0; r < 3; r++) {
        AKV_REAL tv = V[r][i]; V[r][i] = V[r][k]; V[r][k] = tv;
      }
    }
  }
}

static void normalize_gen_evec3(const AKV_REAL N[3][3], AKV_REAL x[3]) {
  AKV_REAL Nx[3];
  mat3_mul_vec(N, x, Nx);
  AKV_REAL nrm2 = vec3_dot(x, Nx);
  if (nrm2 > 0) {
    AKV_REAL inv = 1.0 / sqrt(nrm2);
    x[0] *= inv; x[1] *= inv; x[2] *= inv;
  }
}

// -----------------------------
// Gridpoint mean-zero basis utilities with cached map and buffers
// -----------------------------
typedef struct {
  int N_full;      // N_theta*N_phi
  int N_red;       // N_full - 1
  int ref;         // reference DOF index removed
  AKV_REAL *w_dof; // length N_full, per-DOF weights
  AKV_REAL wref;   // w_dof[ref]
  int *map;        // length N_red, reduced index -> full index (skip ref)

  // Scratch buffers (reused)
  AKV_REAL *tmpN;  // length max(N_full,N_red)
  AKV_REAL *tmpK;  // length max(N_full,N_red)
  AKV_REAL *tmpR;  // length max(N_full,N_red)
} akv_mean0_basis_t;

static void akv_mean0_basis_free(akv_mean0_basis_t *b) {
  if (!b) return;
  free(b->w_dof);
  free(b->map);
  free(b->tmpN);
  free(b->tmpK);
  free(b->tmpR);
  memset(b, 0, sizeof(*b));
}

static akv_error_t akv_mean0_basis_build(const akv_horizon_grid_t *grid, akv_mean0_basis_t *b) {
  if (!grid || !b || !grid->w) return AKV_ERR_NULLPTR;
  const int N_full = grid->N_theta * grid->N_phi;
  if (N_full <= 1) return AKV_ERR_BADPARAM;

  memset(b, 0, sizeof(*b));
  b->N_full = N_full;
  b->N_red = N_full - 1;
  b->ref = -1;

  b->w_dof = (AKV_REAL*)calloc((size_t)N_full, sizeof(AKV_REAL));
  b->map   = (int*)malloc((size_t)b->N_red * sizeof(int));
  if (!b->w_dof || !b->map) { akv_mean0_basis_free(b); return AKV_ERR_ALLOC; }

  for (int i = 0; i < N_full; i++) {
    int p = akv_dof_to_full_p(grid, i);
    b->w_dof[i] = grid->w[p];
  }

  for (int i = 0; i < N_full; i++) {
    if (fabs(b->w_dof[i]) > 0) { b->ref = i; b->wref = b->w_dof[i]; break; }
  }
  if (b->ref < 0 || b->wref == 0.0) { akv_mean0_basis_free(b); return AKV_ERR_BADPARAM; }

  for (int i = 0, k = 0; i < N_full; i++) if (i != b->ref) b->map[k++] = i;

  int scratch_len = (b->N_full > b->N_red) ? b->N_full : b->N_red;
  b->tmpN = (AKV_REAL*)malloc((size_t)scratch_len * sizeof(AKV_REAL));
  b->tmpK = (AKV_REAL*)malloc((size_t)scratch_len * sizeof(AKV_REAL));
  b->tmpR = (AKV_REAL*)malloc((size_t)scratch_len * sizeof(AKV_REAL));
  if (!b->tmpN || !b->tmpK || !b->tmpR) { akv_mean0_basis_free(b); return AKV_ERR_ALLOC; }

  return AKV_SUCCESS;
}

// Transform full-space dense matrix M_full (N_full x N_full) into reduced mean-zero matrix (N_red x N_red):
// M_red = B^T M_full B, where basis columns are b_i = e_i - (w_i/w_ref) e_ref for i != ref.
static void akv_mean0_transform_dense(
    const akv_mean0_basis_t *b,
    const AKV_REAL *M_full, int ld_full,
    AKV_REAL *M_red, int ld_red) {

  const int N = b->N_full;
  const int Nr = b->N_red;
  const int r = b->ref;

  for (int a = 0; a < Nr; a++) {
    const int ia = b->map[a];
    const AKV_REAL alpha_a = -b->w_dof[ia] / b->wref;

    for (int c = 0; c < Nr; c++) {
      const int ic = b->map[c];
      const AKV_REAL alpha_c = -b->w_dof[ic] / b->wref;

      const AKV_REAL M_ia_ic = M_full[ia + (size_t)ld_full * ic];
      const AKV_REAL M_ia_r  = M_full[ia + (size_t)ld_full * r];
      const AKV_REAL M_r_ic  = M_full[r  + (size_t)ld_full * ic];
      const AKV_REAL M_r_r   = M_full[r  + (size_t)ld_full * r];

      (void)N;
      M_red[a + (size_t)ld_red * c] =
          M_ia_ic + alpha_c * M_ia_r + alpha_a * M_r_ic + (alpha_a * alpha_c) * M_r_r;
    }
  }
}

// Lift reduced vector x_red (Nr) to full mean-zero vector x_full (N_full) using the cached basis.
static void akv_mean0_lift_vector(
    const akv_mean0_basis_t *b,
    const AKV_REAL *x_red,
    AKV_REAL *x_full) {

  const int N = b->N_full;
  const int Nr = b->N_red;
  const int r = b->ref;

  memset(x_full, 0, (size_t)N * sizeof(AKV_REAL));

  AKV_REAL xref = 0.0;
  for (int a = 0; a < Nr; a++) {
    const int ia = b->map[a];
    x_full[ia] = x_red[a];
    xref += (-b->w_dof[ia] / b->wref) * x_red[a];
  }
  x_full[r] = xref;
}

// -----------------------------
// Dense helpers (column-major)
// -----------------------------
static void akv_symmetrize_dense(AKV_REAL *M, int N) {
  for (int i = 0; i < N; i++) {
    for (int j = i+1; j < N; j++) {
      AKV_REAL s = 0.5*(M[i + (size_t)N*j] + M[j + (size_t)N*i]);
      M[i + (size_t)N*j] = s;
      M[j + (size_t)N*i] = s;
    }
  }
}

static void akv_dense_add_diag(AKV_REAL *M, int N, AKV_REAL eps) {
  for (int i = 0; i < N; i++) M[i + (size_t)N*i] += eps;
}

static void akv_dense_matvec(const AKV_REAL *M, int N, const AKV_REAL *x, AKV_REAL *y) {
  for (int i = 0; i < N; i++) {
    AKV_REAL s = 0.0;
    const AKV_REAL *col = &M[i]; // M[i + N*j]
    for (int j = 0; j < N; j++) s += col[(size_t)N*j] * x[j];
    y[i] = s;
  }
}

// Residual ratio for dense matrices using caller-provided scratch buffers of length N.
static AKV_REAL akv_residual_ratio_dense_inplace(
    const AKV_REAL *K, const AKV_REAL *B, int N,
    const AKV_REAL *z, AKV_REAL lambda,
    AKV_REAL *Kz, AKV_REAL *Bz, AKV_REAL *r) {

  akv_dense_matvec(K, N, z, Kz);
  akv_dense_matvec(B, N, z, Bz);
  for (int i = 0; i < N; i++) r[i] = Kz[i] - lambda * Bz[i];

  AKV_REAL nr = vecn_norm2(r, N);
  AKV_REAL nK = vecn_norm2(Kz, N);
  AKV_REAL nB = vecn_norm2(Bz, N);
  AKV_REAL denom = nK + nB;
  return (denom > 0) ? (nr / denom) : 0.0;
}

/*
// Residual ratio using operator application (no dense matrices), with caller-provided scratch buffers of length N.
static AKV_REAL akv_residual_ratio_apply_inplace(
    const akv_horizon_grid_t *grid,
    akv_apply_K_f applyK,
    akv_apply_B_f applyB,
    int N,
    const AKV_REAL *z,
    AKV_REAL lambda,
    AKV_REAL *Kz, AKV_REAL *Bz, AKV_REAL *r) {

  applyK(grid, z, Kz);
  applyB(grid, z, Bz);
  for (int i = 0; i < N; i++) r[i] = Kz[i] - lambda * Bz[i];

  AKV_REAL nr = vecn_norm2(r, N);
  AKV_REAL nK = vecn_norm2(Kz, N);
  AKV_REAL nB = vecn_norm2(Bz, N);
  AKV_REAL denom = nK + nB;
  return (denom > 0) ? (nr / denom) : 0.0;
}
*/

// B-inner-product Gram–Schmidt on k vectors (columns) in V (N x k), using dense B (N x N).
// Reuses tmp buffer of length N.
static void akv_b_gram_schmidt_inplace(AKV_REAL *V, int N, int k, const AKV_REAL *B, AKV_REAL *tmp) {
  for (int a = 0; a < k; a++) {
    AKV_REAL *va = &V[(size_t)N * a];

    for (int b = 0; b < a; b++) {
      AKV_REAL *vb = &V[(size_t)N * b];

      akv_dense_matvec(B, N, vb, tmp);
      AKV_REAL proj = vecn_dot(va, tmp, N); // <va,vb>_B
      vecn_axpy(va, vb, N, -proj);
    }

    akv_dense_matvec(B, N, va, tmp);
    AKV_REAL n2 = vecn_dot(va, tmp, N);
    if (n2 > 0) vecn_scale(va, N, 1.0 / sqrt(n2));
  }
}

// -----------------------------
// L1-reduced assembly + solve + postprocess
// -----------------------------
static AKV_REAL default_sign_ref_l1(const AKV_REAL coeffs[3], const AKV_REAL Jm_assembled[3]) {
  return coeffs[0]*Jm_assembled[0] + coeffs[1]*Jm_assembled[1] + coeffs[2]*Jm_assembled[2];
}

static akv_error_t akv_solve_l1_reduced(
    const akv_horizon_grid_t *grid,
    const akv_params_t *pars,
    akv_eval_l1_integrands_f eval_l1,
    akv_eval_sign_l1_f eval_sign_l1,
    akv_map_to_spinvec_f map_to_spinvec,
    akv_diagnostics_t *out) {

  if (!grid || !pars || !eval_l1 || !out) return AKV_ERR_NULLPTR;
  if (!grid->w) return AKV_ERR_NULLPTR;

  AKV_REAL H[3][3], N[3][3], Jm[3];
  mat3_zero(H); mat3_zero(N); vec3_zero(Jm);

#ifdef _OPENMP
  AKV_REAL H_acc[3][3], N_acc[3][3], J_acc[3];
  mat3_zero(H_acc); mat3_zero(N_acc); vec3_zero(J_acc);
#pragma omp parallel
  {
    AKV_REAL Hloc[3][3], Nloc[3][3], Jloc[3];
    mat3_zero(Hloc); mat3_zero(Nloc); vec3_zero(Jloc);

#pragma omp for collapse(2) nowait
    for (int it = grid->NG_theta; it < grid->NG_theta + grid->N_theta; it++) {
      for (int ip = grid->NG_phi; ip < grid->NG_phi + grid->N_phi; ip++) {
        const int p = akv_idx2_full(grid, it, ip);
        const AKV_REAL wp = grid->w[p];

        AKV_REAL Hij[3][3], Nij[3][3], Jpi[3];
        eval_l1(p, grid, Hij, Nij, Jpi);

        mat3_add_inplace(Hloc, Hij, wp);
        mat3_add_inplace(Nloc, Nij, wp);
        vec3_add_inplace(Jloc, Jpi, wp);
      }
    }

#pragma omp critical
    {
      mat3_add_inplace(H_acc, Hloc, 1.0);
      mat3_add_inplace(N_acc, Nloc, 1.0);
      vec3_add_inplace(J_acc, Jloc, 1.0);
    }
  }
  memcpy(H, H_acc, 9*sizeof(AKV_REAL));
  memcpy(N, N_acc, 9*sizeof(AKV_REAL));
  memcpy(Jm, J_acc, 3*sizeof(AKV_REAL));
#else
  for (int it = grid->NG_theta; it < grid->NG_theta + grid->N_theta; it++) {
    for (int ip = grid->NG_phi; ip < grid->NG_phi + grid->N_phi; ip++) {
      const int p = akv_idx2_full(grid, it, ip);
      const AKV_REAL wp = grid->w[p];

      AKV_REAL Hij[3][3], Nij[3][3], Jpi[3];
      eval_l1(p, grid, Hij, Nij, Jpi);

      mat3_add_inplace(H, Hij, wp);
      mat3_add_inplace(N, Nij, wp);
      vec3_add_inplace(Jm, Jpi, wp);
    }
  }
#endif

  for (int i = 0; i < 3; i++) for (int j = i+1; j < 3; j++) {
    H[i][j] = H[j][i] = 0.5*(H[i][j] + H[j][i]);
    N[i][j] = N[j][i] = 0.5*(N[i][j] + N[j][i]);
  }

  if (pars->reg_eps > 0) {
    N[0][0] += pars->reg_eps;
    N[1][1] += pars->reg_eps;
    N[2][2] += pars->reg_eps;
  }

  AKV_REAL L[3][3];
  if (!chol3_lower(N, L)) {
    out->akv_quality_flag |= (1<<3);
    return AKV_ERR_EIGEN_FAIL;
  }

  AKV_REAL C[3][3];
  generalized_to_standard3(H, L, C);

  AKV_REAL eval[3], V[3][3];
  jacobi_eig_sym3(C, eval, V);
  sort_eigs3(eval, V);

  out->akv_lambda[0] = eval[0];
  out->akv_lambda[1] = eval[1];
  out->akv_lambda[2] = eval[2];

  for (int k = 0; k < 3; k++) {
    AKV_REAL u[3] = { V[0][k], V[1][k], V[2][k] };
    AKV_REAL x[3];
    tri3_solve_upper_from_lowerT(L, u, x);
    normalize_gen_evec3(N, x);

    AKV_REAL sref = 0.0;
    if (eval_sign_l1) sref = eval_sign_l1(x, Jm, grid);
    else sref = default_sign_ref_l1(x, Jm);
    if (sref < 0) { x[0] = -x[0]; x[1] = -x[1]; x[2] = -x[2]; }

    out->akv_l1_evecs[0][k] = x[0];
    out->akv_l1_evecs[1][k] = x[1];
    out->akv_l1_evecs[2][k] = x[2];

    AKV_REAL Hx[3], Nx[3], r[3];
    mat3_mul_vec(H, x, Hx);
    mat3_mul_vec(N, x, Nx);
    r[0] = Hx[0] - eval[k]*Nx[0];
    r[1] = Hx[1] - eval[k]*Nx[1];
    r[2] = Hx[2] - eval[k]*Nx[2];

    AKV_REAL nr = sqrt(vec3_dot(r, r));
    AKV_REAL nHx = sqrt(vec3_dot(Hx, Hx));
    AKV_REAL nNx = sqrt(vec3_dot(Nx, Nx));
    AKV_REAL denom = nHx + nNx;
    out->akv_eig_resid[k] = (denom > 0) ? (nr / denom) : 0.0;
    if (pars->eig_tol > 0 && out->akv_eig_resid[k] > pars->eig_tol) out->akv_quality_flag |= (1<<0);

    out->akv_J[k] = vec3_dot(x, Jm) / (8.0 * (AKV_REAL)M_PI);
  }
  akv_compute_dimensionless_spins(out->akv_J, pars->horizon_area, out->akv_a);

  if (pars->build_spin_vector && map_to_spinvec) {
    int k = pars->l1_choose_index;
    if (k < 0) k = 0;
    if (k > 2) k = 2;
    AKV_REAL coeffs[3] = { out->akv_l1_evecs[0][k], out->akv_l1_evecs[1][k], out->akv_l1_evecs[2][k] };
    map_to_spinvec(coeffs, grid, out->akv_spin_vec);
  } else {
    out->akv_spin_vec[0] = 0.0;
    out->akv_spin_vec[1] = 0.0;
    out->akv_spin_vec[2] = 0.0;
  }

  out->akv_eig_gap_43 = 0.0;
  out->method_used = AKV_METHOD_L1_REDUCED;

  return AKV_SUCCESS;
}

// -----------------------------
// Gridpoint dense assembly
// -----------------------------

typedef struct {
  int j;
  AKV_REAL Kv;
  AKV_REAL Bv;
} akv_triplet_t;

static int cmp_triplet_j(const void *a, const void *b) {
  const akv_triplet_t *x = (const akv_triplet_t*)a;
  const akv_triplet_t *y = (const akv_triplet_t*)b;
  return (x->j < y->j) ? -1 : (x->j > y->j);
}

// Assemble dense matrices from stencil/bilinear-form callback.
// For each row i, the callback provides (i,j,K,B) entries; duplicates in j are summed.
static akv_error_t akv_full_assemble_dense_from_stencil(
    const akv_horizon_grid_t *grid,
    akv_eval_row_stencil_f eval_row,
    AKV_REAL *Kmat, AKV_REAL *Bmat, int N_full) {

  if (!grid || !eval_row || !Kmat || !Bmat) return AKV_ERR_NULLPTR;

  const int max_entries = 512;
  int *j_cols = (int*)malloc((size_t)max_entries * sizeof(int));
  AKV_REAL *K_vals = (AKV_REAL*)malloc((size_t)max_entries * sizeof(AKV_REAL));
  AKV_REAL *B_vals = (AKV_REAL*)malloc((size_t)max_entries * sizeof(AKV_REAL));
  akv_triplet_t *t = (akv_triplet_t*)malloc((size_t)max_entries * sizeof(akv_triplet_t));
  if (!j_cols || !K_vals || !B_vals || !t) {
    free(j_cols); free(K_vals); free(B_vals); free(t);
    return AKV_ERR_ALLOC;
  }

  for (int i = 0; i < N_full; i++) {
    int nnz = eval_row(i, grid, max_entries, j_cols, K_vals, B_vals);
    if (nnz < 0) { free(j_cols); free(K_vals); free(B_vals); free(t); return AKV_ERR_BADPARAM; }
    if (nnz > max_entries) nnz = max_entries;

    int m = 0;
    for (int e = 0; e < nnz; e++) {
      int j = j_cols[e];
      if (j < 0 || j >= N_full) continue;
      t[m].j = j;
      t[m].Kv = K_vals[e];
      t[m].Bv = B_vals[e];
      m++;
    }

    if (m == 0) continue;

    qsort(t, (size_t)m, sizeof(akv_triplet_t), cmp_triplet_j);

    int curj = t[0].j;
    AKV_REAL sumK = t[0].Kv;
    AKV_REAL sumB = t[0].Bv;

    for (int e = 1; e < m; e++) {
      if (t[e].j == curj) {
        sumK += t[e].Kv;
        sumB += t[e].Bv;
      } else {
        Kmat[i + (size_t)N_full * curj] = sumK;
        Bmat[i + (size_t)N_full * curj] = sumB;
        curj = t[e].j;
        sumK = t[e].Kv;
        sumB = t[e].Bv;
      }
    }
    Kmat[i + (size_t)N_full * curj] = sumK;
    Bmat[i + (size_t)N_full * curj] = sumB;
  }

  free(j_cols); free(K_vals); free(B_vals); free(t);

  akv_symmetrize_dense(Kmat, N_full);
  akv_symmetrize_dense(Bmat, N_full);
  return AKV_SUCCESS;
}

// Debug-only dense assembly by applying operator to basis vectors (O(N^3) operator work).
static akv_error_t akv_full_assemble_dense_debug_apply(
    const akv_horizon_grid_t *grid,
    akv_apply_K_f applyK,
    akv_apply_B_f applyB,
    AKV_REAL *Kmat, AKV_REAL *Bmat, int N_full) {

  if (!grid || !applyK || !applyB || !Kmat || !Bmat) return AKV_ERR_NULLPTR;

  AKV_REAL *e = (AKV_REAL*)calloc((size_t)N_full, sizeof(AKV_REAL));
  AKV_REAL *Kz = (AKV_REAL*)calloc((size_t)N_full, sizeof(AKV_REAL));
  AKV_REAL *Bz = (AKV_REAL*)calloc((size_t)N_full, sizeof(AKV_REAL));
  if (!e || !Kz || !Bz) {
    free(e); free(Kz); free(Bz);
    return AKV_ERR_ALLOC;
  }

  for (int n = 0; n < N_full; n++) {
    memset(e, 0, (size_t)N_full * sizeof(AKV_REAL));
    e[n] = 1.0;
    applyK(grid, e, Kz);
    applyB(grid, e, Bz);
    for (int m = 0; m < N_full; m++) {
      Kmat[m + (size_t)N_full * n] = Kz[m];
      Bmat[m + (size_t)N_full * n] = Bz[m];
    }
  }

  akv_symmetrize_dense(Kmat, N_full);
  akv_symmetrize_dense(Bmat, N_full);

  free(e); free(Kz); free(Bz);
  return AKV_SUCCESS;
}

// -----------------------------
// Gridpoint dense generalized eigen solve with explicit SPD handling
// -----------------------------

#ifdef USE_LAPACKE
static bool akv_try_cholesky_spd(const AKV_REAL *B, int N) {
  AKV_REAL *tmp = (AKV_REAL*)malloc((size_t)N*(size_t)N*sizeof(AKV_REAL));
  if (!tmp) return false;
  memcpy(tmp, B, (size_t)N*(size_t)N*sizeof(AKV_REAL));
  int info = LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'U', N, tmp, N);
  free(tmp);
  return (info == 0);
}
#endif

static akv_error_t akv_full_solve_dense_generalized_reduced_spd_retry(
    AKV_REAL *Kred, AKV_REAL *Bred, int Nr, int want_neigs,
    const akv_params_t *pars,
    AKV_REAL *eval_out, AKV_REAL *evec_out,
    int *quality_flag_io) {

  if (!Kred || !Bred || !eval_out || !evec_out || !pars || !quality_flag_io) return AKV_ERR_NULLPTR;
  if (Nr <= 0 || want_neigs <= 0 || want_neigs > Nr) return AKV_ERR_BADPARAM;

#ifndef USE_LAPACKE
  (void)Kred; (void)Bred; (void)Nr; (void)want_neigs; (void)pars; (void)eval_out; (void)evec_out;
  *quality_flag_io |= (1<<2);
  return AKV_ERR_NOT_IMPLEMENTED;
#else
  const int max_tries = (pars->reg_max_tries > 0) ? pars->reg_max_tries : 1;
  AKV_REAL reg = (pars->reg_eps > 0) ? pars->reg_eps : 0.0;
  const AKV_REAL reg_max = (pars->reg_eps_max > reg) ? pars->reg_eps_max : reg;

  bool max_eps = false;
  for (int attempt = 0; attempt < max_tries; attempt++) {
    AKV_REAL *A = (AKV_REAL*)malloc((size_t)Nr*(size_t)Nr*sizeof(AKV_REAL));
    AKV_REAL *B = (AKV_REAL*)malloc((size_t)Nr*(size_t)Nr*sizeof(AKV_REAL));
    AKV_REAL *w = (AKV_REAL*)malloc((size_t)Nr*sizeof(AKV_REAL));
    if (!A || !B || !w) { free(A); free(B); free(w); return AKV_ERR_ALLOC; }

    memcpy(A, Kred, (size_t)Nr*(size_t)Nr*sizeof(AKV_REAL));
    memcpy(B, Bred, (size_t)Nr*(size_t)Nr*sizeof(AKV_REAL));

    if (reg > 0) {
      akv_dense_add_diag(B, Nr, reg);
      *quality_flag_io |= (1<<3);
    }

    // Explicit SPD probe (Cholesky). If it fails, increase reg and retry.
    if (!akv_try_cholesky_spd(B, Nr)) {
      free(A); free(B); free(w);
      *quality_flag_io |= (1<<3);
      if (attempt == max_tries - 1 || max_eps) return AKV_ERR_EIGEN_FAIL;
      if (reg == 0.0) reg = 1e-14;
      else reg *= 10.0;
      if (reg_max > 0 && reg > reg_max) { reg = reg_max; max_eps = true; }
      continue;
    }
    
    // Solve A x = λ B x. LAPACKE_dsygvd overwrites A with eigenvectors.
    int info = LAPACKE_dsygvd(LAPACK_COL_MAJOR, 1, 'V', 'U', Nr, A, Nr, B, Nr, w);
    if (info != 0) {
      free(A); free(B); free(w);
      *quality_flag_io |= (1<<3);
      if (attempt == max_tries - 1 || max_eps) return AKV_ERR_EIGEN_FAIL;
      if (reg == 0.0) reg = 1e-14;
      else reg *= 10.0;
      if (reg_max > 0 && reg > reg_max) { reg = reg_max; max_eps = true; }
      continue;
    }

    for (int i = 0; i < want_neigs; i++) {
      eval_out[i] = w[i];
      for (int r = 0; r < Nr; r++) evec_out[r + (size_t)Nr*i] = A[r + (size_t)Nr*i];
    }

    free(A); free(B); free(w);
    return AKV_SUCCESS;
  }

  return AKV_ERR_EIGEN_FAIL;
#endif
}

// -----------------------------
// Gridpoint solve + postprocess + eigenvector storage + J loop scaffold
// -----------------------------
static akv_error_t akv_solve_gridpoint_basis(
    const akv_horizon_grid_t *grid,
    const akv_params_t *pars,
    akv_eval_row_stencil_f eval_row_stencil,
    akv_apply_K_f applyK_debug,
    akv_apply_B_f applyB_debug,
    akv_eval_sign_full_f eval_sign_full,
    akv_eval_J_integrand_full_f eval_J_full,
    akv_diagnostics_t *out) {

  if (!grid || !pars || !out) return AKV_ERR_NULLPTR;
  if (!grid->w) return AKV_ERR_NULLPTR;

  const int N_full = grid->N_theta * grid->N_phi;
  if (N_full <= 1) return AKV_ERR_BADPARAM;

  akv_mean0_basis_t b;
  akv_error_t err = akv_mean0_basis_build(grid, &b);
  if (err != AKV_SUCCESS) return err;

  const int Nr = b.N_red;

  // Dense memory warning: two N_full^2 matrices may be large at typical resolutions.
  AKV_REAL *Kfull = (AKV_REAL*)calloc((size_t)N_full*(size_t)N_full, sizeof(AKV_REAL));
  AKV_REAL *Bfull = (AKV_REAL*)calloc((size_t)N_full*(size_t)N_full, sizeof(AKV_REAL));
  if (!Kfull || !Bfull) {
    free(Kfull); free(Bfull);
    akv_mean0_basis_free(&b);
    return AKV_ERR_ALLOC;
  }

  if (eval_row_stencil) {
    err = akv_full_assemble_dense_from_stencil(grid, eval_row_stencil, Kfull, Bfull, N_full);
    if (err != AKV_SUCCESS) {
      free(Kfull); free(Bfull);
      akv_mean0_basis_free(&b);
      return err;
    }
  } else {
    if (!pars->allow_debug_assembly || !applyK_debug || !applyB_debug) {
      out->akv_quality_flag |= (1<<2);
      free(Kfull); free(Bfull);
      akv_mean0_basis_free(&b);
      return AKV_ERR_BADPARAM;
    }
    err = akv_full_assemble_dense_debug_apply(grid, applyK_debug, applyB_debug, Kfull, Bfull, N_full);
    if (err != AKV_SUCCESS) {
      free(Kfull); free(Bfull);
      akv_mean0_basis_free(&b);
      return err;
    }
  }

  AKV_REAL *Kred = (AKV_REAL*)calloc((size_t)Nr*(size_t)Nr, sizeof(AKV_REAL));
  AKV_REAL *Bred = (AKV_REAL*)calloc((size_t)Nr*(size_t)Nr, sizeof(AKV_REAL));
  if (!Kred || !Bred) {
    free(Kfull); free(Bfull); free(Kred); free(Bred);
    akv_mean0_basis_free(&b);
    return AKV_ERR_ALLOC;
  }

  akv_mean0_transform_dense(&b, Kfull, N_full, Kred, Nr);
  akv_mean0_transform_dense(&b, Bfull, N_full, Bred, Nr);
  akv_symmetrize_dense(Kred, Nr);
  akv_symmetrize_dense(Bred, Nr);

  int want = pars->full_num_eigs;
  if (want < 4) want = 4;
  if (want > Nr) want = Nr;

  AKV_REAL *eval = (AKV_REAL*)malloc((size_t)want * sizeof(AKV_REAL));
  AKV_REAL *evec = (AKV_REAL*)malloc((size_t)Nr*(size_t)want * sizeof(AKV_REAL));
  if (!eval || !evec) {
    free(Kfull); free(Bfull); free(Kred); free(Bred); free(eval); free(evec);
    akv_mean0_basis_free(&b);
    return AKV_ERR_ALLOC;
  }

  err = akv_full_solve_dense_generalized_reduced_spd_retry(Kred, Bred, Nr, want, pars, eval, evec, &out->akv_quality_flag);
  if (err != AKV_SUCCESS) {
    if (err == AKV_ERR_NOT_IMPLEMENTED) out->akv_quality_flag |= (1<<2);
    free(Kfull); free(Bfull); free(Kred); free(Bred); free(eval); free(evec);
    akv_mean0_basis_free(&b);
    return err;
  }

  // B-orthonormalize first three modes using Bred.
  const int kmodes = (want >= 3) ? 3 : want;
  akv_b_gram_schmidt_inplace(evec, Nr, kmodes, Bred, b.tmpN);

  // Gap diagnostic R_{3,4} = λ4/λ3 in mean-zero subspace.
  if (want >= 4 && eval[2] != 0.0) out->akv_eig_gap_43 = eval[3] / eval[2];
  else out->akv_eig_gap_43 = 0.0;

  if (pars->gap_ratio_thresh > 0 && want >= 4 && out->akv_eig_gap_43 > 0 &&
      out->akv_eig_gap_43 < pars->gap_ratio_thresh) {
    out->akv_quality_flag |= (1<<1);
  }

  out->akv_lambda[0] = (want > 0 ? eval[0] : 0.0);
  out->akv_lambda[1] = (want > 1 ? eval[1] : 0.0);
  out->akv_lambda[2] = (want > 2 ? eval[2] : 0.0);

  // Residual ratios for first three modes using dense reduced matrices and reused buffers.
  for (int a = 0; a < 3; a++) {
    if (a >= want) { out->akv_eig_resid[a] = 0.0; continue; }
    const AKV_REAL *za = &evec[(size_t)Nr * a];
    out->akv_eig_resid[a] = akv_residual_ratio_dense_inplace(Kred, Bred, Nr, za, eval[a], b.tmpK, b.tmpN, b.tmpR);
    if (pars->eig_tol > 0 && out->akv_eig_resid[a] > pars->eig_tol) out->akv_quality_flag |= (1<<0);
  }

  // Store full-space eigenvectors (interior DOF space) for downstream postprocessing.
  out->gp_N = N_full;
  for (int a = 0; a < 3; a++) {
    out->gp_z[a] = NULL;
    if (a >= kmodes) continue;
    out->gp_z[a] = (AKV_REAL*)malloc((size_t)N_full * sizeof(AKV_REAL));
    if (!out->gp_z[a]) {
      for (int k = 0; k < a; k++) { free(out->gp_z[k]); out->gp_z[k] = NULL; }
      out->gp_N = 0;
      free(Kfull); free(Bfull); free(Kred); free(Bred); free(eval); free(evec);
      akv_mean0_basis_free(&b);
      return AKV_ERR_ALLOC;
    }
    akv_mean0_lift_vector(&b, &evec[(size_t)Nr * a], out->gp_z[a]);
  }

  // Sign-fix in full space (if provided). Flip stored gp_z and the reduced eigenvector consistently.
  for (int a = 0; a < kmodes; a++) {
    if (!out->gp_z[a]) continue;
    if (!eval_sign_full) continue;

    AKV_REAL sref = eval_sign_full(out->gp_z[a], N_full, grid);
    if (sref < 0) {
      // Flip reduced eigenvector
      for (int i = 0; i < Nr; i++) evec[i + (size_t)Nr * a] = -evec[i + (size_t)Nr * a];
      // Re-lift
      akv_mean0_lift_vector(&b, &evec[(size_t)Nr * a], out->gp_z[a]);
    }
  }

  // Loop scaffold for angular momentum integrals J[a] using per-mode z_full.
  // J integrand callback must not include 1/(8π). This function applies weights and sums.
  for (int a = 0; a < 3; a++) out->akv_J[a] = 0.0;
  if (eval_J_full) {
    for (int a = 0; a < kmodes; a++) {
      const AKV_REAL *z = out->gp_z[a];
      if (!z) continue;

      AKV_REAL sum = 0.0;
      for (int dof = 0; dof < N_full; dof++) {
        int p = akv_dof_to_full_p(grid, dof);
        AKV_REAL wp = grid->w[p];
        AKV_REAL integrand = eval_J_full(p, grid, z, N_full);
        sum += wp * integrand;
      }
      out->akv_J[a] = sum / (8.0 * (AKV_REAL)M_PI);
    }
  }

  // Dimensionless spin components normalized by M_H^2.
  akv_compute_dimensionless_spins(out->akv_J, pars->horizon_area, out->akv_a);

  out->akv_spin_vec[0] = 0.0;
  out->akv_spin_vec[1] = 0.0;
  out->akv_spin_vec[2] = 0.0;

  out->method_used = AKV_METHOD_GRIDPOINT;

  free(Kfull); free(Bfull);
  free(Kred); free(Bred);
  free(eval); free(evec);
  akv_mean0_basis_free(&b);
  return AKV_SUCCESS;
}

// -----------------------------
// Public entry point
// -----------------------------
akv_error_t akv_compute(
    const akv_horizon_grid_t *grid,
    const akv_params_t *pars,
    // L1-reduced callbacks
    akv_eval_l1_integrands_f eval_l1,
    akv_eval_sign_l1_f eval_sign_l1,
    akv_map_to_spinvec_f map_to_spinvec,
    // Gridpoint callbacks
    akv_eval_row_stencil_f eval_row_stencil,
    akv_apply_K_f applyK_debug,
    akv_apply_B_f applyB_debug,
    akv_eval_sign_full_f eval_sign_full,
    akv_eval_J_integrand_full_f eval_J_full,
    akv_diagnostics_t *out) {

  if (!grid || !pars || !out) return AKV_ERR_NULLPTR;
  memset(out, 0, sizeof(*out));

  switch (pars->method) {
    case AKV_METHOD_L1_REDUCED:
      return akv_solve_l1_reduced(grid, pars, eval_l1, eval_sign_l1, map_to_spinvec, out);
    case AKV_METHOD_GRIDPOINT:
      return akv_solve_gridpoint_basis(grid, pars, eval_row_stencil, applyK_debug, applyB_debug,
                                       eval_sign_full, eval_J_full, out);
    default:
      return AKV_ERR_BADPARAM;
  }
}

// End of scaffold


#endif /* BHAHAHA_BHAHAHA_NRPY_AKV_SCAFFOLD_INCLUDED */

#ifdef __GNUC__
#define AKV_WEAK __attribute__((weak))
#else
#define AKV_WEAK
#endif

/* Stub-used mask bits (8 bits). If BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION==1:
 *   bit  BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT      = any optional stub used
 *   bits [BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT .. +7] = packed 8-bit stub mask
 *   (all bit positions are configurable via BHAHAHA_AKV_QUAL_EXT_* macros)
 */
// Stub-usage bits (8-bit mask). Kept as explicit 32-bit values to make shifts/ORs portable.
// Note: bit0 (L1 integrands) is tracked internally but *excluded* from the packed mask/any-stub bit
// so that "any stub used" does not become trivially true for all runs.
#define BHAHAHA_AKV_STUB_USED_L1_INTEGRANDS (UINT32_C(1) << 0)  /* wrapper-provided L1 integrands used */
#define BHAHAHA_AKV_STUB_USED_SIGN_L1       (UINT32_C(1) << 1)
#define BHAHAHA_AKV_STUB_USED_SIGN_FULL     (UINT32_C(1) << 2)
#define BHAHAHA_AKV_STUB_USED_J_FULL        (UINT32_C(1) << 3)
#define BHAHAHA_AKV_STUB_USED_MAP_SPINVEC   (UINT32_C(1) << 4)
#define BHAHAHA_AKV_STUB_USED_ROW_STENCIL   (UINT32_C(1) << 5)
#define BHAHAHA_AKV_STUB_USED_APPLYK        (UINT32_C(1) << 6)
#define BHAHAHA_AKV_STUB_USED_APPLYB        (UINT32_C(1) << 7)

typedef struct {
  commondata_struct *restrict commondata;
  griddata_struct   *restrict griddata;
  int grid;

  // Pointers needed by the generated AKV integrands callback:
  const params_struct *restrict params;
  REAL *restrict xx[3];
  REAL *restrict auxevol_gfs;
  const REAL *restrict in_gfs;

  uint32_t *stub_flags_ptr; /* per-call flags sink */
} bhahaha_akv_context_t;

static inline void bhahaha_akv_mark_stub(const akv_horizon_grid_t *grid, uint32_t bit) {
  /* Prefer per-call accumulation via ctx->stub_flags_ptr. We intentionally do not use
     any global/TLS fallback here to keep the behavior unambiguous and race-free.
     NOTE: stub tracking therefore requires g.geom to be set to a valid ctx. */
  if (!grid || !grid->geom) { (void)bit; return; }
  const bhahaha_akv_context_t *ctx = (const bhahaha_akv_context_t*)grid->geom;
  if (!ctx || !ctx->stub_flags_ptr) { (void)bit; return; }

  {
#ifdef _OPENMP
#pragma omp critical(bhahaha_akv_spin_stubflags_or)
#endif
    { *(ctx->stub_flags_ptr) |= bit; }
  }
}


/* -------------------- Weak stub definitions -------------------- */
static void bhahaha_akv_eval_l1_integrands_default_impl(
    int p, const akv_horizon_grid_t *grid,
    AKV_REAL Hmn[3][3], AKV_REAL Nmn[3][3], AKV_REAL Jm[3]) {
  if (!grid || !grid->geom || !Hmn || !Nmn || !Jm) return;

  const bhahaha_akv_context_t *ctx = (const bhahaha_akv_context_t*)grid->geom;
  const commondata_struct *restrict commondata = ctx->commondata;
  const params_struct *restrict params = ctx->params;
  REAL *restrict auxevol_gfs = ctx->auxevol_gfs;
  const REAL *restrict in_gfs = ctx->in_gfs;

  #include "set_CodeParameters.h"

  // p is the full 2D ghosted horizon-grid index: p = it_full*N_phi_tot + ip_full.
  const int N_phi_tot = grid->N_phi_tot;
  const int ip_full = p % N_phi_tot;            // phi index (full, incl. ghosts)
  const int it_full = (p - ip_full) / N_phi_tot;// theta index (full, incl. ghosts)
  const int i2 = ip_full;                       // AH phi index (full, incl. ghosts)
  const int i1 = it_full;                       // AH theta index (full, incl. ghosts)
  const int i0 = BHAHAHA_AKV_AH_SLAB_I0;         // AH slab index in the 3D grid

  const REAL xx0 = ctx->xx[0][i0];
  const REAL xx1 = ctx->xx[1][i1];
  const REAL xx2 = ctx->xx[2][i2];

  // Human-readable SymPy expressions, CSE'd and FD-optimized by NRPy:
  """
        + akv_codegen
        + r"""

#if (BHAHAHA_AKV_CALLBACK_DENSITY_FREE) && (BHAHAHA_AKV_EXPR_ARE_DENSITIES)
  // Return density-free integrands if expressions include sqrt(q) and w2d includes sqrt(det q).
  const AKV_REAL inv_sqrtq = (sqrt_det_q_cb != 0.0) ? (1.0 / sqrt_det_q_cb) : 0.0;
#else
  const AKV_REAL inv_sqrtq = 1.0;
#endif

  // Fill outputs. Generated symbols are densities (include sqrt(q)).
  Hmn[0][0] = akv_Hraw_0_0 * inv_sqrtq;
  Hmn[0][1] = akv_Hraw_0_1 * inv_sqrtq;
  Hmn[0][2] = akv_Hraw_0_2 * inv_sqrtq;
  Hmn[1][0] = akv_Hraw_1_0 * inv_sqrtq;
  Hmn[1][1] = akv_Hraw_1_1 * inv_sqrtq;
  Hmn[1][2] = akv_Hraw_1_2 * inv_sqrtq;
  Hmn[2][0] = akv_Hraw_2_0 * inv_sqrtq;
  Hmn[2][1] = akv_Hraw_2_1 * inv_sqrtq;
  Hmn[2][2] = akv_Hraw_2_2 * inv_sqrtq;

  Nmn[0][0] = akv_Nraw_0_0 * inv_sqrtq;
  Nmn[0][1] = akv_Nraw_0_1 * inv_sqrtq;
  Nmn[0][2] = akv_Nraw_0_2 * inv_sqrtq;
  Nmn[1][0] = akv_Nraw_1_0 * inv_sqrtq;
  Nmn[1][1] = akv_Nraw_1_1 * inv_sqrtq;
  Nmn[1][2] = akv_Nraw_1_2 * inv_sqrtq;
  Nmn[2][0] = akv_Nraw_2_0 * inv_sqrtq;
  Nmn[2][1] = akv_Nraw_2_1 * inv_sqrtq;
  Nmn[2][2] = akv_Nraw_2_2 * inv_sqrtq;

  // J integrand components without 1/(8π).
  Jm[0] = akv_Jraw_0 * inv_sqrtq;
  Jm[1] = akv_Jraw_1 * inv_sqrtq;
  Jm[2] = akv_Jraw_2 * inv_sqrtq;

  // Ensure symmetry at the callback boundary (robustness; harmless if already symmetric).
  for (int m = 0; m < 3; m++) for (int n = m+1; n < 3; n++) {
    const AKV_REAL hs = 0.5*(Hmn[m][n] + Hmn[n][m]);
    const AKV_REAL ns = 0.5*(Nmn[m][n] + Nmn[n][m]);
    Hmn[m][n] = Hmn[n][m] = hs;
    Nmn[m][n] = Nmn[n][m] = ns;
  }

  (void)commondata; // available for future extensions
}

AKV_WEAK void bhahaha_akv_eval_l1_integrands(
    int p, const akv_horizon_grid_t *grid,
    AKV_REAL Hmn[3][3], AKV_REAL Nmn[3][3], AKV_REAL Jm[3]) {
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_L1_INTEGRANDS);
  bhahaha_akv_eval_l1_integrands_default_impl(p, grid, Hmn, Nmn, Jm);
}

AKV_WEAK AKV_REAL bhahaha_akv_eval_sign_l1(
    const AKV_REAL coeffs[3], const AKV_REAL Jm_assembled[3], const akv_horizon_grid_t *grid) {
  (void)grid;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_SIGN_L1);
  return coeffs[0]*Jm_assembled[0] + coeffs[1]*Jm_assembled[1] + coeffs[2]*Jm_assembled[2];
}

AKV_WEAK void bhahaha_akv_map_to_spinvec(
    const AKV_REAL coeffs[3],
    const akv_horizon_grid_t *grid,
    AKV_REAL S_out[3]) {
  (void)coeffs; (void)grid;
  if (!S_out) return;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_MAP_SPINVEC);
  S_out[0] = 0.0; S_out[1] = 0.0; S_out[2] = 0.0;
}

AKV_WEAK int bhahaha_akv_eval_row_stencil(
    int i_row, const akv_horizon_grid_t *grid, int max_entries,
    int *j_cols, AKV_REAL *K_vals, AKV_REAL *B_vals) {
  (void)grid;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_ROW_STENCIL);
  /* Gridpoint method is not implemented in this wrapper. Returning -1 entries makes failure
   * mode more obvious than silently emitting an identity row.
   */
  (void)i_row;
  (void)max_entries;
  (void)j_cols;
  (void)K_vals;
  (void)B_vals;
  return -1;
}

/* IMPORTANT CONTRACT (AKV scaffold): apply_K/apply_B act on interior-DOF vectors of length N_theta*N_phi. */
AKV_WEAK void bhahaha_akv_apply_K(const akv_horizon_grid_t *grid, const AKV_REAL *z_in, AKV_REAL *Kz_out) {
  (void)z_in;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_APPLYK);
  if (!grid || !Kz_out) return;
  const int N_full = grid->N_theta * grid->N_phi;
  memset(Kz_out, 0, (size_t)N_full * sizeof(AKV_REAL));
}

AKV_WEAK void bhahaha_akv_apply_B(const akv_horizon_grid_t *grid, const AKV_REAL *z_in, AKV_REAL *Bz_out) {
  (void)z_in;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_APPLYB);
  if (!grid || !Bz_out) return;
  const int N_full = grid->N_theta * grid->N_phi;
  memset(Bz_out, 0, (size_t)N_full * sizeof(AKV_REAL));
}

AKV_WEAK AKV_REAL bhahaha_akv_eval_sign_full(
    const AKV_REAL *z_full, int N_full, const akv_horizon_grid_t *grid) {
  (void)z_full; (void)N_full;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_SIGN_FULL);
  return 1.0;
}

AKV_WEAK AKV_REAL bhahaha_akv_eval_J_integrand_full(
    int p, const akv_horizon_grid_t *grid, const AKV_REAL *z_full, int N_full) {
  (void)p; (void)z_full; (void)N_full;
  bhahaha_akv_mark_stub(grid, BHAHAHA_AKV_STUB_USED_J_FULL);
  return 0.0;
}
"""
    )

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "BHaHAHA apparent horizon diagnostics: Approximate Killing Vector (AKV) spin diagnostics."
    cfunc_type = "int"
    name = "diagnostics_approx_killing_vector_spin"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    body = r"""
  const int grid = 0;
  bhahaha_diagnostics_struct *restrict bhahaha_diags = commondata->bhahaha_diagnostics;
  const params_struct *restrict params = &griddata[grid].params;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
  const REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  (void)auxevol_gfs; (void)in_gfs;

  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++) xx[ww] = griddata[grid].xx[ww];
#include "set_CodeParameters.h"

  const REAL *restrict weights_1d;
  int weight_stencil_size;
  bah_diagnostics_integration_weights(Nxx1, Nxx2, &weights_1d, &weight_stencil_size);
  /* NOTE: bah_diagnostics_integration_weights() returns a single 1D repeating stencil table
   * (length = weight_stencil_size) intended to be tiled across BOTH horizon-grid directions.
   * See bah_diagnostics_integration_weights.c and existing diagnostics (e.g., proper circumferences)
   * for the canonical usage pattern: weights_1d[i % weight_stencil_size].
   * The returned coefficients are dimensionless; multiply by coordinate spacings separately if desired.
   */
  if (weight_stencil_size <= 0 || weights_1d == NULL) return (int)AKV_ERR_BADPARAM;

  const int NGt = NGHOSTS;
  const int NGp = NGHOSTS;
  const int N_theta = Nxx1;
  const int N_phi   = Nxx2;
  const int N_theta_tot = N_theta + 2*NGt;
  const int N_phi_tot   = N_phi   + 2*NGp;
  const int N_tot_2d = N_theta_tot * N_phi_tot;

  REAL *restrict w2d = NULL;

#if BHAHAHA_HAVE_TLS
  static BHAHAHA_TLS REAL *w2d_cache = NULL;
  static BHAHAHA_TLS size_t w2d_cache_len = 0;
  if (w2d_cache_len < (size_t)N_tot_2d) {
    REAL *newptr = (REAL*)realloc(w2d_cache, sizeof(REAL) * (size_t)N_tot_2d);
    if (!newptr) return (int)AKV_ERR_ALLOC;
    w2d_cache = newptr;
    w2d_cache_len = (size_t)N_tot_2d;
  }
  w2d = w2d_cache;
#else
  /* No TLS: disable caching to avoid shared static pointer races. */
  w2d = (REAL*)malloc(sizeof(REAL) * (size_t)N_tot_2d);
  if (!w2d) return (int)AKV_ERR_ALLOC;
#endif

  for (int p = 0; p < N_tot_2d; p++) w2d[p] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
  for (int ip_full = NGp; ip_full < NGp + N_phi; ip_full++) {
    for (int it_full = NGt; it_full < NGt + N_theta; it_full++) {
      const int i2 = ip_full;
      const int i1 = it_full;
      const int i0 = BHAHAHA_AKV_AH_SLAB_I0;

      // Map physical-node indices (i1,i2) to 1D quadrature weights.
      // Contract: bah_diagnostics_integration_weights returns a single 1D *stencil table*
      // of length weight_stencil_size that is tiled across both theta and phi directions.
      // Therefore we apply the same periodic tiling in both directions using modulo
      // weight_stencil_size, with positive wrap.
      int it_w = i1 - NGHOSTS;
      int ip_w = i2 - NGHOSTS;
      it_w = (it_w % weight_stencil_size + weight_stencil_size) % weight_stencil_size;
      ip_w = (ip_w % weight_stencil_size + weight_stencil_size) % weight_stencil_size;

      const REAL wtheta = weights_1d[it_w];
      const REAL wphi   = weights_1d[ip_w];


"""
    body += area_codegen + r"""
      const int p2d = it_full * N_phi_tot + ip_full;

#if BHAHAHA_AKV_WEIGHTS_INCLUDE_DX
      w2d[p2d] = sqrt_det_q * wtheta * wphi * params->dxx1 * params->dxx2;
#else
      w2d[p2d] = sqrt_det_q * wtheta * wphi;
#endif
    }
  }

  akv_horizon_grid_t g;
  g.N_theta = N_theta;
  g.N_phi   = N_phi;
  g.NG_theta = NGt;
  g.NG_phi   = NGp;
  g.dtheta = params->dxx1;
  g.dphi   = params->dxx2;
  g.N_theta_tot = N_theta_tot;
  g.N_phi_tot   = N_phi_tot;
  g.N_tot       = N_tot_2d;
  g.w = w2d;

  uint32_t stub_flags = 0u;

  // Inputs required by the L1 integrands callback: reuse params/in_gfs/xx defined above.

  bhahaha_akv_context_t ctx;
  ctx.commondata = commondata;
  ctx.griddata   = griddata;
  ctx.grid       = grid;
  ctx.params     = params;
  ctx.xx[0]      = xx[0];
  ctx.xx[1]      = xx[1];
  ctx.xx[2]      = xx[2];
  ctx.auxevol_gfs = auxevol_gfs;
  ctx.in_gfs     = in_gfs;
  ctx.stub_flags_ptr = &stub_flags;

  g.geom = (const void *)&ctx;

  // Compute horizon area A = sum_p w^{(p)} over interior points:
  REAL A_H = 0.0;
  for (int it_full = NGt; it_full < NGt + N_theta; it_full++) {
    for (int ip_full = NGp; ip_full < NGp + N_phi; ip_full++) {
      const int p2d = it_full * N_phi_tot + ip_full;
      A_H += w2d[p2d];
    }
  }

  akv_params_t pars;
  memset(&pars, 0, sizeof(pars));
  pars.method = AKV_METHOD_L1_REDUCED;
  pars.l1_choose_index = 0;
  pars.full_num_eigs = 8;
  pars.eig_tol = 1e-10;
  pars.reg_eps = 1e-14;
  pars.reg_eps_max = 1e-6;
  pars.reg_max_tries = 6;
  pars.gap_ratio_thresh = 1.5;
  pars.horizon_area = (AKV_REAL)A_H;
  pars.build_spin_vector = (BHAHAHA_AKV_ENABLE_SPIN_VECTOR ? true : false);
  pars.allow_debug_assembly = false;

#ifdef BHAHAHA_AKV_DIAGS_HAVE_PARAMS
  pars.method = (akv_method_t)bhahaha_diags->akv_method;
  pars.l1_choose_index = bhahaha_diags->akv_l1_choose_index;
  pars.full_num_eigs = bhahaha_diags->akv_full_num_eigs;
  pars.eig_tol = bhahaha_diags->akv_eig_tol;
  pars.reg_eps = bhahaha_diags->akv_reg_eps;
  pars.reg_eps_max = bhahaha_diags->akv_reg_eps_max;
  pars.reg_max_tries = bhahaha_diags->akv_reg_max_tries;
  pars.gap_ratio_thresh = bhahaha_diags->akv_gap_ratio_thresh;
  pars.build_spin_vector = (BHAHAHA_AKV_ENABLE_SPIN_VECTOR ? true : false);
  pars.allow_debug_assembly = (bhahaha_diags->akv_allow_debug_assembly != 0);
#endif


#if !BHAHAHA_AKV_ENABLE_GRIDPOINT
  if (pars.method == AKV_METHOD_GRIDPOINT) {
#ifdef BHAHAHA_AKV_DIAGS_HAVE_FIELDS
    bhahaha_diags->akv_quality_flag |= (1<<2);
    bhahaha_diags->akv_method_used = (int)AKV_METHOD_GRIDPOINT;
#endif
#if !BHAHAHA_HAVE_TLS
    free(w2d);
#endif
    return (int)AKV_ERR_NOT_IMPLEMENTED;
  }
#endif

#ifndef USE_LAPACKE
  if (pars.method == AKV_METHOD_GRIDPOINT) {
#ifdef BHAHAHA_AKV_DIAGS_HAVE_FIELDS
    bhahaha_diags->akv_quality_flag |= (1<<2);
    bhahaha_diags->akv_method_used = (int)AKV_METHOD_GRIDPOINT;
#endif
#if !BHAHAHA_HAVE_TLS
    free(w2d);
#endif
    return (int)AKV_ERR_NOT_IMPLEMENTED;
  }
#endif

  akv_diagnostics_t out;
  const akv_error_t err = akv_compute(&g, &pars,
                                      bhahaha_akv_eval_l1_integrands,
                                      bhahaha_akv_eval_sign_l1,
                                      (BHAHAHA_AKV_ENABLE_SPIN_VECTOR ? bhahaha_akv_map_to_spinvec : NULL),
                                      (BHAHAHA_AKV_ENABLE_GRIDPOINT ? bhahaha_akv_eval_row_stencil : NULL),
                                      (BHAHAHA_AKV_ENABLE_GRIDPOINT ? bhahaha_akv_apply_K : NULL),
                                      (BHAHAHA_AKV_ENABLE_GRIDPOINT ? bhahaha_akv_apply_B : NULL),
                                      (BHAHAHA_AKV_ENABLE_GRIDPOINT ? bhahaha_akv_eval_sign_full : NULL),
                                      (BHAHAHA_AKV_ENABLE_GRIDPOINT ? bhahaha_akv_eval_J_integrand_full : NULL),
                                      &out);

#if BHAHAHA_AKV_ENABLE_QUALITY_EXTENSION
  /* Pack an 8-bit stub mask (excluding the always-present L1-integrands bit) so the
   * "any stub used" indicator is not permanently set in normal runs.
   */
  const uint32_t stub_mask8_all = (stub_flags & UINT32_C(0xFF));
  const uint32_t stub_mask8 = (stub_mask8_all & ~BHAHAHA_AKV_STUB_USED_L1_INTEGRANDS);
  if (stub_mask8) {
    uint32_t q = (uint32_t)out.akv_quality_flag;
    q |= (UINT32_C(1) << BHAHAHA_AKV_QUAL_EXT_ANY_STUB_BIT);
    q |= (stub_mask8 << BHAHAHA_AKV_QUAL_EXT_STUBMASK_SHIFT);
    out.akv_quality_flag = (int)q;
  }
#else
  (void)stub_flags;
#endif

#ifdef BHAHAHA_AKV_DIAGS_HAVE_FIELDS
  for (int a = 0; a < 3; a++) {
    bhahaha_diags->akv_lambda[a] = out.akv_lambda[a];
    bhahaha_diags->akv_J[a]      = out.akv_J[a];
    bhahaha_diags->akv_a[a]      = out.akv_a[a];
    bhahaha_diags->akv_spin_vec[a] = out.akv_spin_vec[a];
    bhahaha_diags->akv_eig_resid[a] = out.akv_eig_resid[a];
  }
  bhahaha_diags->akv_eig_gap_43 = out.akv_eig_gap_43;
  bhahaha_diags->akv_quality_flag = out.akv_quality_flag;
  bhahaha_diags->akv_method_used = (int)out.method_used;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      bhahaha_diags->akv_l1_evecs[i][j] = out.akv_l1_evecs[i][j];
#endif

  akv_diagnostics_free(&out);

#if !BHAHAHA_HAVE_TLS
  free(w2d);
#endif

#ifdef BHAHAHA_DIAGNOSTICS_SUCCESS
  return (err == AKV_SUCCESS) ? BHAHAHA_DIAGNOSTICS_SUCCESS : BHAHAHA_DIAGNOSTICS_FAILURE;
#else
  return (int)err;
#endif
"""

    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

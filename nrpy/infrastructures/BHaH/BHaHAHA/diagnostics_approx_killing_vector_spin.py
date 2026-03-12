"""
Scaffold registration and implementation module.

This code registers the AKV C scaffold implements the equations
through SymPy to generate efficient C code,
referencing the callbacks defined in the scaffold.

Author: Wesley Inselman
"""

from __future__ import annotations

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
from nrpy.equations.general_relativity.bhahaha.approx_killing_vector_spin import (
    ApproxKillingSpinClass,
)
from nrpy.infrastructures import BHaH


def register_CFunction_diagnostics_akv_spin(
    enable_fd_functions: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BHaHAHA approximate Killing vector (AKV) spin diagnostic C function.

    This generates C code for the AKV spin diagnostic, including code for the
    horizon surface element and the L1-reduced AKV integrands, then registers
    the resulting 'diagnostics_akv_spin' C function with NRPy. During the
    parallel-codegen registration phase, the function records the call and
    returns 'None' without generating code.

    :param enable_fd_functions: Whether to have NRPy emit finite-difference helper
        functions when generating C code for expressions that require finite
        differencing.

    :return: 'None' during the parallel-codegen registration phase. Otherwise,
        returns the populated NRPy environment after registering the generated
        C function.

    :raises FileNotFoundError: If the required companion C scaffold file
        'bah_diagnostics_akv_spin.c' is not found alongside this module.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    akv_c_path = Path(__file__).with_name("bah_diagnostics_akv_spin.c")
    if not akv_c_path.is_file():
        raise FileNotFoundError(f"Missing required companion file: {akv_c_path}")
    akv_c_code = akv_c_path.read_text(encoding="utf-8")

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

    area_codegen_cb = ccg.c_codegen(
        BHaH.BHaHAHA.area.area3(),
        "const REAL sqrt_det_q_cb",
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
"""
        + "\n"
        + akv_c_code
        + r"""
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
  REAL *restrict *xx = ctx->xx;
  const REAL *restrict in_gfs = ctx->in_gfs;

#include "set_CodeParameters.h"

  // p is the full 2D ghosted horizon-grid index: p = it_full*N_phi_tot + ip_full.
  const int N_phi_tot = grid->N_phi_tot;
  const int ip_full = p % N_phi_tot;            // phi index (full, incl. ghosts)
  const int it_full = (p - ip_full) / N_phi_tot;// theta index (full, incl. ghosts)
  const int i2 = ip_full;                       // AH phi index (full, incl. ghosts)
  const int i1 = it_full;                       // AH theta index (full, incl. ghosts)
  const int i0 = BHAHAHA_AKV_AH_SLAB_I0;         // AH slab index in the 3D grid

  // Compute the surface element using the exact same sqrt(det q) as used in w2d.
  """
        + area_codegen_cb
        + r"""

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
    name = "diagnostics_akv_spin"
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
    body += (
        area_codegen
        + r"""
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
    )

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

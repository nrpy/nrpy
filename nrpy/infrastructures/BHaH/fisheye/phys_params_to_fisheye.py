"""
Register a C helper that converts physical fisheye parameters to internal fisheye parameters.

This module registers a C function that solves for (R_i, s_i) and computes c given physical
transition centers and widths (r_trans_i, w_trans_i) and the outer radius L,
matching the conventions in nrpy/equations/generalrfm/fisheye.py for an
N-transition fisheye (zoom factors a0 -> a1 -> ... -> aN).

The physical parameters are registered as NRPy CodeParameters in commondata_struct.

Physical parameter meanings (physical radius):
- fisheye_phys_L: outer physical boundary radius.
- fisheye_phys_r_trans{i}: physical radius where transition i is centered.
- fisheye_phys_w_trans{i}: physical width of transition i.
"""

from typing import List, Union

import nrpy.params as par
from nrpy.c_function import register_CFunction
from nrpy.equations.generalrfm import fisheye as fisheye_eqs


def _register_physical_fisheye_codeparams(num_transitions: int) -> None:
    """
    Register physical fisheye CodeParameters in commondata_struct.

    :param num_transitions: Number of fisheye transitions.
    """
    _ = par.register_CodeParameter(
        "REAL",
        __name__,
        "fisheye_phys_L",
        10.0,
        commondata=True,
        add_to_parfile=True,
        description="Outer physical boundary radius.",
    )

    r_names = [f"fisheye_phys_r_trans{i + 1}" for i in range(num_transitions)]
    w_names = [f"fisheye_phys_w_trans{i + 1}" for i in range(num_transitions)]

    # Reasonable defaults: increasing centers, moderate widths.
    r_defaults: List[Union[str, int, float]] = [
        2.0 * (i + 1) for i in range(num_transitions)
    ]
    w_defaults: List[Union[str, int, float]] = [1.0 for _ in range(num_transitions)]

    _ = par.register_CodeParameters(
        "REAL",
        __name__,
        r_names,
        r_defaults,
        commondata=True,
        add_to_parfile=True,
        descriptions=[
            f"Physical center radius of transition {i + 1}."
            for i in range(num_transitions)
        ],
    )
    _ = par.register_CodeParameters(
        "REAL",
        __name__,
        w_names,
        w_defaults,
        commondata=True,
        add_to_parfile=True,
        descriptions=[
            f"Physical width of transition {i + 1}." for i in range(num_transitions)
        ],
    )


def _c_array_initializer(values: List[str]) -> str:
    return ", ".join(values)


def build_post_params_struct_set_to_default_hook(
    num_transitions: int,
    compute_griddata: str,
) -> str:
    """
    Build C code for main() hook after params_struct_set_to_default().

    This hook computes internal fisheye parameters from commondata->fisheye_phys_*.

    :param num_transitions: Number of fisheye transitions.
    :param compute_griddata: Griddata array symbol in generated main C code.
    :raises ValueError: If num_transitions is less than 1.
    :return: C code string suitable for post_params_struct_set_to_default hook.
    """
    if num_transitions < 1:
        raise ValueError("num_transitions must be >= 1")

    return f"""for(int grid=0; grid<commondata.NUMGRIDS; grid++) {{
  if (fisheye_params_from_physical_N{num_transitions}(&commondata, &{compute_griddata}[grid].params) != 0) {{
    fprintf(stderr, "Error: fisheye_params_from_physical_N{num_transitions} failed for grid %d. Check fisheye_a* and fisheye_phys_* values.\\n", grid);
    exit(1);
  }}
}}\n"""


def register_CFunction_fisheye_params_from_physical_N(
    num_transitions: int = 2,
) -> None:
    """
    Register a C function that solves for fisheye (R_i, s_i) and computes c from physical parameters.

    The solver uses a damped Newton method with a numerical Jacobian to enforce:
      R_phys(R_i) = r_trans_i,
      R_phys(R_i+s_i) - R_phys(R_i-s_i) = w_trans_i,
    with R_phys(r) = c * rbar_unscaled(r), c = L / rbar_unscaled(L).

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> for n in (2, 3):
    ...     cfc.CFunction_dict.clear()
    ...     register_CFunction_fisheye_params_from_physical_N(n)
    ...     name = f"fisheye_params_from_physical_N{n}"
    ...     generated_str = cfc.CFunction_dict[name].full_function
    ...     validate_strings(generated_str, name, file_ext="c")

    :param num_transitions: Number of fisheye transitions.
    :raises ValueError: If num_transitions is less than 1.
    """
    if num_transitions < 1:
        raise ValueError("num_transitions must be >= 1")

    # Ensure internal fisheye CodeParameters (a_i, R_i, s_i, c) are registered.
    _ = fisheye_eqs.build_fisheye(num_transitions=num_transitions)

    # Register physical parameters in commondata_struct.
    _register_physical_fisheye_codeparams(num_transitions)

    includes = ["<math.h>", "<stdio.h>", "<stdlib.h>"]
    desc = """\
Compute fisheye internal parameters (R_i, s_i) and compute c from physical fisheye inputs
(r_trans_i, w_trans_i, L) for an N-transition fisheye.

Physical parameter meanings (physical radius):
- fisheye_phys_L: outer physical boundary radius.
- fisheye_phys_r_trans{i}: physical radius where transition i is centered.
- fisheye_phys_w_trans{i}: physical width of transition i.
"""
    cfunc_type = "int"
    name = f"fisheye_params_from_physical_N{num_transitions}"
    params = (
        "const commondata_struct *restrict commondata, params_struct *restrict params"
    )

    a_vals = [f"params->fisheye_a{i}" for i in range(num_transitions + 1)]
    r_trans_vals = [
        f"commondata->fisheye_phys_r_trans{i + 1}" for i in range(num_transitions)
    ]
    w_trans_vals = [
        f"commondata->fisheye_phys_w_trans{i + 1}" for i in range(num_transitions)
    ]

    prefunc = rf"""
#ifndef STANDALONE
#include "BHaH_defines.h"
#include "BHaH_function_prototypes.h"
#else
typedef double REAL;
typedef struct __commondata_struct__ {{
  REAL fisheye_phys_L;
  {" ".join([f"REAL fisheye_phys_r_trans{i + 1};" for i in range(num_transitions)])}
  {" ".join([f"REAL fisheye_phys_w_trans{i + 1};" for i in range(num_transitions)])}
}} commondata_struct;
typedef struct __params_struct__ {{
  {" ".join([f"REAL fisheye_a{i};" for i in range(num_transitions + 1)])}
  {" ".join([f"REAL fisheye_R{i + 1};" for i in range(num_transitions)])}
  {" ".join([f"REAL fisheye_s{i + 1};" for i in range(num_transitions)])}
  REAL fisheye_c;
}} params_struct;
#endif

// Numerically stable log(cosh(x)) helper for large |x|.
// log(cosh(x)) = |x| + log(1 + exp(-2|x|)) - log(2)
static inline REAL logcosh_stable(const REAL x) {{
  const REAL ax = fabs(x);
  return ax + log1p(exp((REAL)-2.0 * ax)) - log((REAL)2.0);
}}

static inline REAL log_cosh_ratio(const REAL u, const REAL v) {{
  return logcosh_stable(u) - logcosh_stable(v);
}}

static inline REAL rbar_unscaled(const REAL r,
                                 const REAL a[],
                                 const REAL R[],
                                 const REAL s[],
                                 const int N) {{
  REAL rb = a[N] * r;
  for (int i = 0; i < N; i++) {{
    const REAL delta_a = a[i] - a[i + 1];
    const REAL denom = (REAL)2.0 * tanh(R[i] / s[i]);
    const REAL u = (r + R[i]) / s[i];
    const REAL v = (r - R[i]) / s[i];
    const REAL term = (delta_a * s[i]) / denom * log_cosh_ratio(u, v);
    rb += term;
  }}
  return rb;
}}

// Return 0 on success, nonzero on failure (invalid/non-finite scaling).
static inline int evaluate_constraints(const REAL L,
                                        const REAL r_trans[],
                                        const REAL w_trans[],
                                        const REAL a[],
                                        const REAL R[],
                                        const REAL s[],
                                        const int N,
                                        REAL F[], REAL *c_out) {{
  const REAL rbar_L = rbar_unscaled(L, a, R, s, N);
  if (!(isfinite(rbar_L)) || !(rbar_L > (REAL)0.0)) return 1;
  const REAL c = L / rbar_L;
  if (!(isfinite(c)) || !(c > (REAL)0.0)) return 1;

  for (int i = 0; i < N; i++) {{
    const REAL Rm = R[i] - s[i];
    const REAL Rp = R[i] + s[i];
    const REAL Rphys_R = c * rbar_unscaled(R[i], a, R, s, N);
    const REAL Rphys_Rp = c * rbar_unscaled(Rp, a, R, s, N);
    const REAL Rphys_Rm = c * rbar_unscaled(Rm, a, R, s, N);
    F[2 * i + 0] = Rphys_R - r_trans[i];
    F[2 * i + 1] = (Rphys_Rp - Rphys_Rm) - w_trans[i];
  }}

  if (c_out) *c_out = c;
  return 0;
}}

static inline int solve_linear_system(const int n,
                                      const REAL *A_in,
                                      const REAL *b_in,
                                      REAL *x_out) {{
  // Gaussian elimination with partial pivoting; A_in is row-major n x n.
  enum {{ MAX_LINEAR_N = 64 }};
  if (n <= 0 || n > MAX_LINEAR_N) {{
    fprintf(stderr, "ERROR: fisheye solve_linear_system requires 0 < n <= %d (got n=%d)\\n",
            MAX_LINEAR_N, n);
    return 1;
  }}
  REAL A[MAX_LINEAR_N][MAX_LINEAR_N + 1];
  for (int i = 0; i < n; i++) {{
    for (int j = 0; j < n; j++) A[i][j] = A_in[i * n + j];
    A[i][n] = b_in[i];
  }}

  for (int k = 0; k < n; k++) {{
    int pivot = k;
    REAL maxval = fabs(A[k][k]);
    for (int i = k + 1; i < n; i++) {{
      const REAL val = fabs(A[i][k]);
      if (val > maxval) {{ maxval = val; pivot = i; }}
    }}
    if (!(maxval > (REAL)0.0)) return 1;
    if (pivot != k) {{
      for (int j = k; j < n + 1; j++) {{
        const REAL tmp = A[k][j];
        A[k][j] = A[pivot][j];
        A[pivot][j] = tmp;
      }}
    }}
    const REAL inv_pivot = (REAL)1.0 / A[k][k];
    for (int j = k; j < n + 1; j++) A[k][j] *= inv_pivot;
    for (int i = 0; i < n; i++) {{
      if (i == k) continue;
      const REAL factor = A[i][k];
      for (int j = k; j < n + 1; j++) A[i][j] -= factor * A[k][j];
    }}
  }}
  for (int i = 0; i < n; i++) x_out[i] = A[i][n];
  return 0;
}}
"""

    body = rf"""
  const int NTRANS = {num_transitions};
  const int NUNK = 2 * NTRANS;

  const REAL L = commondata->fisheye_phys_L;
  const REAL a[{num_transitions + 1}] = {{ {_c_array_initializer(a_vals)} }};
  const REAL r_trans[{num_transitions}] = {{ {_c_array_initializer(r_trans_vals)} }};
  const REAL w_trans[{num_transitions}] = {{ {_c_array_initializer(w_trans_vals)} }};

  for (int i = 0; i < NTRANS + 1; i++) {{
    if (!(a[i] > (REAL)0.0)) return 1;
  }}
  if (!(L > (REAL)0.0)) return 1;
  for (int i = 0; i < NTRANS; i++) {{
    if (!(r_trans[i] > (REAL)0.0) || !(w_trans[i] > (REAL)0.0)) return 1;
    if (i > 0 && !(r_trans[i] > r_trans[i - 1])) return 1;
  }}

  REAL x[NUNK];
  for (int i = 0; i < NTRANS; i++) {{
    x[2 * i + 0] = r_trans[i];
    x[2 * i + 1] = (REAL)0.5 * w_trans[i];
  }}

  const REAL tol = (REAL)1e-12;
  const int max_iter = 80;
  const REAL eps = (REAL)1e-6;

  int converged = 0;
  REAL c = (REAL)1.0;

  for (int iter = 0; iter < max_iter; iter++) {{
    // Build R and s arrays from x
    REAL R[NTRANS];
    REAL s[NTRANS];
    for (int i = 0; i < NTRANS; i++) {{
      R[i] = x[2 * i + 0];
      s[i] = x[2 * i + 1];
      if (!(R[i] > (REAL)0.0) || !(s[i] > (REAL)0.0)) return 1;
      if (!(R[i] > s[i])) return 1;
    }}

    REAL F[NUNK];
    if (evaluate_constraints(L, r_trans, w_trans, a, R, s, NTRANS, F, &c)) return 1;

    REAL Fnorm = (REAL)0.0;
    for (int i = 0; i < NUNK; i++) Fnorm += fabs(F[i]);
    if (Fnorm < tol) {{ converged = 1; break; }}

    // Numerical Jacobian
    REAL J[NUNK][NUNK];
    for (int j = 0; j < NUNK; j++) {{
      const REAL xj = x[j];
      const REAL dx = eps * (fabs(xj) + (REAL)1.0);
      x[j] = xj + dx;

      REAL Rp[NTRANS];
      REAL sp[NTRANS];
      for (int i = 0; i < NTRANS; i++) {{
        Rp[i] = x[2 * i + 0];
        sp[i] = x[2 * i + 1];
      }}

      REAL Fp[NUNK];
      if (evaluate_constraints(L, r_trans, w_trans, a, Rp, sp, NTRANS, Fp, NULL)) return 1;
      for (int i = 0; i < NUNK; i++) {{
        J[i][j] = (Fp[i] - F[i]) / dx;
      }}
      x[j] = xj;
    }}

    REAL delta[NUNK];
    REAL minusF[NUNK];
    for (int i = 0; i < NUNK; i++) minusF[i] = -F[i];

    if (solve_linear_system(NUNK, &J[0][0], minusF, delta)) return 1;

    // Damped step to enforce positivity and improve robustness
    REAL alpha = (REAL)1.0;
    int accepted = 0;
    for (int attempt = 0; attempt < 12; attempt++) {{
      int ok = 1;
      for (int i = 0; i < NUNK; i++) {{
        const REAL trial = x[i] + alpha * delta[i];
        if (!(trial > (REAL)0.0)) {{ ok = 0; break; }}
      }}
      if (ok) {{
        for (int i = 0; i < NTRANS; i++) {{
          const REAL Rtrial = x[2 * i + 0] + alpha * delta[2 * i + 0];
          const REAL strial = x[2 * i + 1] + alpha * delta[2 * i + 1];
          if (!(Rtrial > strial)) {{ ok = 0; break; }}
        }}
      }}
      if (ok) {{
        REAL x_trial[NUNK];
        for (int i = 0; i < NUNK; i++) x_trial[i] = x[i] + alpha * delta[i];
        REAL Rt[NTRANS];
        REAL st[NTRANS];
        for (int i = 0; i < NTRANS; i++) {{
          Rt[i] = x_trial[2 * i + 0];
          st[i] = x_trial[2 * i + 1];
        }}
        REAL F_trial[NUNK];
        if (evaluate_constraints(L, r_trans, w_trans, a, Rt, st, NTRANS, F_trial, NULL)) return 1;
        REAL Fnorm_trial = (REAL)0.0;
        for (int i = 0; i < NUNK; i++) Fnorm_trial += fabs(F_trial[i]);
        if (Fnorm_trial < Fnorm) {{
          for (int i = 0; i < NUNK; i++) x[i] = x_trial[i];
          accepted = 1;
          break;
        }}
      }}
      alpha *= (REAL)0.5;
    }}
    if (!accepted) return 1;
  }}

  if (!converged) return 1;

  // Commit results into params
  for (int i = 0; i < NTRANS; i++) {{
    const REAL Ri = x[2 * i + 0];
    const REAL si = x[2 * i + 1];
    switch (i) {{
"""

    # Emit switch cases for R_i and s_i assignments (1-based names).
    for i in range(num_transitions):
        body += rf"""
    case {i}:
      params->fisheye_R{i + 1} = Ri;
      params->fisheye_s{i + 1} = si;
      break;
"""

    body += r"""
    }
  }
  params->fisheye_c = c;
  return 0;
"""

    # Standalone defaults: monotone transition radii and widths; a_i = 2^i.
    standalone_r_trans_vals: List[float] = []
    standalone_w_trans_vals: List[float] = []
    if num_transitions >= 1:
        standalone_r_trans_vals.append(4.5)
        standalone_w_trans_vals.append(1.0)
    if num_transitions >= 2:
        standalone_r_trans_vals.append(8.0)
        standalone_w_trans_vals.append(1.5)
    for i in range(3, num_transitions + 1):
        standalone_r_trans_vals.append(8.0 + 4.0 * (i - 2))
        standalone_w_trans_vals.append(1.5)
    standalone_L = standalone_r_trans_vals[-1] + 2.0
    standalone_a_vals = [float(2**i) for i in range(num_transitions + 1)]

    def _fmt_list(vals: List[float]) -> str:
        return ", ".join(f"{v:.15e}" for v in vals)

    standalone_r_trans = _fmt_list(standalone_r_trans_vals)
    standalone_w_trans = _fmt_list(standalone_w_trans_vals)
    standalone_a = _fmt_list(standalone_a_vals)

    postfunc = rf"""
#ifdef STANDALONE

#define FISHEYE_GRID_N 70

static int write_fisheye_grid_txt(const char *fname) {{
  const int N = FISHEYE_GRID_N;
  const REAL L = {standalone_L:.15e};
  // Physical parameters:
  //   r_trans[i]: physical transition centers; w_trans[i]: physical widths.
  const REAL r_trans[] = {{ {standalone_r_trans} }};
  const REAL w_trans[] = {{ {standalone_w_trans} }};
  const REAL a[] = {{ {standalone_a} }};
  const int NTRANS = (int)(sizeof(r_trans) / sizeof(r_trans[0]));
  const int NUNK = 2 * NTRANS;

  for (int i = 0; i < NTRANS + 1; i++) {{
    if (!(a[i] > (REAL)0.0))
      return 1;
  }}

  REAL x[NUNK];
  for (int i = 0; i < NTRANS; i++) {{
    x[2 * i + 0] = r_trans[i];
    x[2 * i + 1] = (REAL)0.5 * w_trans[i];
  }}
  REAL c = (REAL)1.0;
  const REAL tol = (REAL)1e-12;
  const int max_iter = 60;

  int converged = 0;
  for (int iter = 0; iter < max_iter; iter++) {{
  for (int i = 0; i < NUNK; i++) {{
      if (!(x[i] > (REAL)0.0))
        return 1;
    }}

    REAL R[NTRANS];
    REAL s[NTRANS];
    for (int i = 0; i < NTRANS; i++) {{
      R[i] = x[2 * i + 0];
      s[i] = x[2 * i + 1];
      if (!(R[i] > s[i]))
        return 1;
    }}
    REAL F[NUNK];
    if (evaluate_constraints(L, r_trans, w_trans, a, R, s, NTRANS, F, &c)) return 1;

    REAL Fnorm = (REAL)0.0;
    for (int i = 0; i < NUNK; i++) Fnorm += fabs(F[i]);
    if (Fnorm < tol) {{
      converged = 1;
      break;
    }}

    REAL J[NUNK][NUNK];
    const REAL eps = (REAL)1e-6;
    for (int j = 0; j < NUNK; j++) {{
      REAL xj = x[j];
      const REAL dx = eps * (fabs(xj) + (REAL)1.0);
      x[j] = xj + dx;

      REAL Rp[NTRANS];
      REAL sp[NTRANS];
      for (int i = 0; i < NTRANS; i++) {{
        Rp[i] = x[2 * i + 0];
        sp[i] = x[2 * i + 1];
      }}
      REAL Fp[NUNK];
      if (evaluate_constraints(L, r_trans, w_trans, a, Rp, sp, NTRANS, Fp, NULL)) return 1;
      for (int i = 0; i < NUNK; i++) {{
        J[i][j] = (Fp[i] - F[i]) / dx;
      }}
      x[j] = xj;
    }}

    REAL delta[NUNK];
    REAL minusF[NUNK];
    for (int i = 0; i < NUNK; i++) minusF[i] = -F[i];
    if (solve_linear_system(NUNK, &J[0][0], minusF, delta))
      return 1;

    REAL alpha = (REAL)1.0;
    int accepted = 0;
    for (int attempt = 0; attempt < 12; attempt++) {{
      int ok = 1;
      REAL x_trial[NUNK];
      for (int i = 0; i < NUNK; i++) {{
        x_trial[i] = x[i] + alpha * delta[i];
        if (!(x_trial[i] > (REAL)0.0)) ok = 0;
      }}
      if (ok) {{
        for (int i = 0; i < NTRANS; i++) {{
          const REAL Rtrial = x[2 * i + 0] + alpha * delta[2 * i + 0];
          const REAL strial = x[2 * i + 1] + alpha * delta[2 * i + 1];
          if (!(Rtrial > strial)) ok = 0;
        }}
      }}
      if (ok) {{
        REAL Rt[NTRANS];
        REAL st[NTRANS];
        for (int i = 0; i < NTRANS; i++) {{
          Rt[i] = x_trial[2 * i + 0];
          st[i] = x_trial[2 * i + 1];
        }}
        REAL F_trial[NUNK];
        if (evaluate_constraints(L, r_trans, w_trans, a, Rt, st, NTRANS, F_trial, NULL)) return 1;
        REAL Fnorm_trial = (REAL)0.0;
        for (int i = 0; i < NUNK; i++) Fnorm_trial += fabs(F_trial[i]);
        if (Fnorm_trial < Fnorm) {{
          for (int i = 0; i < NUNK; i++) x[i] = x_trial[i];
          accepted = 1;
          break;
        }}
      }}
      alpha *= (REAL)0.5;
    }}
    if (!accepted)
      return 1;
  }}
  if (!converged)
    return 1;

  REAL R[NTRANS];
  REAL s[NTRANS];
  for (int i = 0; i < NTRANS; i++) {{
    R[i] = x[2 * i + 0];
    s[i] = x[2 * i + 1];
    if (!(R[i] > s[i]))
      return 1;
  }}
  const REAL dx = (REAL)2.0 * L / (REAL)(N - 1);

  FILE *fp = fopen(fname, "w");
  if (!fp)
    return 1;
  for (int i = 0; i < N; i++) {{
    const REAL X = -L + dx * (REAL)i;
    for (int j = 0; j < N; j++) {{
      const REAL Y = -L + dx * (REAL)j;
      const REAL r = sqrt(X * X + Y * Y);
      REAL Xb = 0.0, Yb = 0.0;
      if (r > (REAL)0.0) {{
        const REAL Rphys = c * rbar_unscaled(r, a, R, s, NTRANS);
        const REAL scale = Rphys / r;
        Xb = X * scale;
        Yb = Y * scale;
      }}
      fprintf(fp, "%d %d %.15e %.15e\n", i, j, Xb, Yb);
    }}
  }}
  fclose(fp);
  return 0;
}}

int main(void) {{
  if (write_fisheye_grid_txt("fisheye_grid.txt") != 0) {{
    fprintf(stderr, "ERROR: write_fisheye_grid_txt failed.\n");
    return EXIT_FAILURE;
  }}
  printf("Wrote fisheye_grid.txt\n");
  printf("Plot command:\n");
  printf("python3 -c \"import numpy as np, matplotlib.pyplot as plt; d=np.loadtxt('fisheye_grid.txt'); "
         "N=int(d[:,0].max()+1); Xb=d[:,2].reshape(N,N); Yb=d[:,3].reshape(N,N); "
         "[plt.plot(Xb[i,:],Yb[i,:],c='k',lw=0.45) for i in range(N)]; "
         "[plt.plot(Xb[:,j],Yb[:,j],c='k',lw=0.45) for j in range(N)]; "
         "r_trans=[{standalone_r_trans}]; "
         "theta=np.linspace(0,2*np.pi,400); "
         "[plt.plot(r*np.cos(theta), r*np.sin(theta), c='k', lw=2.5, ls=(0,(2,4))) for r in r_trans]; "
         "plt.gca().set_aspect('equal','box'); plt.axis('off'); plt.show()\"\n");
  return EXIT_SUCCESS;
}}

#endif // STANDALONE
"""

    register_CFunction(
        subdirectory="fisheye",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        postfunc=postfunc,
    )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()
    if results.failed:
        raise SystemExit(1)

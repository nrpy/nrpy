r"""
Reference:
J. A. Faber, T. W. Baumgarte, Z. B. Etienne, S. L. Shapiro, and K. Taniguchi,
"Relativistic hydrodynamics in the presence of puncture black holes,"
Phys. Rev. D 96, 104035 (2017), Appendix: Fisheye Coordinates.

Fisheye radial profile r'(r) and its derivative dr'/dr (equations-style module).

This module provides reusable expressions for the purely *radial* fisheye map:

  - FisheyeRadialProfile.rprime_of_r(r):  r'(r)
  - FisheyeRadialProfile.drprime_dr(r):   d(r')/dr

Multi-transition profile:

  r' = a_n r + sum_{i=1..n} ((a_{i-1} - a_i) s_i / [2 tanh(R_i/s_i)])
       * ln[ cosh((r+R_i)/s_i) / cosh((r-R_i)/s_i) ],

  dr'/dr = a_n + sum_{i=1..n} ((a_{i-1} - a_i) / tanh(R_i/s_i))
           * [ (tanh((r+R_i)/s_i) - tanh((r-R_i)/s_i)) / 2 ].

Author: Nishita Jadoo
"""

from __future__ import annotations
from typing import Optional, Sequence

import sympy as sp  # SymPy: CAS

# Optional read-only dependency; not required for doctests.
try:
    import nrpy.params as par  # NRPy+: parameter interface
except Exception:  # pragma: no cover
    par = None


class FisheyeRadialProfile:
    """
    Build symbolic expressions for the *radial* fisheye profile r'(r) and dr'/dr.

    Parameters
    ----------
    n : int in [0..4]
        Number of fisheye transition zones.
    a, R, s : sequences (length ≥ 5)
        Fisheye coefficients (a), transition centers (R), and widths (s).
        Only entries up to n (or n+1 for a) are used.

    Notes
    -----
    If any parameter is omitted and `nrpy.params` is available, values are taken
    from CodeParameters: "fisheye_n", "fisheye_a", "fisheye_R", "fisheye_s".

    When n=0: r' = a0 * r and dr'/dr = a0.

    Doctests
    --------
    >>> # Analytic sanity check: n=0, a0=2 -> r' = 2 r; dr'/dr = 2
    >>> prof = FisheyeRadialProfile(n=0, a=[2,1,1,1,1], R=[0,0,0,0,0], s=[1,1,1,1,1])
    >>> prof.rprime_of_r(sp.Integer(3))
    6
    >>> sp.simplify(prof.drprime_dr(sp.symbols("r", positive=True)))
    2

    >>> # Symbolic identity (strong): d/dr r'(r) == dr'/dr for n=0
    >>> r = sp.symbols("r", positive=True)
    >>> rp = prof.rprime_of_r(r)
    >>> drp = prof.drprime_dr(r)
    >>> sp.simplify(sp.diff(rp, r) - drp)
    0

    >>> # Symbolic identity for a simple n=1 case
    >>> prof1 = FisheyeRadialProfile(n=1, a=[2,1,1,1,1], R=[3,0,0,0,0], s=[2,1,1,1,1])
    >>> r = sp.symbols("r", positive=True)
    >>> sp.simplify(sp.diff(prof1.rprime_of_r(r), r) - prof1.drprime_dr(r))
    0
    """

    __slots__ = ("n", "a", "R", "s")

    def __init__(
        self,
        n: Optional[int] = None,
        a: Optional[Sequence[object]] = None,
        R: Optional[Sequence[object]] = None,
        s: Optional[Sequence[object]] = None,
    ) -> None:
        # Optionally read from nrpy.params, but never register here.
        if n is None or a is None or R is None or s is None:
            if par is not None:
                if n is None:
                    n = int(par.parval_from_str("fisheye_n"))
                if a is None:
                    a = list(par.parval_from_str("fisheye_a"))
                if R is None:
                    R = list(par.parval_from_str("fisheye_R"))
                if s is None:
                    s = list(par.parval_from_str("fisheye_s"))
            else:
                raise ValueError(
                    "FisheyeRadialProfile requires (n, a, R, s), or provide CodeParameters via nrpy.params."
                )

        if not (0 <= int(n) <= 4):
            raise ValueError("fisheye_n must be in [0, 4].")
        if len(a) < 5 or len(R) < 5 or len(s) < 5:
            raise ValueError("a, R, s must each have length ≥ 5.")

        self.n = int(n)
        self.a = tuple(sp.sympify(val) for val in a[:5])
        self.R = tuple(sp.sympify(val) for val in R[:5])
        self.s = tuple(sp.sympify(val) for val in s[:5])

    # -------------------------
    # Core expressions
    # -------------------------
    def rprime_of_r(self, r: sp.Expr) -> sp.Expr:
        """Return r'(r) for the current radial fisheye settings."""
        rp = self.a[self.n] * r
        for i in range(1, self.n + 1):
            a_im1, a_i = self.a[i - 1], self.a[i]
            R_i, s_i = self.R[i - 1], self.s[i - 1]
            denom = 2 * sp.tanh(R_i / s_i)
            term = sp.log(sp.cosh((r + R_i) / s_i) / sp.cosh((r - R_i) / s_i))
            rp += (a_im1 - a_i) * s_i / denom * term
        return sp.simplify(rp)

    def drprime_dr(self, r: sp.Expr) -> sp.Expr:
        """Return dr'/dr for the current radial fisheye settings."""
        drp = self.a[self.n]
        for i in range(1, self.n + 1):
            a_im1, a_i = self.a[i - 1], self.a[i]
            R_i, s_i = self.R[i - 1], self.s[i - 1]
            coeff = (a_im1 - a_i) / sp.tanh(R_i / s_i)
            bracket = (sp.tanh((r + R_i) / s_i) - sp.tanh((r - R_i) / s_i)) / 2
            drp += coeff * bracket
        return sp.simplify(drp)


# ================================================================
# Doctest + validate_expressions hook (numeric trusted results only)
# ================================================================
if __name__ == "__main__":
    import doctest
    import os
    import sys
    from mpmath import mp, mpf

    # Optional: generate/compare trusted results (if NRPy test harness is available)
    try:
        import nrpy.validate_expressions.validate_expressions as ve
    except Exception:
        ve = None

    # 1) run doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # helper to convert SymPy expr to high-precision mpf
    def to_mpf(expr, prec_dps: int = 80) -> mpf:
        mp.dps = prec_dps
        return mpf(str(sp.N(expr, prec_dps)))

    # ---------- Trusted test: analytic baseline n=0 ----------
    prof0 = FisheyeRadialProfile(n=0, a=[2, 1, 1, 1, 1], R=[0, 0, 0, 0, 0], s=[1, 1, 1, 1, 1])
    r0_a = sp.Integer(3)
    r0_b = sp.Integer(7)
    results_dict = {
        "n0_rprime_at_3": to_mpf(prof0.rprime_of_r(r0_a)),   # = 6
        "n0_drprime_at_7": to_mpf(prof0.drprime_dr(r0_b)),   # = 2
    }

    # ---------- Trusted test: nontrivial two-transition case n=2 ----------
    # Parameters chosen to be smooth and positive, with distinct transitions.
    prof2 = FisheyeRadialProfile(
        n=2,
        a=[2.0, 1.5, 1.0, 1.0, 1.0],
        R=[3.0, 6.0, 0.0, 0.0, 0.0],
        s=[1.0, 1.5, 1.0, 1.0, 1.0],
    )
    r2_a = sp.Integer(4)
    r2_b = sp.Integer(10)
    results_dict.update({
        "n2_rprime_at_4": to_mpf(prof2.rprime_of_r(r2_a)),
        "n2_drprime_at_4": to_mpf(prof2.drprime_dr(r2_a)),
        "n2_rprime_at_10": to_mpf(prof2.rprime_of_r(r2_b)),
        "n2_drprime_at_10": to_mpf(prof2.drprime_dr(r2_b)),
    })

    if ve is not None:
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            os.path.splitext(os.path.basename(__file__))[0],  # => tests/trusted/fisheye_radial_profile.py
            results_dict,
        )

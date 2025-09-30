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
        njadoo **at** uidaho **dot* edu
"""

from __future__ import annotations
from typing import Optional, Sequence, Tuple

import sympy as sp  # SymPy: CAS
import nrpy.params as par  # NRPy+: parameter interface


def _make_symbol_list(prefix: str, indices: Sequence[int]) -> Tuple[sp.Symbol, ...]:
    """Create a tuple of symbols like ('a0','a1',...) or ('R1','R2',...)."""
    return tuple(sp.symbols(" ".join(f"{prefix}{i}" for i in indices)))


class FisheyeRadialProfile:
    r"""
    Build symbolic expressions for the *radial* fisheye profile r'(r) and dr'/dr.

    Parameters
    ----------
    n : int, optional
        Number of fisheye transition zones. If omitted, read from
        `par.parval_from_str("fisheye_n")`.
    a, R, s : sequences of scalars/symbols, optional
        Fisheye coefficients (a), transition centers (R), and widths (s).
        If omitted, **symbolic** sequences of appropriate lengths are created:

          - a = (a0, a1, ..., a_n)           length n+1
          - R = (R1, R2, ..., R_n)           length n
          - s = (s1, s2, ..., s_n)           length n

        You may pass numeric or symbolic entries. Only entries up to n are used.

    Notes
    -----
    When n = 0: r' = a0 * r and dr'/dr = a0.

    Doctests (self-contained: do not use nrpy.params)
    -----------------------------------------------
    >>> # Analytic sanity check: n=0, a0=2 -> r' = 2 r; dr'/dr = 2
    >>> prof = FisheyeRadialProfile(n=0, a=[2])
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

    >>> # Symbolic identity for a simple n=1 case with explicit numbers
    >>> prof1 = FisheyeRadialProfile(n=1, a=[2, 1], R=[3], s=[2])
    >>> r = sp.symbols("r", positive=True)
    >>> sp.simplify(sp.diff(prof1.rprime_of_r(r), r) - prof1.drprime_dr(r))
    0

    >>> # Purely symbolic build via automatic symbols when args omitted:
    >>> # (we pass n explicitly so doctest remains self-contained)
    >>> prof_sym = FisheyeRadialProfile(n=2)
    >>> len(prof_sym.a), len(prof_sym.R), len(prof_sym.s)
    (3, 2, 2)
    >>> r = sp.symbols("r", positive=True)
    >>> sp.simplify(sp.diff(prof_sym.rprime_of_r(r), r) - prof_sym.drprime_dr(r))
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
        # Determine n
        if n is None:
            n = int(par.parval_from_str("fisheye_n"))
        n = int(n)
        if n < 0:
            raise ValueError("fisheye_n must be >= 0.")
        self.n = n

        # Build a, R, s (use user-provided if present; else generate *symbols* of correct lengths)
        if a is None:
            # a0..a_n
            self.a = _make_symbol_list("a", range(0, n + 1))
        else:
            if len(a) < n + 1:
                raise ValueError(f"length of a must be at least n+1 = {n+1}.")
            self.a = tuple(sp.sympify(val) for val in a[: n + 1])

        if n == 0:
            # No transitions; R and s unused.
            self.R = tuple()
            self.s = tuple()
        else:
            if R is None:
                # R1..R_n
                self.R = _make_symbol_list("R", range(1, n + 1))
            else:
                if len(R) < n:
                    raise ValueError(f"length of R must be at least n = {n}.")
                self.R = tuple(sp.sympify(val) for val in R[:n])

            if s is None:
                # s1..s_n
                self.s = _make_symbol_list("s", range(1, n + 1))
            else:
                if len(s) < n:
                    raise ValueError(f"length of s must be at least n = {n}.")
                self.s = tuple(sp.sympify(val) for val in s[:n])

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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

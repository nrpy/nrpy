"""
General N transition fisheye raw to physical (xx -> Cart) mapping and reference metric.

This module implements a radially symmetric, arbitrary N transition fisheye mapping
from raw Cartesian coordinates xx^i to physical Cartesian coordinates Cart^i,
together with the associated flat space reference metric in the raw coordinates and
its first and second derivatives:

* Raw to physical map: xx^i -> Cart^i(xx).
* Reference metric: ghat_ij = delta_mn (d Cart^m / d xx^i) (d Cart^n / d xx^j).
* First derivatives: ghat_ij,k computed using sympy.diff.
* Second derivatives: ghat_ij,kl computed using sympy.diff.

The raw radius is

    r = sqrt( sum_i (xx^i)^2 ).

The N transition fisheye map is defined by

* Plateau stretch (zoom out) factors a_0, ..., a_N.
* Transition centers R_1, ..., R_N.
* Width parameters s_1, ..., s_N.
* Differences Delta a_i = a_{i-1} - a_i.

The single transition kernel is

    G(r; R, s) =
        s / (2 * tanh(R / s)) *
        log( cosh( (r + R) / s ) / cosh( (r - R) / s ) ).

The unscaled radius map is

    rbar_unscaled(r) =
        a_N * r + sum_{i=1}^N (a_{i-1} - a_i) * G(r; R_i, s_i).

A global scale factor c produces the final radius map

    rbar(r) = c * rbar_unscaled(r).

The physical Cartesian coordinates are obtained by a purely radial rescaling,

    Cart^i = (rbar(r) / r) * xx^i,

which leaves the angular coordinates unchanged.

From this mapping, the flat reference metric in raw coordinates is

    ghat_ij = delta_mn (d Cart^m / d xx^i) (d Cart^n / d xx^j),

and its derivatives are constructed via direct differentiation with respect to
the xx^i using sympy.diff only. No subs or simplify calls are used anywhere in
this module. The inverse mapping Cart -> xx is intentionally not provided here;
it is handled by a separate Newton Raphson based module.

Fisheye parameters are stored as NRPy code parameters so that they can be
configured at code generation time.
"""

from typing import List, Optional, Sequence, Union, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par

CodeParameterDefaultList = List[Union[str, int, float]]


class FisheyeGeneralRFM:
    """
    Construct an N transition fisheye raw to physical map and reference metric.

    The mapping is radially symmetric in the raw Cartesian coordinates xx^i and
    is defined by a multi-transition kernel. Given plateau stretch factors a_i,
    transition centers R_i, widths s_i, and a global scale c, the raw radius
    r = sqrt( sum_i (xx^i)^2 ) is mapped to a physical radius rbar(r). The
    physical Cartesian coordinates are then

        Cart^i = (rbar(r) / r) * xx^i,

    for r > 0, with the expression understood in the limiting sense at r = 0.

    From this transformation the flat space reference metric ghat_ij in the raw
    coordinates and its first and second derivatives are constructed by direct
    differentiation with respect to xx^i.

    :ivar num_transitions: Number of fisheye transitions N.
    :ivar a_list: Plateau stretch factors a_0, ..., a_N, stored as NRPy code
        parameters.
    :ivar R_list: Raw transition centers R_1, ..., R_N, stored as NRPy code
        parameters.
    :ivar s_list: Raw transition widths s_1, ..., s_N, stored as NRPy code
        parameters.
    :ivar c: Global scaling factor c, stored as an NRPy code parameter.
    :ivar xx: Raw Cartesian coordinate symbols [xx0, xx1, xx2].
    :ivar r: Raw radius r = sqrt( sum_i xx_i^2 ).
    :ivar rbar_unscaled: Unscaled physical radius rbar_unscaled(r).
    :ivar rbar: Scaled physical radius rbar(r) = c * rbar_unscaled(r).
    :ivar xx_to_CartU: Raw to physical Cartesian map Cart^i(xx^j).
    :ivar dCart_dxxUD: Jacobian d Cart^i / d xx^j.
    :ivar ghatDD: Reference metric ghat_ij in raw coordinates.
    :ivar ghatDDdD: First derivatives ∂_k ghat_ij.
    :ivar ghatDDdDD: Second derivatives ∂_l ∂_k ghat_ij.
    """

    def __init__(
        self,
        num_transitions: int,
        a_default: Optional[Sequence[float]] = None,
        R_default: Optional[Sequence[float]] = None,
        s_default: Optional[Sequence[float]] = None,
        c_default: float = 1.0,
    ) -> None:
        """
        Initialize the N transition fisheye mapping and reference metric in 3D.

        Fisheye parameters are registered as NRPy code parameters so that they
        can be overridden at code generation time.

        :param num_transitions: Number of fisheye transitions N. Must be at least 1.
        :param a_default: Default plateau stretch factors [a0, ..., aN].
                          If None, all entries are set to 1.0. Length must be
                          num_transitions + 1.
        :param R_default: Default raw transition centers [R1, ..., RN].
                          If None, they are set to [1.0, 2.0, ..., float(N)].
                          Length must be num_transitions.
        :param s_default: Default raw transition widths [s1, ..., sN].
                          If None, all entries are set to 0.5. Length must be
                          num_transitions.
        :param c_default: Default global scaling factor c. Defaults to 1.0.

        :raises ValueError: If num_transitions is less than 1.
        :raises ValueError: If any of the default lists have inconsistent lengths.
        """
        if num_transitions < 1:
            raise ValueError(
                f"num_transitions must be >= 1; got num_transitions = {num_transitions}."
            )

        self.num_transitions = num_transitions

        # ---------------------------------------------------------------------
        # Step 1: Set up default parameter values and register NRPy CodeParameters.
        # ---------------------------------------------------------------------
        if a_default is None:
            # Default: all plateaus have unit stretch before global scaling.
            a_default_list: List[float] = [1.0 for _ in range(num_transitions + 1)]
        else:
            a_default_list = [float(val) for val in a_default]
        if len(a_default_list) != num_transitions + 1:
            raise ValueError(
                "a_default must have length num_transitions + 1; "
                f"got len(a_default) = {len(a_default_list)} while "
                f"num_transitions + 1 = {num_transitions + 1}."
            )

        if R_default is None:
            # Simple monotonically increasing centers as a fallback.
            R_default_list: List[float] = [float(i + 1) for i in range(num_transitions)]
        else:
            R_default_list = [float(val) for val in R_default]
        if len(R_default_list) != num_transitions:
            raise ValueError(
                "R_default must have length num_transitions; "
                f"got len(R_default) = {len(R_default_list)} while "
                f"num_transitions = {num_transitions}."
            )

        if s_default is None:
            # Default: moderate transition widths.
            s_default_list: List[float] = [0.5 for _ in range(num_transitions)]
        else:
            s_default_list = [float(val) for val in s_default]
        if len(s_default_list) != num_transitions:
            raise ValueError(
                "s_default must have length num_transitions; "
                f"got len(s_default) = {len(s_default_list)} while "
                f"num_transitions = {num_transitions}."
            )

        # Plateau stretch factors a_0..a_N:
        a_names = [f"fisheye_a{i}" for i in range(num_transitions + 1)]
        self.a_list = list(
            par.register_CodeParameters(
                "REAL",
                __name__,
                a_names,
                cast(CodeParameterDefaultList, a_default_list),
                commondata=True,
            )
        )

        # Raw transition centers R_1..R_N:
        R_names = [f"fisheye_R{i + 1}" for i in range(num_transitions)]
        self.R_list = list(
            par.register_CodeParameters(
                "REAL",
                __name__,
                R_names,
                cast(CodeParameterDefaultList, R_default_list),
                commondata=True,
            )
        )

        # Raw transition widths s_1..s_N:
        s_names = [f"fisheye_s{i + 1}" for i in range(num_transitions)]
        self.s_list = list(
            par.register_CodeParameters(
                "REAL",
                __name__,
                s_names,
                cast(CodeParameterDefaultList, s_default_list),
                commondata=True,
            )
        )

        # Global scale factor c:
        self.c = par.register_CodeParameter(
            "REAL",
            __name__,
            "fisheye_c",
            c_default,
            commondata=True,
        )

        # ---------------------------------------------------------------------
        # Step 2: Define raw coordinates and radial map.
        # ---------------------------------------------------------------------
        # Raw Cartesian coordinates xx^i (3D).
        self.xx = list(ixp.declarerank1("xx", dimension=3))

        # Raw radius r = sqrt(xx^2 + yy^2 + zz^2):
        self.r = sp.sqrt(sum(self.xx[i] ** 2 for i in range(3)))

        # Build the unscaled radius map rbar_unscaled(r) =
        # a_N r + sum_{i=1}^N Delta a_i G(r; R_i, s_i).
        self.rbar_unscaled = self._build_unscaled_radius_map(self.r)

        # Final scaled physical radius: rbar(r) = c * rbar_unscaled(r).
        self.rbar = self.c * self.rbar_unscaled

        # ---------------------------------------------------------------------
        # Step 3: Raw (xx) -> physical Cartesian (Cart) map.
        #         Cart^i = (rbar / r) * xx^i
        # ---------------------------------------------------------------------
        self.xx_to_CartU = list(ixp.zerorank1(dimension=3))
        for i in range(3):
            self.xx_to_CartU[i] = (self.rbar / self.r) * self.xx[i]

        # ---------------------------------------------------------------------
        # Step 4: Jacobian d Cart^i / d xx^j.
        # ---------------------------------------------------------------------
        self.dCart_dxxUD = [[sp.sympify(0) for _ in range(3)] for __ in range(3)]
        for mu in range(3):
            for i in range(3):
                # IMPORTANT: use sp.diff for all derivatives.
                self.dCart_dxxUD[mu][i] = sp.diff(self.xx_to_CartU[mu], self.xx[i])

        # ---------------------------------------------------------------------
        # Step 5: Reference metric ghat_ij = delta_mn (d Cart^m / d xx^i)
        #         (d Cart^n / d xx^j).
        # ---------------------------------------------------------------------
        self.ghatDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                g_ij = sp.sympify(0)
                for mu in range(3):
                    g_ij += self.dCart_dxxUD[mu][i] * self.dCart_dxxUD[mu][j]
                self.ghatDD[i][j] = g_ij

        # ---------------------------------------------------------------------
        # Step 6: First derivatives ghat_ij,k via sp.diff.
        # ---------------------------------------------------------------------
        self.ghatDDdD = ixp.zerorank3(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ghatDDdD[i][j][k] = sp.diff(self.ghatDD[i][j], self.xx[k])

        # ---------------------------------------------------------------------
        # Step 7: Second derivatives ghat_ij,kl via sp.diff.
        # ---------------------------------------------------------------------
        self.ghatDDdDD = ixp.zerorank4(dimension=3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        self.ghatDDdDD[i][j][k][l] = sp.diff(
                            self.ghatDDdD[i][j][k], self.xx[l]
                        )

    # -------------------------------------------------------------------------
    # Internal helper methods.
    # -------------------------------------------------------------------------

    def _G_kernel(self, r: sp.Expr, R: sp.Expr, s: sp.Expr) -> sp.Expr:
        """
        Single transition kernel G(r; R, s) used in the multi transition map.

        The kernel is

            G(r; R, s) =
                s / (2 * tanh(R / s)) *
                log( cosh( (r + R) / s ) / cosh( (r - R) / s ) ).

        :param r: Raw radius r.
        :param R: Raw transition center R.
        :param s: Raw transition width s, assumed to be strictly positive.
        :return: The kernel value G(r; R, s).
        """
        prefactor = s / (2 * sp.tanh(R / s))
        arg_plus = (r + R) / s
        arg_minus = (r - R) / s
        return prefactor * sp.log(sp.cosh(arg_plus) / sp.cosh(arg_minus))

    def _build_unscaled_radius_map(self, r: sp.Expr) -> sp.Expr:
        """
        Construct the unscaled N transition radius map rbar_unscaled(r).

        The map is

            rbar_unscaled(r) =
                a_N * r + sum_{i=1}^N (a_{i-1} - a_i) * G(r; R_i, s_i),

        where G(r; R, s) is the single transition kernel defined in _G_kernel.

        :param r: Raw radius r.
        :return: The unscaled physical radius rbar_unscaled(r).
        """
        a_N = self.a_list[-1]
        rbar_unscaled = a_N * r

        # Delta a_i = a_{i-1} - a_i for i = 1..N
        for i in range(self.num_transitions):
            delta_a_i = self.a_list[i] - self.a_list[i + 1]
            R_i = self.R_list[i]
            s_i = self.s_list[i]
            rbar_unscaled += delta_a_i * self._G_kernel(r, R_i, s_i)

        return rbar_unscaled


def build_fisheye_generalrfm(
    n_fisheye_transitions: int,
    a_default: Optional[Sequence[float]] = None,
    R_default: Optional[Sequence[float]] = None,
    s_default: Optional[Sequence[float]] = None,
    c_default: float = 1.0,
) -> FisheyeGeneralRFM:
    """
    Construct a FisheyeGeneralRFM instance.

    This function instantiates FisheyeGeneralRFM using the input
    n_fisheye_transitions value and passes through the default plateau
    stretches, transition centers, widths, and global scale, which in turn
    define the default values of the underlying NRPy CodeParameters.

    :param n_fisheye_transitions: Number of fisheye transitions N. Must be at least 1.
    :param a_default: Default plateau stretch factors [a0, ..., aN].
                      If None, all entries are set to 1.0. Length must be
                      n_fisheye_transitions + 1.
    :param R_default: Default raw transition centers [R1, ..., RN].
                      If None, they are set to [1.0, 2.0, ..., float(N)].
                      Length must be n_fisheye_transitions.
    :param s_default: Default raw transition widths [s1, ..., sN].
                      If None, all entries are set to 0.5. Length must be
                      n_fisheye_transitions.
    :param c_default: Default global scaling factor c. Defaults to 1.0.
    :return: A newly constructed FisheyeGeneralRFM instance.

    :raises ValueError: If n_fisheye_transitions is less than 1.
    :raises ValueError: If any of the default lists have inconsistent lengths.
    """
    if n_fisheye_transitions < 1:
        raise ValueError(
            "n_fisheye_transitions must be >= 1; "
            f"got n_fisheye_transitions = {n_fisheye_transitions}."
        )

    return FisheyeGeneralRFM(
        num_transitions=n_fisheye_transitions,
        a_default=a_default,
        R_default=R_default,
        s_default=s_default,
        c_default=c_default,
    )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Run doctests for this module (none are defined explicitly, but this keeps
    # behavior consistent with other NRPy modules).
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Basic validation for a small set of configurations (3D only).
    for N in (1, 2):
        print(f"Setting up FisheyeGeneralRFM[N={N}]...")
        fisheye = FisheyeGeneralRFM(num_transitions=N)
        results_dict = ve.process_dictionary_of_expressions(
            fisheye.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_N{N}",
            results_dict,
        )

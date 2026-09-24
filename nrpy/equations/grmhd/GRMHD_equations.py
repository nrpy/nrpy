"""
Construct GRMHD expressions from the GRHD fluid and metric quantities.

The comoving magnetic four-vector and magnetic stress-energy follow Duez et
al., Phys. Rev. D 72, 024028 (2005), https://arxiv.org/abs/astro-ph/0503420v2.
Reference-metric component rescaling follows Jacques et al.,
https://arxiv.org/abs/2412.03659v2. This module does not evolve the magnetic
field.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

# Step 1.a: import all needed modules from NRPy:
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from nrpy.equations.general_relativity.g4munu_conversions import (
    ADM_to_g4DD,
    ADM_to_g4UU,
)
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def compute_smallb4U(
    gammaDD: List[List[sp.Expr]],
    betaU: List[sp.Expr],
    alpha: sp.Expr,
    u4U: List[sp.Expr],
    BmagU: List[sp.Expr],
) -> List[sp.Expr]:
    """
    Compute the comoving magnetic four-vector from the Eulerian field BmagU.

    Duez et al. (2005), Eqs. (23)-(24) and (31). BmagU includes the
    Gaussian-to-Heaviside-Lorentz factor 1/sqrt(4*pi) before this call.
    BmagU, gammaDD, betaU, and u4U must use the same spatial basis.

    :param gammaDD: Covariant spatial metric.
    :param betaU: Shift vector.
    :param alpha: Lapse.
    :param u4U: Fluid four-velocity.
    :param BmagU: Eulerian magnetic field divided by sqrt(4 pi), in the basis of gammaDD.
    :return: Comoving magnetic four-vector.
    """
    g4DD = ADM_to_g4DD(gammaDD, betaU, alpha)
    u_dot_B = sp.sympify(0)
    for i in range(3):
        for mu in range(4):
            u_dot_B += g4DD[i + 1][mu] * u4U[mu] * BmagU[i]

    smallb4U = ixp.zerorank1(dimension=4)
    smallb4U[0] = u_dot_B / alpha
    for i in range(3):
        smallb4U[i + 1] = (BmagU[i] + u_dot_B * u4U[i + 1]) / (alpha * u4U[0])
    return smallb4U


def compute_smallb2(
    gammaDD: List[List[sp.Expr]],
    betaU: List[sp.Expr],
    alpha: sp.Expr,
    smallb4U: List[sp.Expr],
) -> sp.Expr:
    """
    Contract the comoving magnetic four-vector with the spacetime metric.

    The scalar b^2 = b^mu b_mu appears in Duez et al. (2005), Eq. (32).

    :param gammaDD: Covariant spatial metric.
    :param betaU: Shift vector.
    :param alpha: Lapse.
    :param smallb4U: Comoving magnetic four-vector.
    :return: Magnetic scalar b^mu b_mu.
    """
    g4DD = ADM_to_g4DD(gammaDD, betaU, alpha)
    smallb2 = sp.sympify(0)
    for mu in range(4):
        for nu in range(4):
            smallb2 += g4DD[mu][nu] * smallb4U[mu] * smallb4U[nu]
    return smallb2


class GRMHD_Equations(GRHD_Equations):
    """Add magnetic stress-energy to the inherited GRHD equations."""

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
    ) -> None:
        """
        Initialize and set up all GRMHD quantities, storing them within the class object.

        BmagU^i = ReU[i] * rescaledBmagU[i] uses the reference-metric
        scale factors of Jacques et al., Eq. (22).

        :param enable_rfm_precompute: Whether to enable reference-metric
            precomputation, defaults to False.
        :param CoordSystem: The coordinate system being used, defaults
            to "Cartesian".

        """
        super().__init__(CoordSystem, enable_rfm_precompute)
        self.rescaledBmagU = ixp.declarerank1("rescaledBmagU", dimension=3)
        self.BmagU = ixp.zerorank1(dimension=3)
        for i in range(3):
            self.BmagU[i] = self.rescaledBmagU[i] * self.ReU[i]
        self.smallb4U: List[sp.Expr]
        self.smallb2: sp.Expr
        self._fluid_T4UU: Tuple[Tuple[sp.Expr, ...], ...]
        self._fluid_rescaledT4UU: Tuple[Tuple[sp.Expr, ...], ...]

    def compute_T4UU(self) -> None:
        """
        Add magnetic stress-energy from Duez et al. (2005), Eqs. (32)-(33).

        Tensor-component rescaling uses Jacques et al., Eq. (22).
        """
        super().compute_T4UU()
        # Keep the fluid tensors for compute_T4UD and compute_tau_tilde_fluxU,
        #   which apply the GRHD methods to the fluid part only. Tuples keep
        #   them out of the trusted-value dictionaries.
        self._fluid_T4UU = tuple(tuple(row) for row in self.T4UU)
        self._fluid_rescaledT4UU = tuple(tuple(row) for row in self.rescaledT4UU)
        gammaDD = self.gammaDD
        betaU = self.betaU
        alpha = self.alpha
        u4U = self.u4U
        BmagU = self.BmagU

        self.smallb4U = compute_smallb4U(gammaDD, betaU, alpha, u4U, BmagU)
        self.smallb2 = compute_smallb2(gammaDD, betaU, alpha, self.smallb4U)
        smallb4U = self.smallb4U
        smallb2 = self.smallb2

        # define g^{mu nu} in terms of the ADM quantities:
        g4UU = ADM_to_g4UU(gammaDD, betaU, alpha)

        # add the magnetic part of T^{mu nu}
        for mu in range(4):
            for nu in range(4):
                magnetic = (
                    smallb2 * u4U[mu] * u4U[nu]
                    + sp.Rational(1, 2) * smallb2 * g4UU[mu][nu]
                    - smallb4U[mu] * smallb4U[nu]
                )
                self.T4UU[mu][nu] += magnetic
                rescale_mu = self.ReU[mu - 1] if mu else sp.sympify(1)
                rescale_nu = self.ReU[nu - 1] if nu else sp.sympify(1)
                self.rescaledT4UU[mu][nu] += magnetic / (rescale_mu * rescale_nu)

    def compute_T4UD(self) -> None:
        """
        Lower one index of the fluid and magnetic stress-energy tensors.

        The fluid part uses the GRHD method. The magnetic part is lowered
        analytically, b^2 u^mu u_nu + (b^2/2) delta^mu_nu - b^mu b_nu, which
        avoids combining the magnetic terms over a common denominator.
        Mixed-index component rescaling follows Jacques et al., Eqs. (22)-(24).
        """
        T4UU = self.T4UU
        rescaledT4UU = self.rescaledT4UU
        self.T4UU = [list(row) for row in self._fluid_T4UU]
        self.rescaledT4UU = [list(row) for row in self._fluid_rescaledT4UU]
        super().compute_T4UD()
        self.T4UU = T4UU
        self.rescaledT4UU = rescaledT4UU

        gammaDD = self.gammaDD
        betaU = self.betaU
        alpha = self.alpha
        u4U = self.u4U
        smallb4U = self.smallb4U
        smallb2 = self.smallb2

        # we'll need g_{alpha nu} in terms of ADM quantities:
        g4DD = ADM_to_g4DD(gammaDD, betaU, alpha)
        u4D = ixp.zerorank1(dimension=4)
        smallb4D = ixp.zerorank1(dimension=4)
        for mu in range(4):
            for nu in range(4):
                u4D[mu] += g4DD[mu][nu] * u4U[nu]
                smallb4D[mu] += g4DD[mu][nu] * smallb4U[nu]

        # add the magnetic part of T^mu_nu
        for mu in range(4):
            for nu in range(4):
                magnetic = smallb2 * u4U[mu] * u4D[nu] - smallb4U[mu] * smallb4D[nu]
                if mu == nu:
                    magnetic += sp.Rational(1, 2) * smallb2
                self.T4UD[mu][nu] += magnetic
                rescale_mu = self.ReU[mu - 1] if mu else sp.sympify(1)
                rescale_nu = self.ReU[nu - 1] if nu else sp.sympify(1)
                self.rescaledT4UD[mu][nu] += magnetic * rescale_nu / rescale_mu

    def compute_tau_tilde_fluxU(self) -> None:
        """
        Compute the energy flux with the fluid part from the GRHD method.

        The magnetic part alpha^2 e^(6 phi) T_EM^{0j} is added without
        combining it over a common denominator.
        """
        T4UU = self.T4UU
        rescaledT4UU = self.rescaledT4UU
        fluid_T4UU = [list(row) for row in self._fluid_T4UU]
        fluid_rescaledT4UU = [list(row) for row in self._fluid_rescaledT4UU]
        self.T4UU = fluid_T4UU
        self.rescaledT4UU = fluid_rescaledT4UU
        super().compute_tau_tilde_fluxU()
        self.T4UU = T4UU
        self.rescaledT4UU = rescaledT4UU

        alpha = self.alpha
        e6phi = self.e6phi
        for j in range(3):
            self.tau_tilde_fluxU[j] += (
                alpha**2 * e6phi * (T4UU[0][j + 1] - fluid_T4UU[0][j + 1])
            )
            self.rescaled_tau_tilde_fluxU[j] += (
                alpha**2
                * e6phi
                * (rescaledT4UU[0][j + 1] - fluid_rescaledT4UU[0][j + 1])
            )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Step 1: With B^i = 0, every GRHD expression must be reproduced exactly.
    grhd_eqs = GRHD_Equations(CoordSystem="Cartesian")
    grhd_eqs.construct_all_equations()
    grmhd_zero_eqs = GRMHD_Equations(CoordSystem="Cartesian")
    grmhd_zero_eqs.BmagU = ixp.zerorank1(dimension=3)
    grmhd_zero_eqs.construct_all_equations()
    for key, value in grhd_eqs.__dict__.items():
        if isinstance(value, (sp.Basic, list)) and (
            grmhd_zero_eqs.__dict__[key] != value
        ):
            raise AssertionError(f"{key} with B^i = 0 differs from GRHD")

    # Step 2: Compare the magnetic expressions with trusted values.
    for Coord in [
        "Spherical",
        "SinhSpherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        enable_rfm_pre = "rfm_precompute" in Coord
        grmhd_eqs = GRMHD_Equations(
            Coord.replace("_rfm_precompute", ""),
            enable_rfm_precompute=enable_rfm_pre,
        )
        grmhd_eqs.construct_all_equations()
        results_dict = ve.process_dictionary_of_expressions(
            grmhd_eqs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )

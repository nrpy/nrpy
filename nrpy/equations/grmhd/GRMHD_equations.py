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

from typing import List

import sympy as sp

import nrpy.indexedexp as ixp
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


class GRMHDEquations(GRHD_Equations):
    """Add magnetic stress-energy to the inherited GRHD equations."""

    def __init__(
        self, CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
    ) -> None:
        """
        Initialize fluid, metric, and rescaled Eulerian magnetic symbols.

        BmagU^i = ReU[i] * rescaledBmagU[i] uses the reference-metric
        scale factors of Jacques et al., Eq. (22).

        :param CoordSystem: Reference-metric coordinate system.
        :param enable_rfm_precompute: Whether to precompute reference-metric factors.
        """
        super().__init__(CoordSystem, enable_rfm_precompute)
        self.rescaledBmagU = ixp.declarerank1("rescaledBmagU", dimension=3)
        self.BmagU = ixp.zerorank1(dimension=3)
        for i in range(3):
            self.BmagU[i] = self.rescaledBmagU[i] * self.ReU[i]
        self.smallb4U: List[sp.Expr]
        self.smallb2: sp.Expr

    def compute_T4UU(self) -> None:
        """
        Add magnetic stress-energy from Duez et al. (2005), Eqs. (32)-(33).

        Tensor-component rescaling uses Jacques et al., Eq. (22).
        """
        super().compute_T4UU()
        self.smallb4U = compute_smallb4U(
            self.gammaDD, self.betaU, self.alpha, self.u4U, self.BmagU
        )
        self.smallb2 = compute_smallb2(
            self.gammaDD, self.betaU, self.alpha, self.smallb4U
        )
        g4UU = ADM_to_g4UU(self.gammaDD, self.betaU, self.alpha)
        for mu in range(4):
            for nu in range(4):
                magnetic = (
                    self.smallb2 * self.u4U[mu] * self.u4U[nu]
                    + sp.Rational(1, 2) * self.smallb2 * g4UU[mu][nu]
                    - self.smallb4U[mu] * self.smallb4U[nu]
                )
                self.T4UU[mu][nu] += magnetic
                rescale_mu = self.ReU[mu - 1] if mu else sp.sympify(1)
                rescale_nu = self.ReU[nu - 1] if nu else sp.sympify(1)
                self.rescaledT4UU[mu][nu] += magnetic / (rescale_mu * rescale_nu)


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

    # Step 1: Compare every shared expression with GRHD when B^i = 0.
    grhd_eqs = GRHD_Equations(CoordSystem="Cartesian")
    grmhd_zero = GRMHDEquations(CoordSystem="Cartesian")
    grmhd_zero.BmagU = ixp.zerorank1(dimension=3)
    grhd_eqs.construct_all_equations()
    grmhd_zero.construct_all_equations()
    shared_results = (
        "T4UU",
        "rescaledT4UU",
        "T4UD",
        "rescaledT4UD",
        "rho_star",
        "Ye_star",
        "S_star",
        "tau_tilde",
        "S_tildeD",
        "rescaledS_tildeD",
        "rho_star_fluxU",
        "rescaled_rho_star_fluxU",
        "Ye_star_fluxU",
        "rescaled_Ye_star_fluxU",
        "S_star_fluxU",
        "rescaled_S_star_fluxU",
        "tau_tilde_fluxU",
        "rescaled_tau_tilde_fluxU",
        "S_tilde_fluxUD",
        "rescaled_S_tilde_fluxUD",
        "tau_source_term",
        "rho_star_connection_term",
        "Ye_star_connection_term",
        "S_star_connection_term",
        "tau_connection_term",
        "S_tilde_source_termD",
        "S_tilde_connection_termsD",
    )
    for result_name in shared_results:
        if getattr(grmhd_zero, result_name) != getattr(grhd_eqs, result_name):
            raise AssertionError(f"{result_name}: GRMHD with B^i = 0 differs from GRHD")
    if grmhd_zero.smallb2 != 0:
        raise AssertionError("B=0 must give b^2=0")

    # Step 2: Compare Cartesian magnetic expressions with trusted values.
    grmhd_eqs = GRMHDEquations(CoordSystem="Cartesian")
    grmhd_eqs.gammaDD = sp.diag(1, 4, 9).tolist()
    grmhd_eqs.betaU = [sp.Rational(1, 10), sp.Rational(1, 20), -sp.Rational(1, 30)]
    grmhd_eqs.alpha = sp.Rational(5, 4)
    grmhd_eqs.e6phi = sp.sympify(6)
    grmhd_eqs.u4U = [
        sp.sympify(1),
        sp.Rational(13, 20),
        -sp.Rational(1, 20),
        sp.Rational(1, 30),
    ]
    grmhd_eqs.rho_b = sp.sympify(1)
    grmhd_eqs.h = sp.Rational(3, 2)
    grmhd_eqs.P = sp.Rational(1, 3)
    grmhd_eqs.BmagU = ixp.declarerank1("BmagU", dimension=3)
    grmhd_eqs.compute_vU_from_u4U__no_speed_limit()
    grmhd_eqs.VU = grmhd_eqs.VU_from_u4U
    grmhd_eqs.compute_T4UU()
    grmhd_eqs.compute_T4UD()
    grmhd_eqs.compute_rho_star()
    grmhd_eqs.compute_tau_tilde()
    grmhd_eqs.compute_S_tildeD()
    grmhd_eqs.compute_tau_tilde_fluxU()
    grmhd_eqs.compute_S_tilde_fluxUD()
    # The sampled u^mu is normalized, so b^mu u_mu = 0 and the comoving
    # energy density is T^{mu nu} u_mu u_nu = rho_b h - P + b^2/2.
    sample_g4DD = ADM_to_g4DD(grmhd_eqs.gammaDD, grmhd_eqs.betaU, grmhd_eqs.alpha)
    u4D = [
        sum(sample_g4DD[mu][nu] * grmhd_eqs.u4U[nu] for nu in range(4))
        for mu in range(4)
    ]
    b_dot_u = sum(grmhd_eqs.smallb4U[mu] * u4D[mu] for mu in range(4))
    comoving_energy_density = sum(
        grmhd_eqs.T4UU[mu][nu] * u4D[mu] * u4D[nu] for mu in range(4) for nu in range(4)
    )
    if not ve.check_zero(b_dot_u, fixed_mpfs_for_free_symbols=True):
        raise AssertionError("Sampled b^mu u_mu is nonzero")
    if not ve.check_zero(
        comoving_energy_density
        - (grmhd_eqs.rho_b * grmhd_eqs.h - grmhd_eqs.P + grmhd_eqs.smallb2 / 2),
        fixed_mpfs_for_free_symbols=True,
    ):
        raise AssertionError("Sampled T^{mu nu} u_mu u_nu != rho_b h - P + b^2/2")
    expressions = {
        "smallb4U": grmhd_eqs.smallb4U,
        "smallb2": grmhd_eqs.smallb2,
        "T4UU": grmhd_eqs.T4UU,
        "T4UD": grmhd_eqs.T4UD,
        "tau_tilde": grmhd_eqs.tau_tilde,
        "S_tildeD": grmhd_eqs.S_tildeD,
        "tau_tilde_fluxU": grmhd_eqs.tau_tilde_fluxU,
        "S_tilde_fluxUD": grmhd_eqs.S_tilde_fluxUD,
    }
    sampled_results = ve.process_dictionary_of_expressions(
        expressions, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        "GRMHD_equations_Cartesian",
        sampled_results,
    )

    # Step 3: Compare Spherical rescaled fields, stress-energy, and energy flux
    #         with trusted values. The radially moving fluid gives nonzero
    #         magnetic T^{0i}.
    spherical_eqs = GRMHDEquations(CoordSystem="Spherical")
    radius, theta = ixp.declarerank1("xx", dimension=3)[:2]
    spherical_eqs.gammaDD = sp.diag(
        1, radius**2, radius**2 * sp.sin(theta) ** 2
    ).tolist()
    spherical_eqs.betaU = ixp.zerorank1(dimension=3)
    spherical_eqs.alpha = sp.sympify(1)
    spherical_eqs.e6phi = sp.sympify(1)
    spherical_eqs.u4U = [sp.Rational(5, 4), sp.Rational(3, 4), 0, 0]
    spherical_eqs.rho_b = sp.sympify(1)
    spherical_eqs.h = sp.Rational(3, 2)
    spherical_eqs.P = sp.Rational(1, 3)
    spherical_eqs.compute_vU_from_u4U__no_speed_limit()
    spherical_eqs.VU = spherical_eqs.VU_from_u4U
    spherical_eqs.compute_T4UU()
    spherical_eqs.compute_rho_star()
    spherical_eqs.compute_tau_tilde_fluxU()
    spherical_expressions = {
        "BmagU": spherical_eqs.BmagU,
        "smallb4U": spherical_eqs.smallb4U,
        "smallb2": spherical_eqs.smallb2,
        "T4UU": spherical_eqs.T4UU,
        "rescaledT4UU": spherical_eqs.rescaledT4UU,
        "rescaled_tau_tilde_fluxU": spherical_eqs.rescaled_tau_tilde_fluxU,
    }
    spherical_sampled_results = ve.process_dictionary_of_expressions(
        spherical_expressions, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        "GRMHD_equations_Spherical",
        spherical_sampled_results,
    )

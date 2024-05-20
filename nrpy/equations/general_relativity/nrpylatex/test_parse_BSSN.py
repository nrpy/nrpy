"""
Test that handwritten BSSN quantities within equations/general_relativity/*.py match those generated by NRPyLaTeX below.

Author: Ken Sible
Email:  ksible *at* outlook *dot* com
"""

from typing import List, cast

import sympy as sp
from nrpylatex import parse_latex  # type: ignore

import nrpy.params as par
from nrpy.equations.general_relativity import (
    BSSN_constraints,
    BSSN_gauge_RHSs,
    BSSN_quantities,
    BSSN_RHSs,
)
from nrpy.validate_expressions.validate_expressions import assert_equal


def test_example_BSSN() -> bool:
    """
    Function to test the BSSN example.

    :returns: True if test is successful, otherwise False.
    """
    # Initialize to junk values, to make the linter and mypy happy.
    h_rhsDD: List[List[sp.Expr]] = [[]]
    a_rhsDD: List[List[sp.Expr]] = [[]]
    RbarDD: List[List[sp.Expr]] = [[]]
    bet_rhsU: List[sp.Expr] = []
    vet_rhsU: List[sp.Expr] = []
    Lambdabar_rhsU: List[sp.Expr] = []
    MU: List[sp.Expr] = []
    H = cast(sp.Expr, None)
    cf_rhs = cast(sp.Expr, None)
    alpha_rhs = cast(sp.Expr, None)
    trK_rhs = cast(sp.Expr, None)

    parse_latex(
        r"""
        % declare coord x y z
        % ignore "\begin{align}" "\end{align}" "\\%" "\qquad"

        \begin{align}
            % declare metric gammahatDD --zeros --dim 3
            % \hat{\gamma}_{ii} = 1 % noimpsum
            % declare metric gammabarDD hDD --dim 3
            % \bar{\gamma}_{ij} = h_{ij} + \hat{\gamma}_{ij}

            % replace "\beta" -> "\mathrm{vet}"
            % declare vetU --dim 3
            %% upwind pattern inside Lie derivative expansion
            % replace "\mathrm{vet}^{\1*} \partial_{\1*}" -> "\mathrm{vet}^{\1*} % suffix dupD
            \partial_{\1*}"
            %% substitute tensor identity (see appropriate BSSN notebook)
            % replace "\bar{D}_k \mathrm{vet}^k" -> "(\partial_k \mathrm{vet}^k + \frac{\partial_k \mathrm{gammahatdet} \mathrm{vet}^k}{2 \mathrm{gammahatdet}})"

            % replace "\bar{A}" -> "\mathrm{a}"
            % declare aDD --dim 3 --sym sym01
            % replace "\partial_t \bar{\gamma}" -> "\mathrm{h_rhs}"
            \partial_t \bar{\gamma}_{ij} &= \mathcal{L}_\beta \bar{\gamma}_{ij} + \frac{2}{3} \bar{\gamma}_{ij} \left(\alpha \bar{A}^k{}_k - \bar{D}_k \beta^k\right) - 2 \alpha \bar{A}_{ij} \\

            % declare cf trK phi --dim 3
            % replace "K" -> "\mathrm{trK}"
            %% replace 'phi' with conformal factor cf = W = e^{-2\phi}
            % replace "e^{-4\phi}" -> "\mathrm{cf}^2"
            % replace "\partial_t \phi = \1* \\" -> "\mathrm{cf_rhs} = -2 \mathrm{cf} (\1*) \\"
            % replace "\partial_{\1*} \phi" -> "\partial_{\1*} \mathrm{cf} \frac{-1}{2 \mathrm{cf}}"
            \partial_t \phi &= \mathcal{L}_\beta \phi + \frac{1}{6} \left(\bar{D}_k \beta^k - \alpha K \right) \\

            % declare alpha --dim 3
            % replace "\partial_t \mathrm{trK}" -> "\mathrm{trK_rhs}"
            \partial_t K &= \mathcal{L}_\beta K + \frac{1}{3} \alpha K^2 + \alpha \bar{A}_{ij} \bar{A}^{ij}
                - e^{-4\phi} \left(\bar{D}_i \bar{D}^i \alpha + 2 \bar{D}^i \alpha \bar{D}_i \phi\right) \\

            % replace "\bar{\Lambda}" -> "\mathrm{lambda}"
            % declare lambdaU --dim 3
            % \Delta^k_{ij} = \bar{\Gamma}^k_{ij} - \hat{\Gamma}^k_{ij}
            %% assign DeltaUDD --metric gammabar
            % \Delta^k = \bar{\gamma}^{ij} \Delta^k_{ij}
            % replace "\partial_t \mathrm{lambda}" -> "\mathrm{Lambdabar_rhs}"
            \partial_t \bar{\Lambda}^i &= \mathcal{L}_\beta \bar{\Lambda}^i + \bar{\gamma}^{jk} \hat{D}_j \hat{D}_k \beta^i
                + \frac{2}{3} \Delta^i \bar{D}_k \beta^k + \frac{1}{3} \bar{D}^i \bar{D}_k \beta^k \\%
                &\qquad- 2 \bar{A}^{ij} \left(\partial_j \alpha - 6 \alpha \partial_j \phi\right)
                + 2 \alpha \bar{A}^{jk} \Delta^i_{jk} - \frac{4}{3} \alpha \bar{\gamma}^{ij} \partial_j K \\

            % declare RbarDD --dim 3 --sym sym01
            X_{ij} &= -2 \alpha \bar{D}_i \bar{D}_j \phi + 4 \alpha \bar{D}_i \phi \bar{D}_j \phi
                + 2 \bar{D}_i \alpha \bar{D}_j \phi + 2 \bar{D}_j \alpha \bar{D}_i \phi
                - \bar{D}_i \bar{D}_j \alpha + \alpha \bar{R}_{ij} \\
            \hat{X}_{ij} &= X_{ij} - \frac{1}{3} \bar{\gamma}_{ij} \bar{\gamma}^{kl} X_{kl} \\
            % replace "\partial_t \mathrm{a}" -> "\mathrm{a_rhs}"
            \partial_t \bar{A}_{ij} &= \mathcal{L}_\beta \bar{A}_{ij} - \frac{2}{3} \bar{A}_{ij} \bar{D}_k \beta^k
                - 2 \alpha \bar{A}_{ik} \bar{A}^k_j + \alpha \bar{A}_{ij} K + e^{-4\phi} \hat{X}_{ij} \\

            % replace "\partial_t \alpha" -> "\mathrm{alpha_rhs}"
            \partial_t \alpha &= \mathcal{L}_\beta \alpha - 2 \alpha K \\

            % replace "B" -> "\mathrm{bet}"
            % declare betU --dim 3
            % replace "\partial_t \mathrm{vet}" -> "\mathrm{vet_rhs}"
            \partial_t \beta^i &= \left[\beta^j % suffix dupD
            \bar{D}_j \beta^i\right] + B^i \\

            % declare eta --const
            % replace "\partial_t \mathrm{bet}" -> "\mathrm{bet_rhs}"
            \partial_t B^i &= \left[\beta^j % suffix dupD
            \bar{D}_j B^i\right]
                + \frac{3}{4} \left(\partial_t \bar{\Lambda}^i - \left[\beta^j % suffix dupD
                \bar{D}_j \bar{\Lambda}^i\right]\right) - \eta B^i \\

            % \bar{R} = \bar{\gamma}^{ij} \bar{R}_{ij}
            % replace "\bar{D}^2" -> "\bar{D}^i \bar{D}_i"
            % replace "\mathcal{\1*}" -> "\mathrm{\1*}"
            \mathcal{H} &= \frac{2}{3} K^2 - \bar{A}_{ij} \bar{A}^{ij}
                + e^{-4\phi} \left(\bar{R} - 8 \bar{D}^i \phi \bar{D}_i \phi - 8 \bar{D}^2 \phi\right) \\

            \mathcal{M}^i &= e^{-4\phi} \left(\bar{D}_j \bar{A}^{ij} + 6 \bar{A}^{ij} \partial_j \phi
                - \frac{2}{3} \bar{\gamma}^{ij} \partial_j K\right) \\

            \bar{R}_{ij} &= -\frac{1}{2} \bar{\gamma}^{kl} \hat{D}_k \hat{D}_l \bar{\gamma}_{ij}
                + \frac{1}{2} \left(\bar{\gamma}_{ki} \hat{D}_j \bar{\Lambda}^k + \bar{\gamma}_{kj} \hat{D}_i \bar{\Lambda}^k\right)
                + \frac{1}{2} \Delta^k \left(\Delta_{ijk} + \Delta_{jik}\right) \\%
                &\qquad+ \bar{\gamma}^{kl} \left(\Delta^m_{ki} \Delta_{jml} + \Delta^m_{kj} \Delta_{iml} + \Delta^m_{ik} \Delta_{mjl}\right)
        \end{align}
    """
    )
    par.set_parval_from_str("enable_RbarDD_gridfunctions", True)
    rhs = BSSN_RHSs.BSSN_RHSs["Cartesian"]
    (
        trusted_alpha_rhs,
        trusted_vet_rhsU,
        trusted_bet_rhsU,
    ) = BSSN_gauge_RHSs.BSSN_gauge_RHSs()
    bssncon = BSSN_constraints.BSSN_constraints["Cartesian"]
    par.set_parval_from_str("enable_RbarDD_gridfunctions", False)
    # Clear BSSN_quantities.BSSN_quantities["Cartesian"], as it left Ricci symbolic.
    del BSSN_quantities.BSSN_quantities["Cartesian"]
    # Construct full symbolic expression for Ricci (RbarDD)
    Bq = BSSN_quantities.BSSN_quantities["Cartesian"]
    try:
        assert_equal(
            {
                "h_rhsDD": h_rhsDD,
                "cf_rhs": cf_rhs,
                "trK_rhs": trK_rhs,
                "Lambdabar_rhsU": Lambdabar_rhsU,
                "a_rhsDD": a_rhsDD,
                "alpha_rhs": alpha_rhs,
                "vet_rhsU": vet_rhsU,
                "bet_rhsU": bet_rhsU,
                "H": H,
                "MU": MU,
                "RbarDD": RbarDD,
            },
            {
                "h_rhsDD": rhs.h_rhsDD,
                "cf_rhs": rhs.cf_rhs,
                "trK_rhs": rhs.trK_rhs,
                "Lambdabar_rhsU": rhs.Lambdabar_rhsU,
                "a_rhsDD": rhs.a_rhsDD,
                "alpha_rhs": trusted_alpha_rhs,
                "vet_rhsU": trusted_vet_rhsU,
                "bet_rhsU": trusted_bet_rhsU,
                "H": bssncon.H,
                "MU": bssncon.MU,
                "RbarDD": cast(List[List[sp.Expr]], Bq.RbarDD),
            },
            suppress_message=False,
        )
    except AssertionError:
        return False
    return True


if __name__ == "__main__":
    if not test_example_BSSN():
        raise AssertionError(
            "Error: NRPy+ based BSSN expressions (Cartesian, Tmunu=False) "
            "disagree with NRPyLaTeX generated expressions"
        )

"""
Construct symbolic expressions for geodesic equations in various spacetime backgrounds.

This module can provide symbolic Christoffel symbols and geodesic equations of motion for numerical spacetime data.
It also can take a pre-defined analytic spacetime metric and compute the geometric quantities necessary for geodesic
integration, namely the Christoffel symbols and the geodesic equations of motion.
It supports both massive and photon particles.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

# Step 0.a: Import standard Python modules.
import logging
from typing import Dict, List, Optional, Tuple, cast

# Step 0.b: Import third-party modules
import sympy as sp
from typing_extensions import Literal

# Step 0.c: Import NRPy core modules
import nrpy.indexedexp as ixp
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)


class GeodesicEquations:
    """Generate and store symbolic expressions for geodesic motion from a given metric."""

    # mypy --strict requires class attributes to be declared.
    spacetime: str
    particle_type: str
    g4DD: List[List[sp.Expr]]
    xx: List[sp.Symbol]
    g4DD_dD: List[List[List[sp.Expr]]]
    Gamma4UDD: List[List[List[sp.Expr]]]
    geodesic_rhs: List[sp.Expr]
    u0_massive: Optional[sp.Expr]
    p0_photon: Optional[sp.Expr]
    norm_constraint_expr: sp.Expr
    Gamma4UDD_from_metric_recipe: List[List[List[sp.Expr]]]

    def __init__(self, spacetime: str, particle_type: str = "massive") -> None:
        """
        Initialize and compute all geodesic equation quantities.

        This constructor orchestrates the acquisition of the metric, the computation
        of its derivatives and Christoffel symbols, and the final generation of
        the geodesic equations of motion and Hamiltonian constraint for the
        specified particle type.

        :param spacetime: The spacetime to use (e.g., "KerrSchild_Cartesian").
        :param particle_type: The type of particle, either "massive" or "photon".
        :raises ValueError: If the particle type is not supported.
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Initialize optional attributes to None to satisfy type safety
        self.u0_massive = None
        self.p0_photon = None

        # Step 1: Acquire the metric and its associated coordinates.
        metric = Analytic_Spacetimes[spacetime]
        self.g4DD = metric.g4DD
        self.xx = metric.xx

        # Step 2: Compute geometric quantities derived from the metric.
        self.g4DD_dD = self.derivative_g4DD()
        self.Gamma4UDD = self.four_connections()

        # Step 3: Generate the right-hand-side of the geodesic ODEs and the
        #         Hamiltonian constraint equation based on the particle type.
        if self.particle_type == "massive":
            self.geodesic_rhs = self.geodesic_eom_rhs_massive()
            self.u0_massive = self.hamiltonian_constraint_massive()
        elif self.particle_type == "photon":
            self.geodesic_rhs = self.geodesic_eom_rhs_photon()
            self.p0_photon = self.hamiltonian_constraint_photon()
        else:
            raise ValueError(f"Particle type '{self.particle_type}' is not supported.")

        # Step 4: Reuse the precomputed symbolic recipe for Christoffel symbols
        #         built directly from generic metric data.
        self.Gamma4UDD_from_metric_recipe = _GAMMA4UDD_METRIC_RECIPE

        # Step 5: Generate generic normalization constraint for validation
        self.norm_constraint_expr = self.normalization_constraint()

    def derivative_g4DD(self) -> List[List[List[sp.Expr]]]:
        r"""
        Compute symbolic first derivatives of the metric tensor.

        This function calculates g_{\mu\nu,\alpha} = \partial_\alpha g_{\mu\nu}.
        It exploits the symmetry g_{\mu\nu} = g_{\nu\mu} by computing only the
        upper-triangular components of the derivative (where \nu >= \mu) and
        copying them to the lower-triangular part.

        :return: A 4x4x4 rank-3 tensor containing the derivatives.
        """
        g4DD_dD = ixp.zerorank3(dimension=4)
        # Loop over each coordinate to take the derivative
        for alpha in range(4):
            # Exploit symmetry by looping over the upper triangle of the metric
            for mu in range(4):
                for nu in range(mu, 4):
                    deriv = sp.diff(self.g4DD[mu][nu], self.xx[alpha])
                    g4DD_dD[mu][nu][alpha] = deriv
                    g4DD_dD[nu][mu][alpha] = deriv  # Set the symmetric component
        return g4DD_dD

    def four_connections(self) -> List[List[List[sp.Expr]]]:
        r"""
        Compute and simplify Christoffel symbols of the second kind, \Gamma^\alpha_{\mu\nu}.

        This function calculates the Christoffel symbols using the formula:
        \Gamma^\alpha_{\mu\nu} = (1/2) g^{\alpha\beta} (g_{\beta\mu,\nu} + g_{\beta\nu,\mu} - g_{\mu\nu,\beta})
        It exploits the symmetry \Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu} by computing
        only the upper-triangular components (where \nu >= \mu) and copying them.

        Reference:
        Wikipedia: Christoffel symbols
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Christoffel_symbols&oldid=1291002340#General_definition
        (See last equation in Section : General definition - Christoffel symbols of the second kind (symmetric definition))

        :return: A 4x4x4 rank-3 tensor containing the Christoffel symbols.
        """
        g4UU, _ = ixp.symm_matrix_inverter4x4(self.g4DD)
        Gamma4UDD = ixp.zerorank3(dimension=4)
        # Loop over the contravariant index alpha
        for alpha in range(4):
            # Exploit symmetry by looping over the upper triangle of the lower indices
            for mu in range(4):
                for nu in range(mu, 4):
                    term = sp.sympify(0)
                    # Sum over the dummy index beta
                    for beta in range(4):
                        term += (
                            sp.Rational(1, 2)
                            * g4UU[alpha][beta]
                            * (
                                self.g4DD_dD[nu][beta][mu]
                                + self.g4DD_dD[mu][beta][nu]
                                - self.g4DD_dD[mu][nu][beta]
                            )
                        )
                    Gamma4UDD[alpha][mu][nu] = term
                    Gamma4UDD[alpha][nu][mu] = term
        return Gamma4UDD

    def geodesic_eom_rhs_massive(self) -> List[sp.Expr]:
        r"""
        Generate the symbolic right-hand-side for the 8 massive geodesic ODEs.

        The equations of motion are: d(u^\alpha)/d(\tau) = -\Gamma^\alpha_{\mu\nu} u^\mu u^\nu
        This function exploits the symmetry of the Christoffel symbols in their
        lower indices to optimize the summation, doubling off-diagonal terms.


        Reference:
        Wikipedia: Geodesics in general relativity
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Geodesics_in_general_relativity&oldid=1320086873
        (See first equation in Section: Mathematical expression)

        :return: A list of 8 SymPy expressions for the RHS of the ODEs.
        """
        uU = ixp.declarerank1("uU", dimension=4)
        pos_rhs = [uU[0], uU[1], uU[2], uU[3]]

        Gamma4UDD = ixp.declarerank3("conn_Gamma4UDD", dimension=4, symmetry="sym12")
        vel_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            # Exploit symmetry by summing over the upper triangle of (mu, nu)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * uU[mu] * uU[nu]
                    if mu != nu:
                        term *= 2  # Double the off-diagonal terms
                    sum_term += term
            vel_rhs[alpha] = -sum_term

        return pos_rhs + vel_rhs

    def hamiltonian_constraint_massive(self) -> sp.Expr:
        r"""
        Symbolically derive u^0 from the Hamiltonian constraint for a massive particle.

        For this function geometrized units are used, i.e. G=c=1.
        This function solves the 4-velocity normalization condition:
        g_{\mu\nu} u^\mu u^\nu = -1
        Assuming the usual (-,+,+,+) signature with g_{00} < 0, this returns the positive u^0 root, conventional for forward time evolution.

        Reference:
        Wikipedia: Line element
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Line_element&oldid=1325490955
        (See last equation in Section: General formulation - Identification of the square of the line element with the metric tensor;
        Note: In the referenced equation, a parameter $g$ encodes the causal character of the curve; for massive particles (timelike curves) one has $g = -1$.)

                .. warning::
           This function assumes the particle is initialized in a region where g_{00} < 0 (e.g., outside the ergosphere).
           It cannot be used to set initial conditions inside the ergosphere where g_{00} > 0.
           Note: This limitation applies only to initialization. The time-evolution integrator handles transitions across the ergosphere correctly.

        :return: A SymPy expression for the positive root of u^0.
        """
        uU = ixp.declarerank1("uU", dimension=4)

        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)

        # The constraint is g_4DD00 (u^0)^2 + 2 g_4DD0i u^0 u^i + g_4DDij u^i u^j = -1
        # This is a quadratic equation: A (u^0)^2 + B u^0 + C = 0

        # A = g4DD_00
        A = g4DD[0][0]

        # B = 2 * sum(g_4DD0i * u^i) for i in {1,2,3}
        B = 2 * (g4DD[0][1] * uU[1] + g4DD[0][2] * uU[2] + g4DD[0][3] * uU[3])

        # C includes the +1 from moving the -1 to the LHS
        # C = 1 + sum(g_4DDij * u^i * u^j)
        # We manually unroll the loops to ensure we only access the upper triangle (i <= j)

        # Diagonal terms (i == j)
        term_diagonals = (
            g4DD[1][1] * uU[1] ** 2 + g4DD[2][2] * uU[2] ** 2 + g4DD[3][3] * uU[3] ** 2
        )

        # Off-diagonal terms (i != j) multiplied by 2 because g_4DDij is symmetric
        term_off_diagonals = 2 * (
            g4DD[1][2] * uU[1] * uU[2]
            + g4DD[1][3] * uU[1] * uU[3]
            + g4DD[2][3] * uU[2] * uU[3]
        )

        C = sp.sympify(1) + term_diagonals + term_off_diagonals

        # Quadratic formula: (-B +/- sqrt(B^2 - 4AC)) / 2A
        discriminant = sp.sqrt(B**2 - 4 * A * C)
        sol1 = (-B + discriminant) / (2 * A)
        sol2 = (-B - discriminant) / (2 * A)
        solutions = [sol1, sol2]
        return cast(sp.Expr, solutions[1])

    def geodesic_eom_rhs_photon_christoffel(self) -> List[sp.Expr]:
        r"""
        Generate the symbolic right-hand-side for the 9 photon geodesic ODEs.

        The equations of motion are: d(p^\alpha)/d(\kappa) = -\Gamma^\alpha_{\mu\nu} p^\mu p^\nu.
        The ninth ODE provides a diagnostic for the spatial distance traveled by the photon,
        as measured by a local Eulerian observer.

        **Equation 9: Eulerian Distance**
        The rate of change of the proper spatial distance L_Euler with respect
        to the affine parameter \lambda is given by:

            dL_Euler/d\lambda = \alpha p^0 = p^0 / sqrt(-g^00)

        **Physical Interpretation:**
        This value represents the "lab frame" distance measured by an **Eulerian Observer**
        (also known as a Normal Observer). This observer has a 4-velocity n^\mu normal
        to the spatial hypersurface (slice of constant coordinate time t).

        * **Observer:** The Eulerian observer is "hovering" at a specific spatial coordinate,
            moving orthogonal to the spatial slice. In the ADM formalism, their 4-velocity
            is n_\mu = (-\alpha, 0, 0, 0).
        * **Measurement:** It answers the question: "How much distance does the local
            Eulerian observer see the photon cover in this instant?"
        * **Invariance:** This scalar is spatially invariant (valid under spatial rotations)
            but slicing dependent (depends on the definition of time t).

        **Variables:**
        * p^0: Time component of the photon 4-momentum (dt/d\lambda).
        * g^00: Time-time component of the **inverse** metric tensor.
        * \alpha: The Lapse function, defined as \alpha = 1/sqrt(-g^00) (assuming g^00 < 0).

        .. note::
            This formulation assumes a metric signature where g^00 < 0 (e.g., -+++).
            It is well-behaved in horizon-penetrating coordinates (like Kerr-Schild) where
            g^00 remains finite, but may diverge in static coordinates (like Schwarzschild)
            at the horizon where the Eulerian observer cannot exist.

        Reference:
        Wikipedia: Geodesics in general relativity
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Geodesics_in_general_relativity&oldid=1320086873
        (See first equation in Section: Mathematical expression)

        :return: A list of 9 SymPy expressions for the RHS of the ODEs.
        """
        pU = ixp.declarerank1("pU", dimension=4)
        pos_rhs = [pU[0], pU[1], pU[2], pU[3]]

        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)

        # We need the inverse metric g^UU to compute g^00 for the lapse function.
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)

        Gamma4UDD = ixp.declarerank3("conn_Gamma4UDD", dimension=4, symmetry="sym12")
        mom_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            # Exploit symmetry by summing over the upper triangle of (mu, nu)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * pU[mu] * pU[nu]
                    if mu != nu:
                        term *= 2
                    sum_term += term
            mom_rhs[alpha] = -sum_term

        # Ninth Equation: Eulerian Distance Evolution
        # dL_Euler/dlambda = p^0 / sqrt(-g^00) = alpha * p^0
        # Note: L_Euler is a *signed* distance. In reverse ray tracing, p^0 < 0
        # generally causes L_Euler to decrease, except inside an ergosphere
        # where p^0 > 0 causes it to increase.
        path_len_rhs = [pU[0] / sp.sqrt(-g4UU[0][0])]

        return pos_rhs + mom_rhs + path_len_rhs

    def geodesic_eom_rhs_photon(self) -> List[sp.Expr]:
        r"""
        Generate the symbolic right-hand-side for the 9 photon geodesic ODEs.

        The equations of motion are:

            d(x^\alpha)/d(\lambda) = p^\alpha

            d(p^\alpha)/d(\lambda)
                = -\Gamma^\alpha_{\mu\nu} p^\mu p^\nu

        where the Christoffel symbols are computed from the 4-metric and its
        partial derivatives:

            \Gamma^\alpha_{\mu\nu}
                = 1/2 g^{\alpha\beta}
                  (g_{\beta\nu,\mu} + g_{\beta\mu,\nu} - g_{\mu\nu,\beta})

        The ninth ODE provides a diagnostic for the spatial distance traveled
        by the photon, as measured by a local Eulerian observer:

            dL_Euler/d\lambda = p^0 / sqrt(-g^00)

        :return: A list of 9 SymPy expressions for the RHS of the ODEs.
        """
        pU = ixp.declarerank1("pU", dimension=4)
        pos_rhs = [pU[0], pU[1], pU[2], pU[3]]

        g4DD = ixp.declarerank2(
            "metric_g4DD",
            symmetry="sym01",
            dimension=4,
        )

        # g4DD_dD[mu][nu][sigma] = partial_sigma g_{mu nu}
        g4DD_dD = ixp.declarerank3(
            "metric_g4DD_dD",
            symmetry="sym01",
            dimension=4,
        )

        # Compute the inverse metric g^{mu nu}.
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)

        # Compute the Christoffel symbols from g4DD and g4DD_dD.
        Gamma4UDD = ixp.zerorank3(dimension=4)

        for alpha in range(4):
            for mu in range(4):
                for nu in range(4):
                    sum_term = sp.sympify(0)

                    for beta in range(4):
                        sum_term += g4UU[alpha][beta] * (
                            g4DD_dD[beta][nu][mu]
                            + g4DD_dD[beta][mu][nu]
                            - g4DD_dD[mu][nu][beta]
                        )

                    Gamma4UDD[alpha][mu][nu] = sp.Rational(1, 2) * sum_term

        # Compute the four momentum equations.
        mom_rhs = ixp.zerorank1(dimension=4)

        for alpha in range(4):
            sum_term = sp.sympify(0)

            # Exploit symmetry by summing over the upper triangle of (mu, nu).
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * pU[mu] * pU[nu]

                    if mu != nu:
                        term *= 2

                    sum_term += term

            mom_rhs[alpha] = -sum_term

        # Ninth Equation: Eulerian Distance Evolution
        # dL_Euler/dlambda = p^0 / sqrt(-g^00) = alpha * p^0
        path_len_rhs = [pU[0] / sp.sqrt(-g4UU[0][0])]

        return pos_rhs + mom_rhs + path_len_rhs

    def geodesic_eom_rhs_photon_normalized(self) -> List[sp.Expr]:
        r"""
        Generate the normalized 9-component photon RHS used by the numerical build.

        This formulation is a normalized photon evolution
        system, using coordinate time ``t`` as the independent variable and the
        state

        ``(t, x, y, z, u, \Pi_1, \Pi_2, \Pi_3, \lambda)``,

        where

        ``u = \ln|\alpha p^0|``

        and

        ``\Pi_i = p_i / (\alpha p^0)``.

        The ninth equation evolves the affine parameter along the past-directed
        branch conventional for reverse ray tracing:

        ``d\lambda / dt = -\alpha e^{-u}``.

        All geometric quantities are constructed directly from the generic
        placeholder four-metric ``metric_g4DD`` and its first derivatives
        ``metric_g4DD_dD``.

        Reference:
        Bohn et al., "What does a binary black hole merger look like?"
        (See Eqs. (4) and (5))

        :return: A list of 9 SymPy expressions for the normalized photon RHS.
        """
        g4DD = ixp.declarerank2(
            "metric_g4DD",
            symmetry="sym01",
            dimension=4,
        )

        # g4DD_dD[mu][nu][sigma] = partial_sigma g_{mu nu}
        g4DD_dD = ixp.declarerank3(
            "metric_g4DD_dD",
            symmetry="sym01",
            dimension=4,
        )

        # The normalized variables use the lapse, shift, and inverse spatial
        # metric. These are algebraic combinations of the supplied four-metric,
        # not additional numerical payload fields.
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)
        alpha = sp.sympify(1) / sp.sqrt(-g4UU[0][0])

        gammaDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                gammaDD[i][j] = g4DD[i + 1][j + 1]
        gammaUU, _ = ixp.symm_matrix_inverter3x3(gammaDD)

        betaU = ixp.zerorank1(dimension=3)
        for i in range(3):
            betaU[i] = -g4UU[0][i + 1] / g4UU[0][0]

        # Evolved variables for the normalized photon formulation:
        # f[0] = t
        # f[1:4] = x^i
        # f[4] = u = ln|alpha p^0|
        # f[5:8] = Pi_i = p_i / (alpha p^0)
        # f[8] = lambda
        u = sp.Symbol("u", real=True)
        PiD = ixp.declarerank1("PiD", dimension=3)

        # Raise the covariant normalized momentum:
        # Pi^i = gamma^ij Pi_j
        PiU = ixp.zerorank1(dimension=3)
        for i in range(3):
            PiU[i] = sp.sympify(0)
            for j in range(3):
                PiU[i] += gammaUU[i][j] * PiD[j]

        # Coordinate-time tangent:
        # V^mu = dx^mu / dt = p^mu / p^0
        # V^0 = 1
        # V^i = alpha Pi^i - beta^i
        coordinate_velocityU = ixp.zerorank1(dimension=4)
        coordinate_velocityU[0] = sp.sympify(1)

        pos_rhs = []
        for i in range(3):
            coordinate_velocityU[i + 1] = alpha * PiU[i] - betaU[i]
            pos_rhs.append(coordinate_velocityU[i + 1])

        # Normalized contravariant momentum:
        # p^mu / (alpha p^0) = V^mu / alpha
        normalized_pU = ixp.zerorank1(dimension=4)
        for mu in range(4):
            normalized_pU[mu] = coordinate_velocityU[mu] / alpha

        # Differentiate alpha = (-g^00)^(-1/2) directly from the metric:
        # alpha_,sigma = (1/2) alpha^3 (g^00)_,sigma
        # (g^00)_,sigma = -g^0mu g^0nu g_munu,sigma
        alpha_dD = ixp.zerorank1(dimension=4)
        for sigma in range(4):
            g4UU00_dD = sp.sympify(0)
            for mu in range(4):
                for nu in range(4):
                    g4UU00_dD -= g4UU[0][mu] * g4UU[0][nu] * g4DD_dD[mu][nu][sigma]
            alpha_dD[sigma] = sp.Rational(1, 2) * alpha**3 * g4UU00_dD

        # The time component of the contravariant geodesic equation gives
        # d/dt ln|p^0| = -Gamma^0_munu V^mu V^nu. Construct only the required
        # contraction directly from g_munu and its first derivatives.
        Gamma0_contract = sp.sympify(0)
        for mu in range(4):
            for nu in range(mu, 4):
                Gamma0_munu = sp.sympify(0)
                for rho in range(4):
                    Gamma0_munu += (
                        sp.Rational(1, 2)
                        * g4UU[0][rho]
                        * (
                            g4DD_dD[rho][nu][mu]
                            + g4DD_dD[rho][mu][nu]
                            - g4DD_dD[mu][nu][rho]
                        )
                    )

                term = Gamma0_munu * coordinate_velocityU[mu] * coordinate_velocityU[nu]
                if mu != nu:
                    term *= 2
                Gamma0_contract += term

        # u = ln|alpha p^0|, so
        # du/dt = V^mu alpha_,mu / alpha
        #         - Gamma^0_munu V^mu V^nu
        alpha_transport = sp.sympify(0)
        for mu in range(4):
            alpha_transport += alpha_dD[mu] * coordinate_velocityU[mu]

        u_rhs_expr = alpha_transport / alpha - Gamma0_contract
        u_rhs = [u_rhs_expr]

        # The covariant geodesic equation can be written as
        # dp_i/dlambda = (1/2) g_munu,i p^mu p^nu. Therefore,
        # dPi_i/dt = (alpha/2) g_munu,i
        #             [p^mu/(alpha p^0)] [p^nu/(alpha p^0)]
        #           - (du/dt) Pi_i.
        Pi_rhs = ixp.zerorank1(dimension=3)
        for i in range(3):
            metric_derivative_contract = sp.sympify(0)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = (
                        g4DD_dD[mu][nu][i + 1] * normalized_pU[mu] * normalized_pU[nu]
                    )
                    if mu != nu:
                        term *= 2
                    metric_derivative_contract += term

            Pi_rhs[i] = (
                sp.Rational(1, 2) * alpha * metric_derivative_contract
                - u_rhs_expr * PiD[i]
            )

        # First Equation: Coordinate-Time Evolution
        time_rhs = [sp.sympify(1)]

        # Ninth Equation: Affine Parameter Evolution
        # Reverse ray tracing uses the past-directed branch:
        # dlambda/dt = -alpha exp(-u)
        lambda_rhs = [-alpha * sp.exp(-u)]

        return time_rhs + pos_rhs + u_rhs + Pi_rhs + lambda_rhs

    def normalization_constraint_photon_normalized(self) -> sp.Expr:
        r"""
        Generate the normalized photon null constraint from the four-metric.

        The normalized photon formulation evolves the covariant spatial momentum
        ``\Pi_i`` and should satisfy

        ``\gamma^{ij} \Pi_i \Pi_j = 1``,

        where ``\gamma^{ij}`` is the inverse of the spatial block ``g_ij`` of
        the generic placeholder four-metric ``metric_g4DD``.

        This helper returns the left-hand side of that relation. A normalized-
        photon diagnostic should therefore compare the returned scalar against
        ``1``.

        :return: A SymPy expression for ``\gamma^{ij} \Pi_i \Pi_j``.
        """
        g4DD = ixp.declarerank2(
            "metric_g4DD",
            symmetry="sym01",
            dimension=4,
        )

        # The spatial metric is the spatial block of the four-metric:
        # gamma_ij = g_ij
        gammaDD = ixp.zerorank2(dimension=3)
        for i in range(3):
            for j in range(3):
                gammaDD[i][j] = g4DD[i + 1][j + 1]
        gammaUU, _ = ixp.symm_matrix_inverter3x3(gammaDD)

        PiD = ixp.declarerank1("PiD", dimension=3)

        constraint = sp.sympify(0)
        for i in range(3):
            for j in range(i, 3):
                term = gammaUU[i][j] * PiD[i] * PiD[j]
                if i != j:
                    term *= 2
                constraint += term

        return constraint

    def hamiltonian_constraint_photon(self) -> sp.Expr:
        r"""
        Symbolically derive p^0 from the Hamiltonian constraint for a photon particle.

        This function solves the null geodesic condition: g_{\mu\nu} p^\mu p^\nu = 0.
        Assuming the usual (-,+,+,+) signature with g_{00} < 0, this returns the negative p^0 root, conventional for reverse ray tracing.

        Reference:
        Wikipedia: Line element
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Line_element&oldid=1325490955
        (See last equation in Section: General formulation - Identification of the square of the line element with the metric tensor;
         Note: In the referenced equation, null (lightlike) curves correspond to $g = 0$, which is the case for photon particles.)

        .. warning::
           This function assumes the particle is initialized in a region where g_{00} < 0 (e.g., outside the ergosphere).
           It cannot be used to set initial conditions inside the ergosphere where g_{00} > 0.
           Note: This limitation applies only to initialization. The time-evolution integrator handles transitions across the ergosphere correctly.

        :return: A SymPy expression for the negative root of p^0.
        """
        pU = ixp.declarerank1("pU", dimension=4)

        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)

        # The constraint is g_00 (p^0)^2 + 2 g_0i p^0 p^i + g_ij p^i p^j = 0
        # This is a quadratic equation: A (p^0)^2 + B p^0 + C = 0

        # A = g4DD_00
        A = g4DD[0][0]

        # B = 2 * sum(g_4DD0i * p^i) for i in {1,2,3}
        B = 2 * (g4DD[0][1] * pU[1] + g4DD[0][2] * pU[2] + g4DD[0][3] * pU[3])

        # C = sum(g_4DDij * p^i * p^j)
        # We manually unroll the loops to ensure we only access the upper triangle (i <= j)

        # Diagonal terms (i == j)
        term_diagonals = (
            g4DD[1][1] * pU[1] ** 2 + g4DD[2][2] * pU[2] ** 2 + g4DD[3][3] * pU[3] ** 2
        )

        # Off-diagonal terms (i != j) multiplied by 2 because g_4DDij is symmetric
        term_off_diagonals = 2 * (
            g4DD[1][2] * pU[1] * pU[2]
            + g4DD[1][3] * pU[1] * pU[3]
            + g4DD[2][3] * pU[2] * pU[3]
        )

        C = term_diagonals + term_off_diagonals

        # Quadratic formula: (-B +/- sqrt(B^2 - 4AC)) / 2A
        discriminant = sp.sqrt(B**2 - 4 * A * C)
        sol1 = (-B + discriminant) / (2 * A)
        sol2 = (-B - discriminant) / (2 * A)
        solutions = [sol1, sol2]
        return cast(sp.Expr, solutions[0])

    def normalization_constraint(self) -> sp.Expr:
        r"""
        Generate the symbolic expression for the metric normalization constraint.

        This computes the scalar invariant: C = g_{\mu\nu} v^\mu v^\nu.

        The variable 'vU' represents the tangent vector to the geodesic curve:
        - For massive particles, vU corresponds to 4-velocity (u^\mu).
            Expected Result: -1 (in -+++ signature).
        - For photon particles, vU corresponds to 4-momentum (p^\mu).
            Expected Result: 0.

        Reference:
        Wikipedia: Line element
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Line_element&oldid=1325490955
        (See last equation in Section: General formulation - Identification of the square of the line element with the metric tensor;
            Note: In the referenced equation, null (lightlike) curves correspond to $g = 0$, which is the case for photon particles.
            Note: In the referenced equation for massive particles (timelike curves) $g = -1$.)

        .. warning:: This function is meant for validation purposes

        :return: A SymPy expression for the contraction g_{\mu\nu} v^\mu v^\nu.
        """
        # Generic tangent vector v^mu (vU0, vU1, vU2, vU3)
        vU = ixp.declarerank1("vU", dimension=4)

        # Generic metric g_mu_nu
        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)

        constraint = sp.sympify(0)

        # Loop over indices, exploiting g_{mu,nu} = g_{nu,mu} symmetry
        for mu in range(4):
            for nu in range(mu, 4):
                term = g4DD[mu][nu] * vU[mu] * vU[nu]

                # Double the off-diagonal terms (mu != nu)
                if mu != nu:
                    term *= 2

                constraint += term

        return constraint

    @staticmethod
    def photon_momentum_to_normalized_quantities() -> Tuple[sp.Expr, List[sp.Expr]]:
        r"""
        Convert direct photon momentum variables to normalized photon quantities.

        This helper is intended for the numerical photon initialization pipeline,
        where a one-time conversion is needed from the legacy direct four-momentum
        state ``(p^0, p^1, p^2, p^3)`` to the normalized variables used by the
        Bohn-style coordinate-time evolution system:

        ``u = \ln|\alpha p^0|``

        ``\Pi_i = p_i / (\alpha p^0)``

        The conversion uses only the generic placeholder metric ``metric_g4DD``
        and direct four-momentum ``pU`` so it can be consumed by a dedicated
        photon-side infrastructure helper that operates after the initial camera
        geometry has been seeded and the local metric has been interpolated.

        :return: Tuple containing ``u`` and the 3-component covariant normalized
            spatial momentum ``PiD``.
        """
        pU = ixp.declarerank1("pU", dimension=4)
        g4DD = ixp.declarerank2("metric_g4DD", symmetry="sym01", dimension=4)
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)

        # ADM lapse from the inverse four-metric:
        # alpha = 1 / sqrt(-g^00)
        alpha = sp.sympify(1) / sp.sqrt(-g4UU[0][0])

        # Lower the spatial momentum index:
        # p_i = g_i0 p^0 + g_ij p^j
        pD = ixp.zerorank1(dimension=3)
        for i in range(3):
            pD[i] = g4DD[i + 1][0] * pU[0]
            for j in range(3):
                pD[i] += g4DD[i + 1][j + 1] * pU[j + 1]

        # The normalization denominator shared by u and Pi_i:
        # alpha p^0
        alpha_p0 = alpha * pU[0]

        # Normalized variables:
        # u = ln|alpha p^0|
        # Pi_i = p_i / (alpha p^0)
        u = sp.log(sp.Abs(alpha_p0))
        PiD = ixp.zerorank1(dimension=3)
        for i in range(3):
            PiD[i] = pD[i] / alpha_p0

        return u, PiD

    @staticmethod
    def symbolic_g4DD_recipe_from_bssn_grid_basis(
        bssn_coord_system: str,
        target_basis: Literal["Cartesian", "Spherical"] = "Cartesian",
        enable_bssn_rfm_precompute: bool = False,
    ) -> List[List[sp.Expr]]:
        r"""
        Generate a transformed covariant four-metric recipe from BSSN grid-basis data.

        This routine specializes to the NRPy BSSN/reference-metric data model.
        The spacetime map is assumed to be a time-independent spatial map lifted
        into 4D with identity time:

        ``x^0 = u^0``

        ``x^i = X^i(u^1, u^2, u^3)``

        Here ``u^1``, ``u^2``, and ``u^3`` are the spatial grid coordinates
        stored in ``ReferenceMetric.xx``. The target spatial map ``X^i`` is
        selected by ``target_basis`` from the matching ``ReferenceMetric``.
        This routine does not rebuild Jacobians from ``xx_to_Cart`` or
        ``xxSph``. Instead, it consumes the canonical inverse-Jacobian objects
        stored on ``ReferenceMetric``.

        Since the four-metric is a covariant rank-2 tensor, the transformation is

        ``g'_{\mu\nu} = K^a{}_\mu K^b{}_\nu g_{ab}``,

        where ``K4UD[a][mu] = \partial u^a / \partial x^\mu`` is the inverse
        Jacobian from target-basis coordinates back to grid-basis coordinates.
        All mixed time-space Jacobian entries vanish except ``K4UD[0][0] = 1``.
        Unlike the Christoffel transformation, no forward Jacobian derivative or
        inhomogeneous connection term appears.

        :param bssn_coord_system: BSSN coordinate system used to construct the
                                  source grid-basis four-metric.
        :param target_basis: Target spatial basis. Supported values are
                             ``"Cartesian"`` and ``"Spherical"``.
        :param enable_bssn_rfm_precompute: Whether to enable reference-metric
                                           precomputation inside the BSSN helper.
        :return: A 4x4 list of SymPy expressions for the transformed covariant
                 four-metric.
        :raises ValueError: If ``bssn_coord_system`` is not a string.
        :raises ValueError: If ``target_basis`` is unsupported.
        :raises ValueError: If ``enable_bssn_rfm_precompute`` is not a bool.
        :raises ValueError: If ``bssn_coord_system`` is unsupported.
        :raises ValueError: If the matching reference metric is not 3D.
        """
        # Step 1.a: Validate the public API inputs before triggering cached lookups.
        if not isinstance(bssn_coord_system, str):
            raise ValueError(
                "bssn_coord_system must be a string, "
                f"got {type(bssn_coord_system).__name__}"
            )
        if target_basis not in ("Cartesian", "Spherical"):
            raise ValueError(
                f"Unsupported target_basis '{target_basis}'. "
                "Supported values are 'Cartesian' and 'Spherical'."
            )
        if not isinstance(enable_bssn_rfm_precompute, bool):
            raise ValueError(
                "enable_bssn_rfm_precompute must be a bool, "
                f"got {type(enable_bssn_rfm_precompute).__name__}"
            )

        # Step 1.b: Import the BSSN and reference-metric helpers lazily so
        #           analytic-only users do not pull in the full BSSN stack at import time.
        # pylint: disable-next=import-outside-toplevel,redefined-outer-name
        import nrpy.reference_metric as refmetric

        # pylint: disable-next=import-outside-toplevel
        from nrpy.equations.general_relativity.g4munu_conversions import (
            BSSN_to_g4DD,
        )

        if bssn_coord_system not in refmetric.supported_CoordSystems:
            raise ValueError(
                f"Unsupported CoordSystem '{bssn_coord_system}'. "
                f"Supported coordinate systems: {refmetric.supported_CoordSystems}"
            )

        # Step 2: Construct the source-basis four-metric g_{ab} directly from
        #         the BSSN variables native to the requested grid basis.
        grid_g4DD = BSSN_to_g4DD(
            CoordSystem=bssn_coord_system,
            enable_rfm_precompute=enable_bssn_rfm_precompute,
        )

        # Step 3: Fetch the matching ReferenceMetric object so the spatial
        #         basis transformation uses the same coordinate model as the
        #         underlying BSSN quantities.
        coord_system_key = bssn_coord_system + (
            "_rfm_precompute" if enable_bssn_rfm_precompute else ""
        )
        rfm = refmetric.reference_metric[coord_system_key]
        if len(rfm.xx) != 3:
            raise ValueError(
                "ReferenceMetric coordinates are inconsistent with 3D space."
            )

        # Step 4: Select the canonical inverse spatial Jacobian
        #         K^j_i = partial u^j / partial x^i for the requested target basis.
        #         Since g_{mu nu} is a covariant tensor, the forward Jacobian and
        #         its derivatives are not needed.
        if target_basis == "Cartesian":
            spatial_K3UD = rfm.Jac_dUrfm_dDCartUD
        else:
            spatial_K3UD = rfm.Jac_dUrfm_dDSphUD

        # Step 5: Lift the time-independent inverse spatial Jacobian into the
        #         4D map x^0 = u^0, x^i = X^i(u^j). Thus K^0_0 = 1 and all
        #         mixed time-space entries vanish.
        K4UD = ixp.zerorank2(dimension=4)
        K4UD[0][0] = sp.sympify(1)
        for spatial_grid_index in range(3):
            for spatial_target_index in range(3):
                K4UD[spatial_grid_index + 1][spatial_target_index + 1] = spatial_K3UD[
                    spatial_grid_index
                ][spatial_target_index]

        # Step 6: Apply the covariant tensor transformation law:
        #         g'_{mu nu} = K^a_mu K^b_nu g_{ab}.
        return GeodesicEquations._transform_covariant_metric_from_grid_basis_data(
            grid_g4DD=grid_g4DD, K4UD=K4UD
        )

    @staticmethod
    def symbolic_numerical_christoffel_recipe() -> List[List[List[sp.Expr]]]:
        r"""
        Generate a symbolic recipe for Christoffel symbols from generic metric data.

        This function exploits the symmetry \Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu}
        by computing only the upper-triangular components (where \nu >= \mu).

        Reference:
        Wikipedia: Christoffel symbols
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Christoffel_symbols&oldid=1291002340#General_definition
        (See last equation in Section : General definition - Christoffel symbols of the second kind (symmetric definition))

        :return: A list of lists of lists of SymPy expressions representing the recipe.
        """
        g4DD_sym = ixp.declarerank2("g4DD_sym", dimension=4)
        g4DD_dD_sym = ixp.declarerank3("g4DD_dD_sym", dimension=4)
        g4UU_sym, _ = ixp.symm_matrix_inverter4x4(g4DD_sym)

        Gamma4UDD_recipe = ixp.zerorank3(dimension=4)
        for alpha in range(4):
            for mu in range(4):
                for nu in range(mu, 4):
                    term = sp.sympify(0)
                    for beta in range(4):
                        term += (
                            sp.Rational(1, 2)
                            * g4UU_sym[alpha][beta]
                            * (
                                g4DD_dD_sym[beta][mu][nu]
                                + g4DD_dD_sym[beta][nu][mu]
                                - g4DD_dD_sym[mu][nu][beta]
                            )
                        )
                    Gamma4UDD_recipe[alpha][mu][nu] = term
                    Gamma4UDD_recipe[alpha][nu][mu] = term
        return Gamma4UDD_recipe

    @staticmethod
    def _transform_covariant_metric_from_grid_basis_data(
        grid_g4DD: List[List[sp.Expr]], K4UD: List[List[sp.Expr]]
    ) -> List[List[sp.Expr]]:
        r"""
        Transform a covariant four-metric using a supplied inverse Jacobian recipe.

        This routine implements the coordinate-transformation law

        ``g'_{mu nu} = K^a{}_{mu} K^b{}_{nu} g_{ab}``,

        where ``K4UD[a][mu] = \partial u^a / \partial x^\mu`` is the inverse
        Jacobian from target-basis coordinates back to the source grid basis.

        :param grid_g4DD: Covariant four-metric in the source grid basis.
        :param K4UD: Inverse Jacobian from target basis to grid basis.
        :return: A 4x4 list of SymPy expressions for the transformed covariant
                 four-metric.
        """
        # Step 1: Loop over the target-basis covariant indices. Since the
        #         output metric is symmetric, compute only the upper triangle
        #         and mirror each entry to the lower triangle.
        target_g4DD = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(mu, 4):
                term = sp.sympify(0)

                # Step 2: Contract the source-basis metric with one inverse
                #         Jacobian factor for each covariant index:
                #         g'_{mu nu} = K^a_mu K^b_nu g_{ab}.
                #         The dummy indices grid_a and grid_b run over the
                #         source-basis spacetime coordinates.
                for grid_a in range(4):
                    for grid_b in range(4):
                        term += (
                            K4UD[grid_a][mu]
                            * K4UD[grid_b][nu]
                            * grid_g4DD[grid_a][grid_b]
                        )

                # Step 3: Store the transformed component and copy it across
                #         the diagonal to preserve explicit symmetry.
                target_g4DD[mu][nu] = term
                target_g4DD[nu][mu] = term
        return target_g4DD

    @staticmethod
    def _transform_christoffel_recipe_from_grid_basis_data(
        grid_Gamma4UDD: List[List[List[sp.Expr]]],
        J4UD: List[List[sp.Expr]],
        J4UD_dD: List[List[List[sp.Expr]]],
        K4UD: Optional[List[List[sp.Expr]]] = None,
    ) -> List[List[List[sp.Expr]]]:
        r"""
        Transform grid-basis Christoffel symbols using a supplied Jacobian recipe.

        This routine implements the coordinate-transformation law

        \Gamma^\alpha{}_{\mu\nu}
            =
            J^\alpha{}_a K^b{}_\mu K^c{}_\nu \Gamma^a{}_{bc}
            -
            K^b{}_\mu K^c{}_\nu \partial_b J^\alpha{}_c,

        where ``J4UD[alpha][a] = \partial x^\alpha / \partial u^a`` is the
        grid-to-target Jacobian and ``K4UD[a][mu] = \partial u^a / \partial x^\mu``
        is its inverse.

        The derivative tensor is indexed as
        ``J4UD_dD[alpha][c][b] = \partial_b J^\alpha{}_c``.

        If ``K4UD`` is supplied, this routine consumes it directly as the
        inverse Jacobian associated with ``J4UD``. As an internal helper, it
        assumes the Christoffel and Jacobian inputs from module-local callers
        already satisfy the expected 4D tensor structure.

        :param grid_Gamma4UDD: Christoffel symbols in the source grid basis.
        :param J4UD: Grid-to-target Jacobian.
        :param J4UD_dD: Derivatives of the grid-to-target Jacobian.
        :param K4UD: Inverse Jacobian. If ``None``, compute it from ``J4UD``.
        :return: A 4x4x4 list of SymPy expressions for the transformed Christoffel
                 symbols.
        """
        # Step 1: Invert the non-symmetric grid-to-target Jacobian if needed.
        if K4UD is None:
            K4UD, _ = ixp.generic_matrix_inverter4x4(J4UD)

        # Step 2: Construct the target-basis Christoffels from the standard
        #         coordinate-transformation law.
        Gamma4UDD_recipe = ixp.zerorank3(dimension=4)
        for alpha in range(4):
            for mu in range(4):
                for nu in range(mu, 4):
                    term = sp.sympify(0)

                    # Step 2.a: Add the tensorial pullback term
                    #           J^alpha_a K^b_mu K^c_nu Gamma^a_bc.
                    for grid_a in range(4):
                        for grid_b in range(4):
                            for grid_c in range(4):
                                term += (
                                    J4UD[alpha][grid_a]
                                    * K4UD[grid_b][mu]
                                    * K4UD[grid_c][nu]
                                    * grid_Gamma4UDD[grid_a][grid_b][grid_c]
                                )

                    # Step 2.b: Subtract the inhomogeneous connection term
                    #           K^b_mu K^c_nu partial_b J^alpha_c.
                    for grid_b in range(4):
                        for grid_c in range(4):
                            term -= (
                                K4UD[grid_b][mu]
                                * K4UD[grid_c][nu]
                                * J4UD_dD[alpha][grid_c][grid_b]
                            )

                    Gamma4UDD_recipe[alpha][mu][nu] = term
                    Gamma4UDD_recipe[alpha][nu][mu] = term
        return Gamma4UDD_recipe

    @staticmethod
    def symbolic_christoffel_recipe_from_grid_basis() -> List[List[List[sp.Expr]]]:
        r"""
        Generate the fully generic symbolic transformed-Christoffel recipe.

        This routine declares symbolic grid-basis Christoffel symbols
        ``grid_Gamma4UDD``, a symbolic 4D Jacobian ``J4UD``, and a symbolic
        Hessian ``J4UD_dD`` for an arbitrary smooth spacetime coordinate map
        ``x^\alpha(u^a)``. It therefore returns a purely symbolic recipe for
        the transformed connection, without assuming any particular spacetime,
        coordinate system, or reference metric.

        :return: A 4x4x4 list of SymPy expressions for the transformed Christoffel
                 symbols.
        """
        # Step 1: Declare symbolic source-basis Christoffels
        #         Gamma^a_bc in an arbitrary smooth grid basis u^a.
        grid_Gamma4UDD = ixp.declarerank3(
            "grid_Gamma4UDD", dimension=4, symmetry="sym12"
        )

        # Step 2: Declare the fully generic Jacobian J^alpha_a =
        #         partial x^alpha / partial u^a.
        J4UD = ixp.declarerank2("J4UD", dimension=4, symmetry="nosym")

        # Step 3: Declare the Jacobian derivatives
        #         partial_b J^alpha_c = partial^2 x^alpha / partial u^b partial u^c.
        #         Because these are mixed partials of one smooth coordinate map,
        #         the last two indices are symmetric.
        J4UD_dD = ixp.declarerank3("J4UD_dD", dimension=4, symmetry="sym12")

        # Step 4: Apply the general connection-transformation formula.
        return GeodesicEquations._transform_christoffel_recipe_from_grid_basis_data(
            grid_Gamma4UDD=grid_Gamma4UDD,
            J4UD=J4UD,
            J4UD_dD=J4UD_dD,
        )

    @staticmethod
    def symbolic_christoffel_recipe_from_bssn_grid_basis(
        bssn_coord_system: str,
        target_basis: Literal["Cartesian", "Spherical"] = "Cartesian",
        enable_bssn_rfm_precompute: bool = False,
    ) -> List[List[List[sp.Expr]]]:
        r"""
        Generate a transformed-Christoffel recipe from BSSN/reference-metric data.

        This routine specializes to the NRPy BSSN/reference-metric data model.
        The spacetime map is assumed to be a time-independent spatial map lifted
        into 4D with identity time:

        ``x^0 = u^0``

        ``x^i = X^i(u^1, u^2, u^3)``

        Here ``u^1``, ``u^2``, and ``u^3`` are the spatial grid coordinates
        stored in ``ReferenceMetric.xx``. The target spatial map ``X^i`` is
        selected by ``target_basis`` from the matching ``ReferenceMetric``.
        This routine does not rebuild those Jacobians from ``xx_to_Cart`` or
        ``xxSph``. Instead, it consumes the public canonical Jacobian objects
        already stored on ``ReferenceMetric``:

        - ``"Cartesian"`` uses the public Jacobians
          ``rfm.Jac_dUCart_dDrfmUD`` and ``rfm.Jac_dUrfm_dDCartUD``
        - ``"Spherical"`` uses the public Jacobians
          ``rfm.Jac_dUSph_dDrfmUD`` and ``rfm.Jac_dUrfm_dDSphUD``

        All mixed time-space Jacobian entries vanish except ``J4UD[0][0] = 1``,
        and all second derivatives involving time vanish as well. This branch
        therefore supports only time-independent spatial maps lifted into
        spacetime by identity time. It does not represent moving grids,
        time-dependent spacetime maps, or general non-coordinate basis
        transformations.

        :param bssn_coord_system: BSSN coordinate system used to construct the
                                  grid-basis Christoffel symbols.
        :param target_basis: Target spatial basis. Supported values are
                             ``"Cartesian"`` and ``"Spherical"``.
        :param enable_bssn_rfm_precompute: Whether to enable reference-metric
                                           precomputation inside the BSSN
                                           Christoffel helper.
        :return: A 4x4x4 list of SymPy expressions for the transformed Christoffel
                 symbols.
        :raises ValueError: If ``bssn_coord_system`` is not a string.
        :raises ValueError: If ``target_basis`` is unsupported.
        :raises ValueError: If ``enable_bssn_rfm_precompute`` is not a bool.
        """
        # Step 1.a: Validate the public API inputs before triggering cached lookups.
        if not isinstance(bssn_coord_system, str):
            raise ValueError(
                "bssn_coord_system must be a string, "
                f"got {type(bssn_coord_system).__name__}"
            )
        if target_basis not in ("Cartesian", "Spherical"):
            raise ValueError(
                f"Unsupported target_basis '{target_basis}'. "
                "Supported values are 'Cartesian' and 'Spherical'."
            )
        if not isinstance(enable_bssn_rfm_precompute, bool):
            raise ValueError(
                "enable_bssn_rfm_precompute must be a bool, "
                f"got {type(enable_bssn_rfm_precompute).__name__}"
            )

        # Step 1.b: Import the BSSN and reference-metric helpers lazily so
        #         analytic-only users do not pull in the full BSSN stack at import time.
        # pylint: disable-next=import-outside-toplevel,redefined-outer-name
        import nrpy.reference_metric as refmetric

        # pylint: disable-next=import-outside-toplevel
        from nrpy.equations.general_relativity.BSSN_to_g4Christoffel import (
            BSSN_to_g4Christoffel,
        )

        if bssn_coord_system not in refmetric.supported_CoordSystems:
            raise ValueError(
                f"Unsupported CoordSystem '{bssn_coord_system}'. "
                f"Supported coordinate systems: {refmetric.supported_CoordSystems}"
            )

        # Step 2: Construct the source-basis four-Christoffels Gamma^a_bc from
        #         the BSSN variables native to the requested grid basis.
        bssn_to_g4christoffel = BSSN_to_g4Christoffel(
            CoordSystem=bssn_coord_system,
            enable_rfm_precompute=enable_bssn_rfm_precompute,
        )
        grid_Gamma4UDD = bssn_to_g4christoffel.Gamma4UDD

        # Step 3: Fetch the matching ReferenceMetric object so the spatial
        #         basis transformation uses the same coordinate model as the
        #         underlying BSSN quantities.
        coord_system_key = bssn_coord_system + (
            "_rfm_precompute" if enable_bssn_rfm_precompute else ""
        )
        rfm = refmetric.reference_metric[coord_system_key]
        if len(rfm.xx) != 3:
            raise ValueError(
                "ReferenceMetric coordinates are inconsistent with 3D space."
            )

        # Step 4: Select the canonical spatial Jacobian pair
        #         J^i_j = partial x^i / partial u^j and
        #         K^j_i = partial u^j / partial x^i
        #         for the requested target basis.
        if target_basis == "Cartesian":
            spatial_J4UD = rfm.Jac_dUCart_dDrfmUD
            spatial_K4UD = rfm.Jac_dUrfm_dDCartUD
        else:
            spatial_J4UD = rfm.Jac_dUSph_dDrfmUD
            spatial_K4UD = rfm.Jac_dUrfm_dDSphUD

        # Step 5: Lift the time-independent spatial Jacobians into the 4D map
        #         x^0 = u^0, x^i = X^i(u^j), so only the purely spatial block
        #         differs from the identity.
        J4UD = ixp.zerorank2(dimension=4)
        K4UD = ixp.zerorank2(dimension=4)
        J4UD[0][0] = sp.sympify(1)
        K4UD[0][0] = sp.sympify(1)
        for spatial_target_index in range(3):
            for spatial_grid_index in range(3):
                J4UD[spatial_target_index + 1][spatial_grid_index + 1] = spatial_J4UD[
                    spatial_target_index
                ][spatial_grid_index]
                K4UD[spatial_grid_index + 1][spatial_target_index + 1] = spatial_K4UD[
                    spatial_grid_index
                ][spatial_target_index]

        # Step 6: Construct the lifted Hessian
        #         partial_k J^i_j = partial^2 X^i / partial u^k partial u^j.
        #         All entries involving time vanish because the lifted map is
        #         explicitly time independent.
        J4UD_dD = ixp.zerorank3(dimension=4)
        for spatial_target_i in range(3):
            for spatial_grid_j in range(3):
                for spatial_grid_k in range(3):
                    J4UD_dD[spatial_target_i + 1][spatial_grid_j + 1][
                        spatial_grid_k + 1
                    ] = sp.diff(
                        spatial_J4UD[spatial_target_i][spatial_grid_j],
                        rfm.xx[spatial_grid_k],
                    )

        # Step 7: Apply the connection-transformation formula to expose the
        #         BSSN Christoffels in the requested target basis.
        return GeodesicEquations._transform_christoffel_recipe_from_grid_basis_data(
            grid_Gamma4UDD=grid_Gamma4UDD,
            J4UD=J4UD,
            J4UD_dD=J4UD_dD,
            K4UD=K4UD,
        )


# Compute the metric-based Christoffel recipe once at module scope to avoid recomputation.
_GAMMA4UDD_METRIC_RECIPE = GeodesicEquations.symbolic_numerical_christoffel_recipe()
# Module-level convenience alias for callers that need only the transformed
# Christoffel recipe and do not need to instantiate GeodesicEquations.
symbolic_christoffel_recipe_from_grid_basis = (
    GeodesicEquations.symbolic_christoffel_recipe_from_grid_basis
)
symbolic_g4DD_recipe_from_bssn_grid_basis = (
    GeodesicEquations.symbolic_g4DD_recipe_from_bssn_grid_basis
)
symbolic_christoffel_recipe_from_bssn_grid_basis = (
    GeodesicEquations.symbolic_christoffel_recipe_from_bssn_grid_basis
)


class GeodesicEquations_dict(Dict[str, "GeodesicEquations"]):
    """A caching dictionary for GeodesicEquations instances."""

    def __getitem__(self, key: str) -> "GeodesicEquations":
        """
        Get or create a GeodesicEquations instance for a given configuration.

        :param key: A string key in the format "Spacetime_ParticleType".
        :return: A GeodesicEquations instance for the specified configuration.
        """
        if key not in self:
            parts = key.split("_")
            particle_type = parts[-1]
            spacetime = "_".join(parts[:-1])
            logging.getLogger(__name__).info(
                "Setting up geodesic equations for spacetime: '%s', particle_type: '%s'...",
                spacetime,
                particle_type,
            )
            self[key] = GeodesicEquations(
                spacetime=spacetime, particle_type=particle_type
            )
        return super().__getitem__(key)


Geodesic_Equations = GeodesicEquations_dict()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    import os

    import nrpy.validate_expressions.validate_expressions as ve

    # Configure logging to output to the console for direct script execution.
    logging.basicConfig(level=logging.INFO)

    # Step 1: Perform validation checks for a sample spacetime.
    print("-" * 40)
    print("Performing symmetry and sample point validation...")

    # We use Kerr-Schild (Cartesian) for basic validation in this configuration.
    geo_eqs_sample = Geodesic_Equations["KerrSchild_Cartesian_massive"]

    # Check that the metric is symmetric, g_{\mu\nu} = g_{\nu\mu}.
    for val_i in range(4):
        for val_j in range(val_i + 1, 4):
            diff = sp.simplify(
                geo_eqs_sample.g4DD[val_i][val_j] - geo_eqs_sample.g4DD[val_j][val_i]
            )
            if diff != 0:
                raise ValueError(
                    f"Metric g4DD is not symmetric: g[{val_i}][{val_j}] != g[{val_j}][{val_i}]"
                )

    # Check that the Christoffel symbols are symmetric in their lower indices
    for val_a in range(4):
        for val_i in range(4):
            for val_j in range(val_i + 1, 4):
                diff = sp.simplify(
                    geo_eqs_sample.Gamma4UDD[val_a][val_i][val_j]
                    - geo_eqs_sample.Gamma4UDD[val_a][val_j][val_i]
                )
                if diff != 0:
                    raise ValueError(
                        f"Christoffel symbols not symmetric: Gamma[{val_a}][{val_i}][{val_j}] != Gamma[{val_a}][{val_j}][{val_i}]"
                    )
    print("Symmetry checks passed.")
    print("Validation successful.")
    print("-" * 40)

    # Step 2: Validate the transformed-Christoffel recipe helpers using
    #         numeric substitution. First evaluate the fully generic
    #         grid-basis recipe on cylindrical Minkowski data and confirm the
    #         transformed Cartesian Christoffels vanish. Then validate the
    #         BSSN-specialized wrapper in an identity-map case where it must
    #         reproduce the underlying BSSN four-Christoffels exactly.
    print("Checking transformed-Christoffel recipes with numeric substitution...")
    import nrpy.equations.general_relativity.BSSN_to_g4Christoffel
    import nrpy.reference_metric as refmetric

    cylindrical_rfm = refmetric.reference_metric["Cylindrical"]
    cylindrical_xx = cylindrical_rfm.xx
    sample_rho = sp.Rational(5, 4)
    sample_phi = sp.Rational(2, 5)
    sample_z = sp.Rational(7, 10)
    cylindrical_sample_values: Dict[sp.Basic, sp.Basic] = {
        cylindrical_xx[0]: sample_rho,
        cylindrical_xx[1]: sample_phi,
        cylindrical_xx[2]: sample_z,
    }

    # Step 2.a: Supply the 4D cylindrical-to-Cartesian Jacobian and its
    #           second derivatives at a fixed nonsingular sample point.
    transformed_recipe = symbolic_christoffel_recipe_from_grid_basis()
    transformed_grid_Gamma4UDD = ixp.declarerank3(
        "grid_Gamma4UDD", dimension=4, symmetry="sym12"
    )
    transformed_J4UD = ixp.declarerank2("J4UD", dimension=4, symmetry="nosym")
    transformed_J4UD_dD = ixp.declarerank3("J4UD_dD", dimension=4, symmetry="sym12")

    transform_values: Dict[sp.Basic, sp.Basic] = {}
    for jac_alpha in range(4):
        for jac_beta in range(4):
            if jac_alpha == 0 and jac_beta == 0:
                transform_values[transformed_J4UD[jac_alpha][jac_beta]] = sp.sympify(1)
            elif jac_alpha == 0 or jac_beta == 0:
                transform_values[transformed_J4UD[jac_alpha][jac_beta]] = sp.sympify(0)
            else:
                transform_values[transformed_J4UD[jac_alpha][jac_beta]] = (
                    cylindrical_rfm.xx_to_Cart[jac_alpha - 1]
                    .diff(cylindrical_rfm.xx[jac_beta - 1])
                    .xreplace(cylindrical_sample_values)
                )

    for hess_alpha in range(4):
        for hess_nu in range(4):
            for hess_mu in range(4):
                if hess_alpha == 0 or hess_nu == 0 or hess_mu == 0:
                    transform_values[
                        transformed_J4UD_dD[hess_alpha][hess_nu][hess_mu]
                    ] = sp.sympify(0)
                else:
                    transform_values[
                        transformed_J4UD_dD[hess_alpha][hess_nu][hess_mu]
                    ] = (
                        cylindrical_rfm.xx_to_Cart[hess_alpha - 1]
                        .diff(cylindrical_rfm.xx[hess_nu - 1])
                        .diff(cylindrical_rfm.xx[hess_mu - 1])
                        .xreplace(cylindrical_sample_values)
                    )

    # Step 2.b: Evaluate the metric-based Christoffel recipe on the
    #           cylindrical Minkowski metric at the same sample point.
    cylindrical_flat_g4DD = ixp.zerorank2(dimension=4)
    cylindrical_flat_g4DD[0][0] = sp.sympify(-1)
    cylindrical_flat_g4DD[1][1] = sp.sympify(1)
    cylindrical_flat_g4DD[2][2] = cylindrical_xx[0] ** 2
    cylindrical_flat_g4DD[3][3] = sp.sympify(1)

    cylindrical_flat_g4DD_dD = ixp.zerorank3(dimension=4)
    for cyl_mu in range(4):
        for cyl_nu in range(cyl_mu, 4):
            for cyl_rho in range(4):
                cylindrical_deriv = (
                    sp.diff(
                        cylindrical_flat_g4DD[cyl_mu][cyl_nu],
                        cylindrical_xx[cyl_rho - 1],
                    ).xreplace(cylindrical_sample_values)
                    if cyl_rho > 0
                    else sp.sympify(0)
                )
                cylindrical_flat_g4DD_dD[cyl_mu][cyl_nu][cyl_rho] = cylindrical_deriv
                cylindrical_flat_g4DD_dD[cyl_nu][cyl_mu][cyl_rho] = cylindrical_deriv

    metric_recipe_g4DD_sym = ixp.declarerank2("g4DD_sym", dimension=4)
    metric_recipe_g4DD_dD_sym = ixp.declarerank3("g4DD_dD_sym", dimension=4)
    metric_recipe_values: Dict[sp.Basic, sp.Basic] = {}
    for sym_mu in range(4):
        for sym_nu in range(4):
            metric_recipe_values[metric_recipe_g4DD_sym[sym_mu][sym_nu]] = (
                cylindrical_flat_g4DD[sym_mu][sym_nu].xreplace(
                    cylindrical_sample_values
                )
            )
            for sym_rho in range(4):
                metric_recipe_values[
                    metric_recipe_g4DD_dD_sym[sym_mu][sym_nu][sym_rho]
                ] = cylindrical_flat_g4DD_dD[sym_mu][sym_nu][sym_rho]

    cylindrical_flat_gamma = GeodesicEquations.symbolic_numerical_christoffel_recipe()
    for gamma_alpha in range(4):
        for gamma_mu in range(4):
            for gamma_nu in range(4):
                cylindrical_flat_gamma[gamma_alpha][gamma_mu][gamma_nu] = (
                    cylindrical_flat_gamma[gamma_alpha][gamma_mu][gamma_nu].xreplace(
                        metric_recipe_values
                    )
                )

    # Step 2.c: Feed the numerically evaluated cylindrical Christoffels into
    #           the coordinate-transformation recipe.
    for grid_alpha_val in range(4):
        for grid_mu_val in range(4):
            for grid_nu_val in range(4):
                transform_values[
                    transformed_grid_Gamma4UDD[grid_alpha_val][grid_mu_val][grid_nu_val]
                ] = cylindrical_flat_gamma[grid_alpha_val][grid_mu_val][grid_nu_val]

    # Step 2.d: Evaluate the transformed recipe numerically. Because the
    #           target coordinates are Cartesian flat spacetime, every
    #           transformed Christoffel component should vanish.
    for check_alpha in range(4):
        for check_mu in range(4):
            for check_nu in range(check_mu, 4):
                transformed_diff = transformed_recipe[check_alpha][check_mu][
                    check_nu
                ].xreplace(transform_values)
                if not ve.check_zero(
                    transformed_diff,
                    fixed_mpfs_for_free_symbols=True,
                    verbose=False,
                ):
                    raise ValueError(
                        "Nontrivial Christoffel transform check failed for "
                        f"Gamma[{check_alpha}][{check_mu}][{check_nu}]"
                    )
    print("Generic transformed-Christoffel recipe check passed.")

    # Step 2.e: Validate the BSSN-specialized wrapper in the Cartesian
    #           identity-map case, where the transformed Christoffels must
    #           reproduce the source-basis BSSN Christoffels exactly.
    print(
        "Checking BSSN-specialized transformed-Christoffel recipe in the "
        "Cartesian identity-map case..."
    )
    bssn_identity_helper = (
        nrpy.equations.general_relativity.BSSN_to_g4Christoffel.BSSN_to_g4Christoffel(
            CoordSystem="Cartesian"
        )
    )
    bssn_identity_recipe = symbolic_christoffel_recipe_from_bssn_grid_basis(
        bssn_coord_system="Cartesian", target_basis="Cartesian"
    )
    for check_alpha in range(4):
        for check_mu in range(4):
            for check_nu in range(check_mu, 4):
                identity_diff = (
                    bssn_identity_recipe[check_alpha][check_mu][check_nu]
                    - bssn_identity_helper.Gamma4UDD[check_alpha][check_mu][check_nu]
                )
                if not ve.check_zero(
                    identity_diff,
                    fixed_mpfs_for_free_symbols=True,
                    verbose=False,
                ):
                    raise ValueError(
                        "BSSN transformed-Christoffel identity-map check failed for "
                        f"Gamma[{check_alpha}][{check_mu}][{check_nu}]"
                    )
    print("BSSN-specialized transformed-Christoffel recipe check passed.")

    # Step 2.f: Validate the BSSN-specialized transformed-four-metric recipe
    #           in the Cartesian identity-map case, where the transformed
    #           four-metric must reproduce the source-basis BSSN four-metric.
    print(
        "Checking BSSN-specialized transformed-four-metric recipe in the "
        "Cartesian identity-map case..."
    )
    bssn_identity_g4_recipe = symbolic_g4DD_recipe_from_bssn_grid_basis(
        bssn_coord_system="Cartesian", target_basis="Cartesian"
    )
    for check_mu in range(4):
        for check_nu in range(check_mu, 4):
            identity_g4_diff = (
                bssn_identity_g4_recipe[check_mu][check_nu]
                - bssn_identity_helper.g4DD[check_mu][check_nu]
            )
            if not ve.check_zero(
                identity_g4_diff,
                fixed_mpfs_for_free_symbols=True,
                verbose=False,
            ):
                raise ValueError(
                    "BSSN transformed-four-metric identity-map check failed for "
                    f"g4DD[{check_mu}][{check_nu}]"
                )
    print("BSSN-specialized transformed-four-metric identity-map check passed.")

    # Step 2.g: Validate the spherical-to-Cartesian transformed-four-metric
    #           recipe by plugging in flat-space spherical BSSN data at a
    #           fixed nonsingular sample point. The transformed metric must
    #           reduce numerically to Cartesian Minkowski spacetime.
    print(
        "Checking BSSN-specialized transformed-four-metric recipe in the "
        "Spherical-to-Cartesian case..."
    )
    spherical_recipe_g4DD = symbolic_g4DD_recipe_from_bssn_grid_basis(
        bssn_coord_system="Spherical", target_basis="Cartesian"
    )
    spherical_symbol_by_name: Dict[str, sp.Basic] = {}
    for metric_row in spherical_recipe_g4DD:
        for metric_component in metric_row:
            for free_symbol in metric_component.free_symbols:
                spherical_symbol_by_name[str(free_symbol)] = free_symbol
    spherical_numeric_values: Dict[sp.Basic, sp.Basic] = {
        spherical_symbol_by_name["alpha"]: sp.sympify(1),
        spherical_symbol_by_name["cf"]: sp.sympify(1),
        spherical_symbol_by_name["hDD00"]: sp.sympify(0),
        spherical_symbol_by_name["hDD01"]: sp.sympify(0),
        spherical_symbol_by_name["hDD02"]: sp.sympify(0),
        spherical_symbol_by_name["hDD11"]: sp.sympify(0),
        spherical_symbol_by_name["hDD12"]: sp.sympify(0),
        spherical_symbol_by_name["hDD22"]: sp.sympify(0),
        spherical_symbol_by_name["vetU0"]: sp.sympify(0),
        spherical_symbol_by_name["vetU1"]: sp.sympify(0),
        spherical_symbol_by_name["vetU2"]: sp.sympify(0),
        spherical_symbol_by_name["xx0"]: sp.Rational(5, 4),
        spherical_symbol_by_name["xx1"]: sp.Rational(2, 5),
    }
    spherical_expected_g4DD = ixp.zerorank2(dimension=4)
    spherical_expected_g4DD[0][0] = sp.sympify(-1)
    spherical_expected_g4DD[1][1] = sp.sympify(1)
    spherical_expected_g4DD[2][2] = sp.sympify(1)
    spherical_expected_g4DD[3][3] = sp.sympify(1)
    for check_mu in range(4):
        for check_nu in range(check_mu, 4):
            spherical_g4_diff = sp.simplify(
                spherical_recipe_g4DD[check_mu][check_nu].xreplace(
                    spherical_numeric_values
                )
                - spherical_expected_g4DD[check_mu][check_nu]
            )
            if spherical_g4_diff != sp.sympify(0):
                raise ValueError(
                    "BSSN transformed-four-metric spherical-to-Cartesian check "
                    f"failed for g4DD[{check_mu}][{check_nu}]"
                )
    print("BSSN-specialized transformed-four-metric spherical check passed.")

    # Step 3: Generate trusted results for all configurations.
    # This loop ensures that the __init__ logic (including the solver) works for all spacetimes.
    for config_key in ["KerrSchild_Cartesian_massive", "KerrSchild_Cartesian_photon"]:
        print(f"Processing configuration: {config_key}...")
        geodesic_eqs = Geodesic_Equations[config_key]
        results_dict = ve.process_dictionary_of_expressions(
            geodesic_eqs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{config_key}",
            results_dict,
        )

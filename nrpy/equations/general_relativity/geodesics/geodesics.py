"""
Construct symbolic expressions for geodesic equations in various spacetime backgrounds.

This module can provide symbolic Christoffel symbols and geodesic equations of motion for numerical spacetime data.
It also can take a pre-defined analytic spacetime metric and compute the geometric quantities necessary for geodesic
integration, namely the Christoffel symbols and the geodesic equations of motion.
It supports both massive and photon particles.

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
import doctest
import logging
import os
import sys
from typing import Dict, List, Optional, cast

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.indexedexp as ixp
import nrpy.validate_expressions.validate_expressions as ve
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
    Gamma4UDD_from_generic_metric: List[List[List[sp.Expr]]]

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

        # Step 4: Reuse the precomputed symbolic recipe for Christoffel symbols from a generic metric.
        self.Gamma4UDD_from_generic_metric = _GAMMA4UDD_GENERIC_RECIPE

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

        Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
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

        g4DD = ixp.declarerank2("metric->g4DD", sym="sym01", dimension=4)

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

    def geodesic_eom_rhs_photon(self) -> List[sp.Expr]:
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

        g4DD = ixp.declarerank2("metric->g4DD", sym="sym01", dimension=4)

        # We need the inverse metric g^UU to compute g^00 for the lapse function.
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)

        Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
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
        # dL/dlambda = p^0 / sqrt(-g^00)
        # We use g4UU (inverse metric) to match the definition of lapse alpha.
        path_len_rhs = [pU[0] / sp.sqrt(-g4UU[0][0])]

        return pos_rhs + mom_rhs + path_len_rhs

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

        g4DD = ixp.declarerank2("metric->g4DD", sym="sym01", dimension=4)

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
        g4DD = ixp.declarerank2("metric->g4DD", sym="sym01", dimension=4)

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


# Compute the generic Christoffel recipe once at module scope to avoid recomputation.
_GAMMA4UDD_GENERIC_RECIPE = GeodesicEquations.symbolic_numerical_christoffel_recipe()


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

    # Step 2: Run doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

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

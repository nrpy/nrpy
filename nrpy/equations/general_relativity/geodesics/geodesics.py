"""
Construct symbolic expressions for geodesic equations in various spacetime backgrounds.

This module acts as a "metric processor." It takes a pre-defined analytic
spacetime metric and computes the geometric quantities necessary for geodesic
integration, namely the Christoffel symbols and the geodesic equations of motion.
It supports both massive and massless (photon) particles.
It adheres to the nrpy "gold standard" for equations modules.

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
import doctest
import logging
import os
import sys
from typing import Dict, List

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
    ut_from_vel_expr: sp.Expr
    p0_expr: sp.Expr
    Gamma4UDD_from_generic_metric: List[List[List[sp.Expr]]]

    def __init__(self, spacetime: str, particle_type: str = "massive") -> None:
        """
        Initialize and compute all geodesic equation quantities.

        This constructor orchestrates the acquisition of the metric, the computation
        of its derivatives and Christoffel symbols, and the final generation of
        the geodesic equations of motion and Hamiltonian constraint for the
        specified particle type.

        :param spacetime: The spacetime to use (e.g., "KerrSchild_Cartesian").
        :param particle_type: The type of particle, either "massive" or "massless".
        :raises ValueError: If the particle type is not supported.
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Step 1: Acquire the metric and its associated coordinates.
        metric = Analytic_Spacetimes[spacetime]
        self.g4DD = metric.g4DD
        self.xx = metric.xx

        # Step 2: Compute geometric quantities derived from the metric.
        self.g4DD_dD = self._derivative_g4DD()
        self.Gamma4UDD = self._four_connections()

        # Step 3: Generate the right-hand-side of the geodesic ODEs and the
        #         Hamiltonian constraint equation based on the particle type.
        if self.particle_type == "massive":
            self.geodesic_rhs = self._geodesic_eom_rhs_massive()
            self.ut_from_vel_expr = self._hamiltonian_constraint_massive()
        elif self.particle_type == "massless":
            self.geodesic_rhs = self._geodesic_eom_rhs_massless()
            self.p0_expr = self._hamiltonian_constraint_massless()
        else:
            raise ValueError(f"Particle type '{self.particle_type}' is not supported.")

        # Step 4: Generate the symbolic recipe for Christoffel symbols from a generic metric.
        self.Gamma4UDD_from_generic_metric = (
            self._symbolic_numerical_christoffel_recipe()
        )

    def _derivative_g4DD(self) -> List[List[List[sp.Expr]]]:
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

    def _four_connections(self) -> List[List[List[sp.Expr]]]:
        r"""
        Compute and simplify Christoffel symbols of the second kind, \Gamma^\alpha_{\mu\nu}.

        This function calculates the Christoffel symbols using the formula:
        \Gamma^\alpha_{\mu\nu} = (1/2) g^{\alpha\beta} (g_{\beta\mu,\nu} + g_{\beta\nu,\mu} - g_{\mu\nu,\beta})
        It exploits the symmetry \Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu} by computing
        only the upper-triangular components (where \nu >= \mu) and copying them.

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
                    # Use sp.cancel for faster simplification of rational functions
                    # simplified_term = sp.cancel(term)
                    simplified_term = term
                    Gamma4UDD[alpha][mu][nu] = simplified_term
                    Gamma4UDD[alpha][nu][
                        mu
                    ] = simplified_term  # Set the symmetric component
        return Gamma4UDD

    def _geodesic_eom_rhs_massive(self) -> List[sp.Expr]:
        r"""
        Generate the symbolic right-hand-side for the 8 massive geodesic ODEs.

        The equations of motion are: d(u^\alpha)/d(\tau) = -\Gamma^\alpha_{\mu\nu} u^\mu u^\nu
        This function exploits the symmetry of the Christoffel symbols in their
        lower indices to optimize the summation, doubling off-diagonal terms.

        :return: A list of 8 SymPy expressions for the RHS of the ODEs.
        """
        uU = ixp.declarerank1("uU", dimension=4)
        pos_rhs = [uU[0], uU[1], uU[2], uU[3]]

        vel_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            # Exploit symmetry by summing over the upper triangle of (mu, nu)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = self.Gamma4UDD[alpha][mu][nu] * uU[mu] * uU[nu]
                    if mu != nu:
                        term *= 2  # Double the off-diagonal terms
                    sum_term += term
            vel_rhs[alpha] = -sum_term

        return pos_rhs + vel_rhs

    def _hamiltonian_constraint_massive(self) -> sp.Expr:
        r"""
        Symbolically derive u^0 from the Hamiltonian constraint for a massive particle.

        This function solves the 4-velocity normalization condition:
        g_{\mu\nu} u^\mu u^\nu = -1
        It returns the negative root solution, conventional for reverse ray tracing.

        :return: A SymPy expression for the negative root of u^0.
        """
        uU = ixp.declarerank1("uU", dimension=4)
        g4DD = self.g4DD

        hamiltonian = sp.sympify(1)
        # Exploit symmetry of g_{\mu\nu} to optimize summation
        for mu in range(4):
            for nu in range(mu, 4):
                term = g4DD[mu][nu] * uU[mu] * uU[nu]
                if mu != nu:
                    term *= 2
                hamiltonian += term

        solutions = sp.solve(hamiltonian, uU[0])
        return solutions[0] if len(solutions) > 1 else solutions[0]

    def _geodesic_eom_rhs_massless(self) -> List[sp.Expr]:
        r"""
        Generate the symbolic right-hand-side for the 9 massless geodesic ODEs.

        The equations of motion are: d(p^\alpha)/d(\kappa) = -\Gamma^\alpha_{\mu\nu} p^\mu p^\nu
        This function exploits the symmetry of the Christoffel symbols and the
        metric to optimize the summations.

        :return: A list of 9 SymPy expressions for the RHS of the ODEs.
        """
        pU = ixp.declarerank1("pU", dimension=4)
        pos_rhs = [pU[0], pU[1], pU[2], pU[3]]

        mom_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            # Exploit symmetry by summing over the upper triangle of (mu, nu)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = self.Gamma4UDD[alpha][mu][nu] * pU[mu] * pU[nu]
                    if mu != nu:
                        term *= 2
                    sum_term += term
            mom_rhs[alpha] = -sum_term

        path_len_sum = sp.sympify(0)
        # Exploit symmetry of g_{ij} to optimize summation
        for i in range(1, 4):
            for j in range(i, 4):
                term = self.g4DD[i][j] * pU[i] * pU[j]
                if i != j:
                    term *= 2
                path_len_sum += term
        path_len_rhs = [sp.sqrt(path_len_sum)]

        return pos_rhs + mom_rhs + path_len_rhs

    def _hamiltonian_constraint_massless(self) -> sp.Expr:
        r"""
        Symbolically derive p^0 from the Hamiltonian constraint for a massless particle.

        This function solves the null geodesic condition: g_{\mu\nu} p^\mu p^\nu = 0.

        :return: A SymPy expression for p^0.
        """
        pU = ixp.declarerank1("pU", dimension=4)
        g4DD = self.g4DD

        hamiltonian = sp.sympify(0)
        # Exploit symmetry of g_{\mu\nu} to optimize summation
        for mu in range(4):
            for nu in range(mu, 4):
                term = g4DD[mu][nu] * pU[mu] * pU[nu]
                if mu != nu:
                    term *= 2
                hamiltonian += term

        solutions = sp.solve(hamiltonian, pU[0])
        return solutions[1] if len(solutions) > 1 else solutions[0]

    @staticmethod
    def _symbolic_numerical_christoffel_recipe() -> List[List[List[sp.Expr]]]:
        r"""
        Generate a symbolic recipe for Christoffel symbols from generic metric data.

        This function exploits the symmetry \Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu}
        by computing only the upper-triangular components (where \nu >= \mu).

        :return: A list of lists of lists of SymPy expressions representing the recipe.
        """
        g4DD_sym = ixp.declarerank2("g4DD_sym", dimension=4)
        g4DD_dD_sym = ixp.declarerank3("g4DD_dD_sym", dimension=4)
        g4UU_sym, _ = ixp.symm_matrix_inverter4x4(g4DD_sym)  # type: ignore

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
    # Initialize a sample spacetime to perform checks on its symbolic expressions.
    geo_eqs_sample = Geodesic_Equations["Schwarzschild_Cartesian_Isotropic_massive"]

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

    # Check that the Christoffel symbols are symmetric in their lower indices, \Gamma^\alpha_{\mu\nu} = \Gamma^\alpha_{\nu\mu}.
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
    print("Sample point checks passed (no sample points defined for Isotropic yet).")
    print("Validation successful.")
    print("-" * 40)

    # Step 2: Run doctests and generate trusted results for all configurations.
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for config_key in [
        "KerrSchild_Cartesian_massive",
        "Schwarzschild_Cartesian_Isotropic_massive",
        "KerrSchild_Cartesian_massless",
        "Schwarzschild_Cartesian_Isotropic_massless",
    ]:
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

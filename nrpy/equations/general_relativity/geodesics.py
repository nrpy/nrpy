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
import os
import sys
from typing import Dict, List, Tuple

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy+ core modules
import nrpy.indexedexp as ixp
import nrpy.validate_expressions.validate_expressions as ve

# Step 0.d: Import the analytic spacetimes module to get metric definitions.
from nrpy.equations.general_relativity import analytic_spacetimes as anasp


class GeodesicEquations:
    """
    Generate and store symbolic expressions for geodesic motion from a given metric.
    """

    def __init__(self, spacetime: str, particle_type: str = "massive") -> None:
        """
        Initialize and compute core geodesic equations from a given metric.

        Args:
            spacetime: The spacetime to use (e.g., "KerrSchild").
            particle_type: The type of particle ("massive" or "massless").
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Step 1: Get the symbolic metric tensor for the specified spacetime.
        metric = anasp.Analytic_Spacetimes[spacetime]
        self.g4DD = metric.g4DD
        self.xx = metric.xx

        # Step 2: Compute metric derivatives and Christoffel symbols from the provided metric.
        self.g4DD_dD = self._derivative_g4DD(self.g4DD, self.xx)
        self.Gamma4UDD = self._four_connections(self.g4DD, self.g4DD_dD)

        # Step 3: Generate geodesic equations of motion based on particle type.
        if self.particle_type == "massive":
            self.geodesic_rhs = self._geodesic_eom_rhs_massive()
            self.ut_from_vel_expr = self._ut_massive()
        elif self.particle_type == "massless":
            self.geodesic_rhs = self._geodesic_eom_rhs_massless()
            self.p0_expr = self._p0_massless()
            self.Gamma4UDD_num_recipe = self._symbolic_numerical_christoffel_recipe()
        else:
            raise ValueError(f"Particle type '{self.particle_type}' is not supported.")

    @staticmethod
    def _derivative_g4DD(
        g4DD: List[List[sp.Expr]], xx: List[sp.Symbol]
    ) -> List[List[List[sp.Expr]]]:
        """Compute symbolic first derivatives of the metric tensor, g_{munu,alpha}."""
        g4DD_dD = ixp.zerorank3(dimension=4)
        for mu in range(4):
            for nu in range(4):
                for alpha in range(4):
                    g4DD_dD[mu][nu][alpha] = sp.diff(g4DD[mu][nu], xx[alpha])
        return g4DD_dD

    @staticmethod
    def _four_connections(
        g4DD: List[List[sp.Expr]], g4DD_dD: List[List[List[sp.Expr]]]
    ) -> List[List[List[sp.Expr]]]:
        """Compute and simplify Christoffel symbols of the second kind, Gamma^delta_{mu,nu}."""
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)
        Gamma4UDD = ixp.zerorank3(dimension=4)
        for delta in range(4):
            for mu in range(4):
                for nu in range(4):
                    term = sp.sympify(0)
                    for alpha in range(4):
                        term += sp.Rational(1, 2) * g4UU[delta][alpha] * (
                            g4DD_dD[nu][alpha][mu]
                            + g4DD_dD[mu][alpha][nu]
                            - g4DD_dD[mu][nu][alpha]
                        )
                    Gamma4UDD[delta][mu][nu] = sp.trigsimp(term)
        return Gamma4UDD

    @staticmethod
    def _geodesic_eom_rhs_massive() -> List[sp.Expr]:
        """Generate symbolic RHS for the 8 massive geodesic ODEs."""
        ut, ux, uy, uz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pos_rhs = [ut, ux, uy, uz]
        Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
        uU = [ut, ux, uy, uz]
        vel_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * uU[mu] * uU[nu]
                    if mu != nu:
                        term *= 2
                    sum_term += term
            vel_rhs[alpha] = -sum_term
        return pos_rhs + vel_rhs

    @staticmethod
    def _ut_massive() -> sp.Expr:
        """Symbolically derive u^t for a massive particle from g_munu u^mu u^nu = -1."""
        u0, u1, u2, u3 = sp.symbols("u0 u1 u2 u3", real=True)
        uU = [u0, u1, u2, u3]
        g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)
        sum_g0i_ui = sp.sympify(0)
        for i in range(1, 4):
            sum_g0i_ui += g4DD[0][i] * uU[i]
        sum_gij_ui_uj = sp.sympify(0)
        for i in range(1, 4):
            for j in range(i, 4):
                term = g4DD[i][j] * uU[i] * uU[j]
                if i != j:
                    term *= 2
                sum_gij_ui_uj += term
        discriminant = (2 * sum_g0i_ui) ** 2 - 4 * g4DD[0][0] * (sum_gij_ui_uj + 1)
        answer = (-2 * sum_g0i_ui - sp.sqrt(discriminant)) / (2 * g4DD[0][0])
        return answer

    @staticmethod
    def _geodesic_eom_rhs_massless() -> List[sp.Expr]:
        """Generate symbolic RHS for the 9 massless geodesic ODEs."""
        # Position and Momentum ODEs
        pt, pr, pth, pph = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pos_rhs = [pt, pr, pth, pph]
        Gamma4UDD = ixp.declarerank3("conn->Gamma4UDD", dimension=4, sym="sym12")
        pU = [pt, pr, pth, pph]
        mom_rhs = ixp.zerorank1(dimension=4)
        for alpha in range(4):
            sum_term = sp.sympify(0)
            for mu in range(4):
                for nu in range(mu, 4):
                    term = Gamma4UDD[alpha][mu][nu] * pU[mu] * pU[nu]
                    if mu != nu:
                        term *= 2
                    sum_term += term
            mom_rhs[alpha] = -sum_term

        # Path Length ODE: dL/dκ = sqrt(γ_ij p^i p^j)
        g4DD = ixp.declarerank2("metric->g", dimension=4, sym="sym01")
        path_len_sum = sp.sympify(0)
        for i in range(1, 4):
            for j in range(i, 4):
                term = g4DD[i][j] * pU[i] * pU[j]
                if i != j:
                    term *= 2
                path_len_sum += term
        path_len_rhs = [sp.sqrt(path_len_sum)]

        return pos_rhs + mom_rhs + path_len_rhs

    @staticmethod
    def _p0_massless() -> sp.Expr:
        """Symbolically derive p^0 for a massless particle from g_munu p^mu p^nu = 0."""
        p0, p1, p2, p3 = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pU = [p0, p1, p2, p3]
        g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)
        sum_g0i_pi = sp.sympify(0)
        for i in range(1, 4):
            sum_g0i_pi += g4DD[0][i] * pU[i]
        sum_gij_pi_pj = sp.sympify(0)
        for i in range(1, 4):
            for j in range(i, 4):
                term = g4DD[i][j] * pU[i] * pU[j]
                if i != j:
                    term *= 2
                sum_gij_pi_pj += term
        discriminant = sum_g0i_pi**2 - g4DD[0][0] * sum_gij_pi_pj
        answer = (-sum_g0i_pi + sp.sqrt(discriminant)) / g4DD[0][0]
        return answer

    @staticmethod
    def _symbolic_numerical_christoffel_recipe() -> List[List[List[sp.Expr]]]:
        """Generate the symbolic recipe for Christoffels from numerical metric data."""
        g4DD = ixp.declarerank2("g4DD", symmetry="sym01", dimension=4)
        g4DDdD = ixp.zerorank3(dimension=4)
        for i in range(4):
            for j in range(i, 4):
                for k in range(4):
                    symbol_name = f"g4DDdD{i}{j}d{k}"
                    g4DDdD[i][j][k] = sp.Symbol(symbol_name)
                    if i != j:
                        g4DDdD[j][i][k] = g4DDdD[i][j][k]
        g4UU, _ = ixp.symm_matrix_inverter4x4(g4DD)
        Gamma4UDD_num_recipe = ixp.zerorank3(dimension=4)
        for alpha in range(4):
            for mu in range(4):
                for nu in range(mu, 4):
                    for delta in range(4):
                        Gamma4UDD_num_recipe[alpha][mu][nu] += sp.Rational(1, 2) * g4UU[alpha][delta] * (
                            g4DDdD[nu][delta][mu]
                            + g4DDdD[mu][delta][nu]
                            - g4DDdD[mu][nu][delta]
                        )
        return Gamma4UDD_num_recipe


class GeodesicEquations_dict(Dict[str, GeodesicEquations]):
    """A caching dictionary for GeodesicEquations instances."""

    def __getitem__(self, key: str) -> GeodesicEquations:
        """
        Get or create a GeodesicEquations instance for a given configuration.

        Args:
            key: A string key in the format "Spacetime_ParticleType".

        Returns:
            A GeodesicEquations instance for the specified configuration.
        """
        if key not in self:
            parts = key.split("_")
            if len(parts) < 2:
                raise ValueError("Invalid key format. Expected 'Spacetime_ParticleType'.")
            particle_type = parts[-1]
            spacetime = "_".join(parts[:-1])

            print(
                f"Setting up GeodesicEquations for spacetime='{spacetime}', particle='{particle_type}'..."
            )
            self[key] = GeodesicEquations(
                spacetime=spacetime, particle_type=particle_type
            )
        return super().__getitem__(key)


Geodesic_Equations = GeodesicEquations_dict()


if __name__ == "__main__":
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Update validation loop to include massless configurations.
    for config in [
        "KerrSchild_massive",
        "Schwarzschild_Cartesian_massive",
        "KerrSchild_massless",
        "Schwarzschild_Cartesian_massless",
    ]:
        geo_eqs = Geodesic_Equations[config]
        results_dict = ve.process_dictionary_of_expressions(
            geo_eqs.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{config}",
            results_dict,
        )

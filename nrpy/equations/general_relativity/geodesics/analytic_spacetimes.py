"""
Construct symbolic expressions for various analytic spacetime metrics.

This module provides a class-based structure for generating the symbolic
metric tensor for several common analytic solutions to Einstein's equations,
such as Kerr-Schild and Schwarzschild. It is designed to integrate with
nrpy's CodeParameter system.

Author: Dalton J. Moone
"""

# Step 0.a: Import standard Python modules
from typing import Dict, List, Tuple

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import NRPy core modules
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.validate_expressions.validate_expressions as ve

# Step 0.d: Define global CodeParameters for physical parameters.
# These are special sympy symbols that the C code generator recognizes.
M_scale = par.register_CodeParameter("REAL", __name__, "M_scale", 1.0)
a_spin = par.register_CodeParameter("REAL", __name__, "a_spin", 0.0)


class AnalyticSpacetimes:
    """
    Generate and store symbolic expressions for analytic spacetime metrics.

    This class is instantiated with a specific spacetime name. It then calls
    the appropriate recipe to generate the 4-metric g_munu and the underlying
    coordinate system symbols, storing them as instance attributes.
    """

    # mypy --strict requires class attributes to be declared.
    spacetime_name: str
    g4DD: List[List[sp.Expr]]
    xx: List[sp.Symbol]

    def __init__(self, spacetime_name: str) -> None:
        """
        Initialize and generate the symbolic metric for a given spacetime.

        :param spacetime_name: The name of the spacetime to generate
                               (e.g., "KerrSchild", "Schwarzschild_Cartesian").
        :raises ValueError: If the requested spacetime is not supported.
        """
        self.spacetime_name = spacetime_name

        if self.spacetime_name == "KerrSchild":
            self.g4DD, self.xx = self._define_kerr_metric_Cartesian_Kerr_Schild()
        elif self.spacetime_name == "Schwarzschild_Cartesian":
            self.g4DD, self.xx = self._define_schwarzschild_metric_cartesian()
        else:
            raise ValueError(f"Spacetime '{self.spacetime_name}' is not supported.")

    @staticmethod
    def _define_kerr_metric_Cartesian_Kerr_Schild() -> (
        Tuple[List[List[sp.Expr]], List[sp.Symbol]]
    ):
        """
        Define the Kerr metric in Cartesian Kerr-Schild coordinates.

        The metric is constructed as g_munu = eta_munu + 2H * l_mu * l_nu.

        :return: A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
                 and xx is the list of symbolic coordinate variables (x, y, z).
        """
        # Step 1: Define generic symbolic coordinates.
        x, y, z = sp.symbols("x y z", real=True)
        xx = [x, y, z]

        # Step 2: Define intermediate geometric quantities.
        r2 = x**2 + y**2 + z**2
        r = sp.sqrt(r2)

        # Step 3: Define the Kerr-Schild null vector l_mu.
        # Formula: l_mu = (1, (r*x + a*y)/(r^2+a^2), (r*y - a*x)/(r^2+a^2), z/r)
        l_down = ixp.zerorank1(dimension=4)
        l_down[0] = sp.sympify(1)
        l_down[1] = (r * x + a_spin * y) / (r2 + a_spin**2)
        l_down[2] = (r * y - a_spin * x) / (r2 + a_spin**2)
        l_down[3] = z / r

        # Step 4: Define the scalar function H.
        # Formula: H = (M*r^3) / (r^4 + a^2*z^2)
        H = (M_scale * r**3) / (r**4 + a_spin**2 * z**2)

        # Step 5: Construct the Kerr-Schild metric g_munu = eta_munu + 2H * l_mu * l_nu.
        eta4DD = ixp.zerorank2(dimension=4)
        eta4DD[0][0] = -1
        eta4DD[1][1] = eta4DD[2][2] = eta4DD[3][3] = 1
        g4DD = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(4):
                g4DD[mu][nu] = eta4DD[mu][nu] + 2 * H * l_down[mu] * l_down[nu]

        return g4DD, xx

    @staticmethod
    def _define_schwarzschild_metric_cartesian() -> (
        Tuple[List[List[sp.Expr]], List[sp.Symbol]]
    ):
        """
        Define the Schwarzschild metric in standard Cartesian coordinates.

        This recipe uses the global CodeParameter symbol for M.

        :return: A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
                 and xx is the list of symbolic coordinate variables.
        """
        # Step 1: Define generic symbolic coordinates and radial distance.
        x, y, z = sp.symbols("x y z", real=True)
        xx = [x, y, z]
        r = sp.sqrt(x**2 + y**2 + z**2)

        g4DD = ixp.zerorank2(dimension=4)
        # Step 2: Set the time-time component.
        # Formula: g_tt = -(1 - 2M/r).
        g4DD[0][0] = -(1 - 2 * M_scale / r)

        # Step 3: Set the space-space components.
        # Formula: g_ij = delta_ij + (2M/r) * (x_i * x_j / r^2).
        for i in range(3):
            for j in range(3):
                delta_ij = sp.sympify(1) if i == j else sp.sympify(0)
                g4DD[i + 1][j + 1] = delta_ij + (2 * M_scale / r) * (
                    xx[i] * xx[j] / (r**2)
                )
        return g4DD, xx


class AnalyticSpacetimes_dict(Dict[str, "AnalyticSpacetimes"]):
    """A caching dictionary for AnalyticSpacetimes instances."""

    def __getitem__(self, key: str) -> "AnalyticSpacetimes":
        """
        Get or create an AnalyticSpacetimes instance for a given configuration.

        :param key: A string key identifying the spacetime (e.g., "KerrSchild").
        :return: An AnalyticSpacetimes instance for the specified configuration.
        """
        if key not in self:
            print(f"Setting up analytic spacetime: '{key}'...")
            self[key] = AnalyticSpacetimes(spacetime_name=key)
        return super().__getitem__(key)


Analytic_Spacetimes = AnalyticSpacetimes_dict()


if __name__ == "__main__":
    # Move imports here to avoid pylint unused-import errors.
    import doctest
    import os
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Use a distinct loop variable name to avoid pylint redefined-outer-name warnings.
    for spacetime_name_str in ["KerrSchild", "Schwarzschild_Cartesian"]:
        spacetimes = Analytic_Spacetimes[spacetime_name_str]
        results_dict = ve.process_dictionary_of_expressions(
            spacetimes.__dict__, fixed_mpfs_for_free_symbols=True
        )
        # Break long line to satisfy pylint line-too-long.
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{spacetime_name_str}",
            results_dict,
        )

"""
Construct symbolic expressions for analytic spacetime metrics.

This module provides a class-based structure for generating the symbolic
metric tensor for supported analytic or analytic-like spacetime recipes.
It is designed to integrate with nrpy's CodeParameter system.

For Kerr, ``a_spin`` is the dimensional Kerr parameter ``a = J/M`` in
geometric units, with the same length units as ``M_scale``. The dimensionless
spin is ``a_spin / M_scale``. Analytic geodesic examples use ``M_scale = 1``.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com

"""

# Step 0.a: Import standard Python modules
import logging
from typing import Dict, List, Tuple

# Step 0.b: Import third-party modules
import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.validate_expressions.validate_expressions as ve


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
    inertial_coordinates: List[sp.Expr]
    inertial_xx: List[sp.Symbol]
    computational_coordinates_from_inertial: List[sp.Expr]
    flat_a: sp.Expr
    flat_b: sp.Expr
    flat_psi: sp.Expr
    flat_kappa: sp.Expr

    def __init__(self, spacetime_name: str) -> None:
        """
        Initialize and generate the symbolic metric for a given spacetime.

        :param spacetime_name: The name of the spacetime to generate
                               (e.g., "KerrSchild_Cartesian").
        :raises ValueError: If the requested spacetime is not supported.
        """
        self.spacetime_name = spacetime_name

        if self.spacetime_name == "KerrSchild_Cartesian":
            self.g4DD, self.xx = self._define_kerr_metric_Cartesian_Kerr_Schild()
        elif self.spacetime_name == "BrillLindquist_InitialData_Static_Cartesian":
            # Static continuation of coincident Brill-Lindquist spatial data:
            # gamma_ij = psi^4 delta_ij, alpha = psi^-2, and zero shift.
            t, x, y, z = sp.symbols("t x y z", real=True)
            xx = [t, x, y, z]

            M_total = par.register_CodeParameter(
                "REAL", __name__, "M_total", 1.0, commondata=True
            )
            isotropic_radius = sp.sqrt(x**2 + y**2 + z**2)
            psi = sp.sympify(1) + M_total / (2 * isotropic_radius)

            g4DD = ixp.zerorank2(dimension=4)
            g4DD[0][0] = -(psi ** (-4))
            for spatial_index in range(1, 4):
                g4DD[spatial_index][spatial_index] = psi**4
            self.g4DD, self.xx = g4DD, xx
        elif self.spacetime_name == "Minkowski_TimeDependentAxisymmetric_Cartesian":
            # Pull back the inertial flat metric through the time-dependent
            # scale, rotation, and twist map. The inverse maps photon positions.
            t, x, y, z = sp.symbols("t x y z", real=True)
            self.xx = [t, x, y, z]
            parameter_defaults = (
                ("flat_epsilon_a", 0.0),
                ("flat_omega_a", 1.0),
                ("flat_epsilon_b", 0.0),
                ("flat_omega_b", 1.0),
                ("flat_epsilon_Omega", 0.0),
                ("flat_omega_Omega", 1.0),
                ("flat_epsilon_kappa", 0.0),
                ("flat_omega_kappa", 1.0),
                ("flat_L", 1.0),
            )
            parameters = {
                name: par.register_CodeParameter(
                    "REAL", __name__, name, default, commondata=True
                )
                for name, default in parameter_defaults
            }
            a = 1 + parameters["flat_epsilon_a"] * sp.sin(
                parameters["flat_omega_a"] * t
            )
            b = 1 + parameters["flat_epsilon_b"] * sp.sin(
                parameters["flat_omega_b"] * t
            )
            Omega = parameters["flat_epsilon_Omega"] * sp.sin(
                parameters["flat_omega_Omega"] * t
            )
            kappa = (
                parameters["flat_epsilon_kappa"]
                * sp.sin(parameters["flat_omega_kappa"] * t)
                / parameters["flat_L"]
            )
            psi = Omega + kappa * z
            self.flat_a, self.flat_b = a, b
            self.flat_psi, self.flat_kappa = psi, kappa
            cos_psi, sin_psi = sp.cos(psi), sp.sin(psi)
            self.inertial_coordinates = [
                t,
                a * (x * cos_psi - y * sin_psi),
                a * (x * sin_psi + y * cos_psi),
                b * z,
            ]
            T, X, Y, Z = sp.symbols("T X Y Z", real=True)
            self.inertial_xx = [T, X, Y, Z]
            z_from_inertial = Z / b.subs(t, T)
            psi_from_inertial = Omega.subs(t, T) + kappa.subs(t, T) * z_from_inertial
            self.computational_coordinates_from_inertial = [
                T,
                (X * sp.cos(psi_from_inertial) + Y * sp.sin(psi_from_inertial))
                / a.subs(t, T),
                (-X * sp.sin(psi_from_inertial) + Y * sp.cos(psi_from_inertial))
                / a.subs(t, T),
                z_from_inertial,
            ]

            jacobian = [
                [sp.diff(inertial_coordinate, coordinate) for coordinate in self.xx]
                for inertial_coordinate in self.inertial_coordinates
            ]
            self.g4DD = ixp.zerorank2(dimension=4)
            for mu in range(4):
                for nu in range(mu, 4):
                    metric_component = -jacobian[0][mu] * jacobian[0][nu] + sum(
                        jacobian[A][mu] * jacobian[A][nu] for A in range(1, 4)
                    )
                    self.g4DD[mu][nu] = self.g4DD[nu][mu] = metric_component
        else:
            raise ValueError(f"Spacetime '{self.spacetime_name}' is not supported.")

    @staticmethod
    def _define_kerr_metric_Cartesian_Kerr_Schild() -> (
        Tuple[List[List[sp.Expr]], List[sp.Symbol]]
    ):
        """
        Define the Kerr metric in Cartesian Kerr-Schild coordinates.

        The metric is constructed as g_munu = eta_munu + 2H * l_mu * l_nu.
        This form is regular everywhere, including the horizon.

        Reference:
        Wikipedia: Kerr-Schild coordinates
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Kerr_metric&oldid=1318460406
        (See section on Kerr–Schild coordinates)

        :return: A tuple (g4DD, xx), where g4DD is the symbolic 4x4 metric tensor
                 and xx is the list of symbolic coordinate variables (t, x, y, z).
        """
        # Step 1.a: Define generic symbolic coordinates.
        t, x, y, z = sp.symbols("t x y z", real=True)
        xx = [t, x, y, z]

        # Step 1.b: Register physical parameters in geometric units (G=c=1).
        # M_scale is ADM mass; a_spin is dimensional Kerr a=J/M. The
        # dimensionless spin is a_spin / M_scale.
        M_scale = par.register_CodeParameter(
            "REAL", __name__, "M_scale", 1.0, commondata=True
        )
        a_spin = par.register_CodeParameter(
            "REAL", __name__, "a_spin", 0.0, commondata=True
        )

        # Step 2: Define intermediate geometric quantities.
        # The Kerr-Schild radius 'r' is not the Euclidean radius. It is solved
        # for implicitly from the Cartesian coordinates (x, y, z) and spin a.
        # rho2 is the squared Euclidean distance from the origin.
        rho2 = x**2 + y**2 + z**2
        a_spin2 = a_spin**2

        # This is the solution to the quartic equation for r:
        # r^4 - (rho^2 - a^2)r^2 - a^2 z^2 = 0
        r2 = sp.Rational(1, 2) * (
            rho2 - a_spin2 + sp.sqrt((rho2 - a_spin2) ** 2 + 4 * a_spin2 * z**2)
        )
        r = sp.sqrt(r2)

        # Step 3: Define the Kerr-Schild null vector l_mu.
        l_down = ixp.zerorank1(dimension=4)
        l_down[0] = sp.sympify(1)
        l_down[1] = (r * x + a_spin * y) / (r2 + a_spin2)
        l_down[2] = (r * y - a_spin * x) / (r2 + a_spin2)
        l_down[3] = z / r

        # Step 4: Define the scalar function H.
        H = (M_scale * r**3) / (r**4 + a_spin2 * z**2)

        # Step 5: Construct the Kerr-Schild metric g_munu = eta_munu + 2H * l_mu * l_nu.
        eta4DD = ixp.zerorank2(dimension=4)
        eta4DD[0][0] = sp.sympify(-1)
        eta4DD[1][1] = eta4DD[2][2] = eta4DD[3][3] = sp.sympify(1)
        g4DD = ixp.zerorank2(dimension=4)
        for mu in range(4):
            for nu in range(4):
                g4DD[mu][nu] = eta4DD[mu][nu] + 2 * H * l_down[mu] * l_down[nu]

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
            # If the spacetime is not cached, generate it and add it to the cache.
            logging.getLogger(__name__).info(
                "Setting up analytic spacetime: '%s'...", key
            )
            self[key] = AnalyticSpacetimes(spacetime_name=key)
        return super().__getitem__(key)


Analytic_Spacetimes = AnalyticSpacetimes_dict()


if __name__ == "__main__":
    import doctest
    import os
    import sys

    # Configure logging to output to the console.
    logging.basicConfig(level=logging.INFO)

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Use a distinct loop variable name to avoid pylint redefined-outer-name warnings.
    for spacetime_name_str in [
        "KerrSchild_Cartesian",
        "BrillLindquist_InitialData_Static_Cartesian",
    ]:
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

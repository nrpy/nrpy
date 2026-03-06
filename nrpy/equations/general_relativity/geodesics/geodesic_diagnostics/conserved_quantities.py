"""
Construct symbolic expressions for diagnostic quantities along geodesics.

When the spacetime admits the relevant symmetries, this includes constants
of motion such as Energy, L_z, and (in Kerr) the Carter constant. Other
returned quantities may be non-conserved diagnostics.
It supports both massive and photon particles.

Author: Dalton J. Moone

--------------------------------------------------------------------------------
List of Assumptions
--------------------------------------------------------------------------------
Physical Assumptions:
1.  **Metric Signature:** The code calculates Energy as E = -p_0. This assumes a
    metric signature of (-, +, +, +). If the signature were (+, -, -, -), Energy
    would be p_0.
2.  **Stationarity:** The calculation of Energy assumes the spacetime possesses
    a timelike Killing vector xi = d/dt. This is applied to all spacetimes.
3.  **Axial Symmetry:** The calculation of the conserved Angular Momentum component
    (L_z) assumes the spacetime possesses an axial Killing vector xi = d/dphi.
4.  **Mass Normalization (m=1):** In the massive particle case for the Carter
    Constant, the code uses the term (1 - E**2). This implies the code assumes
    particle mass m=1 (or affine parameter normalization p_mu p^mu = -1).
5.  **Photon Particles:** For the photon case, the code assumes m=0.
6.  **Geometric Units:** The code assumes geometric units where G = c = 1 (implied
    by the usage of M_scale and a_spin without constants).

Coordinate & Geometric Assumptions:
7.  **Cartesian Coordinates for Angular Momentum:** The code assumes that if a
    spacetime name ends with "Cartesian", the coordinates at indices [1, 2, 3]
    correspond to standard Cartesian [x, y, z] to satisfy the cross-product
    relation L = x X p.
8.  **Kerr-Schild Geometry:** The Carter Constant logic assumes specific
    geometric identities (e.g., implicit r definition, oblate spheroidal
    Jacobians). It assumes the standard transformation between Cartesian
    Kerr-Schild and Boyer-Lindquist coordinates.

Implementation & Logic Assumptions:
9.  **Spacetime Specificity for Q:** The Carter Constant (Q) is **only**
    calculated if the spacetime is exactly "KerrSchild_Cartesian". For all
    other spacetimes (including Numerical), Q is set to None.
10. **Spacetime Specificity for L:** Angular momentum vectors are **only**
    calculated if the spacetime name ends with "Cartesian".
11. **Numerical Fallback:** If the spacetime is set to "Numerical", the code
    assumes a generic 4D metric and generic coordinates [x0, x1, x2, x3]
    without loading analytic expressions.
12. **Contravariant Input:** The input pU is assumed to be the contravariant
    4-momentum vector p^mu.
--------------------------------------------------------------------------------
"""

# Step 0.a: Import standard Python modules
import doctest
import logging
import os
import sys
from typing import Any, Dict, List, Optional, Tuple, cast

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import nrpy core modules
import nrpy.indexedexp as ixp
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
)


class GeodesicDiagnostics:
    """
    Generate and store symbolic expressions for conserved quantities.

    This class is instantiated with a specific spacetime and particle type.
    It computes the symbolic expressions for Energy (E), Angular Momentum (L),
    and the Carter Constant (Q) using separate, modular methods.
    """

    # mypy --strict requires class attributes to be declared.
    spacetime: str
    particle_type: str
    g4DD: List[List[sp.Expr]]
    xx: List[sp.Symbol]
    E_expr: Optional[sp.Expr]
    L_exprs: List[sp.Expr]
    Q_expr: Optional[sp.Expr]

    def __init__(
        self,
        spacetime: str,
        particle_type: str = "massive",
    ) -> None:
        """
        Initialize and generate all symbolic diagnostic expressions.

        :param spacetime: The spacetime to use (e.g., "KerrSchild_Cartesian") or "Numerical".
        :param particle_type: The type of particle ("massive" or "photon").
        :raises ValueError: If the particle type is not supported or spacetime is invalid.
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Initialize optional attributes
        self.E_expr = None
        self.L_exprs = []
        self.Q_expr = None

        if self.particle_type not in ["massive", "photon"]:
            raise ValueError(f"Unsupported particle_type: {self.particle_type}")

        # Step 1: Set up coordinates and metric
        if self.spacetime == "Numerical":
            # Generate generic symbolic metric and coordinates
            self.xx = [sp.Symbol(f"x{i}", real=True) for i in range(4)]
            # Suppress mypy complaint about Sequence vs List
            self.g4DD = ixp.declarerank2("g4DD", sym="sym01", dimension=4)
        else:
            # Try to load from analytic module.
            try:
                metric_data = Analytic_Spacetimes[self.spacetime]
                self.xx = metric_data.xx
                self.g4DD = metric_data.g4DD
            except ValueError as exc:
                # AnalyticSpacetimes raises ValueError if the key is not supported internally
                raise ValueError(
                    f"Spacetime '{self.spacetime}' is not supported. "
                    "Please check analytic_spacetimes.py for available metrics "
                    "or use 'Numerical' for generic derivation."
                ) from exc

        # Define 4-momentum vector pU = [p0, p1, p2, p3]
        pU = [sp.Symbol(f"p{i}", real=True) for i in range(4)]

        # Step 2: Compute quantities using modular functions.
        # Energy-like quantity: E := -p_0. It is conserved only if the spacetime is stationary
        # (i.e., admits a timelike Killing vector aligned with the chosen time coordinate).
        self.E_expr = self.compute_energy(pU)

        # Angular Momentum is computed based on coordinate basis.
        if self.spacetime.endswith("Cartesian"):
            self.L_exprs = self.compute_angular_momentum_cartesian(pU)

        # Step 3: Compute Spacetime-Specific Quantities (Carter Constant)
        # This formula is specific to Kerr geometry in Cartesian coordinates.
        if spacetime == "KerrSchild_Cartesian":
            self.Q_expr = self.compute_carter_constant_KerrSchild_Cartesian(
                self.E_expr, self.L_exprs, pU
            )
        else:
            # For spacetimes where Q is not defined or not yet implemented
            self.Q_expr = None

    def compute_energy(self, pU: List[sp.Symbol]) -> sp.Expr:
        """
        Compute the conserved energy E = -p_0.

        This assumes the spacetime is stationary (has a timelike Killing vector d/dt).

        Reference:
        Wikipedia: Kerr Metric
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Kerr_metric&oldid=1323832196#Boyer%E2%80%93Lindquist_coordinates
        (See second equation in Section: Trajectory Equations)

        :param pU: The symbolic 4-momentum vector.
        :return: Symbolic expression for Energy.
        """
        p_0 = sp.sympify(0)
        for mu in range(4):
            # Exploit symmetry: g_{0,mu} = g_{mu,0}. Use upper triangle.
            p_0 += self.g4DD[0][mu] * pU[mu]
        return cast(sp.Expr, -p_0)

    def compute_angular_momentum_cartesian(self, pU: List[sp.Symbol]) -> List[sp.Expr]:
        """
        Compute the angular momentum vector components L_i.

        Note: In non-spherically symmetric spacetimes like Kerr, only the
        component aligned with the symmetry axis (usually L_z) is conserved.

        Note: This implementation uses the cross-product formula L = x X p.
        This assumes the input coordinates `self.xx` are Cartesian-like, where
        xx[1]=x, xx[2]=y, xx[3]=z.

        Reference:
        Wikipedia: Angular momentum
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Angular_momentum&oldid=1326178215
        (See fourth equation in Section: Orbital angular momentum in three dimensions)


        :param pU: The symbolic 4-momentum vector.
        :return: A list of symbolic expressions [L_x, L_y, L_z].
        """
        # First, compute covariant spatial momentum p_k = g_{k,mu} p^mu.
        pD = ixp.zerorank1(dimension=4)
        for k in range(1, 4):
            for mu in range(4):
                if k <= mu:
                    pD[k] += self.g4DD[k][mu] * pU[mu]
                else:
                    pD[k] += self.g4DD[mu][k] * pU[mu]

        # Map inputs to Cartesian logic for cross product
        # xx[0] is time, so x,y,z are indices 1,2,3
        x, y, z = self.xx[1], self.xx[2], self.xx[3]
        p_x, p_y, p_z = pD[1], pD[2], pD[3]

        # L = x × p  (Cartesian cross product)
        L_x = y * p_z - z * p_y
        L_y = z * p_x - x * p_z
        L_z = x * p_y - y * p_x

        return [cast(sp.Expr, L_x), cast(sp.Expr, L_y), cast(sp.Expr, L_z)]

    def compute_carter_constant_KerrSchild_Cartesian(
        self,
        E: sp.Expr,
        L_exprs: List[sp.Expr],
        pU: List[sp.Symbol],
    ) -> sp.Expr:
        """
        Compute the Carter Constant Q.

        This function handles the specific geometric identities required to
        compute Q in Kerr-Schild (Cartesian) coordinates.

        Reference:
        Wikipedia: Kerr Metric
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Kerr_metric&oldid=1323832196#Boyer%E2%80%93Lindquist_coordinates
        (See fourth equation in Section: Trajectory Equations)

        :param E: Symbolic Energy expression.
        :param L_exprs: List of Angular Momentum components [Lx, Ly, Lz].
        :param pU: Symbolic 4-momentum vector.
        :return: Symbolic expression for Q.
        """
        L_z = L_exprs[2]
        # Define a_spin locally (G=c=1)
        a_spin = sp.Symbol("a_spin", real=True)

        # Kerr spacetime logic (Cartesian Kerr-Schild)
        # xx[0] is time, so x,y,z are indices 1,2,3
        x, y, z = self.xx[1], self.xx[2], self.xx[3]

        # Re-compute p_down locally to avoid dependency on previous method return values
        p_down = ixp.zerorank1(dimension=4)
        for k in range(1, 4):
            for mu in range(4):
                if k <= mu:
                    p_down[k] += self.g4DD[k][mu] * pU[mu]
                else:
                    p_down[k] += self.g4DD[mu][k] * pU[mu]
        p_x, p_y, p_z = p_down[1], p_down[2], p_down[3]

        # Calculate r^2 implicitly for Kerr-Schild
        Sigma = x**2 + y**2 + z**2 - a_spin**2
        r_sq = (Sigma + sp.sqrt(Sigma**2 + 4 * a_spin**2 * z**2)) / 2

        # Define r and rho for p_theta calculation
        r = sp.sqrt(r_sq)
        rho_sq = x**2 + y**2
        rho = sp.sqrt(rho_sq)

        # Calculate p_theta correctly by projecting 4-momentum onto theta basis vector
        # p_theta = (x*p_x + y*p_y) * cot(theta) - r*p_z*sin(theta)
        # Using identities for Kerr-Schild coordinates:
        # cot(theta) = (z * sqrt(r^2 + a^2)) / (r * rho)
        # sin(theta) = rho / sqrt(r^2 + a^2)

        sqrt_r2_plus_a2 = sp.sqrt(r_sq + a_spin**2)

        # Factor for (x*px + y*py) term
        factor_1 = (z * sqrt_r2_plus_a2) / (r * rho)

        # Factor for p_z term
        factor_2 = (r * rho) / sqrt_r2_plus_a2

        p_theta = factor_1 * (x * p_x + y * p_y) - factor_2 * p_z
        p_theta_sq = p_theta**2

        # Geometric identity: 1/sin^2(theta) = (r^2 + a^2) / rho^2
        inv_sin_theta_sq = (r_sq + a_spin**2) / rho_sq

        if self.particle_type == "massive":
            # Q = p_theta^2 + cos^2(theta) * (a^2(1 - E^2) + L_z^2/sin^2(theta))
            # cos^2(theta) = z^2/r^2
            # Assumes mass m=1 (4-velocity normalization)
            second_term = (z**2 / r_sq) * (
                a_spin**2 * (1 - E**2) + L_z**2 * inv_sin_theta_sq
            )
        else:
            # Q = p_theta^2 + cos^2(theta) * (-a^2 E^2 + L_z^2/sin^2(theta))
            second_term = (z**2 / r_sq) * (
                -(a_spin**2) * E**2 + L_z**2 * inv_sin_theta_sq
            )

        Q_formula = p_theta_sq + second_term

        return cast(sp.Expr, Q_formula)


class GeodesicDiagnostics_dict(Dict[str, GeodesicDiagnostics]):
    """A caching dictionary for GeodesicDiagnostics instances."""

    def __getitem__(self, key: str) -> GeodesicDiagnostics:
        """
        Get or create a GeodesicDiagnostics instance for a given configuration.

        :param key: A string key in the format "Spacetime_ParticleType".
        :return: A GeodesicDiagnostics instance for the specified configuration.
        :raises ValueError: If the key format is incorrect.
        """
        if key not in self:
            parts = key.split("_")
            if len(parts) < 2:
                raise ValueError(
                    "Invalid key format. Expected 'Spacetime_ParticleType'."
                )
            particle_type = parts[-1]
            spacetime = "_".join(parts[:-1])

            logging.getLogger(__name__).info(
                "Setting up GeodesicDiagnostics for spacetime='%s', particle='%s'...",
                spacetime,
                particle_type,
            )

            # The class now handles metric lookup internally
            self[key] = GeodesicDiagnostics(
                spacetime=spacetime,
                particle_type=particle_type,
            )
        return super().__getitem__(key)


Geodesic_Diagnostics = GeodesicDiagnostics_dict()


if __name__ == "__main__":
    # Configure logging to output to the console for direct script execution.
    logging.basicConfig(level=logging.INFO)

    # Step 1: Run doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Step 2: Perform Physics Validation Checks
    print("-" * 40)
    print("Performing physics validation checks...")

    # Check 1: In Kerr (limit a->0), Q must reduce to L^2 - L_z^2.
    # Note: The general Kerr Q formula (Carter constant) reduces to L^2 - L_z^2
    # in the spherical limit.
    # We verify that Q_kerr(a=0) + L_z^2 == L^2.
    kerr_diag = Geodesic_Diagnostics["KerrSchild_Cartesian_massive"]

    L_sq_kerr = (
        kerr_diag.L_exprs[0] ** 2
        + kerr_diag.L_exprs[1] ** 2
        + kerr_diag.L_exprs[2] ** 2
    )
    L_z_sq_kerr = kerr_diag.L_exprs[2] ** 2

    # Extract the main formula from the diagnostic expression.
    assert kerr_diag.Q_expr is not None
    Q_expr = kerr_diag.Q_expr

    if isinstance(Q_expr, sp.Piecewise):
        # FIX: Explicitly cast .args to a list of tuples (Expression, Condition)
        # generic SymPy .args is Tuple[Basic, ...], which isn't unpackable.
        piecewise_args = cast(Tuple[Tuple[sp.Expr, Any], ...], Q_expr.args)

        Q_kerr_formula = next(
            expr
            for (expr, cond) in piecewise_args
            if (cond is True) or (cond == sp.true)
        )
    else:
        Q_kerr_formula = Q_expr

    # Substitute a_spin -> 0 in the formula and L^2
    a_spin_sym = sp.Symbol("a_spin", real=True)
    Q_kerr_a0 = Q_kerr_formula.subs(a_spin_sym, 0)
    L_sq_kerr_a0 = L_sq_kerr.subs(a_spin_sym, 0)
    L_z_sq_kerr_a0 = L_z_sq_kerr.subs(a_spin_sym, 0)

    # Check identity: Q(a=0) + L_z^2 - L^2 == 0
    kerr_diff = sp.simplify(Q_kerr_a0 + L_z_sq_kerr_a0 - L_sq_kerr_a0)

    if kerr_diff != 0:
        raise ValueError(
            f"Error: In Kerr (limit a->0), Q + L_z^2 != L^2. Diff: {kerr_diff}"
        )
    print("PASS: Kerr (a->0) limit Q + L_z^2 == L^2 identity verified.")
    print("-" * 40)

    # Step 3: Generate trusted results
    for config in [
        "KerrSchild_Cartesian_massive",
        "KerrSchild_Cartesian_photon",
    ]:
        geo_diags = Geodesic_Diagnostics[config]
        results_dict = ve.process_dictionary_of_expressions(
            geo_diags.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{config}",
            results_dict,
        )

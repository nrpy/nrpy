"""
Construct symbolic expressions for conserved quantities along geodesics.

This module provides a class-based structure for generating the symbolic
expressions for quantities that are conserved along geodesics in certain
spacetimes, such as energy, angular momentum, and the Carter constant.
It supports both massive and massless (photon) particles.

Author: Dalton J. Moone

--------------------------------------------------------------------------------
List of Assumptions
--------------------------------------------------------------------------------
Physical Assumptions:
1.  **Metric Signature:** The code calculates Energy as E = -p_0. This assumes a
    metric signature of (-, +, +, +). If the signature were (+, -, -, -), Energy
    would be p_0.
2.  **Stationarity:** The calculation of Energy assumes the spacetime possesses
    a timelike Killing vector xi = d/dt, meaning the metric components are
    independent of time t.
3.  **Axial Symmetry:** The calculation of Angular Momentum (L_z) assumes the
    spacetime possesses an axial Killing vector xi = d/dphi, meaning the metric
    components are independent of phi.
4.  **Mass Normalization (m=1):** In the massive particle case, the code uses
    the term (1 - E**2). The general form is (m^2 - E^2). This implies the code
    assumes the particle mass m=1 (or that the affine parameter is normalized
    such that p_mu p^mu = -1).
5.  **Massless Particles:** For the massless case, the code assumes m=0,
    resulting in the term -E**2 (derived from 0 - E^2).

Coordinate & Geometric Assumptions:
6.  **Kerr-Schild / Boyer-Lindquist Correspondence:** The code assumes that the
    r and theta appearing in the Carter Constant formula are the standard
    Boyer-Lindquist radial and polar coordinates. It assumes the standard
    transformation between Cartesian Kerr-Schild.
7.  **Cartesian Input:** The code assumes self.xx corresponds to [t, x, y, z]
    where x, y, z are Cartesian-like coordinates.
8.  **Spin Parameter:** The symbol a_spin is assumed to be the angular momentum
    per unit mass (J/M) of the black hole.

Mathematical & Numerical Assumptions:
9.  **Contravariant Input Momentum:** The input pU is assumed to be the
    contravariant 4-momentum vector p^mu.
10. **Covariant Lowering:** The code assumes self.g4DD is the correct metric
    tensor for lowering indices (p_mu = g_mu_nu p^nu).
11. **Axis Regularization:** The code uses sp.Piecewise to return 0 when
    rho_sq < 1e-12. This assumes that on the z-axis (x=y=0), the computation
    of Q should default to 0 to avoid coordinate singularities.
12. **Definition of Q in Schwarzschild vs. Kerr:** The code assumes different
    definitions for Q depending on the spacetime tag. For "Schwarzschild", Q is
    defined as the total angular momentum squared (L^2). For "Kerr", Q is
    defined as the standard Boyer-Lindquist separation constant, which reduces
    to L^2 - L_z^2 in the non-spinning limit.
--------------------------------------------------------------------------------
"""

# Step 0.a: Import standard Python modules
import doctest
import logging
import os
import sys
from typing import Dict, List, Optional

# Step 0.b: Import third-party modules
import sympy as sp

# Step 0.c: Import nrpy core modules
import nrpy.indexedexp as ixp
import nrpy.validate_expressions.validate_expressions as ve
from nrpy.equations.general_relativity.geodesics.analytic_spacetimes import (
    Analytic_Spacetimes,
    a_spin,
)

# Step 0.d: Define global symbolic parameters for spin.
a_spin = sp.Symbol("a_spin", real=True)


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
        xx: Optional[List[sp.Symbol]] = None,
        g4DD: Optional[List[List[sp.Expr]]] = None,
    ) -> None:
        """
        Initialize and generate all symbolic diagnostic expressions.

        :param spacetime: The spacetime to use (e.g., "KerrSchild").
        :param particle_type: The type of particle ("massive" or "massless").
        :param xx: (Optional) List of symbolic coordinates [t, x, y, z].
                   Defaults to [x0, x1, x2, x3].
        :param g4DD: (Optional) Symbolic 4-metric tensor.
                     Defaults to a generic symmetric rank-2 tensor "g4DD".
        :raises ValueError: If the particle type is not supported.
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Initialize optional attributes
        self.E_expr = None
        self.L_exprs = []
        self.Q_expr = None

        if self.particle_type not in ["massive", "massless"]:
            raise ValueError(f"Unsupported particle_type: {self.particle_type}")

        # Step 1: Set up coordinates and metric (Defaults if not provided)
        if xx is None:
            self.xx = [sp.Symbol(f"x{i}", real=True) for i in range(4)]
        else:
            self.xx = xx

        if g4DD is None:
            self.g4DD = ixp.declarerank2("g4DD", symmetry="sym01", dimension=4)  # type: ignore
        else:
            self.g4DD = g4DD

        # Define 4-momentum vector pU = [p0, p1, p2, p3]
        pU = [sp.Symbol(f"p{i}", real=True) for i in range(4)]

        # Step 2: Compute quantities using modular functions.
        # Energy and Angular Momentum are computed regardless of particle type.
        self.E_expr = self.compute_energy(pU)
        self.L_exprs = self.compute_angular_momentum_cartesian(pU)

        # Carter Constant depends on particle type and spacetime.
        self.Q_expr = self.compute_carter_constant_cartesian(
            self.E_expr, self.L_exprs, pU
        )

    def compute_energy(self, pU: List[sp.Symbol]) -> sp.Expr:
        """
        Compute the conserved energy E = -p_0.

        This assumes the spacetime is stationary (has a timelike Killing vector d/dt).

        Reference:
        Wikipedia: Geodesics in general relativity
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Geodesics_in_general_relativity&oldid=1320086873

        :param pU: The symbolic 4-momentum vector.
        :return: Symbolic expression for Energy.
        """
        p_0 = sp.sympify(0)
        for mu in range(4):
            # Exploit symmetry: g_{0,mu} = g_{mu,0}. Use upper triangle.
            p_0 += self.g4DD[0][mu] * pU[mu]
        return -p_0

    def compute_angular_momentum_cartesian(self, pU: List[sp.Symbol]) -> List[sp.Expr]:
        """
        Compute the angular momentum vector L_i.

        Note: This implementation uses the cross-product formula L = x X p.
        This assumes the input coordinates `self.xx` are Cartesian-like, where
        xx[1]=x, xx[2]=y, xx[3]=z.

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

        # L = x x p
        L_x = y * p_z - z * p_y
        L_y = z * p_x - x * p_z
        L_z = x * p_y - y * p_x

        return [L_x, L_y, L_z]

    def compute_carter_constant_cartesian(
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
        Wikipedia: Carter constant
        Permanent Link: https://en.wikipedia.org/w/index.php?title=Carter_constant&oldid=1295996328

        :param E: Symbolic Energy expression.
        :param L_exprs: List of Angular Momentum components [Lx, Ly, Lz].
        :param pU: Symbolic 4-momentum vector.
        :return: Symbolic expression for Q.
        """
        L_x, L_y, L_z = L_exprs[0], L_exprs[1], L_exprs[2]

        if "Schwarzschild" in self.spacetime:
            # For Schwarzschild (a=0), Q simplifies to L^2.
            return L_x**2 + L_y**2 + L_z**2

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
        # Handle axial singularity
        return sp.Piecewise((sp.sympify(0), rho_sq < 1e-12), (Q_formula, True))


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

            # Retrieve analytic metric and coordinates
            metric_data = Analytic_Spacetimes[spacetime]

            self[key] = GeodesicDiagnostics(
                spacetime=spacetime,
                particle_type=particle_type,
                xx=metric_data.xx,
                g4DD=metric_data.g4DD,
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

    # Check 1: In Schwarzschild, Q must equal L^2 exactly.
    sch_diag = Geodesic_Diagnostics["Schwarzschild_Cartesian_Isotropic_massive"]
    L_sq_sch = (
        sch_diag.L_exprs[0] ** 2 + sch_diag.L_exprs[1] ** 2 + sch_diag.L_exprs[2] ** 2
    )
    sch_diff = sp.simplify(sch_diag.Q_expr - L_sq_sch)
    if sch_diff != 0:
        raise ValueError(f"Error: In Schwarzschild, Q != L^2. Diff: {sch_diff}")
    print("PASS: Schwarzschild Q == L^2 identity verified.")

    # Check 2: In Kerr (limit a->0), Q must reduce to L^2 - L_z^2.
    # Note: The general Kerr Q formula (Carter constant) reduces to L^2 - L_z^2
    # in the spherical limit, whereas the Schwarzschild case explicitly returns L^2.
    # We verify that Q_kerr(a=0) + L_z^2 == L^2.
    kerr_diag = Geodesic_Diagnostics["KerrSchild_Cartesian_massive"]

    L_sq_kerr = (
        kerr_diag.L_exprs[0] ** 2
        + kerr_diag.L_exprs[1] ** 2
        + kerr_diag.L_exprs[2] ** 2
    )
    L_z_sq_kerr = kerr_diag.L_exprs[2] ** 2

    # Extract the main formula from the Piecewise object to avoid coordinate substitution.
    # The Piecewise is ((0, rho < epsilon), (formula, True)). We want the formula.
    assert kerr_diag.Q_expr is not None
    Q_kerr_formula = kerr_diag.Q_expr.args[1][0]

    # Substitute a_spin -> 0 in the formula and L^2
    Q_kerr_a0 = Q_kerr_formula.subs(a_spin, 0)
    L_sq_kerr_a0 = L_sq_kerr.subs(a_spin, 0)
    L_z_sq_kerr_a0 = L_z_sq_kerr.subs(a_spin, 0)

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
        "Schwarzschild_Cartesian_Isotropic_massive",
        "KerrSchild_Cartesian_massless",
        "Schwarzschild_Cartesian_Isotropic_massless",
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

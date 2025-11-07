"""
Construct symbolic expressions for conserved quantities along geodesics.

This module provides a class-based structure for generating the symbolic
expressions for quantities that are conserved along geodesics in certain
spacetimes, such as energy, angular momentum, and the Carter constant.
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

# Step 0.d: Define global symbolic parameters for mass and spin.
M_scale = sp.Symbol("M_scale", real=True)
a_spin = sp.Symbol("a_spin", real=True)


class GeodesicDiagnostics:
    """
    Generate and store symbolic expressions for conserved quantities.

    This class is instantiated with a specific spacetime and particle type.
    It then computes the symbolic expressions for Energy (E), Angular Momentum (L),
    and the Carter Constant (Q), storing them as instance attributes.
    """

    def __init__(self, spacetime: str, particle_type: str = "massive") -> None:
        """
        Initialize and generate all symbolic diagnostic expressions.

        Args:
            spacetime: The spacetime to use (e.g., "KerrSchild").
            particle_type: The type of particle ("massive" or "massless").
        """
        self.spacetime = spacetime
        self.particle_type = particle_type

        # Step 1: Call the unified symbolic generation function.
        # This function contains the conditional logic for particle type.
        if self.particle_type not in ["massive", "massless"]:
            raise ValueError(f"Unsupported particle_type: {self.particle_type}")

        (
            self.E_expr,
            self.L_exprs,
            self.Q_expr,
        ) = self._symbolic_conserved_quantities(particle_type=self.particle_type)

    def _symbolic_conserved_quantities(
        self, particle_type: str
    ) -> Tuple[sp.Expr, List[sp.Expr], sp.Expr]:
        """
        Generate symbolic recipes for conserved quantities E, L_i, and Q.

        This unified function handles both massive and massless particles by
        implementing conditional logic for the Carter Constant Q.

        Args:
            particle_type: The type of particle ("massive" or "massless").

        Returns:
            A tuple (E_expr, L_exprs, Q_expr).
        """
        # Step 1: Define symbolic state vector and metric placeholder.
        # These are identical for both massive and massless cases.
        x, y, z = sp.symbols("y[1] y[2] y[3]", real=True)
        pt, px, py, pz = sp.symbols("y[4] y[5] y[6] y[7]", real=True)
        pU = [pt, px, py, pz]
        g4DD = ixp.declarerank2("metric->g", sym="sym01", dimension=4)

        # Step 2: Conserved Energy, E = -p_t = -g_{t,mu} p^mu.
        # This formula is the same for both particle types.
        p_t = sp.sympify(0)
        for mu in range(4):
            # Exploit symmetry: g_{0,mu} = g_{mu,0}. Always use upper triangle index order.
            p_t += g4DD[0][mu] * pU[mu]
        E_expr = -p_t

        # Step 3: Angular Momentum, L_i = epsilon_{ijk} x^j p_k.
        # This formula is the same for both particle types.
        # Step 3.a: First, compute covariant momentum p_k = g_{k,mu} p^mu.
        p_down = ixp.zerorank1(dimension=4)
        for k in range(1, 4):
            for mu in range(4):
                # *** THE FIX IS HERE ***
                # CORRECTED: Exploit metric symmetry for correctness.
                # Always access the metric with the smaller index first to ensure
                # we only reference components that exist in the C struct (g01, not g10).
                if k <= mu:
                    p_down[k] += g4DD[k][mu] * pU[mu]
                else:
                    p_down[k] += g4DD[mu][k] * pU[mu]

        p_x, p_y, p_z = p_down[1], p_down[2], p_down[3]
        # Step 3.b: Compute the components of L = x x p.
        L_x, L_y, L_z = y * p_z - z * p_y, z * p_x - x * p_z, x * p_y - y * p_x
        L_exprs = [L_x, L_y, L_z]

        # Step 4: Carter Constant, Q.
        # This is the only quantity with a particle-type-dependent formula.
        if "Schwarzschild" in self.spacetime:
            # For Schwarzschild (a=0), Q simplifies to L^2 for both particle types.
            Q_expr = L_x**2 + L_y**2 + L_z**2
        else:  # Kerr spacetime (a != 0)
            # The Boyer-Lindquist r coordinate is needed for the Carter constant formula.
            # It is found by solving r^4 - (x^2+y^2+z^2-a^2)r^2 - a^2 z^2 = 0 for r^2.
            Sigma = x**2 + y**2 + z**2 - a_spin**2
            r_sq = (Sigma + sp.sqrt(Sigma**2 + 4 * a_spin**2 * z**2)) / 2

            # p_θ^2 in Cartesian coordinates
            rho_sq = x**2 + y**2
            xpx_plus_ypy = x * p_x + y * p_y
            p_theta_sq = (
                (z**2 * xpx_plus_ypy**2 / rho_sq)
                - (2 * z * p_z * xpx_plus_ypy)
                + (rho_sq * p_z**2)
            )

            # --- Main Conditional Logic for Q ---
            if particle_type == "massive":
                # Carter Constant for massive particles (mass m=1):
                # Q = p_θ² + cos²θ * (a²(m² - E²) + L_z²/sin²θ)
                # In Cartesian-like coords, cos²θ = z²/r², sin²θ = ρ²/r²
                second_term = (z**2 / r_sq) * (
                    a_spin**2 * (1 - E_expr**2) + L_z**2 * (r_sq / rho_sq)
                )
            elif particle_type == "massless":
                # Carter Constant for massless particles (photons):
                # Q = p_θ² + cos²θ * (a² * E² - L_z²/sin²θ)
                second_term = (z**2 / r_sq) * (
                    a_spin**2 * E_expr**2 - L_z**2 * (r_sq / rho_sq)
                )
            else:
                # This should not be reached due to the check in __init__
                raise ValueError(f"Unsupported particle_type for Carter Constant: {particle_type}")

            Q_formula = p_theta_sq + second_term
            # Handle the axial singularity where rho_sq -> 0
            Q_expr = sp.Piecewise((sp.sympify(0), rho_sq < 1e-12), (Q_formula, True))

        return E_expr, L_exprs, Q_expr


class GeodesicDiagnostics_dict(Dict[str, GeodesicDiagnostics]):
    """A caching dictionary for GeodesicDiagnostics instances."""

    def __getitem__(self, key: str) -> GeodesicDiagnostics:
        """
        Get or create a GeodesicDiagnostics instance for a given configuration.

        Args:
            key: A string key in the format "Spacetime_ParticleType".

        Returns:
            A GeodesicDiagnostics instance for the specified configuration.
        """
        if key not in self:
            parts = key.split("_")
            if len(parts) < 2:
                raise ValueError("Invalid key format. Expected 'Spacetime_ParticleType'.")
            particle_type = parts[-1]
            spacetime = "_".join(parts[:-1])

            print(
                f"Setting up GeodesicDiagnostics for spacetime='{spacetime}', particle='{particle_type}'..."
            )
            self[key] = GeodesicDiagnostics(
                spacetime=spacetime, particle_type=particle_type
            )
        return super().__getitem__(key)


Geodesic_Diagnostics = GeodesicDiagnostics_dict()


if __name__ == "__main__":
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for config in [
        "KerrSchild_massive",
        "Schwarzschild_Cartesian_massive",
        "KerrSchild_massless",
        "Schwarzschild_Cartesian_massless",
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
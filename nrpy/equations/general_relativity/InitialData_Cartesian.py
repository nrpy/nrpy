"""
Set up initial data for solving Einstein's equations of general relativity, for data most naturally specified in Cartesian coordinates.

Outputs ADM quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp

# NRPy: This module depends on the parameter EvolvedConformalFactor_cf,
#        which is defined in BSSN.BSSN_quantities
import nrpy.equations.general_relativity.BSSN_quantities  # pylint: disable=unused-import
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN


class InitialData_Cartesian:
    """Construct and store Cartesian initial data for Einstein's equations of general relativity, as ADM quantities."""

    def __init__(self, IDtype: str, override_gauge_with_standard: bool = False, **kwargs) -> None:
        """
        Initialize Cartesian initial data.

        :param IDtype: The initial data type.
        :param override_gauge_with_standard: Whether to override gauge with standard values or not.
        :param kwargs: Keyword arguments for specific ID types.

        :raises ValueError: If IDtype not a supported option or kwargs are missing.
        """
        self.IDtype = IDtype

        self.gammaDD = ixp.zerorank2()
        self.KDD = ixp.zerorank2()
        self.alpha = sp.sympify(0)
        self.betaU = ixp.zerorank1()
        self.BU = ixp.zerorank1()

        self.x, self.y, self.z = sp.symbols("x y z", real=True)

        if IDtype == "BrillLindquist":
            ID_defines_gauge_quantities = False
            self.gammaDD, self.KDD = self.BrillLindquist()
        elif IDtype == "MassiveParticle_StableCircularOrbit":
            ID_defines_gauge_quantities = True
            if "M_scale" not in kwargs or "a_spin" not in kwargs:
                raise ValueError(
                    "M_scale and a_spin must be provided for MassiveParticle_StableCircularOrbit ID."
                )
            self.massive_particle_stable_circular_orbit(
                kwargs["M_scale"], kwargs["a_spin"]
            )
        else:
            raise ValueError(f"IDtype = {IDtype} is not supported.")

        if not ID_defines_gauge_quantities or override_gauge_with_standard:
            # Set shift betaU and time derivative of shift BU to zero -- the standard choice.
            self.betaU = ixp.zerorank1()
            self.BU = ixp.zerorank1()
            # Next compute alpha. Standard choice is alpha = 1/psi**2 (psi = BSSN conformal factor),
            #   where psi = exp(phi); chi = 1/psi**4; W = 1/psi**2
            adm2bssn = ADM_to_BSSN(
                self.gammaDD,
                self.KDD,
                self.betaU,
                self.BU,
                "Cartesian",
                compute_cf_only=True,
            )
            EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
            cf = adm2bssn.cf
            if EvolvedConformalFactor_cf == "phi":
                self.alpha = sp.exp(-2 * cf)
            elif EvolvedConformalFactor_cf == "chi":
                self.alpha = sp.sqrt(cf)
            elif EvolvedConformalFactor_cf == "W":
                self.alpha = cf
            else:
                raise ValueError(
                    f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
                )

    # fmt: off
    def BrillLindquist(self) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
        """
        Calculate Brill-Lindquist initial data for a pair of black holes.

        :return: Tuple containing gammaDD and KDD.
        """
        BH1_posn_x, BH1_posn_y, BH1_posn_z = par.register_CodeParameters(
            "REAL", __name__, ["BH1_posn_x", "BH1_posn_y", "BH1_posn_z"], [0.0, 0.0, +0.5], commondata=True
        )
        BH1_mass = par.register_CodeParameter("REAL", __name__, "BH1_mass", 0.5, commondata=True)
        BH2_posn_x, BH2_posn_y, BH2_posn_z = par.register_CodeParameters(
            "REAL", __name__, ["BH2_posn_x", "BH2_posn_y", "BH2_posn_z"], [0.0, 0.0, -0.5], commondata=True
        )
        BH2_mass = par.register_CodeParameter("REAL", __name__, "BH2_mass", 0.5, commondata=True)

        x = self.x
        y = self.y
        z = self.z

        # Step 2.b: Set psi, the conformal factor:
        psi = sp.sympify(1)
        psi += BH1_mass / ( 2 * sp.sqrt((x-BH1_posn_x)**2 + (y-BH1_posn_y)**2 + (z-BH1_posn_z)**2) )
        psi += BH2_mass / ( 2 * sp.sqrt((x-BH2_posn_x)**2 + (y-BH2_posn_y)**2 + (z-BH2_posn_z)**2) )

        # Step 2.c: Set all needed ADM variables in Cartesian coordinates
        gammaDD = ixp.zerorank2()
        KDD     = ixp.zerorank2() # K_{ij} = 0 for these initial data
        for i in range(3):
            gammaDD[i][i] = psi**4

        return gammaDD, KDD
    # fmt: on

    def massive_particle_stable_circular_orbit(
        self, M_scale: sp.Expr, a_spin: sp.Expr
    ) -> None:
        """
        Generate symbolic expressions for stable, circular, equatorial orbits.

        This recipe is stored as instance attributes ut_stable_expr and uphi_stable_expr.
        It provides the 4-velocity components u^t and u^phi needed to place a
        massive particle on a stable circular orbit at a given radius 'r_initial'
        in a Kerr spacetime. It uses a numerically stable three-step method.

        Args:
            M_scale: The symbolic parameter for the black hole's mass.
            a_spin: The symbolic parameter for the black hole's dimensionless spin.
        """
        r_initial = sp.Symbol("r_initial", real=True)

        # Step 1: Angular velocity Omega = d(phi)/dt for stable prograde orbits.
        Omega = sp.sqrt(M_scale) / (
            r_initial ** sp.Rational(3, 2) + a_spin * sp.sqrt(M_scale)
        )

        # Step 2: Define required Boyer-Lindquist metric components in the equatorial plane.
        g_tt = -(1 - 2 * M_scale / r_initial)
        g_tphi = -2 * a_spin * M_scale / r_initial
        g_phiphi = (
            r_initial**2 + a_spin**2 + (2 * M_scale * a_spin**2) / r_initial
        )

        # Step 3: Solve for u^t from the normalization condition g_munu u^mu u^nu = -1.
        ut_squared_inv_denom = g_tt + 2 * g_tphi * Omega + g_phiphi * Omega**2
        ut_squared = -1 / ut_squared_inv_denom
        self.ut_stable_expr = sp.sqrt(ut_squared)

        # Step 4: Calculate u^phi = Omega * u^t.
        self.uphi_stable_expr = Omega * self.ut_stable_expr


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    for ID_type in [
        "BrillLindquist",
        "MassiveParticle_StableCircularOrbit",
    ]:
        if ID_type == "MassiveParticle_StableCircularOrbit":
            # Provide the required symbolic parameters for this specific test case
            M_param = sp.Symbol("M_scale", real=True)
            a_param = sp.Symbol("a_spin", real=True)
            ID = InitialData_Cartesian(ID_type, M_scale=M_param, a_spin=a_param)
        else:
            ID = InitialData_Cartesian(ID_type)
        results_dict = ve.process_dictionary_of_expressions(
            ID.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{ID_type}",
            results_dict,
        )

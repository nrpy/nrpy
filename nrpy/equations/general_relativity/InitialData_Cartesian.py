"""
Set up initial data for solving Einstein's equations of general relativity, for data most naturally specified in Cartesian coordinates.

Outputs ADM quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Tuple, List
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN

# NRPy+: This module depends on the parameter EvolvedConformalFactor_cf,
#        which is defined in BSSN.BSSN_quantities
import nrpy.equations.general_relativity.BSSN_quantities  # pylint: disable=unused-import


class InitialData_Cartesian:
    """Construct and store Cartesian initial data for Einstein's equations of general relativity, as ADM quantities."""

    def __init__(self, IDtype: str, override_gauge_with_standard: bool = False) -> None:
        """
        Initialize Cartesian initial data.

        :param IDtype: The initial data type.
        :param override_gauge_with_standard: Whether to override gauge with standard values or not.
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
    ]:
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

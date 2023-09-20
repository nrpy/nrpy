"""
This module constructs expressions for ADM
  quantities in terms of BSSN quantities.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1.a: import all needed modules from NRPy+:
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities


class BSSN_to_ADM:
    """Class for setting up and storing expressions for ADM quantities in terms of BSSN quantities."""

    def __init__(
        self, CoordSystem: str = "Cartesian", enable_rfm_precompute: bool = False
    ):
        """
        Initialize the BSSN_to_ADM class, sets: gammaDD, gammaDDdD, gammaDDdDD, gammaUU, detgamma, GammaUDD, KDD, KDDdD

        :param CoordSystem: String describing the coordinate system of the inputs.
        :param enable_rfm_precompute: Boolean flag to enable reference metric precomputation
        """
        # Step 1.b: Set EvolvedConformalFactor_cf from NRPy+ parameter
        EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")

        # Step 1.c: Import all needed basic (unrescaled) BSSN scalars & tensors from BSSN_quantities
        Bq = BSSN_quantities[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        gammabarDD = Bq.gammabarDD
        cf = Bq.cf
        AbarDD = Bq.AbarDD
        trK = Bq.trK
        gammabarDD_dD = Bq.gammabarDD_dD
        gammabarDD_dDD = Bq.gammabarDD_dDD
        AbarDD_dD = Bq.AbarDD_dD

        # Step 2: The ADM three-metric gammaDD and its
        #         derivatives in terms of BSSN quantities.
        self.gammaDD = ixp.zerorank2()

        exp4phi = sp.sympify(0)
        if EvolvedConformalFactor_cf == "phi":
            exp4phi = sp.exp(4 * cf)
        elif EvolvedConformalFactor_cf == "chi":
            exp4phi = 1 / cf
        elif EvolvedConformalFactor_cf == "W":
            exp4phi = 1 / cf**2
        else:
            raise ValueError(
                f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
            )

        for i in range(3):
            for j in range(3):
                self.gammaDD[i][j] = exp4phi * gammabarDD[i][j]

        # Step 2.a: Derivatives of $e^{4\phi}$
        phidD = ixp.zerorank1()
        phidDD = ixp.zerorank2()
        cf_dD = ixp.declarerank1("cf_dD")
        cf_dDD = ixp.declarerank2("cf_dDD", symmetry="sym01")
        if EvolvedConformalFactor_cf == "phi":
            for i in range(3):
                phidD[i] = cf_dD[i]
                for j in range(3):
                    phidDD[i][j] = cf_dDD[i][j]
        elif EvolvedConformalFactor_cf == "chi":
            for i in range(3):
                phidD[i] = -sp.Rational(1, 4) * exp4phi * cf_dD[i]
                for j in range(3):
                    phidDD[i][j] = sp.Rational(1, 4) * (
                        exp4phi**2 * cf_dD[i] * cf_dD[j] - exp4phi * cf_dDD[i][j]
                    )
        elif EvolvedConformalFactor_cf == "W":
            exp2phi = 1 / cf
            for i in range(3):
                phidD[i] = -sp.Rational(1, 2) * exp2phi * cf_dD[i]
                for j in range(3):
                    phidDD[i][j] = sp.Rational(1, 2) * (
                        exp4phi * cf_dD[i] * cf_dD[j] - exp2phi * cf_dDD[i][j]
                    )
        else:
            raise ValueError(
                f"Error EvolvedConformalFactor_cf type = {EvolvedConformalFactor_cf} unknown."
            )

        exp4phidD = ixp.zerorank1()
        exp4phidDD = ixp.zerorank2()
        for i in range(3):
            exp4phidD[i] = 4 * exp4phi * phidD[i]
            for j in range(3):
                exp4phidDD[i][j] = (
                    16 * exp4phi * phidD[i] * phidD[j] + 4 * exp4phi * phidDD[i][j]
                )

        # Step 2.b: Derivatives of self.gammaDD, the ADM three-metric
        self.gammaDDdD = ixp.zerorank3()
        self.gammaDDdDD = ixp.zerorank4()

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.gammaDDdD[i][j][k] = (
                        exp4phidD[k] * gammabarDD[i][j]
                        + exp4phi * gammabarDD_dD[i][j][k]
                    )
                    for l in range(3):
                        self.gammaDDdDD[i][j][k][l] = (
                            exp4phidDD[k][l] * gammabarDD[i][j]
                            + exp4phidD[k] * gammabarDD_dD[i][j][l]
                            + exp4phidD[l] * gammabarDD_dD[i][j][k]
                            + exp4phi * gammabarDD_dDD[i][j][k][l]
                        )

        # Step 2.c: 3-Christoffel symbols associated with ADM 3-metric gammaDD
        # Step 2.c.i: First compute the inverse 3-metric gammaUU:
        self.gammaUU, self.detgamma = ixp.symm_matrix_inverter3x3(self.gammaDD)

        self.GammaUDD = ixp.zerorank3()

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        self.GammaUDD[i][j][k] += (
                            sp.Rational(1, 2)
                            * self.gammaUU[i][l]
                            * (
                                self.gammaDDdD[l][j][k]
                                + self.gammaDDdD[l][k][j]
                                - self.gammaDDdD[j][k][l]
                            )
                        )

        # Step 3: Define ADM extrinsic curvature KDD and
        #         its first spatial derivatives KDDdD
        #         in terms of BSSN quantities
        self.KDD = ixp.zerorank2()

        for i in range(3):
            for j in range(3):
                self.KDD[i][j] = (
                    exp4phi * AbarDD[i][j]
                    + sp.Rational(1, 3) * self.gammaDD[i][j] * trK
                )

        self.KDDdD = ixp.zerorank3()
        trK_dD = ixp.declarerank1("trK_dD")
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.KDDdD[i][j][k] = (
                        exp4phidD[k] * AbarDD[i][j]
                        + exp4phi * AbarDD_dD[i][j][k]
                        + sp.Rational(1, 3)
                        * (
                            self.gammaDDdD[i][j][k] * trK
                            + self.gammaDD[i][j] * trK_dD[k]
                        )
                    )


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

    for Coord in [
        "Spherical",
        "SinhSpherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        enable_rfm_pre = "rfm_precompute" in Coord
        bssn2adm = BSSN_to_ADM(Coord.replace("_rfm_precompute", "", enable_rfm_pre))
        results_dict = ve.process_dictionary_of_expressions(
            bssn2adm.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )

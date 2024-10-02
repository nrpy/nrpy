"""
Compute the Expansion Function Theta in Spherical Coordinates.

This module calculates the Expansion Function Theta, defined as:

Theta = Dbar_i s^i - K + s^i s^j K_ij,

where:
- Dbar_i is the covariant derivative associated with the physical 3-metric gammabar_jk,
- K_ij is the extrinsic curvature,
- K = gammabar^jk K_jk is the trace of the extrinsic curvature,
- s^i is the unit normal to the horizon surface defined by the level-set function F(r, theta, phi) = r - h(theta, phi) = 0.

At a marginally trapped surface, Theta = 0.

Note that only the definition of F(r,theta,phi) is in Spherical coordinates; adjusting only
this should make the construction of Theta fully covariant.

This implementation follows the notation and methodology primarily from:
- Thornburg (https://arxiv.org/pdf/gr-qc/9508014)
- Gundlach (https://arxiv.org/pdf/gr-qc/9707050)
- Hui & Lin (https://arxiv.org/pdf/2404.16511)

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from collections import OrderedDict
from typing import Dict, List, cast

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

import nrpy.grid as gri
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities


class ExpansionFunctionThetaClass:
    """Compute and store the Expansion Function Theta in Spherical coordinates."""

    def __init__(
        self,
        CoordSystem: str = "Spherical",
        enable_rfm_precompute: bool = False,
    ):
        """
        Initialize and set up all necessary quantities for Theta computation in Spherical coordinates.

        :param CoordSystem: The coordinate system being used; must be "Spherical".
        :param enable_rfm_precompute: Whether to enable reference-metric precomputation; defaults to False.

        :raises ValueError: If CoordSystem is not set to "Spherical".
        """
        # Validate CoordSystem
        if CoordSystem != "Spherical":
            raise ValueError(
                f"Unsupported CoordSystem '{CoordSystem}'. Currently, only 'Spherical' is supported."
            )
        self.CoordSystem = CoordSystem

        # Setup reference metric based on the coordinate system and precompute flag
        self.rfm = refmetric.reference_metric[
            (
                self.CoordSystem + "_rfm_precompute"
                if enable_rfm_precompute
                else self.CoordSystem
            )
        ]

        # Import BSSN quantities specific to the coordinate system and precompute flag
        self.Bq = BSSN_quantities[
            self.CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
        ]

        # Compute the Expansion Function Theta
        self.Theta = self._compute_theta()

        # Store the result in an ordered dictionary for easy access and extension
        self.Theta_varname_to_expr_dict: Dict[str, sp.Expr] = OrderedDict()
        self.Theta_varname_to_expr_dict["Theta"] = self.Theta

    def _compute_theta(self) -> sp.Expr:
        """
        Compute the Expansion Function Theta based on the defined equations.

        The computation follows these steps:
        0. Classify h, gammabarDDdD, and KDD as gridfunctions.
        1. Define the level-set function F and its first and second derivatives.
        2. Compute the derivatives of the inverse 3-metric gammabar^{ij}.
        3. Compute the unnormalized normal vector s_i and normalize it to obtain s^i.
        4. Compute the derivative of the normalization factor lambda.
        5. Compute the divergence of the unit normal vector s^i.
        6. Assemble the final expression for Theta.

        :return: The symbolic expression for Theta.
        """
        # Step 0: Classify h, gammabarDDdD, and KDD as gridfunctions.
        # Register gridfunctions for 'hh', 'gammabarDDdD' (partial_k gammabar_mn), and
        #   KDD (extrinsic curvature tensor K_ij), if not already registered. Otherwise
        #   just declare the variables.
        if "hh" not in gri.glb_gridfcs_dict:
            self.h = gri.register_gridfunctions(
                "hh", gf_array_name="evol_gfs", wavespeed=1.0
            )
            self.partialrhDD = gri.register_gridfunctions_for_single_rank2(
                "partialrhDD", symmetry="sym01", gf_array_name="auxevol_gfs"
            )
            self.KDD = gri.register_gridfunctions_for_single_rank2(
                "KDD", symmetry="sym01", gf_array_name="auxevol_gfs"
            )
        else:
            self.h = sp.symbols("hh", real=True)
            self.partialrhDD = ixp.declarerank2("partialrhDD", symmetry="sym01")
            self.KDD = ixp.declarerank2("KDD", symmetry="sym01")

        hDD_dD = cast(
            List[List[List[sp.Expr]]], ixp.declarerank3("hDD_dD", symmetry="sym01")
        )
        # Overwrite what hDD_dD[i][j][0] means:
        for i in range(3):
            for j in range(3):
                # MYPY ERROR HERE:
                hDD_dD[i][j][0] = self.partialrhDD[i][j]

        # Stolen from BSSN_quantities.py:
        self.gammabarDDdD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.gammabarDDdD[i][j][k] = (
                        self.rfm.ghatDDdD[i][j][k]
                        + hDD_dD[i][j][k] * self.rfm.ReDD[i][j]
                        + self.Bq.hDD[i][j] * self.rfm.ReDDdD[i][j][k]
                    )

        # Step 1: Define the level-set function F and its derivatives
        # Declare variables that will be computed from finite differencing.
        self.h_dD = ixp.declarerank1("hh_dD")
        self.h_dDD = ixp.declarerank2("hh_dDD", symmetry="sym01")
        # F(r, theta, phi) = r - h(theta, phi)
        # Partial derivatives of F with respect to coordinates: F_dD = [1, -h_theta, -h_phi]
        F_dD = [sp.sympify(1), -self.h_dD[1], -self.h_dD[2]]  # Partial derivatives of F
        F_dDD = ixp.zerorank2()  # Second partial derivatives of F
        # F_dDD[i][j] = second derivative of F with respect to coordinates i and j
        F_dDD[0][0] = sp.sympify(0)
        F_dDD[0][1] = F_dDD[1][0] = sp.sympify(0)
        F_dDD[0][2] = F_dDD[2][0] = sp.sympify(0)
        F_dDD[1][1] = -self.h_dDD[1][1]
        F_dDD[1][2] = F_dDD[2][1] = -self.h_dDD[1][2]
        F_dDD[2][2] = -self.h_dDD[2][2]

        # Step 2: Compute derivatives of the inverse 3-metric gammabar^{ij}
        # Using the identity: gammabar^{ij}_{,k} = -gammabar^{im} gammabar^{jn} gammabar_{mn,k}
        gammabarUU_dD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        for n in range(3):
                            gammabarUU_dD[i][j][k] += (
                                -self.Bq.gammabarUU[i][m]
                                * self.Bq.gammabarUU[j][n]
                                * self.gammabarDDdD[m][n][k]
                            )

        # Step 3: Compute the unnormalized normal vector s_i = Dbar_i F
        # Given F = r - h(theta, phi), we have Dbar_i F = [1, -h_theta, -h_phi]
        # Compute s^i = gammabar^{ij} F_j
        unnormalized_sU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                unnormalized_sU[i] += self.Bq.gammabarUU[i][j] * F_dD[j]

        # Step 4: Compute the normalization factor \lambda = (gammabar^{ij} \partial_i F\, \partial_j F)^{1/2}
        lambda_squared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                lambda_squared += self.Bq.gammabarUU[i][j] * F_dD[i] * F_dD[j]
        lamb = sp.sqrt(lambda_squared)

        # Step 5: Compute the unit normal vector s^i = \frac{gammabar^{ij} \partial_j F}{\lambda}
        sU = ixp.zerorank1()
        for i in range(3):
            sU[i] = unnormalized_sU[i] / lamb

        # Step 6: Compute the derivative of \lambda, i.e., \partial_i \lambda
        # \partial_i \lambda = \frac{1}{2\lambda} (gammabar^{km} \partial_k F \partial_m F \partial_i gammabar^{km} + 2 gammabar^{km} \partial_m F \partial_i \partial_k F)
        #                      ^^^ PREFACTOR ^^^^  ^^^^^^^^^^^^^^^^^^^^^^^^ TERM 1 ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^    ^^^^^^^^^^^^^^^^^^^^^^^^ TERM 2 ^^^^^^^^^^^^^^^^^^^^
        lamb_dD = ixp.zerorank1()
        # TERM 1:
        for i in range(3):
            for k in range(3):
                for m in range(3):
                    lamb_dD[i] += F_dD[k] * F_dD[m] * gammabarUU_dD[k][m][i]
        # TERM 2:
        for i in range(3):
            for k in range(3):
                for m in range(3):
                    lamb_dD[i] += 2 * self.Bq.gammabarUU[k][m] * F_dD[m] * F_dDD[i][k]
        # PREFACTOR
        for i in range(3):
            lamb_dD[i] *= 1 / (2 * lamb)

        # Step 7: Compute the divergence of the unit normal vector \partial_i s^i
        # \partial_i s^i = \frac{1}{\lambda} (\partial_j F \partial_i gammabar^{ij} + gammabar^{ij} \partial_i \partial_j F) - \frac{1}{\lambda^2} \partial_i \lambda gammabar^{ij} \partial_j F
        #                  ^^ PREFACTOR 1 ^^ ^^^^ FIRST PARENTHETICAL TERM ^^^^^^^^   ^^^^ SECOND PARENTHETICAL TERM ^^^^^^  ^^^^ PREFACTOR 2 ^^^  ^^^^^^^^^ TERM AFTER PARENTHETICAL ^^^^^^^^^^
        partial_i_si_parenthetical_term = sp.sympify(0)
        # FIRST PARENTHETICAL TERM
        for i in range(3):
            for j in range(3):
                partial_i_si_parenthetical_term += F_dD[j] * gammabarUU_dD[i][j][i]
        # SECOND PARENTHETICAL TERM
        for i in range(3):
            for j in range(3):
                partial_i_si_parenthetical_term += (
                    self.Bq.gammabarUU[i][j] * F_dDD[i][j]
                )
        # PREFACTOR 1
        partial_i_si_parenthetical_term /= lamb

        partial_i_si_term_after_parenthetical = sp.sympify(0)
        # TERM AFTER PARENTHETICAL
        for i in range(3):
            for j in range(3):
                partial_i_si_term_after_parenthetical += (
                    lamb_dD[i] * self.Bq.gammabarUU[i][j] * F_dD[j]
                )
        # PREFACTOR 2
        partial_i_si_term_after_parenthetical /= -lamb * lamb

        # OVERALL SUM:
        partial_i_si = (
            partial_i_si_parenthetical_term + partial_i_si_term_after_parenthetical
        )

        # Step 8: Compute the covariant divergence of s^i
        # Dbar_i s^i = \frac{1}{2 gammabar} s^i \partial_i gammabar + \partial_i s^i
        covariant_divergence_of_s_i = sp.sympify(0)
        for i in range(3):
            covariant_divergence_of_s_i += sU[i] * self.Bq.detgammabar_dD[i]
        covariant_divergence_of_s_i *= 1 / (2 * self.Bq.detgammabar)
        covariant_divergence_of_s_i += partial_i_si

        # Step 9: Assemble the final expression for Theta
        # Theta = Dbar_i s^i - K + s^i s^j K_{ij}
        Theta = covariant_divergence_of_s_i
        for i in range(3):
            for j in range(3):
                Theta += -self.Bq.gammabarUU[i][j] * self.KDD[i][j]
        for i in range(3):
            for j in range(3):
                Theta += sU[i] * sU[j] * self.KDD[i][j]

        return Theta


# Class to manage different configurations of ExpansionFunctionThetaClass
class ExpansionFunctionThetaClass_dict(Dict[str, ExpansionFunctionThetaClass]):
    """Custom dictionary for storing ExpansionFunctionThetaClass objects."""

    def __getitem__(self, key: str) -> ExpansionFunctionThetaClass:
        """
        Retrieve an ExpansionFunctionThetaClass object based on the provided key.

        Supported keys:
        - "Spherical": Spherical coordinates without reference-metric precomputation.
        - "Spherical_rfm_precompute": Spherical coordinates with reference-metric precomputation enabled.

        :param key: The key representing the desired configuration.
        :return: The corresponding ExpansionFunctionThetaClass object.
        :raises KeyError: If an unsupported key is provided.
        """
        if key not in self:
            # Determine if reference-metric precompute is enabled based on the key
            if key == "Spherical":
                enable_rfm_precompute = False
            elif key == "Spherical_rfm_precompute":
                enable_rfm_precompute = True
            else:
                raise KeyError(
                    f"Unsupported key: '{key}'. Supported keys are 'Spherical' and 'Spherical_rfm_precompute'."
                )

            print(f"Setting up ExpansionFunctionThetaClass[{key}]...")
            # Create a new instance of ExpansionFunctionThetaClass with the specified settings
            self.__setitem__(
                key,
                ExpansionFunctionThetaClass(
                    CoordSystem="Spherical",
                    enable_rfm_precompute=enable_rfm_precompute,
                ),
            )
        return dict.__getitem__(self, key)

    def __setitem__(self, key: str, value: ExpansionFunctionThetaClass) -> None:
        """
        Set an ExpansionFunctionThetaClass object with the provided key.

        Only "Spherical" and "Spherical_rfm_precompute" are supported.

        :param key: The key representing the desired configuration.
        :param value: The ExpansionFunctionThetaClass object to store.
        :raises KeyError: If an unsupported key is provided.
        """
        if key not in ["Spherical", "Spherical_rfm_precompute"]:
            raise KeyError(
                f"Unsupported key: '{key}'. Supported keys are 'Spherical' and 'Spherical_rfm_precompute'."
            )
        dict.__setitem__(self, key, value)

    def __delitem__(self, key: str) -> None:
        """
        Delete an ExpansionFunctionThetaClass object based on the provided key.

        :param key: The key representing the configuration to delete.
        """
        dict.__delitem__(self, key)


# Instantiate the custom dictionary for ExpansionFunctionThetaClass objects
ExpansionFunctionTheta = ExpansionFunctionThetaClass_dict()

# Restore the doctest at the bottom
if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve

    # Run doctests
    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Validate expressions for the two specialized cases
    for validation_key in ["Spherical", "Spherical_rfm_precompute"]:
        theta_calc = ExpansionFunctionTheta[validation_key]
        # Process the expressions in the theta_calc object
        results_dict = ve.process_dictionary_of_expressions(
            theta_calc.__dict__, fixed_mpfs_for_free_symbols=True
        )
        # Compare or generate trusted results
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            # trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{validation_key}",
            results_dict,
        )

"""
Compute the Expansion Function Theta in Spherical Coordinates.

Calculates the Expansion Function Theta, defined as:

Theta = D_i s^i - K + s^i s^j K_ij,

where:
- D_i is the covariant derivative associated with the physical 3-metric gamma_jk,
- K_ij is the extrinsic curvature,
- K = gamma^jk K_jk is the trace of the extrinsic curvature,
- s^i is the unit normal to the horizon surface defined by the level-set function F(r, theta, phi) = r - h(theta, phi) = 0.

At a marginally trapped surface, Theta = 0.

Note that only the definition of F(r,theta,phi) and derivatives of hDD are in Spherical coordinates; adjusting
these should make the construction of Theta fully covariant.

Follows the notation and methodology primarily from:
- Thornburg (https://arxiv.org/pdf/gr-qc/9508014)
- Gundlach (https://arxiv.org/pdf/gr-qc/9707050)
- Hui & Lin (https://arxiv.org/pdf/2404.16511)

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List, cast

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression support
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support


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

        # Compute the Expansion Function Theta
        self.Theta = self._compute_theta()

    def _compute_theta(self) -> sp.Expr:
        """
        Compute the Expansion Function Theta based on the defined equations.

        The computation follows these steps:
        0. Register as gridfunctions:
           * h ("hh")
           * h_{ij} ("hDD")
           * W conformal factor ("WW"), W_{,k} ("partial_D_WW")
           * h_{ij,k} ("partial_D_hDD"), and
           * K_{ij} ("KDD").
        1. Define the level-set function F and its first and second derivatives.
        2. Compute needed metric quantities, including gamma_{ij}, gamma^{ij}, and their derivatives.
        3. Compute derivatives of the inverse 3-metric gamma^{ij}.
        4. Compute the unnormalized normal vector s^i = gamma^{ij} partial_j F.
        5. Compute the normalization factor lambda = sqrt(gamma^{ij} partial_i F partial_j F).
        6. Normalize the normal vector s^i by dividing by lambda.
        7. Compute the derivative of lambda, i.e., partial_i lambda.
        8. Compute the divergence of the unit normal vector partial_i s^i.
        9. Compute the covariant divergence of s^i.
        10. Compute K_{ij} from a_{ij} and trK.
        11. Assemble the final expression for Theta.

        :return: The symbolic expression for Theta.
        """
        # Step 0: Register h, h_{ij}, h_{ij,k}, W, W_{,k}, and K_{ij} as gridfunctions, if not already
        #         registered. If already registered, just declare as NRPy indexed expressions.

        # Step 1: Define the level-set function F and its derivatives
        # Declare variables that will be computed from finite differencing.
        h = sp.Symbol("hh", real=True)
        h_dD = ixp.declarerank1("hh_dD")
        h_dDD = ixp.declarerank2("hh_dDD", symmetry="sym01")
        # F(r, theta, phi) = r - h(theta, phi)
        # Partial derivatives of F with respect to coordinates: F_dD = [1, -h_theta, -h_phi]
        F_dD = [sp.sympify(1), -h_dD[1], -h_dD[2]]  # Partial derivatives of F
        F_dDD = ixp.zerorank2()  # Second partial derivatives of F
        # F_dDD[i][j] = second derivative of F with respect to coordinates i and j
        F_dDD[0][0] = sp.sympify(0)
        F_dDD[0][1] = F_dDD[1][0] = sp.sympify(0)
        F_dDD[0][2] = F_dDD[2][0] = sp.sympify(0)
        F_dDD[1][1] = -h_dDD[1][1]
        F_dDD[1][2] = F_dDD[2][1] = -h_dDD[1][2]
        F_dDD[2][2] = -h_dDD[2][2]

        # Step 2: Next compute needed metric quantities
        # Stolen from BSSN_quantities.py:
        # gammabar_{ij}  = h_{ij}*ReDD[i][j] + gammahat_{ij}
        gammabarDD = ixp.zerorank2()
        hDD = ixp.declarerank2("hDD", symmetry="sym01")
        for i in range(3):
            for j in range(3):
                gammabarDD[i][j] = (
                    hDD[i][j] * self.rfm.ReDD[i][j] + self.rfm.ghatDD[i][j]
                )
        # gamma_{ij} = e^{4 phi} gammabar_{ij}, but
        #   W = e^{-2 phi}
        # -> gamma_{ij} = 1/W^2 gammabar_{ij}
        self.gammaDD = ixp.zerorank2()
        W = sp.symbols("WW", real=True)
        for i in range(3):
            for j in range(3):
                self.gammaDD[i][j] = 1 / W**2 * gammabarDD[i][j]
        gammabarDDdD = ixp.zerorank3()
        # IMPORTANT: partial_D_hDD[i][j][k] = h_{jk,i}, so sym12 is correct here:
        partial_D_hDD = cast(
            List[List[List[sp.Expr]]],
            ixp.declarerank3("partial_D_hDD", symmetry="sym12"),
        )
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    # Stolen from BSSN_quantities.py:
                    # gammabar_{ij,k} = gammahat_{ij,k} + h_{ij,k} ReDD[i][j] + h_{ij} ReDDdD[i][j][k]
                    gammabarDDdD[i][j][k] = (
                        self.rfm.ghatDDdD[i][j][k]
                        + partial_D_hDD[k][i][j] * self.rfm.ReDD[i][j]
                        + hDD[i][j] * self.rfm.ReDDdD[i][j][k]
                    )

        # W = e^{-2 phi}
        # -> gamma_{ij} = 1/W^2 gammabar_{ij}
        # -> gamma_{ij,k} = -2 1/W^3 gammabar_{ij} W_{,k} + 1/W^2 gammabar_{ij,k}
        gammaDDdD = ixp.zerorank3()
        WdD = ixp.declarerank1("partial_D_WW")
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    gammaDDdD[i][j][k] = (
                        -2 / W**3 * gammabarDD[i][j] * WdD[k]
                        + gammabarDDdD[i][j][k] / W**2
                    )

        # Compute derivatives of det(gamma) using Jacobi's formula:
        # detgamma_{,k} = detgamma * gamma^{ij} gamma_{ij,k}
        self.gammaUU, self.detgamma = ixp.symm_matrix_inverter3x3(self.gammaDD)
        self.detgamma_dD = ixp.zerorank1()
        for k in range(3):
            for j in range(3):
                for i in range(3):
                    self.detgamma_dD[k] += (
                        self.detgamma * self.gammaUU[i][j] * gammaDDdD[i][j][k]
                    )

        # Step 3: Compute derivatives of the inverse 3-metric gamma^{ij}
        # Using the identity: gamma^{ij}_{,k} = -gamma^{im} gamma^{jn} gamma_{mn,k}
        self.gammaUUdD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for m in range(3):
                        for n in range(3):
                            self.gammaUUdD[i][j][k] += (
                                -self.gammaUU[i][m]
                                * self.gammaUU[j][n]
                                * gammaDDdD[m][n][k]
                            )

        # Step 4: Compute the unnormalized normal vector s_i = D_i F
        # Given F = r - h(theta, phi), we have D_i F = [1, -h_theta, -h_phi]
        # Compute s^i = gamma^{ij} F_j
        unnormalized_sU = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                unnormalized_sU[i] += self.gammaUU[i][j] * F_dD[j]

        # Step 5: Compute the normalization factor lambda = (gamma^{ij} partial_i F partial_j F)^{1/2}
        lambda_squared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                #                 gamma^{ij}    partial_i F partial_j F
                lambda_squared += self.gammaUU[i][j] * F_dD[i] * F_dD[j]
        lamb = sp.sqrt(lambda_squared)

        # Step 6: Compute the unit normal vector s^i = (gamma^{ij} partial_j F) / (lambda)
        self.sU = ixp.zerorank1()
        for i in range(3):
            self.sU[i] = unnormalized_sU[i] / lamb

        # Step 7: Compute the derivative of lambda, i.e., partial_i lambda
        # partial_i lambda = 1/(2 lambda) [partial_k F partial_m F partial_i gamma^{km} + 2 gamma^{km} partial_m F partial_i partial_k F]
        #                      PREFACTOR   ^^^^^^^^^^^^^^^^^^ TERM 1 ^^^^^^^^^^^^^^^^^^   ^^^^^^^^^^^^^^^^^^^ TERM 2 ^^^^^^^^^^^^^^^^^^^
        lamb_dD = ixp.zerorank1()
        # TERM 1:
        for i in range(3):
            for k in range(3):
                for m in range(3):
                    #           partial_k F partial_m F partial_i gamma^{km}
                    lamb_dD[i] += F_dD[k] * F_dD[m] * self.gammaUUdD[k][m][i]
        # TERM 2:
        for i in range(3):
            for k in range(3):
                for m in range(3):
                    #             2    gamma^{km} partial_m F partial_i partial_k F
                    lamb_dD[i] += 2 * self.gammaUU[k][m] * F_dD[m] * F_dDD[i][k]
        # PREFACTOR
        for i in range(3):
            #              1/(2 lambda)
            lamb_dD[i] *= 1 / (2 * lamb)

        # Step 8: Compute the divergence of the unit normal vector partial_i s^i
        # partial_i s^i = (1 / lambda) (partial_j F partial_i gamma^{ij} + gamma^{ij} partial_i partial_j F) - (1/ lambda^2) partial_i lambda gamma^{ij} partial_j F
        #                  PREFACTOR 1   ^^^ FIRST PARENTHETICAL TERM ^^^   ^^ SECOND PARENTHETICAL TERM ^^^     PREFACTOR 2  ^^^^^^^ TERM AFTER PARENTHETICAL ^^^^^^
        partial_i_si_parenthetical_term = sp.sympify(0)
        # FIRST PARENTHETICAL TERM
        for i in range(3):
            for j in range(3):
                #                              partial_j F partial_i gamma^{ij}
                partial_i_si_parenthetical_term += F_dD[j] * self.gammaUUdD[i][j][i]
        # SECOND PARENTHETICAL TERM
        for i in range(3):
            for j in range(3):
                #                                    gamma^{ij} partial_i partial_j F
                partial_i_si_parenthetical_term += self.gammaUU[i][j] * F_dDD[i][j]
        # PREFACTOR 1
        #                             (1 / lambda)
        partial_i_si_parenthetical_term /= lamb

        partial_i_si_term_after_parenthetical = sp.sympify(0)
        # TERM AFTER PARENTHETICAL
        for i in range(3):
            for j in range(3):
                # partial_i lambda gamma^{ij} partial_j F
                partial_i_si_term_after_parenthetical += (
                    lamb_dD[i] * self.gammaUU[i][j] * F_dD[j]
                )
        # PREFACTOR 2
        #                                     - (1/ lambda^2)
        partial_i_si_term_after_parenthetical /= -lamb * lamb

        # OVERALL SUM
        partial_i_si = (
            partial_i_si_parenthetical_term + partial_i_si_term_after_parenthetical
        )

        # Step 9: Compute the covariant divergence of s^i
        # D_i s^i = 1/(2 detgamma) s^i partial_i detgamma + partial_i s^i
        self.covariant_divergence_of_s_i = sp.sympify(0)
        for i in range(3):
            self.covariant_divergence_of_s_i += self.sU[i] * self.detgamma_dD[i]
        self.covariant_divergence_of_s_i *= 1 / (2 * self.detgamma)
        self.covariant_divergence_of_s_i += partial_i_si

        # Step 10: Compute K_{ij} from a_{ij} and trK
        aDD = ixp.declarerank2("aDD", symmetry="sym01")
        trK = sp.symbols("trK", real=True)
        self.AbarDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                # Abar_{ij}      = a_{ij}*ReDD[i][j]
                self.AbarDD[i][j] = aDD[i][j] * self.rfm.ReDD[i][j]
        self.KDD = ixp.zerorank2()
        exp4phi = 1 / W**2
        for i in range(3):
            for j in range(3):
                self.KDD[i][j] = (
                    exp4phi * self.AbarDD[i][j]
                    + sp.Rational(1, 3) * self.gammaDD[i][j] * trK
                )

        # Step 11: Assemble the final expression for Theta
        # Theta = D_i s^i - K + s^i s^j K_{ij}
        Theta = self.covariant_divergence_of_s_i
        # Subtract the trace of the extrinsic curvature K = gamma^{ij} K_{ij}
        Theta -= trK
        # Add s^i s^j K_{ij}
        for i in range(3):
            for j in range(3):
                Theta += self.sU[i] * self.sU[j] * self.KDD[i][j]

        return Theta.subs(self.rfm.xx[0], h).subs(sp.sympify("f0_of_xx0"), h)


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

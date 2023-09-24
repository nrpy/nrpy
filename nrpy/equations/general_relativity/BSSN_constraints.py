"""
This module constructs expressions for the
  BSSN Hamiltonian and momentum constraint equations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1: Initialize needed Python/NRPy+ modules
from typing import Dict
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support

# NRPy+: Computes useful BSSN quantities
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity import T4munu

par.register_param(bool, __name__, "register_MU_gridfunctions", False)


class BSSNconstraints:
    """Sets up and stores expressions for BSSN constraints."""

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
    ) -> None:
        register_MU_gridfunctions = par.parval_from_str("register_MU_gridfunctions")

        # Step 1.b: Given the chosen coordinate system, set up
        #           corresponding reference metric and needed
        #           reference metric quantities
        # The following function call sets up the reference metric
        #    and related quantities, including rescaling matrices ReDD,
        #    ReU, and hatted quantities.
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        Bq = BSSN_quantities[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]

        # Step 1.c: Register H and MU gridfunctions for
        #           Hamiltonian & momentum constraints,
        #           respectively.
        if "H" not in gri.glb_gridfcs_dict:
            _ = gri.register_gridfunctions(
                "H", group="AUX", gf_array_name="diagnostic_output_gfs"
            )
            _ = gri.register_gridfunctions(
                "MSQUARED", group="AUX", gf_array_name="diagnostic_output_gfs"
            )
        if register_MU_gridfunctions and "MU" not in gri.glb_gridfcs_dict:
            _ = gri.register_gridfunctions_for_single_rank1(
                "MU", group="AUX", gf_array_name="diagnostic_output_gfs"
            )

        # Step 2: Hamiltonian constraint.
        #################################
        # -={ HAMILTONIAN CONSTRAINT }=-
        #################################

        # Term 1: 2/3 K^2
        self.H = sp.Rational(2, 3) * Bq.trK**2

        # Term 2: -A_{ij} A^{ij}
        for i in range(3):
            for j in range(3):
                self.H += -Bq.AbarDD[i][j] * Bq.AbarUU[i][j]

        # Term 3a: trace(Rbar)
        Rbartrace = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                Rbartrace += Bq.gammabarUU[i][j] * Bq.RbarDD[i][j]

        # Term 3b: -8 \bar{\gamma}^{ij} \bar{D}_i \phi \bar{D}_j \phi = -8*phi_dBar_times_phi_dBar
        # Term 3c: -8 \bar{\gamma}^{ij} \bar{D}_i \bar{D}_j \phi      = -8*phi_dBarDD_contraction
        phi_dBar_times_phi_dBar = sp.sympify(0)  # Term 3b
        phi_dBarDD_contraction = sp.sympify(0)  # Term 3c
        for i in range(3):
            for j in range(3):
                phi_dBar_times_phi_dBar += (
                    Bq.gammabarUU[i][j] * Bq.phi_dBarD[i] * Bq.phi_dBarD[j]
                )
                phi_dBarDD_contraction += Bq.gammabarUU[i][j] * Bq.phi_dBarDD[i][j]

        # Add Term 3:
        self.H += Bq.exp_m4phi * (
            Rbartrace - 8 * (phi_dBar_times_phi_dBar + phi_dBarDD_contraction)
        )

        # Step 3: M^i, the momentum constraint
        ##############################
        # -={ MOMENTUM CONSTRAINT }=-
        ##############################
        self.MU = ixp.zerorank1()

        # Term 2: 6 A^{ij} \partial_j \phi:
        for i in range(3):
            for j in range(3):
                self.MU[i] += 6 * Bq.AbarUU[i][j] * Bq.phi_dD[j]

        # Term 3: -2/3 \bar{\gamma}^{ij} K_{,j}
        trK_dD = ixp.declarerank1(
            "trK_dD"
        )  # Not defined in BSSN_constraints; only trK_dupD is defined there.
        for i in range(3):
            for j in range(3):
                self.MU[i] += -sp.Rational(2, 3) * Bq.gammabarUU[i][j] * trK_dD[j]

        # Next evaluate the conformal covariant derivative \bar{D}_j \bar{A}_{lm}
        AbarDD_dBarD = ixp.zerorank3()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    AbarDD_dBarD[i][j][k] = Bq.AbarDD_dD[i][j][k]
                    for l in range(3):
                        AbarDD_dBarD[i][j][k] += (
                            -Bq.GammabarUDD[l][k][i] * Bq.AbarDD[l][j]
                        )
                        AbarDD_dBarD[i][j][k] += (
                            -Bq.GammabarUDD[l][k][j] * Bq.AbarDD[i][l]
                        )

        # Term 1: Contract twice with the metric to make \bar{D}_{j} \bar{A}^{ij}
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        self.MU[i] += (
                            Bq.gammabarUU[i][k]
                            * Bq.gammabarUU[j][l]
                            * AbarDD_dBarD[k][l][j]
                        )

        # Finally, we multiply by e^{-4 phi} and rescale the momentum constraint:
        for i in range(3):
            self.MU[i] *= Bq.exp_m4phi

        # Next add T4UUmunu source terms if desired.
        if par.parval_from_str("enable_T4munu"):
            if "T4UU00" not in gri.glb_gridfcs_dict:
                _ = gri.register_gridfunctions_for_single_rank2(
                    "T4UU",
                    symmetry="sym01",
                    dimension=4,
                    group="AUXEVOL",
                )

            sourceterm_H, sourceterm_MU = T4munu.BSSN_constraints_T4UU_source_terms(
                CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
            )
            self.H += sourceterm_H
            for i in range(3):
                self.MU[i] += sourceterm_MU[i]

        # Then compute M^2 = M^i M_i
        self.Msquared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.Msquared += Bq.gammabarDD[i][j] * self.MU[i] * self.MU[j]

        # Finally construct the rescaled MU:
        self.mU = ixp.zerorank1()
        for i in range(3):
            self.mU[i] = self.MU[i] / rfm.ReU[i]


class BSSNconstraints_dict(Dict[str, BSSNconstraints]):
    """Custom dictionary for storing BSSNconstraints objects."""

    def __getitem__(self, CoordSystem_in: str) -> BSSNconstraints:
        if CoordSystem_in not in self:
            # In case [CoordSystem]_rfm_precompute is passed:
            CoordSystem = CoordSystem_in.replace("_rfm_precompute", "")
            enable_T4munu = par.parval_from_str("enable_T4munu")
            print(
                f"Setting up BSSN_constraints for CoordSystem = {CoordSystem}, enable_T4munu={enable_T4munu}."
            )
            self.__setitem__(
                CoordSystem, BSSNconstraints(CoordSystem, enable_rfm_precompute=False)
            )
            self.__setitem__(
                CoordSystem + "_rfm_precompute",
                BSSNconstraints(CoordSystem, enable_rfm_precompute=True),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: BSSNconstraints) -> None:
        dict.__setitem__(self, CoordSystem, value)

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)


BSSN_constraints = BSSNconstraints_dict()


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

    # enable T4munu to ensure a maximally comprehensive test.
    par.set_parval_from_str("enable_T4munu", True)

    for Coord in [
        "Spherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        bconstraints = BSSN_constraints[Coord]
        results_dict = ve.process_dictionary_of_expressions(
            bconstraints.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )

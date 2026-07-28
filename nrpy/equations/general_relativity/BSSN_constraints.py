# nrpy/equations/general_relativity/BSSN_constraints.py
"""
Construct expressions for the BSSN Hamiltonian, momentum, and conformal connection constraint equations.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Step 1: Initialize needed Python/NRPy modules
from typing import Dict, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.grid as gri  # NRPy: Functions having to do with numerical grids
import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.params as par  # NRPy: Parameter interface
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.equations.general_relativity import T4munu

# NRPy: Computes useful BSSN quantities
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities

par.register_param(bool, __name__, "register_MU_gridfunctions", False)


class BSSNconstraints:
    """Set up and store Hamiltonian, momentum, and Lambda constraints."""

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        enable_RbarDD_gridfunctions: bool = False,
        enable_T4munu: bool = False,
    ) -> None:
        """
        Set up symbolic expressions for BSSN constraints, storing them within the class object.

        :param CoordSystem: The coordinate system being used, defaults to "Cartesian".
        :param enable_rfm_precompute: Whether to enable reference-metric precomputation, defaults to False.
        :param enable_RbarDD_gridfunctions: Whether to enable RbarDD gridfunctions, defaults to False.
        :param enable_T4munu: Whether to enable T4munu (stress-energy terms), defaults to False.
        """
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
            CoordSystem
            + ("_rfm_precompute" if enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
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
        if register_MU_gridfunctions and "MU0" not in gri.glb_gridfcs_dict:
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
        if enable_T4munu:
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

        # Then compute M^2 = gamma_{ij} M^i M^j, where
        # gamma_{ij} = gammabar_{ij} / exp(-4 phi).
        self.Msquared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.Msquared += Bq.gammabarDD[i][j] * self.MU[i] * self.MU[j]
        self.Msquared /= Bq.exp_m4phi

        # Finally construct the rescaled MU:
        self.mU = ixp.zerorank1()
        for i in range(3):
            self.mU[i] = self.MU[i] / rfm.ReU[i]

        # Step 4: Covariant conformal connection constraint.
        self.LambdaConstraintU = ixp.zerorank1()
        for i in range(3):
            self.LambdaConstraintU[i] = Bq.LambdabarU[i] - Bq.DGammaU[i]

        # Contract with the conformal metric.
        self.LambdaConstraintSquared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.LambdaConstraintSquared += (
                    Bq.gammabarDD[i][j]
                    * self.LambdaConstraintU[i]
                    * self.LambdaConstraintU[j]
                )
        self.LambdaConstraintMagnitude = sp.sqrt(self.LambdaConstraintSquared)


class BSSNconstraints_dict(Dict[str, BSSNconstraints]):
    """Custom dictionary for storing BSSNconstraints objects."""

    def __init__(self) -> None:
        """Initialize an empty cache and its construction-parameter metadata."""
        super().__init__()
        self._construction_parameters: Dict[str, Tuple[str, bool, bool]] = {}

    def __getitem__(self, CoordSystem_in: str) -> BSSNconstraints:
        """
        Return constraints built for the current construction parameters.

        :param CoordSystem_in: Coordinate-system and option cache key.
        :return: Cached or newly built BSSN constraints.

        Doctests:
        >>> import contextlib
        >>> import io
        >>> import nrpy.validate_expressions.validate_expressions as ve
        >>> from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
        >>> original_cf = par.parval_from_str("EvolvedConformalFactor_cf")
        >>> original_register_MU = par.parval_from_str("register_MU_gridfunctions")
        >>> original_gridfunctions = gri.glb_gridfcs_dict.copy()
        >>> original_parameters = par.glb_params_dict.copy()
        >>> original_code_parameters = par.glb_code_params_dict.copy()
        >>> original_reference_metrics = refmetric.reference_metric.copy()
        >>> original_quantities = dict(BSSN_quantities)
        >>> original_quantity_parameters = (
        ...     BSSN_quantities._construction_parameters.copy()
        ... )
        >>> cache = BSSNconstraints_dict()
        >>> objects = []
        >>> invalid_cf_raised = False
        >>> try:
        ...     for name in ("MU0", "MU1", "MU2"):
        ...         _ = gri.glb_gridfcs_dict.pop(name, None)
        ...     par.set_parval_from_str("register_MU_gridfunctions", False)
        ...     with contextlib.redirect_stdout(io.StringIO()):
        ...         for cf_choice in ("W", "chi", "phi"):
        ...             par.set_parval_from_str(
        ...                 "EvolvedConformalFactor_cf", cf_choice
        ...             )
        ...             constraint = cache["Cartesian"]
        ...             objects.append(constraint)
        ...             adm = BSSN_to_ADM("Cartesian")
        ...             physical_norm = sum(
        ...                 adm.gammaDD[i][j]
        ...                 * constraint.MU[i]
        ...                 * constraint.MU[j]
        ...                 for i in range(3)
        ...                 for j in range(3)
        ...             )
        ...             ve.assert_equal(
        ...                 constraint.Msquared,
        ...                 physical_norm,
        ...                 suppress_message=True,
        ...             )
        ...         par.set_parval_from_str("register_MU_gridfunctions", True)
        ...         registered_object = cache["Cartesian"]
        ...         same_parameters_reused = cache["Cartesian"] is registered_object
        ...         cartesian_quantities = BSSN_quantities["Cartesian"]
        ...         spherical_object = cache["Spherical"]
        ...         spherical_quantities = BSSN_quantities["Spherical"]
        ...         coordinate_contraction = sum(
        ...             spherical_quantities.gammabarDD[i][j]
        ...             * spherical_object.LambdaConstraintU[i]
        ...             * spherical_object.LambdaConstraintU[j]
        ...             for i in range(3)
        ...             for j in range(3)
        ...         )
        ...         lambda_magnitude_matches = (
        ...             registered_object.LambdaConstraintMagnitude
        ...             == sp.sqrt(registered_object.LambdaConstraintSquared)
        ...         )
        ...         controlled_values = {
        ...             "hDD00": sp.Integer(1),
        ...             "hDD01": sp.Integer(0),
        ...             "hDD02": sp.Integer(0),
        ...             "hDD11": -sp.Rational(1, 2),
        ...             "hDD12": sp.Integer(0),
        ...             "hDD22": sp.Integer(0),
        ...             "lambdaU0": sp.Integer(5),
        ...             "lambdaU1": sp.Integer(6),
        ...             "lambdaU2": sp.Integer(7),
        ...         }
        ...         for i in range(3):
        ...             for j in range(i, 3):
        ...                 for k in range(3):
        ...                     controlled_values[f"hDD_dD{i}{j}{k}"] = sp.Integer(0)
        ...         controlled_values["hDD_dD000"] = sp.Integer(3)
        ...         controlled_values["hDD_dD110"] = -sp.Rational(3, 4)
        ...         controlled_symbols = set().union(
        ...             *(
        ...                 expression.free_symbols
        ...                 for expression in registered_object.LambdaConstraintU
        ...             )
        ...         )
        ...         controlled_mapping = {
        ...             symbol: controlled_values[str(symbol)]
        ...             for symbol in controlled_symbols
        ...         }
        ...         cartesian_divergence = ixp.zerorank1()
        ...         for i in range(3):
        ...             for j in range(3):
        ...                 for k in range(3):
        ...                     for l in range(3):
        ...                         cartesian_divergence[i] += (
        ...                             -cartesian_quantities.gammabarUU[i][k]
        ...                             * cartesian_quantities.gammabarUU[j][l]
        ...                             * cartesian_quantities.gammabarDD_dD[k][l][j]
        ...                         )
        ...         actual_dict = {
        ...             "LambdaConstraintU": registered_object.LambdaConstraintU,
        ...             "LambdaConstraintSquared": (
        ...                 spherical_object.LambdaConstraintSquared
        ...             ),
        ...             "ReferenceDataCartesian": [
        ...                 expression.xreplace(
        ...                     {
        ...                         symbol: sp.sympify(0)
        ...                         for symbol in expression.free_symbols
        ...                         if str(symbol).startswith(("hDD", "lambdaU"))
        ...                     }
        ...                 )
        ...                 for expression in registered_object.LambdaConstraintU
        ...             ],
        ...             "ReferenceDataSpherical": [
        ...                 expression.xreplace(
        ...                     {
        ...                         symbol: sp.sympify(0)
        ...                         for symbol in expression.free_symbols
        ...                         if str(symbol).startswith(("hDD", "lambdaU"))
        ...                     }
        ...                 )
        ...                 for expression in spherical_object.LambdaConstraintU
        ...             ],
        ...             "CartesianGammaReduction": [
        ...                 expression.xreplace(controlled_mapping)
        ...                 for expression in registered_object.LambdaConstraintU
        ...             ],
        ...         }
        ...         expected_dict = {
        ...             "LambdaConstraintU": [
        ...                 cartesian_quantities.LambdabarU[i]
        ...                 - cartesian_quantities.DGammaU[i]
        ...                 for i in range(3)
        ...             ],
        ...             "LambdaConstraintSquared": coordinate_contraction,
        ...             "ReferenceDataCartesian": ixp.zerorank1(),
        ...             "ReferenceDataSpherical": ixp.zerorank1(),
        ...             "CartesianGammaReduction": [
        ...                 (
        ...                     cartesian_quantities.LambdabarU[i]
        ...                     + cartesian_divergence[i]
        ...                 ).xreplace(controlled_mapping)
        ...                 for i in range(3)
        ...             ],
        ...         }
        ...         ve.assert_equal(
        ...             actual_dict,
        ...             expected_dict,
        ...             suppress_message=True,
        ...         )
        ...         controlled_squared = (
        ...             registered_object.LambdaConstraintSquared.xreplace(
        ...                 controlled_mapping
        ...             )
        ...         )
        ...         controlled_magnitude = (
        ...             registered_object.LambdaConstraintMagnitude.xreplace(
        ...                 controlled_mapping
        ...             )
        ...         )
        ...         positive_magnitude_matches = (
        ...             controlled_squared.is_positive is True
        ...             and controlled_magnitude == sp.sqrt(controlled_squared)
        ...         )
        ...         par.set_parval_from_str("EvolvedConformalFactor_cf", "invalid")
        ...         try:
        ...             _ = cache["Cartesian"]
        ...         except ValueError:
        ...             invalid_cf_raised = True
        ...     mu_registered = all(
        ...         f"MU{i}" in gri.glb_gridfcs_dict for i in range(3)
        ...     )
        ... finally:
        ...     par.set_parval_from_str("EvolvedConformalFactor_cf", original_cf)
        ...     par.set_parval_from_str(
        ...         "register_MU_gridfunctions", original_register_MU
        ...     )
        ...     dict.clear(BSSN_quantities)
        ...     dict.update(BSSN_quantities, original_quantities)
        ...     BSSN_quantities._construction_parameters.clear()
        ...     BSSN_quantities._construction_parameters.update(
        ...         original_quantity_parameters
        ...     )
        ...     gri.glb_gridfcs_dict.clear()
        ...     gri.glb_gridfcs_dict.update(original_gridfunctions)
        ...     par.glb_params_dict.clear()
        ...     par.glb_params_dict.update(original_parameters)
        ...     par.glb_code_params_dict.clear()
        ...     par.glb_code_params_dict.update(original_code_parameters)
        ...     refmetric.reference_metric.clear()
        ...     refmetric.reference_metric.update(original_reference_metrics)
        >>> (
        ...     objects[0].H != objects[1].H,
        ...     objects[0] is not objects[1],
        ...     registered_object is not objects[-1],
        ...     mu_registered,
        ...     same_parameters_reused,
        ...     invalid_cf_raised,
        ...     lambda_magnitude_matches,
        ...     positive_magnitude_matches,
        ... )
        (True, True, True, True, True, True, True, True)
        """
        construction_parameters = (
            par.parval_from_str("EvolvedConformalFactor_cf"),
            par.parval_from_str("detgbarOverdetghat_equals_one"),
            par.parval_from_str("register_MU_gridfunctions"),
        )
        if (
            CoordSystem_in not in self
            or self._construction_parameters.get(CoordSystem_in)
            != construction_parameters
        ):
            # In case e.g., [CoordSystem]_rfm_precompute_T4munu or [CoordSystem]_rfm_precompute are passed:
            CoordSystem = (
                CoordSystem_in.replace("_rfm_precompute", "")
                .replace("_RbarDD_gridfunctions", "")
                .replace("_T4munu", "")
            )
            enable_rfm_precompute = "_rfm_precompute" in CoordSystem_in
            enable_RbarDD_gridfunctions = "_RbarDD_gridfunctions" in CoordSystem_in
            enable_T4munu = "_T4munu" in CoordSystem_in

            print(
                f"Setting up BSSN_constraints for CoordSystem = {CoordSystem}, enable_T4munu={enable_T4munu}, "
                f"rfm_precompute={enable_rfm_precompute}, Rij gridfuncs={enable_RbarDD_gridfunctions}."
            )
            self.__setitem__(
                CoordSystem_in,
                BSSNconstraints(
                    CoordSystem,
                    enable_rfm_precompute=enable_rfm_precompute,
                    enable_RbarDD_gridfunctions=enable_RbarDD_gridfunctions,
                    enable_T4munu=enable_T4munu,
                ),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: BSSNconstraints) -> None:
        dict.__setitem__(self, CoordSystem, value)
        self._construction_parameters[CoordSystem] = (
            par.parval_from_str("EvolvedConformalFactor_cf"),
            par.parval_from_str("detgbarOverdetghat_equals_one"),
            par.parval_from_str("register_MU_gridfunctions"),
        )

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)
        self._construction_parameters.pop(CoordSystem, None)

    def clear(self) -> None:
        """Remove all cached constraints and construction metadata."""
        dict.clear(self)
        self._construction_parameters.clear()


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

    for Coord in [
        "Spherical",
        "SinhSpherical_rfm_precompute",
        "Cartesian",
        "SinhCartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        bconstraints = BSSN_constraints[Coord + "_T4munu"]
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

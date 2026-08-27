"""
Construct fully covariant conformal Z4 (fCCZ4) evolution equations.

This module implements the three-dimensional reference-metric formulation of
Mewes et al. (2020), with the Lagrangian conformal-determinant condition and
``kappa1`` damping terms that contain no additional lapse factor.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Validation code below intentionally mirrors equation-construction indices.
# pylint: disable=redefined-outer-name

from collections import OrderedDict
from typing import Dict, Tuple

import sympy as sp

import nrpy.equations.general_relativity.fCCZ4_constraints as fccz4_constraints
import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_RHSs import BSSN_RHSs


class FCCZ4RHSs:
    """Construct and store fCCZ4 right-hand sides and Z4 constraints."""

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        enable_RbarDD_gridfunctions: bool = False,
        enable_T4munu: bool = False,
    ) -> None:
        """
        Build fCCZ4 expressions for one reference-metric configuration.

        :param CoordSystem: Reference-metric coordinate system.
        :param enable_rfm_precompute: Use precomputed reference-metric symbols.
        :param enable_RbarDD_gridfunctions: Read the BSSN-shaped conformal
            Ricci components, with ``lambdaU`` interpreted as
            ``LambdatildeU``, from auxiliary gridfunctions.
        :param enable_T4munu: Include stress-energy source terms.
        :raises ValueError: If existing ``Theta_fCCZ4`` storage is not evolved
            state or, for BHaH, is not backed by ``in_gfs``.

        """
        theta_gf = gri.glb_gridfcs_dict.get("Theta_fCCZ4")
        if theta_gf is not None:
            theta_array_name = getattr(theta_gf, "gf_array_name", None)
            if theta_gf.group != "EVOL" or theta_array_name not in (None, "in_gfs"):
                raise ValueError(
                    "Theta_fCCZ4 must be an EVOL gridfunction stored in in_gfs."
                )

        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        suffix = (
            CoordSystem
            + ("_rfm_precompute" if enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        )
        Bq = BSSN_quantities[suffix]
        Brhs = BSSN_RHSs[suffix + ("_T4munu" if enable_T4munu else "")]

        # Register damping coefficients only while constructing equations.  In
        # the 2020 convention kappa1 has inverse-length units; kappa2 is
        # dimensionless.  kappa3 is the literal one for full covariance.
        self.kappa1, self.kappa2 = par.register_CodeParameters(
            "REAL",
            __name__,
            ["kappa1", "kappa2"],
            [0.1, 0.0],
            commondata=True,
            descriptions=[
                "fCCZ4 damping coefficient with inverse-length units",
                "dimensionless fCCZ4 damping coefficient",
            ],
        )

        if theta_gf is not None:
            self.Theta = sp.symbols("Theta_fCCZ4", real=True)
        else:
            self.Theta = gri.register_gridfunctions(
                "Theta_fCCZ4", f_infinity=0.0, wavespeed=1.0, is_basename=False
            )[0]
        self.Theta_dD = ixp.declarerank1("Theta_fCCZ4_dD")
        self.Theta_dupD = ixp.declarerank1("Theta_fCCZ4_dupD")

        constraints = fccz4_constraints.fCCZ4_constraints[
            suffix + ("_T4munu" if enable_T4munu else "")
        ]

        # Copy cached tensor aggregates before adding fCCZ4 evolution terms.
        self.LambdatildeU = [constraints.LambdatildeU[i] for i in range(3)]
        self.LambdatildeU_dupD = ixp.zerorank2()
        self.Z4constraintU = ixp.zerorank1()
        self.ZbarU = ixp.zerorank1()
        self.ZU = ixp.zerorank1()
        self.ZD = ixp.zerorank1()
        self.RbarZ4DD_delta = ixp.zerorank2()
        self.RbarZ4DD = ixp.zerorank2()
        for i in range(3):
            self.Z4constraintU[i] = constraints.Z4constraintU[i]
            self.ZbarU[i] = constraints.ZbarU[i]
            self.ZU[i] = constraints.ZU[i]
            self.ZD[i] = constraints.ZD[i]
            for j in range(3):
                self.LambdatildeU_dupD[i][j] = Brhs.LambdabarU_dupD[i][j]
                self.RbarZ4DD_delta[i][j] = constraints.RbarZ4DD_delta[i][j]
                self.RbarZ4DD[i][j] = constraints.RbarZ4DD[i][j]

        self.RbarZ4 = constraints.RbarZ4
        self.H_Z4 = constraints.H_Z4

        self.cf_rhs = Brhs.cf_rhs
        self.trK_rhs = (
            Brhs.trK_rhs
            + Bq.alpha * (self.H_Z4 - 2 * self.Theta * Bq.trK)
            - 3 * self.kappa1 * (1 + self.kappa2) * self.Theta
        )
        self.Theta_rhs = (
            Bq.alpha * sp.Rational(1, 2) * (self.H_Z4 - 2 * self.Theta * Bq.trK)
            - self.kappa1 * (2 + self.kappa2) * self.Theta
        )
        alpha_dD = ixp.declarerank1("alpha_dD")
        for i in range(3):
            self.Theta_rhs += (
                Bq.betaU[i] * self.Theta_dupD[i] - self.ZU[i] * alpha_dD[i]
            )

        self.h_rhsDD = ixp.zerorank2()
        self.a_rhsDD = ixp.zerorank2()
        delta_trace = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                delta_trace += Bq.gammabarUU[i][j] * self.RbarZ4DD_delta[i][j]
        for i in range(3):
            for j in range(3):
                self.h_rhsDD[i][j] = Brhs.h_rhsDD[i][j]
                self.a_rhsDD[i][j] = (
                    Brhs.a_rhsDD[i][j]
                    + (
                        -2 * Bq.alpha * self.Theta * Bq.AbarDD[i][j]
                        + Bq.alpha
                        * Bq.exp_m4phi
                        * (
                            self.RbarZ4DD_delta[i][j]
                            - sp.Rational(1, 3) * Bq.gammabarDD[i][j] * delta_trace
                        )
                    )
                    / rfm.ReDD[i][j]
                )

        Dhat_div_beta = sp.sympify(0)
        for k in range(3):
            Dhat_div_beta += Bq.betaU_dD[k][k]
            for ell in range(3):
                Dhat_div_beta += rfm.GammahatUDD[k][ell][k] * Bq.betaU[ell]
        self.Lambdatilde_rhsU_delta = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                self.Lambdatilde_rhsU_delta[i] += (
                    2
                    * Bq.gammabarUU[i][j]
                    * (Bq.alpha * self.Theta_dD[j] - self.Theta * alpha_dD[j])
                )
            self.Lambdatilde_rhsU_delta[i] += (
                -sp.Rational(2, 3) * Bq.alpha * Bq.trK * self.Z4constraintU[i]
                - self.kappa1 * self.Z4constraintU[i]
                + sp.Rational(2, 3) * self.Z4constraintU[i] * Dhat_div_beta
            )
            for k in range(3):
                Dhat_beta = Bq.betaU_dD[i][k]
                for ell in range(3):
                    Dhat_beta += rfm.GammahatUDD[i][ell][k] * Bq.betaU[ell]
                self.Lambdatilde_rhsU_delta[i] += -self.Z4constraintU[k] * Dhat_beta

        self.Lambdatilde_rhsU = ixp.zerorank1()
        self.lambda_rhsU = ixp.zerorank1()
        for i in range(3):
            self.Lambdatilde_rhsU[i] = (
                Brhs.Lambdabar_rhsU[i] + self.Lambdatilde_rhsU_delta[i]
            )
            self.lambda_rhsU[i] = self.Lambdatilde_rhsU[i] / rfm.ReU[i]

        self.fCCZ4_RHSs_varname_to_expr_dict: Dict[str, sp.Expr] = OrderedDict()
        self.fCCZ4_RHSs_varname_to_expr_dict["cf_rhs"] = self.cf_rhs
        self.fCCZ4_RHSs_varname_to_expr_dict["Theta_fCCZ4_rhs"] = self.Theta_rhs
        self.fCCZ4_RHSs_varname_to_expr_dict["trK_rhs"] = self.trK_rhs
        for i in range(3):
            self.fCCZ4_RHSs_varname_to_expr_dict[f"lambda_rhsU{i}"] = self.lambda_rhsU[
                i
            ]
            for j in range(i, 3):
                self.fCCZ4_RHSs_varname_to_expr_dict[f"a_rhsDD{i}{j}"] = self.a_rhsDD[
                    i
                ][j]
                self.fCCZ4_RHSs_varname_to_expr_dict[f"h_rhsDD{i}{j}"] = self.h_rhsDD[
                    i
                ][j]
        self.fCCZ4_RHSs_varname_to_expr_dict = OrderedDict(
            sorted(self.fCCZ4_RHSs_varname_to_expr_dict.items())
        )


class FCCZ4RHSsDict(Dict[str, FCCZ4RHSs]):
    """Cache fCCZ4 expressions by coordinate system and construction options."""

    def __init__(self) -> None:
        """Initialize empty expression cache and parameter metadata."""
        super().__init__()
        self._construction_parameters: Dict[str, Tuple[str, bool]] = {}

    def __getitem__(self, CoordSystem_in: str) -> FCCZ4RHSs:
        """
        Return cached expressions matching current global parameters.

        :param CoordSystem_in: Coordinate-system and option cache key.
        :return: Cached or newly built fCCZ4 expressions.

        """
        construction_parameters = (
            par.parval_from_str("EvolvedConformalFactor_cf"),
            par.parval_from_str("detgbarOverdetghat_equals_one"),
        )
        if (
            CoordSystem_in not in self
            or self._construction_parameters.get(CoordSystem_in)
            != construction_parameters
        ):
            CoordSystem = (
                CoordSystem_in.replace("_rfm_precompute", "")
                .replace("_RbarDD_gridfunctions", "")
                .replace("_T4munu", "")
            )
            print(f"Setting up fCCZ4_RHSs[{CoordSystem_in}]...")
            self.__setitem__(
                CoordSystem_in,
                FCCZ4RHSs(
                    CoordSystem=CoordSystem,
                    enable_rfm_precompute="_rfm_precompute" in CoordSystem_in,
                    enable_RbarDD_gridfunctions=(
                        "_RbarDD_gridfunctions" in CoordSystem_in
                    ),
                    enable_T4munu="_T4munu" in CoordSystem_in,
                ),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: FCCZ4RHSs) -> None:
        """
        Store expressions and current construction metadata.

        :param CoordSystem: Cache key.
        :param value: Constructed fCCZ4 expressions.
        """
        dict.__setitem__(self, CoordSystem, value)
        self._construction_parameters[CoordSystem] = (
            par.parval_from_str("EvolvedConformalFactor_cf"),
            par.parval_from_str("detgbarOverdetghat_equals_one"),
        )

    def __delitem__(self, CoordSystem: str) -> None:
        """
        Delete expressions and matching construction metadata.

        :param CoordSystem: Cache key to delete.
        """
        dict.__delitem__(self, CoordSystem)
        self._construction_parameters.pop(CoordSystem, None)

    def clear(self) -> None:
        """Remove all cached expressions and construction metadata."""
        dict.clear(self)
        self._construction_parameters.clear()


fCCZ4_RHSs = FCCZ4RHSsDict()


if __name__ == "__main__":
    import doctest
    import os
    import sys
    from typing import Mapping, cast  # pylint: disable=ungrouped-imports

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    # Reject incompatible storage before constructing caches or installing
    # fCCZ4-owned damping parameters.
    gri.glb_gridfcs_dict.pop("Theta_fCCZ4", None)
    validation_caches: Tuple[Mapping[str, object], ...] = (
        cast(Mapping[str, object], refmetric.reference_metric),
        cast(Mapping[str, object], BSSN_quantities),
        cast(Mapping[str, object], BSSN_RHSs),
        cast(Mapping[str, object], fccz4_constraints.fCCZ4_constraints),
    )
    validation_cache_snapshots = tuple(dict(cache) for cache in validation_caches)
    validation_kappa_snapshots = {
        name: par.glb_code_params_dict.get(name) for name in ("kappa1", "kappa2")
    }
    gri.register_gridfunctions(
        "Theta_fCCZ4",
        group="AUX",
        gf_array_name="diagnostic_output_gfs",
        is_basename=False,
    )
    try:
        FCCZ4RHSs()
    except ValueError as exc:
        if str(exc) != "Theta_fCCZ4 must be an EVOL gridfunction stored in in_gfs.":
            raise
    else:
        raise AssertionError("malformed Theta_fCCZ4 storage must be rejected")
    for cache, snapshot in zip(validation_caches, validation_cache_snapshots):
        if cache.keys() != snapshot.keys() or any(
            cache[key] is not value for key, value in snapshot.items()
        ):
            raise AssertionError(
                "rejected fCCZ4 construction mutated construction caches"
            )
    if any(
        par.glb_code_params_dict.get(name) is not value
        for name, value in validation_kappa_snapshots.items()
    ):
        raise AssertionError("rejected fCCZ4 construction changed damping parameters")
    gri.glb_gridfcs_dict.pop("Theta_fCCZ4", None)

    # BHaHAHA owns an unrelated AUX diagnostic named Theta. Verify both
    # registration orders retain distinct storage identities.
    gri.glb_gridfcs_dict.pop("Theta", None)
    gri.glb_gridfcs_dict.pop("Theta_fCCZ4", None)
    gri.register_gridfunctions(
        "Theta", group="AUX", gf_array_name="diagnostic_output_gfs"
    )
    diagnostic_first_rhs = FCCZ4RHSs()
    if (
        diagnostic_first_rhs.Theta.name != "Theta_fCCZ4"
        or gri.glb_gridfcs_dict["Theta_fCCZ4"].group != "EVOL"
        or gri.glb_gridfcs_dict["Theta"].group != "AUX"
    ):
        raise AssertionError("BHaHAHA-first Theta registration changed fCCZ4 storage")
    gri.glb_gridfcs_dict.pop("Theta", None)
    gri.glb_gridfcs_dict.pop("Theta_fCCZ4", None)
    fccz4_first_rhs = FCCZ4RHSs()
    gri.register_gridfunctions(
        "Theta", group="AUX", gf_array_name="diagnostic_output_gfs"
    )
    if (
        fccz4_first_rhs.Theta.name != "Theta_fCCZ4"
        or gri.glb_gridfcs_dict["Theta_fCCZ4"].group != "EVOL"
        or gri.glb_gridfcs_dict["Theta"].group != "AUX"
    ):
        raise AssertionError("fCCZ4-first Theta registration changed BHaHAHA storage")
    gri.glb_gridfcs_dict.pop("Theta_fCCZ4", None)

    cases = (
        ("Cartesian", False, False, False, "Cartesian"),
        (
            "Cartesian",
            False,
            True,
            False,
            "Cartesian_RbarDD_gridfunctions",
        ),
        ("SinhSpherical", False, False, False, "SinhSpherical"),
        (
            "SinhSpherical",
            False,
            True,
            False,
            "SinhSpherical_RbarDD_gridfunctions",
        ),
        (
            "SinhSpherical",
            True,
            False,
            True,
            "SinhSpherical_rfm_precompute_T4munu",
        ),
    )
    for (
        case_coord,
        case_enable_rfm_precompute,
        case_enable_RbarDD_gridfunctions,
        case_enable_T4munu,
        trusted_suffix,
    ) in cases:
        rhs = FCCZ4RHSs(
            CoordSystem=case_coord,
            enable_rfm_precompute=case_enable_rfm_precompute,
            enable_RbarDD_gridfunctions=case_enable_RbarDD_gridfunctions,
            enable_T4munu=case_enable_T4munu,
        )
        validation_suffix = (
            case_coord
            + ("_rfm_precompute" if case_enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if case_enable_RbarDD_gridfunctions else "")
        )
        canonical = fccz4_constraints.fCCZ4_constraints[
            validation_suffix + ("_T4munu" if case_enable_T4munu else "")
        ]
        validation_Bq = BSSN_quantities[validation_suffix]
        validation_rfm = refmetric.reference_metric[
            case_coord + ("_rfm_precompute" if case_enable_rfm_precompute else "")
        ]

        # Rebuild the literal Mewes kappa3=1 shift sector independently,
        # retaining the off-constraint vector C^i = Z4constraintU^i.
        validation_Dhat_div_beta = sp.sympify(0)
        validation_Dhat_betaUD = ixp.zerorank2()
        validation_alpha_dD = ixp.declarerank1("alpha_dD")
        for k in range(3):
            validation_Dhat_div_beta += validation_Bq.betaU_dD[k][k]
            for ell in range(3):
                validation_Dhat_div_beta += (
                    validation_rfm.GammahatUDD[k][ell][k] * validation_Bq.betaU[ell]
                )
        for i in range(3):
            for k in range(3):
                validation_Dhat_betaUD[i][k] = validation_Bq.betaU_dD[i][k]
                for ell in range(3):
                    validation_Dhat_betaUD[i][k] += (
                        validation_rfm.GammahatUDD[i][ell][k] * validation_Bq.betaU[ell]
                    )
        expected_delta = ixp.zerorank1()
        expected_mewes_shift_delta = ixp.zerorank1()
        expected_alic_shift_delta = ixp.zerorank1()
        expected_mewes_minus_alic = ixp.zerorank1()
        for i in range(3):
            for j in range(3):
                expected_delta[i] += (
                    2
                    * validation_Bq.gammabarUU[i][j]
                    * (
                        validation_Bq.alpha * rhs.Theta_dD[j]
                        - rhs.Theta * validation_alpha_dD[j]
                    )
                )
            expected_delta[i] += (
                -sp.Rational(2, 3)
                * validation_Bq.alpha
                * validation_Bq.trK
                * rhs.Z4constraintU[i]
                - rhs.kappa1 * rhs.Z4constraintU[i]
            )
            # Relative to the reference-metric BSSN shift sector, literal
            # Mewes kappa3=1 contributes 2 C^i D/3 - C^k D_k beta^i.
            expected_mewes_shift_delta[i] += (
                sp.Rational(2, 3) * rhs.Z4constraintU[i] * validation_Dhat_div_beta
            )
            for k in range(3):
                shift_gradient_term = (
                    -rhs.Z4constraintU[k] * validation_Dhat_betaUD[i][k]
                )
                expected_mewes_shift_delta[i] += shift_gradient_term
            expected_delta[i] += expected_mewes_shift_delta[i]

            # Alic Eq. (19) is Cartesian. Its kappa3-independent shift terms
            # use the geometric contracted connection, so at kappa3=1 its
            # delta relative to the same BSSN base is only 2 C^i D/3.
            if case_coord == "Cartesian":
                expected_alic_shift_delta[i] = (
                    sp.Rational(2, 3) * rhs.Z4constraintU[i] * validation_Dhat_div_beta
                )
                for k in range(3):
                    expected_mewes_minus_alic[i] += (
                        -rhs.Z4constraintU[k] * validation_Dhat_betaUD[i][k]
                    )

        actual_reconstructions: Dict[str, ve.ExpressionValue] = {
            "Mewes_delta": rhs.Lambdatilde_rhsU_delta
        }
        expected_reconstructions: Dict[str, ve.ExpressionValue] = {
            "Mewes_delta": expected_delta
        }
        if case_coord == "Cartesian":
            actual_reconstructions["Mewes_minus_Alic"] = [
                expected_mewes_shift_delta[i] - expected_alic_shift_delta[i]
                for i in range(3)
            ]
            expected_reconstructions["Mewes_minus_Alic"] = expected_mewes_minus_alic
        ve.assert_equal(
            actual_reconstructions,
            expected_reconstructions,
            suppress_message=True,
        )

        # Cartesian coefficients distinguish the Mewes and Alic/Sanchis-Gual
        # off-constraint equations without a same-constructor snapshot.
        if case_coord == "Cartesian":
            shift_gradient_coefficients = []
            expected_shift_gradient_coefficients = []
            for i in range(3):
                for a in range(3):
                    for b in range(3):
                        shift_gradient_coefficient = sp.diff(
                            rhs.Lambdatilde_rhsU_delta[i],
                            validation_Bq.betaU_dD[a][b],
                        )
                        expected_shift_gradient_coefficient = sp.Rational(
                            2, 3
                        ) * rhs.Z4constraintU[i] * sp.KroneckerDelta(
                            a, b
                        ) - rhs.Z4constraintU[
                            b
                        ] * sp.KroneckerDelta(
                            i, a
                        )
                        shift_gradient_coefficients.append(shift_gradient_coefficient)
                        expected_shift_gradient_coefficients.append(
                            expected_shift_gradient_coefficient
                        )
            ve.assert_equal(
                {"shift_gradient_coefficients": shift_gradient_coefficients},
                {"shift_gradient_coefficients": (expected_shift_gradient_coefficients)},
                suppress_message=True,
            )
            aggregate_shift_gradient_coefficients = []
            expected_aggregate_shift_gradient_coefficients = []
            for i in range(3):
                for a in range(3):
                    for b in range(3):
                        aggregate_shift_gradient_coefficient = sp.diff(
                            rhs.Lambdatilde_rhsU[i],
                            validation_Bq.betaU_dD[a][b],
                        )
                        expected_aggregate_shift_gradient_coefficient = -(
                            rhs.LambdatildeU[b] + rhs.Z4constraintU[b]
                        ) * sp.KroneckerDelta(i, a) + sp.Rational(
                            2, 3
                        ) * rhs.LambdatildeU[
                            i
                        ] * sp.KroneckerDelta(
                            a, b
                        )
                        aggregate_shift_gradient_coefficients.append(
                            aggregate_shift_gradient_coefficient
                        )
                        expected_aggregate_shift_gradient_coefficients.append(
                            expected_aggregate_shift_gradient_coefficient
                        )
            ve.assert_equal(
                {
                    "aggregate_shift_gradient_coefficients": (
                        aggregate_shift_gradient_coefficients
                    )
                },
                {
                    "aggregate_shift_gradient_coefficients": (
                        expected_aggregate_shift_gradient_coefficients
                    )
                },
                suppress_message=True,
            )

        aggregate_Theta_rhs = rhs.Theta_rhs
        aggregate_cf_rhs = rhs.cf_rhs
        aggregate_trK_rhs = rhs.trK_rhs
        expression_dict: Dict[str, ve.ExpressionValue] = dict(
            rhs.fCCZ4_RHSs_varname_to_expr_dict
        )
        expression_dict.update(
            {
                "Theta": rhs.Theta,
                "Theta_dD": rhs.Theta_dD,
                "Theta_dupD": rhs.Theta_dupD,
                "H_Z4": rhs.H_Z4,
                "LambdatildeU": rhs.LambdatildeU,
                "LambdatildeU_dupD": rhs.LambdatildeU_dupD,
                "Lambdatilde_rhsU": rhs.Lambdatilde_rhsU,
                "Lambdatilde_rhsU_delta": rhs.Lambdatilde_rhsU_delta,
                "RbarZ4": rhs.RbarZ4,
                "RbarZ4DD": rhs.RbarZ4DD,
                "RbarZ4DD_delta": rhs.RbarZ4DD_delta,
                "Z4constraintU": rhs.Z4constraintU,
                "ZD": rhs.ZD,
                "ZU": rhs.ZU,
                "ZbarU": rhs.ZbarU,
                "aggregate_Theta_rhs": aggregate_Theta_rhs,
                "aggregate_a_rhsDD": rhs.a_rhsDD,
                "aggregate_cf_rhs": aggregate_cf_rhs,
                "aggregate_h_rhsDD": rhs.h_rhsDD,
                "aggregate_lambda_rhsU": rhs.lambda_rhsU,
                "aggregate_trK_rhs": aggregate_trK_rhs,
                "kappa1": rhs.kappa1,
                "kappa2": rhs.kappa2,
                "canonical_H_Z4_residual": rhs.H_Z4 - canonical.H_Z4,
                "canonical_LambdatildeU_residual": [
                    rhs.LambdatildeU[i] - canonical.LambdatildeU[i] for i in range(3)
                ],
                "canonical_RbarZ4DD_delta_residual": [
                    [
                        rhs.RbarZ4DD_delta[i][j] - canonical.RbarZ4DD_delta[i][j]
                        for j in range(3)
                    ]
                    for i in range(3)
                ],
                "canonical_RbarZ4DD_residual": [
                    [rhs.RbarZ4DD[i][j] - canonical.RbarZ4DD[i][j] for j in range(3)]
                    for i in range(3)
                ],
                "canonical_RbarZ4_residual": rhs.RbarZ4 - canonical.RbarZ4,
                "canonical_Z4constraintU_residual": [
                    rhs.Z4constraintU[i] - canonical.Z4constraintU[i] for i in range(3)
                ],
                "canonical_ZD_residual": [
                    rhs.ZD[i] - canonical.ZD[i] for i in range(3)
                ],
                "canonical_ZU_residual": [
                    rhs.ZU[i] - canonical.ZU[i] for i in range(3)
                ],
                "canonical_ZbarU_residual": [
                    rhs.ZbarU[i] - canonical.ZbarU[i] for i in range(3)
                ],
                "canonical_copy_residual": sp.Integer(
                    0
                    if rhs.LambdatildeU is not canonical.LambdatildeU
                    and rhs.Z4constraintU is not canonical.Z4constraintU
                    and rhs.ZbarU is not canonical.ZbarU
                    and rhs.ZU is not canonical.ZU
                    and rhs.ZD is not canonical.ZD
                    and rhs.RbarZ4DD_delta is not canonical.RbarZ4DD_delta
                    and rhs.RbarZ4DD is not canonical.RbarZ4DD
                    else 1
                ),
                "lambda_rescaling_residual": [
                    rhs.lambda_rhsU[i] - rhs.Lambdatilde_rhsU[i] / validation_rfm.ReU[i]
                    for i in range(3)
                ],
                "mapped_output_count_residual": sp.Integer(
                    len(rhs.fCCZ4_RHSs_varname_to_expr_dict) - 18
                ),
                "mapped_output_order_residual": sp.Integer(
                    0
                    if list(rhs.fCCZ4_RHSs_varname_to_expr_dict)
                    == sorted(rhs.fCCZ4_RHSs_varname_to_expr_dict)
                    else 1
                ),
                "RbarZ4_trace_residual": rhs.RbarZ4
                - sum(
                    validation_Bq.gammabarUU[i][j] * rhs.RbarZ4DD[i][j]
                    for i in range(3)
                    for j in range(3)
                ),
                "Z4constraintU_definition_residual": [
                    rhs.Z4constraintU[i]
                    - (rhs.LambdatildeU[i] - validation_Bq.DGammaU[i])
                    for i in range(3)
                ],
                "ZD_definition_residual": [
                    rhs.ZD[i]
                    - sp.Rational(1, 2)
                    * sum(
                        validation_Bq.gammabarDD[i][j] * rhs.Z4constraintU[j]
                        for j in range(3)
                    )
                    for i in range(3)
                ],
                "ZU_definition_residual": [
                    rhs.ZU[i] - validation_Bq.exp_m4phi * rhs.ZbarU[i] for i in range(3)
                ],
                "ZbarU_definition_residual": [
                    rhs.ZbarU[i] - sp.Rational(1, 2) * rhs.Z4constraintU[i]
                    for i in range(3)
                ],
            }
        )
        cache_key = validation_suffix + ("_T4munu" if case_enable_T4munu else "")
        cached_rhs = fCCZ4_RHSs[cache_key]
        expression_dict["cache_reuse_residual"] = sp.Integer(
            0 if fCCZ4_RHSs[cache_key] is cached_rhs else 1
        )
        processed = ve.process_dictionary_of_expressions(
            expression_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{trusted_suffix}",
            processed,
        )

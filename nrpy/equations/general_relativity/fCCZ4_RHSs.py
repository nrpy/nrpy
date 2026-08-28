"""
Construct fully covariant conformal Z4 (fCCZ4) evolution equations.

This module implements the three-dimensional reference-metric formulation of
Mewes et al. (2020), with the Lagrangian conformal-determinant condition and
``kappa1`` damping terms that contain no additional lapse factor.

Mewes et al. appear to contain a convention error in the shift sector.  They
define ``partial_0 = partial_t - L_beta``, whose full vector Lie derivative
stretches the evolved connection by
``-Lambdatilde^k Dhat_k beta^i``.  Since
``Lambdatilde^i = DeltaGamma^i + C^i``, that stretch already contains
``-C^k Dhat_k beta^i``, but the printed ``kappa3=1`` bracket adds this term
again.  NRPy corrects that apparent paper error, not a discretionary
formulation departure.  It retains ``+(2/3) C^i Dhat_k beta^k`` because the
reused BSSN base divergence contains geometric ``DeltaGamma`` and needs
promotion to evolved ``Lambdatilde``.

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

        # Brhs already applies the full vector Lie derivative to its evolved
        # connection slot, which stores LambdatildeU here.  Mewes et al.
        # define partial_0 = partial_t - L_beta but then print an additional
        # -C^k Dhat_k beta^i in the kappa3=1 bracket, duplicating that stretch
        # for C^i = LambdatildeU^i - DeltaGamma^i.  Omit the apparent paper
        # error.  Keep +2 C^i Dhat_k beta^k/3 because the BSSN base divergence
        # contains geometric DeltaGamma and needs promotion to LambdatildeU.
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

    import nrpy.validate_expressions.validate_expressions as ve

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    cases = (
        ("SinhCartesian", False, False, False, "SinhCartesian"),
        (
            "SinhCartesian",
            False,
            True,
            False,
            "SinhCartesian_RbarDD_gridfunctions",
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
        ("Cartesian", False, False, False, "Cartesian"),
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
        processed = ve.process_dictionary_of_expressions(
            rhs.fCCZ4_RHSs_varname_to_expr_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{trusted_suffix}",
            processed,
        )

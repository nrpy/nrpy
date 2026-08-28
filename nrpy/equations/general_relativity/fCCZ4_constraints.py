"""
Construct fully covariant conformal Z4 (fCCZ4) constraint expressions.

This module owns the spatial Z4 connection constraint, its associated
conformal Ricci tensor, and the fCCZ4 Hamiltonian expression.  The Hamiltonian
excludes ``-2 Theta K``; evolution modules append that term where required.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, Tuple

import sympy as sp

import nrpy.grid as gri
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.equations.general_relativity import T4munu
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities


class FCCZ4Constraints:
    """
    Construct and store fCCZ4 connection and Hamiltonian constraints.

    Instances expose the conformal connection constraint, three spatial Z4
    vectors, the Z4-corrected conformal Ricci tensor and trace, and ``H_Z4``.
    Construction reuses parameter-aware BSSN quantities without mutating them.
    """

    def __init__(
        self,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        enable_RbarDD_gridfunctions: bool = False,
        enable_T4munu: bool = False,
    ) -> None:
        r"""
        Build fCCZ4 constraint expressions for one metric configuration.

        ``Z4constraintU`` is
        :math:`C^i=\widetilde{\Lambda}^i-\Delta\Gamma^i`.  ``RbarZ4DD``
        uses the evolved connection in the linear Ricci slot and adds the
        derivative-free spatial-Z4 correction.  ``H_Z4`` includes the
        established stress-energy Hamiltonian source only when requested.

        :param CoordSystem: Reference-metric coordinate system.
        :param enable_rfm_precompute: Use precomputed reference-metric symbols.
        :param enable_RbarDD_gridfunctions: Read the BSSN-shaped conformal
            Ricci components, with ``lambdaU`` interpreted as
            ``LambdatildeU``, from auxiliary gridfunctions.
        :param enable_T4munu: Include the Hamiltonian stress-energy source.

        """
        suffix = (
            CoordSystem
            + ("_rfm_precompute" if enable_rfm_precompute else "")
            + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        )
        Bq = BSSN_quantities[suffix]

        # The evolved BSSN connection slot stores tilde Lambda in fCCZ4.
        self.LambdatildeU = [Bq.LambdabarU[i] for i in range(3)]
        self.Z4constraintU = ixp.zerorank1()
        Z4constraintD = ixp.zerorank1()
        self.ZbarU = ixp.zerorank1()
        self.ZU = ixp.zerorank1()
        self.ZD = ixp.zerorank1()
        for i in range(3):
            self.Z4constraintU[i] = self.LambdatildeU[i] - Bq.DGammaU[i]
            self.ZbarU[i] = sp.Rational(1, 2) * self.Z4constraintU[i]
            self.ZU[i] = Bq.exp_m4phi * self.ZbarU[i]
        for i in range(3):
            for j in range(3):
                Z4constraintD[i] += Bq.gammabarDD[i][j] * self.Z4constraintU[j]
            self.ZD[i] = sp.Rational(1, 2) * Z4constraintD[i]

        # Bq.RbarDD uses tilde Lambda in its linear connection-derivative slot
        # and geometric DGammaU/DGammaUDD in its nonlinear terms.  Add the
        # remaining derivative-free Z contribution from Mewes et al. Eq. (32).
        self.RbarZ4DD_delta = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                self.RbarZ4DD_delta[i][j] = -2 * (
                    Z4constraintD[i] * Bq.phi_dD[j] + Z4constraintD[j] * Bq.phi_dD[i]
                )
                for ell in range(3):
                    self.RbarZ4DD_delta[i][j] += (
                        2
                        * Bq.gammabarDD[i][j]
                        * self.Z4constraintU[ell]
                        * Bq.phi_dD[ell]
                    )
                    for k in range(3):
                        self.RbarZ4DD_delta[i][j] += (
                            sp.Rational(1, 2)
                            * self.Z4constraintU[ell]
                            * (
                                Bq.gammabarDD[k][i] * Bq.DGammaUDD[k][j][ell]
                                + Bq.gammabarDD[k][j] * Bq.DGammaUDD[k][i][ell]
                            )
                        )

        self.RbarZ4DD = ixp.zerorank2()
        self.RbarZ4 = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                self.RbarZ4DD[i][j] = Bq.RbarDD[i][j] + self.RbarZ4DD_delta[i][j]
                self.RbarZ4 += Bq.gammabarUU[i][j] * self.RbarZ4DD[i][j]

        gradphi_squared = sp.sympify(0)
        laplacianphi = sp.sympify(0)
        Abar_squared = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                gradphi_squared += (
                    Bq.gammabarUU[i][j] * Bq.phi_dBarD[i] * Bq.phi_dBarD[j]
                )
                laplacianphi += Bq.gammabarUU[i][j] * Bq.phi_dBarDD[i][j]
                Abar_squared += Bq.AbarDD[i][j] * Bq.AbarUU[i][j]
        self.H_Z4 = (
            Bq.exp_m4phi * (self.RbarZ4 - 8 * gradphi_squared - 8 * laplacianphi)
            + sp.Rational(2, 3) * Bq.trK**2
            - Abar_squared
        )
        if enable_T4munu:
            if "T4UU00" not in gri.glb_gridfcs_dict:
                _ = gri.register_gridfunctions_for_single_rank2(
                    "T4UU",
                    symmetry="sym01",
                    dimension=4,
                    group="AUXEVOL",
                )
            source_H, _source_MU = T4munu.BSSN_constraints_T4UU_source_terms(
                CoordSystem=CoordSystem,
                enable_rfm_precompute=enable_rfm_precompute,
            )
            self.H_Z4 += source_H


class FCCZ4ConstraintsDict(Dict[str, FCCZ4Constraints]):
    """
    Cache fCCZ4 constraints by coordinate system and construction options.

    Entries rebuild when the conformal-factor or conformal-determinant
    parameters change.
    """

    def __init__(self) -> None:
        """Initialize an empty constraint cache and parameter metadata."""
        super().__init__()
        self._construction_parameters: Dict[str, Tuple[str, bool]] = {}

    def __getitem__(self, CoordSystem_in: str) -> FCCZ4Constraints:
        """
        Return constraints matching current expression-affecting parameters.

        :param CoordSystem_in: Coordinate-system and option cache key.
        :return: Cached or newly built fCCZ4 constraints.

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
            print(f"Setting up fCCZ4_constraints[{CoordSystem_in}]...")
            self.__setitem__(
                CoordSystem_in,
                FCCZ4Constraints(
                    CoordSystem=CoordSystem,
                    enable_rfm_precompute="_rfm_precompute" in CoordSystem_in,
                    enable_RbarDD_gridfunctions=(
                        "_RbarDD_gridfunctions" in CoordSystem_in
                    ),
                    enable_T4munu="_T4munu" in CoordSystem_in,
                ),
            )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: FCCZ4Constraints) -> None:
        """
        Store constraints with current construction metadata.

        :param CoordSystem: Cache key.
        :param value: Constructed fCCZ4 constraints.
        """
        dict.__setitem__(self, CoordSystem, value)
        self._construction_parameters[CoordSystem] = (
            par.parval_from_str("EvolvedConformalFactor_cf"),
            par.parval_from_str("detgbarOverdetghat_equals_one"),
        )

    def __delitem__(self, CoordSystem: str) -> None:
        """
        Delete constraints and matching construction metadata.

        :param CoordSystem: Cache key to delete.
        """
        dict.__delitem__(self, CoordSystem)
        self._construction_parameters.pop(CoordSystem, None)

    def clear(self) -> None:
        """Remove all cached constraints and construction metadata."""
        dict.clear(self)
        self._construction_parameters.clear()


fCCZ4_constraints = FCCZ4ConstraintsDict()


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

    original_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    cases = (
        ("SinhCartesian", False, False, False, "W", "SinhCartesian"),
        ("SinhCartesian", False, False, False, "phi", "SinhCartesian_phi"),
        ("SinhCartesian", False, False, False, "chi", "SinhCartesian_chi"),
        (
            "SinhCartesian",
            False,
            True,
            False,
            "W",
            "SinhCartesian_RbarDD_gridfunctions",
        ),
        ("SinhSpherical", False, False, False, "W", "SinhSpherical"),
        (
            "SinhSpherical",
            False,
            False,
            False,
            "phi",
            "SinhSpherical_phi",
        ),
        (
            "SinhSpherical",
            False,
            False,
            False,
            "chi",
            "SinhSpherical_chi",
        ),
        (
            "SinhSpherical",
            False,
            True,
            False,
            "W",
            "SinhSpherical_RbarDD_gridfunctions",
        ),
        (
            "SinhSpherical",
            True,
            False,
            True,
            "W",
            "SinhSpherical_rfm_precompute_T4munu",
        ),
        ("Cartesian", False, False, False, "W", "Cartesian"),
    )
    for (
        case_coord,
        case_enable_rfm_precompute,
        case_enable_RbarDD_gridfunctions,
        case_enable_T4munu,
        case_conformal_factor,
        trusted_suffix,
    ) in cases:
        par.set_parval_from_str("EvolvedConformalFactor_cf", case_conformal_factor)
        try:
            constraints = FCCZ4Constraints(
                CoordSystem=case_coord,
                enable_rfm_precompute=case_enable_rfm_precompute,
                enable_RbarDD_gridfunctions=case_enable_RbarDD_gridfunctions,
                enable_T4munu=case_enable_T4munu,
            )
        finally:
            par.set_parval_from_str("EvolvedConformalFactor_cf", original_cf)
        processed = ve.process_dictionary_of_expressions(
            {"H_Z4": constraints.H_Z4}, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{trusted_suffix}",
            processed,
        )

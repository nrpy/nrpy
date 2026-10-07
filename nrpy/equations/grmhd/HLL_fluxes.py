"""
Construct GRMHD fluid fluxes from the HLL approximate Riemann solver at cell interfaces.

Magnetic stress-energy follows Duez et al., Phys. Rev. D 72, 024028
(2005), Eqs. (32)-(33); its HLL flux is Eq. (48),
https://arxiv.org/abs/astro-ph/0503420v2.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List, Tuple, Union

import sympy as sp

import nrpy.indexedexp as ixp
from nrpy.equations.grhd.HLL_fluxes import (
    HLL_solver,
    calculate_Tmunu_and_contractions_from_equations,
)
from nrpy.equations.grmhd.characteristic_speeds import find_cmax_cmin
from nrpy.equations.grmhd.GRMHD_equations import GRMHD_Equations


def calculate_HLL_fluxes(
    flux_dirn: int,
    alpha_face: sp.Expr,
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    e6phi_face: sp.Expr,
    u4rU: List[sp.Expr],
    u4lU: List[sp.Expr],
    BmagrU: List[sp.Expr],
    BmaglU: List[sp.Expr],
    rho_b_r: sp.Expr,
    rho_b_l: sp.Expr,
    Ye_r: sp.Expr,
    Ye_l: sp.Expr,
    S_r: sp.Expr,
    S_l: sp.Expr,
    P_r: sp.Expr,
    P_l: sp.Expr,
    h_r: sp.Expr,
    h_l: sp.Expr,
    cs2_r: sp.Expr,
    cs2_l: sp.Expr,
) -> Tuple[sp.Expr, sp.Expr, sp.Expr, sp.Expr, List[sp.Expr]]:
    """
    Calculate symbolic HLL fluxes for the GRMHD fluid evolution system.

    The HLL combination is Duez et al. (2005), Eq. (48); see HLL_solver.
    For each face state, a Cartesian GRMHD_Equations object supplies the
    magnetic stress-energy to the GRHD contraction function. All face inputs,
    including BmagrU and BmaglU, must use the basis of gamma_faceDD; for
    reference-metric evolutions this is the rescaled basis, and the returned
    fluxes are rescaled fluxes.

    :param flux_dirn: Flux direction.
    :param alpha_face: Lapse on the cell face.
    :param gamma_faceDD: Spatial metric on the cell face.
    :param beta_faceU: Shift vector on the cell face.
    :param e6phi_face: Face reference-metric volume factor e^(6 phi), equal to sqrt(gamma/gammahat) when det(gammabar) = det(gammahat); for rescaled face data it equals sqrt(det gamma_faceDD).
    :param u4rU: Four-velocity reconstructed to the right side.
    :param u4lU: Four-velocity reconstructed to the left side.
    :param BmagrU: Eulerian magnetic field scaled by 1/sqrt(4 pi) on the right side, in the basis of gamma_faceDD.
    :param BmaglU: Eulerian magnetic field scaled by 1/sqrt(4 pi) on the left side, in the basis of gamma_faceDD.
    :param rho_b_r: Density on the right side.
    :param rho_b_l: Density on the left side.
    :param Ye_r: Electron fraction on the right side.
    :param Ye_l: Electron fraction on the left side.
    :param S_r: Primitive entropy variable on the right side.
    :param S_l: Primitive entropy variable on the left side.
    :param P_r: Pressure on the right side.
    :param P_l: Pressure on the left side.
    :param h_r: Specific enthalpy on the right side.
    :param h_l: Specific enthalpy on the left side.
    :param cs2_r: Sound speed squared on the right side.
    :param cs2_l: Sound speed squared on the left side.
    :return: HLL fluxes for density, electron fraction, entropy, energy, and
        momentum.
    """
    # Step 1: Compute the conserved variables and physical fluxes on each side.
    grmhd_eqs_r = GRMHD_Equations(CoordSystem="Cartesian", enable_rfm_precompute=False)
    grmhd_eqs_r.BmagU = BmagrU.copy()
    (
        U_rho_star_r,
        F_rho_star_r,
        U_Ye_star_r,
        F_Ye_star_r,
        U_S_star_r,
        F_S_star_r,
        U_tau_tilde_r,
        F_tau_tilde_r,
        U_S_tilde_rD,
        F_S_tilde_rD,
    ) = calculate_Tmunu_and_contractions_from_equations(
        grmhd_eqs_r,
        flux_dirn,
        gamma_faceDD,
        beta_faceU,
        alpha_face,
        e6phi_face,
        rho_b_r,
        Ye_r,
        S_r,
        P_r,
        h_r,
        u4rU,
    )

    grmhd_eqs_l = GRMHD_Equations(CoordSystem="Cartesian", enable_rfm_precompute=False)
    grmhd_eqs_l.BmagU = BmaglU.copy()
    (
        U_rho_star_l,
        F_rho_star_l,
        U_Ye_star_l,
        F_Ye_star_l,
        U_S_star_l,
        F_S_star_l,
        U_tau_tilde_l,
        F_tau_tilde_l,
        U_S_tilde_lD,
        F_S_tilde_lD,
    ) = calculate_Tmunu_and_contractions_from_equations(
        grmhd_eqs_l,
        flux_dirn,
        gamma_faceDD,
        beta_faceU,
        alpha_face,
        e6phi_face,
        rho_b_l,
        Ye_l,
        S_l,
        P_l,
        h_l,
        u4lU,
    )

    # Step 2: Compute the fastest left- and right-going signal speeds.
    cmin, cmax = find_cmax_cmin(
        flux_dirn,
        gamma_faceDD,
        beta_faceU,
        alpha_face,
        u4rU,
        u4lU,
        BmagrU,
        BmaglU,
        rho_b_r,
        rho_b_l,
        h_r,
        h_l,
        cs2_r,
        cs2_l,
    )

    # Step 3: Assemble the HLL fluxes.
    rho_star_HLL_flux = HLL_solver(
        cmax, cmin, F_rho_star_r, F_rho_star_l, U_rho_star_r, U_rho_star_l
    )

    Ye_star_HLL_flux = HLL_solver(
        cmax, cmin, F_Ye_star_r, F_Ye_star_l, U_Ye_star_r, U_Ye_star_l
    )

    S_star_HLL_flux = HLL_solver(
        cmax, cmin, F_S_star_r, F_S_star_l, U_S_star_r, U_S_star_l
    )

    tau_tilde_HLL_flux = HLL_solver(
        cmax, cmin, F_tau_tilde_r, F_tau_tilde_l, U_tau_tilde_r, U_tau_tilde_l
    )

    Stilde_flux_HLLD = ixp.zerorank1()
    for mom_comp in range(3):
        Stilde_flux_HLLD[mom_comp] = HLL_solver(
            cmax,
            cmin,
            F_S_tilde_rD[mom_comp],
            F_S_tilde_lD[mom_comp],
            U_S_tilde_rD[mom_comp],
            U_S_tilde_lD[mom_comp],
        )

    return (
        rho_star_HLL_flux,
        Ye_star_HLL_flux,
        S_star_HLL_flux,
        tau_tilde_HLL_flux,
        Stilde_flux_HLLD,
    )


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.reference_metric as refmetric
    import nrpy.validate_expressions.validate_expressions as ve
    from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
    from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
    from nrpy.equations.grhd.HLL_fluxes import (
        calculate_HLL_fluxes as calculate_grhd_HLL_fluxes,
    )
    from nrpy.equations.grmhd.characteristic_speeds import _nrpyAbs_to_Abs

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    rfm = refmetric.reference_metric["Cartesian"]

    alpha_face_test = sp.symbols("alpha_face", real=True)
    cf_face = sp.symbols("cf_face", real=True)
    h_faceDD = ixp.declarerank2("h_faceDD", symmetry="sym01", dimension=3)
    vet_faceU = ixp.declarerank1("vet_faceU", dimension=3)

    # Step 1: Define symbolic right and left interface states.
    rescaledvrU = ixp.declarerank1("rescaledvrU", dimension=3)
    rescaledvlU = ixp.declarerank1("rescaledvlU", dimension=3)
    BmagrU_test = ixp.declarerank1("BmagrU", dimension=3)
    BmaglU_test = ixp.declarerank1("BmaglU", dimension=3)

    VrU = ixp.zerorank1()
    VlU = ixp.zerorank1()

    u4rUt_test = sp.symbols("u4rUt", real=True)
    u4lUt_test = sp.symbols("u4lUt", real=True)
    u4rU_test = ixp.zerorank1(dimension=4)
    u4rU_test[0] = u4rUt_test
    u4lU_test = ixp.zerorank1(dimension=4)
    u4lU_test[0] = u4lUt_test

    beta_faceU_test = ixp.zerorank1()
    for i in range(3):
        VrU[i] = rescaledvrU[i] * rfm.ReU[i]
        VlU[i] = rescaledvlU[i] * rfm.ReU[i]
        u4rU_test[i + 1] = VrU[i] * u4rU_test[0]
        u4lU_test[i + 1] = VlU[i] * u4lU_test[0]
        beta_faceU_test[i] = vet_faceU[i] * rfm.ReU[i]

    rho_b_r_test = sp.symbols("rhob_r", real=True)
    rho_b_l_test = sp.symbols("rhob_l", real=True)

    Ye_r_test = sp.symbols("Ye_r", real=True)
    Ye_l_test = sp.symbols("Ye_l", real=True)

    S_r_test = sp.symbols("S_r", real=True)
    S_l_test = sp.symbols("S_l", real=True)

    P_r_test = sp.symbols("P_r", real=True)
    P_l_test = sp.symbols("P_l", real=True)

    h_r_test = sp.symbols("h_r", real=True)
    h_l_test = sp.symbols("h_l", real=True)

    cs2_r_test = sp.symbols("cs2_r", real=True)
    cs2_l_test = sp.symbols("cs2_l", real=True)

    # Step 2: Build face-centered metric quantities in Cartesian form.
    AitoB = BSSN_to_ADM(CoordSystem="Cartesian")
    Bq = BSSN_quantities["Cartesian"]

    e6phi_face_test = (Bq.exp_m4phi ** sp.Rational(-3, 2)).subs(Bq.cf, cf_face)

    gamma_faceDD_test = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gamma_faceDD_test[i][j] = (
                AitoB.gammaDD[i][j]
                .subs(Bq.hDD[i][j], h_faceDD[i][j])
                .subs(Bq.cf, cf_face)
            )

    # Step 3: Evaluate conserved variables and fluxes symbolically.
    grmhd_eqs_test = GRMHD_Equations(
        CoordSystem="Cartesian", enable_rfm_precompute=False
    )
    grmhd_eqs_test.BmagU = BmagrU_test.copy()
    (
        U_rho_star_test,
        F_rho_star_test,
        U_Ye_star_test,
        F_Ye_star_test,
        U_S_star_test,
        F_S_star_test,
        U_tau_tilde_test,
        F_tau_tilde_test,
        U_S_tildeD_test,
        F_S_tildeD_test,
    ) = calculate_Tmunu_and_contractions_from_equations(
        grmhd_eqs_test,
        2,
        gamma_faceDD_test,
        beta_faceU_test,
        alpha_face_test,
        e6phi_face_test,
        rho_b_r_test,
        Ye_r_test,
        S_r_test,
        P_r_test,
        h_r_test,
        u4rU_test,
    )

    cmin_test, cmax_test = find_cmax_cmin(
        1,
        gamma_faceDD_test,
        beta_faceU_test,
        alpha_face_test,
        u4rU_test,
        u4lU_test,
        BmagrU_test,
        BmaglU_test,
        rho_b_r_test,
        rho_b_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    )

    HLL_test = HLL_solver(
        _nrpyAbs_to_Abs(cmax_test),
        _nrpyAbs_to_Abs(cmin_test),
        F_rho_star_test,
        F_Ye_star_test,
        U_rho_star_test,
        U_Ye_star_test,
    )

    # Step 4: Evaluate the magnetic HLL flux expressions. In flux direction 1,
    #         both speed bounds are nonzero at the trusted sample values, so
    #         both face states and the magnetic signal speeds enter the fluxes.
    (
        rho_star_HLL_flux_test,
        Ye_star_HLL_flux_test,
        S_star_HLL_flux_test,
        tau_tilde_HLL_flux_test,
        Stilde_flux_HLLD_test,
    ) = calculate_HLL_fluxes(
        1,
        alpha_face_test,
        gamma_faceDD_test,
        beta_faceU_test,
        e6phi_face_test,
        u4rU_test,
        u4lU_test,
        BmagrU_test,
        BmaglU_test,
        rho_b_r_test,
        rho_b_l_test,
        Ye_r_test,
        Ye_l_test,
        S_r_test,
        S_l_test,
        P_r_test,
        P_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    )

    # At the sample values the left state sets both speed bounds, so repeat the
    #   evaluation with the face states exchanged to test the right-state inputs.
    (
        exchanged_rho_star_HLL_flux_test,
        exchanged_Ye_star_HLL_flux_test,
        exchanged_S_star_HLL_flux_test,
        exchanged_tau_tilde_HLL_flux_test,
        exchanged_Stilde_flux_HLLD_test,
    ) = calculate_HLL_fluxes(
        1,
        alpha_face_test,
        gamma_faceDD_test,
        beta_faceU_test,
        e6phi_face_test,
        u4lU_test,
        u4rU_test,
        BmaglU_test,
        BmagrU_test,
        rho_b_l_test,
        rho_b_r_test,
        Ye_l_test,
        Ye_r_test,
        S_l_test,
        S_r_test,
        P_l_test,
        P_r_test,
        h_l_test,
        h_r_test,
        cs2_l_test,
        cs2_r_test,
    )

    exprs_dict: Dict[str, Union[sp.Expr, List[sp.Expr]]] = {
        "U_rho_star": U_rho_star_test,
        "F_rho_star": F_rho_star_test,
        "U_Ye_star": U_Ye_star_test,
        "F_Ye_star": F_Ye_star_test,
        "U_S_star": U_S_star_test,
        "F_S_star": F_S_star_test,
        "U_tau_tilde": U_tau_tilde_test,
        "F_tau_tilde": F_tau_tilde_test,
        "U_S_tildeD": U_S_tildeD_test,
        "F_S_tildeD": F_S_tildeD_test,
        "HLL_test": HLL_test,
        "rho_star_HLL_flux": _nrpyAbs_to_Abs(rho_star_HLL_flux_test),
        "Ye_star_HLL_flux": _nrpyAbs_to_Abs(Ye_star_HLL_flux_test),
        "S_star_HLL_flux": _nrpyAbs_to_Abs(S_star_HLL_flux_test),
        "tau_tilde_HLL_flux": _nrpyAbs_to_Abs(tau_tilde_HLL_flux_test),
        "Stilde_flux_HLLD": [
            _nrpyAbs_to_Abs(component) for component in Stilde_flux_HLLD_test
        ],
        "exchanged_rho_star_HLL_flux": _nrpyAbs_to_Abs(
            exchanged_rho_star_HLL_flux_test
        ),
        "exchanged_Ye_star_HLL_flux": _nrpyAbs_to_Abs(exchanged_Ye_star_HLL_flux_test),
        "exchanged_S_star_HLL_flux": _nrpyAbs_to_Abs(exchanged_S_star_HLL_flux_test),
        "exchanged_tau_tilde_HLL_flux": _nrpyAbs_to_Abs(
            exchanged_tau_tilde_HLL_flux_test
        ),
        "exchanged_Stilde_flux_HLLD": [
            _nrpyAbs_to_Abs(component) for component in exchanged_Stilde_flux_HLLD_test
        ],
    }

    # Step 5: With B^i = 0, the HLL fluxes must equal the GRHD fluxes. Only the
    #         selected component index depends on flux_dirn, and a flat face
    #         metric keeps these exact comparisons small.
    zero_BmagU = ixp.zerorank1(dimension=3)
    flat_gammaDD = ixp.zerorank2(dimension=3)
    for i in range(3):
        flat_gammaDD[i][i] = sp.sympify(1)
    zero_betaU = ixp.zerorank1(dimension=3)
    unit_alpha = sp.sympify(1)
    if calculate_HLL_fluxes(
        0,
        unit_alpha,
        flat_gammaDD,
        zero_betaU,
        unit_alpha,
        u4rU_test,
        u4lU_test,
        zero_BmagU,
        zero_BmagU,
        rho_b_r_test,
        rho_b_l_test,
        Ye_r_test,
        Ye_l_test,
        S_r_test,
        S_l_test,
        P_r_test,
        P_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    ) != calculate_grhd_HLL_fluxes(
        0,
        unit_alpha,
        flat_gammaDD,
        zero_betaU,
        unit_alpha,
        u4rU_test,
        u4lU_test,
        rho_b_r_test,
        rho_b_l_test,
        Ye_r_test,
        Ye_l_test,
        S_r_test,
        S_l_test,
        P_r_test,
        P_l_test,
        h_r_test,
        h_l_test,
        cs2_r_test,
        cs2_l_test,
    ):
        raise AssertionError("B=0 HLL fluxes differ from GRHD")

    results_dict = ve.process_dictionary_of_expressions(
        exprs_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )

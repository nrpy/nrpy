"""
Construct fluxes from the HLL approximate Riemann solver at cell interfaces.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from typing import Any, Dict, List, Tuple, cast

import sympy as sp

import nrpy.indexedexp as ixp
from nrpy.equations.grhd.characteristic_speeds import find_cmax_cmin
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def calculate_GRHD_Tmunu_and_contractions(
    flux_dirn: int,
    gammaDD: List[List[sp.Expr]],
    betaU: List[sp.Expr],
    alpha: sp.Expr,
    e6phi: sp.Expr,
    rho_b: sp.Expr,
    Ye: sp.Expr,
    S: sp.Expr,
    P: sp.Expr,
    h: sp.Expr,
    u4U: List[sp.Expr],
) -> Tuple[
    sp.Expr,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    sp.Expr,
    List[sp.Expr],
    List[sp.Expr],
]:
    """
    Compute conserved quantities and fluxes entering the HLL solver.

    :param flux_dirn: Flux direction.
    :param gammaDD: Spatial metric.
    :param betaU: Shift vector.
    :param alpha: Lapse function.
    :param e6phi: Exponential conformal factor.
    :param rho_b: Baryon density.
    :param Ye: Electron fraction.
    :param S: Specific entropy.
    :param P: Pressure.
    :param h: Specific enthalpy.
    :param u4U: Four-velocity.
    :return: Conserved quantities and fluxes in the requested direction.

    Note: Written in terms of rescaled quantities, the rescaled fluxes have the
    same mathematical form as the Cartesian expressions.
    """
    # Step 1: Initialize the GRHD equations in Cartesian-equivalent form.
    grhd_eqs = GRHD_Equations(CoordSystem="Cartesian", enable_rfm_precompute=False)

    grhd_eqs.gammaDD = gammaDD.copy()
    grhd_eqs.betaU = betaU.copy()
    grhd_eqs.u4U = u4U.copy()
    grhd_eqs.alpha = alpha
    grhd_eqs.e6phi = e6phi
    grhd_eqs.rho_b = rho_b
    grhd_eqs.Ye = Ye
    grhd_eqs.S = S
    grhd_eqs.P = P
    grhd_eqs.h = h

    # Step 2: Recover the Valencia velocity and stress-energy tensor pieces.
    grhd_eqs.compute_vU_from_u4U__no_speed_limit()
    grhd_eqs.VU = grhd_eqs.VU_from_u4U
    grhd_eqs.compute_T4UU()
    grhd_eqs.compute_T4UD()

    # Step 3: Compute conserved variables.
    grhd_eqs.compute_rho_star()
    grhd_eqs.compute_Ye_star()
    grhd_eqs.compute_S_star()
    grhd_eqs.compute_tau_tilde()
    grhd_eqs.compute_S_tildeD()

    # Step 4: Compute their fluxes.
    grhd_eqs.compute_rho_star_fluxU()
    grhd_eqs.compute_Ye_star_fluxU()
    grhd_eqs.compute_S_star_fluxU()
    grhd_eqs.compute_tau_tilde_fluxU()
    grhd_eqs.compute_S_tilde_fluxUD()

    U_rho_star = grhd_eqs.rho_star
    F_rho_star = grhd_eqs.rho_star_fluxU[flux_dirn]

    U_Ye_star = grhd_eqs.Ye_star
    F_Ye_star = grhd_eqs.Ye_star_fluxU[flux_dirn]

    U_S_star = grhd_eqs.S_star
    F_S_star = grhd_eqs.S_star_fluxU[flux_dirn]

    U_tau_tilde = grhd_eqs.tau_tilde
    F_tau_tilde = grhd_eqs.tau_tilde_fluxU[flux_dirn]

    U_S_tildeD = grhd_eqs.S_tildeD.copy()
    F_S_tildeD = grhd_eqs.S_tilde_fluxUD[flux_dirn].copy()

    return (
        U_rho_star,
        F_rho_star,
        U_Ye_star,
        F_Ye_star,
        U_S_star,
        F_S_star,
        U_tau_tilde,
        F_tau_tilde,
        U_S_tildeD,
        F_S_tildeD,
    )


def HLL_solver(
    cmax: sp.Expr,
    cmin: sp.Expr,
    Fr: sp.Expr,
    Fl: sp.Expr,
    Ur: sp.Expr,
    Ul: sp.Expr,
) -> sp.Expr:
    """
    Solve the one-dimensional Riemann problem using the HLL algorithm.

    :param cmax: Maximum characteristic speed.
    :param cmin: Minimum characteristic speed.
    :param Fr: Hydrodynamic flux at the right state.
    :param Fl: Hydrodynamic flux at the left state.
    :param Ur: Conserved variable at the right state.
    :param Ul: Conserved variable at the left state.
    :return: HLL flux at the interface.
    """
    return cast(
        sp.Expr, (cmin * Fr + cmax * Fl - cmin * cmax * (Ur - Ul)) / (cmax + cmin)
    )


def calculate_HLL_fluxes(
    flux_dirn: int,
    alpha_face: sp.Expr,
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    e6phi_face: sp.Expr,
    u4rU: List[sp.Expr],
    u4lU: List[sp.Expr],
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
    Calculate symbolic HLL fluxes for the GRHD evolution system.

    :param flux_dirn: Flux direction.
    :param alpha_face: Lapse on the cell face.
    :param gamma_faceDD: Spatial metric on the cell face.
    :param beta_faceU: Shift vector on the cell face.
    :param e6phi_face: Exponential conformal factor on the cell face.
    :param u4rU: Four-velocity reconstructed to the right side.
    :param u4lU: Four-velocity reconstructed to the left side.
    :param rho_b_r: Density on the right side.
    :param rho_b_l: Density on the left side.
    :param Ye_r: Electron fraction on the right side.
    :param Ye_l: Electron fraction on the left side.
    :param S_r: Specific entropy on the right side.
    :param S_l: Specific entropy on the left side.
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
    ) = calculate_GRHD_Tmunu_and_contractions(
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
    ) = calculate_GRHD_Tmunu_and_contractions(
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
        flux_dirn, gamma_faceDD, beta_faceU, alpha_face, u4rU, u4lU, cs2_r, cs2_l
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

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    exprs_dict: Dict[str, Any] = {}

    rfm = refmetric.reference_metric["Cartesian"]

    alpha_face_test = sp.symbols("alpha_face", real=True)
    cf_face = sp.symbols("cf_face", real=True)
    h_faceDD = ixp.declarerank2("h_faceDD", symmetry="sym01", dimension=3)
    vet_faceU = ixp.declarerank1("vet_faceU", dimension=3)

    # Step 1: Define symbolic right and left interface states.
    rescaledvrU = ixp.declarerank1("rescaledvrU", dimension=3)
    rescaledvlU = ixp.declarerank1("rescaledvlU", dimension=3)

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
    (
        exprs_dict["U_rho_star"],
        exprs_dict["F_rho_star"],
        exprs_dict["U_Ye_star"],
        exprs_dict["F_Ye_star"],
        exprs_dict["U_S_star"],
        exprs_dict["F_S_star"],
        exprs_dict["U_tau_tilde"],
        exprs_dict["F_tau_tilde"],
        exprs_dict["U_S_tildeD"],
        exprs_dict["F_S_tildeD"],
    ) = calculate_GRHD_Tmunu_and_contractions(
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
        cs2_r_test,
        cs2_l_test,
    )

    cmin_test = cmin_test.subs(sp.Function("nrpyAbs"), sp.Abs)
    cmax_test = cmax_test.subs(sp.Function("nrpyAbs"), sp.Abs)

    exprs_dict["HLL_test"] = HLL_solver(
        cmax_test,
        cmin_test,
        exprs_dict["F_rho_star"],
        exprs_dict["F_Ye_star"],
        exprs_dict["U_rho_star"],
        exprs_dict["U_Ye_star"],
    )

    # Step 4: Evaluate the final HLL flux expressions.
    (
        exprs_dict["rho_star_HLL_flux"],
        exprs_dict["Ye_star_HLL_flux"],
        exprs_dict["S_star_HLL_flux"],
        exprs_dict["tau_tilde_HLL_flux"],
        exprs_dict["Stilde_flux_HLLD"],
    ) = calculate_HLL_fluxes(
        0,
        alpha_face_test,
        gamma_faceDD_test,
        beta_faceU_test,
        e6phi_face_test,
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
    )

    exprs_dict["rho_star_HLL_flux"] = exprs_dict["rho_star_HLL_flux"].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["Ye_star_HLL_flux"] = exprs_dict["Ye_star_HLL_flux"].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["S_star_HLL_flux"] = exprs_dict["S_star_HLL_flux"].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["tau_tilde_HLL_flux"] = exprs_dict["tau_tilde_HLL_flux"].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["Stilde_flux_HLLD"][0] = exprs_dict["Stilde_flux_HLLD"][0].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["Stilde_flux_HLLD"][1] = exprs_dict["Stilde_flux_HLLD"][1].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )
    exprs_dict["Stilde_flux_HLLD"][2] = exprs_dict["Stilde_flux_HLLD"][2].subs(
        sp.Function("nrpyAbs"), sp.Abs
    )

    results_dict = ve.process_dictionary_of_expressions(
        exprs_dict, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__),
        os.getcwd(),
        f"{os.path.splitext(os.path.basename(__file__))[0]}",
        results_dict,
    )

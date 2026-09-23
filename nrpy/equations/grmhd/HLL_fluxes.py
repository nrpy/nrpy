"""
Construct GRMHD fluid fluxes with the HLL approximate Riemann solver.

Magnetic stress-energy follows Duez et al., Phys. Rev. D 72, 024028
(2005), Eqs. (32)-(33); its HLL flux is Eq. (48),
https://arxiv.org/abs/astro-ph/0503420v2.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp

import nrpy.indexedexp as ixp
from nrpy.equations.grhd.HLL_fluxes import (
    HLL_solver,
    calculate_Tmunu_and_contractions_from_equations,
)
from nrpy.equations.grmhd.characteristic_speeds import find_cmax_cmin
from nrpy.equations.grmhd.GRMHD_equations import GRMHDEquations


def calculate_HLL_fluxes(
    flux_dirn: int,
    alpha_face: sp.Expr,
    gamma_faceDD: List[List[sp.Expr]],
    beta_faceU: List[sp.Expr],
    e6phi_face: sp.Expr,
    u4U_r: List[sp.Expr],
    u4U_l: List[sp.Expr],
    BmagU_r: List[sp.Expr],
    BmagU_l: List[sp.Expr],
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
    Calculate GRMHD fluid HLL fluxes across one face.

    Applies Duez et al. (2005), Eq. (48), to both magnetic face states.
    For each state, a Cartesian GRMHDEquations object supplies the
    magnetic stress-energy to the GRHD contraction function. Mass, momentum,
    and energy variables and fluxes use Duez et al. (2005), Eqs. (34)-(38),
    in reference-metric form. Electron fraction follows Jacques et al.,
    Eq. (13); entropy transport is a GRHayL choice. All face inputs,
    including BmagU_r and BmagU_l, must use the basis of gamma_faceDD; for
    reference-metric evolutions this is the rescaled basis, and the returned
    fluxes are rescaled fluxes.

    :param flux_dirn: Flux direction.
    :param alpha_face: Face lapse.
    :param gamma_faceDD: Covariant face metric.
    :param beta_faceU: Face shift vector.
    :param e6phi_face: Face reference-metric volume factor e^(6 phi), equal to sqrt(gamma/gammahat) when det(gammabar) = det(gammahat); for rescaled face data it equals sqrt(det gamma_faceDD).
    :param u4U_r: Right fluid four-velocity.
    :param u4U_l: Left fluid four-velocity.
    :param BmagU_r: Right Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param BmagU_l: Left Eulerian magnetic field scaled by 1/sqrt(4 pi), in the basis of gamma_faceDD.
    :param rho_b_r: Right baryon density.
    :param rho_b_l: Left baryon density.
    :param Ye_r: Right electron fraction.
    :param Ye_l: Left electron fraction.
    :param S_r: Right primitive entropy variable.
    :param S_l: Left primitive entropy variable.
    :param P_r: Right pressure.
    :param P_l: Left pressure.
    :param h_r: Right specific enthalpy.
    :param h_l: Left specific enthalpy.
    :param cs2_r: Right sound speed squared.
    :param cs2_l: Left sound speed squared.
    :return: Density, electron fraction, entropy, energy, and momentum HLL fluxes.
    """
    # Step 1: Compute conserved variables and physical fluxes for each face state.
    face_results = []
    for u4U, BmagU, rho_b, Ye, S, P, h in (
        (u4U_r, BmagU_r, rho_b_r, Ye_r, S_r, P_r, h_r),
        (u4U_l, BmagU_l, rho_b_l, Ye_l, S_l, P_l, h_l),
    ):
        grmhd_eqs = GRMHDEquations(CoordSystem="Cartesian")
        grmhd_eqs.BmagU = BmagU.copy()
        face_results.append(
            calculate_Tmunu_and_contractions_from_equations(
                grmhd_eqs,
                flux_dirn,
                gamma_faceDD,
                beta_faceU,
                alpha_face,
                e6phi_face,
                rho_b,
                Ye,
                S,
                P,
                h,
                u4U,
            )
        )
    (
        U_rho_r,
        F_rho_r,
        U_Ye_r,
        F_Ye_r,
        U_S_r,
        F_S_r,
        U_tau_r,
        F_tau_r,
        U_S_tilde_rD,
        F_S_tilde_rD,
    ) = face_results[0]
    (
        U_rho_l,
        F_rho_l,
        U_Ye_l,
        F_Ye_l,
        U_S_l,
        F_S_l,
        U_tau_l,
        F_tau_l,
        U_S_tilde_lD,
        F_S_tilde_lD,
    ) = face_results[1]
    # Step 2: Combine both face states with the HLL solver.
    cmin, cmax = find_cmax_cmin(
        flux_dirn,
        gamma_faceDD,
        beta_faceU,
        alpha_face,
        u4U_r,
        u4U_l,
        BmagU_r,
        BmagU_l,
        rho_b_r,
        rho_b_l,
        h_r,
        h_l,
        cs2_r,
        cs2_l,
    )
    rho_flux = HLL_solver(cmax, cmin, F_rho_r, F_rho_l, U_rho_r, U_rho_l)
    Ye_flux = HLL_solver(cmax, cmin, F_Ye_r, F_Ye_l, U_Ye_r, U_Ye_l)
    S_flux = HLL_solver(cmax, cmin, F_S_r, F_S_l, U_S_r, U_S_l)
    tau_flux = HLL_solver(cmax, cmin, F_tau_r, F_tau_l, U_tau_r, U_tau_l)
    Stilde_flux_HLLD = ixp.zerorank1(dimension=3)
    for i in range(3):
        Stilde_flux_HLLD[i] = HLL_solver(
            cmax,
            cmin,
            F_S_tilde_rD[i],
            F_S_tilde_lD[i],
            U_S_tilde_rD[i],
            U_S_tilde_lD[i],
        )
    return rho_flux, Ye_flux, S_flux, tau_flux, Stilde_flux_HLLD


if __name__ == "__main__":
    import doctest
    import os
    import sys

    import nrpy.validate_expressions.validate_expressions as ve
    from nrpy.equations.grhd.HLL_fluxes import (
        calculate_HLL_fluxes as calculate_GRHD_HLL_fluxes,
    )

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    gamma_faceDD_test = sp.eye(3).tolist()
    beta_faceU_test = ixp.zerorank1(dimension=3)
    alpha_face_test = sp.sympify(1)
    e6phi_face_test = sp.sympify(1)
    u4U_r_test = [sp.Rational(5, 4), sp.Rational(3, 4), 0, 0]
    u4U_l_test = [sp.Rational(13, 12), -sp.Rational(5, 12), 0, 0]
    rho_b_r_test = sp.sympify(1)
    rho_b_l_test = sp.Rational(4, 5)
    Ye_r_test = sp.Rational(1, 4)
    Ye_l_test = sp.Rational(1, 5)
    S_r_test = sp.Rational(3, 4)
    S_l_test = sp.Rational(4, 5)
    P_r_test = sp.Rational(1, 3)
    P_l_test = sp.Rational(1, 4)
    h_r_test = sp.Rational(3, 2)
    h_l_test = sp.Rational(7, 5)
    cs2_r_test = sp.Rational(1, 4)
    cs2_l_test = sp.Rational(1, 5)

    # Step 1: Compare B=0 HLL fluxes with GRHD for arbitrary face states in every direction.
    zero_BmagU = ixp.zerorank1(dimension=3)
    symbolic_face = (
        sp.symbols("alpha_face", real=True),
        ixp.declarerank2("gamma_faceDD", symmetry="sym01", dimension=3),
        ixp.declarerank1("beta_faceU", dimension=3),
        sp.symbols("e6phi_face", real=True),
        ixp.declarerank1("u4U_r", dimension=4),
        ixp.declarerank1("u4U_l", dimension=4),
    )
    symbolic_fluid = sp.symbols(
        "rho_b_r rho_b_l Ye_r Ye_l S_r S_l P_r P_l h_r h_l cs2_r cs2_l", real=True
    )
    for zero_field_dirn in range(3):
        grmhd_zero_fluxes = calculate_HLL_fluxes(
            zero_field_dirn, *symbolic_face, zero_BmagU, zero_BmagU, *symbolic_fluid
        )
        grhd_fluxes = calculate_GRHD_HLL_fluxes(
            zero_field_dirn, *symbolic_face, *symbolic_fluid
        )
        if grmhd_zero_fluxes != grhd_fluxes:
            raise AssertionError(
                f"B=0 HLL fluxes differ from GRHD in flux direction {zero_field_dirn}"
            )

    shared_state = (
        0,
        alpha_face_test,
        gamma_faceDD_test,
        beta_faceU_test,
        e6phi_face_test,
        u4U_r_test,
        u4U_l_test,
    )
    fluid_state = (
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

    # Step 2: Compare magnetic HLL fluxes with trusted values.
    BmagU_r_test = [sp.symbols("Bmag_r", real=True), sp.Rational(1, 3), 0]
    BmagU_l_test = [sp.symbols("Bmag_l", real=True), 0, sp.Rational(1, 4)]
    (
        sampled_rho_flux,
        sampled_Ye_flux,
        sampled_S_flux,
        sampled_tau_flux,
        sampled_Stilde_flux_HLLD,
    ) = calculate_HLL_fluxes(*shared_state, BmagU_r_test, BmagU_l_test, *fluid_state)
    expressions = {
        "rho_star_HLL_flux": sampled_rho_flux,
        "Ye_star_HLL_flux": sampled_Ye_flux,
        "S_star_HLL_flux": sampled_S_flux,
        "tau_tilde_HLL_flux": sampled_tau_flux,
        "Stilde_flux_HLLD": sampled_Stilde_flux_HLLD,
    }
    expressions = {
        key: (
            [component.subs(sp.Function("nrpyAbs"), sp.Abs) for component in value]
            if isinstance(value, list)
            else value.subs(sp.Function("nrpyAbs"), sp.Abs)
        )
        for key, value in expressions.items()
    }
    sampled_results = ve.process_dictionary_of_expressions(
        expressions, fixed_mpfs_for_free_symbols=True
    )
    ve.compare_or_generate_trusted_results(
        os.path.abspath(__file__), os.getcwd(), "HLL_fluxes", sampled_results
    )

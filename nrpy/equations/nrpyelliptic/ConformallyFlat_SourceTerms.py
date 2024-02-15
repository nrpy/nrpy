"""
Generating source terms of Hamiltonian constraint equation using a reference-metric formalism.

Author: Thiago Assumpção
        assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step Import needed modules:
from typing import cast, List, Tuple, Union
import sympy as sp  # For symbolic computations
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.nrpyelliptic.CommonParams import (
    bare_mass_0,
    bare_mass_1,
    zPunc,
    P0U,
    P1U,
    S0U,
    S1U,
)


def compute_psi_background_and_ADD_times_AUU(
    CoordSystem: str,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute conformal factor and contraction of conformal extrinsic curvature.

    :param CoordSystem: The coordinate system being used.

    :return: psi_background, ADD_times_AUU
    """
    # Step 1: Set up the reference metric
    rfm = refmetric.reference_metric[CoordSystem]

    # Step 2: Define Cartesian coordinates
    xxCart = [rfm.Cartx, rfm.Carty, rfm.Cartz]

    # Step 3: Set relative positions of each puncture

    # Step 3.a: First copy array with Cartesian coordinates
    x0U = xxCart[:]
    x1U = xxCart[:]
    # Step 3.b: Then shift z coordinates
    x0U[2] -= zPunc
    x1U[2] += zPunc

    # Step 4: Compute psi_background
    psi_background_cart = psi_background_cartesian(bare_mass_0, x0U, bare_mass_1, x1U)
    psi_background = replace_cart_coord_by_xx(rfm, psi_background_cart)

    # Step 5: Define Bowen-York vector field V_i
    VU_cart = VU_cart_two_punctures(x0U, P0U, S0U, x1U, P1U, S1U)

    # Step 6: Define extrinsic curvature
    ADD_conf_cart = ADD_conf_cartesian(xxCart, VU_cart)

    # Step 7: Compute the scalar ADD*AUU in Cartesian coordinates
    ADD_times_AUU_conf_cart = ADD_times_AUU_conf_cartesian(ADD_conf_cart)

    # Step 8: Compute ADD_times_AUU in xx coordinates
    ADD_times_AUU = replace_cart_coord_by_xx(rfm, ADD_times_AUU_conf_cart)

    return psi_background, ADD_times_AUU


# FIXME: replace with xx_to_Cart() in C code
def replace_cart_coord_by_xx(rfm: refmetric.ReferenceMetric, expr: sp.Expr) -> sp.Expr:
    """
    Replace Cartx, Carty, and Cartz by their definitions in terms of the rfm coordinates.

    :param rfm: Computed reference metric quantities
    :param expr: Original SymPy expression

    :return: The modified SymPy expression

    """
    Cart = [rfm.Cartx, rfm.Carty, rfm.Cartz]
    for k in range(3):
        expr = expr.subs(Cart[k], rfm.xx_to_Cart[k])
    return expr


def VU_cart_two_punctures(
    _x0U: List[sp.Symbol],
    _P0U: List[sp.Symbol],
    _S0U: List[sp.Symbol],
    _x1U: List[sp.Symbol],
    _P1U: List[sp.Symbol],
    _S1U: List[sp.Symbol],
) -> List[sp.Expr]:
    """
    Compute Bowen-York vector for two punctures in Cartesian coordinates.

    :param _x0U: position of first puncture
    :param _P0U: linear momentum of first puncture
    :param _S0U: spin of first puncture
    :param _x1U: position of second puncture
    :param _P1U: linear momentum of second puncture
    :param _S1U: spin of second puncture

    :return: Bowen-York vector VU in Cartesian coordinates

    """
    # Set spatial dimension
    dimension = 3

    def dot(
        vec1: Union[List[sp.Symbol], List[sp.Expr]],
        vec2: Union[List[sp.Symbol], List[sp.Expr]],
    ) -> sp.Expr:
        """Compute dot product of two vectors in 3D space, assuming Cartesian coordinates."""
        vec1_dot_vec2 = sp.sympify(0)
        for i in range(dimension):
            vec1_dot_vec2 += vec1[i] * vec2[i]
        return cast(sp.Expr, vec1_dot_vec2)

    def cross(
        vec1: Union[List[sp.Symbol], List[sp.Expr]],
        vec2: Union[List[sp.Symbol], List[sp.Expr]],
    ) -> List[sp.Expr]:
        """Compute cross product of two vectors in 3D space, assuming Cartesian coordinates."""
        vec1_cross_vec2 = ixp.zerorank1()
        LeviCivitaSymbol = ixp.LeviCivitaSymbol_dim3_rank3()
        for i in range(dimension):
            for j in range(dimension):
                for k in range(dimension):
                    vec1_cross_vec2[i] += LeviCivitaSymbol[i][j][k] * vec1[j] * vec2[k]
        return vec1_cross_vec2

    def VU_cart_single_puncture(
        xU: List[sp.Symbol], PU: List[sp.Symbol], SU: List[sp.Symbol]
    ) -> List[sp.Expr]:
        """Compute Bowen-York vector for a single puncture."""
        r = sp.sympify(0)
        for i in range(dimension):
            r += sp.Pow(xU[i], 2)
        r = sp.sqrt(r)

        VU_cart = ixp.zerorank1()
        for i in range(dimension):
            VU_cart[i] += (
                sp.Rational(-7, 4) * PU[i] / r
                + sp.Rational(-1, 4) * dot(xU, PU) * xU[i] / sp.Pow(r, 3)
                + cross(xU, SU)[i] / sp.Pow(r, 3)
            )
        return VU_cart

    # Compute Bowen-York vector for each puncture
    V0U_cart = VU_cart_single_puncture(_x0U, _P0U, _S0U)
    V1U_cart = VU_cart_single_puncture(_x1U, _P1U, _S1U)

    # Compute superposition of the two Bowen-York vectors
    VU_cart = ixp.zerorank1()
    for i in range(dimension):
        VU_cart[i] = V0U_cart[i] + V1U_cart[i]

    return VU_cart


def ADD_conf_cartesian(
    xxCart: List[sp.Expr], VU_cart: List[sp.Expr]
) -> List[List[sp.Expr]]:
    """
    Compute conformal, tracefree extrinsic curvature in Cartesian coordinates.

    :param xxCart: Vector of Cartesian coordinates
    :param VU_cart: Bowen-York vector in Cartesian coordinates

    :return: Conformal extrinsic curvature with two lower indices, in Cartesian coordinates
    """
    # Set spatial dimension
    dimension = 3

    # Define vector field V_i
    VD_cart = ixp.zerorank1()
    for i in range(dimension):
        VD_cart[i] = VU_cart[i]  # In Cartesian coordinates, V_i = V^i

    # Compute partial derivatives of the Bowen-York vector
    VD_cart_dD = ixp.zerorank2()
    for i in range(dimension):
        for j in range(dimension):
            VD_cart_dD[i][j] = sp.diff(VD_cart[i], xxCart[j])

    #  Compute the divergence
    V_cart_div = sp.sympify(0)
    for i in range(dimension):
        V_cart_div += VD_cart_dD[i][i]

    # Define Kronecker delta symbol as a function
    def kronecker_delta(i: int, j: int) -> sp.Expr:
        """Compute Kronecker delta symbol."""
        delta_ij = sp.sympify(0) if i == j else sp.sympify(1)
        return cast(sp.Expr, delta_ij)

    #  Compute the components of the conformal extrinsic curvature
    ADD_conf_cart = ixp.zerorank2()
    for i in range(dimension):
        for j in range(dimension):
            ADD_conf_cart[i][j] = (
                VD_cart_dD[i][j]
                + VD_cart_dD[j][i]
                - sp.Rational(2, 3) * kronecker_delta(i, j) * V_cart_div
            )
    return ADD_conf_cart


def ADD_times_AUU_conf_cartesian(
    ADD_conf_cart: List[List[sp.Expr]],
) -> sp.Expr:
    """
    Compute contraction of conformal extrinsic curvature in Cartesian coordinates.

    :param ADD_conf_cart: Conformal extrinsic curvature in Cartesian coordinates with two lower indices

    :return: contraction ADD_conf_cart[i][j]*ADD_conf_cart[i][j]
    """
    # Set spatial dimension
    dimension = 3

    # Compute contraction assuming Cartesian coordinates
    ADD_times_AUU_conf_cart = sp.sympify(0)
    for i in range(dimension):
        for j in range(dimension):
            ADD_times_AUU_conf_cart += ADD_conf_cart[i][j] * ADD_conf_cart[i][j]

    return cast(sp.Expr, ADD_times_AUU_conf_cart)


def psi_background_cartesian(
    _bare_mass_0: sp.Expr,
    _xU0: List[sp.Expr],
    _bare_mass_1: sp.Expr,
    _xU1: List[sp.Expr],
) -> sp.Expr:
    """
    Compute singular piece of conformal factor for two punctures in Cartesian coordinates.

    :param _bare_mass_0: bare mass of first puncture
    :param _xU0: relative position first puncture
    :param _bare_mass_1: bare mass second puncture
    :param _xU1: relative position second puncture

    :return: background conformal factor in Cartesian coordinates
    """
    # Set spatial dimension
    dimension = 3

    def psi_background_cartesian_single_puncture(
        bare_mass: sp.Expr, xU: List[sp.Expr]
    ) -> sp.Expr:
        """Compute conformal factor for a single puncture."""
        # Compute distance of puncture from the origin
        r = sp.sympify(0)
        for i in range(dimension):
            r += sp.Pow(xU[i], 2)
        r = sp.sqrt(r)

        # Compute contribution to conformal factor
        psi_single = bare_mass / (2 * r)
        return cast(sp.Expr, psi_single)

    psi_puncture_0 = psi_background_cartesian_single_puncture(_bare_mass_0, _xU0)
    psi_puncture_1 = psi_background_cartesian_single_puncture(_bare_mass_1, _xU1)

    # Compute superposition of conformal factors, adding 1 for the asymptotic limit
    psi = sp.sympify(1) + psi_puncture_0 + psi_puncture_1
    return cast(sp.Expr, psi)


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
        "SinhSpherical",
        "Cartesian",
        "SinhCylindrical",
        "SinhSymTP",
    ]:
        # Compute source terms
        _psi_background, _ADD_times_AUU = compute_psi_background_and_ADD_times_AUU(
            CoordSystem=Coord
        )

        # Declare dictionary to store expressions
        input_dict = {}

        # Extend input_dict with the symbolic expressions for the source terms
        input_dict[f"psi_background_{Coord}"] = _psi_background
        input_dict[f"ADD_times_AUU_{Coord}"] = _ADD_times_AUU

        results_dict = ve.process_dictionary_of_expressions(
            input_dict, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )

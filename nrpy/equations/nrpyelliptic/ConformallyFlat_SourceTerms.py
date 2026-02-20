"""
Generating source terms of Hamiltonian constraint equation using a reference-metric formalism.

This module constructs the initial data source terms for the Hamiltonian constraint
equation for a binary black hole system (TwoPunctures) using the Bowen-York extrinsic
curvature ansatz.

It calculates:
1. The background conformal factor, psi_background, which captures the singular
   behavior of the punctures (1/r terms).
2. The contraction of the conformal extrinsic curvature, Abar_{ij} Abar^{ij},
   derived from the Bowen-York vector field.

Author: Thiago Assumpção
        assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step Import needed modules:
from typing import List, Tuple, Union, cast

import sympy as sp  # For symbolic computations

import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.params as par
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support


def compute_psi_background_and_ADD_times_AUU(
    CoordSystem: str,
) -> Tuple[sp.Expr, sp.Expr]:
    """
    Compute the background conformal factor and the squared extrinsic curvature contraction.

    This function sets up the "Two Punctures" scenario by registering the necessary
    NRPy parameters (masses, positions, linear momenta, and spins for two black holes).
    It then computes the Bowen-York vector field and the resulting source terms required
    to solve the Hamiltonian constraint.

    :param CoordSystem: The coordinate system being used (e.g., "Spherical", "Cartesian").

    :return: A tuple containing:
             1. psi_background (psi_BG): The singular part of the conformal factor.
             2. ADD_times_AUU (Abar_{ij} Abar^{ij}): The contraction of the conformal
                extrinsic curvature tensor.
    """
    # Step 0: Define puncture parameters

    # Step 0.a: bare masses
    bare_mass_0 = par.register_CodeParameter(
        "REAL", __name__, "bare_mass_0", 0.0, commondata=True
    )
    bare_mass_1 = par.register_CodeParameter(
        "REAL", __name__, "bare_mass_1", 0.0, commondata=True
    )

    # Step 0.b: position of the punctures in the z axis
    zPunc = par.register_CodeParameter("REAL", __name__, "zPunc", 5.0, commondata=True)

    # Step 0.c.1: linear momentum 0
    P0U = par.register_CodeParameters(
        "REAL", __name__, ["P0_x", "P0_y", "P0_z"], 0.0, commondata=True
    )
    # Step 0.c.2: linear momentum 1
    P1U = par.register_CodeParameters(
        "REAL", __name__, ["P1_x", "P1_y", "P1_z"], 0.0, commondata=True
    )

    # Step 0.d.1: angular momentum 0
    S0U = par.register_CodeParameters(
        "REAL", __name__, ["S0_x", "S0_y", "S0_z"], 0.0, commondata=True
    )
    # Step 0.d.2: angular momentum 1
    S1U = par.register_CodeParameters(
        "REAL", __name__, ["S1_x", "S1_y", "S1_z"], 0.0, commondata=True
    )

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
    Replace Cartesian coordinate symbols with their definitions in the chosen reference metric.

    This maps {x, y, z} to the curvilinear coordinates {xx_0, xx_1, xx_2}
    defined in the rfm object (e.g., x = r sin(theta) cos(phi) for spherical).

    :param rfm: The ReferenceMetric object containing the coordinate mapping xx_to_Cart.
    :param expr: The original SymPy expression containing Cartesian symbols.

    :return: The modified SymPy expression in terms of reference metric coordinates.
    """
    Cart = [rfm.Cartx, rfm.Carty, rfm.Cartz]
    for k in range(3):
        expr = expr.subs(Cart[k], rfm.xx_to_Cart[k])
    return expr


def VU_cart_two_punctures(
    _x0U: List[sp.Expr],
    _P0U: List[sp.Expr],
    _S0U: List[sp.Expr],
    _x1U: List[sp.Expr],
    _P1U: List[sp.Expr],
    _S1U: List[sp.Expr],
) -> List[sp.Expr]:
    """
    Compute the superposition of Bowen-York vectors for two punctures in Cartesian coordinates.

    The Bowen-York extrinsic curvature is defined via a vector potential V^i.
    Because the momentum constraint is linear in a conformally flat background,
    the vector potential for two punctures is simply the sum of individual potentials:
    V^i_total = V^i_puncture1 + V^i_puncture2

    :param _x0U: Cartesian position vector of the first puncture relative to current point.
    :param _P0U: Linear momentum vector of the first puncture.
    :param _S0U: Spin (angular momentum) vector of the first puncture.
    :param _x1U: Cartesian position vector of the second puncture relative to current point.
    :param _P1U: Linear momentum vector of the second puncture.
    :param _S1U: Spin (angular momentum) vector of the second puncture.

    :return: The total Bowen-York vector V^i in Cartesian coordinates.
    """
    # Set spatial dimension
    dimension = 3

    def dot(
        vec1: Union[List[sp.Symbol], List[sp.Expr]],
        vec2: Union[List[sp.Symbol], List[sp.Expr]],
    ) -> sp.Expr:
        """
        Compute the Euclidean dot product of two vectors.

        vec(v)_1 . vec(v)_2 = sum_{i=1}^3 v_1^i v_2^i

        :param vec1: First 3D vector.
        :param vec2: Second 3D vector.
        :return: Scalar dot product as a sympy expression.
        """
        vec1_dot_vec2 = sp.sympify(0)
        for i in range(dimension):
            vec1_dot_vec2 += vec1[i] * vec2[i]
        return cast(sp.Expr, vec1_dot_vec2)

    def cross(
        vec1: Union[List[sp.Symbol], List[sp.Expr]],
        vec2: Union[List[sp.Symbol], List[sp.Expr]],
    ) -> List[sp.Expr]:
        """
        Compute the Euclidean cross product of two vectors.

        (vec(v)_1 x vec(v)_2)^i = epsilon^{ijk} v_{1j} v_{2k}

        :param vec1: First 3D vector.
        :param vec2: Second 3D vector.
        :return: 3D vector cross product.
        """
        vec1_cross_vec2 = ixp.zerorank1()
        LeviCivitaSymbol = ixp.LeviCivitaSymbol_dim3_rank3()
        for i in range(dimension):
            for j in range(dimension):
                for k in range(dimension):
                    vec1_cross_vec2[i] += LeviCivitaSymbol[i][j][k] * vec1[j] * vec2[k]
        return vec1_cross_vec2

    def VU_cart_single_puncture(
        xU: List[sp.Expr], PU: List[sp.Expr], SU: List[sp.Expr]
    ) -> List[sp.Expr]:
        """
        Compute the Bowen-York vector term for a single puncture.

        The formula for the vector potential V^i is:
        V^i = -7/(4r) * P^i - 1/(4r^3) * (n_j P^j) x^i + 1/r^3 * epsilon^{ijk} x_j S_k
        where r is the distance from the puncture and n^i = x^i/r.

        :param xU: Relative Cartesian position x^i.
        :param PU: Linear momentum P^i.
        :param SU: Spin S^i.
        :return: The vector V^i components.
        """
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
    Compute the conformal, tracefree extrinsic curvature Abar_{ij}.

    The physical extrinsic curvature K_{ij} is related to Abar_{ij} by the conformal factor psi:
    K_{ij} = psi^{-2} Abar_{ij}
    In the Bowen-York formalism on a conformally flat background, Abar_{ij} is constructed
    from the vector potential V^i via the conformal Killing operator:
    Abar_{ij} = partial_i V_j + partial_j V_i - (2/3) delta_{ij} (partial_k V^k)

    :param xxCart: Vector of Cartesian coordinates used for differentiation.
    :param VU_cart: Bowen-York vector potential V^i (which equals V_i in Cartesian).

    :return: Rank-2 tensor Abar_{ij}.
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
        """
        Compute the Kronecker delta symbol delta_{ij}.

        :param i: The first index.
        :param j: The second index.
        :return: 1 if i == j, else 0.
        """
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
    Compute the scalar contraction of the conformal extrinsic curvature.

    This computes the term required for the Hamiltonian constraint source:
    Abar_{ij} Abar^{ij}
    Note that in Cartesian coordinates on a flat background, indices are raised/lowered with delta_{ij},
    so this is equivalent to the sum of squares of components.

    :param ADD_conf_cart: Conformal extrinsic curvature Abar_{ij}.

    :return: The scalar contraction Abar_{ij} Abar^{ij}.
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
    Compute the "background" singular conformal factor for two punctures.

    The conformal factor psi is decomposed as psi = psi_BG + u, where
    u is non-singular and solved for numerically. psi_BG handles the
    singularities at the puncture locations analytically:
    psi_BG = 1 + sum_n (m_n / (2 |vec(r) - vec(r)_n|))

    :param _bare_mass_0: Bare mass parameter m_0 of the first puncture.
    :param _xU0: Relative Cartesian position vector of the first puncture.
    :param _bare_mass_1: Bare mass parameter m_1 of the second puncture.
    :param _xU1: Relative Cartesian position vector of the second puncture.

    :return: The scalar field psi_BG.
    """
    # Set spatial dimension
    dimension = 3

    def psi_background_cartesian_single_puncture(
        bare_mass: sp.Expr, xU: List[sp.Expr]
    ) -> sp.Expr:
        """
        Compute the m / (2r) term for a single puncture.

        :param bare_mass: The bare mass parameter.
        :param xU: Relative position vector.
        :return: The scalar contribution to the conformal factor.
        """
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

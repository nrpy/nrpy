"""
This module constructs the right-hand sides (RHSs)
  expressions for the time evolution equations of the
  BSSN gauge quantities alpha and beta^i (i.e., the
  lapse and shift). Several popular gauge conditions are
  implemented.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

# Step 1: Import all needed modules from NRPy+:
from typing import Tuple, Sequence
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import (
    BSSN_quantities,
)  # NRPy+: Computes useful BSSN quantities
from nrpy.equations.general_relativity.BSSN_RHSs import (
    BSSN_RHSs,
)  # NRPy+: Constructs BSSN right-hand-side expressions

# Step 1.a: Declare/initialize parameters for this module


def BSSN_gauge_RHSs(
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
    LapseEvolutionOption: str = "OnePlusLog",
    ShiftEvolutionOption: str = "GammaDriving2ndOrder_Covariant",
) -> Tuple[sp.Expr, Sequence[sp.Expr], Sequence[sp.Expr]]:
    """
    Core gauge evolution equation right-hand-side expression generation function.

    :param CoordSystem: Specifies the coordinate system.
    :param enable_rfm_precompute: Flag for enabling precomputation of reference metric quantities.
    :param LapseEvolutionOption: Specifies the lapse condition to use.
    :param ShiftEvolutionOption: Specifies the shift condition to use.

    :return: Returns a tuple of sympy expressions for the right-hand-side of the gauge evolution equations.
    """
    # Step 1.b: Given the chosen coordinate system, set up
    #           corresponding reference metric and needed
    #           reference metric quantities
    # The following function call sets up the reference metric
    #    and related quantities, including rescaling matrices ReDD,
    #    ReU, and hatted quantities.
    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Step 1.c: Define needed BSSN quantities:
    # Declare scalars & tensors (in terms of rescaled BSSN quantities)
    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    # Step 1.d: Declare BSSN_RHSs (excluding the time evolution equations for the gauge conditions),
    #    if they haven't already been declared.
    Brhs = BSSN_RHSs[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    ########################################
    # Step 2: Lapse conditions
    # Step 2.a: The 1+log lapse condition:
    #   \partial_t \alpha = \beta^i \alpha_{,i} - 2*\alpha*K
    # First import expressions from BSSN_quantities
    cf = Bq.cf
    trK = Bq.trK
    alpha = Bq.alpha
    betaU = Bq.betaU

    # Implement the 1+log lapse condition
    EvolvedConformalFactor_cf = par.parval_from_str("EvolvedConformalFactor_cf")
    alpha_rhs = sp.sympify(0)
    if LapseEvolutionOption == "OnePlusLog":
        alpha_rhs = -2 * alpha * trK
        alpha_dupD = ixp.declarerank1("alpha_dupD")
        for i in range(3):
            alpha_rhs += betaU[i] * alpha_dupD[i]

    # Step 2.b: Implement the harmonic slicing lapse condition
    elif LapseEvolutionOption == "HarmonicSlicing":
        if EvolvedConformalFactor_cf == "W":
            alpha_rhs = -3 * cf ** (-4) * Brhs.cf_rhs
        elif LapseEvolutionOption == "phi":
            alpha_rhs = 6 * sp.exp(6 * cf) * Brhs.cf_rhs
        else:
            raise ValueError(
                "LapseEvolutionOption==HarmonicSlicing unsupported for EvolvedConformalFactor_cf!=(W or phi)"
            )

    # Step 2.c: Frozen lapse
    #    \partial_t \alpha = 0
    elif LapseEvolutionOption == "Frozen":
        alpha_rhs = sp.sympify(0)

    # Step 2.d: Alternative 1+log lapse condition:
    #   \partial_t \alpha = \beta^i \alpha_{,i} -\alpha*(1 - \alpha)*K
    elif LapseEvolutionOption == "OnePlusLogAlt":
        alpha_rhs = -alpha * (1 - alpha) * trK
        alpha_dupD = ixp.declarerank1("alpha_dupD")
        for i in range(3):
            alpha_rhs += betaU[i] * alpha_dupD[i]

    else:
        raise ValueError(
            f"Error: LapseEvolutionOption == {LapseEvolutionOption} not supported!"
        )
    ########################################

    ########################################
    # Step 3: Shift conditions
    # Step 3.a: Set \partial_t \beta^i
    # First check that ShiftEvolutionOption parameter choice is supported.
    if ShiftEvolutionOption not in (
        "Frozen",
        "GammaDriving2ndOrder_NoCovariant",
        "GammaDriving2ndOrder_Covariant",
        "GammaDriving2ndOrder_Covariant__Hatted",
        "GammaDriving1stOrder_Covariant",
        "GammaDriving1stOrder_Covariant__Hatted",
        "NonAdvectingGammaDriving",
    ):
        raise ValueError(
            f"Error: ShiftEvolutionOption == {ShiftEvolutionOption} unsupported!"
        )

    # Next import expressions from BSSN_quantities
    BU = Bq.BU
    betU = Bq.betU
    betaU_dupD = Bq.betaU_dupD
    # Define needed quantities
    beta_rhsU = ixp.zerorank1()
    B_rhsU = ixp.zerorank1()

    # In the case of Frozen shift condition, we
    #    explicitly set the betaU and BU RHS's to zero
    #    instead of relying on the ixp.zerorank1()'s above,
    #    for safety.
    if ShiftEvolutionOption == "Frozen":
        for i in range(3):
            beta_rhsU[i] = sp.sympify(0)
            BU[i] = sp.sympify(0)

    if ShiftEvolutionOption == "GammaDriving2ndOrder_NoCovariant":
        # Step 3.a.i: Compute right-hand side of beta^i
        # *  \partial_t \beta^i = \beta^j \beta^i_{,j} + B^i
        for i in range(3):
            beta_rhsU[i] += BU[i]
            for j in range(3):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]
        # Compute right-hand side of B^i:
        eta = par.register_CodeParameter("REAL", __name__, "eta", 2.0, commondata=True)

        # Step 3.a.ii: Compute right-hand side of B^i
        # *  \partial_t B^i     = \beta^j B^i_{,j} + 3/4 * \partial_0 \Lambda^i - eta B^i
        # Step 3.a.iii: Define BU_dupD, in terms of derivative of rescaled variable \bet^i
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD", symmetry="nosym")
        for i in range(3):
            for j in range(3):
                BU_dupD[i][j] = betU_dupD[i][j] * rfm.ReU[i] + betU[i] * rfm.ReUdD[i][j]

        # Step 3.a.iv: Compute \partial_0 \bar{\Lambda}^i = (\partial_t - \beta^i \partial_i) \bar{\Lambda}^j
        Lambdabar_partial0 = ixp.zerorank1()
        for i in range(3):
            Lambdabar_partial0[i] = Brhs.Lambdabar_rhsU[i]
        for i in range(3):
            for j in range(3):
                Lambdabar_partial0[j] += -betaU[i] * Brhs.LambdabarU_dupD[j][i]

        # Step 3.a.v: Evaluate RHS of B^i:
        for i in range(3):
            B_rhsU[i] += sp.Rational(3, 4) * Lambdabar_partial0[i] - eta * BU[i]
            for j in range(3):
                B_rhsU[i] += betaU[j] * BU_dupD[i][j]

    # Step 3.b: The right-hand side of the \partial_t \beta^i equation
    if "GammaDriving2ndOrder_Covariant" in ShiftEvolutionOption:
        # Step 3.b Option 2: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + B^{i}
        ConnectionUDD = Bq.GammabarUDD
        # If instead we wish to use the Hatted covariant derivative, we replace
        #    ConnectionUDD with GammahatUDD:
        if ShiftEvolutionOption == "GammaDriving2ndOrder_Covariant__Hatted":
            ConnectionUDD = rfm.GammahatUDD
        # Then compute right-hand side:
        # Term 1: \beta^j \beta^i_{,j}
        for i in range(3):
            for j in range(3):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    beta_rhsU[i] += betaU[j] * ConnectionUDD[i][m][j] * betaU[m]
        # Term 3: B^i
        for i in range(3):
            beta_rhsU[i] += BU[i]

    if "GammaDriving2ndOrder_Covariant" in ShiftEvolutionOption:
        ConnectionUDD = Bq.GammabarUDD
        # If instead we wish to use the Hatted covariant derivative, we replace
        #    ConnectionUDD with GammahatUDD:
        if ShiftEvolutionOption == "GammaDriving2ndOrder_Covariant__Hatted":
            ConnectionUDD = rfm.GammahatUDD

        # Step 3.c: Covariant option:
        #  \partial_t B^i = \beta^j \bar{D}_j B^i
        #               + \frac{3}{4} ( \partial_t \bar{\Lambda}^{i} - \beta^j \bar{D}_j \bar{\Lambda}^{i} )
        #               - \eta B^{i}
        #                 = \beta^j B^i_{,j} + \beta^j \bar{\Gamma}^i_{mj} B^m
        #               + \frac{3}{4}[ \partial_t \bar{\Lambda}^{i}
        #                            - \beta^j (\bar{\Lambda}^i_{,j} + \bar{\Gamma}^i_{mj} \bar{\Lambda}^m)]
        #               - \eta B^{i}
        # Term 1, part a: First compute B^i_{,j} using upwinded derivative
        BU_dupD = ixp.zerorank2()
        betU_dupD = ixp.declarerank2("betU_dupD", symmetry="nosym")
        for i in range(3):
            for j in range(3):
                BU_dupD[i][j] = betU_dupD[i][j] * rfm.ReU[i] + betU[i] * rfm.ReUdD[i][j]
        # Term 1: \beta^j B^i_{,j}
        for i in range(3):
            for j in range(3):
                B_rhsU[i] += betaU[j] * BU_dupD[i][j]
        # Term 2: \beta^j \bar{\Gamma}^i_{mj} B^m
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    B_rhsU[i] += betaU[j] * ConnectionUDD[i][m][j] * BU[m]
        # Term 3: \frac{3}{4}\partial_t \bar{\Lambda}^{i}
        for i in range(3):
            B_rhsU[i] += sp.Rational(3, 4) * Brhs.Lambdabar_rhsU[i]
        # Term 4: -\frac{3}{4}\beta^j \bar{\Lambda}^i_{,j}
        for i in range(3):
            for j in range(3):
                B_rhsU[i] += -sp.Rational(3, 4) * betaU[j] * Brhs.LambdabarU_dupD[i][j]
        # Term 5: -\frac{3}{4}\beta^j \bar{\Gamma}^i_{mj} \bar{\Lambda}^m
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    B_rhsU[i] += (
                        -sp.Rational(3, 4)
                        * betaU[j]
                        * ConnectionUDD[i][m][j]
                        * Bq.LambdabarU[m]
                    )
        # Term 6: - \eta B^i
        # eta is a free parameter; we declare it here:
        eta = par.register_CodeParameter("REAL", __name__, "eta", 2.0, commondata=True)
        for i in range(3):
            B_rhsU[i] += -eta * BU[i]

    if "GammaDriving1stOrder_Covariant" in ShiftEvolutionOption:
        # Step 3.c: \partial_t \beta^i = \left[\beta^j \bar{D}_j \beta^i\right] + 3/4 Lambdabar^i - eta*beta^i

        # First set \partial_t B^i = 0:
        B_rhsU = ixp.zerorank1()  # \partial_t B^i = 0

        # Second, set \partial_t beta^i RHS:

        # Compute covariant advection term:
        ConnectionUDD = Bq.GammabarUDD
        # If instead we wish to use the Hatted covariant derivative, we replace
        #    ConnectionUDD with GammahatUDD:
        if ShiftEvolutionOption == "GammaDriving1stOrder_Covariant__Hatted":
            ConnectionUDD = rfm.GammahatUDD

        # Term 1: \beta^j \beta^i_{,j}
        for i in range(3):
            for j in range(3):
                beta_rhsU[i] += betaU[j] * betaU_dupD[i][j]

        # Term 2: \beta^j \bar{\Gamma}^i_{mj} \beta^m
        for i in range(3):
            for j in range(3):
                for m in range(3):
                    beta_rhsU[i] += betaU[j] * ConnectionUDD[i][m][j] * betaU[m]

        # Term 3: 3/4 Lambdabar^i - eta*beta^i
        eta = par.register_CodeParameter("REAL", __name__, "eta", 2.0, commondata=True)
        for i in range(3):
            beta_rhsU[i] += sp.Rational(3, 4) * Bq.LambdabarU[i] - eta * betaU[i]

    if ShiftEvolutionOption == "NonAdvectingGammaDriving":
        # Step 3.c.i: Compute right-hand side of beta^i
        # *  \partial_t \beta^i = B^i
        for i in range(3):
            beta_rhsU[i] += BU[i]

        # Compute right-hand side of B^i:
        eta = par.register_CodeParameter("REAL", __name__, "eta", 2.0, commondata=True)

        # Step 3.c.ii: Compute right-hand side of B^i
        # *  \partial_t B^i     = 3/4 * \partial_t \Lambda^i - eta B^i
        # Step 3.c.iii: Evaluate RHS of B^i:
        for i in range(3):
            B_rhsU[i] += sp.Rational(3, 4) * Brhs.Lambdabar_rhsU[i] - eta * BU[i]

    # Step 4: Rescale the BSSN gauge RHS quantities so that the evolved
    #         variables may remain smooth across coord singularities
    vet_rhsU = ixp.zerorank1()
    bet_rhsU = ixp.zerorank1()
    for i in range(3):
        vet_rhsU[i] = beta_rhsU[i] / rfm.ReU[i]
        bet_rhsU[i] = B_rhsU[i] / rfm.ReU[i]
    # Mathematica validation
    # print(str(Abar_rhsDD[2][2]).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("sin(2*x2)","Sin[2*x2]").replace("cos(x2)","Cos[x2]").replace("detgbaroverdetghat","detg"))
    # print(str(Dbarbetacontraction).replace("**","^").replace("_","").replace("xx","x").replace("sin(x2)","Sin[x2]").replace("detgbaroverdetghat","detg"))
    # print(betaU_dD)
    # print(str(trK_rhs).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))
    # print(str(bet_rhsU[0]).replace("xx2","xx3").replace("xx1","xx2").replace("xx0","xx1").replace("**","^").replace("_","").replace("sin(xx2)","Sinx2").replace("xx","x").replace("sin(2*x2)","Sin2x2").replace("cos(x2)","Cosx2").replace("detgbaroverdetghat","detg"))

    return alpha_rhs, vet_rhsU, bet_rhsU


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

    for LapseEvolOption in ["OnePlusLog", "HarmonicSlicing", "Frozen", "OnePlusLogAlt"]:
        for ShiftEvolOption in [
            "Frozen",
            "GammaDriving2ndOrder_NoCovariant",
            "GammaDriving2ndOrder_Covariant",
            "GammaDriving2ndOrder_Covariant__Hatted",
            "GammaDriving1stOrder_Covariant",
            "GammaDriving1stOrder_Covariant__Hatted",
            "NonAdvectingGammaDriving",
        ]:
            for Coord in [
                "SinhSpherical",
                "SinhSpherical_rfm_precompute",
                "SinhCartesian",
                "SinhSymTP",
            ]:
                enable_rfm_pre = "rfm_precompute" in Coord
                base_CoordSystem = Coord.replace("_rfm_precompute", "")
                inalpha_rhs, invet_rhsU, inbet_rhsU = BSSN_gauge_RHSs(
                    base_CoordSystem,
                    enable_rfm_pre,
                    LapseEvolOption,
                    ShiftEvolOption,
                )
                expr_dict = {
                    "alpha_rhs": inalpha_rhs,
                    "vet_rhsU0": invet_rhsU[0],
                    "vet_rhsU1": invet_rhsU[1],
                    "vet_rhsU2": invet_rhsU[2],
                    "bet_rhsU0": inbet_rhsU[0],
                    "bet_rhsU1": inbet_rhsU[1],
                    "bet_rhsU2": inbet_rhsU[2],
                }
                results_dict = ve.process_dictionary_of_expressions(
                    expr_dict,
                    fixed_mpfs_for_free_symbols=True,
                )
                ve.compare_or_generate_trusted_results(
                    os.path.abspath(__file__),
                    os.getcwd(),
                    # File basename. If this is set to "trusted_module_test1", then
                    #   trusted results_dict will be stored in tests/trusted_module_test1.py
                    f"{os.path.splitext(os.path.basename(__file__))[0]}_{LapseEvolOption}_{ShiftEvolOption}_{Coord}",
                    results_dict,
                )

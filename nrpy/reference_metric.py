# nrpy/reference_metric.py
"""
Define the ReferenceMetric class and its related functionalities.

This module defines the ReferenceMetric class which handles various
functionalities related to the reference metric.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         SinhSymTP fixes:
         Thiago Assumpção; assumpcaothiago **at** gmail **dot** com
"""

from typing import Any, Dict, List, Tuple, cast

import sympy as sp

import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.helpers.cached_functions import cached_simplify
from nrpy.helpers.generic import superfast_uniq

par.register_param(str, __name__, "CoordSystem_to_register_CodeParameters", "All")


class ReferenceMetric:
    """
    Handle computations and storage related to the reference metric.

    This class computes and stores quantities related to the reference
    metric for a given coordinate system.
    """

    def __init__(
        self,
        CoordSystem: str,
        enable_rfm_precompute: bool,
        SymPySimplifyExpressions: bool = False,
    ) -> None:
        """
        Initialize reference metric for the given coordinate system.

        :param CoordSystem: The coordinate system for the reference metric, such as 'Cartesian', 'Spherical', etc.
        :param enable_rfm_precompute: A Boolean indicating whether to enable precomputation of reference metric quantities.
        :param SymPySimplifyExpressions: A boolean indicating whether to simplify expressions using SymPy. Default is False.
        :raises ValueError: If an unrecognized coordinate system is specified.
        """
        self.CoordSystem = CoordSystem
        self.SymPySimplifyExpressions = SymPySimplifyExpressions

        # grid_physical_size is set based entirely on CoordSystem. So it is a rfm parameter not a grid parameter.
        par.register_param(bool, __name__, "enable_grid_physical_size", True)
        # For multipatch runs, set to "All". For single grid runs, set to the grid's CoordSystem.
        #    Otherwise you'll find this pollutes the CodeParameter namespace.
        par.register_CodeParameter(
            "REAL",
            __name__,
            "grid_physical_size",
            defaultvalue=10.0,
            add_to_glb_code_params_dict=True,
        )
        par.register_CodeParameter(
            "REAL",
            __name__,
            "grid_hole_radius",
            defaultvalue=2.0,
            add_to_glb_code_params_dict=True,
        )
        par.register_CodeParameter(
            "char[100]",
            __name__,
            "CoordSystemName",
            "unset",
            add_to_glb_code_params_dict=True,
            add_to_parfile=False,
            add_to_set_CodeParameters_h=True,
        )

        # Lower corner of the grid:
        # Annotated as List[Any] because it can hold Integers, Exprs, or Strings ("xmin")
        self.xxmin: List[Any] = [sp.sympify(0)] * 3
        # Upper corner of the grid:
        self.xxmax: List[Any] = [sp.sympify(0)] * 3
        # Physical extent of the grid;
        # grid_physical_size_dict automatically sets rfm parameters related
        #    to domain size as functions of grid_physical_size CodeParameter
        #    (usually grid_physical_size or -grid_physical_size)
        self.grid_physical_size_dict: Dict[str, str] = {}
        # grid_hole_radius_dict automatically sets rfm parameters related
        #    to the size of the hole in *Holey* CoordSystems, as functions
        #    of holey_radius.
        self.grid_hole_radius_dict: Dict[str, str] = {}
        self.add_rfm_params_to_parfile = not par.parval_from_str(
            "enable_grid_physical_size"
        )
        # Grid coordinates. In Cartesian self.xx[0],xx[1],xx[2] = x,y,z; in Spherical r, theta, phi, etc.
        self.xx = ixp.declarerank1("xx", dimension=3)
        # Cartesian coordinate; will only be a linear function of self.xx if CoordSystem==Cartesian.
        self.Cartx, self.Carty, self.Cartz = sp.symbols("Cartx Carty Cartz")
        # self.xx_to_Cart must be set as a function of (self.xx[0],xx[1],xx[2])
        self.xx_to_Cart: List[sp.Expr] = [sp.sympify(0)] * 3
        # self.Cart_to_xx must be set as a function of (Cartx, Carty, Cartz)
        # Annotated as List[Any] because it may hold "NewtonRaphson" string
        self.Cart_to_xx: List[Any] = [sp.sympify(0)] * 3
        # self.xxSph must be set as a function of (self.xx[0],xx[1],xx[2])
        self.xxSph: List[sp.Expr] = [sp.sympify(0)] * 3
        # self.scalefactor_orthog must be set as a function of (self.xx[0],xx[1],xx[2])
        self.scalefactor_orthog: List[sp.Expr] = [sp.sympify(0)] * 3
        # UnitVectors must be set as a function of (self.xx[0],xx[1],xx[2])
        self.UnitVectors: List[List[sp.Expr]] = ixp.zerorank2(dimension=3)
        # Non-angular coordinate directions
        self.radial_like_dirns: List[int] = []
        # module name for CodeParameters
        self.CodeParam_modulename = f"{__name__}"  # We don't add CoordSystem here,
        #                                            since many CoordSystems share the same parameter; e.g., AMPL.
        self.add_CodeParams_to_glb_code_params_dict = False
        self.add_CodeParams_to_glb_code_params_dict = (
            self.CoordSystem
            == par.parval_from_str("CoordSystem_to_register_CodeParameters")
        ) or (par.parval_from_str("CoordSystem_to_register_CodeParameters") == "All")
        self.requires_NewtonRaphson_for_Cart_to_xx = False
        self.NewtonRaphson_f_of_xx = [sp.sympify(0)] * 3

        # START: RFM PRECOMPUTE STUFF
        # Must be set in terms of generic functions of xx[]s
        self.scalefactor_orthog_funcform = [sp.sympify(0)] * 3

        self.f0_of_xx0_funcform = sp.Function("f0_of_xx0_funcform")(self.xx[0])
        self.f1_of_xx1_funcform = sp.Function("f1_of_xx1_funcform")(self.xx[1])
        self.f2_of_xx0_funcform = sp.Function("f2_of_xx0_funcform")(self.xx[0])
        self.f3_of_xx2_funcform = sp.Function("f3_of_xx2_funcform")(self.xx[2])
        self.f4_of_xx1_funcform = sp.Function("f4_of_xx1_funcform")(self.xx[1])

        # Annotate these as Expr because they are initialized as Symbols but later assigned full Expressions
        (
            self.f0_of_xx0,
            self.f1_of_xx1,
            self.f2_of_xx0,
            self.f3_of_xx2,
            self.f4_of_xx1,
        ) = par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            [
                "f0_of_xx0",
                "f1_of_xx1",
                "f2_of_xx0",
                "f3_of_xx2",
                "f4_of_xx1",
            ],
            add_to_parfile=False,
            add_to_glb_code_params_dict=False,
        )
        # END: RFM PRECOMPUTE STUFF

        if CoordSystem.startswith("GeneralRFM"):
            # Fully general reference metric coordinate system; coordinates are
            # treated as Cartesian-labeled with identity maps to/from Cartesian.
            self.EigenCoord = "Cartesian"
            self.general_rfm_like()
        elif "Cartesian" in CoordSystem:
            self.EigenCoord = "Cartesian"
            self.cartesian_like()
        elif "Spherical" in CoordSystem:
            self.EigenCoord = "Spherical"
            self.spherical_like()
        elif "Wedge" in CoordSystem:
            # Wedges behave the same as Cartesian -- no inner boundaries; every face is an "outer boundary":
            self.EigenCoord = "Cartesian"
            self.spherical_wedge_like()
        elif "SymTP" in CoordSystem:
            self.EigenCoord = "SymTP"
            self.prolate_spheroidal_like()
        elif "Cylindrical" in CoordSystem:
            self.EigenCoord = "Cylindrical"
            self.cylindrical_like()
        else:
            raise ValueError(
                f"Error: CoordSystem = {CoordSystem} unrecognized. Please check for typos, or add this CoordSystem support to reference_metric"
            )

        # to/from Cartesian coordinates
        def compute_Jacobian_and_inverseJacobian_tofrom_Cartesian(
            self: ReferenceMetric,
        ) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
            # Step 2.a: First construct Jacobian matrix:

            Jac_dUCart_dDrfmUD: List[List[sp.Expr]] = [
                [sp.sympify(0) for _ in range(3)] for _ in range(3)
            ]
            for i in range(3):
                for j in range(3):
                    Jac_dUCart_dDrfmUD[i][j] = sp.diff(self.xx_to_Cart[i], self.xx[j])
            Jac_dUrfm_dDCartUD, dummyDET = ixp.generic_matrix_inverter3x3(
                Jac_dUCart_dDrfmUD
            )
            return Jac_dUCart_dDrfmUD, Jac_dUrfm_dDCartUD

        (
            self.Jac_dUCart_dDrfmUD,
            self.Jac_dUrfm_dDCartUD,
        ) = compute_Jacobian_and_inverseJacobian_tofrom_Cartesian(self)

        # to/from Spherical coordinates
        def compute_Jacobian_and_inverseJacobian_tofrom_Spherical(
            self: ReferenceMetric,
        ) -> Tuple[List[List[sp.Expr]], List[List[sp.Expr]]]:
            # Step 2.a: First construct Jacobian matrix:
            Jac_dUSph_dDrfmUD = ixp.zerorank2()
            for i in range(3):
                for j in range(3):
                    Jac_dUSph_dDrfmUD[i][j] = sp.diff(self.xxSph[i], self.xx[j])
            Jac_dUrfm_dDSphUD, dummyDET = ixp.generic_matrix_inverter3x3(
                Jac_dUSph_dDrfmUD
            )
            return Jac_dUSph_dDrfmUD, Jac_dUrfm_dDSphUD

        (
            self.Jac_dUSph_dDrfmUD,
            self.Jac_dUrfm_dDSphUD,
        ) = compute_Jacobian_and_inverseJacobian_tofrom_Spherical(self)

        # Step 1: Compute ghatDD (reference metric), ghatUU
        #         (inverse reference metric), as well as
        #         rescaling vector ReU & rescaling matrix ReDD
        self.ReU = ixp.zerorank1()
        self.ReD = ixp.zerorank1()
        self.ReDD = ixp.zerorank2()
        self.ghatDD = ixp.zerorank2()
        if self.CoordSystem.startswith("GeneralRFM"):
            # GeneralRFM* reference metrics are symmetric rank-2 tensors in 3 dimensions.
            # Unlike all other reference metrics, it is *not* diagonal, and unlike
            # *some* other reference metrics, it is *nonsingular*. As such the
            # rescaling infrastructure (used to combat singular behavior) is *not needed*,
            # so we set ReU = ReD = 1, and ReDD = 1 for all components.
            for i in range(3):
                self.ReU[i] = sp.sympify(1)
                self.ReD[i] = sp.sympify(1)
                for j in range(3):
                    self.ReDD[i][j] = sp.sympify(1)

            # Treat the six independent components of ghatDD[i][j] as primitive
            # symmetric symbols.
            self.ghatDD = ixp.declarerank2("ghatDD", dimension=3, symmetry="sym01")
        elif not enable_rfm_precompute:
            for i in range(3):
                self.scalefactor_orthog[i] = sp.sympify(self.scalefactor_orthog[i])
                self.ghatDD[i][i] = self.scalefactor_orthog[i] ** 2
                self.ReU[i] = 1 / self.scalefactor_orthog[i]
                self.ReD[i] = self.scalefactor_orthog[i]
                for j in range(3):
                    self.ReDD[i][j] = (
                        self.scalefactor_orthog[i] * self.scalefactor_orthog[j]
                    )
        else:
            for i in range(3):
                self.scalefactor_orthog[i] = sp.sympify(
                    self.scalefactor_orthog_funcform[i]
                )
                self.ghatDD[i][i] = self.scalefactor_orthog_funcform[i] ** 2
                self.ReU[i] = 1 / self.scalefactor_orthog_funcform[i]
                self.ReD[i] = self.scalefactor_orthog_funcform[i]
                for j in range(3):
                    self.ReDD[i][j] = (
                        self.scalefactor_orthog_funcform[i]
                        * self.scalefactor_orthog_funcform[j]
                    )
        # Step 1b: Compute ghatUU and detgammahat
        self.ghatUU = ixp.zerorank2()
        self.ghatUU, self.detgammahat = ixp.symm_matrix_inverter3x3(self.ghatDD)

        # Step 1c: Sanity check: verify that ReDD, ghatDD,
        #          and ghatUU are all symmetric rank-2:
        for i in range(3):
            for j in range(3):
                if self.ReDD[i][j] != self.ReDD[j][i]:
                    raise ValueError(
                        f"Error: ReDD[{i}][{j}] != ReDD[{j}][{i}]: {self.ReDD[i][j]}!={self.ReDD[j][i]}"
                    )
                if self.ghatDD[i][j] != self.ghatDD[j][i]:
                    raise ValueError(
                        f"Error: ghatDD[{i}][{j}] != ghatDD[{j}][{i}]: {self.ghatDD[i][j]}!={self.ghatDD[j][i]}"
                    )
                if self.ghatUU[i][j] != self.ghatUU[j][i]:
                    raise ValueError(
                        f"Error: ghatUU[{i}][{j}] != ghatUU[{j}][{i}]: {self.ghatUU[i][j]}!={self.ghatUU[j][i]}"
                    )

        # Step 2: Compute det(ghat) and its 1st & 2nd derivatives
        self.detgammahatdD = ixp.zerorank1(3)
        self.detgammahatdDD = ixp.zerorank2(3)

        # Step 3a: Compute 1st & 2nd derivatives of rescaling vectors.
        #          (E.g., needed in BSSN for betaUdDD computation)
        self.ReUdD = ixp.zerorank2(3)
        self.ReUdDD = ixp.zerorank3(3)
        self.ReDdD = ixp.zerorank2(3)
        self.ReDdDD = ixp.zerorank3(3)

        # Step 3b: Compute 1st & 2nd derivatives of rescaling matrix.
        self.ReDDdD = ixp.zerorank3(3)
        self.ReDDdDD = ixp.zerorank4(3)

        # Step 3c: Compute 1st & 2nd derivatives of reference metric.
        self.ghatDDdD = ixp.zerorank3(3)
        self.ghatDDdDD = ixp.zerorank4(3)

        if self.CoordSystem.startswith("GeneralRFM"):
            # In GeneralRFM, the first and second derivatives of the reference
            # metric are treated as *primitive* tensors, never obtained by
            # differentiating ghatDD inside ReferenceMetric.
            # First derivatives: symmetric in (i,j).
            self.ghatDDdD = ixp.declarerank3("ghatDDdD", dimension=3, symmetry="sym01")
            self.ghatDDdDD = ixp.declarerank4(
                "ghatDDdDD", dimension=3, symmetry="sym01sym23"
            )
        else:
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        if self.SymPySimplifyExpressions:
                            self.ghatDDdD[i][j][k] = cached_simplify(
                                sp.diff(self.ghatDD[i][j], self.xx[k])
                            )
                        else:
                            self.ghatDDdD[i][j][k] = sp.diff(
                                self.ghatDD[i][j], self.xx[k]
                            )
                        for l in range(3):
                            self.ghatDDdDD[i][j][k][l] = sp.diff(
                                self.ghatDDdD[i][j][k], self.xx[l]
                            )

        # Algebraic computation of detgammahatdD using Jacobi's formula:
        # D_{,k} = D * g^{ij} * g_{ij,k}
        self.detgammahatdD = ixp.zerorank1(dimension=3)
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    self.detgammahatdD[k] += (
                        self.detgammahat * self.ghatUU[i][j] * self.ghatDDdD[i][j][k]
                    )

        # Algebraic computation of detgammahatdDD:
        # D_{,kl} = D * [ g^{mn} g_{mn,l} * g^{ij} g_{ij,k} <- TERM 1
        #                - g^{ip} g^{jq} g_{pq,l} g_{ij,k}  <- TERM 2
        #                + g^{ij} g_{ij,kl} ]               <- TERM 3
        self.detgammahatdDD = ixp.zerorank2(dimension=3)
        # TERM 1:
        # D * ( g^{mn} g_{mn,l} * g^{ij} g_{ij,k} )
        for k in range(3):
            for l in range(3):
                for i in range(3):
                    for j in range(3):
                        for m in range(3):
                            for n in range(3):
                                self.detgammahatdDD[k][l] += self.detgammahat * (
                                    self.ghatUU[m][n]
                                    * self.ghatDDdD[m][n][l]
                                    * self.ghatUU[i][j]
                                    * self.ghatDDdD[i][j][k]
                                )
        # TERM 2:
        # D * ( - g^{ip} g^{jq} g_{pq,l} g_{ij,k} )
        for k in range(3):
            for l in range(3):
                for i in range(3):
                    for j in range(3):
                        for p in range(3):
                            for q in range(3):
                                self.detgammahatdDD[k][l] += self.detgammahat * (
                                    -self.ghatUU[i][p]
                                    * self.ghatUU[j][q]
                                    * self.ghatDDdD[p][q][l]
                                    * self.ghatDDdD[i][j][k]
                                )
        # TERM 3:
        # D * ( g^{ij} g_{ij,kl} )
        for k in range(3):
            for l in range(3):
                for i in range(3):
                    for j in range(3):
                        self.detgammahatdDD[k][l] += self.detgammahat * (
                            self.ghatUU[i][j] * self.ghatDDdDD[i][j][k][l]
                        )

        for i in range(3):
            for j in range(3):
                self.ReUdD[i][j] = sp.diff(self.ReU[i], self.xx[j])
                self.ReDdD[i][j] = sp.diff(self.ReD[i], self.xx[j])
                for k in range(3):
                    self.ReUdDD[i][j][k] = sp.diff(self.ReUdD[i][j], self.xx[k])
                    self.ReDdDD[i][j][k] = sp.diff(self.ReDdD[i][j], self.xx[k])

        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ReDDdD[i][j][k] = sp.diff(self.ReDD[i][j], self.xx[k])
                    for l in range(3):
                        # Simplifying this doesn't appear to help overall NRPy run time.
                        self.ReDDdDD[i][j][k][l] = sp.diff(
                            self.ReDDdD[i][j][k], self.xx[l]
                        )

        # Step 4a: Compute Christoffel symbols of reference metric.
        # (Standard algebraic construction works for all systems including GeneralRFM)
        self.GammahatUDD = ixp.zerorank3(3)
        for i in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        self.GammahatUDD[i][k][l] += (
                            (sp.Rational(1, 2))
                            * self.ghatUU[i][m]
                            * (
                                self.ghatDDdD[m][k][l]
                                + self.ghatDDdD[m][l][k]
                                - self.ghatDDdD[k][l][m]
                            )
                        )

        # Step 4b: Compute derivs of Christoffel symbols of reference metric.
        self.GammahatUDDdD = ixp.zerorank4(3)
        # Algebraic construction of GammahatUDDdD from ghatUU, ghatDDdD, ghatDDdDD.
        # Gamma^i_{jk,l} = [g^{mi} Gamma_{mjk}],l
        #                = g^{mi}_{,l} Gamma_{mjk} + g^{mi} Gamma_{mjk,l},
        #  where Gamma_{cab} = 1/2 ( g_{ca,b} + g_{cb,a} - g_{ab,c} )
        GammahatDDD = ixp.zerorank3()
        for c in range(3):
            for a in range(3):
                for b in range(3):
                    GammahatDDD[c][a][b] += sp.Rational(1, 2) * (
                        self.ghatDDdD[c][a][b]
                        + self.ghatDDdD[c][b][a]
                        - self.ghatDDdD[a][b][c]
                    )
        #  Then, Gamma_{cab,d} = 1/2 ( g_{ca,bd} + g_{cb,ad} - g_{ab,cd} )
        GammahatDDDdD = ixp.zerorank4()
        for c in range(3):
            for a in range(3):
                for b in range(3):
                    for d in range(3):
                        GammahatDDDdD[c][a][b][d] += sp.Rational(1, 2) * (
                            self.ghatDDdDD[c][a][b][d]
                            + self.ghatDDdDD[c][b][a][d]
                            - self.ghatDDdDD[a][b][c][d]
                        )
        # So we have
        # Gamma^i_{jk,l} = [g^{mi} Gamma_{mjk}],l
        #                 = g^{mi}_{,l} Gamma_{mjk} + g^{mi} Gamma_{mjk,l},
        # where we can use the identity g^{mi}_{,l} = -g^{ip} g^{mq} g_{pq,l}
        ghatUUdD = ixp.zerorank3()
        for i in range(3):
            for l in range(3):
                for m in range(3):
                    for p in range(3):
                        for q in range(3):
                            ghatUUdD[m][i][l] += (
                                -self.ghatUU[i][p]
                                * self.ghatUU[m][q]
                                * self.ghatDDdD[p][q][l]
                            )
        # Gamma^i_{jk,l} = [g^{mi} Gamma_{mjk}],l
        #                 = g^{mi}_{,l} Gamma_{mjk} + g^{mi} Gamma_{mjk,l},
        self.GammahatUDDdD = ixp.zerorank4()
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        for m in range(3):
                            self.GammahatUDDdD[i][j][k][l] += (
                                ghatUUdD[m][i][l] * GammahatDDD[m][j][k]
                                + self.ghatUU[m][i] * GammahatDDDdD[m][j][k][l]
                            )

        # Step 4c: If rfm_precompute is disabled, then we are finished with this function.
        #          Otherwise continue to Step 5.
        if not enable_rfm_precompute or CoordSystem.startswith("GeneralRFM"):
            return

        # Step 5: Now that all hatted quantities are written in terms of generic SymPy functions,
        #         we will now replace SymPy functions with simple variables using rigid NRPy syntax,
        #         and store these variables to globals defined above.
        def make_replacements(expr: sp.Expr) -> Any:
            """
            Replace SymPy functions with simple variables using rigid NRPy syntax.

            :param expr: SymPy expression to be replaced
            :return: Expression with replaced variables
            """
            sp_version = (
                sp.__version__.lower().replace("rc", ".").replace("b", ".")
            )  # turn 1.2rc1 -> 1.2.1, 1.11b1 -> 1.11.1
            parts = sp_version.split(".")
            major = int(parts[0])
            minor = int(parts[1])
            is_old_sympy_version = (major, minor) < (1, 2)
            # The derivative representation changed with SymPy 1.2, forcing version-dependent behavior.

            # Example: Derivative(f0_of_xx0_funcform(xx0)(xx0), (xx0, 2)) >> f0_of_xx0__DD00
            rule = {}  # replacement dictionary
            for item in sp.preorder_traversal(expr):
                if item.func == sp.Derivative:
                    # extract function name before '_funcform'
                    strfunc = str(item.args[0]).split("_funcform(", 1)[0]
                    if is_old_sympy_version:
                        # extract differentiation variable and derivative order (SymPy <= 1.1)
                        var, order = str(item.args[1])[2:], len(item.args) - 1
                    else:
                        # extract differentiation variable and derivative order (SymPy >= 1.2)
                        args1 = cast(Tuple[Any, ...], item.args[1])
                        var, order = str(args1[0])[2:], args1[1]
                    # build derivative operator with format: __DD...D(var)(var)...(var) where
                    # D and (var) are repeated for every derivative order
                    oper = "__D" + "D" * (order - 1) + var * order
                    # add replacement rule to dictionary
                    rule[item] = sp.sympify(strfunc + oper)
            expr = expr.xreplace(rule)
            rule = {}

            # Example: f0_of_xx0_funcform(xx0)(xx0) >> f0_of_xx0
            for item in sp.preorder_traversal(expr):
                if "_funcform" in str(item.func):
                    # extract function name before '_funcform'
                    strfunc = str(item.func).split("_funcform", 1)[0]
                    # add replacement rule to dictionary
                    rule[item] = sp.sympify(strfunc)
            return expr.xreplace(rule)

        # Step 6: enable_rfm_precompute: precompute and store in memory
        #    expressions related to the reference metric (a.k.a.,
        #    "hatted quantities"). These expressions may involve
        #    transcendental functions, which are expensive to compute
        #    within nested loops in C. Hence we precompute them and
        #    store the result.

        # The precomputed "hatted quantity" expressions will be stored in
        #    a C struct called rfmstruct. As these expressions generally
        #    involve computationally expensive transcendental functions
        #    of xx0,xx1,or xx2, and xx0,xx1, and xx2 remain fixed across
        #    most (if not all) of a given simulation, setting up the
        #    rfmstruct can greatly improve performance.

        # The core challenge in setting up the rfmstruct is collecting
        #    all the information needed to automatically generate it.
        # Step 5 and onwards implements this algorithm, using the
        #    *generic functional form* of the hatted quantities (as
        #    opposed to the exact closed-form expressions of the
        #    hatted quantities) computed above.
        self.detgammahat = make_replacements(self.detgammahat)
        for i in range(3):
            self.ReU[i] = make_replacements(self.ReU[i])
            self.ReD[i] = make_replacements(self.ReD[i])
            self.detgammahatdD[i] = make_replacements(self.detgammahatdD[i])
            for j in range(3):
                self.ReDD[i][j] = make_replacements(self.ReDD[i][j])
                self.ReUdD[i][j] = make_replacements(self.ReUdD[i][j])
                self.ReDdD[i][j] = make_replacements(self.ReDdD[i][j])
                self.ghatDD[i][j] = make_replacements(self.ghatDD[i][j])
                self.ghatUU[i][j] = make_replacements(self.ghatUU[i][j])
                self.detgammahatdDD[i][j] = make_replacements(self.detgammahatdDD[i][j])
                for k in range(3):
                    self.ReDDdD[i][j][k] = make_replacements(self.ReDDdD[i][j][k])
                    self.ReUdDD[i][j][k] = make_replacements(self.ReUdDD[i][j][k])
                    self.ReDdDD[i][j][k] = make_replacements(self.ReDdDD[i][j][k])
                    self.ghatDDdD[i][j][k] = make_replacements(self.ghatDDdD[i][j][k])
                    self.GammahatUDD[i][j][k] = make_replacements(
                        self.GammahatUDD[i][j][k]
                    )
                    for l in range(3):
                        self.ReDDdDD[i][j][k][l] = make_replacements(
                            self.ReDDdDD[i][j][k][l]
                        )
                        self.ghatDDdDD[i][j][k][l] = make_replacements(
                            self.ghatDDdDD[i][j][k][l]
                        )
                        self.GammahatUDDdD[i][j][k][l] = make_replacements(
                            self.GammahatUDDdD[i][j][k][l]
                        )

        # Step 6: At this point, each expression is written in terms of the generic functions
        #         of xx0, xx1, and/or xx2 and their derivatives. Depending on the functions, some
        #         of these derivatives may be zero. In this step, we'll evaluate the function
        #         derivatives exactly and set the expressions to zero. Otherwise in the C code
        #         we'd be storing performing arithmetic with zeros -- wasteful!

        # Step 6.a: Construct the full list of *unique* NRPy variables representing the
        #           SymPy functions and derivatives, so that all zero derivatives can be
        #           computed.
        freevars: List[sp.Basic] = []
        freevars.extend(self.detgammahat.free_symbols)
        for i in range(3):
            freevars.extend(self.ReU[i].free_symbols)
            freevars.extend(self.ReD[i].free_symbols)
            freevars.extend(self.detgammahatdD[i].free_symbols)
            for j in range(3):
                freevars.extend(self.ReDD[i][j].free_symbols)
                freevars.extend(self.ReUdD[i][j].free_symbols)
                freevars.extend(self.ReDdD[i][j].free_symbols)
                freevars.extend(self.ghatDD[i][j].free_symbols)
                freevars.extend(self.ghatUU[i][j].free_symbols)
                freevars.extend(self.detgammahatdDD[i][j].free_symbols)
                for k in range(3):
                    freevars.extend(self.ReDDdD[i][j][k].free_symbols)
                    freevars.extend(self.ReUdDD[i][j][k].free_symbols)
                    freevars.extend(self.ReDdDD[i][j][k].free_symbols)
                    freevars.extend(self.ghatDDdD[i][j][k].free_symbols)
                    freevars.extend(self.GammahatUDD[i][j][k].free_symbols)
                    for l in range(3):
                        freevars.extend(self.ReDDdDD[i][j][k][l].free_symbols)
                        freevars.extend(self.ghatDDdDD[i][j][k][l].free_symbols)
                        freevars.extend(self.GammahatUDDdD[i][j][k][l].free_symbols)

        freevars_uniq = superfast_uniq(list(freevars))

        self.freevars_uniq_xx_indep = []
        for freevar in freevars_uniq:
            self.freevars_uniq_xx_indep.append(freevar)

        # Step 6.b: Using the expressions f?_of_xx? set in reference_metric(),
        #           evaluate each needed derivative and, in the case it is zero,
        #           set the corresponding "freevar" variable to zero.
        self.freevars_uniq_vals = []
        for i, var in enumerate(freevars_uniq):
            basename = str(var).split("__")[0].replace("_funcform", "")
            derivatv = ""
            if "__" in str(var):
                derivatv = str(var).split("__")[1].replace("_funcform", "")
            if basename == "f0_of_xx0":
                basefunc = self.f0_of_xx0
            elif basename == "f1_of_xx1":
                basefunc = self.f1_of_xx1
            elif basename == "f2_of_xx0":
                basefunc = self.f2_of_xx0
            elif basename == "f3_of_xx2":
                basefunc = self.f3_of_xx2
            elif basename == "f4_of_xx1":
                basefunc = self.f4_of_xx1
            else:
                raise ValueError(f"Error: function inside {str(var)} undefined.")
            diff_result = basefunc
            if derivatv == "":
                pass
            else:
                derivorder = (
                    derivatv.replace("d", "")
                    .replace("D", "")
                    .replace("0", "0 ")
                    .replace("1", "1 ")
                    .replace("2", "2 ")
                    .split(" ")
                )
                for derivdirn in derivorder:
                    if derivdirn != "":
                        derivwrt = self.xx[int(derivdirn)]
                        diff_result = sp.diff(diff_result, derivwrt)
            self.freevars_uniq_vals.append(diff_result)

            frees_uniq = superfast_uniq(list(diff_result.free_symbols))
            has_xx_dependence: bool = False
            for dirn in range(3):
                if self.xx[dirn] in frees_uniq:
                    has_xx_dependence = True
            if not has_xx_dependence:
                self.freevars_uniq_xx_indep[i] = diff_result

        # Step 6.c: Finally, substitute integers for all functions & derivatives that evaluate to integers
        for varidx, freevar in enumerate(freevars_uniq):
            self.detgammahat = self.detgammahat.subs(
                freevar, self.freevars_uniq_xx_indep[varidx]
            )
            for i in range(3):
                self.ReU[i] = self.ReU[i].subs(
                    freevar, self.freevars_uniq_xx_indep[varidx]
                )
                self.ReD[i] = self.ReD[i].subs(
                    freevar, self.freevars_uniq_xx_indep[varidx]
                )
                self.detgammahatdD[i] = self.detgammahatdD[i].subs(
                    freevar, self.freevars_uniq_xx_indep[varidx]
                )
                for j in range(3):
                    self.ReDD[i][j] = self.ReDD[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    self.ReUdD[i][j] = self.ReUdD[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    self.ReDdD[i][j] = self.ReDdD[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    self.ghatDD[i][j] = self.ghatDD[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    self.ghatUU[i][j] = self.ghatUU[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    self.detgammahatdDD[i][j] = self.detgammahatdDD[i][j].subs(
                        freevar, self.freevars_uniq_xx_indep[varidx]
                    )
                    for k in range(3):
                        self.ReDDdD[i][j][k] = self.ReDDdD[i][j][k].subs(
                            freevar, self.freevars_uniq_xx_indep[varidx]
                        )
                        self.ReUdDD[i][j][k] = self.ReUdDD[i][j][k].subs(
                            freevar, self.freevars_uniq_xx_indep[varidx]
                        )
                        self.ReDdDD[i][j][k] = self.ReDdDD[i][j][k].subs(
                            freevar, self.freevars_uniq_xx_indep[varidx]
                        )
                        self.ghatDDdD[i][j][k] = self.ghatDDdD[i][j][k].subs(
                            freevar, self.freevars_uniq_xx_indep[varidx]
                        )
                        self.GammahatUDD[i][j][k] = self.GammahatUDD[i][j][k].subs(
                            freevar, self.freevars_uniq_xx_indep[varidx]
                        )
                        for l in range(3):
                            self.ReDDdDD[i][j][k][l] = self.ReDDdDD[i][j][k][l].subs(
                                freevar, self.freevars_uniq_xx_indep[varidx]
                            )
                            self.ghatDDdDD[i][j][k][l] = self.ghatDDdDD[i][j][k][
                                l
                            ].subs(freevar, self.freevars_uniq_xx_indep[varidx])
                            self.GammahatUDDdD[i][j][k][l] = self.GammahatUDDdD[i][j][
                                k
                            ][l].subs(freevar, self.freevars_uniq_xx_indep[varidx])

    def Sinhv1(self, x: sp.Expr, AMPL: sp.Expr, SINHW: sp.Expr) -> Any:
        """
        Set the sinh transformation used by SinhSpherical, SinhCylindrical, and SinhCartesian.

        :param x: The input symbol for the transformation.
        :param AMPL: The amplitude of the sinh transformation.
        :param SINHW: The width of the sinh transformation.
        :return: The transformed value using the sinh function.
        """
        return (
            AMPL
            * (sp.exp(x / SINHW) - sp.exp(-x / SINHW))
            / (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
        )

    def Sinhv2(self, x: sp.Expr, AMPL: sp.Expr, SINHW: sp.Expr, slope: sp.Expr) -> Any:
        """
        Set the sinh transformation used by SinhSphericalv2*, SinhCylindricalv2, and SinhCartesianv2 (future).

        Implements the Sinhv2* coordinate transformation.

        This function computes a modified hyperbolic sine transformation used in several coordinate systems,
        including Sinhv2*, SinhCylindricalv2, and (in the future) SinhCartesianv2. It builds upon
        the base transformation defined in Sinhv1 (see the Sinhv1 function above), extending it by adding a
        polynomial prefactor and a linear slope modifier. These adjustments are designed to:
          - Suppress the linear term inherent in the sinh Taylor expansion for small x,
          - Maintain the odd-function property (i.e. r(-x) = -r(x)) necessary for symmetric grid mappings,
          - Ensure that r(1) exactly equals AMPL despite the addition of the slope term.

        Transformation Definition:
            r(x) = (AMPL - slope) * x^n * sinh(x / SINHW) / sinh(1 / SINHW) + slope * x

        where:
            - x:      The normalized coordinate, typically in the interval [0, 1].
            - AMPL:   The amplitude of the transformation (ensuring that r(1) equals AMPL).
            - SINHW:  The width parameter controlling the steepness of the hyperbolic sine.
            - slope:  A small linear modifier added to adjust the transformation for small x.
            - n:      An even integer (extracted from self.CoordSystem) that scales the sinh component by x^n.
                      This polynomial factor delays the onset of nonlinearity, making the transformation
                      approximate a simple linear behavior (slope * x) for small x.

        Design and Rationale:
            - The base transformation Sinhv1 is defined as:
                  Sinhv1(x, AMPL, SINHW) = AMPL * (exp(x/SINHW) - exp(-x/SINHW)) /
                                            (exp(1/SINHW) - exp(-1/SINHW))
              which naturally satisfies r(1) = AMPL.
            - Adding slope*x directly (i.e., the *original* sinhv2 approach) causes r(1) to become AMPL + slope.
              To correct this, the amplitude of the sinh term is reduced to (AMPL - slope) so that:
                  r(1) = (AMPL - slope) + slope = AMPL.
            - The multiplication by x^n (with n an even integer) ensures that the overall function remains
              odd and that for small x, the nonlinear sinh contribution is sufficiently suppressed.
            - When n = 2 or 4, the transformation follows the linear behavior (slope * x) for a longer portion
              of the domain before transitioning to the nonlinear sinh-dominated behavior.

        Coordinate System Naming and Parameter Extraction:
            This function uses the attribute 'self.CoordSystem', which must follow the naming convention:
                "Sinh...v2n{even_integer}"
            Specifically:
                - The string must start with "Sinh" and contain "v2n".
                - The suffix after "v2n" must be a digit representing an even integer (e.g., "2" or "4").
                  This integer is parsed and used as the exponent n (referred to as power_n in the code).
            If these conditions are not met, a ValueError is raised.

        :param x: The input symbol for the transformation.
        :param AMPL: The amplitude of the sinh transformation.
        :param SINHW: The width of the sinh transformation.
        :param slope: The slope modifier for the transformation.
        :return: The transformed value using the sinh function with an additional slope adjustment.
        :raises ValueError: If the coordinate system prefix does not start with "Sinh+EigenCoord+v2n" or if the suffix is not an even integer.
        """
        # Extract the CoordSystem prefix.
        # Check if 'v2n' is part of the string
        v2n_index = self.CoordSystem.find("v2n")
        if v2n_index != -1 and self.CoordSystem.startswith("Sinh"):
            # Capture the prefix up to "v2n"
            prefix = self.CoordSystem[: v2n_index + 3]  # '+3' to include "v2n"
        else:
            raise ValueError(
                f'Sinhv2() called for CoordSystem {self.CoordSystem}, which does not start with "Sinh+EigenCoord+v2n".'
            )
        # Extract the suffix after the prefix
        suffix = self.CoordSystem[len(prefix) :]

        # Check if the suffix is a positive even integer
        if suffix.isdigit() and int(suffix) >= 2 and int(suffix) % 2 == 0:
            power_n = int(suffix)  # Convert the suffix to an integer
        else:
            raise ValueError(
                f'What comes after "{prefix}" must be an even integer, so that r(xx0) is an odd function.'
            )
        return (
            (AMPL - slope)  # Technically this is slope*1 [1 in units of x]
            * x**power_n
            * (sp.exp(x / SINHW) - sp.exp(-x / SINHW))
            / (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
        ) + slope * x

    def cartesian_like(self) -> None:
        """Initialize class for Cartesian-like coordinate systems."""
        if self.CoordSystem == "Cartesian":
            # return values of par.Cparameters() in the following line of code are unused, so we ignore them.
            par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"],
                [-10.0, 10.0, -10.0, 10.0, -10.0, 10.0],
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = ["xmin", "ymin", "zmin"]
            self.xxmax = ["xmax", "ymax", "zmax"]
            self.grid_physical_size_dict = {
                "xmin": "-grid_physical_size",
                "ymin": "-grid_physical_size",
                "zmin": "-grid_physical_size",
                "xmax": "grid_physical_size",
                "ymax": "grid_physical_size",
                "zmax": "grid_physical_size",
            }

            self.xx_to_Cart[0] = self.xx[0]
            self.xx_to_Cart[1] = self.xx[1]
            self.xx_to_Cart[2] = self.xx[2]

            self.xxSph[0] = sp.sqrt(self.xx[0] ** 2 + self.xx[1] ** 2 + self.xx[2] ** 2)
            self.xxSph[1] = sp.acos(self.xx[2] / self.xxSph[0])
            self.xxSph[2] = sp.atan2(self.xx[1], self.xx[0])

            self.Cart_to_xx[0] = self.Cartx
            self.Cart_to_xx[1] = self.Carty
            self.Cart_to_xx[2] = self.Cartz

            self.scalefactor_orthog[0] = sp.sympify(1)
            self.scalefactor_orthog[1] = sp.sympify(1)
            self.scalefactor_orthog[2] = sp.sympify(1)

            self.scalefactor_orthog_funcform[0] = sp.sympify(1)
            self.scalefactor_orthog_funcform[1] = sp.sympify(1)
            self.scalefactor_orthog_funcform[2] = sp.sympify(1)

        elif self.CoordSystem == "SinhCartesian":
            # SinhCartesian coordinates allows us to push the outer boundary of the
            # computational domain a lot further away, while keeping reasonably high
            # resolution towards the center of the computational grid.
            # Declare basic parameters of the coordinate system and their default values
            AMPLXYZ = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "AMPLXYZ",
                10.0,
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            SINHWXYZ = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "SINHWXYZ",
                0.2,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.grid_physical_size_dict = {"AMPLXYZ": "grid_physical_size"}

            # Set default values for min and max (x,y,z)
            self.xxmin = [sp.sympify(-1), sp.sympify(-1), sp.sympify(-1)]
            self.xxmax = [sp.sympify(+1), sp.sympify(+1), sp.sympify(+1)]

            # Compute (xx_to_Cart0,xx_to_Cart1,xx_to_Cart2) from (xx0,xx1,xx2)
            for ii in [0, 1, 2]:
                self.xx_to_Cart[ii] = self.Sinhv1(self.xx[ii], AMPLXYZ, SINHWXYZ)

            # Compute (r,th,ph) from (xx_to_Cart2,xx_to_Cart1,xx_to_Cart2)
            self.xxSph[0] = sp.sqrt(
                self.xx_to_Cart[0] ** 2
                + self.xx_to_Cart[1] ** 2
                + self.xx_to_Cart[2] ** 2
            )
            self.xxSph[1] = sp.acos(self.xx_to_Cart[2] / self.xxSph[0])
            self.xxSph[2] = sp.atan2(self.xx_to_Cart[1], self.xx_to_Cart[0])

            # Compute (xx0,xx1,xx2) from (Cartx,Carty,Cartz)
            self.Cart_to_xx[0] = SINHWXYZ * sp.asinh(
                self.Cartx * sp.sinh(1 / SINHWXYZ) / AMPLXYZ
            )
            self.Cart_to_xx[1] = SINHWXYZ * sp.asinh(
                self.Carty * sp.sinh(1 / SINHWXYZ) / AMPLXYZ
            )
            self.Cart_to_xx[2] = SINHWXYZ * sp.asinh(
                self.Cartz * sp.sinh(1 / SINHWXYZ) / AMPLXYZ
            )

            # Compute scale factors
            self.scalefactor_orthog[0] = sp.diff(self.xx_to_Cart[0], self.xx[0])
            self.scalefactor_orthog[1] = sp.diff(self.xx_to_Cart[1], self.xx[1])
            self.scalefactor_orthog[2] = sp.diff(self.xx_to_Cart[2], self.xx[2])

            self.f0_of_xx0 = sp.diff(self.xx_to_Cart[0], self.xx[0])
            self.f1_of_xx1 = sp.diff(self.xx_to_Cart[1], self.xx[1])
            self.f3_of_xx2 = sp.diff(self.xx_to_Cart[2], self.xx[2])

            self.scalefactor_orthog_funcform[0] = self.f0_of_xx0_funcform
            self.scalefactor_orthog_funcform[1] = self.f1_of_xx1_funcform
            self.scalefactor_orthog_funcform[2] = self.f3_of_xx2_funcform

        # Set the transpose of the matrix of unit vectors for all Cartesian-like coordinate systems.
        self.UnitVectors = [
            [sp.sympify(1), sp.sympify(0), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
        ]
        # All Cartesian directions are radial-like; none are angular-like:
        self.radial_like_dirns = [0, 1, 2]

    def spherical_wedge_like(self) -> None:
        """Initialize class for Spherical wedge-like coordinate systems."""
        # Definitions:
        # * UWHSinhSpherical = Upper-wedge HoleySinhSpherical
        # * LWHSinhSpherical = Lower-wedge HoleySinhSpherical

        PI = self.register_pi()

        AMPL, SINHW, RMIN = par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["AMPL", "SINHW", "RMIN"],
            [10.0, 0.2, 2.0],
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        self.grid_hole_radius_dict["RMIN"] = "grid_hole_radius"
        self.grid_physical_size_dict["AMPL"] = "grid_physical_size"

        # RMIN = AMPL * sinh(x0min/SINHW) / sinh(1/SINHW)
        # -> x0min = SINHW * asinh(RMIN * sinh(1/SINHW) / AMPL)
        self.xxmin = [
            SINHW * sp.asinh(RMIN * sp.sinh(1 / SINHW) / AMPL),
            sp.Rational(1, 4) * PI,
            -PI * sp.Rational(1, 4),
        ]
        # Note that RMAX = AMPL:
        self.xxmax = [
            SINHW * sp.asinh(AMPL * sp.sinh(1 / SINHW) / AMPL),
            sp.Rational(3, 4) * PI,
            +PI * sp.Rational(1, 4),
        ]

        # BEGIN: Set universal attributes for all spherical-like coordinate systems:
        self.xxSph[0] = self.Sinhv1(self.xx[0], AMPL, SINHW)
        self.xxSph[1] = self.xx[1]
        self.xxSph[2] = self.xx[2]

        # In spherical-like coords, only the zeroth direction is radial-like; the rest are angular-like:
        self.radial_like_dirns = [0]

        # Wedges take a chunk of a spherical grid and rotate it.
        #   Scale factors are invariant under such rotations.
        self.scalefactor_orthog[0] = sp.diff(self.xxSph[0], self.xx[0])
        self.scalefactor_orthog[1] = self.xxSph[0]
        self.scalefactor_orthog[2] = self.xxSph[0] * sp.sin(self.xxSph[1])

        # fmt: off
        self.f0_of_xx0 = self.xxSph[0]
        self.f1_of_xx1 = sp.sin(self.xxSph[1])
        self.scalefactor_orthog_funcform[0] = sp.diff(self.f0_of_xx0_funcform, self.xx[0])
        self.scalefactor_orthog_funcform[1] = self.f0_of_xx0_funcform
        self.scalefactor_orthog_funcform[2] = self.f0_of_xx0_funcform * self.f1_of_xx1_funcform
        # fmt: on

        r_ito_Cart = sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
        if self.CoordSystem == "UWedgeHSinhSph":
            # Upper-wedge HoleySinhSpherical
            # This rotates the wedge 90-degrees counter-clockwise along the y-axis, causing x_old=z_new and z_old=-x_new
            #   Thus, x_new = -z_old, y_new = y_old, z_new = x_old

            # Upper wedge: x_new = -z_old, y_new = y_old, z_new = x_old
            self.Cart_to_xx[0] = SINHW * sp.asinh(
                r_ito_Cart * sp.sinh(1 / SINHW) / AMPL
            )
            # Cart_to_xx[1] = sp.acos(Cartz / rr)
            self.Cart_to_xx[1] = sp.acos(-self.Cartx / r_ito_Cart)
            # Cart_to_xx[2] = sp.atan2(Carty, Cartx)
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartz)

            # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
            #   Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
            # orig:
            # xx_to_Cart[0] = xxSph[0] * sp.sin(xxSph[1]) * sp.cos(xxSph[2])
            # xx_to_Cart[1] = xxSph[0] * sp.sin(xxSph[1]) * sp.sin(xxSph[2])
            # xx_to_Cart[2] = xxSph[0] * sp.cos(xxSph[1])

            # Upper wedge: x_new = -z_old, y_new = y_old, z_new = x_old
            # fmt: off
            self.xx_to_Cart[0] = -self.xxSph[0] * sp.cos(self.xxSph[1])
            self.xx_to_Cart[1] = +self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.sin(self.xxSph[2])
            self.xx_to_Cart[2] = +self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.cos(self.xxSph[2])
            # fmt: on

            # Set the unit vectors
            # orig:
            # UnitVectors = [[sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
            #                [sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
            #                [                -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]

            # Upper wedge: x_new = -z_old, y_new = y_old, z_new = x_old
            xxS = self.xxSph
            # fmt: off
            self.UnitVectors = [[-sp.cos(xxS[1]), sp.sin(xxS[1]) * sp.sin(xxS[2]), sp.sin(xxS[1]) * sp.cos(xxS[2])],
                                [+sp.sin(xxS[1]), sp.cos(xxS[1]) * sp.sin(xxS[2]), sp.cos(xxS[1]) * sp.cos(xxS[2])],
                                [-sp.sympify(0), sp.cos(xxS[2]), -sp.sin(xxS[2])]]
            # fmt: on
        if self.CoordSystem == "LWedgeHSinhSph":
            # Lower-wedge HoleySinhSpherical
            # This rotates the wedge 90-degrees clockwise along the y-axis, causing x_old=-z_new and z_old=x_new
            #   Thus, x_new = z_old, y_new = y_old, z_new = -x_old

            # Lower wedge: x_new = z_old, y_new = y_old, z_new = -x_old
            self.Cart_to_xx[0] = SINHW * sp.asinh(
                r_ito_Cart * sp.sinh(1 / SINHW) / AMPL
            )
            # Cart_to_xx[1] = sp.acos(Cartz / rr)
            self.Cart_to_xx[1] = sp.acos(self.Cartx / r_ito_Cart)
            # Cart_to_xx[2] = sp.atan2(Carty, Cartx)
            self.Cart_to_xx[2] = sp.atan2(self.Carty, -self.Cartz)

            # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
            #   Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
            # orig:
            # xx_to_Cart[0] = xxSph[0] * sp.sin(xxSph[1]) * sp.cos(xxSph[2])
            # xx_to_Cart[1] = xxSph[0] * sp.sin(xxSph[1]) * sp.sin(xxSph[2])
            # xx_to_Cart[2] = xxSph[0] * sp.cos(xxSph[1])

            # Lower wedge: x_new = z_old, y_new = y_old, z_new = -x_old
            # fmt: off
            self.xx_to_Cart[0] = +self.xxSph[0] * sp.cos(self.xxSph[1])
            self.xx_to_Cart[1] = +self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.sin(self.xxSph[2])
            self.xx_to_Cart[2] = -self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.cos(self.xxSph[2])

            # Set the unit vectors
            # orig:
            # UnitVectors = [[sp.sin(xxSph[1])*sp.cos(xxSph[2]), sp.sin(xxSph[1])*sp.sin(xxSph[2]),  sp.cos(xxSph[1])],
            #                [sp.cos(xxSph[1])*sp.cos(xxSph[2]), sp.cos(xxSph[1])*sp.sin(xxSph[2]), -sp.sin(xxSph[1])],
            #                [                -sp.sin(xxSph[2]),                  sp.cos(xxSph[2]),  sp.sympify(0)   ]]

            # Lower wedge: x_new = z_old, y_new = y_old, z_new = -x_old
            xxS = self.xxSph
            self.UnitVectors = [[+sp.cos(xxS[1]), sp.sin(xxS[1]) * sp.sin(xxS[2]), -sp.sin(xxS[1]) * sp.cos(xxS[2])],
                                [-sp.sin(xxS[1]), sp.cos(xxS[1]) * sp.sin(xxS[2]), -sp.cos(xxS[1]) * sp.cos(xxS[2])],
                                [sp.sympify(0), sp.cos(xxS[2]), +sp.sin(xxS[2])]]
            # fmt: on

    def spherical_like(self) -> None:
        """
        Initialize class for Spherical-like coordinate systems.

        :raises ValueError: If an unrecognized Spherical-like coordinate system is specified.
        """
        PI = self.register_pi()

        if self.CoordSystem == "Spherical":
            RMAX = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "RMAX",
                10.0,
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = [sp.sympify(0), sp.sympify(0), -PI]
            self.xxmax = [RMAX, PI, PI]
            self.grid_physical_size_dict = {"RMAX": "grid_physical_size"}

            r = self.xx[0]
            th = self.xx[1]
            ph = self.xx[2]

            self.Cart_to_xx[0] = sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
            self.Cart_to_xx[1] = sp.acos(self.Cartz / self.Cart_to_xx[0])
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)
        elif self.CoordSystem.endswith("SinhSpherical"):
            AMPL = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "AMPL",
                10.0,
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            SINHW = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "SINHW",
                0.2,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            if "HoleySinhSpherical" in self.CoordSystem:
                RMIN = par.register_CodeParameter(
                    "REAL",
                    self.CodeParam_modulename,
                    "RMIN",
                    2.0,
                    add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
                )
                # RMIN = AMPL * sinh(x0min/SINHW) / sinh(1/SINHW)
                # -> x0min = SINHW * asinh(RMIN * sinh(1/SINHW) / AMPL)
                # fmt: off
                if self.CoordSystem == "HoleySinhSpherical":
                    self.xxmin = [SINHW * sp.asinh(RMIN * sp.sinh(1 / SINHW) / AMPL), sp.sympify(0), -PI]
                    # Note that RMAX = AMPL:
                    self.xxmax = [SINHW * sp.asinh(AMPL * sp.sinh(1 / SINHW) / AMPL), PI, PI]
                elif self.CoordSystem == "RingHoleySinhSpherical":
                    self.xxmin = [SINHW * sp.asinh(RMIN * sp.sinh(1 / SINHW) / AMPL), sp.Rational(1, 4) * PI, -PI]
                    # Note that RMAX = AMPL:
                    self.xxmax = [SINHW * sp.asinh(AMPL * sp.sinh(1 / SINHW) / AMPL), sp.Rational(3, 4) * PI, PI]
                else:
                    raise ValueError(f"SinhSpherical coordinate system {self.CoordSystem} not supported!")
                # fmt: on
                self.grid_hole_radius_dict = {"RMIN": "grid_hole_radius"}
            else:
                self.xxmin = [sp.sympify(0), sp.sympify(0), -PI]
                self.xxmax = [sp.sympify(1), PI, PI]
            self.grid_physical_size_dict = {"AMPL": "grid_physical_size"}

            # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2n*".
            r = self.Sinhv1(self.xx[0], AMPL, SINHW)
            th = self.xx[1]
            ph = self.xx[2]

            rCart = sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
            self.Cart_to_xx[0] = SINHW * sp.asinh(rCart * sp.sinh(1 / SINHW) / AMPL)
            self.Cart_to_xx[1] = sp.acos(self.Cartz / rCart)
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)

        # SinhSphericalv2n[even integer] adds the parameter "r_slope", which allows for a region near self.xx[0]=0 to
        # have near-constant grid spacing Deltar. If the grid extent is Nx0, then near x0 = 0:
        # dr/dx0 \approx r_slope
        # -> dr ~ r_slope dx0
        # -> Deltar ~ r_slope * Deltax0 ~ r_slope * AMPL/Nx0
        elif self.CoordSystem.startswith("SinhSphericalv2n"):
            self.requires_NewtonRaphson_for_Cart_to_xx = True

            AMPL = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "AMPL",
                10.0,
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            SINHW = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "SINHW",
                0.2,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = [sp.sympify(0), sp.sympify(0), -PI]
            self.xxmax = [sp.sympify(1), PI, PI]
            self.grid_physical_size_dict = {"AMPL": "grid_physical_size"}
            r_slope = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "r_slope",
                0.0625,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            r = self.Sinhv2(self.xx[0], AMPL, SINHW, r_slope)
            th = self.xx[1]
            ph = self.xx[2]

            # NO CLOSED-FORM EXPRESSION FOR RADIAL INVERSION.
            self.Cart_to_xx[0] = "NewtonRaphson"
            self.NewtonRaphson_f_of_xx[0] = r - sp.sqrt(
                self.Cartx**2 + self.Carty**2 + self.Cartz**2
            )
            self.Cart_to_xx[1] = sp.acos(
                self.Cartz / sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
            )
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)
        else:
            raise ValueError(
                f"Spherical-like CoordSystem == {self.CoordSystem} unrecognized"
            )
        # BEGIN: Set universal attributes for all spherical-like coordinate systems:
        self.xxSph[0] = r
        self.xxSph[1] = th
        self.xxSph[2] = ph

        # In spherical-like coords, only the zeroth direction is radial-like; the rest are angular-like:
        self.radial_like_dirns = [0]

        # Now define xCart, yCart, and zCart in terms of x0,xx[1],xx[2].
        #   Note that the relation between r and x0 is not necessarily trivial in SinhSpherical coordinates. See above.
        self.xx_to_Cart[0] = (
            self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.cos(self.xxSph[2])
        )
        self.xx_to_Cart[1] = (
            self.xxSph[0] * sp.sin(self.xxSph[1]) * sp.sin(self.xxSph[2])
        )
        self.xx_to_Cart[2] = self.xxSph[0] * sp.cos(self.xxSph[1])

        self.scalefactor_orthog[0] = sp.diff(self.xxSph[0], self.xx[0])
        self.scalefactor_orthog[1] = self.xxSph[0]
        self.scalefactor_orthog[2] = self.xxSph[0] * sp.sin(self.xxSph[1])

        self.f0_of_xx0 = self.xxSph[0]
        self.f1_of_xx1 = sp.sin(self.xxSph[1])
        self.scalefactor_orthog_funcform[0] = sp.diff(
            self.f0_of_xx0_funcform, self.xx[0]
        )
        self.scalefactor_orthog_funcform[1] = self.f0_of_xx0_funcform
        self.scalefactor_orthog_funcform[2] = (
            self.f0_of_xx0_funcform * self.f1_of_xx1_funcform
        )

        # fmt: off
        # Set the transpose of the matrix of unit vectors
        xxS = self.xxSph
        self.UnitVectors = [[sp.sin(xxS[1]) * sp.cos(xxS[2]), sp.sin(xxS[1]) * sp.sin(xxS[2]), sp.cos(xxS[1])],
                            [sp.cos(xxS[1]) * sp.cos(xxS[2]), sp.cos(xxS[1]) * sp.sin(xxS[2]), -sp.sin(xxS[1])],
                            [-sp.sin(xxS[2]), sp.cos(xxS[2]), sp.sympify(0)]]
        # fmt: on
        # END: Set universal attributes for all spherical-like coordinate systems:

    def prolate_spheroidal_like(self) -> None:
        """
        Initialize class for Prolate spheroidal (SymTP)-like coordinate systems.

        :raises ValueError: If an unrecognized Prolate spheroidal (SymTP)-like coordinate system is specified.
        """
        PI = self.register_pi()
        SQRT1_2 = self.register_sqrt1_2()
        AMAX = par.register_CodeParameter(
            "REAL",
            self.CodeParam_modulename,
            "AMAX",
            10.0,
            add_to_parfile=self.add_rfm_params_to_parfile,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        bScale, SINHWAA = par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["bScale", "SINHWAA"],
            [0.5, 0.2],
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )

        self.xxmin = [sp.sympify(0), sp.sympify(0), -PI]
        self.xxmax = [AMAX, PI, PI]
        self.grid_physical_size_dict = {"AMAX": "grid_physical_size"}

        AA = self.xx[0]

        if self.CoordSystem == "SinhSymTP":
            self.xxmax[0] = sp.sympify(1)
            # With xxmax[0] = 1, sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and 1.
            #   Then AA = AMAX * sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and AMAX.
            AA = self.Sinhv1(self.xx[0], AMAX, SINHWAA)

        var1 = sp.sqrt(AA**2 + (bScale * sp.sin(self.xx[1])) ** 2)
        var2 = sp.sqrt(AA**2 + bScale**2)

        RHOSYMTP = AA * sp.sin(self.xx[1])
        PHSYMTP = self.xx[2]
        ZSYMTP = var2 * sp.cos(self.xx[1])

        self.xx_to_Cart[0] = AA * sp.sin(self.xx[1]) * sp.cos(self.xx[2])
        self.xx_to_Cart[1] = AA * sp.sin(self.xx[1]) * sp.sin(self.xx[2])
        self.xx_to_Cart[2] = ZSYMTP

        self.xxSph[0] = sp.sqrt(RHOSYMTP**2 + ZSYMTP**2)
        self.xxSph[1] = sp.acos(ZSYMTP / self.xxSph[0])
        self.xxSph[2] = PHSYMTP

        if self.CoordSystem == "SymTP":
            rSph = sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
            thSph = sp.acos(self.Cartz / rSph)
            phSph = sp.atan2(self.Carty, self.Cartx)

            # Mathematica script to compute Cart_to_xx[]
            #             AA = x1;
            #             var2 = Sqrt[AA^2 + bScale^2];
            #             RHOSYMTP = AA*Sin[x2];
            #             ZSYMTP = var2*Cos[x2];
            #             Solve[{rSph == Sqrt[RHOSYMTP^2 + ZSYMTP^2],
            #                    thSph == ArcCos[ZSYMTP/Sqrt[RHOSYMTP^2 + ZSYMTP^2]],
            #                    phSph == x3},
            #                   {x1, x2, x3}]
            self.Cart_to_xx[0] = (
                sp.sqrt(
                    -(bScale**2)
                    + rSph**2
                    + sp.sqrt(
                        bScale**4
                        + 2 * bScale**2 * rSph**2
                        + rSph**4
                        - 4 * bScale**2 * rSph**2 * sp.cos(thSph) ** 2
                    )
                )
                * SQRT1_2
            )  # SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            # The sign() function in the following expression ensures the correct root is taken.
            self.Cart_to_xx[1] = sp.acos(
                sp.sign(self.Cartz)
                * (
                    sp.sqrt(
                        1
                        + rSph**2 / bScale**2
                        - sp.sqrt(
                            bScale**4
                            + 2 * bScale**2 * rSph**2
                            + rSph**4
                            - 4 * bScale**2 * rSph**2 * sp.cos(thSph) ** 2
                        )
                        / bScale**2
                    )
                    * SQRT1_2
                )
            )  # SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            self.Cart_to_xx[2] = phSph

        elif self.CoordSystem == "SinhSymTP":
            # These expressions come from a Mathematica notebook shared by Thiago on Feb 26, 2025 in #nrpy2.
            #  They are validated at the bottom of this file.

            # Define intermediate terms based on Cartx, Carty, Cartz for clarity
            # Note: The sqrt_term_inner and denom_sqrt terms match the validation script.
            sqrt_term_inner = sp.sqrt(
                (self.Cartx**2 + self.Carty**2 + (bScale - self.Cartz) ** 2)
                * (self.Cartx**2 + self.Carty**2 + (bScale + self.Cartz) ** 2)
            )
            denom_sqrt_x0 = sp.sqrt(
                -(bScale**2)
                + self.Cartx**2
                + self.Carty**2
                + self.Cartz**2
                + sqrt_term_inner
            )
            denom_sqrt_x1 = sp.sqrt(
                bScale**2
                + self.Cartx**2
                + self.Carty**2
                + self.Cartz**2
                + sqrt_term_inner
            )

            # Inverse transformation derived from Mathematica (see validation script)
            def acsch(x: sp.Expr) -> sp.Expr:
                return cast(sp.Expr, sp.log(sp.sqrt(1 + x ** (-2)) + 1 / x))

            def csch(x: sp.Expr) -> sp.Expr:
                return cast(sp.Expr, 2 / (sp.exp(x) - sp.exp(-x)))

            self.Cart_to_xx[0] = SINHWAA * acsch(
                sp.sqrt(sp.Integer(2)) * AMAX * csch(1 / SINHWAA) / (denom_sqrt_x0)
            )
            self.Cart_to_xx[1] = sp.acos(sp.sqrt(2) * self.Cartz / denom_sqrt_x1)
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)
        else:
            raise ValueError(
                f"prolate-spheroidal-like CoordSystem == {self.CoordSystem} unrecognized"
            )

        # BEGIN: Set universal attributes for all prolate-spheroidal-like coordinate systems:
        self.scalefactor_orthog[0] = sp.diff(AA, self.xx[0]) * var1 / var2
        self.scalefactor_orthog[1] = var1
        self.scalefactor_orthog[2] = AA * sp.sin(self.xx[1])

        # In (Sinh)SymTP coords:
        #   xx[2] (phi) is an azimuthal angle everywhere.
        #   xx[1] (theta) approaches the usual polar angle for AA >> bScale, but as AA -> 0 it
        #   primarily parametrizes the focal line segment via z -> bScale*cos(theta). Some
        #   boundary/regularity logic therefore treats xx[1] as "radial-like" near AA ~ 0.
        self.radial_like_dirns = [0, 1]

        self.f0_of_xx0 = AA
        self.f1_of_xx1 = sp.sin(self.xx[1])
        self.f2_of_xx0 = var2
        self.f4_of_xx1 = bScale * sp.sin(self.xx[1])

        # var1 = sp.sqrt(AA**2 + (bScale * sp.sin(self.xx[1])) ** 2)
        #      = sp.sqrt(self.f0_of_xx0_funcform**2 + (bScale * self.f1_of_xx1_funcform) ** 2)

        var1_funcform = sp.sqrt(self.f0_of_xx0_funcform**2 + self.f4_of_xx1_funcform**2)
        self.scalefactor_orthog_funcform[0] = (
            sp.diff(self.f0_of_xx0_funcform, self.xx[0])
            * var1_funcform
            / self.f2_of_xx0_funcform
        )
        self.scalefactor_orthog_funcform[1] = var1_funcform
        self.scalefactor_orthog_funcform[2] = (
            self.f0_of_xx0_funcform * self.f1_of_xx1_funcform
        )

        # Set the transpose of the matrix of unit vectors
        self.UnitVectors = [
            [
                sp.sin(self.xx[1]) * sp.cos(self.xx[2]) * var2 / var1,
                sp.sin(self.xx[1]) * sp.sin(self.xx[2]) * var2 / var1,
                AA * sp.cos(self.xx[1]) / var1,
            ],
            [
                AA * sp.cos(self.xx[1]) * sp.cos(self.xx[2]) / var1,
                AA * sp.cos(self.xx[1]) * sp.sin(self.xx[2]) / var1,
                -sp.sin(self.xx[1]) * var2 / var1,
            ],
            [-sp.sin(self.xx[2]), sp.cos(self.xx[2]), sp.sympify(0)],
        ]
        # END: Set universal attributes for all prolate-spheroidal-like coordinate systems:

    def cylindrical_like(self) -> None:
        """
        Initialize class for Cylindrical-like coordinate systems.

        :raises ValueError: If an unrecognized Cylindrical-like coordinate system is specified.
        """
        PI = self.register_pi()

        # Assuming the cylindrical radial coordinate
        #   is positive makes nice simplifications of
        #   unit vectors possible.
        if self.CoordSystem == "Cylindrical":
            RHOMAX, ZMIN, ZMAX = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["RHOMAX", "ZMIN", "ZMAX"],
                [10.0, -10.0, 10.0],
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = [sp.sympify(0), -PI, ZMIN]
            self.xxmax = [RHOMAX, PI, ZMAX]
            self.grid_physical_size_dict = {
                "RHOMAX": "grid_physical_size",
                "ZMIN": "-grid_physical_size",
                "ZMAX": "grid_physical_size",
            }

            RHOCYL = self.xx[0]
            PHICYL = self.xx[1]
            ZCYL = self.xx[2]

            self.Cart_to_xx[0] = sp.sqrt(self.Cartx**2 + self.Carty**2)
            self.Cart_to_xx[1] = sp.atan2(self.Carty, self.Cartx)
            self.Cart_to_xx[2] = self.Cartz

        elif self.CoordSystem == "SinhCylindrical":
            AMPLRHO, AMPLZ = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["AMPLRHO", "AMPLZ"],
                [10.0, 10.0],
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            SINHWRHO, SINHWZ = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["SINHWRHO", "SINHWZ"],
                [0.2, 0.2],
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )

            self.xxmin = [sp.sympify(0), -PI, sp.sympify(-1)]
            self.xxmax = [sp.sympify(1), PI, sp.sympify(+1)]
            self.grid_physical_size_dict = {
                "AMPLRHO": "grid_physical_size",
                "AMPLZ": "grid_physical_size",
            }

            # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
            RHOCYL = self.Sinhv1(self.xx[0], AMPLRHO, SINHWRHO)
            # phi coordinate remains unchanged.
            PHICYL = self.xx[1]
            ZCYL = self.Sinhv1(self.xx[2], AMPLZ, SINHWZ)
            self.Cart_to_xx[0] = SINHWRHO * sp.asinh(
                sp.sqrt(self.Cartx**2 + self.Carty**2) * sp.sinh(1 / SINHWRHO) / AMPLRHO
            )
            self.Cart_to_xx[1] = sp.atan2(self.Carty, self.Cartx)
            self.Cart_to_xx[2] = SINHWZ * sp.asinh(
                self.Cartz * sp.sinh(1 / SINHWZ) / AMPLZ
            )

        # SinhCylindricalv2n[even integer] adds the parameters "rho_slope" & "z_slope", which allows for a region
        # near self.xx[0]=0 and self.xx[2]=0 to
        # have near-constant grid spacing Deltarho and Deltaz, respectively.
        # For example, if the grid extent in rho is Nx0, then near x0 = 0:
        # drho/dx0 \approx rho_slope
        # -> dr ~ rho_slope dx0
        # -> Deltarho ~ rho_slope * Deltax0 ~ rho_slope * AMPLRHO/Nx0
        elif self.CoordSystem.startswith("SinhCylindricalv2n"):
            self.requires_NewtonRaphson_for_Cart_to_xx = True

            AMPLRHO, AMPLZ = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["AMPLRHO", "AMPLZ"],
                [10.0, 10.0],
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            SINHWRHO, SINHWZ = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["SINHWRHO", "SINHWZ"],
                [0.2, 0.2],
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = [sp.sympify(0), -PI, sp.sympify(-1)]
            self.xxmax = [sp.sympify(1), PI, sp.sympify(+1)]
            self.grid_physical_size_dict = {
                "AMPLRHO": "grid_physical_size",
                "AMPLZ": "grid_physical_size",
            }
            rho_slope, z_slope = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["rho_slope", "z_slope"],
                [0.0625, 0.0625],
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )

            RHOCYL = self.Sinhv2(self.xx[0], AMPLRHO, SINHWRHO, rho_slope)
            PHICYL = self.xx[1]
            ZCYL = self.Sinhv2(self.xx[2], AMPLZ, SINHWZ, z_slope)

            # NO CLOSED-FORM EXPRESSION FOR RADIAL OR Z INVERSION.
            self.Cart_to_xx[0] = "NewtonRaphson"
            self.NewtonRaphson_f_of_xx[0] = RHOCYL - sp.sqrt(
                self.Cartx**2 + self.Carty**2
            )
            self.Cart_to_xx[1] = sp.atan2(self.Carty, self.Cartx)
            self.Cart_to_xx[2] = "NewtonRaphson"
            self.NewtonRaphson_f_of_xx[2] = ZCYL - self.Cartz
        else:
            raise ValueError(
                f"Cylindrical-like CoordSystem == {self.CoordSystem} unrecognized"
            )

        # START: Set universal attributes for all Cylindrical-like coordinate systems:
        self.xx_to_Cart[0] = RHOCYL * sp.cos(PHICYL)
        self.xx_to_Cart[1] = RHOCYL * sp.sin(PHICYL)
        self.xx_to_Cart[2] = ZCYL

        # In cylindrical-like coords, only the 1st coordinate direction is angular; the rest are radial-like:
        self.radial_like_dirns = [0, 2]

        self.xxSph[0] = sp.sqrt(RHOCYL**2 + ZCYL**2)
        self.xxSph[1] = sp.acos(ZCYL / self.xxSph[0])
        self.xxSph[2] = PHICYL

        self.scalefactor_orthog[0] = sp.diff(RHOCYL, self.xx[0])
        self.scalefactor_orthog[1] = RHOCYL
        self.scalefactor_orthog[2] = sp.diff(ZCYL, self.xx[2])

        self.f0_of_xx0 = RHOCYL
        self.f3_of_xx2 = sp.diff(ZCYL, self.xx[2])

        self.scalefactor_orthog_funcform[0] = sp.diff(
            self.f0_of_xx0_funcform, self.xx[0]
        )
        self.scalefactor_orthog_funcform[1] = self.f0_of_xx0_funcform
        self.scalefactor_orthog_funcform[2] = self.f3_of_xx2_funcform

        # Set the transpose of the matrix of unit vectors
        self.UnitVectors = [
            [sp.cos(PHICYL), sp.sin(PHICYL), sp.sympify(0)],
            [-sp.sin(PHICYL), sp.cos(PHICYL), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
        ]
        # END: Set universal attributes for all Cylindrical-like coordinate systems:

    def register_pi(self) -> sp.Expr:
        """
        Register the mathematical constant pi as a code parameter and return PI as a sympy symbol.

        :return: The sympy symbol "PI".
        """
        return par.register_CodeParameter(
            "REAL",
            __name__,
            "PI",
            "3.14159265358979323846264338327950288",
            add_to_parfile=False,
            add_to_set_CodeParameters_h=True,
            add_to_glb_code_params_dict=True,
        )

    def register_sqrt1_2(self) -> sp.Expr:
        """
        Register the mathematical constant sqrt(1/2) as a code parameter and return SQRT1_2 as a sympy symbol.

        :return: The sympy symbol "SQRT1_2".
        """
        return par.register_CodeParameter(
            "REAL",
            __name__,
            "SQRT1_2",
            "0.707106781186547524400844362105",
            add_to_parfile=False,
            add_to_set_CodeParameters_h=True,
            add_to_glb_code_params_dict=True,
        )

    def general_rfm_like(self) -> None:
        """
        Initialize class attributes for the GeneralRFM coordinate system.

        The coordinates xx^i are treated as Cartesian-labeled coordinates
        with identity maps to and from the Cartesian variables (Cartx,
        Carty, Cartz). The reference metric itself is specified separately
        and may be fully general and non-orthogonal.
        """
        # Neutral default domain; users are free to override these bounds.
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["xmin", "xmax", "ymin", "ymax", "zmin", "zmax"],
            [-1.0, 1.0, -1.0, 1.0, -1.0, 1.0],
            add_to_parfile=self.add_rfm_params_to_parfile,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        self.xxmin = ["xmin", "ymin", "zmin"]
        self.xxmax = ["xmax", "ymax", "zmax"]
        self.grid_physical_size_dict = {
            "xmin": "-grid_physical_size",
            "ymin": "-grid_physical_size",
            "zmin": "-grid_physical_size",
            "xmax": "grid_physical_size",
            "ymax": "grid_physical_size",
            "zmax": "grid_physical_size",
        }

        # Identity map between reference-metric coordinates and Cartesian.
        self.xx_to_Cart[0] = self.xx[0]
        self.xx_to_Cart[1] = self.xx[1]
        self.xx_to_Cart[2] = self.xx[2]

        self.Cart_to_xx[0] = self.Cartx
        self.Cart_to_xx[1] = self.Carty
        self.Cart_to_xx[2] = self.Cartz

        # Define spherical coordinates from these Cartesian-like coordinates
        # (for ancillary purposes only; not used for metric construction).
        self.xxSph[0] = sp.sqrt(self.xx[0] ** 2 + self.xx[1] ** 2 + self.xx[2] ** 2)
        self.xxSph[1] = sp.acos(self.xx[2] / self.xxSph[0])
        self.xxSph[2] = sp.atan2(self.xx[1], self.xx[0])

        # Trivial scale factors (placeholders; not used to build metric in GeneralRFM).
        self.scalefactor_orthog[0] = sp.sympify(1)
        self.scalefactor_orthog[1] = sp.sympify(1)
        self.scalefactor_orthog[2] = sp.sympify(1)

        # Unit vectors: identity; all directions are "radial-like".
        self.UnitVectors = [
            [sp.sympify(1), sp.sympify(0), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
        ]
        self.radial_like_dirns = [0, 1, 2]


class rfm_dict(Dict[str, ReferenceMetric]):
    """Custom dictionary for storing ReferenceMetric objects."""

    def __getitem__(self, CoordSystem_in: str) -> ReferenceMetric:
        if CoordSystem_in not in self:
            print(f"Setting up reference_metric[{CoordSystem_in}]...")

            # In case [CoordSystem]_rfm_precompute is passed:
            CoordSystem = CoordSystem_in.replace("_rfm_precompute", "")

            # Base ReferenceMetric without precomputation:
            base_rfm = ReferenceMetric(CoordSystem, enable_rfm_precompute=False)
            self.__setitem__(CoordSystem, base_rfm)

            # For GeneralRFM, rfm_precompute is enabled by design; map the
            # _rfm_precompute variant to the same underlying object.
            if CoordSystem.startswith("GeneralRFM"):
                self.__setitem__(CoordSystem + "_rfm_precompute", base_rfm)
            else:
                self.__setitem__(
                    CoordSystem + "_rfm_precompute",
                    ReferenceMetric(CoordSystem, enable_rfm_precompute=True),
                )
        return dict.__getitem__(self, CoordSystem_in)

    def __setitem__(self, CoordSystem: str, value: ReferenceMetric) -> None:
        dict.__setitem__(self, CoordSystem, value)

    def __delitem__(self, CoordSystem: str) -> None:
        dict.__delitem__(self, CoordSystem)


reference_metric = rfm_dict()

supported_CoordSystems = [
    "Spherical",
    "SinhSpherical",
    "SinhSphericalv2n2",
    "Cartesian",
    "SinhCartesian",
    "Cylindrical",
    "SinhCylindrical",
    "SinhCylindricalv2n2",
    "SymTP",
    "SinhSymTP",
    "LWedgeHSinhSph",
    "UWedgeHSinhSph",
    "RingHoleySinhSpherical",
    "HoleySinhSpherical",
    "GeneralRFM_fisheyeN2",
]

unittest_CoordSystems = [
    "SinhSymTP",
    "HoleySinhSpherical",
    "Cartesian",
    "SinhCylindricalv2n2",
]

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

    for Coord in supported_CoordSystems:
        rfm = reference_metric[Coord]
        results_dict = ve.process_dictionary_of_expressions(
            rfm.__dict__, fixed_mpfs_for_free_symbols=True
        )
        ve.compare_or_generate_trusted_results(
            os.path.abspath(__file__),
            os.getcwd(),
            # File basename. If this is set to "trusted_module_test1", then
            #   trusted results_dict will be stored in tests/trusted_module_test1.py
            f"{os.path.splitext(os.path.basename(__file__))[0]}_{Coord}",
            results_dict,
        )

    # --------------------------------------------------------------------------
    # Add specific validation test for SinhSymTP coordinate inversion here
    # --------------------------------------------------------------------------
    # Get the ReferenceMetric object for SinhSymTP
    rfm_sinhsymtp = reference_metric["SinhSymTP"]

    # Define the substitution dictionary: xx[0], xx[1], xx[2] mapped to Cart_to_xx expressions
    # Typing as Dict[Any, Any] resolves the variance issue with sp.Expr.subs() expecting strict Mappings
    sub_dict_backward_in_forward: Dict[Any, Any] = {
        rfm_sinhsymtp.xx[0]: rfm_sinhsymtp.Cart_to_xx[0],
        rfm_sinhsymtp.xx[1]: rfm_sinhsymtp.Cart_to_xx[1],
        rfm_sinhsymtp.xx[2]: rfm_sinhsymtp.Cart_to_xx[2],
    }

    # Substitute backward into forward for each Cartesian component
    Cartx_check = rfm_sinhsymtp.xx_to_Cart[0].subs(sub_dict_backward_in_forward)
    Carty_check = rfm_sinhsymtp.xx_to_Cart[1].subs(sub_dict_backward_in_forward)
    Cartz_check = rfm_sinhsymtp.xx_to_Cart[2].subs(sub_dict_backward_in_forward)

    # Calculate the differences (should be zero if inversion is correct)
    diff_Cartx = Cartx_check - rfm_sinhsymtp.Cartx
    diff_Carty = Carty_check - rfm_sinhsymtp.Carty
    diff_Cartz = Cartz_check - rfm_sinhsymtp.Cartz

    is_zero_Cartx = ve.check_zero(
        diff_Cartx, verbose=False, fixed_mpfs_for_free_symbols=True, hex_offset=1
    )
    is_zero_Carty = ve.check_zero(
        diff_Carty, verbose=False, fixed_mpfs_for_free_symbols=True, hex_offset=0
    )
    is_zero_Cartz = ve.check_zero(
        diff_Cartz, verbose=False, fixed_mpfs_for_free_symbols=True, hex_offset=0
    )

    if is_zero_Cartx and is_zero_Carty and is_zero_Cartz:
        print("PASS: SinhSymTP coordinate inversion validation (Cart -> xx -> Cart)")
    else:
        print(
            "FAILURE: SinhSymTP coordinate inversion validation (Cart -> xx -> Cart) FAILED."
        )
        print(is_zero_Cartx, is_zero_Carty, is_zero_Cartz)
        sys.exit(1)

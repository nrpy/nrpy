"""
Define the ReferenceMetric class and its related functionalities.

This module defines the ReferenceMetric class which handles various
functionalities related to the reference metric.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Dict, List, Any, Tuple, cast
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.helpers.generic import superfast_uniq
from nrpy.helpers.cached_functions import cached_simplify

# grid_physical_size is set based entirely on CoordSystem. So it is a rfm parameter not a grid parameter.
par.register_param(bool, __name__, "enable_grid_physical_size", True)
# For multipatch runs, set to "All". For single grid runs, set to the grid's CoordSystem.
#    Otherwise you'll find this pollutes the CodeParameter namespace.
par.register_param(str, __name__, "CoordSystem_to_register_CodeParameters", "All")
par.register_CodeParameter(
    "REAL",
    __name__,
    "grid_physical_size",
    defaultvalue=10.0,
    add_to_glb_code_params_dict=True,
)
par.register_CodeParameter(
    "char[50]",
    __name__,
    "CoordSystemName",
    "must set",
    add_to_glb_code_params_dict=True,
    add_to_parfile=False,
    add_to_set_CodeParameters_h=True,
)


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
        """
        self.CoordSystem = CoordSystem
        self.SymPySimplifyExpressions = SymPySimplifyExpressions

        # Lower corner of the grid:
        self.xxmin = [sp.sympify(0)] * 3
        # Upper corner of the grid:
        self.xxmax = [sp.sympify(0)] * 3
        # Physical extent of the grid;
        # grid_physical_size_dict automatically sets rfm parameters related
        #    to domain size as functions of grid_physical_size CodeParameter
        #    (usually grid_physical_size or -grid_physical_size)
        self.grid_physical_size_dict: Dict[str, str] = {}
        self.add_rfm_params_to_parfile = not par.parval_from_str(
            "enable_grid_physical_size"
        )
        # Grid coordinates. In Cartesian self.xx[0],xx[1],xx[2] = x,y,z; in Spherical r, theta, phi, etc.
        self.xx = cast(List[sp.Symbol], ixp.declarerank1("xx", dimension=3))
        # Cartesian coordinate; will only be a linear function of self.xx if CoordSystem==Cartesian.
        self.Cartx, self.Carty, self.Cartz = sp.symbols("Cartx Carty Cartz")
        # self.xx_to_Cart must be set as a function of (self.xx[0],xx[1],xx[2])
        self.xx_to_Cart = [sp.sympify(0)] * 3
        # self.Cart_to_xx must be set as a function of (Cartx, Carty, Cartz)
        self.Cart_to_xx = [sp.sympify(0)] * 3
        # self.xxSph must be set as a function of (self.xx[0],xx[1],xx[2])
        self.xxSph = [sp.sympify(0)] * 3
        # self.scalefactor_orthog must be set as a function of (self.xx[0],xx[1],xx[2])
        self.scalefactor_orthog = [sp.sympify(0)] * 3
        # UnitVectors must be set as a function of (self.xx[0],xx[1],xx[2])
        self.UnitVectors = ixp.zerorank2(dimension=3)
        # module name for CodeParameters
        self.CodeParam_modulename = f"{__name__}_{CoordSystem}"
        self.add_CodeParams_to_glb_code_params_dict = False
        self.add_CodeParams_to_glb_code_params_dict = (
            self.CoordSystem
            == par.parval_from_str("CoordSystem_to_register_CodeParameters")
        ) or (par.parval_from_str("CoordSystem_to_register_CodeParameters") == "All")

        # START: RFM PRECOMPUTE STUFF
        # Must be set in terms of generic functions of xx[]s
        self.scalefactor_orthog_funcform = [sp.sympify(0)] * 3

        self.f0_of_xx0_funcform = sp.Function("f0_of_xx0_funcform")(self.xx[0])
        self.f1_of_xx1_funcform = sp.Function("f1_of_xx1_funcform")(self.xx[1])
        self.f2_of_xx0_funcform = sp.Function("f2_of_xx0_funcform")(self.xx[0])
        self.f3_of_xx2_funcform = sp.Function("f3_of_xx2_funcform")(self.xx[2])
        self.f4_of_xx1_funcform = sp.Function("f4_of_xx1_funcform")(self.xx[1])
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
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )

        # return values of register_rfm_CodeParameters() in the following code block are unused, so we ignore them.
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["f0_of_xx0__D0", "f0_of_xx0__DD00", "f0_of_xx0__DDD000"],
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["f1_of_xx1__D1", "f1_of_xx1__DD11", "f1_of_xx1__DDD111"],
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["f2_of_xx0__D0", "f2_of_xx0__DD00"],
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["f3_of_xx2__D2", "f3_of_xx2__DD22"],
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        par.register_CodeParameters(
            "REAL",
            self.CodeParam_modulename,
            ["f4_of_xx1__D1", "f4_of_xx1__DD11", "f4_of_xx1__DDD111"],
            add_to_parfile=False,
            add_to_set_CodeParameters_h=False,
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        # END: RFM PRECOMPUTE STUFF

        if "Cartesian" in CoordSystem:
            self.EigenCoord = "Cartesian"
            self.cartesian_like()
        elif "Spherical" in CoordSystem:
            self.EigenCoord = "Spherical"
            self.spherical_like()
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

            Jac_dUCart_dDrfmUD = [[sp.sympify(0) for _ in range(3)] for _ in range(3)]
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
        if not enable_rfm_precompute:
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
        for i in range(3):
            self.detgammahatdD[i] = sp.diff(self.detgammahat, self.xx[i])
            for j in range(3):
                self.detgammahatdDD[i][j] = sp.diff(self.detgammahatdD[i], self.xx[j])

        # Step 3a: Compute 1st & 2nd derivatives of rescaling vectors.
        #          (E.g., needed in BSSN for betaUdDD computation)
        self.ReUdD = ixp.zerorank2(3)
        self.ReUdDD = ixp.zerorank3(3)
        self.ReDdD = ixp.zerorank2(3)
        self.ReDdDD = ixp.zerorank3(3)
        for i in range(3):
            for j in range(3):
                self.ReUdD[i][j] = sp.diff(self.ReU[i], self.xx[j])
                self.ReDdD[i][j] = sp.diff(self.ReD[i], self.xx[j])
                for k in range(3):
                    self.ReUdDD[i][j][k] = sp.diff(self.ReUdD[i][j], self.xx[k])
                    self.ReDdDD[i][j][k] = sp.diff(self.ReDdD[i][j], self.xx[k])

        # Step 3b: Compute 1st & 2nd derivatives of rescaling matrix.
        self.ReDDdD = ixp.zerorank3(3)
        self.ReDDdDD = ixp.zerorank4(3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    self.ReDDdD[i][j][k] = sp.diff(self.ReDD[i][j], self.xx[k])
                    for l in range(3):
                        # Simplifying this doesn't appear to help overall NRPy run time.
                        self.ReDDdDD[i][j][k][l] = sp.diff(
                            self.ReDDdD[i][j][k], self.xx[l]
                        )

        # Step 3c: Compute 1st & 2nd derivatives of reference metric.
        self.ghatDDdD = ixp.zerorank3(3)
        self.ghatDDdDD = ixp.zerorank4(3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    if self.SymPySimplifyExpressions:
                        #                    ghatDDdD[i][j][k] = sp.trigsimp(sp.diff(ghatDD[i][j], xx[k])) # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                        self.ghatDDdD[i][j][k] = cached_simplify(
                            sp.diff(self.ghatDD[i][j], self.xx[k])
                        )  # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                    else:
                        self.ghatDDdD[i][j][k] = sp.diff(
                            self.ghatDD[i][j], self.xx[k]
                        )  # FIXME: BAD: MUST BE SIMPLIFIED OR ANSWER IS INCORRECT! Must be some bug in sympy...
                    for l in range(3):
                        self.ghatDDdDD[i][j][k][l] = sp.diff(
                            self.ghatDDdD[i][j][k], self.xx[l]
                        )

        # Step 4a: Compute Christoffel symbols of reference metric.
        self.GammahatUDD = ixp.zerorank3(3)
        for i in range(3):
            for k in range(3):
                for l in range(3):
                    for m in range(3):
                        #                    GammahatUDD[i][k][l] += sp.trigsimp((sp.Rational(1,2))*ghatUU[i][m]*\
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
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    for l in range(3):
                        self.GammahatUDDdD[i][j][k][l] = sp.diff(
                            self.GammahatUDD[i][j][k], self.xx[l]
                        )

        # Step 4c: If rfm_precompute is disabled, then we are finished with this function.
        #          Otherwise continue to Step 5.
        if not enable_rfm_precompute:
            return

        # Step 5: Now that all hatted quantities are written in terms of generic SymPy functions,
        #         we will now replace SymPy functions with simple variables using rigid NRPy+ syntax,
        #         and store these variables to globals defined above.
        def make_replacements(expr: sp.Expr) -> Any:
            """
            Replace SymPy functions with simple variables using rigid NRPy+ syntax.

            :param expr: SymPy expression to be replaced
            :return: Expression with replaced variables
            """
            sympy_version = sp.__version__.replace("rc", "...").replace(
                "b", "..."
            )  # Ignore the rc's and b's
            # (release candidates & betas).
            sympy_version_decimal = float(
                int(sympy_version.split(".")[0])
                + int(sympy_version.split(".")[1]) / 10.0
            )
            is_old_sympy_version = sympy_version_decimal < 1.2
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
                        var, order = str(item.args[1][0])[2:], item.args[1][1]
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

        # Step 6.a: Construct the full list of *unique* NRPy+ variables representing the
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

        freevars_uniq = superfast_uniq(freevars)

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

            frees_uniq = superfast_uniq(diff_result.free_symbols)
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
                self.xx_to_Cart[ii] = (
                    AMPLXYZ
                    * (sp.exp(self.xx[ii] / SINHWXYZ) - sp.exp(-self.xx[ii] / SINHWXYZ))
                    / (sp.exp(1 / SINHWXYZ) - sp.exp(-1 / SINHWXYZ))
                )

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

        # Set the transpose of the matrix of unit vectors
        self.UnitVectors = [
            [sp.sympify(1), sp.sympify(0), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(1), sp.sympify(0)],
            [sp.sympify(0), sp.sympify(0), sp.sympify(1)],
        ]

    def spherical_like(self) -> None:
        """Initialize class for Spherical-like coordinate systems."""
        M_PI = par.register_CodeParameter(
            "#define",
            self.CodeParam_modulename,
            "M_PI",
            "",
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
        if self.CoordSystem == "Spherical":
            RMAX = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "RMAX",
                10.0,
                add_to_parfile=self.add_rfm_params_to_parfile,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            self.xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            self.xxmax = [RMAX, M_PI, M_PI]
            self.grid_physical_size_dict = {"RMAX": "grid_physical_size"}

            r = self.xx[0]
            th = self.xx[1]
            ph = self.xx[2]

            self.Cart_to_xx[0] = sp.sqrt(
                self.Cartx**2 + self.Carty**2 + self.Cartz**2
            )
            self.Cart_to_xx[1] = sp.acos(self.Cartz / self.Cart_to_xx[0])
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)
        elif self.CoordSystem == "SinhSpherical":
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
            self.xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            self.xxmax = [sp.sympify(1), M_PI, M_PI]
            self.grid_physical_size_dict = {"AMPL": "grid_physical_size"}

            # Set SinhSpherical radial coordinate by default; overwrite later if CoordSystem == "SinhSphericalv2".
            r = (
                AMPL
                * (sp.exp(self.xx[0] / SINHW) - sp.exp(-self.xx[0] / SINHW))
                / (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
            )
            th = self.xx[1]
            ph = self.xx[2]

            rCart = sp.sqrt(self.Cartx**2 + self.Carty**2 + self.Cartz**2)
            self.Cart_to_xx[0] = SINHW * sp.asinh(rCart * sp.sinh(1 / SINHW) / AMPL)
            self.Cart_to_xx[1] = sp.acos(self.Cartz / rCart)
            self.Cart_to_xx[2] = sp.atan2(self.Carty, self.Cartx)

        # SinhSphericalv2 adds the parameter "const_dr", which allows for a region near self.xx[0]=0 to have
        # constant radial resolution of const_dr, provided the sinh() term does not dominate near self.xx[0]=0.
        elif self.CoordSystem == "SinhSphericalv2":
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
            self.xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
            self.xxmax = [sp.sympify(1), M_PI, M_PI]
            self.grid_physical_size_dict = {"AMPL": "grid_physical_size"}

            const_dr = par.register_CodeParameter(
                "REAL",
                self.CodeParam_modulename,
                "const_dr",
                0.0625,
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )
            r = AMPL * (
                const_dr * self.xx[0]
                + (sp.exp(self.xx[0] / SINHW) - sp.exp(-self.xx[0] / SINHW))
                / (sp.exp(1 / SINHW) - sp.exp(-1 / SINHW))
            )
            th = self.xx[1]
            ph = self.xx[2]

            # NO CLOSED-FORM EXPRESSION FOR RADIAL INVERSION.
            # self.Cart_to_xx[0] = "NewtonRaphson"
            # self.Cart_to_xx[1] = sp.acos(Cartz / sp.sqrt(Cartx ** 2 + Carty ** 2 + Cartz ** 2))
            # self.Cart_to_xx[2] = sp.atan2(Carty, Cartx)
        else:
            raise ValueError(
                f"Spherical-like CoordSystem == {self.CoordSystem} unrecognized"
            )

        # BEGIN: Set universal attributes for all spherical-like coordinate systems:
        self.xxSph[0] = r
        self.xxSph[1] = th
        self.xxSph[2] = ph

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
        self.UnitVectors = [[sp.sin(xxS[1])*sp.cos(xxS[2]), sp.sin(xxS[1])*sp.sin(xxS[2]),  sp.cos(xxS[1])],
                            [sp.cos(xxS[1])*sp.cos(xxS[2]), sp.cos(xxS[1])*sp.sin(xxS[2]), -sp.sin(xxS[1])],
                            [              -sp.sin(xxS[2]),                sp.cos(xxS[2]),  sp.sympify(0) ]]
        # fmt: on
        # END: Set universal attributes for all spherical-like coordinate systems:

    def prolate_spheroidal_like(self) -> None:
        """Initialize class for Prolate spheroidal (SymTP)-like coordinate systems."""
        M_PI, M_SQRT1_2 = par.register_CodeParameters(
            "#define",
            self.CodeParam_modulename,
            ["M_PI", "M_SQRT1_2"],
            "",
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )
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

        self.xxmin = [sp.sympify(0), sp.sympify(0), -M_PI]
        self.xxmax = [AMAX, M_PI, M_PI]
        self.grid_physical_size_dict = {"AMAX": "grid_physical_size"}

        AA = self.xx[0]

        if self.CoordSystem == "SinhSymTP":
            self.xxmax[0] = sp.sympify(1)
            # With xxmax[0] = 1, sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and 1.
            #   Then AA = AMAX * sinh(xx0/SINHWAA) / sinh(1/SINHWAA) will evaluate to a number between 0 and AMAX.
            AA = (
                AMAX
                * (sp.exp(self.xx[0] / SINHWAA) - sp.exp(-self.xx[0] / SINHWAA))
                / (sp.exp(1 / SINHWAA) - sp.exp(-1 / SINHWAA))
            )

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
                * M_SQRT1_2
            )  # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

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
                    * M_SQRT1_2
                )
            )  # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            self.Cart_to_xx[2] = phSph

        elif self.CoordSystem == "SinhSymTP":
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
                * M_SQRT1_2
            )  # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

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
                    * M_SQRT1_2
                )
            )  # M_SQRT1_2 = 1/sqrt(2); define this way for UnitTesting

            self.Cart_to_xx[2] = phSph
        else:
            raise ValueError(
                f"prolate-spheroidal-like CoordSystem == {self.CoordSystem} unrecognized"
            )

        # BEGIN: Set universal attributes for all prolate-spheroidal-like coordinate systems:
        self.scalefactor_orthog[0] = sp.diff(AA, self.xx[0]) * var1 / var2
        self.scalefactor_orthog[1] = var1
        self.scalefactor_orthog[2] = AA * sp.sin(self.xx[1])

        self.f0_of_xx0 = AA
        self.f1_of_xx1 = sp.sin(self.xx[1])
        self.f2_of_xx0 = var2
        self.f4_of_xx1 = bScale * sp.sin(self.xx[1])

        # var1 = sp.sqrt(AA**2 + (bScale * sp.sin(self.xx[1])) ** 2)
        #      = sp.sqrt(self.f0_of_xx0_funcform**2 + (bScale * self.f1_of_xx1_funcform) ** 2)

        var1_funcform = sp.sqrt(
            self.f0_of_xx0_funcform**2 + self.f4_of_xx1_funcform**2
        )
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
        """Initialize class for Cylindrical-like coordinate systems."""
        M_PI = par.register_CodeParameter(
            "#define",
            self.CodeParam_modulename,
            "M_PI",
            "",
            add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
        )

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
            self.xxmin = [sp.sympify(0), -M_PI, ZMIN]
            self.xxmax = [RHOMAX, M_PI, ZMAX]
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

            self.xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
            self.xxmax = [sp.sympify(1), M_PI, sp.sympify(+1)]
            self.grid_physical_size_dict = {
                "AMPLRHO": "grid_physical_size",
                "AMPLZ": "grid_physical_size",
            }

            # Set SinhCylindrical radial & z coordinates by default; overwrite later if CoordSystem == "SinhCylindricalv2".
            RHOCYL = (
                AMPLRHO
                * (sp.exp(self.xx[0] / SINHWRHO) - sp.exp(-self.xx[0] / SINHWRHO))
                / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO))
            )
            # phi coordinate remains unchanged.
            PHICYL = self.xx[1]
            ZCYL = (
                AMPLZ
                * (sp.exp(self.xx[2] / SINHWZ) - sp.exp(-self.xx[2] / SINHWZ))
                / (sp.exp(1 / SINHWZ) - sp.exp(-1 / SINHWZ))
            )
            self.Cart_to_xx[0] = SINHWRHO * sp.asinh(
                sp.sqrt(self.Cartx**2 + self.Carty**2)
                * sp.sinh(1 / SINHWRHO)
                / AMPLRHO
            )
            self.Cart_to_xx[1] = sp.atan2(self.Carty, self.Cartx)
            self.Cart_to_xx[2] = SINHWZ * sp.asinh(
                self.Cartz * sp.sinh(1 / SINHWZ) / AMPLZ
            )

        # SinhCylindricalv2 adds the parameters "const_drho", "const_dz", which allows for regions near xx[0]=0
        # and xx[2]=0 to have constant rho and z resolution of const_drho and const_dz, provided the sinh() terms
        # do not dominate near xx[0]=0 and xx[2]=0.
        elif self.CoordSystem == "SinhCylindricalv2":
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
            self.xxmin = [sp.sympify(0), -M_PI, sp.sympify(-1)]
            self.xxmax = [sp.sympify(1), M_PI, sp.sympify(+1)]
            self.grid_physical_size_dict = {
                "AMPLRHO": "grid_physical_size",
                "AMPLZ": "grid_physical_size",
            }
            const_drho, const_dz = par.register_CodeParameters(
                "REAL",
                self.CodeParam_modulename,
                ["const_drho", "const_dz"],
                [0.0625, 0.0625],
                add_to_glb_code_params_dict=self.add_CodeParams_to_glb_code_params_dict,
            )

            RHOCYL = AMPLRHO * (
                const_drho * self.xx[0]
                + (sp.exp(self.xx[0] / SINHWRHO) - sp.exp(-self.xx[0] / SINHWRHO))
                / (sp.exp(1 / SINHWRHO) - sp.exp(-1 / SINHWRHO))
            )
            PHICYL = self.xx[1]
            ZCYL = AMPLZ * (
                const_dz * self.xx[2]
                + (sp.exp(self.xx[2] / SINHWZ) - sp.exp(-self.xx[2] / SINHWZ))
                / (sp.exp(1 / SINHWZ) - sp.exp(-1 / SINHWZ))
            )

            # NO CLOSED-FORM EXPRESSION FOR RADIAL OR Z INVERSION.
            # Cart_to_xx[0] = "NewtonRaphson"
            # Cart_to_xx[1] = sp.atan2(Carty, Cartx)
            # Cart_to_xx[2] = "NewtonRaphson"
        else:
            raise ValueError(
                f"Cylindrical-like CoordSystem == {self.CoordSystem} unrecognized"
            )

        # START: Set universal attributes for all Cylindrical-like coordinate systems:
        self.xx_to_Cart[0] = RHOCYL * sp.cos(PHICYL)
        self.xx_to_Cart[1] = RHOCYL * sp.sin(PHICYL)
        self.xx_to_Cart[2] = ZCYL

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


class rfm_dict(Dict[str, ReferenceMetric]):
    """Custom dictionary for storing ReferenceMetric objects."""

    def __getitem__(self, CoordSystem_in: str) -> ReferenceMetric:
        if CoordSystem_in not in self:
            # In case [CoordSystem]_rfm_precompute is passed:
            CoordSystem = CoordSystem_in.replace("_rfm_precompute", "")
            print(f"Setting up reference metric for CoordSystem = {CoordSystem}.")
            self.__setitem__(
                CoordSystem, ReferenceMetric(CoordSystem, enable_rfm_precompute=False)
            )
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


if __name__ == "__main__":
    import doctest
    import sys
    import os
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
        "SinhSphericalv2",
        "Cartesian",
        "SinhCartesian",
        "Cylindrical",
        "SinhCylindrical",
        "SymTP",
        "SinhSymTP",
    ]:
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

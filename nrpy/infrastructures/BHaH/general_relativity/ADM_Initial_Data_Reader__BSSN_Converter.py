"""
Performs the conversion from ADM spacetime variables
 in Spherical or Cartesian coordinates,
 to BSSN quantities in any basis supported by NRPy+.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""


# Initialize core Python/NRPy+ modules
# Step 1: Initialize core Python/NRPy+ modules
from typing import List, cast, Sequence, Optional
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import nrpy.c_function as cfc  # NRPy+: C function registration
import nrpy.c_codegen as ccg  # NRPy+: C code generation
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
import nrpy.helpers.jacobians as jac

# NRPy+: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN

import nrpy.infrastructures.BHaH.simple_loop as lp


def register_CFunction_exact_ADM_ID_function(
    IDCoordSystem: str,
    IDtype: str,
    alpha: sp.Expr,
    betaU: List[sp.Expr],
    BU: List[sp.Expr],
    gammaDD: List[List[sp.Expr]],
    KDD: List[List[sp.Expr]],
) -> None:
    """
    Adds a C function for exact ADM initial data of a given ID type.

    :param IDCoordSystem: The ID coordinate system.
    :param IDtype: The ID type.
    :param alpha: The lapse function.
    :param betaU: The beta upper indices.
    :param BU: partial_t beta upper indices.
    :param gammaDD: The 3-metric with lower indices.
    :param KDD: The extrinsic curvature with lower indices.
    """

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"{IDtype} initial data"
    c_type = "void"
    name = IDtype
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"""
    body = ""

    if IDCoordSystem == "Spherical":
        body += r"""
  REAL xx0,xx1,xx2 __attribute__((unused));
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
  const REAL r  = xx0;
  const REAL th = xx1;
  const REAL ph = xx2;
"""
    elif IDCoordSystem == "Cartesian":
        body += r"""  const REAL x=xCart[0], y=xCart[1], z=xCart[2];
"""
    else:
        raise ValueError(
            f"register_CFunction_exact_ADM_ID_function() Error: IDCoordSystem == {IDCoordSystem} unsupported."
        )

    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["initial_data->alpha"]

    for i in range(3):
        list_of_output_exprs += [betaU[i]]
        list_of_output_varnames += [f"initial_data->betaSphorCartU{i}"]
        list_of_output_exprs += [BU[i]]
        list_of_output_varnames += [f"initial_data->BSphorCartU{i}"]

        for j in range(i, 3):
            list_of_output_exprs += [gammaDD[i][j]]
            list_of_output_varnames += [f"initial_data->gammaSphorCartDD{i}{j}"]
            list_of_output_exprs += [KDD[i][j]]
            list_of_output_varnames += [f"initial_data->KSphorCartDD{i}{j}"]

    # Sort the outputs before calling outputC()
    list_of_output_varnames, list_of_output_exprs = (
        list(t)
        for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs)))
    )

    body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )

    cfc.register_CFunction(
        includes=includes,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    # return pickle_NRPy_env()


def Cfunction_ADM_SphorCart_to_Cart(
    IDCoordSystem: str = "Spherical",
    include_T4UU: bool = False,
) -> str:
    """
    Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis.

    :param IDCoordSystem: The input coordinate system. Defaults to "Spherical".
    :param include_T4UU: Whether to include the stress-energy tensor. Defaults to False.

    :return: The name of the generated C function.
    """

    desc = "Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis"
    c_type = "static void"
    name = "ADM_SphorCart_to_Cart"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xCart[3], const initial_data_struct *restrict initial_data,
    ADM_Cart_basis_struct *restrict ADM_Cart_basis"""

    body = r"""// Unpack initial_data for ADM vectors/tensors
"""
    for var in ["betaSphorCartU", "BSphorCartU"]:
        for j in range(3):
            body += f"const REAL {var}{j} = initial_data->{var}{j};\n"
        body += "\n"
    for var in ["gammaSphorCartDD", "KSphorCartDD"]:
        for j in range(3):
            for k in range(j, 3):
                body += f"  const REAL {var}{j}{k} = initial_data->{var}{j}{k};\n"
        body += "\n"
    # Read stress-energy tensor in spherical or Cartesian basis if desired.
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                body += f"const REAL T4SphorCartUU{mu}{nu} = initial_data->T4SphorCartUU{mu}{nu};\n"
        body += "\n"

    # Set reference_metric to the IDCoordSystem
    rfm = refmetric.reference_metric[IDCoordSystem]

    # Cartesian -> Cartesian does not require setting xx0,xx1,xx2, as all Jacobians are just multiplying by 1
    if IDCoordSystem != "Cartesian":
        body += rf"""
      // Perform the basis transform on ADM vectors/tensors from {IDCoordSystem} to Cartesian:

      // Set destination xx[3] based on desired xCart[3]
      REAL xx0,xx1,xx2;
      """ + ccg.c_codegen(
            rfm.Cart_to_xx,
            ["xx0", "xx1", "xx2"],
            include_braces=True,
        ).replace(
            "Cartx", "xCart[0]"
        ).replace(
            "Carty", "xCart[1]"
        ).replace(
            "Cartz", "xCart[2]"
        )

    # Define the input variables:
    gammaSphorCartDD = ixp.declarerank2("gammaSphorCartDD", symmetry="sym01")
    KSphorCartDD = ixp.declarerank2("KSphorCartDD", symmetry="sym01")
    betaSphorCartU = ixp.declarerank1("betaSphorCartU")
    BSphorCartU = ixp.declarerank1("BSphorCartU")
    T4SphorCartUU = ixp.declarerank2("T4SphorCartUU", symmetry="sym01", dimension=4)

    # Compute Jacobian to convert to Cartesian coordinates
    gammaCartDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        IDCoordSystem, gammaSphorCartDD
    )

    KCartDD = jac.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        IDCoordSystem, KSphorCartDD
    )
    betaCartU = jac.basis_transform_vectorU_from_rfmbasis_to_Cartesian(
        IDCoordSystem, betaSphorCartU
    )
    BCartU = jac.basis_transform_vectorU_from_rfmbasis_to_Cartesian(
        IDCoordSystem, BSphorCartU
    )
    T4CartUU = cast(Sequence[Sequence[sp.Expr]], ixp.zerorank2(dimension=4))
    if include_T4UU:
        T4CartUU = jac.basis_transform_4tensorUU_from_time_indep_rfmbasis_to_Cartesian(
            IDCoordSystem, T4SphorCartUU
        )

    alpha = sp.symbols("initial_data->alpha", real=True)
    list_of_output_exprs = [alpha]
    list_of_output_varnames = ["ADM_Cart_basis->alpha"]
    for i in range(3):
        list_of_output_exprs += [betaCartU[i]]
        list_of_output_varnames += [f"ADM_Cart_basis->betaU{i}"]
        list_of_output_exprs += [BCartU[i]]
        list_of_output_varnames += [f"ADM_Cart_basis->BU{i}"]
        for j in range(i, 3):
            list_of_output_exprs += [gammaCartDD[i][j]]
            list_of_output_varnames += [f"ADM_Cart_basis->gammaDD{i}{j}"]
            list_of_output_exprs += [KCartDD[i][j]]
            list_of_output_varnames += [f"ADM_Cart_basis->KDD{i}{j}"]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += [f"ADM_Cart_basis->T4UU{mu}{nu}"]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (
        list(t)
        for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs)))
    )

    body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )

    return cfc.CFunction(
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


def Cfunction_ADM_Cart_to_BSSN_Cart(include_T4UU: bool = False) -> str:
    """
    Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.

    :param include_T4UU: Boolean flag to indicate whether to include T4UU or not.

    :return: A string representing the full C function.
    """
    desc = "Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis"
    c_type = "static void"
    name = "ADM_Cart_to_BSSN_Cart"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xCart[3], const ADM_Cart_basis_struct *restrict ADM_Cart_basis,
    BSSN_Cart_basis_struct *restrict BSSN_Cart_basis"""

    gammaCartDD = ixp.declarerank2("ADM_Cart_basis->gammaDD", symmetry="sym01")
    KCartDD = ixp.declarerank2("ADM_Cart_basis->KDD", symmetry="sym01")

    dummy1U = ixp.zerorank1()
    dummy2U = ixp.zerorank1()
    adm2bssn = ADM_to_BSSN(gammaCartDD, KCartDD, dummy1U, dummy2U, "Cartesian")

    body = r"""
  // *In the Cartesian basis*, convert ADM quantities gammaDD & KDD
  //   into BSSN gammabarDD, AbarDD, cf, and trK.
  BSSN_Cart_basis->alpha = ADM_Cart_basis->alpha;
  BSSN_Cart_basis->betaU0 = ADM_Cart_basis->betaU0;
  BSSN_Cart_basis->betaU1 = ADM_Cart_basis->betaU1;
  BSSN_Cart_basis->betaU2 = ADM_Cart_basis->betaU2;
  BSSN_Cart_basis->BU0 = ADM_Cart_basis->BU0;
  BSSN_Cart_basis->BU1 = ADM_Cart_basis->BU1;
  BSSN_Cart_basis->BU2 = ADM_Cart_basis->BU2;
"""
    list_of_output_exprs = [adm2bssn.cf, adm2bssn.trK]
    list_of_output_varnames = ["BSSN_Cart_basis->cf", "BSSN_Cart_basis->trK"]
    for i in range(3):
        for j in range(i, 3):
            list_of_output_exprs += [adm2bssn.gammabarDD[i][j]]
            list_of_output_varnames += [f"BSSN_Cart_basis->gammabarDD{i}{j}"]
            list_of_output_exprs += [adm2bssn.AbarDD[i][j]]
            list_of_output_varnames += [f"BSSN_Cart_basis->AbarDD{i}{j}"]
    if include_T4UU:
        T4CartUU = ixp.declarerank2(
            "ADM_Cart_basis->T4UU", symmetry="sym01", dimension=4
        )
        for mu in range(4):
            for nu in range(mu, 4):
                list_of_output_exprs += [T4CartUU[mu][nu]]
                list_of_output_varnames += [f"BSSN_Cart_basis->T4UU{mu}{nu}"]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (
        list(t)
        for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs)))
    )
    body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )

    return cfc.CFunction(
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


def Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(
    CoordSystem: str, include_T4UU: bool = False
) -> str:
    """
    Convert Cartesian-basis BSSN vectors/tensors *except* lambda^i,
    to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities.

    :param CoordSystem: Coordinate system to which the variables are to be transformed
    :param include_T4UU: Whether to include T4UU tensor in the transformation
    :return: Returns the generated C code as a string.
    """

    desc = r"""Convert Cartesian-basis BSSN vectors/tensors *except* lambda^i,
to the basis specified by `reference_metric::CoordSystem`, then rescale these BSSN quantities"""
    c_type = "static void"
    name = "BSSN_Cart_to_rescaled_BSSN_rfm"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis"""

    body = r"""
  REAL xx0,xx1,xx2 __attribute__((unused));  // xx2 might be unused in the case of axisymmetric initial data.
  {
    int unused_Cart_to_i0i1i2[3];
    REAL xx[3];
    Cart_to_xx_and_nearest_i0i1i2(commondata, params, xCart, xx, unused_Cart_to_i0i1i2);
    xx0=xx[0];  xx1=xx[1];  xx2=xx[2];
  }
"""

    # Define the input variables:
    gammabarCartDD = ixp.declarerank2("BSSN_Cart_basis->gammabarDD", symmetry="sym01")
    AbarCartDD = ixp.declarerank2("BSSN_Cart_basis->AbarDD", symmetry="sym01")
    betaCartU = ixp.declarerank1("BSSN_Cart_basis->betaU")
    BCartU = ixp.declarerank1("BSSN_Cart_basis->BU")

    # Compute Jacobian to convert to Cartesian coordinates
    gammabarDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        CoordSystem, gammabarCartDD
    )
    AbarDD = jac.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        CoordSystem, AbarCartDD
    )
    betaU = jac.basis_transform_vectorU_from_Cartesian_to_rfmbasis(
        CoordSystem, betaCartU
    )
    BU = jac.basis_transform_vectorU_from_Cartesian_to_rfmbasis(CoordSystem, BCartU)
    if include_T4UU:
        T4CartUU = ixp.declarerank2(
            "BSSN_Cart_basis->T4UU", symmetry="sym01", dimension=4
        )
        T4UU = jac.basis_transform_4tensorUU_from_Cartesian_to_time_indep_rfmbasis(
            CoordSystem, T4CartUU
        )

    # Next rescale:
    rfm = refmetric.reference_metric[CoordSystem]
    vetU = ixp.zerorank1()
    betU = ixp.zerorank1()
    hDD = ixp.zerorank2()
    aDD = ixp.zerorank2()
    for i in range(3):
        vetU[i] = betaU[i] / rfm.ReU[i]
        betU[i] = BU[i] / rfm.ReU[i]
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]
            aDD[i][j] = AbarDD[i][j] / rfm.ReDD[i][j]

    alpha, cf, trK = sp.symbols(
        "BSSN_Cart_basis->alpha BSSN_Cart_basis->cf BSSN_Cart_basis->trK", real=True
    )

    list_of_output_exprs = [alpha, cf, trK]
    list_of_output_varnames = [
        "rescaled_BSSN_rfm_basis->alpha",
        "rescaled_BSSN_rfm_basis->cf",
        "rescaled_BSSN_rfm_basis->trK",
    ]
    for i in range(3):
        list_of_output_exprs += [vetU[i]]
        list_of_output_varnames += [f"rescaled_BSSN_rfm_basis->vetU{i}"]
        list_of_output_exprs += [betU[i]]
        list_of_output_varnames += [f"rescaled_BSSN_rfm_basis->betU{i}"]
        for j in range(i, 3):
            list_of_output_exprs += [hDD[i][j]]
            list_of_output_varnames += [f"rescaled_BSSN_rfm_basis->hDD{i}{j}"]
            list_of_output_exprs += [aDD[i][j]]
            list_of_output_varnames += [f"rescaled_BSSN_rfm_basis->aDD{i}{j}"]
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                # T4UU IS ASSUMED NOT RESCALED; RESCALINGS ARE HANDLED WITHIN BSSN RHSs, etc.
                list_of_output_exprs += [T4UU[mu][nu]]
                list_of_output_varnames += [f"rescaled_BSSN_rfm_basis->T4UU{mu}{nu}"]

    # Sort the outputs before calling outputC()
    # https://stackoverflow.com/questions/9764298/is-it-possible-to-sort-two-listswhich-reference-each-other-in-the-exact-same-w
    list_of_output_varnames, list_of_output_exprs = (
        list(t)
        for t in zip(*sorted(zip(list_of_output_varnames, list_of_output_exprs)))
    )

    body += ccg.c_codegen(
        list_of_output_exprs,
        list_of_output_varnames,
        verbose=False,
        include_braces=False,
    )

    return cfc.CFunction(
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


# Cfunction_initial_data_lambdaU_grid_interior() computes lambdaU from
# finite-difference derivatives of rescaled metric quantities
def Cfunction_initial_data_lambdaU_grid_interior(CoordSystem: str) -> str:
    """
    Compute lambdaU in the specified coordinate system.

    :param CoordSystem: The coordinate system to be used.
    :return: The full function generated for computing lambdaU.
    """

    c_type = "static void"

    desc = f"Compute lambdaU in {CoordSystem} coordinates"
    name = "initial_data_lambdaU_grid_interior"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], REAL *restrict in_gfs"""
    # Step 7: Compute $\bar{\Lambda}^i$ from finite-difference derivatives of rescaled metric quantities

    # We will need all BSSN gridfunctions to be defined, as well as
    # expressions for gammabarDD_dD in terms of exact derivatives of
    # the rescaling matrix and finite-difference derivatives of
    # hDD's. This functionality is provided by BSSN.BSSN_unrescaled_and_barred_vars,
    # which we call here to overwrite above definitions of gammabarDD,gammabarUU, etc.
    Bq = BSSN_quantities[CoordSystem]
    gammabarUU = Bq.gammabarUU
    GammabarUDD = Bq.GammabarUDD

    # Next evaluate \bar{\Lambda}^i, based on GammabarUDD above and GammahatUDD
    # (from the reference metric):
    rfm = refmetric.reference_metric[CoordSystem]
    LambdabarU = ixp.zerorank1()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                LambdabarU[i] += gammabarUU[j][k] * (
                    GammabarUDD[i][j][k] - rfm.GammahatUDD[i][j][k]
                )

    # Finally apply rescaling:
    # lambda^i = Lambdabar^i/\text{ReU[i]}
    lambdaU = ixp.zerorank1()
    for i in range(3):
        lambdaU[i] = LambdabarU[i] / rfm.ReU[i]

    body = lp.simple_loop(
        ccg.c_codegen(
            lambdaU,
            [
                gri.BHaHGridFunction.access_gf("lambdaU0"),
                gri.BHaHGridFunction.access_gf("lambdaU1"),
                gri.BHaHGridFunction.access_gf("lambdaU2"),
            ],
            verbose=False,
            include_braces=False,
            enable_fd_codegen=True,
        ),
        loop_region="interior",
        read_xxs=True,
    )

    return cfc.CFunction(
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


def register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
    CoordSystem: str,
    addl_includes: Optional[List[str]] = None,
    IDCoordSystem: str = "Spherical",
    include_T4UU: bool = False,
    enable_fd_functions: bool = False,
) -> None:
    """
    Register the CFunction for converting initial ADM data to BSSN variables.

    :param CoordSystem: Coordinate system for output BSSN variables.
    :param addl_includes: Additional header files to include.
    :param IDCoordSystem: Coordinate system for input ADM variables. Defaults to "Spherical".
    :param include_T4UU: Whether to include stress-energy tensor components.
    :param enable_fd_functions: Whether to enable finite-difference functions.

    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_fd_functions:
        includes += ["finite_difference_functions.h"]

    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            raise ValueError("Error: addl_includes must be a list.")
        includes += addl_includes

    def T4UU_prettyprint() -> str:
        """
        Returns a pretty-printed string for T4UU variables in C code.

        :return: A string containing the C declarations for T4UU variables.
        """
        return r"""
  REAL T4UU00,T4UU01,T4UU02,T4UU03;
  REAL        T4UU11,T4UU12,T4UU13;
  REAL               T4UU22,T4UU23;
  REAL                      T4UU33;
"""

    prefunc = """
// ADM variables in the Cartesian basis:
typedef struct __ADM_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL gammaDD00,gammaDD01,gammaDD02,gammaDD11,gammaDD12,gammaDD22;
  REAL KDD00,KDD01,KDD02,KDD11,KDD12,KDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} ADM_Cart_basis_struct;\n"
    ##############
    prefunc += """
// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL cf, trK;
  REAL gammabarDD00,gammabarDD01,gammabarDD02,gammabarDD11,gammabarDD12,gammabarDD22;
  REAL AbarDD00,AbarDD01,AbarDD02,AbarDD11,AbarDD12,AbarDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} BSSN_Cart_basis_struct;\n"
    ##############
    prefunc += """
// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0,vetU1,vetU2, betU0,betU1,betU2;
  REAL cf, trK;
  REAL hDD00,hDD01,hDD02,hDD11,hDD12,hDD22;
  REAL aDD00,aDD01,aDD02,aDD11,aDD12,aDD22;
"""
    if include_T4UU:
        prefunc += T4UU_prettyprint()
    prefunc += "} rescaled_BSSN_rfm_basis_struct;\n"
    ##############
    ##############
    prefunc += Cfunction_ADM_SphorCart_to_Cart(
        IDCoordSystem=IDCoordSystem, include_T4UU=include_T4UU
    )
    prefunc += Cfunction_ADM_Cart_to_BSSN_Cart(include_T4UU=include_T4UU)
    prefunc += Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(
        CoordSystem=CoordSystem, include_T4UU=include_T4UU
    )
    prefunc += Cfunction_initial_data_lambdaU_grid_interior(CoordSystem=CoordSystem)

    desc = f"Read ADM data in the {IDCoordSystem} basis, and output rescaled BSSN data in the {CoordSystem} basis"
    c_type = "void"
    name = f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN"
    params = """const commondata_struct *restrict commondata, griddata_struct *restrict griddata, ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist,
                     initial_data_struct *restrict initial_data)"""

    body = r"""
  const int Nxx_plus_2NGHOSTS0 = griddata->params.Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = griddata->params.Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = griddata->params.Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for", i0,0,Nxx_plus_2NGHOSTS0, i1,0,Nxx_plus_2NGHOSTS1, i2,0,Nxx_plus_2NGHOSTS2) {
    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];  xx_to_Cart(commondata, &griddata->params, griddata->xx, i0,i1,i2, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(commondata, &griddata->params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(commondata, &griddata->params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(commondata, &griddata->params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(commondata, &griddata->params, xCart, &BSSN_Cart_basis, &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3(i0,i1,i2);

    // Output data to gridfunctions
"""
    gf_list = ["alpha", "trK", "cf"]
    for i in range(3):
        gf_list += [f"vetU{i}", f"betU{i}"]
        for j in range(i, 3):
            gf_list += [f"hDD{i}{j}", f"aDD{i}{j}"]
    for gf in sorted(gf_list):
        body += f"griddata->gridfuncs.y_n_gfs[IDX4pt({gf.upper()}GF, idx3)] = rescaled_BSSN_rfm_basis.{gf};\n"
    if include_T4UU:
        for mu in range(4):
            for nu in range(mu, 4):
                gf = f"T4UU{mu}{nu}"
                body += f"griddata->gridfuncs.auxevol_gfs[IDX4pt({gf.upper()}GF, idx3)] = rescaled_BSSN_rfm_basis.{gf};\n"
    body += """
  } // END LOOP over all gridpoints on given grid

  initial_data_lambdaU_grid_interior(commondata, &griddata->params, griddata->xx, griddata->gridfuncs.y_n_gfs);
"""

    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    # return pickle_NRPy_env()


def generate_BHaH_defines_contribution(
    ID_persist_struct_contents_str: str = "", include_T4UU: bool = False
) -> str:
    """
    Generates C code for the BHaH initial data definitions and contributions.

    :param ID_persist_struct_contents_str: String containing the content of the ID_persist_struct structure.
    :param include_T4UU: Boolean indicating whether to include the T4UU terms.
    :return: A string containing the generated C code for the initial data structures.

    >>> generate_BHaH_defines_contribution("", False).startswith("typedef struct __initial_data_struct__ {")
    True
    >>> "T4SphorCartUU00" in generate_BHaH_defines_contribution("", True)
    True
    """

    Nbd = r"""typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
"""
    if include_T4UU:
        Nbd += """
  REAL T4SphorCartUU00,T4SphorCartUU01,T4SphorCartUU02,T4SphorCartUU03;
  REAL                 T4SphorCartUU11,T4SphorCartUU12,T4SphorCartUU13;
  REAL                                 T4SphorCartUU22,T4SphorCartUU23;
  REAL                                                 T4SphorCartUU33;
"""
    Nbd += """
} initial_data_struct;
"""

    Nbd += "typedef struct __ID_persist_struct__ {\n"
    Nbd += ID_persist_struct_contents_str + "\n"
    Nbd += "} ID_persist_struct;\n"
    return Nbd


if __name__ == "__main__":
    import sys
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

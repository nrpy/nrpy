"""
Perform the conversion from ADM quantities in spherical/Cartesian coordinates, to BSSN quantities in any NRPy-supported basis.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

# Initialize core Python/NRPy modules
# Step 1: Initialize core Python/NRPy modules
from typing import List, Optional, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.c_codegen as ccg  # NRPy: C code generation
import nrpy.c_function as cfc  # NRPy: C function registration
import nrpy.grid as gri  # NRPy: Functions having to do with numerical grids
import nrpy.helpers.jacobians as jac
import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.params as par  # NRPy: Parameter interface
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN

# NRPy: Computes useful BSSN quantities; e.g., gammabarUU & GammabarUDD needed below
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)
from nrpy.helpers.parallelization import utilities
from nrpy.infrastructures import BHaH


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
    Register C function for exact ADM initial data of a given ID type.

    :param IDCoordSystem: The ID coordinate system.
    :param IDtype: The ID type.
    :param alpha: The lapse function.
    :param betaU: The beta upper indices.
    :param BU: partial_t beta upper indices.
    :param gammaDD: The 3-metric with lower indices.
    :param KDD: The extrinsic curvature with lower indices.

    :raises ValueError: If an unsupported coordinate system is specified, ensuring that the function generation is restricted to supported coordinate systems.
    """
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = f"{IDtype} initial data"
    cfunc_type = "void"
    name = IDtype
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL xCart[3], const ID_persist_struct *restrict ID_persist, initial_data_struct *restrict initial_data"""
    body = ""

    if IDCoordSystem == "Spherical":
        rfm_Spherical = refmetric.reference_metric["Spherical"]
        body += "const REAL Cartx=xCart[0],Carty=xCart[1],Cartz=xCart[2];\n"
        body += "REAL r, th, ph;\n"
        body += ccg.c_codegen(
            rfm_Spherical.Cart_to_xx,
            ["r", "th", "ph"],
            verbose=False,
            include_braces=True,
        )
        body += "MAYBE_UNUSED const REAL xx0=r, xx1=th, xx2=ph;\n"
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
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )


def Cfunction_ADM_SphorCart_to_Cart(
    IDCoordSystem: str = "Spherical",
    enable_T4munu: bool = False,
) -> str:
    """
    Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis.

    :param IDCoordSystem: The input coordinate system. Defaults to "Spherical".
    :param enable_T4munu: Whether to include the stress-energy tensor. Defaults to False.

    :return: The name of the generated C function.
    """
    desc = "Convert ADM variables from the spherical or Cartesian basis to the Cartesian basis"
    cfunc_type = "static void"
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
    if enable_T4munu:
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
        ).replace("Cartx", "xCart[0]").replace("Carty", "xCart[1]").replace(
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
    T4CartUU = ixp.zerorank2(dimension=4)
    if enable_T4munu:
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
    if enable_T4munu:
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
        subdirectory=IDCoordSystem,  # Probably not needed
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


def Cfunction_ADM_Cart_to_BSSN_Cart(
    enable_T4munu: bool = False,
) -> str:
    """
    Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis.

    :param enable_T4munu: Boolean flag to indicate whether to include T4UU or not.

    :return: A string representing the full C function.
    """
    desc = "Convert ADM variables in the Cartesian basis to BSSN variables in the Cartesian basis"
    cfunc_type = "static void"
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
    list_of_output_exprs: List[sp.Expr] = [adm2bssn.cf, adm2bssn.trK]
    list_of_output_varnames = ["BSSN_Cart_basis->cf", "BSSN_Cart_basis->trK"]
    for i in range(3):
        for j in range(i, 3):
            list_of_output_exprs += [adm2bssn.gammabarDD[i][j]]
            list_of_output_varnames += [f"BSSN_Cart_basis->gammabarDD{i}{j}"]
            list_of_output_exprs += [adm2bssn.AbarDD[i][j]]
            list_of_output_varnames += [f"BSSN_Cart_basis->AbarDD{i}{j}"]
    if enable_T4munu:
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
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
    ).full_function


def Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(
    CoordSystem: str,
    enable_T4munu: bool = False,
) -> str:
    """
    Convert Cartesian-basis BSSN vectors/tensors (except lambda^i) to CoordSystem basis, then rescale these BSSN quantities.

    :param CoordSystem: Coordinate system to which the variables are to be transformed.
    :param enable_T4munu: Whether to include T4UU tensor in the transformation.
    :return: Returns the generated C code as a string.
    """
    desc = rf"""Cartesian -> {CoordSystem} basis transformation of BSSN vectors/tensors *except* lambda^i.
After the basis transform, all BSSN quantities are rescaled."""
    cfunc_type = "static void"
    name = "BSSN_Cart_to_rescaled_BSSN_rfm"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xxL[3],
                                           const BSSN_Cart_basis_struct *restrict BSSN_Cart_basis,
                                           rescaled_BSSN_rfm_basis_struct *restrict rescaled_BSSN_rfm_basis"""

    body = ""
    if CoordSystem != "Cartesian":
        body += r"""
  const REAL xx0=xxL[0], xx1=xxL[1], xx2=xxL[2];
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
    if enable_T4munu:
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
    if enable_T4munu:
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
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


# Cfunction_initial_data_lambdaU_grid_interior() computes lambdaU from
# finite-difference derivatives of rescaled metric quantities
def Cfunction_initial_data_lambdaU_grid_interior(
    CoordSystem: str,
) -> Tuple[str, str]:
    """
    Compute lambdaU in the specified coordinate system.

    :param CoordSystem: The coordinate system to be used.
    :return: The full function generated for computing lambdaU.
    """
    cfunc_type = "static void"

    desc = f"Compute lambdaU in {CoordSystem} coordinates"
    name = "initial_data_lambdaU_grid_interior"
    parallelization = par.parval_from_str("parallelization")
    arg_dict_cuda = {
        "x0": "const REAL *restrict",
        "x1": "const REAL *restrict",
        "x2": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        "xx[3]": "const REAL *restrict",
        "in_gfs": "REAL *restrict",
    }
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

    enable_intrinsics = (
        False  # FAILS FOR float: True if parallelization != "openmp" else False
    )
    kernel_body = BHaH.simple_loop.simple_loop(
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
            enable_simd=enable_intrinsics,
        ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD"),
        loop_region="interior",
        read_xxs=True,
        enable_intrinsics=enable_intrinsics,
    )
    loop_params = utilities.get_loop_parameters(
        parallelization, enable_intrinsics=enable_intrinsics
    )
    param_symbols, _ = get_params_commondata_symbols_from_expr_list(lambdaU)
    params_definitions = generate_definition_header(
        param_symbols,
        enable_intrinsics=enable_intrinsics,
        var_access=utilities.get_params_access(parallelization),
    ).replace("SIMD", "CUDA" if parallelization == "cuda" else "SIMD")

    kernel_body = f"{loop_params}\n{params_definitions}\n{kernel_body}"

    prefunc, launch_body = utilities.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=desc,
        cfunc_type=cfunc_type,
        launchblock_with_braces=False,
        thread_tiling_macro_suffix="ID_LAMBDAU",
    )

    launch_body = launch_body.replace("in_gfs", "gridfuncs->y_n_gfs")
    for i in range(3):
        launch_body = (
            launch_body.replace(f"x{i},", f"xx[{i}],")
            if parallelization == "cuda"
            else launch_body.replace("xx[3],", "xx,")
        )

    return prefunc, launch_body


def register_BHaH_defines_h(
    ID_persist_struct_str: str,
    enable_T4munu: bool = False,
) -> None:
    """
    Register the initial_data_struct and ID_persist_struct contribution to BHaH_defines.h.

    This function constructs the C typedefs for initial_data_struct—
    including optional stress-energy (T4) fields—and for ID_persist_struct,
    then registers them via register_BHaH_defines.

    :param ID_persist_struct_str: body of the ID_persist_struct (C code fragment)
    :param enable_T4munu:        if True, include the T4SphorCartUU block
    """
    # build the initial_data_struct
    BHd_str = r"""typedef struct __initial_data_struct__ {
  REAL alpha;

  REAL betaSphorCartU0, betaSphorCartU1, betaSphorCartU2;
  REAL BSphorCartU0, BSphorCartU1, BSphorCartU2;

  REAL gammaSphorCartDD00, gammaSphorCartDD01, gammaSphorCartDD02;
  REAL gammaSphorCartDD11, gammaSphorCartDD12, gammaSphorCartDD22;

  REAL KSphorCartDD00, KSphorCartDD01, KSphorCartDD02;
  REAL KSphorCartDD11, KSphorCartDD12, KSphorCartDD22;
"""
    if enable_T4munu:
        BHd_str += r"""
  REAL T4SphorCartUU00,T4SphorCartUU01,T4SphorCartUU02,T4SphorCartUU03;
  REAL                 T4SphorCartUU11,T4SphorCartUU12,T4SphorCartUU13;
  REAL                                 T4SphorCartUU22,T4SphorCartUU23;
  REAL                                                 T4SphorCartUU33;
"""

    BHd_str += r"""
} initial_data_struct;
"""

    # append the ID_persist_struct definition
    BHd_str += "typedef struct __ID_persist_struct__ {\n"
    BHd_str += ID_persist_struct_str + "\n"
    BHd_str += "} ID_persist_struct;\n"

    # register into BHaH_defines.h
    BHaH.BHaH_defines_h.register_BHaH_defines(__name__, BHd_str)


def generate_ADM_Initial_Data_Reader_prefunc_and_lambdaU_launch(
    enable_T4munu: bool, CoordSystem: str, IDCoordSystem: str = "Spherical"
) -> Tuple[str, str]:
    """
    Generate the C “prefunc” string and the lambdaU launch snippet for the initial-data reader converting ADM→BSSN.

    :param enable_T4munu: whether to include the T4UU stress-energy block
    :param CoordSystem:   coordinate system for the BSSN conversion CFunctions
    :param IDCoordSystem: coordinate system for the ADM→Cart conversion (default "Spherical")
    :return: a tuple (prefunc, lambdaU_launch) where
             - prefunc is the C code defining the ADM, BSSN, and rescaled-BSSN structs
               and the conversion CFunction bodies plus lambdaU prefuc
             - lambdaU_launch is the C code snippet to launch the lambdaU grid interior
    """

    def T4UU_prettyprint() -> str:
        """
        Return a pretty-printed string for T4UU variables in C code.

        :return: The pretty-printed T4UU variables declaration block.
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
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} ADM_Cart_basis_struct;\n"

    prefunc += """
// BSSN variables in the Cartesian basis:
typedef struct __BSSN_Cart_basis_struct__ {
  REAL alpha, betaU0,betaU1,betaU2, BU0,BU1,BU2;
  REAL cf, trK;
  REAL gammabarDD00,gammabarDD01,gammabarDD02,gammabarDD11,gammabarDD12,gammabarDD22;
  REAL AbarDD00,AbarDD01,AbarDD02,AbarDD11,AbarDD12,AbarDD22;
"""
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} BSSN_Cart_basis_struct;\n"

    prefunc += """
// Rescaled BSSN variables in the rfm basis:
typedef struct __rescaled_BSSN_rfm_basis_struct__ {
  REAL alpha, vetU0,vetU1,vetU2, betU0,betU1,betU2;
  REAL cf, trK;
  REAL hDD00,hDD01,hDD02,hDD11,hDD12,hDD22;
  REAL aDD00,aDD01,aDD02,aDD11,aDD12,aDD22;
"""
    if enable_T4munu:
        prefunc += T4UU_prettyprint()
    prefunc += "} rescaled_BSSN_rfm_basis_struct;\n"

    prefunc += Cfunction_ADM_SphorCart_to_Cart(
        IDCoordSystem=IDCoordSystem,
        enable_T4munu=enable_T4munu,
    )
    prefunc += Cfunction_ADM_Cart_to_BSSN_Cart(enable_T4munu=enable_T4munu)
    prefunc += Cfunction_BSSN_Cart_to_rescaled_BSSN_rfm(
        CoordSystem=CoordSystem,
        enable_T4munu=enable_T4munu,
    )

    lambdaU_prefunc, lambdaU_launch = Cfunction_initial_data_lambdaU_grid_interior(
        CoordSystem=CoordSystem
    )
    prefunc += lambdaU_prefunc

    return prefunc, lambdaU_launch


def setup_ADM_initial_data_reader(
    ID_persist_struct_str: str,
    enable_T4munu: bool,
    enable_fd_functions: bool,
    addl_includes: Optional[List[str]],
    CoordSystem: str,
    IDCoordSystem: str = "Spherical",
) -> Tuple[List[str], str, str]:
    """
    Perform Steps 1–3 for the ADM initial-data reader registration.

    1. Register the initial_data_struct and ID_persist_struct via BHaH_defines_h.
    2. Assemble the include list.
    3. Generate the prefunc string and lambdaU launch snippet.

    :param ID_persist_struct_str:  String for persistent ID structure.
    :param enable_T4munu:         Whether to include T4UU blocks.
    :param enable_fd_functions:   Whether to add finite-difference headers.
    :param addl_includes:         Additional headers to include.
    :param CoordSystem:           Target coordinate system for CFunctions.
    :param IDCoordSystem:         Input ADM coordinate system (default "Spherical").
    :return:                       A tuple (includes, prefunc, lambdaU_launch).
    :raises ValueError:           If `addl_includes` is provided but not a list.
    """
    # Step 1: register BHaH_defines.h contribution
    register_BHaH_defines_h(
        ID_persist_struct_str=ID_persist_struct_str,
        enable_T4munu=enable_T4munu,
    )

    # Step 2: build include list
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_fd_functions:
        includes += ["finite_difference_functions.h"]

    if addl_includes is not None:
        if not isinstance(addl_includes, list):
            raise ValueError("Error: addl_includes must be a list.")
        includes += addl_includes

    # Step 3: generate prefunc and lambdaU_launch
    prefunc, lambdaU_launch = (
        generate_ADM_Initial_Data_Reader_prefunc_and_lambdaU_launch(
            enable_T4munu=enable_T4munu,
            CoordSystem=CoordSystem,
            IDCoordSystem=IDCoordSystem,
        )
    )

    return includes, prefunc, lambdaU_launch


def build_initial_data_conversion_loop(enable_T4munu: bool) -> str:
    """
    Generate the string for the initial data conversion loop.

      1) Declare the three Nxx_plus_2NGHOSTS constants
      2) Open the `LOOP_OMP("omp parallel for", i0,i1,i2)` triple‐loop
      3) Compute xxL, xCart, call ID_function
      4) Convert through ADM→Cart, Cart→BSSN, BSSN→rescaled_RFM
      5) Set idx3 = IDX3(i0,i1,i2)
      6) Write out every BSSN field into gridfuncs->y_n_gfs, and
         (if enable_T4munu) T4UU fields into auxevol_gfs.

    :param enable_T4munu: whether to include the stress‐energy (T4UU) fields
    :returns: a raw string containing the entire loop + assignments
    """
    header = r"""
  const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

  LOOP_OMP("omp parallel for",
           i0, 0, Nxx_plus_2NGHOSTS0,
           i1, 0, Nxx_plus_2NGHOSTS1,
           i2, 0, Nxx_plus_2NGHOSTS2) {
    // xxL are the local coordinates on the destination grid
    const REAL xxL[3] = { xx[0][i0], xx[1][i1], xx[2][i2] };

    // xCart is the global Cartesian coordinate, which accounts for any grid offsets from the origin.
    REAL xCart[3];
    REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
    xx_to_Cart(params, xOrig, xCart);

    // Read or compute initial data at destination point xCart
    initial_data_struct initial_data;
    ID_function(commondata, params, xCart, ID_persist, &initial_data);

    ADM_Cart_basis_struct ADM_Cart_basis;
    ADM_SphorCart_to_Cart(commondata, params, xCart, &initial_data, &ADM_Cart_basis);

    BSSN_Cart_basis_struct BSSN_Cart_basis;
    ADM_Cart_to_BSSN_Cart(commondata, params, xCart, &ADM_Cart_basis, &BSSN_Cart_basis);

    rescaled_BSSN_rfm_basis_struct rescaled_BSSN_rfm_basis;
    BSSN_Cart_to_rescaled_BSSN_rfm(commondata, params,
                                   xxL,
                                   &BSSN_Cart_basis,
                                   &rescaled_BSSN_rfm_basis);

    const int idx3 = IDX3(i0, i1, i2);
"""

    # 2) build the list of fields to write out
    gf_list: List[str] = ["alpha", "trK", "cf"]
    for i in range(3):
        gf_list += [f"vetU{i}", f"betU{i}"]
        for j in range(i, 3):
            gf_list += [f"hDD{i}{j}", f"aDD{i}{j}"]
    if enable_T4munu:
        for mu in range(4):
            for nu in range(mu, 4):
                gf_list.append(f"T4UU{mu}{nu}")

    # 3) emit each assignment line
    lines: List[str] = []
    for field in sorted(gf_list):
        target = "auxevol_gfs" if field.startswith("T4UU") else "y_n_gfs"
        lines.append(
            f"    gridfuncs->{target}[IDX4pt({field.upper()}GF, idx3)]"
            f" = rescaled_BSSN_rfm_basis.{field};"
        )

    # 4) combine and return
    return header + "\n".join(lines) + "\n"


def build_lambdaU_zeroing_block() -> str:
    """
    Build the C code snippet that zeros out lambdaU and closes the OpenMP loop over all gridpoints.

    :return: The formatted C code block for lambda^i initialization
             and loop termination.
    """
    return r"""
    // Initialize lambdaU to zero
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU0GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU1GF, idx3)] = 0.0;
    gridfuncs->y_n_gfs[IDX4pt(LAMBDAU2GF, idx3)] = 0.0;
  } // END LOOP over all gridpoints on given grid
"""


def build_apply_inner_bcs_block(parallelization: Optional[str] = None) -> str:
    """
    Build the C code block that explains why inner boundary conditions must be applied and then invokes apply_bcs_inner_only on the gridfunctions.

    :param parallelization: The current parallelization mode (e.g., "cuda"), or None.
    :return: The formatted C code block for inner‐BC comments and call,
             with pointer adjustments for CUDA if needed.
    """
    code_block = r"""
  // Now we've set all but lambda^i, which will be computed via a finite-difference of hDD.
  //    However, hDD is not correctly set in inner boundary points so we apply inner bcs first.

  // Apply inner bcs to get correct values of all tensor quantities across symmetry boundaries;
  //    BSSN_Cart_to_rescaled_BSSN_rfm() converts each xCart->xx, which guarantees a mapping
  //    to the grid interior. It therefore does not account for parity conditions across
  //    symmetry boundaries being correct.
  apply_bcs_inner_only(commondata, params, bcstruct, gridfuncs->y_n_gfs);
"""
    if parallelization == "cuda":
        code_block = code_block.replace("gridfuncs->", "d_gridfuncs->").replace(
            " xx", " d_xx"
        )
    return code_block


def register_CFunction_initial_data_reader__convert_ADM_Sph_or_Cart_to_BSSN(
    CoordSystem: str,
    addl_includes: Optional[List[str]] = None,
    IDCoordSystem: str = "Spherical",
    enable_T4munu: bool = False,
    enable_fd_functions: bool = False,
    ID_persist_struct_str: str = "",
) -> None:
    """
    Register the CFunction for converting initial ADM data to BSSN variables.

    :param CoordSystem: Coordinate system for output BSSN variables.
    :param addl_includes: Additional header files to include.
    :param IDCoordSystem: Coordinate system for input ADM variables. Defaults to "Spherical".
    :param enable_T4munu: Whether to include stress-energy tensor components.
    :param enable_fd_functions: Whether to enable finite-difference functions.
    :param ID_persist_struct_str: String for persistent ID structure.
    """
    parallelization = par.parval_from_str("parallelization")

    includes, prefunc, lambdaU_launch = setup_ADM_initial_data_reader(
        ID_persist_struct_str=ID_persist_struct_str,
        enable_T4munu=enable_T4munu,
        enable_fd_functions=enable_fd_functions,
        addl_includes=addl_includes,
        CoordSystem=CoordSystem,
        IDCoordSystem=IDCoordSystem,
    )

    desc = f"Read ADM data in the {IDCoordSystem} basis, and output rescaled BSSN data in the {CoordSystem} basis"
    cfunc_type = "void"
    name = f"initial_data_reader__convert_ADM_{IDCoordSystem}_to_BSSN"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL *restrict xx[3], bc_struct *restrict bcstruct, MoL_gridfunctions_struct *restrict gridfuncs,
    ID_persist_struct *restrict ID_persist,
    void ID_function(const commondata_struct *restrict commondata, const params_struct *restrict params, const REAL xCart[3],
                     const ID_persist_struct *restrict ID_persist,
                     initial_data_struct *restrict initial_data)""".replace(
        "gridfuncs,",
        (
            "gridfuncs, MoL_gridfunctions_struct *restrict d_gridfuncs,"
            if parallelization in ["cuda"]
            else "gridfuncs,"
        ),
    ).replace(
        "const REAL *restrict xx[3],",
        (
            "const REAL *restrict xx[3], const REAL *restrict d_xx[3],"
            if parallelization in ["cuda"]
            else "const REAL *restrict xx[3],"
        ),
    )

    body = build_initial_data_conversion_loop(enable_T4munu)

    post_initial_data_call: str = ""
    if parallelization in ["cuda"]:
        post_initial_data_call = """
    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        cpyHosttoDevice__gf(commondata, params, gridfuncs->y_n_gfs, d_gridfuncs->y_n_gfs, which_gf, which_gf, params->grid_idx % NUM_STREAMS);
    }"""
        if enable_T4munu:
            post_initial_data_call += """
    for (int which_gf = 0; which_gf < NUM_AUXEVOL_GFS; which_gf++) {
        cpyHosttoDevice__gf(commondata, params, gridfuncs->auxevol_gfs, d_gridfuncs->auxevol_gfs, which_gf, which_gf, params->grid_idx % NUM_STREAMS);
    }
    """
        post_initial_data_call += "BHAH_DEVICE_SYNC();\n"

    body += build_lambdaU_zeroing_block()

    body += post_initial_data_call

    body += build_apply_inner_bcs_block(parallelization)

    if parallelization in ["cuda"]:
        lambdaU_launch = lambdaU_launch.replace("gridfuncs->", "d_gridfuncs->").replace(
            " xx", " d_xx"
        )
    body += lambdaU_launch

    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

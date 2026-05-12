"""
C function to save data to GRHayL structs from NRPy gridfunctions, and do basis transforms when needed.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.basis_transforms.jacobians as bt
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM


def register_CFunction_basis_transform_rfm_basis_to_Cartesian__read_cons_only(
    CoordSystem: str,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function for converting conservatives and metric quantities to Cartesian coordinates, and loading the data to structs to be used by GRHayL.

    This function performs the necessary coordinate transformations (basis changes)
    to convert data from the simulation's native coordinate system (e.g., Spherical)
    into the Cartesian basis expected by GRHayL.

    :param CoordSystem: The coordinate system.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: whether we're using grhayl to evolve temperature or not
    :param evolving_entropy: whether we're using grhayl to evolve entropy or not

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    # if enable_intrinsics:
    #     includes += [str(Path("simd") / "simd_intrinsics.h")]
    cfunc_type = "void"
    desc = f"""Convert conservative GRHD variables and ADM metric quantities from {CoordSystem} coordinates to Cartesian coordinates for GRHayL.

        @param[in] commondata Common simulation data.
        @param[in] params Grid-local runtime parameters.
        @param[out] cons Conservative variables in Cartesian basis.
        @param[out] metric ADM metric quantities in Cartesian basis.
        @param[in] i0 x0 grid index.
        @param[in] i1 x1 grid index.
        @param[in] i2 x2 grid index.
        @param[in] xx Reference-metric coordinate arrays.
        @param[in] auxevol_gfs Auxiliary/evolved primitive gridfunctions.
        @param[in] in_gfs Conservative gridfunctions.
        """

    name = "basis_transform_rfm_basis_to_Cartesian__read_cons_only"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, "
    params += "ghl_conservative_quantities *restrict cons, ghl_metric_quantities *restrict metric, "
    params += "const int i0, const int i1, const int i2, REAL *restrict xx[3], "
    params += "const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs"

    # We enforce precompute to False here
    enable_rfm_precompute = False

    # Step 3: Initialize Reference Metric and BSSN/ADM conversions
    # ADM in terms of BSSN
    AitoB = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    # Import all basic (unrescaled) BSSN scalars & tensors
    Bq = BSSN_quantities[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]

    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Step 4: Define Symbolic Variables for Scalar Quantities
    rho_star, tau_tilde = sp.symbols("rho_star tau_tilde", real=True)

    if evolving_temperature:
        Ye_star = sp.symbols("Ye_star", real=True)

    if evolving_entropy:
        S_star = sp.symbols("S_star", real=True)

    # Step 5: Process Vector Quantities
    # Declare rank-1 variables for momentum
    rescaledstildeD = ixp.declarerank1("rescaledstildeD")

    StildeD = ixp.zerorank1()
    for i in range(3):
        # Rescale: Unrescale the grid functions to get physical components in the local basis
        # S_i = S_i_{rescaled} / ReU[i] (Covariant vector rescaling)
        StildeD[i] = rescaledstildeD[i] / rfm.ReU[i]

    # Step 6: Perform Basis Transformations to Cartesian
    # Use Jacobians to transform Tensors and Vectors from the simulation basis (e.g., Spherical) to Cartesian.
    basis_transforms = bt.basis_transforms[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Transform Metric: gamma_{ij}
    gammaDD_cart = basis_transforms.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        AitoB.gammaDD,
    )

    # Transform Shift: beta^i
    betaU_Cart = basis_transforms.basis_transform_vectorU_from_rfmbasis_to_Cartesian(
        AitoB.betaU,
    )

    # Transform Momentum: S_i
    StildeD_Cart = basis_transforms.basis_transform_vectorD_from_rfmbasis_to_Cartesian(
        StildeD,
    )

    # Step 7: Prepare Code Generation
    # Read local coordinates to evaluate metric terms
    read_rfm_xx_arrays = r"""
    const REAL xx0 = xx[0][i0];
    const REAL xx1 = xx[1][i1];
    const REAL xx2 = xx[2][i2];
    """

    pre_body = read_rfm_xx_arrays

    # Define Left-Hand Side (Struct fields)
    lhs = [
        "cons->rho",
        "cons->tau",
        "cons->SD[0]",
        "cons->SD[1]",
        "cons->SD[2]",
        "metric->lapse",
        "metric->betaU[0]",
        "metric->betaU[1]",
        "metric->betaU[2]",
        "metric->gammaDD[0][0]",
        "metric->gammaDD[0][1]",
        "metric->gammaDD[0][2]",
        "metric->gammaDD[1][1]",
        "metric->gammaDD[1][2]",
        "metric->gammaDD[2][2]",
    ]

    if evolving_temperature:
        lhs.append("cons->Y_e")

    if evolving_entropy:
        lhs.append("cons->entropy")

    # Define Right-Hand Side (Symbolic Expressions in Cartesian basis)
    rhs = [
        rho_star,
        tau_tilde,
        StildeD_Cart[0],
        StildeD_Cart[1],
        StildeD_Cart[2],
        Bq.alpha,
        betaU_Cart[0],
        betaU_Cart[1],
        betaU_Cart[2],
        gammaDD_cart[0][0],
        gammaDD_cart[0][1],
        gammaDD_cart[0][2],
        gammaDD_cart[1][1],
        gammaDD_cart[1][2],
        gammaDD_cart[2][2],
    ]

    if evolving_temperature:
        rhs.append(Ye_star)

    if evolving_entropy:
        rhs.append(S_star)

    # Step 8: Generate C Code for Assignments
    expr_body = ccg.c_codegen(
        rhs,
        lhs,
        enable_simd=False,
        enable_fd_codegen=True,
        enable_GoldenKernels=enable_GoldenKernels,
    )

    # Step 9: Finalize Struct Initialization
    # Add a C call to ghl_initialize_metric so derived metric fields remain consistent
    post_body = r"""

       ghl_initialize_metric(
        metric->lapse,
        metric->betaU[0],
        metric->betaU[1],
        metric->betaU[2],
        metric->gammaDD[0][0],
        metric->gammaDD[0][1],
        metric->gammaDD[0][2],
        metric->gammaDD[1][1],
        metric->gammaDD[1][2],
        metric->gammaDD[2][2],
        metric);
"""

    # Step 10: Register the final C function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=pre_body + expr_body + post_body,
        enable_simd=False,
    )
    return pcg.NRPyEnv()

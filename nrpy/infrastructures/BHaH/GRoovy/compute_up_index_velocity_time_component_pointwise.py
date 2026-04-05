"""
C function to compute the up-index time component of the four velocity, U4Ut.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions as noif
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM


def register_CFunction_compute_up_index_velocity_time_component_pointwise(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_GoldenKernels: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to compute the time-component of the four-velocity, in a point-wise fashion.

    This function calculates u^0 based on the 3-velocity and metric. It also enforces
    a physical speed limit (gamma cap) by rescaling the velocity if the Lorentz factor
    exceeds a maximum value.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers and signature
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]

    desc = r"""
compute time component of four velocity, via

// Derivation of equation:
// \gamma_{ij} (v^i + \beta^i)(v^j + \beta^j)/(\alpha)^2
//   = \gamma_{ij} 1/(u^0)^2 ( \gamma^{ik} u_k \gamma^{jl} u_l /(\alpha)^2 <- Using Eq. 53 of arXiv:astro-ph/0503420
//   = 1/(u^0 \alpha)^2 u_j u_l \gamma^{jl}  <- Since \gamma_{ij} \gamma^{ik} = \delta^k_j
//   = 1/(u^0 \alpha)^2 ( (u^0 \alpha)^2 - 1 ) <- Using Eq. 56 of arXiv:astro-ph/0503420
//   = 1 - 1/(u^0 \alpha)^2 <= 1
"""
    cfunc_type = "void"
    name = "compute_up_index_velocity_time_component_pointwise"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const ghl_parameters *restrict ghl_params, "
    params += "const double alpha, const double vetU0, const double vetU1, const double vetU2, "
    params += "const double hDD00, const double hDD01, const double hDD02, "
    params += (
        "const double hDD11, const double hDD12, const double hDD22, const double cf, "
    )
    params += "double *restrict rescaledvU0, double *restrict rescaledvU1, double *restrict rescaledvU2, double *restrict u4Ut"

    # Step 3: Define C code for setup (before symbolic expressions)
    #         Calculates the speed limit threshold based on W_max (max Lorentz factor)
    body = r"""
    // 1 - W_max^{-2}
    const double one_minus_one_over_W_max_squared = 1.0 - ghl_params->inv_sq_max_Lorentz_factor;    
    
    const double rescaledvU_old0 = *rescaledvU0;
    const double rescaledvU_old1 = *rescaledvU1;
    const double rescaledvU_old2 = *rescaledvU2;
"""

    # Step 4: Initialize BSSN and ADM quantities
    # ADM in terms of BSSN
    AitoB = BSSN_to_ADM(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    # Import basic BSSN scalars & tensors
    Bq = BSSN_quantities[
        CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    ]

    rfm = refmetric.reference_metric[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]

    # Step 5: physical velocity vectors
    rescaledvU = ixp.declarerank1("rescaledvU_old")

    VU = ixp.zerorank1()
    utildeU = ixp.zerorank1()
    for i in range(3):
        # Compute physical v^i from rescaled input
        VU[i] = rescaledvU[i] * rfm.ReU[i]
        # utildeU^i = v^i + beta^i. This vector helps compute the Lorentz factor.
        utildeU[i] = VU[i] + Bq.betaU[i]

    # Step 6: Compute the squared norm used for Lorentz factor calculation
    # We calculate: (1 - 1/(alpha*u^0)^2) = gamma_{ij} (v^i + beta^i)(v^j + beta^j) / alpha^2
    # NOTE: We type hint as sp.Expr to prevent mypy from narrowing the type to Float,
    # which causes errors when assigning the clamped result later.
    one_minus_one_over_alpha_u0_squared: sp.Expr = sp.sympify(0)
    for i in range(3):
        for j in range(3):
            one_minus_one_over_alpha_u0_squared += (
                AitoB.gammaDD[i][j] * (VU[i] + Bq.betaU[i]) * (VU[j] + Bq.betaU[j])
            )

    one_minus_one_over_alpha_u0_squared /= Bq.alpha**2

    one_minus_one_over_alpha_u0_squared = sp.simplify(
        one_minus_one_over_alpha_u0_squared
    )

    # Step 7: Enforce Speed Limit
    # We check if the calculated squared norm exceeds the maximum allowed by W_max.
    # If it does, we must rescale the velocities.

    one_minus_one_over_W_max_squared = sp.symbols("one_minus_one_over_W_max_squared")

    # prevent divide by zero in cse
    TINYDOUBLE = sp.symbols("TINYDOUBLE")
    correction_fac = sp.sqrt(
        one_minus_one_over_W_max_squared
        / (one_minus_one_over_alpha_u0_squared + TINYDOUBLE)
    )

    # Use branchless logic (Min_Max_and_Piecewise_Expressions) to determine if we are violating bounds.
    # coord_geq_bound(x,xstar) returns 1.0 if x >= xstar, 0.0 otherwise.
    check_too_fast = noif.coord_geq_bound(
        one_minus_one_over_alpha_u0_squared, one_minus_one_over_W_max_squared
    )
    # coord_less_bound(x,xstar) returns 1.0 if x < xstar, 0.0 otherwise.
    check_normal = noif.coord_less_bound(
        one_minus_one_over_alpha_u0_squared, one_minus_one_over_W_max_squared
    )

    # Construct the new velocity vector using the boolean-like masks:
    # If too fast: Apply correction factor.
    # If normal: Keep original value.
    rescaledvU_new = ixp.zerorank1()
    for i in range(3):
        rescaledvU_new[i] = (
            check_too_fast
            * (
                sp.together(utildeU[i] / rfm.ReU[i]) * correction_fac
                - sp.together(Bq.betaU[i] / rfm.ReU[i])
            )
            + check_normal * rescaledvU[i]
        )

    # Clamp the squared norm to the maximum allowed value
    one_minus_one_over_alpha_u0_squared = noif.min_noif(
        one_minus_one_over_W_max_squared, one_minus_one_over_alpha_u0_squared
    )

    # Step 8: Compute final u^0
    # alpha * u^0 = 1 / sqrt(1 - (1 - 1/(alpha*u^0)^2))
    alpha_u0 = 1.0 / sp.sqrt(1.0 - one_minus_one_over_alpha_u0_squared)

    # u^0 = (alpha * u^0) / alpha
    u4Ut = alpha_u0 / Bq.alpha

    # Step 9: Generate C code
    # This outputs the calculated u^0 and the potentially rescaled velocities back to pointers.
    body += ccg.c_codegen(
        [u4Ut, rescaledvU_new[0], rescaledvU_new[1], rescaledvU_new[2]],
        ["*u4Ut", "*rescaledvU0", "*rescaledvU1", "*rescaledvU2"],
        verbose=False,
        include_braces=False,
        enable_fd_codegen=False,
        enable_simd=False,
        enable_GoldenKernels=enable_GoldenKernels,
    )

    # Step 10: Register the final C Function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
    )
    return pcg.NRPyEnv()

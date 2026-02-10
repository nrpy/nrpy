"""
C function to compute the stress energy tensor, T4UU, which is used to evolve the spacetime via the BSSN equations.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def register_CFunction_compute_stress_energy_tensor(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to compute components of the stress-energy tensor T^{mu nu}.

    This function handles the symbolic calculation of T^{mu nu} from GRHD primitives
    and BSSN metric quantities, then generates the corresponding C code using
    NRPy's infrastructure.

    :param CoordSystem: The coordinate system (e.g., "Cartesian", "Spherical").
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable SIMD intrinsics for performance.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels (optimized kernels).
    :param evolving_temperature: Whether we're using GRHayL to evolve temperature or not.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers and includes
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_intrinsics:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = "Compute the stress-energy tensor T4UU"
    cfunc_type = "void"
    name = "compute_stress_energy_tensor"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], "
    params += "const ghl_eos_parameters *restrict eos, "
    params += "const REAL *restrict in_gfs, REAL *restrict auxevol_gfs"

    # Adjust parameters if using reference metric precomputation
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )

    # To compute the densitized conserved momentum, we need the enthalpy h,
    # so we use GRHayL's ghl_compute_h_and_cs2 function, filling in prims struct to do so.
    pre_body_str = ""

    if evolving_temperature:
        pre_body_str += r"""
ghl_primitive_quantities prims;
ghl_initialize_primitives(auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)], 
                          auxevol_gfs[IDX4(PGF, i0, i1, i2)], 
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          auxevol_gfs[IDX4(YEGF, i0, i1, i2)], 
                          auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)], 
                          &prims);

// We must now compute eps and T
ghl_tabulated_enforce_bounds_rho_Ye_P(eos, &prims.rho, &prims.Y_e, &prims.press);
ghl_tabulated_compute_eps_T_from_P(eos, prims.rho, prims.Y_e, prims.press,
                                 &prims.eps, &prims.temperature);

auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)] = prims.temperature;

double h, cs2;
ghl_compute_h_and_cs2(eos, &prims, &h, &cs2);

"""

    else:
        pre_body_str += r"""
ghl_primitive_quantities prims;
ghl_initialize_primitives(auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)], 
                          auxevol_gfs[IDX4(PGF, i0, i1, i2)], 
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          NAN,
                          &prims);
double h, cs2;
ghl_compute_h_and_cs2(eos, &prims, &h, &cs2);

"""

    # Step 4: Initialize GRHD equations and BSSN quantities
    grhd_eqs = GRHD_Equations(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    rescaledvU = ixp.declarerank1("rescaledvU")

    # Access ADM variables (gamma_{ij}, alpha, beta^i) derived from BSSN variables
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

    # Step 5: Compute 4-velocity u^mu
    # Step 5a: Convert rescaled velocity (v^i / Re[i]) to physical 3-velocity v^i
    VU = ixp.zerorank1()
    for i in range(3):
        VU[i] = rescaledvU[i] * rfm.ReU[i]

    # Step 5b: Compute the Lorentz factor term W = alpha * u^0.
    # We use the normalization: -1 = g_{mu nu} u^mu u^nu
    # => -1 = -(alpha u^0)^2 + gamma_{ij} (u^i + beta^i u^0) (u^j + beta^j u^0)
    one_minus_one_over_alpha_u0_squared = sp.sympify(0.0)
    for i in range(3):
        for j in range(3):
            one_minus_one_over_alpha_u0_squared += (
                AitoB.gammaDD[i][j] * (VU[i] + Bq.betaU[i]) * (VU[j] + Bq.betaU[j])
            )

    one_minus_one_over_alpha_u0_squared /= Bq.alpha**2

    # W = alpha * u^0
    alpha_u0 = 1.0 / sp.sqrt(1.0 - one_minus_one_over_alpha_u0_squared)
    # u^0 = W / alpha
    u0 = alpha_u0 / Bq.alpha

    # Step 5c: Construct the full 4-velocity u^mu = (u^0, u^i)
    u4U = ixp.zerorank1(dimension=4)
    u4U[0] = u0
    for i in range(3):
        u4U[i + 1] = VU[i] * u0

    # Step 6: Pass computed velocities to GRHD class and compute T^{mu nu}
    grhd_eqs.VU = VU.copy()
    grhd_eqs.u4U = u4U.copy()

    grhd_eqs.compute_T4UU()

    # Step 7: Select the 10 independent components of the symmetric Stress-Energy Tensor
    output_vars = [
        grhd_eqs.T4UU[0][0],
        grhd_eqs.T4UU[0][1],
        grhd_eqs.T4UU[0][2],
        grhd_eqs.T4UU[0][3],
        grhd_eqs.T4UU[1][1],
        grhd_eqs.T4UU[1][2],
        grhd_eqs.T4UU[1][3],
        grhd_eqs.T4UU[2][2],
        grhd_eqs.T4UU[2][3],
        grhd_eqs.T4UU[3][3],
    ]

    # Step 8: Map symbolic outputs to specific grid functions in auxevol_gfs
    vars_grid_access = [
        gri.BHaHGridFunction.access_gf("T4UU00", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU01", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU02", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU03", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU11", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU12", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU13", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU22", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU23", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("T4UU33", gf_array_name="auxevol_gfs"),
    ]

    # Step 9: Generate C code (with SIMD/Golden Kernel support) and wrap in a loop
    body = lp.simple_loop(
        loop_body=pre_body_str
        + ccg.c_codegen(
            output_vars,
            vars_grid_access,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            enable_GoldenKernels=enable_GoldenKernels,
        ),
        loop_region="all points",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    # Step 10: Register the final C function
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=enable_intrinsics,
    )
    return pcg.NRPyEnv()

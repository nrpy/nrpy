"""
C function to compute source and connection-term contributions to the GRHD RHSs.

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
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def register_CFunction_calculate_all_source_terms(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register C function to compute the source-term contributions to the GRHD right-hand sides.

    This routine adds the geometric source and connection terms for the
    conserved GRHD variables. Optional electron-fraction and entropy channels
    are included when those variables are evolved.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric
        precomputation.
    :param enable_intrinsics: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: Whether temperature and electron fraction are
        evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, signature, and parameters.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_intrinsics:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = (
        "Add source and connection terms to rho_star, tau_tilde, "
        "rescaledStildeD, and optional composition variables"
    )
    cfunc_type = "void"
    name = "calculate_all_source_terms"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "REAL *restrict xx[3], "
        "const ghl_eos_parameters *restrict eos, "
        "REAL *restrict auxevol_gfs, "
        "const REAL *restrict in_gfs, "
        "REAL *restrict rhs_gfs"
    )
    if enable_rfm_precompute:
        params = params.replace(
            "REAL *restrict xx[3]", "const rfm_struct *restrict rfmstruct"
        )

    # Step 3: Initialize the primitive state needed to compute h and c_s^2.
    entropy_arg = "auxevol_gfs[IDX4(SGF, i0, i1, i2)]" if evolving_entropy else "NAN"
    ye_arg = "auxevol_gfs[IDX4(YEGF, i0, i1, i2)]" if evolving_temperature else "NAN"
    temperature_arg = (
        "auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]"
        if evolving_temperature
        else "NAN"
    )
    pre_body = f"""
ghl_primitive_quantities prims;
ghl_initialize_primitives(
    auxevol_gfs[IDX4(RHOBGF, i0, i1, i2)],
    auxevol_gfs[IDX4(PGF, i0, i1, i2)],
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    {entropy_arg},
    {ye_arg},
    {temperature_arg},
    &prims);

double h, cs2;
ghl_compute_h_and_cs2(eos, &prims, &h, &cs2);
"""

    # Step 4: Build the symbolic source terms.
    grhd_eqs = GRHD_Equations(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    grhd_eqs.construct_all_equations()

    rho_star_rhs = -grhd_eqs.rho_star_connection_term
    tau_tilde_rhs = grhd_eqs.tau_source_term - grhd_eqs.tau_connection_term

    StildeD_rhs = ixp.zerorank1()
    rescaledStildeD_rhs = ixp.zerorank1()
    for i in range(3):
        StildeD_rhs[i] = (
            grhd_eqs.S_tilde_source_termD[i] - grhd_eqs.S_tilde_connection_termsD[i]
        )
        rescaledStildeD_rhs[i] = sp.together(StildeD_rhs[i] * grhd_eqs.ReU[i])

    output_vars = [
        rho_star_rhs,
        tau_tilde_rhs,
        rescaledStildeD_rhs[0],
        rescaledStildeD_rhs[1],
        rescaledStildeD_rhs[2],
    ]
    vars_grid_access = [
        gri.BHaHGridFunction.access_gf("rho_star", gf_array_name="rhs_gfs"),
        gri.BHaHGridFunction.access_gf("tau_tilde", gf_array_name="rhs_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD0", gf_array_name="rhs_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD1", gf_array_name="rhs_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD2", gf_array_name="rhs_gfs"),
    ]

    if evolving_temperature:
        output_vars.append(-grhd_eqs.Ye_star_connection_term)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf("Ye_star", gf_array_name="rhs_gfs")
        )

    if evolving_entropy:
        output_vars.append(-grhd_eqs.S_star_connection_term)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf("S_star", gf_array_name="rhs_gfs")
        )

    # Step 5: Wrap the pointwise update in an interior loop.
    body = lp.simple_loop(
        loop_body=pre_body
        + ccg.c_codegen(
            output_vars,
            vars_grid_access,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            enable_GoldenKernels=enable_GoldenKernels,
        ),
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    # Step 6: Register the final C function.
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

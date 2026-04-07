"""
C function to add divergence of cell-centered HLL GRHD fluxes to RHSs.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def register_CFunction_calculate_flux_divergences(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int,
    enable_intrinsics: bool = False,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the flux-divergence contribution to the GRHD right-hand sides.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric
        precomputation.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_intrinsics: Whether to enable SIMD.
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
        "Compute flux divergences for rho_star, tau_tilde, rescaledStildeD, "
        "and optional composition variables."
    )
    cfunc_type = "void"
    name = "calculate_flux_divergences"
    params = (
        "const int flux_dirn, "
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

    # Step 3: Build the symbolic reference-metric rescaling factors.
    rfm_key = CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    rfm = refmetric.reference_metric[rfm_key]

    ReUD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            ReUD[i][j] = rfm.ReU[i] / rfm.ReU[j]

    ReUDdD = ixp.zerorank3()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                ReUDdD[i][j][k] = sp.diff(ReUD[i][j], rfm.xx[k])

    # Step 4: Build the symbolic rescaled fluxes needed for the divergence.
    grhd_eqs = GRHD_Equations(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )
    grhd_eqs.compute_T4UU()
    grhd_eqs.compute_T4UD()
    grhd_eqs.compute_rho_star()
    grhd_eqs.compute_rho_star_fluxU()

    if evolving_temperature:
        grhd_eqs.compute_Ye_star()
        grhd_eqs.compute_Ye_star_fluxU()

    if evolving_entropy:
        grhd_eqs.compute_S_star()
        grhd_eqs.compute_S_star_fluxU()

    grhd_eqs.compute_tau_tilde()
    grhd_eqs.compute_tau_tilde_fluxU()
    grhd_eqs.compute_S_tildeD()
    grhd_eqs.compute_S_tilde_fluxUD()

    c_codegen_lhs: List[str] = []
    c_codegen_rhs: List[sp.Expr] = []
    for i in range(3):
        c_codegen_lhs.append(f"ReU[{i}]")
        c_codegen_rhs.append(rfm.ReU[i])
        for j in range(3):
            c_codegen_lhs.append(f"ReUdD[{i}][{j}]")
            c_codegen_rhs.append(rfm.ReUdD[i][j])
            c_codegen_lhs.append(f"ReUD[{i}][{j}]")
            c_codegen_rhs.append(ReUD[i][j])
            for k in range(3):
                c_codegen_lhs.append(f"ReUDdD[{i}][{j}][{k}]")
                c_codegen_rhs.append(ReUDdD[i][j][k])

    for i in range(3):
        c_codegen_lhs.append(f"rescaledtau_tildefluxU[{i}]")
        c_codegen_rhs.append(grhd_eqs.rescaled_tau_tilde_fluxU[i])
        c_codegen_lhs.append(f"rescaledrho_starfluxU[{i}]")
        c_codegen_rhs.append(grhd_eqs.rescaled_rho_star_fluxU[i])
        for j in range(3):
            c_codegen_lhs.append(f"rescaledStildefluxUD[{i}][{j}]")
            c_codegen_rhs.append(grhd_eqs.rescaled_S_tilde_fluxUD[i][j])

    if evolving_temperature:
        for i in range(3):
            c_codegen_lhs.append(f"rescaledYe_starfluxU[{i}]")
            c_codegen_rhs.append(grhd_eqs.rescaled_Ye_star_fluxU[i])

    if evolving_entropy:
        for i in range(3):
            c_codegen_lhs.append(f"rescaledS_starfluxU[{i}]")
            c_codegen_rhs.append(grhd_eqs.rescaled_S_star_fluxU[i])

    pre_loop = r"""
    const REAL invdxi[3] = {invdxx0, invdxx1, invdxx2};
    const REAL invdx = invdxi[flux_dirn];

    const int xdir = (flux_dirn == 0);
    const int ydir = (flux_dirn == 1);
    const int zdir = (flux_dirn == 2);
"""

    # Step 5: Build the pointwise C code that evaluates the symbolic fluxes.
    entropy_arg = "auxevol_gfs[IDX4(SGF, i0, i1, i2)]" if evolving_entropy else "NAN"
    ye_arg = "auxevol_gfs[IDX4(YEGF, i0, i1, i2)]" if evolving_temperature else "NAN"
    temperature_arg = (
        "auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)]"
        if evolving_temperature
        else "NAN"
    )

    pre_body = r"""
    const int index = IDX3(i0, i1, i2);
    const int indexp1 = IDX3(i0 + xdir, i1 + ydir, i2 + zdir);

    REAL ReU[3];
    REAL ReUdD[3][3];
    REAL ReUD[3][3];
    REAL ReUDdD[3][3][3];

    REAL rescaledStildefluxUD[3][3];
    REAL rescaledtau_tildefluxU[3];
    REAL rescaledrho_starfluxU[3];
"""

    if evolving_temperature:
        pre_body += r"""
    REAL rescaledYe_starfluxU[3];
"""

    if evolving_entropy:
        pre_body += r"""
    REAL rescaledS_starfluxU[3];
"""

    pre_body += f"""
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

    generated_flux_code = ccg.c_codegen(
        c_codegen_rhs,
        c_codegen_lhs,
        verbose=False,
        include_braces=False,
        enable_fd_codegen=True,
        enable_simd=enable_intrinsics,
        enable_GoldenKernels=enable_GoldenKernels,
    )

    full_body = pre_body + generated_flux_code + r"""
    rhs_gfs[IDX4pt(RHO_STARGF, index)] +=
        -rescaledrho_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn]
        + ReU[flux_dirn]
              * (auxevol_gfs[IDX4pt(RHO_STAR_HLL_FLUXGF, index)]
                 - auxevol_gfs[IDX4pt(RHO_STAR_HLL_FLUXGF, indexp1)])
              * invdx;

    rhs_gfs[IDX4pt(TAU_TILDEGF, index)] +=
        -rescaledtau_tildefluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn]
        + ReU[flux_dirn]
              * (auxevol_gfs[IDX4pt(TAU_TILDE_HLL_FLUXGF, index)]
                 - auxevol_gfs[IDX4pt(TAU_TILDE_HLL_FLUXGF, indexp1)])
              * invdx;

    REAL StildeD0_rhs;
    REAL StildeD1_rhs;
    REAL StildeD2_rhs;

    StildeD0_rhs =
        -rescaledStildefluxUD[flux_dirn][0] * ReUDdD[flux_dirn][0][flux_dirn]
        + ReUD[flux_dirn][0]
              * (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD0GF, index)]
                 - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD0GF, indexp1)])
              * invdx;
    StildeD1_rhs =
        -rescaledStildefluxUD[flux_dirn][1] * ReUDdD[flux_dirn][1][flux_dirn]
        + ReUD[flux_dirn][1]
              * (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD1GF, index)]
                 - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD1GF, indexp1)])
              * invdx;
    StildeD2_rhs =
        -rescaledStildefluxUD[flux_dirn][2] * ReUDdD[flux_dirn][2][flux_dirn]
        + ReUD[flux_dirn][2]
              * (auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD2GF, index)]
                 - auxevol_gfs[IDX4pt(RESCALEDSTILDE_FLUX_HLLD2GF, indexp1)])
              * invdx;

    // Rescale the momentum RHSs back to the evolved variables.
    rhs_gfs[IDX4pt(RESCALEDSTILDED0GF, index)] += StildeD0_rhs * ReU[0];
    rhs_gfs[IDX4pt(RESCALEDSTILDED1GF, index)] += StildeD1_rhs * ReU[1];
    rhs_gfs[IDX4pt(RESCALEDSTILDED2GF, index)] += StildeD2_rhs * ReU[2];
"""

    if evolving_temperature:
        full_body += r"""
    rhs_gfs[IDX4pt(YE_STARGF, index)] +=
        -rescaledYe_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn]
        + ReU[flux_dirn]
              * (auxevol_gfs[IDX4pt(YE_STAR_HLL_FLUXGF, index)]
                 - auxevol_gfs[IDX4pt(YE_STAR_HLL_FLUXGF, indexp1)])
              * invdx;
"""

    if evolving_entropy:
        full_body += r"""
    rhs_gfs[IDX4pt(S_STARGF, index)] +=
        -rescaledS_starfluxU[flux_dirn] * ReUdD[flux_dirn][flux_dirn]
        + ReU[flux_dirn]
              * (auxevol_gfs[IDX4pt(S_STAR_HLL_FLUXGF, index)]
                 - auxevol_gfs[IDX4pt(S_STAR_HLL_FLUXGF, indexp1)])
              * invdx;
"""

    # Step 6: Wrap the pointwise divergence update in an interior loop.
    loop_body = lp.simple_loop(
        loop_body=full_body,
        loop_region="interior",
        enable_intrinsics=enable_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )

    # Step 7: Register the final C function.
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=pre_loop + loop_body,
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

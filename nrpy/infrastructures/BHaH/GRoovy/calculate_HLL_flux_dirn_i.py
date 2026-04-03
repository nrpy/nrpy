"""
C function to generate HLL fluxes along a chosen coordinate direction.

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
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.equations.grhd.HLL_fluxes import calculate_HLL_fluxes


def register_CFunction_calculate_HLL_flux_dirn_i(
    flux_dirn: int,
    enable_intrinsics: bool,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the HLL flux routine for one flux direction.

    :param flux_dirn: Flux direction.
    :param enable_intrinsics: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: Whether temperature and electron fraction are
        evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings, clang_format
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp"]
    >>> name = "calculate_HLL_fluxes_direction_0"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    cfc.CFunction_dict.clear()
    ...    _ = register_CFunction_calculate_HLL_flux_dirn_i(0, False, 1)
    ...    generated_str = clang_format(cfc.CFunction_dict[name].full_function)
    ...    validation_desc = f"{name}__{parallelization}"
    ...    _ = validate_strings(generated_str, validation_desc, file_ext="c")
    Setting up reference_metric[Cartesian]...
    Setting up BSSN_Quantities[Cartesian]...
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up C function headers, description, and parameters.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_intrinsics:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = (
        "Compute HLL fluxes on the "
        + str(flux_dirn)
        + "-face for the GRHD evolution variables."
    )
    cfunc_type = "void"
    name = "calculate_HLL_fluxes_direction_" + str(flux_dirn)
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "const ghl_eos_parameters *restrict eos, "
        "REAL *restrict auxevol_gfs"
    )

    # Step 3: Initialize the right and left primitive states needed by the EOS.
    entropy_r_arg = (
        "auxevol_gfs[IDX4(S_RGF, i0, i1, i2)]" if evolving_entropy else "NAN"
    )
    entropy_l_arg = (
        "auxevol_gfs[IDX4(S_LGF, i0, i1, i2)]" if evolving_entropy else "NAN"
    )
    ye_r_arg = (
        "auxevol_gfs[IDX4(YE_RGF, i0, i1, i2)]" if evolving_temperature else "NAN"
    )
    ye_l_arg = (
        "auxevol_gfs[IDX4(YE_LGF, i0, i1, i2)]" if evolving_temperature else "NAN"
    )
    temperature_r_arg = (
        "auxevol_gfs[IDX4(TEMPERATURE_RGF, i0, i1, i2)]"
        if evolving_temperature
        else "NAN"
    )
    temperature_l_arg = (
        "auxevol_gfs[IDX4(TEMPERATURE_LGF, i0, i1, i2)]"
        if evolving_temperature
        else "NAN"
    )
    pre_body = f"""
// Initialize GRHayL primitive states on the right and left sides.
ghl_primitive_quantities prims_r, prims_l;
ghl_initialize_primitives(
    auxevol_gfs[IDX4(RHOB_RGF, i0, i1, i2)],
    auxevol_gfs[IDX4(P_RGF, i0, i1, i2)],
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    {entropy_r_arg},
    {ye_r_arg},
    {temperature_r_arg},
    &prims_r);

ghl_initialize_primitives(
    auxevol_gfs[IDX4(RHOB_LGF, i0, i1, i2)],
    auxevol_gfs[IDX4(P_LGF, i0, i1, i2)],
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    NAN,
    {entropy_l_arg},
    {ye_l_arg},
    {temperature_l_arg},
    &prims_l);

double h_r, h_l, cs2_r, cs2_l;
ghl_compute_h_and_cs2(eos, &prims_r, &h_r, &cs2_r);
ghl_compute_h_and_cs2(eos, &prims_l, &h_l, &cs2_l);
"""

    # Step 4: Build the symbolic face-centered metric and fluid states.
    rfm = refmetric.reference_metric["Cartesian"]

    alpha_face = sp.symbols("alpha_face", real=True)
    cf_face = sp.symbols("cf_face", real=True)
    h_faceDD = ixp.declarerank2("h_faceDD", symmetry="sym01", dimension=3)
    vet_faceU = ixp.declarerank1("vet_faceU", dimension=3)

    rescaledvrU = ixp.declarerank1("rescaledvrU", dimension=3)
    rescaledvlU = ixp.declarerank1("rescaledvlU", dimension=3)

    VrU = ixp.zerorank1()
    VlU = ixp.zerorank1()

    u4rUt = sp.symbols("u4rUt", real=True)
    u4lUt = sp.symbols("u4lUt", real=True)
    u4rU = ixp.zerorank1(dimension=4)
    u4rU[0] = u4rUt
    u4lU = ixp.zerorank1(dimension=4)
    u4lU[0] = u4lUt

    beta_faceU = ixp.zerorank1()
    for i in range(3):
        VrU[i] = rescaledvrU[i] * rfm.ReU[i]
        VlU[i] = rescaledvlU[i] * rfm.ReU[i]
        u4rU[i + 1] = VrU[i] * u4rU[0]
        u4lU[i + 1] = VlU[i] * u4lU[0]
        beta_faceU[i] = vet_faceU[i] * rfm.ReU[i]

    rho_b_r = sp.symbols("rhob_r", real=True)
    rho_b_l = sp.symbols("rhob_l", real=True)
    Ye_r = sp.symbols("Ye_r", real=True)
    Ye_l = sp.symbols("Ye_l", real=True)
    S_r = sp.symbols("S_r", real=True)
    S_l = sp.symbols("S_l", real=True)
    P_r = sp.symbols("P_r", real=True)
    P_l = sp.symbols("P_l", real=True)
    h_r = sp.symbols("h_r", real=True)
    h_l = sp.symbols("h_l", real=True)
    cs2_r = sp.symbols("cs2_r", real=True)
    cs2_l = sp.symbols("cs2_l", real=True)

    AitoB = BSSN_to_ADM(CoordSystem="Cartesian")
    Bq = BSSN_quantities["Cartesian"]
    e6phi_face = ((Bq.exp_m4phi ** sp.Rational(1, 2)) ** (-3)).subs(Bq.cf, cf_face)

    gamma_faceDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gamma_faceDD[i][j] = (
                AitoB.gammaDD[i][j]
                .subs(Bq.hDD[i][j], h_faceDD[i][j])
                .subs(Bq.cf, cf_face)
            )

    # Step 5: Build the symbolic HLL flux expressions.
    (
        rho_star_HLL_flux,
        Ye_star_HLL_flux,
        S_star_HLL_flux,
        tau_tilde_HLL_flux,
        rescaledStilde_flux_HLLD,
    ) = calculate_HLL_fluxes(
        flux_dirn,
        alpha_face,
        gamma_faceDD,
        beta_faceU,
        e6phi_face,
        u4rU,
        u4lU,
        rho_b_r,
        rho_b_l,
        Ye_r,
        Ye_l,
        S_r,
        S_l,
        P_r,
        P_l,
        h_r,
        h_l,
        cs2_r,
        cs2_l,
    )

    output_vars = [
        rho_star_HLL_flux,
        tau_tilde_HLL_flux,
        rescaledStilde_flux_HLLD[0],
        rescaledStilde_flux_HLLD[1],
        rescaledStilde_flux_HLLD[2],
    ]
    vars_grid_access = [
        gri.BHaHGridFunction.access_gf(
            "rho_star_HLL_flux", gf_array_name="auxevol_gfs"
        ),
        gri.BHaHGridFunction.access_gf(
            "tau_tilde_HLL_flux", gf_array_name="auxevol_gfs"
        ),
        gri.BHaHGridFunction.access_gf(
            "rescaledstilde_flux_HLLD0", gf_array_name="auxevol_gfs"
        ),
        gri.BHaHGridFunction.access_gf(
            "rescaledstilde_flux_HLLD1", gf_array_name="auxevol_gfs"
        ),
        gri.BHaHGridFunction.access_gf(
            "rescaledstilde_flux_HLLD2", gf_array_name="auxevol_gfs"
        ),
    ]

    if evolving_temperature:
        output_vars.append(Ye_star_HLL_flux)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf(
                "Ye_star_HLL_flux", gf_array_name="auxevol_gfs"
            )
        )

    if evolving_entropy:
        output_vars.append(S_star_HLL_flux)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf(
                "S_star_HLL_flux", gf_array_name="auxevol_gfs"
            )
        )

    # Step 6: Wrap the pointwise HLL solve in the BH@H loop.
    body = lp.simple_loop(
        loop_body=pre_body
        + ccg.c_codegen(
            output_vars,
            vars_grid_access,
            enable_fd_codegen=True,
            enable_simd=enable_intrinsics,
            enable_GoldenKernels=enable_GoldenKernels,
        ),
        loop_region="interior plus one upper",
        enable_intrinsics=enable_intrinsics,
        enable_rfm_precompute=False,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )

    # Step 7: Register the final C function.
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
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

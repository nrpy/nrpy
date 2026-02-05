"""
C function to compute conservative variables from primitive variables.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from pathlib import Path
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri  # NRPy+: Functionality for handling numerical grids
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.BHaH.simple_loop as lp
from nrpy.equations.grhd.GRHD_equations import GRHD_Equations


def register_CFunction_primitives_to_conservatives_routine(
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_intrinsics: bool,
    OMP_collapse: int,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register function to compute evolved, conserved variables from primitives.

    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_intrinsics: Whether to enable SIMD.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: If True, assume temperature and electron fraction are evolved primitives and also compute the corresponding densitized Ye_star.
    :param evolving_entropy: If True, assume entropy is an evolved primitive and also compute the corresponding densitized S_star.

    :return: An NRPyEnv_type object if registration is successful, otherwise None.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    if enable_intrinsics:
        includes += [str(Path("simd") / "simd_intrinsics.h")]
    desc = "Primitives to Conservatives Routine"
    cfunc_type = "void"
    name = "primitives_to_conservatives_routine"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], "
    params += "const ghl_eos_parameters *restrict eos, "
    params += "REAL *restrict auxevol_gfs, REAL *restrict in_gfs"

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
ghl_tabulated_compute_eps_T_from_P(eos, prims.rho, prims.Y_e, prims.press, &prims.eps, &prims.temperature);

auxevol_gfs[IDX4(TEMPERATUREGF, i0, i1, i2)] = prims.temperature;

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

"""
    pre_body_str += r"""
double h, cs2;
ghl_compute_h_and_cs2(eos, &prims, &h, &cs2);
"""

    grhd_eqs = GRHD_Equations(
        CoordSystem=CoordSystem, enable_rfm_precompute=enable_rfm_precompute
    )

    grhd_eqs.construct_all_equations()

    output_vars = [
        grhd_eqs.rho_star,
        grhd_eqs.tau_tilde,
        grhd_eqs.rescaledS_tildeD[0],
        grhd_eqs.rescaledS_tildeD[1],
        grhd_eqs.rescaledS_tildeD[2],
    ]

    vars_grid_access = [
        gri.BHaHGridFunction.access_gf("rho_star", gf_array_name="in_gfs"),
        gri.BHaHGridFunction.access_gf("tau_tilde", gf_array_name="in_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD0", gf_array_name="in_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD1", gf_array_name="in_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledstildeD2", gf_array_name="in_gfs"),
    ]

    if evolving_temperature:
        output_vars.append(grhd_eqs.Ye_star)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf("Ye_star", gf_array_name="in_gfs")
        )

    if evolving_entropy:
        output_vars.append(grhd_eqs.S_star)
        vars_grid_access.append(
            gri.BHaHGridFunction.access_gf("S_star", gf_array_name="in_gfs")
        )

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

"""
C function registration for Hamiltonian-constraint residual evaluation for hyperbolic relaxation.

This module constructs and registers the C helper routine "residual_H_compute_all_points".
The generated C function evaluates the Hamiltonian-constraint residual across the interior
of the grid for the selected coordinate system and writes the result into a destination
gridfunction array. The loop nest is specialized to the coordinate system and parallelized
with OpenMP as requested. At present, reference-metric precomputation is disabled at
code-generation time, so the generated C signature accepts xx coordinate arrays.

Function
--------
register_CFunction_residual_H_compute_all_points
    Construct and register the "residual_H_compute_all_points" C function.

Authors: Thiago Assumpcao
         assumpcaothiago **at** gmail **dot** com
         Zachariah B. Etienne
         zachetie **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.equations.nrpyelliptic.ConformallyFlat_RHSs import (
    HyperbolicRelaxationCurvilinearRHSs,
)
from nrpy.infrastructures import BHaH


# Define function to compute residual the solution
def register_CFunction_residual_H_compute_all_points(
    CoordSystem: str,
    enable_simd_intrinsics: bool,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a C function that computes the Hamiltonian-constraint residual on all interior points.

    This function generates and registers the C helper "residual_H_compute_all_points", which loops over the
    interior of the grid and computes the Hamiltonian-constraint residual used by the hyperbolic relaxation
    scheme. The loop structure is specialized to the chosen coordinate system and parallelized using OpenMP
    with the requested collapse level. Note that, in the current implementation, reference-metric precomputation
    is forced off at code-generation time, so the emitted C API always accepts xx coordinate arrays and does not
    take an rfm_struct argument. In the generated kernel, commondata is read-only simulation metadata and may be
    unused by this routine.

    :param CoordSystem: Name of the coordinate system that specializes the generated loop bounds and index macros.
    :param enable_simd_intrinsics: Whether to enable SIMD intrinsics.
    :param OMP_collapse: OpenMP collapse level for the nested loops.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> name = "residual_H_compute_all_points"
    >>> for CoordSystem in unittest_CoordSystems:
    ...     cfc.CFunction_dict.clear()
    ...     _ = register_CFunction_residual_H_compute_all_points(CoordSystem, True, True)
    ...     generated_str = cfc.CFunction_dict[f"{name}__rfm__{CoordSystem}"].full_function
    ...     validation_desc = f"{name}__{CoordSystem}"
    ...     validate_strings(generated_str, validation_desc)
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    Setting up reference_metric[HoleySinhSpherical_rfm_precompute]...
    Setting up reference_metric[Cartesian_rfm_precompute]...
    Setting up reference_metric[SinhCylindricalv2n2_rfm_precompute]...
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    includes = ["BHaH_defines.h"]
    if enable_simd_intrinsics:
        includes += ["intrinsics/simd_intrinsics.h"]
    desc = """
 * @file residual_H_compute_all_points.c
 * @brief Compute the Hamiltonian-constraint residual on all interior grid points and store the result.
 *
 * The generated function "residual_H_compute_all_points" iterates over the interior region of the grid
 * for the selected coordinate system and evaluates the Hamiltonian-constraint residual used by the
 * hyperbolic relaxation scheme. Results are written into the destination gridfunction buffer.
 *
 * Reference-metric precomputation is currently disabled at code-generation time, so this routine
 * accepts xx coordinate arrays and does not take an rfm_struct parameter.
 *
 * If a user-editable block is present in the implementation, users may insert custom logic such as
 * additional diagnostics or instrumentation without changing the function interface.
 *
 * @param[in]  commondata        Pointer to read-only global simulation metadata (e.g., time, step counters); may be unused.
 * @param[in]  params            Pointer to read-only per-grid parameters (sizes, ghost zones, strides, names).
 * @param[in]  xx                Array of three coordinate arrays used for coordinate-dependent operations.
 * @param[in]  auxevol_gfs       Pointer to read-only auxiliary evolution gridfunctions required by the residual.
 * @param[in]  in_gfs            Pointer to read-only input gridfunctions (e.g., current solution fields).
 * @param[out] dest_gf_address   Pointer to the destination gridfunction buffer where the residual is stored.
 *
 * @return     void.
"""
    cfunc_type = "void"
    name = "residual_H_compute_all_points"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                const rfm_struct *restrict rfmstruct, const REAL *restrict auxevol_gfs, const REAL *restrict in_gfs,
                REAL *restrict dest_gf_address"""

    # Populate residual_H
    rhs = HyperbolicRelaxationCurvilinearRHSs(CoordSystem, enable_rfm_precompute=True)
    orig_par = par.parval_from_str("parallelization")
    par.set_parval_from_str("parallelization", "openmp")
    body = BHaH.simple_loop.simple_loop(
        loop_body=ccg.c_codegen(
            [rhs.residual],
            ["dest_gf_address[IDX3(i0,i1,i2)]"],
            enable_fd_codegen=True,
            enable_simd=enable_simd_intrinsics,
        ),
        loop_region="interior",
        enable_intrinsics=enable_simd_intrinsics,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=True,
        read_xxs=False,
        OMP_collapse=OMP_collapse,
    )
    par.set_parval_from_str("parallelization", orig_par)

    cfc.register_CFunction(
        subdirectory="diagnostics",
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


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

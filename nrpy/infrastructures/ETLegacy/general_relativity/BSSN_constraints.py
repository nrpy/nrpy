"""
Generates function to compute the constraints H, MU, and MSQUARED.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.ETLegacy.simple_loop as lp
import nrpy.params as par
from nrpy.equations.general_relativity.BSSN_constraints import BSSN_constraints
from nrpy.helpers.expression_utils import get_params_commondata_symbols_from_expr_list
from nrpy.infrastructures.ETLegacy.CodeParameters import read_CodeParameters
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)


def register_CFunction_BSSN_constraints(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    enable_T4munu: bool,
    enable_simd: bool,
    fd_order: int,
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the BSSN constraints evaluation function.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param enable_T4munu: Whether to include the stress-energy tensor.
    :param enable_simd: Whether to enable SIMD instructions.
    :param fd_order: Order of finite difference method
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> vacuum_thorn = "DoctestBaikalVacuumConstraints"
    >>> matter_thorn = "DoctestBaikalConstraints"
    >>> cfunction_names = [
    ...     f"{vacuum_thorn}_BSSN_constraints_order_4",
    ...     f"{matter_thorn}_BSSN_constraints_order_4",
    ... ]
    >>> original_infrastructure = par.parval_from_str("Infrastructure")
    >>> original_parallel_codegen = par.parval_from_str("enable_parallel_codegen")
    >>> original_fd_order = par.parval_from_str("fd_order")
    >>> original_c_codegen = ccg.c_codegen
    >>> cfunctions = []
    >>> try:
    ...     par.set_parval_from_str("Infrastructure", "ETLegacy")
    ...     par.set_parval_from_str("enable_parallel_codegen", False)
    ...     ccg.c_codegen = lambda *_args, **_kwargs: ""
    ...     with contextlib.redirect_stdout(io.StringIO()):
    ...         for thorn_name, enable_t4munu in [
    ...             (vacuum_thorn, False),
    ...             (matter_thorn, True),
    ...         ]:
    ...             _ = register_CFunction_BSSN_constraints(
    ...                 thorn_name=thorn_name,
    ...                 CoordSystem="Cartesian",
    ...                 enable_rfm_precompute=False,
    ...                 enable_T4munu=enable_t4munu,
    ...                 enable_simd=True,
    ...                 fd_order=4,
    ...             )
    ...     cfunctions = [cfc.CFunction_dict[name] for name in cfunction_names]
    ... finally:
    ...     ccg.c_codegen = original_c_codegen
    ...     par.set_parval_from_str("Infrastructure", original_infrastructure)
    ...     par.set_parval_from_str(
    ...         "enable_parallel_codegen", original_parallel_codegen
    ...     )
    ...     par.set_parval_from_str("fd_order", original_fd_order)
    ...     for name in cfunction_names:
    ...         _ = cfc.CFunction_dict.pop(name, None)
    >>> vacuum_pi_lookup = f'CCTK_ParameterGet("PI", "{vacuum_thorn}", NULL)'
    >>> matter_pi_lookup = f'CCTK_ParameterGet("PI", "{matter_thorn}", NULL)'
    >>> [
    ...     vacuum_pi_lookup in cfunctions[0].body,
    ...     "PI" in (cfunctions[0].ET_current_thorn_CodeParams_used or []),
    ...     matter_pi_lookup in cfunctions[1].body,
    ...     "PI" in (cfunctions[1].ET_current_thorn_CodeParams_used or []),
    ... ]
    [False, False, True, True]
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    old_fd_order = par.parval_from_str("fd_order")
    # Set this because parallel codegen needs the correct local values
    par.set_parval_from_str("fd_order", fd_order)
    enable_RbarDD_gridfunctions = False

    includes = define_standard_includes()
    if enable_simd:
        includes += ["./simd/simd_intrinsics.h"]
    desc = r"""Evaluate BSSN constraints."""
    name = f"{thorn_name}_BSSN_constraints_order_{fd_order}"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
"""
    if enable_simd:
        body += r"""
  const REAL_SIMD_ARRAY invdxx0 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(0));
  const REAL_SIMD_ARRAY invdxx1 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(1));
  const REAL_SIMD_ARRAY invdxx2 CCTK_ATTRIBUTE_UNUSED = ConstSIMD(1.0/CCTK_DELTA_SPACE(2));
"""
    else:
        body += """  DECLARE_CCTK_PARAMETERS;

"""

    Bcon = BSSN_constraints[
        CoordSystem
        + ("_rfm_precompute" if enable_rfm_precompute else "")
        + ("_RbarDD_gridfunctions" if enable_RbarDD_gridfunctions else "")
        + ("_T4munu" if enable_T4munu else "")
    ]

    list_of_output_exprs = [Bcon.H]
    Constraints_access_gfs: List[str] = [
        gri.ETLegacyGridFunction.access_gf(gf_name="H")
    ]

    for index in range(3):
        list_of_output_exprs += [Bcon.MU[index]]
        Constraints_access_gfs += [
            gri.ETLegacyGridFunction.access_gf(gf_name="MU" + str(index))
        ]

    list_of_output_exprs += [Bcon.Msquared]
    Constraints_access_gfs += [gri.ETLegacyGridFunction.access_gf(gf_name="MSQUARED")]

    param_symbols, _ = get_params_commondata_symbols_from_expr_list(
        list_of_output_exprs
    )
    if enable_simd:
        body += read_CodeParameters(
            [(thorn_name, symbol) for symbol in param_symbols],
            enable_simd=True,
            declare_invdxxs=False,
        )

    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            list_of_output_exprs,
            Constraints_access_gfs,
            enable_fd_codegen=True,
            enable_simd=enable_simd,
            enable_fd_functions=True,
            enable_GoldenKernels=True,
        ),
        loop_region="interior",
        enable_simd=enable_simd,
        OMP_collapse=OMP_collapse,
    )

    schedule = f"""
if(fd_order == {fd_order}) {{
  schedule FUNC_NAME in MoL_PseudoEvolution as {thorn_name}_BSSN_constraints
  {{
    LANG: C
    READS:  aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF,
            hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
            trKGF, cfGF, lambdaU0GF, lambdaU1GF, lambdaU2GF"""
    if enable_T4munu:
        schedule += """,
            alphaGF, vetU0GF, vetU1GF, vetU2GF,
            T4UU00GF, T4UU01GF, T4UU02GF, T4UU03GF"""
    schedule += f"""
    WRITES: aux_variables
  }} "Compute BSSN (Hamiltonian and momentum) constraints, at finite-differencing order {fd_order}"
}}
"""

    params = param_symbols or None

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=includes,
        desc=desc,
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        prefunc=fin.construct_FD_functions_prefunc().replace(
            "NO_INLINE", "CCTK_ATTRIBUTE_NOINLINE"
        ),  # This prevents a hang when compiling higher-order FD kernels with certain versions of GCC. I'd prefer not adjusting construct_FD_functions_prefunc() for just this infrastructure.
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[("MoL_PseudoEvolution", schedule)],
        ET_current_thorn_CodeParams_used=params,
    )

    # Reset to the initial values
    par.set_parval_from_str("fd_order", old_fd_order)

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

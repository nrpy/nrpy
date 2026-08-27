"""
Generate function enforcing conformal metric determinant and A trace constraints.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.infrastructures.ETLegacy.simple_loop as lp
from nrpy.equations.general_relativity.BSSN_algebraic_constraints import (
    BSSN_algebraic_constraints,
)
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)


def register_CFunction_enforce_detgbar_equals_detghat_trAzero(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    fd_orders: List[int],
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register combined determinant and trace-free conformal-tensor projection.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param fd_orders: Finite-difference orders of registered ADM-to-BSSN producers.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    desc = r"""Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0."""
    name = f"{thorn_name}_enforce_detgbar_equals_detghat_trAzero"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    hprimeDD, aprimeDD = BSSN_algebraic_constraints(CoordSystem, enable_rfm_precompute)

    access_gfs = [
        gri.ETLegacyGridFunction.access_gf(gf_name=f"{basename}{i}{j}")
        for basename in ("hDD", "aDD")
        for i in range(3)
        for j in range(i, 3)
    ]
    projected_exprs = [
        expressions[i][j]
        for expressions in (hprimeDD, aprimeDD)
        for i in range(3)
        for j in range(i, 3)
    ]

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            projected_exprs,
            access_gfs,
            enable_simd=False,
            automatically_read_gf_data_from_memory=True,
            enable_fd_codegen=True,
            enable_fd_functions=True,
        ),
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    reads = """hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF"""
    writes = """hDD00GF(everywhere), hDD01GF(everywhere), hDD02GF(everywhere),
          hDD11GF(everywhere), hDD12GF(everywhere), hDD22GF(everywhere),
          aDD00GF(everywhere), aDD01GF(everywhere), aDD02GF(everywhere),
          aDD11GF(everywhere), aDD12GF(everywhere), aDD22GF(everywhere)"""
    initial_producer_names = [
        f"{thorn_name}_ADM_to_BSSN_order_{fd_order}" for fd_order in fd_orders
    ]
    # Cactus requires a parenthesized whitespace-separated list when an AFTER
    # clause names multiple schedule items. A bare list fails CST parsing.
    initial_producers = " ".join(initial_producer_names)
    if len(initial_producer_names) > 1:
        initial_producers = f"({initial_producers})"
    schedule_initial = (
        "CCTK_INITIAL",
        f"""
schedule FUNC_NAME in CCTK_INITIAL after {initial_producers}
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project freshly initialized or restored BSSN/fCCZ4 state"
""",
    )
    schedule_poststep = (
        "MoL_PostStep",
        f"""
schedule FUNC_NAME in MoL_PostStep after {thorn_name}_evol_ApplyBCs
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project every updated BSSN/fCCZ4 state before its next RHS"
""",
    )
    schedule_pseudoevolution = (
        "MoL_PseudoEvolution",
        f"""
schedule FUNC_NAME in MoL_PseudoEvolution before {thorn_name}_BSSN_constraints
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project BSSN/fCCZ4 state before pseudo-evolution consumers"
""",
    )

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=define_standard_includes(),
        desc=desc,
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[
            schedule_initial,
            schedule_poststep,
            schedule_pseudoevolution,
        ],
    )
    return pcg.NRPyEnv()

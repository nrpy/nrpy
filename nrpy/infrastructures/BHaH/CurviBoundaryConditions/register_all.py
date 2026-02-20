# nrpy/infrastructures/BHaH/CurviBoundaryConditions/register_all.py
"""
Registers all C functions and data for curvilinear boundary conditions.

This module centralizes the registration of all C code components required
for handling curvilinear boundary conditions within the NRPy BHaH
infrastructure. The main function `register_C_functions` sets up:

- Functions to apply various boundary conditions (radiation, extrapolation, inner).
- The `bc_struct` data structure to manage boundary data at runtime.
- C preprocessor definitions for grid function parity.

This process is documented in the tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

from typing import Set

import nrpy.params as par
from nrpy.infrastructures import BHaH


def register_C_functions(
    set_of_CoordSystems: Set[str],
    radiation_BC_fd_order: int = 2,
    set_parity_on_aux: bool = False,
    set_parity_on_auxevol: bool = False,
    enable_masks: bool = False,
) -> None:
    """
    Register various C functions responsible for handling boundary conditions.

    :param set_of_CoordSystems: Set of coordinate systems to use.
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param set_parity_on_aux: If True, set parity on auxiliary grid functions.
    :param set_parity_on_auxevol: If True, set parity on auxiliary evolution grid functions.
    :param enable_masks: If True, make bcstruct algorithm mask-aware.
    """
    _ = par.CodeParameter(
        "char[50]", __name__, "outer_bc_type", "radiation", commondata=True
    )

    if par.parval_from_str("parallelization") == "cuda":
        BHaH.parallelization.cuda_utilities.register_CFunction_cpyHosttoDevice_bc_struct()

    for CoordSystem in set_of_CoordSystems:
        # Register C function to set up the boundary condition struct.
        BHaH.CurviBoundaryConditions.bcstruct_set_up.register_CFunction_bcstruct_set_up(
            CoordSystem=CoordSystem,
            enable_masks=enable_masks,
        )

        # Register C function to apply boundary conditions to both pure outer and inner boundary points.
        BHaH.CurviBoundaryConditions.apply_bcs_outerradiation_and_inner.register_CFunction_apply_bcs_outerradiation_and_inner(
            CoordSystem=CoordSystem,
            radiation_BC_fd_order=radiation_BC_fd_order,
        )

    # Register C function to apply boundary conditions to inner-only boundary points.
    BHaH.CurviBoundaryConditions.apply_bcs_inner_only.register_CFunction_apply_bcs_inner_only()
    BHaH.CurviBoundaryConditions.apply_bcs_inner_only_specific_gfs.register_CFunction_apply_bcs_inner_only_specific_gfs()

    # Register C function to apply boundary conditions to outer-extrapolated and inner boundary points.
    BHaH.CurviBoundaryConditions.apply_bcs_outerextrap_and_inner.register_CFunction_apply_bcs_outerextrap_and_inner()

    # Register bcstruct's contribution to griddata_struct and commondata_struct:
    BHaH.griddata_commondata.register_griddata_commondata(
        __name__,
        "bc_struct bcstruct",
        "all data needed to apply boundary conditions in curvilinear coordinates",
    )

    # Register bcstruct's contribution to BHaH_defines.h:
    BHaH.CurviBoundaryConditions.BHaH_defines.register_BHaH_defines_h(
        set_parity_on_aux=set_parity_on_aux, set_parity_on_auxevol=set_parity_on_auxevol
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

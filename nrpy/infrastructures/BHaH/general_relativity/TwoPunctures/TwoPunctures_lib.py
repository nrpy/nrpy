"""
Register all TwoPunctures functions within NRPy's CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures import BHaH


def register_C_functions(enable_xy_plane: bool = False) -> None:
    """
    Register all C functions needed for the TwoPunctures solve/interpolation path.

    With ``enable_xy_plane=True``, the physical binary stays in the fixed-frame ``xy``
    orientation expected by NRPyPN/TwoPunctures. ``TP_Interp`` then emits native
    fixed-basis components for fixed grids and the legacy startup component ordering
    for rotating startup grids.
    """
    BHaH.general_relativity.TwoPunctures.ID_persist_struct.register_CFunction_initialize_ID_persist_struct(
        enable_xy_plane=enable_xy_plane
    )

    # Original TwoPunctures functions:
    BHaH.general_relativity.TwoPunctures.CoordTransf.register_CFunction_TP_CoordTransf()
    BHaH.general_relativity.TwoPunctures.Equations.register_CFunction_TP_Equations()
    BHaH.general_relativity.TwoPunctures.FuncAndJacobian.register_CFunction_TP_FuncAndJacobian()
    BHaH.general_relativity.TwoPunctures.Newton.register_CFunction_TP_Newton()
    BHaH.general_relativity.TwoPunctures.TP_interp.register_CFunction_TP_Interp(
        enable_xy_plane=enable_xy_plane
    )
    BHaH.general_relativity.TwoPunctures.TP_solve.register_CFunction_TP_solve()
    BHaH.general_relativity.TwoPunctures.TP_utilities.register_CFunction_TP_utilities()

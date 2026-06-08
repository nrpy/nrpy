"""
Register all TwoPunctures functions within NRPy's CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures import BHaH


def register_C_functions_explicit(tp_orientation: str) -> None:
    """
    Register all C functions needed for the TwoPunctures solve/interpolation path.

    This explicit entry point exposes the real shared TwoPunctures orientation
    choices without overloading them with rotating-grid behavior.

    :param tp_orientation: TwoPunctures orientation choice, either
        ``"native_cartesian_xy_plane"`` or ``"legacy_swap_xz"``.
    """
    BHaH.general_relativity.TwoPunctures.ID_persist_struct.register_CFunction_initialize_ID_persist_struct(
        tp_orientation=tp_orientation
    )

    # Original TwoPunctures functions:
    BHaH.general_relativity.TwoPunctures.CoordTransf.register_CFunction_TP_CoordTransf()
    BHaH.general_relativity.TwoPunctures.Equations.register_CFunction_TP_Equations()
    BHaH.general_relativity.TwoPunctures.FuncAndJacobian.register_CFunction_TP_FuncAndJacobian()
    BHaH.general_relativity.TwoPunctures.Newton.register_CFunction_TP_Newton()
    BHaH.general_relativity.TwoPunctures.TP_interp.register_CFunction_TP_Interp(
        tp_orientation=tp_orientation
    )
    BHaH.general_relativity.TwoPunctures.TP_solve.register_CFunction_TP_solve()
    BHaH.general_relativity.TwoPunctures.TP_utilities.register_CFunction_TP_utilities()


def register_C_functions(enable_xy_plane: bool = False) -> None:
    """
    Register all C functions needed for the TwoPunctures solve/interpolation path.

    This compatibility shim maps the older boolean interface onto the explicit
    shared TwoPunctures orientation names.

    :param enable_xy_plane: Whether to keep the physical binary in the native fixed-frame
        xy-plane orientation instead of using the legacy ``swap_xz`` convention.
    """
    register_C_functions_explicit(
        tp_orientation=(
            "native_cartesian_xy_plane" if enable_xy_plane else "legacy_swap_xz"
        )
    )

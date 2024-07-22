"""
Register all TwoPunctures functions within NRPy+'s CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import (
    CoordTransf,
    Equations,
    FuncAndJacobian,
    ID_persist_struct,
    Newton,
    TP_interp,
    TP_solve,
    TP_utilities,
)


def register_C_functions() -> None:
    """Register all C functions needed for TwoPunctures solve."""
    ID_persist_struct.register_CFunction_initialize_ID_persist_struct()

    # Original TwoPunctures functions:
    CoordTransf.register_CFunction_TP_CoordTransf()
    Equations.register_CFunction_TP_Equations()
    FuncAndJacobian.register_CFunction_TP_FuncAndJacobian()
    Newton.register_CFunction_TP_Newton()
    TP_interp.register_CFunction_TP_Interp()
    TP_solve.register_CFunction_TP_solve()
    TP_utilities.register_CFunction_TP_utilities()

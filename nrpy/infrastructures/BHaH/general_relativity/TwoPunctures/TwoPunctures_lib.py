"""
Register all TwoPunctures functions within NRPy+'s CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from pathlib import Path
import shutil

from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import ID_persist_struct
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import CoordTransf
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Equations
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import FuncAndJacobian
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import Newton
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_interp
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_solve
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import TP_utilities


def copy_TwoPunctures_header_files(TwoPunctures_Path: Path) -> None:
    """
    Copy TwoPunctures.h and TP_utilities.h into project directory.

    :param project_Path: The path of the project directory where the file will be copied.
    """
    TwoPunctures_Path.mkdir(parents=True, exist_ok=True)

    try:
        # only Python 3.7+ has importlib.resources
        from importlib import resources  # pylint: disable=E1101,C0415

        for header_file in ["TwoPunctures.h", "TP_utilities.h"]:
            source_path = (
                # pylint: disable=E1101
                resources.files(
                    "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures"
                )
                / header_file
            )
            shutil.copy(str(source_path), str(TwoPunctures_Path / header_file))
    except ImportError:  # Fallback to resource_filename for older Python versions
        # pylint: disable=E1101,C0415
        from pkg_resources import resource_filename  # type: ignore

        for header_file in ["TwoPunctures.h", "TP_utilities.h"]:
            source_path = resource_filename(
                "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures",
                header_file,
            )
            shutil.copy(source_path, str(TwoPunctures_Path))  # type: ignore


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

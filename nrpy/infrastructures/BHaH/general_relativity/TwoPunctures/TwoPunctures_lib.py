"""
Register all TwoPunctures functions within NRPy+'s CFunction dictionary.

License: Lesser GNU Public License, version 2.0+

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path
import shutil

# Attempt to determine the correct files function to use based on Python version
try:
    from importlib.resources import files as resource_files  # Python 3.9 and newer
except ImportError:
    # Fallback for older Python versions: use the backport
    from importlib_resources import files as resource_files

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
    Copy 'TwoPunctures.h' and 'TP_utilities.h' into the specified directory.

    This function uses the appropriate method based on Python version to copy
    the header files from the 'nrpy.infrastructures.BHaH.general_relativity.TwoPunctures'
    package to a specified path.

    :param TwoPunctures_Path: The path of the directory where the files will be copied.
    """
    # Ensure the target directory exists
    TwoPunctures_Path.mkdir(parents=True, exist_ok=True)

    header_files = ["TwoPunctures.h", "TP_utilities.h"]
    package = "nrpy.infrastructures.BHaH.general_relativity.TwoPunctures"

    for header_file in header_files:
        # Use the previously determined files function for resource access
        source_path = resource_files(package) / header_file
        shutil.copy(str(source_path), str(TwoPunctures_Path / header_file))


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

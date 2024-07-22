"""
Copy 'superB.h' into the specified directory.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import shutil
from pathlib import Path

# Try to import the 'files' function from 'importlib.resources' for Python 3.9 and newer versions.
# This provides a consistent API for accessing package resources.
try:
    from importlib.resources import files as resource_files  # Python 3.9 and newer
except ImportError:
    # Fallback to 'importlib_resources' for older Python versions (pre-3.9) to maintain compatibility.
    from importlib_resources import files as resource_files


def copy_superB_header_files(superB_Path: Path) -> None:
    """
    Copy superB.h into the specified directory.

    This function uses the appropriate method based on Python version to copy
    the superB.h file from the 'nrpy.infrastructures.superB.superB' package to a specified path.

    :param superB_Path: The path of the directory where the file will be copied.
    """
    # Ensure the target directory exists
    superB_Path.mkdir(parents=True, exist_ok=True)

    header_file = "superB.h"
    package = "nrpy.infrastructures.superB.superB"

    # Use the previously determined files function for resource access
    source_path = resource_files(package) / header_file
    shutil.copy(str(source_path), str(superB_Path / header_file))

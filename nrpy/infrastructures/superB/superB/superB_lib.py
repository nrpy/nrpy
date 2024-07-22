"""
Copy 'superB.h' into the specified directory.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import shutil
import sys
from pathlib import Path


def copy_superB_header_files(superB_Path: Path) -> None:
    """
    Copy superB.h into the specified directory.

    :param superB_Path: The path of the directory where the files will be copied.
    :raises ImportError: If the importlib.resources module is not found, indicating
                         a Python version earlier than 3.7 where this module was introduced,
                         or if pkg_resources is not available.
    """
    superB_Path.mkdir(parents=True, exist_ok=True)

    try:
        # Only Python 3.9+ has importlib.resources.files()
        if sys.version_info >= (3, 9):
            import importlib.resources  # pylint: disable=E1101,C0415
            from importlib.abc import Traversable  # pylint: disable=C0415,W4904

            for header_file in ["superB.h"]:
                source_path: Traversable = (
                    importlib.resources.files("nrpy.infrastructures.superB.superB")
                    / header_file
                )
                shutil.copy(str(source_path), str(superB_Path / header_file))
        else:
            raise ImportError("Using fallback for importlib.resources")
    except (ImportError, AttributeError):
        # Fallback for versions without importlib.resources or importlib.resources.files()
        try:
            from pkg_resources import resource_filename  # pylint: disable=C0415

            for header_file in ["superB.h"]:
                source_path_str: str = resource_filename(
                    "nrpy.infrastructures.superB.superB", header_file
                )
                shutil.copy(source_path_str, str(superB_Path / header_file))
        except ImportError as e:
            raise ImportError("pkg_resources is not available") from e

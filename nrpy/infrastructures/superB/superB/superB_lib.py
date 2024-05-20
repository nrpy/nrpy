"""
Copy 'superB.h' into the specified directory.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* edu
"""

import shutil
from pathlib import Path


def copy_superB_header_files(superB_Path: Path) -> None:
    """
    Copy superB.h.

    :param superB_Path: The path of the directory where the files will be copied.
    """
    superB_Path.mkdir(parents=True, exist_ok=True)

    try:
        # only Python 3.7+ has importlib.resources
        from importlib import resources  # pylint: disable=E1101,C0415

        for header_file in ["superB.h"]:
            source_path = (
                # pylint: disable=E1101
                resources.files("nrpy.infrastructures.superB.superB")
                / header_file
            )
            shutil.copy(str(source_path), str(superB_Path / header_file))
    except ImportError:  # Fallback to resource_filename for older Python versions
        # pylint: disable=E1101,C0415
        from pkg_resources import resource_filename  # type: ignore

        for header_file in ["superB.h"]:
            source_path = resource_filename(
                "nrpy.infrastructures.superB.superB",
                header_file,
            )
            print("Source path:", source_path)
            shutil.copy(source_path, str(superB_Path))  # type: ignore

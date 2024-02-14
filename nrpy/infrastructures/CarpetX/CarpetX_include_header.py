"""
Define the standard list of headers to include for the CarpetX infrastructure.

Author: Samuel Cupp
"""

from typing import List


def define_standard_includes() -> List[str]:
    """Define the standard list of headers to include for the CarpetX infrastructure."""
    return [
        "loop_device.hxx",
        "math.h",
        "cctk.h",
        "cctk_Arguments.h",
        "cctk_Parameters.h",
    ]

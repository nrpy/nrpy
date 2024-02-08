"""
Define the standard list of headers to include for the ETLegacy infrastructure.

Author: Samuel Cupp
"""

from typing import List


def define_standard_includes() -> List[str]:
    """Define the standard list of headers to include for the ETLegacy infrastructure."""
    return ["math.h", "cctk.h", "cctk_Arguments.h", "cctk_Parameters.h"]

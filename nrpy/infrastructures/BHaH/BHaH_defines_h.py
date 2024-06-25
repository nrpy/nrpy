"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
<<<<<<< HEAD
import nrpy.params as par
from nrpy.infrastructures.BHaH import griddata_commondata

# The ordering of core_modules_list is based largely on data structure dependencies.
# E.g., griddata_struct contains bc_struct.  Module names are parellel scheme dependent
core_modules_list = [
    "general",
    "nrpy.infrastructures.BHaH.diagnostics.progress_indicator",
    "commondata_struct",
    "params_struct",
    "finite_difference",
    "reference_metric",
    "nrpy.infrastructures.BHaH.CurviBoundaryConditions.base_CurviBoundaryConditions",
    "nrpy.infrastructures.BHaH.MoLtimestepping.base_MoL",
    "nrpy.infrastructures.BHaH.interpolation.interpolation",
    "grid",
]
=======
from pathlib import Path
from typing import Dict, List, Optional

import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import griddata_commondata
>>>>>>> main


def register_griddata_struct_and_return_griddata_struct_str(
    enable_rfm_precompute: bool = True,
) -> str:
    """
    Register contributions to the griddata struct from different modules and constructs the griddata_struct string.

    :param enable_rfm_precompute: A boolean value reflecting whether reference metric precomputation is enabled.
    :return: A string representing the typedef structure for grid data, including contributions from BHaH modules,
             reference_metric, CurviBoundaryConditions, and masking, if applicable.
    """
    # Step 1: Register griddata_struct contributions from all BHaH modules:
    griddata_commondata.register_griddata_commondata(
        "params",
        "params_struct params",
        "BHaH parameters, generated from NRPy+'s CodeParameters",
    )
    if any("reference_metric" in key for key in sys.modules):
        griddata_commondata.register_griddata_commondata(
            "reference_metric",
            "char CoordSystemname[100]",
            "the name of the CoordSystem (from reference_metric)",
        )
        griddata_commondata.register_griddata_commondata(
            "reference_metric",
            "char gridname[100]",
            "a user-defined alias for describing the grid",
        )
        if enable_rfm_precompute:
            griddata_commondata.register_griddata_commondata(
                "reference_metric",
                "rfm_struct rfmstruct",
                "includes e.g., 1D arrays of reference metric quantities",
            )

    griddata_struct_def = r"""
typedef struct __griddata__ {
  // griddata_struct stores data needed on each grid
  // xx[3] stores the uniform grid coordinates.
  REAL *restrict xx[3];
"""
    for module, item_list in par.glb_extras_dict["griddata_struct"].items():
        griddata_struct_def += f"  // NRPy+ MODULE: {module}\n"
        for item in item_list:
            griddata_struct_def += f"  {item.c_declaration};"
            if item.description != "":
                griddata_struct_def += f"// <- {item.description}"
            griddata_struct_def += "\n"
    griddata_struct_def += "} griddata_struct;\n"

    return griddata_struct_def


def register_BHaH_defines(module: str, BHaH_defines: str) -> None:
    """
    Register to par.glb_extras_dict["BHaH_defines"] contributions from a given module to BHaH_defines.h.

    :param module: The name of the module for which the defines are being registered.
    :param BHaH_defines: The contribution (string) to BHaH_defines.h.
    """
    if "BHaH_defines" not in par.glb_extras_dict:
        par.glb_extras_dict["BHaH_defines"] = {}

    if module not in par.glb_extras_dict["BHaH_defines"].keys():
        par.glb_extras_dict["BHaH_defines"][module] = BHaH_defines


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

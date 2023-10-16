"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
import sys
from pathlib import Path
from typing import Optional, Dict, List

import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format


def register_griddata_struct_and_return_griddata_struct_str() -> str:
    """
    Register contributions to the griddata struct from different modules and constructs the griddata_struct string.

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


def output_BHaH_defines_h(
    project_dir: str,
    additional_includes: Optional[List[str]] = None,
    REAL_means: str = "double",
    enable_simd: bool = True,
    fin_NGHOSTS_add_one_for_upwinding: bool = False,
    supplemental_defines_dict: Optional[Dict[str, str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> None:
    r"""
    Output C code header file with macro definitions and other configurations for the project.

    :param project_dir: Directory where the project C code is output
    :param additional_includes: Additional header files to be included in the output
    :param REAL_means: The floating-point type to be used in the C code (default is "double")
    :param enable_simd: Flag to enable Single Instruction Multiple Data (SIMD) optimizations
    :param fin_NGHOSTS_add_one_for_upwinding: Option to add one extra ghost zone for upwinding
    :param supplemental_defines_dict: Additional key-value pairs to be included in the output file
    :param clang_format_options: Options for clang formatting.

    :raises ValueError: If the project directory does not exist

    >>> from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> MoL.register_CFunctions(register_MoL_step_forward_in_time=False)
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> output_BHaH_defines_h(project_dir=str(project_dir))
    >>> expected_string = decompress_base64_to_string("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")
    >>> returned_string = (project_dir / "BHaH_defines.h").read_text()
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted BHaH_defines.h string changed!\n"
    ...    error_message += "Here's the diff:\n" + diff_strings(expected_string, returned_string) + "\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    ###############################
    # GENERALLY USEFUL DEFINITIONS
    gen_BHd_str = """#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <errno.h>
#include <ctype.h>
#include <time.h>
"""
    if enable_simd:
        gen_BHd_str += "// output_BHaH_defines_h(...,enable_simd=True) was called so we #include SIMD intrinsics:\n"
        gen_BHd_str += """#include "simd/simd_intrinsics.h"\n"""
    if additional_includes:
        for include in additional_includes:
            gen_BHd_str += f'#include "{include}"\n'
    gen_BHd_str += rf"""#define REAL {REAL_means}\

#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
"""
    register_BHaH_defines("general", gen_BHd_str)

    #####################################
    # PARAMS_STRUCT AND COMMONDATA_STRUCT
    # Generate C code to declare C params_struct & commondata_struct;
    #         output to "declare_Cparameters_struct.h"
    #         We want the elements of this struct to be *sorted*,
    #         to ensure that the struct is consistently ordered
    #         for checkpointing purposes.
    par_BHd_str = "typedef struct __params_struct__ {\n"
    commondata_BHd_str = "typedef struct __commondata_struct__ {\n"
    CCodelines_params_struct: List[str] = []
    CCodelines_commondata_struct: List[str] = []

    # Add all CodeParameters
    # Iterate through the global code parameters dictionary
    for CPname, CodeParam in par.glb_code_params_dict.items():
        CPtype: str = CodeParam.c_type_alias
        if CPtype != "#define":
            comment: str = f"  // {CodeParam.module}::{CPname}"
            c_output = f"  {CPtype} {CPname};{comment}\n"
            if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                chararray_size = CPtype.split("[")[1].replace("]", "")
                c_output = f"char {CPname}[{chararray_size}];{comment}\n"

            # Concatenate module and CPname for the comment
            if CodeParam.commondata:
                CCodelines_commondata_struct.append(c_output)
            else:
                CCodelines_params_struct.append(c_output)

    # Sort CCodelines_params_struct and append them to the par_BHd_str
    for line in sorted(CCodelines_params_struct):
        par_BHd_str += line
    # Sort CCodelines_commondata_struct and append them to the commondata_BHd_str
    for line in sorted(CCodelines_commondata_struct):
        commondata_BHd_str += line

    par_BHd_str += "} params_struct;\n"
    commondata_BHd_str += "} commondata_struct;\n"

    register_BHaH_defines("params_struct", par_BHd_str)
    register_BHaH_defines("commondata_struct", commondata_BHd_str)

    ###############################
    # FINITE DIFFERENCE
    # First register C functions needed by finite_difference

    # Then set up the dictionary entry for finite_difference in BHaH_defines
    if any("finite_difference" in key for key in sys.modules):
        NGHOSTS = int(par.parval_from_str("finite_difference::fd_order") / 2)
        if fin_NGHOSTS_add_one_for_upwinding:
            NGHOSTS += 1
        fin_BHd_str = f"""
// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS {NGHOSTS}
"""
        if not enable_simd:
            fin_BHd_str += """
// When enable_simd = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0\n"""
        register_BHaH_defines("finite_difference", fin_BHd_str)

    ###############################
    # GRID, etc.
    # Then set up the dictionary entry for grid in BHaH_defines
    gri_BHd_str = gri.BHaHGridFunction.gridfunction_defines()
    gri_BHd_str += (
        r"""
// Declare the IDX4(gf,i,j,k) macro, which enables us to store 4-dimensions of
//   data in a 1D array. In this case, consecutive values of "i"
//   (all other indices held to a fixed value) are consecutive in memory, where
//   consecutive values of "j" (fixing all other indices) are separated by
//   Nxx_plus_2NGHOSTS0 elements in memory. Similarly, consecutive values of
//   "k" are separated by Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1 in memory, etc.
#define IDX4(g,i,j,k)                                                  \
  ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) + Nxx_plus_2NGHOSTS2 * (g) ) ) )
#define IDX4pt(g,idx) ( (idx) + (Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2) * (g) )
#define IDX3(i,j,k) ( (i) + Nxx_plus_2NGHOSTS0 * ( (j) + Nxx_plus_2NGHOSTS1 * ( (k) ) ) )
#define LOOP_REGION(i0min,i0max, i1min,i1max, i2min,i2max)              \
  for(int i2=i2min;i2<i2max;i2++) for(int i1=i1min;i1<i1max;i1++) for(int i0=i0min;i0<i0max;i0++)
#define LOOP_OMP(__OMP_PRAGMA__, i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max) \
_Pragma(__OMP_PRAGMA__)  \
    for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_NOOMP(i0,i0min,i0max, i1,i1min,i1max, i2,i2min,i2max)      \
  for(int (i2)=(i2min);(i2)<(i2max);(i2)++) for(int (i1)=(i1min);(i1)<(i1max);(i1)++) for(int (i0)=(i0min);(i0)<(i0max);(i0)++)
#define LOOP_BREAKOUT(i0,i1,i2, i0max,i1max,i2max) { i0=(i0max); i1=(i1max); i2=(i2max); break; }
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2, NG) \
  ( i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0)-(NG) &&       \
    i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1)-(NG) &&       \
    i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2)-(NG) )
"""
        + register_griddata_struct_and_return_griddata_struct_str()
    )
    register_BHaH_defines("grid", gri_BHd_str)

    def output_key(key_name: str, item_name: str) -> str:
        return f"""
//********************************************
// Basic definitions for module {key_name}:\n{item_name}"""

    file_output_str = """// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""
    # The ordering here is based largely on data structure dependencies. E.g., griddata_struct contains bc_struct.
    core_modules_list = [
        "general",
        "nrpy.infrastructures.BHaH.diagnostics.progress_indicator",
        "commondata_struct",
        "params_struct",
        "finite_difference",
        "reference_metric",
        "nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions",
        "nrpy.infrastructures.BHaH.MoLtimestepping.MoL",
        "nrpy.infrastructures.BHaH.interpolate",
        "grid",  # griddata struct depends upon other core modules
    ]
    # Populate BHaH_defines.h with core modules
    for core_module in core_modules_list:
        if core_module in par.glb_extras_dict["BHaH_defines"].keys():
            file_output_str += output_key(
                core_module, par.glb_extras_dict["BHaH_defines"][core_module]
            )

    # Populate BHaH_defines.h with non-core modules
    for key, item in sorted(par.glb_extras_dict["BHaH_defines"].items()):
        if key not in core_modules_list:
            print(f"Outputting non-core modules key = {key} to BHaH_defines.h")
            file_output_str += output_key(key, item)

    # Populate BHaH_defines.h with whatever else is desired.
    if supplemental_defines_dict:
        for key in supplemental_defines_dict:
            file_output_str += output_key(key, supplemental_defines_dict[key])

    bhah_defines_file = project_Path / "BHaH_defines.h"
    with bhah_defines_file.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(file_output_str, clang_format_options=clang_format_options)
        )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

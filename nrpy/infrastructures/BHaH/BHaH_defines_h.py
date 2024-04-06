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


def output_BHaH_defines_h(
    project_dir: str,
    additional_includes: Optional[List[str]] = None,
    REAL_means: str = "double",
    enable_simd: bool = True,
    enable_rfm_precompute: bool = True,
    fin_NGHOSTS_add_one_for_upwinding_or_KO: bool = False,
    supplemental_defines_dict: Optional[Dict[str, str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    r"""
    Output C code header file with macro definitions and other configurations for the project.

    :param project_dir: Directory where the project C code is output
    :param additional_includes: Additional header files to be included in the output
    :param REAL_means: The floating-point type to be used in the C code (default is "double")
    :param enable_simd: Flag to enable Single Instruction Multiple Data (SIMD) optimizations
    :param enable_rfm_precompute: A boolean value reflecting whether reference metric precomputation is enabled.
    :param fin_NGHOSTS_add_one_for_upwinding_or_KO: Option to add one extra ghost zone for upwinding
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
    gen_BHd_str = """#include <ctype.h>   // Character type functions, such as isdigit, isalpha, etc.
#include <errno.h>   // Error number definitions
#include <math.h>    // Transcendental functions, etc.
#include <stdbool.h> // bool-typed variables
#include <stdint.h>  // int8_t-typed variables
#include <stdio.h>   // Basic input/output functions, such as *printf, fopen, fwrite, etc.
#include <stdlib.h>  // malloc/free, etc.
#include <string.h>  // String handling functions, such as strlen, strcmp, etc.
#include <time.h>    // Time-related functions and types, such as time(), clock(),
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
#define SQR(A) ((A) * (A))

// START: CodeParameters declared as #define.
"""
    for CPname, CodeParam in par.glb_code_params_dict.items():
        if CodeParam.cparam_type == "#define":
            gen_BHd_str += f"""#ifndef {CPname}
#define {CPname} {CodeParam.defaultvalue}
#endif
// END: CodeParameters declared as #define.
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
        CPtype = CodeParam.cparam_type
        if CPtype != "#define":
            comment = f"  // {CodeParam.module}::{CPname}"
            c_output = f"  {CPtype} {CPname};{comment}\n"
            if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                chararray_size = CPtype.split("[")[1].replace("]", "")
                c_output = f"char {CPname}[{chararray_size}];{comment}\n"

            # Concatenate module and CPname for the comment
            if CodeParam.commondata:
                CCodelines_commondata_struct.append(c_output)
            else:
                CCodelines_params_struct.append(c_output)

    if "commondata_struct" in par.glb_extras_dict:
        for module, item_list in par.glb_extras_dict["commondata_struct"].items():
            for item in item_list:
                c_code_line = f"  {item.c_declaration};"
                if item.description != "":
                    c_code_line += f"// <- {module}: {item.description}"
                CCodelines_commondata_struct.append(c_code_line + "\n")

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
        if fin_NGHOSTS_add_one_for_upwinding_or_KO:
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
    gri_BHd_str += r"""
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
""" + register_griddata_struct_and_return_griddata_struct_str(
        enable_rfm_precompute=enable_rfm_precompute
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
        "nrpy.infrastructures.BHaH.interpolation.interpolation",
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

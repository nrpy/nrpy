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
"""
    code_params_includes_define_type = False
    for CPname, CodeParam in par.glb_code_params_dict.items():
        if CodeParam.cparam_type == "#define":
            if not code_params_includes_define_type:
                code_params_includes_define_type = True
                gen_BHd_str += "// START: CodeParameters declared as #define.\n"
            gen_BHd_str += f"""#ifndef {CPname}
#define {CPname} {CodeParam.defaultvalue} // {CodeParam.module}
#endif
"""
    if code_params_includes_define_type:
        gen_BHd_str += "// END: CodeParameters declared as #define.\n"
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

    file_output_str = """#ifndef __BHAH_DEFINES_H__
#define __BHAH_DEFINES_H__
// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
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

    file_output_str += """
#endif
"""
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

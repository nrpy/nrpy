"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Union

import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH.BHaH_defines_h import (
    register_BHaH_defines,
    register_griddata_struct_and_return_griddata_struct_str,
)

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


class base_output_BHaH_defines_h:
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

    """

    def __init__(
        self,
        project_dir: str,
        additional_includes: Optional[List[str]] = None,
        REAL_means: str = "double",
        enable_simd: bool = True,
        enable_rfm_precompute: bool = True,
        fin_NGHOSTS_add_one_for_upwinding_or_KO: bool = False,
        supplemental_defines_dict: Optional[Dict[str, str]] = None,
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    ) -> None:
        self.project_dir = project_dir
        self.additional_includes = additional_includes
        self.REAL_means = REAL_means
        self.enable_simd = enable_simd
        self.enable_rfm_precompute = enable_rfm_precompute
        self.fin_NGHOSTS_add_one_for_upwinding_or_KO = (
            fin_NGHOSTS_add_one_for_upwinding_or_KO
        )
        self.supplemental_defines_dict = supplemental_defines_dict
        self.clang_format_options = clang_format_options

        self.project_Path = Path(project_dir)
        self.project_Path.mkdir(parents=True, exist_ok=True)
        self.file_output_str: str = ""

        ###############################
        # GENERALLY USEFUL DEFINITIONS
        self.BHd_include_str = """#include <ctype.h>   // Character type functions, such as isdigit, isalpha, etc.
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
            self.BHd_include_str += "// output_BHaH_defines_h(...,enable_simd=True) was called so we #include SIMD intrinsics:\n"
            self.BHd_include_str += """#include "simd/simd_intrinsics.h"\n"""
        if additional_includes:
            for include in additional_includes:
                self.BHd_include_str += f'#include "{include}"\n'

        self.BHd_definitions_str = rf"""#define REAL {REAL_means}\

// These macros for MIN(), MAX(), and SQR() ensure that if the arguments inside
//   are a function/complex expression, the function/expression is evaluated
//   *only once* per argument. See https://lwn.net/Articles/983965/ for details.
// They are improvements over the original implementations:
// #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
// #define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
// #define SQR(A) ((A) * (A))
#define MIN(A, B) ({{ \
    __typeof__(A) _a = (A); \
    __typeof__(B) _b = (B); \
    _a < _b ? _a : _b; \
}})
#define MAX(A, B) ({{ \
    __typeof__(A) _a = (A); \
    __typeof__(B) _b = (B); \
    _a > _b ? _a : _b; \
}})
#define SQR(A) ({{ \
    __typeof__(A) _a = (A); \
    _a * _a; \
}})
"""
        code_params_includes_define_type = False
        for CPname, CodeParam in par.glb_code_params_dict.items():
            if CodeParam.cparam_type == "#define":
                if not code_params_includes_define_type:
                    code_params_includes_define_type = True
                    self.BHd_definitions_str += (
                        "// START: CodeParameters declared as #define.\n"
                    )
                self.BHd_definitions_str += f"""#ifndef {CPname}
#define {CPname} {CodeParam.defaultvalue} // {CodeParam.module}
#endif
"""
        if code_params_includes_define_type:
            self.BHd_definitions_str += "// END: CodeParameters declared as #define.\n"

        #####################################
        # PARAMS_STRUCT AND COMMONDATA_STRUCT
        # Generate C code to declare C params_struct & commondata_struct;
        #         output to "declare_Cparameters_struct.h"
        #         We want the elements of this struct to be *sorted*,
        #         to ensure that the struct is consistently ordered
        #         for checkpointing purposes.
        self.par_BHd_str = "typedef struct __params_struct__ {\n"
        self.commondata_BHd_str = "typedef struct __commondata_struct__ {\n"
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
            self.par_BHd_str += line
        # Sort CCodelines_commondata_struct and append them to the commondata_BHd_str
        for line in sorted(CCodelines_commondata_struct):
            self.commondata_BHd_str += line

        self.par_BHd_str += "} params_struct;\n"
        self.commondata_BHd_str += "} commondata_struct;\n"

        ###############################
        # FINITE DIFFERENCE
        # First register C functions needed by finite_difference

        # Then set up the dictionary entry for finite_difference in BHaH_defines
        self.NGHOSTS: Union[int, None] = None
        if any("finite_difference" in key for key in sys.modules):
            self.NGHOSTS = int(par.parval_from_str("finite_difference::fd_order") / 2)
            if self.fin_NGHOSTS_add_one_for_upwinding_or_KO:
                self.NGHOSTS += 1
            self.fin_BHd_str = f"""
// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS {self.NGHOSTS}
#define NO_INLINE // Account for cases where NO_INLINE might appear in codegen
"""
            if not enable_simd:
                self.fin_BHd_str += """
// When enable_simd = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0\n"""

        ###############################
        # GRID, etc.
        # Then set up the dictionary entry for grid in BHaH_defines
        self.gri_BHd_str = gri.BHaHGridFunction.gridfunction_defines()
        self.gri_BHd_str += r"""
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

        self.gri_BHd_struct_str = (
            register_griddata_struct_and_return_griddata_struct_str(
                enable_rfm_precompute=enable_rfm_precompute
            )
        )

    def register_define_blocks(self) -> None:
        """Register definition blocks to BHaH_defines."""
        general_str = self.BHd_include_str + self.BHd_definitions_str
        register_BHaH_defines("general", general_str)

        register_BHaH_defines("params_struct", self.par_BHd_str)
        register_BHaH_defines("commondata_struct", self.commondata_BHd_str)

        if any("finite_difference" in key for key in sys.modules):
            register_BHaH_defines("finite_difference", self.fin_BHd_str)

        grid_str = self.gri_BHd_str + self.gri_BHd_struct_str
        register_BHaH_defines("grid", grid_str)

    def generate_output_str(self) -> None:
        """Generate final str."""

        def output_key(key_name: str, item_name: str) -> str:
            return f"""
//********************************************
// Basic definitions for module {key_name}:\n{item_name}"""

        self.file_output_str = """#ifndef __BHAH_DEFINES_H__
#define __BHAH_DEFINES_H__
// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""

        # Populate BHaH_defines.h with core modules
        for core_module in core_modules_list:
            if core_module in par.glb_extras_dict["BHaH_defines"].keys():
                self.file_output_str += output_key(
                    core_module, par.glb_extras_dict["BHaH_defines"][core_module]
                )

        # Populate BHaH_defines.h with non-core modules
        for key, item in sorted(par.glb_extras_dict["BHaH_defines"].items()):
            if key not in core_modules_list:
                print(f"Outputting non-core modules key = {key} to BHaH_defines.h")
                self.file_output_str += output_key(key, item)

        # Populate BHaH_defines.h with whatever else is desired.
        if self.supplemental_defines_dict:
            for key in self.supplemental_defines_dict:
                self.file_output_str += output_key(
                    key, self.supplemental_defines_dict[key]
                )
        self.file_output_str += """
#endif
"""

    # Overload this if you need to structure things differently
    def write_to_file(self) -> None:
        """Write final str to header file."""
        bhah_defines_file = self.project_Path / "BHaH_defines.h"
        with bhah_defines_file.open("w", encoding="utf-8") as file:
            file.write(
                clang_format(
                    self.file_output_str, clang_format_options=self.clang_format_options
                )
            )


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

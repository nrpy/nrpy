"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List, Union

import nrpy.params as par
import nrpy.grid as gri
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH.BHaH_defines_h import (
    core_modules_list,
    register_griddata_struct_and_return_griddata_struct_str,
    register_BHaH_defines,
)


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

    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.openmp import MoL
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> MoL.register_CFunctions(register_MoL_step_forward_in_time=False)
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> C=output_BHaH_defines_h(project_dir=str(project_dir))
    >>> expected_string = decompress_base64_to_string("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")
    >>> returned_string = (project_dir / "BHaH_defines.h").read_text()
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted BHaH_defines.h string changed!\n"
    ...    error_message += "Here's the diff:\n" + diff_strings(expected_string, returned_string) + "\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
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

#define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
#define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
#define SQR(A) ((A) * (A))
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

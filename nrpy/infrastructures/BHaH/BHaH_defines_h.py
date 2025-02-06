"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional

import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import griddata_commondata

# The ordering of core_modules_list is based largely on data structure dependencies.
# E.g., griddata_struct contains bc_struct.  Module names are parellel scheme dependent
# Note: These are here such that they can be modified as needed by relevant modules (e.g. MoL)
core_modules_list = [
    "general",
    "after_general",
    "nrpy.infrastructures.BHaH.diagnostics.progress_indicator",
    "commondata_struct",
    "params_struct",
    "finite_difference",
    "reference_metric",
    "nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions",
    "nrpy.infrastructures.BHaH.MoLtimestepping.MoL_register_all",
    "nrpy.infrastructures.BHaH.interpolation.interpolation",
    "grid",
]


def register_griddata_struct_and_return_griddata_struct_str(
    enable_rfm_precompute: bool = True,
) -> str:
    """
    Register contributions to the griddata struct from different modules and construct the griddata_struct string.

    :param enable_rfm_precompute: A boolean value reflecting whether reference metric precomputation is enabled.
    :return: A string representing the typedef structure for grid data, including contributions from BHaH modules,
             reference_metric, CurviBoundaryConditions, and masking, if applicable.

    DocTests:
    >>> result = register_griddata_struct_and_return_griddata_struct_str()
    >>> isinstance(result, str)
    True
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
                "rfm_struct* rfmstruct",
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
    Register contributions from a given module to par.glb_extras_dict["BHaH_defines"] for BHaH_defines.h.

    :param module: The name of the module for which the defines are being registered.
    :param BHaH_defines: The contribution (string) to BHaH_defines.h.

    DocTests:
    >>> par.glb_extras_dict = {}
    >>> register_BHaH_defines('test_module', '#define TEST_MACRO 1')
    >>> 'test_module' in par.glb_extras_dict['BHaH_defines']
    True
    """
    if "BHaH_defines" not in par.glb_extras_dict:
        par.glb_extras_dict["BHaH_defines"] = {}

    if module not in par.glb_extras_dict["BHaH_defines"].keys():
        par.glb_extras_dict["BHaH_defines"][module] = BHaH_defines
    else:
        par.glb_extras_dict["BHaH_defines"][module] += BHaH_defines


def output_BHaH_defines_h(
    project_dir: str,
    additional_includes: Optional[List[str]] = None,
    restrict_pointer_type: str = "*restrict",
    enable_intrinsics: bool = True,
    intrinsics_header_lst: List[str] = ["simd_intrinsics.h"],
    define_no_simd_UPWIND_ALG: bool = True,
    enable_rfm_precompute: bool = True,
    fin_NGHOSTS_add_one_for_upwinding_or_KO: bool = False,
    DOUBLE_means: str = "double",
    supplemental_defines_dict: Optional[Dict[str, str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    r"""
    Output C code header file with macro definitions and other configurations for the project.

    :param project_dir: Directory where the project C code is output
    :param additional_includes: Additional header files to be included in the output
    :param restrict_pointer_type: Allow modifications of restrict pointer type.  Default is *restrict
    :param enable_intrinsics: Flag to enable hardware intrinsics
    :param intrinsics_header_lst: List of intrinsics header files
    :param define_no_simd_UPWIND_ALG: Flag to #define a SIMD-less UPWIND_ALG. No need to define this if UPWIND_ALG() unused.
    :param enable_rfm_precompute: A boolean value reflecting whether reference metric precomputation is enabled.
    :param fin_NGHOSTS_add_one_for_upwinding_or_KO: Option to add one extra ghost zone for upwinding
    :param DOUBLE_means: Overload DOUBLE macro type for specific calculations that require higher than single precision
    :param supplemental_defines_dict: Additional key-value pairs to be included in the output file
    :param clang_format_options: Options for clang formatting.

    DocTests:
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping import MoL_register_all
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import validate_strings
    >>> MoL_register_all.register_CFunctions(register_MoL_step_forward_in_time=False)
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> output_BHaH_defines_h(project_dir=str(project_dir))
    >>> validate_strings((project_dir / "BHaH_defines.h").read_text(), "BHaH_defines")
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
    if enable_intrinsics:
        gen_BHd_str += "// output_BHaH_defines_h(...,enable_intrinsics=True) was called so we intrinsics headers:\n"
        for intr_header in intrinsics_header_lst:
            gen_BHd_str += f"""#include "intrinsics/{intr_header}"\n"""
    if additional_includes:
        for include in additional_includes:
            gen_BHd_str += f'#include "{include}"\n'
    REAL_means = par.parval_from_str("fp_type")
    gen_BHd_str += f"""#define REAL {REAL_means}
#define DOUBLE {DOUBLE_means}

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
#ifndef MAYBE_UNUSED
#if __cplusplus >= 201703L
#define MAYBE_UNUSED [[maybe_unused]]
#elif defined(__GNUC__) || defined(__clang__) || defined(__NVCC__)
#define MAYBE_UNUSED __attribute__((unused))
#else
#define MAYBE_UNUSED
#endif // END check for GCC, Clang, or NVCC
#endif // END MAYBE_UNUSED
"""

    code_params_includes_define_type = False
    for CPname, CodeParam in par.glb_code_params_dict.items():
        CPtype = CodeParam.cparam_type
        if CPtype == "#define":
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

    # Helper function to format C declarations without using regex
    def format_c_declaration(
        cp_type: str, var_name: str, module: str, description: str
    ) -> str:
        """
        Given a CodeParameter type, variable name, module, and description, return the correct C declaration string.
        Handles both scalar and array types using simple string operations.

        :param cp_type: The type of the C parameter (e.g., "REAL" or "int[8]").
        :param var_name: The name of the variable.
        :param module: The module name to be included in the comment.
        :param description: A description of the parameter to be appended to the comment.
                            If empty, only the module and variable name are included in the comment.

        :return: A formatted C declaration string with an inline comment.
        """
        if "[" in cp_type and "]" in cp_type:
            base_type, size_with_bracket = cp_type.split("[", 1)
            size = size_with_bracket.split("]", 1)[0]
            base_type = base_type.strip()
            size = size.strip()
            if base_type.startswith("char"):
                # Handle char arrays
                decl = f"  char {var_name}[{size}];"
            else:
                decl = f"  {base_type} {var_name}[{size}];"
        else:
            decl = f"  {cp_type} {var_name};"

        # Conditional comment based on description
        if description:
            comment = f" // {description} ({module})"
        else:
            comment = f" // ({module})"

        return decl + comment + "\n"

    # Add all CodeParameters
    # Iterate through the global code parameters dictionary
    for CPname, CodeParam in par.glb_code_params_dict.items():
        CPtype = CodeParam.cparam_type
        if CPtype != "#define":
            # All CodeParams have a description field
            description = CodeParam.description.strip()
            module = CodeParam.module
            c_declaration = format_c_declaration(CPtype, CPname, module, description)
            # Append c_declaration to the appropriate structure: commondata or params
            if CodeParam.commondata:
                CCodelines_commondata_struct.append(c_declaration)
            else:
                CCodelines_params_struct.append(c_declaration)

    if "commondata_struct" in par.glb_extras_dict:
        for module, item_list in par.glb_extras_dict["commondata_struct"].items():
            for item in item_list:
                c_code_line = f"  {item.c_declaration};"
                if item.description != "":
                    c_code_line += f" // <- {module}: {item.description}"
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

// Declare NO_INLINE macro, used in FD functions. GCC v10+ compilations hang on complex RHS expressions (like BSSN) without this.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
    #define NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
    #define NO_INLINE __declspec(noinline)
#else
    #define NO_INLINE // Fallback for unknown compilers
#endif
"""
        if not enable_intrinsics and define_no_simd_UPWIND_ALG:
            fin_BHd_str += """
// When enable_intrinsics = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
"""
        register_BHaH_defines("finite_difference", fin_BHd_str)

    ###############################
    # GRID, etc.
    # Then set up the dictionary entry for grid in BHaH_defines
    gri_BHd_str = gri.BHaHGridFunction.gridfunction_defines()
    gri_BHd_str += r"""
// ----------------------------
// Indexing macros
// ----------------------------
// IDX4: Converts 4D grid indices (gf, i, j, k) into a 1D array index using the strides
//       Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, and Nxx_plus_2NGHOSTS2. This macro assumes
//       that the "i" index varies fastest in memory.
#define IDX4(gf, i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k) + Nxx_plus_2NGHOSTS2 * (gf))))
// IDX4P: Similar to IDX4, but retrieves grid dimensions from the provided parameter structure
//        "params" instead of using global variables.
#define IDX4P(params, gf, i, j, k)                                                                                                                   \
  ((i) + (params)->Nxx_plus_2NGHOSTS0 * ((j) + (params)->Nxx_plus_2NGHOSTS1 * ((k) + (params)->Nxx_plus_2NGHOSTS2 * (gf))))
// IDX4pt: Computes the 1D index offset for a given grid function index (gf) based on an existing index (idx)
//         by using the total number of elements in one grid function, defined as the product of the grid strides.
#define IDX4pt(gf, idx) ((idx) + (Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2) * (gf))
// IDX3: Converts 3D grid indices (i, j, k) into a 1D array index using the strides Nxx_plus_2NGHOSTS0
//       and Nxx_plus_2NGHOSTS1. Like IDX4, this macro assumes the "i" index varies fastest.
#define IDX3(i, j, k) ((i) + Nxx_plus_2NGHOSTS0 * ((j) + Nxx_plus_2NGHOSTS1 * ((k))))
// IDX3P: Similar to IDX3, but retrieves grid dimensions from the provided parameter structure "params".
#define IDX3P(params, i, j, k) ((i) + (params)->Nxx_plus_2NGHOSTS0 * ((j) + (params)->Nxx_plus_2NGHOSTS1 * ((k))))
// END: Indexing macros

// ----------------------------
// Loop-related macros
// ----------------------------
// SET_NXX_PLUS_2NGHOSTS_VARS: Declares local constants for the grid dimensions (including ghost zones) by extracting
// the values from griddata[whichgrid].params.
#define SET_NXX_PLUS_2NGHOSTS_VARS(whichgrid)                                                                                                        \
  const int Nxx_plus_2NGHOSTS0 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS0;                                                                      \
  const int Nxx_plus_2NGHOSTS1 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS1;                                                                      \
  const int Nxx_plus_2NGHOSTS2 = griddata[whichgrid].params.Nxx_plus_2NGHOSTS2;
// LOOP_REGION: Iterates over a 3D region defined by the inclusive lower bounds (i0min, i1min, i2min)
// and exclusive upper bounds (i0max, i1max, i2max) for each dimension.
#define LOOP_REGION(i0min, i0max, i1min, i1max, i2min, i2max)                                                                                        \
  for (int i2 = i2min; i2 < i2max; i2++)                                                                                                             \
    for (int i1 = i1min; i1 < i1max; i1++)                                                                                                           \
      for (int i0 = i0min; i0 < i0max; i0++)
// LOOP_OMP: Similar to LOOP_REGION but inserts an OpenMP pragma (via __OMP_PRAGMA__) for parallelization.
#define LOOP_OMP(__OMP_PRAGMA__, i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                               \
  _Pragma(__OMP_PRAGMA__) for (int(i2) = (i2min); (i2) < (i2max); (i2)++) for (int(i1) = (i1min); (i1) < (i1max);                                    \
                                                                              (i1)++) for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
// LOOP_NOOMP: A non-parallel version of the 3D loop, identical in structure to LOOP_REGION.
#define LOOP_NOOMP(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                                             \
  for (int(i2) = (i2min); (i2) < (i2max); (i2)++)                                                                                                    \
    for (int(i1) = (i1min); (i1) < (i1max); (i1)++)                                                                                                  \
      for (int(i0) = (i0min); (i0) < (i0max); (i0)++)
// LOOP_BREAKOUT: Forces an exit from the nested loops by setting the loop indices to their maximum values and executing a break.
#define LOOP_BREAKOUT(i0, i1, i2, i0max, i1max, i2max)                                                                                               \
  {                                                                                                                                                  \
    i0 = (i0max);                                                                                                                                    \
    i1 = (i1max);                                                                                                                                    \
    i2 = (i2max);                                                                                                                                    \
    break;                                                                                                                                           \
  }
// IS_IN_GRID_INTERIOR: Checks whether the provided 3D index array (i0i1i2) lies within the grid interior,
// defined as the region excluding NG ghost cells on each boundary.
#define IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NG)                                                  \
  (i0i1i2[0] >= (NG) && i0i1i2[0] < (Nxx_plus_2NGHOSTS0) - (NG) && i0i1i2[1] >= (NG) && i0i1i2[1] < (Nxx_plus_2NGHOSTS1) - (NG) &&                   \
   i0i1i2[2] >= (NG) && i0i1i2[2] < (Nxx_plus_2NGHOSTS2) - (NG))

// ----------------------------
// Define griddata struct
// ----------------------------"""
    griddata_struct_def = register_griddata_struct_and_return_griddata_struct_str(
        enable_rfm_precompute=enable_rfm_precompute
    )
    # Add a newline as needed:
    if not griddata_struct_def.startswith("\n"):
        griddata_struct_def = "\n" + griddata_struct_def
    gri_BHd_str += griddata_struct_def

    register_BHaH_defines("grid", gri_BHd_str)

    def output_key(key_name: str, item_name: str) -> str:
        """
        Format a string for outputting the definitions for a specific module.

        :param key_name: The name of the module or key.
        :param item_name: The definitions or content associated with the module.
        :return: A formatted string containing the module name and its content.
        """
        return f"""
// ----------------------------
// Basic definitions for module
// {key_name}:
// ----------------------------
{item_name}"""

    file_output_str = """#ifndef __BHAH_DEFINES_H__
#define __BHAH_DEFINES_H__
// BHaH core header file, automatically generated from output_BHaH_defines_h within BHaH_defines_h.py,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""

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

    file_output_str = file_output_str.replace("*restrict", restrict_pointer_type)

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

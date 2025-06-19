"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import nrpy.grid as gri
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH import griddata_commondata

# The ordering of core_modules_list is based largely on data structure dependencies.
# E.g., griddata_struct contains bc_struct. Module names are parallel scheme dependent.
# Note: These are here such that they can be modified as needed by relevant modules (e.g., MoL).
core_modules_list = [
    "general",
    "after_general",
    "nrpy.infrastructures.BHaH.diagnostics.progress_indicator",
    "commondata_struct",
    "params_struct",
    "finite_difference",
    "reference_metric",
    "nrpy.infrastructures.BHaH.CurviBoundaryConditions.BHaH_defines",
    "nrpy.infrastructures.BHaH.MoLtimestepping.BHaH_defines",
    "nrpydev.infrastructures.BHaH.interpolation.interpolation",
    "grid",
]


def register_griddata_struct_and_return_griddata_struct_str(
    enable_rfm_precompute: bool = True,
) -> str:
    """
    Register contributions to the griddata struct and construct its C-string definition.

    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :return: A string for the typedef of griddata_struct.

    Doctests:
    >>> result = register_griddata_struct_and_return_griddata_struct_str()
    >>> isinstance(result, str)
    True
    """
    # Register griddata_struct contributions from BHaH modules.
    griddata_commondata.register_griddata_commondata(
        "params",
        "params_struct params",
        "BHaH parameters, generated from NRPy+'s CodeParameters",
    )
    if any("reference_metric" in key for key in sys.modules):
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
            if item.description:
                griddata_struct_def += f"// <- {item.description}"
            griddata_struct_def += "\n"
    griddata_struct_def += "} griddata_struct;\n"

    return griddata_struct_def


def parse_cparam_type(cparam_type: str) -> Tuple[str, Optional[str], bool]:
    """
    Parse a cparam_type string into its base type, size, and array status.

    :param cparam_type: The raw CParam type string, e.g., "REAL[8]".
    :return: A tuple (base, size, is_array).
    :raises ValueError: If the array size is not a numeric string.

    Doctests:
    >>> parse_cparam_type("int")
    ('int', None, False)
    >>> parse_cparam_type("REAL[8]")
    ('REAL', '8', True)
    >>> parse_cparam_type("char[100]")
    ('char', '100', True)
    >>> parse_cparam_type("  REAL[ 16 ] ")
    ('REAL', '16', True)
    >>> parse_cparam_type("char[NAME]")  # doctest: +IGNORE_EXCEPTION_DETAIL
    Traceback (most recent call last):
    ...
    ValueError: Invalid array size 'NAME'
    """
    if "[" not in cparam_type or "]" not in cparam_type:
        return cparam_type.strip(), None, False

    base, after = cparam_type.split("[", 1)
    size_str, _ = after.split("]", 1)
    size = size_str.strip()

    if not size.isdigit():
        raise ValueError(f"Invalid array size '{size}'")

    return base.strip(), size, True


def register_BHaH_defines(module: str, bhah_defines_str: str) -> None:
    """
    Register C-code definitions for a given module.

    :param module: The name of the module registering the definitions.
    :param bhah_defines_str: The string containing the C-code definitions.

    Doctests:
    >>> par.glb_extras_dict = {}
    >>> register_BHaH_defines('test_module', '#define TEST_MACRO 1')
    >>> 'test_module' in par.glb_extras_dict['BHaH_defines']
    True
    """
    bhah_defines_dict = par.glb_extras_dict.setdefault("BHaH_defines", {})
    bhah_defines_dict[module] = bhah_defines_dict.get(module, "") + bhah_defines_str


def _register_general_defines(
    additional_includes: Optional[List[str]],
    enable_intrinsics: bool,
    intrinsics_header_lst: List[str],
    double_means: str,
) -> None:
    """
    Register general-purpose headers and macros.

    :param additional_includes: List of additional user-specified include files.
    :param enable_intrinsics: Flag to enable hardware intrinsics.
    :param intrinsics_header_lst: List of intrinsics header files.
    :param double_means: C type for the DOUBLE macro.
    """
    general_defines_str = """#include <ctype.h>   // Character type functions, such as isdigit, isalpha, etc.
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
        general_defines_str += "// output_BHaH_defines_h(...,enable_intrinsics=True) was called so we intrinsics headers:\n"
        for intr_header in intrinsics_header_lst:
            general_defines_str += f'#include "intrinsics/{intr_header}"\n'
    if additional_includes:
        for include in additional_includes:
            general_defines_str += f'#include "{include}"\n'

    real_means = par.parval_from_str("fp_type")
    general_defines_str += f"""#define REAL {real_means}
#define DOUBLE {double_means}

// These macros for MIN(), MAX(), and SQR() ensure that if the arguments inside
//   are a function/complex expression, the function/expression is evaluated
//   *only once* per argument. See https://lwn.net/Articles/983965/ for details.
// They are improvements over the original implementations:
// #define MIN(A, B) ( ((A) < (B)) ? (A) : (B) )
// #define MAX(A, B) ( ((A) > (B)) ? (A) : (B) )
// #define SQR(A) ((A) * (A))
#define MIN(A, B) ({{ \\
    __typeof__(A) _a = (A); \\
    __typeof__(B) _b = (B); \\
    _a < _b ? _a : _b; \\
}})
#define MAX(A, B) ({{ \\
    __typeof__(A) _a = (A); \\
    __typeof__(B) _b = (B); \\
    _a > _b ? _a : _b; \\
}})
#define SQR(A) ({{ \\
    __typeof__(A) _a = (A); \\
    _a * _a; \\
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
    has_define_type_params = False
    for cp_name, code_param in par.glb_code_params_dict.items():
        if code_param.cparam_type == "#define":
            if not has_define_type_params:
                has_define_type_params = True
                general_defines_str += "// START: CodeParameters declared as #define.\n"
            general_defines_str += f"""#ifndef {cp_name}
#define {cp_name} {code_param.defaultvalue} // {code_param.module}
#endif
"""
    if has_define_type_params:
        general_defines_str += "// END: CodeParameters declared as #define.\n"
    register_BHaH_defines("general", general_defines_str)


def _register_param_structs() -> None:
    """Register `params_struct` and `commondata_struct` definitions."""

    def format_c_declaration(
        c_type: str, var_name: str, module: str, description: str
    ) -> str:
        """
        Format a C declaration for a struct member.

        :param c_type: The C type of the variable (e.g., "REAL" or "int[8]").
        :param var_name: The name of the variable.
        :param module: The module defining the variable.
        :param description: The description of the variable.
        :return: A formatted C declaration string with comments.
        """
        base, size, is_array = parse_cparam_type(c_type)
        decl = (
            f"  char {var_name}[{size}];"
            if is_array and base.startswith("char")
            else (
                f"  {base} {var_name}[{size}];" if is_array else f"  {base} {var_name};"
            )
        )
        comment = f" // {description} ({module})" if description else f" // ({module})"
        return f"{decl}{comment}\n"

    params_struct_lines: List[str] = []
    commondata_struct_lines: List[str] = []

    for cp_name, code_param in par.glb_code_params_dict.items():
        if code_param.cparam_type != "#define":
            c_declaration = format_c_declaration(
                code_param.cparam_type,
                cp_name,
                code_param.module,
                code_param.description.strip(),
            )
            if code_param.commondata:
                commondata_struct_lines.append(c_declaration)
            else:
                params_struct_lines.append(c_declaration)

    if "commondata_struct" in par.glb_extras_dict:
        for module, item_list in par.glb_extras_dict["commondata_struct"].items():
            for item in item_list:
                c_code_line = f"  {item.c_declaration};"
                if item.description:
                    c_code_line += f" // <- {module}: {item.description}"
                commondata_struct_lines.append(f"{c_code_line}\n")

    params_struct_str = "typedef struct __params_struct__ {\n"
    params_struct_str += "".join(sorted(params_struct_lines))
    params_struct_str += "} params_struct;\n"

    commondata_struct_str = "typedef struct __commondata_struct__ {\n"
    commondata_struct_str += "".join(sorted(commondata_struct_lines))
    commondata_struct_str += "} commondata_struct;\n"

    register_BHaH_defines("params_struct", params_struct_str)
    register_BHaH_defines("commondata_struct", commondata_struct_str)


def _register_finite_difference_defines(
    add_one_for_upwinding: bool, enable_intrinsics: bool, define_no_simd_upwind: bool
) -> None:
    """
    Register finite difference C-code definitions.

    :param add_one_for_upwinding: Flag to add an extra ghost zone for upwinding.
    :param enable_intrinsics: Flag to enable hardware intrinsics.
    :param define_no_simd_upwind: Flag to define a non-SIMD UPWIND_ALG macro.
    """
    if any("finite_difference" in key for key in sys.modules):
        nghosts = int(par.parval_from_str("finite_difference::fd_order") / 2)
        if add_one_for_upwinding:
            nghosts += 1
        fd_defines_str = f"""
// Set the number of ghost zones
// Note that upwinding in e.g., BSSN requires that NGHOSTS = fd_order/2 + 1 <- Notice the +1.
#define NGHOSTS {nghosts}

// Declare NO_INLINE macro, used in FD functions. GCC v10+ compilations hang on complex RHS expressions (like BSSN) without this.
#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
    #define NO_INLINE __attribute__((noinline))
#elif defined(_MSC_VER)
    #define NO_INLINE __declspec(noinline)
#else
    #define NO_INLINE // Fallback for unknown compilers
#endif
"""
        if not enable_intrinsics and define_no_simd_upwind:
            fd_defines_str += """
// When enable_intrinsics = False, this is the UPWIND_ALG() macro:
#define UPWIND_ALG(UpwindVecU) UpwindVecU > 0.0 ? 1.0 : 0.0
"""
        register_BHaH_defines("finite_difference", fd_defines_str)


def _register_grid_defines(enable_rfm_precompute: bool) -> None:
    """
    Register grid-related C-code definitions.

    :param enable_rfm_precompute: Flag for reference metric precomputation.
    """
    grid_defines_str = gri.BHaHGridFunction.gridfunction_defines()
    grid_defines_str += r"""
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
// IDX4ptP: Similar to IDX4pt, but retrieves grid dimensions from the provided parameter structure
//        "params" instead of using global variables.
#define IDX4Ppt(params, gf, idx) \
  ((idx) + ((params)->Nxx_plus_2NGHOSTS0 * (params)->Nxx_plus_2NGHOSTS1 * (params)->Nxx_plus_2NGHOSTS2) * (gf))
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
    if not griddata_struct_def.startswith("\n"):
        griddata_struct_def = "\n" + griddata_struct_def
    grid_defines_str += griddata_struct_def

    register_BHaH_defines("grid", grid_defines_str)


def output_BHaH_defines_h(
    project_dir: str,
    additional_includes: Optional[List[str]] = None,
    restrict_pointer_type: str = "*restrict",
    enable_intrinsics: bool = True,
    intrinsics_header_lst: Optional[List[str]] = None,
    define_no_simd_UPWIND_ALG: bool = True,
    enable_rfm_precompute: bool = True,
    fin_NGHOSTS_add_one_for_upwinding_or_KO: bool = False,
    DOUBLE_means: str = "double",
    supplemental_defines_dict: Optional[Dict[str, str]] = None,
) -> None:
    r"""
    Output C code header file with macro definitions and other configurations.

    :param project_dir: Directory where the project C code is output.
    :param additional_includes: Additional header files to include.
    :param restrict_pointer_type: Pointer type for restrict, e.g., "*restrict".
    :param enable_intrinsics: Flag to enable hardware intrinsics.
    :param intrinsics_header_lst: List of intrinsics header files.
    :param define_no_simd_UPWIND_ALG: Flag to define a non-SIMD UPWIND_ALG.
    :param enable_rfm_precompute: Flag for reference metric precomputation.
    :param fin_NGHOSTS_add_one_for_upwinding_or_KO: Add ghost zone for upwinding.
    :param DOUBLE_means: C type for the DOUBLE macro.
    :param supplemental_defines_dict: Additional key-value pairs for defines.

    Doctests:
    >>> from nrpy.infrastructures.BHaH import MoLtimestepping
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import validate_strings
    >>> MoLtimestepping.register_all.register_CFunctions(register_MoL_step_forward_in_time=False)
    >>> project_path = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_path.mkdir(parents=True, exist_ok=True)
    >>> output_BHaH_defines_h(project_dir=str(project_path))
    >>> validate_strings((project_path / "BHaH_defines.h").read_text(), "BHaH_defines")
    """
    if intrinsics_header_lst is None:
        intrinsics_header_lst = ["simd_intrinsics.h"]

    project_path = Path(project_dir)
    project_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Register all defines by calling helper functions.
    _register_general_defines(
        additional_includes, enable_intrinsics, intrinsics_header_lst, DOUBLE_means
    )
    _register_param_structs()
    _register_finite_difference_defines(
        fin_NGHOSTS_add_one_for_upwinding_or_KO,
        enable_intrinsics,
        define_no_simd_UPWIND_ALG,
    )
    _register_grid_defines(enable_rfm_precompute)

    # Step 2: Assemble the final header file content.
    def output_key(key_name: str, item_name: str) -> str:
        """
        Format a C-code block for a given module.

        :param key_name: The name of the module.
        :param item_name: The C-code content for the module.
        :return: A formatted string containing the module's definitions.
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

    # for item in par.glb_extras_dict["BHaH_defines"]:
    #     print(item)
    registered_defines = par.glb_extras_dict.get("BHaH_defines", {})
    for core_module in core_modules_list:
        if core_module in registered_defines:
            # print(core_module)
            file_output_str += output_key(core_module, registered_defines[core_module])

    for key, item in sorted(registered_defines.items()):
        if key not in core_modules_list:
            print(f"Outputting non-core modules key = {key} to BHaH_defines.h")
            file_output_str += output_key(key, item)

    if supplemental_defines_dict:
        for key, value in supplemental_defines_dict.items():
            file_output_str += output_key(key, value)

    file_output_str += """
#ifndef BHAH_TYPEOF
#if __cplusplus >= 2000707L
#define BHAH_TYPEOF(a) decltype(a)
#elif defined(__GNUC__) || defined(__clang__) || defined(__NVCC__)
#define BHAH_TYPEOF(a) __typeof__(a)
#else
#define BHAH_TYPEOF(a)
#endif // END check for GCC, Clang, or C++
#endif // END BHAH_TYPEOF

#define BHAH_MALLOC(a, sz) \
do { \
    a = (BHAH_TYPEOF(a)) malloc(sz); \
} while(0);
#define BHAH_MALLOC__PtrMember(a, b, sz) \
do { \
    if (a) { \
        BHAH_MALLOC(a->b, sz); \
    } \
} while(0);

#define BHAH_FREE(a) \
do { \
    if (a) { \
        free((void*)(a)); \
        (a) = NULL; \
    } \
} while (0);
#define BHAH_FREE__PtrMember(a, b) \
do { \
    if (a) { \
        BHAH_FREE(a->b); \
    } \
} while(0);
"""
    parallelization = par.parval_from_str("parallelization")

    if parallelization != "openmp":
        malloc_func = parallel_utils.get_memory_malloc_function(parallelization)
        check_err_malloc = parallel_utils.get_check_errors_str(
            parallelization, malloc_func, opt_msg='Malloc: "#a" failed'
        )
        free_func = parallel_utils.get_memory_free_function(parallelization)
        check_err_free = parallel_utils.get_check_errors_str(
            parallelization, free_func, opt_msg='Free: "#a" failed'
        )
        file_output_str += rf"""
#define BHAH_MALLOC_DEVICE(a, sz) \
do {{ \
    {malloc_func}(&a, sz); \
    {check_err_malloc} \
}} while(0);
#define BHAH_FREE_DEVICE(a) \
do {{ \
    if (a) {{ \
        {free_func}((void*)(a)); \
        {check_err_free} \
        (a) = nullptr; \
    }} \
}} while (0);
"""
    if parallelization == "cuda":
        check_err_malloc_host = parallel_utils.get_check_errors_str(
            parallelization, "cudaMallocHost", opt_msg='Malloc: "#a" failed'
        )
        check_err_free_host = parallel_utils.get_check_errors_str(
            parallelization, "cudaFreeHost", opt_msg='Free: "#a" failed'
        )
        file_output_str += rf"""
#define BHAH_MALLOC_PINNED(a, sz) \
do {{ \
    cudaMallocHost((void**)&a, sz); \
    {check_err_malloc_host} \
}} while (0);

#define BHAH_FREE_PINNED(a) \
do {{ \
    if (a) {{ \
        cudaFreeHost((void*)(a)); \
        {check_err_free_host} \
        (a) = nullptr; \
    }} \
}} while (0);
#define BHAH_FREE_DEVICE__PtrMember(a, b) \
do {{ \
    if (a) {{ \
        decltype(a->b) tmp_ptr_##b = nullptr; \
        cudaMemcpy(&tmp_ptr_##b, &a->b, sizeof(void *), cudaMemcpyDeviceToHost); \
        if(tmp_ptr_##b) {{ \
            BHAH_FREE_DEVICE(tmp_ptr_##b); \
            cudaMemcpy(&a->b, &tmp_ptr_##b, sizeof(void *), cudaMemcpyHostToDevice); \
        }}\
    }} \
}} while(0);
#define BHAH_MALLOC_DEVICE__PtrMember(a, b, sz) \
do {{ \
    if (a) {{ \
        decltype(a->b) tmp_ptr_##b = nullptr; \
        BHAH_MALLOC_DEVICE(tmp_ptr_##b, sz); \
        cudaMemcpy(&a->b, &tmp_ptr_##b, sizeof(void *), cudaMemcpyHostToDevice); \
    }} \
}} while(0);
"""
    if parallelization in ["cuda"]:
        sync_func = parallel_utils.get_device_sync_function(parallelization)
        file_output_str += f"#define BHAH_DEVICE_SYNC() {sync_func}\n"
    file_output_str += r"#endif"

    file_output_str = file_output_str.replace("*restrict", restrict_pointer_type)

    bhah_defines_file = project_path / "BHaH_defines.h"
    with bhah_defines_file.open("w", encoding="utf-8") as file:
        file.write(clang_format(file_output_str))


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

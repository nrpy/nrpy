"""
Register wrapper functions for routines that span multiple coordinate systems.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import List, Dict
import hashlib
import re

import nrpy.c_function as cfc
from nrpy.infrastructures.BHaH.BHaH_defines_h import register_BHaH_defines


def get_CoordSystem_hash(CoordSystem: str) -> int:
    """
    Compute a 32-bit integer hash for the given coordinate system string using MD5.

    :param CoordSystem: The string representation of a coordinate system.
    :return: A 32-bit integer hash.

    >>> get_CoordSystem_hash("Cartesian")
    1445542812
    >>> get_CoordSystem_hash("SinhSpherical")
    579707936
    >>> get_CoordSystem_hash("Spherical")
    1307621703
    """
    # Convert CoordSystem to an md5 sum, then take the last 32 bits.
    # To ensure the result is a positive 32-bit signed integer, modulo by 2**31.
    return int(hashlib.md5(CoordSystem.encode()).hexdigest(), 16) % (2**31)


def register_CFunctions_CoordSystem_wrapper_funcs() -> None:
    """
    Register wrapper functions for C functions with respect to different coordinate systems.

    :raises ValueError: If a CFunction's return type is not void.
    """
    wrapper_func_list: List[str] = []
    base_CFunc_list: List[cfc.CFunction] = []
    for name, CFunction in cfc.CFunction_dict.items():
        CoordSystem = CFunction.CoordSystem_for_wrapper_func
        if CoordSystem:
            wrapper_func_name = CFunction.name.replace(f"__rfm__{CoordSystem}", "")
            if wrapper_func_name not in wrapper_func_list:
                wrapper_func_list += [wrapper_func_name]
                base_CFunc_list += [CFunction]
                if CFunction.c_type != "void":
                    raise ValueError(
                        f"Error creating wrapper function for CFunction {CFunction.name}: the return type {CFunction.c_type} is not void."
                    )

    CoordSystem_hash_dict: Dict[str, int] = {}
    for i, wrapper_func_name in enumerate(wrapper_func_list):
        base_CFunc = base_CFunc_list[i]
        list_of_CoordSystems: List[str] = []
        wrapper_subdir = ""
        for name, CFunc in cfc.CFunction_dict.items():
            if wrapper_func_name == name.replace(
                f"__rfm__{CFunc.CoordSystem_for_wrapper_func}", ""
            ):
                CoordSystem = CFunc.CoordSystem_for_wrapper_func
                CoordSystem_hash_dict[CoordSystem] = get_CoordSystem_hash(CoordSystem)
                list_of_CoordSystems += [CoordSystem]
                wrapper_subdir = CFunc.subdirectory.replace(CoordSystem, "")

        # Construct function call:
        # Remove linebreaks from params
        params = base_CFunc.params.replace("\n", "")
        # Remove multiple spaces to tidy up the line a little.
        params = " ".join(params.split())
        # Case: functions passed to function.
        # Remove all function arguments. E.g., void f(REAL x) -> void f.
        params = re.sub("\\(.*\\)", "", params)
        # Case: fixed-size arrays.
        # Remove size. E.g., const REAL xx[3] -> const REAL xx.
        params = re.sub("\\[[0-9A-Z_]*\\]", "", params)
        # Remove pointers, which might brush against variable name.
        params_sanitized = params.replace("*", " ")
        params = ""
        # Split sanitized params string according to commas:
        for param in params_sanitized.split(","):
            # Remove type, keep only what's after the comma.
            params += param.strip().split(" ")[-1] + ", "

        wrapper_body = "switch (params->CoordSystem_hash) {\n"
        for CoordSystem in sorted(list_of_CoordSystems):
            wrapper_body += f"case {CoordSystem.upper()}:\n"
            wrapper_body += (
                f"  {wrapper_func_name}__rfm__{CoordSystem}({params[:-2]});\n"
            )
            wrapper_body += "  break;\n"
        wrapper_body += rf"""default:
  fprintf(stderr, "ERROR in {wrapper_func_name}(): CoordSystem hash = %d not #define'd!\n", params->CoordSystem_hash);
  exit(1);
}}"""
        cfc.register_CFunction(
            subdirectory=wrapper_subdir,
            enable_simd=False,
            includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
            prefunc="",
            desc=base_CFunc.desc,
            c_type=base_CFunc.c_type,
            name=wrapper_func_name,
            CoordSystem_for_wrapper_func="",
            params=base_CFunc.params,
            include_CodeParameters_h=False,
            body=wrapper_body,
            clang_format_options=base_CFunc.clang_format_options,
        )

    BHd_str = ""
    for key, item in CoordSystem_hash_dict.items():
        BHd_str += f"#define {key.upper()} {item}\n"
    register_BHaH_defines(__name__, BHd_str)


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

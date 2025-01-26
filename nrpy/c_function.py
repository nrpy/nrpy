"""
Provide classes and functions for managing and registering C functions.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
"""

import os
from typing import Dict, List, Optional, Tuple

import nrpy.params as par
from nrpy.helpers.generic import clang_format, prefix_with_star


class CFunction:
    r"""
    Stores a C function and associated properties.

    :param subdirectory: Path from the root source directory to this C function. Defaults to the current directory.
    :param enable_simd: Boolean to enable SIMD. Default is False.
    :param includes: A list of strings representing include files.
    :param prefunc: A string containing code above the core function declaration. Defaults to an empty string.
    :param desc: A description of the function.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param cfunc_type: The C type of the function (e.g., void, int). Default is "void".
    :param name: The name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param include_CodeParameters_h: Boolean to enable C parameters. Default is False.
    :param body: The body of the function.
    :param postfunc: A string containing code below the core function definition. Defaults to an empty string.
    :param CoordSystem_for_wrapper_func: (BHaH only) Coordinate system for the wrapper function. E.g., if set to Cartesian -> create subdirectory/name() wrapper function and subdirectory/Cartesian/name__rfm__Cartesian(). Defaults to an empty string.
    :param ET_thorn_name: (ET only) Thorn home for this function.
    :param ET_schedule_bins_entries: (ET only) List of (schedule bin, schedule entry) tuples for Einstein Toolkit schedule.ccl.
    :param ET_current_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *this thorn's* param.ccl.
    :param ET_other_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *other thorn's* param.ccl.
    :param clang_format_options: Options for the clang-format tool. Defaults to "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".

    DocTests:
    >>> func = CFunction(desc="just a test... testing 1,2,3", name="main", params="", body="return 0;")
    >>> print(func.full_function)
    /**
     * just a test... testing 1,2,3
     */
    void main() { return 0; } // END FUNCTION main
    <BLANKLINE>
    >>> print(func.function_prototype)
    void main();
    >>> func2 = CFunction()  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: Error in CFunction: 'name' attribute must be set....
    """

    def __init__(
        self,
        subdirectory: str = os.path.join("."),
        enable_simd: bool = False,
        includes: Optional[List[str]] = None,
        prefunc: str = "",
        desc: str = "",
        cfunc_decorators: str = "",
        cfunc_type: str = "void",
        name: str = "",
        params: str = "",
        include_CodeParameters_h: bool = False,
        body: str = "",
        postfunc: str = "",
        CoordSystem_for_wrapper_func: str = "",
        ET_thorn_name: str = "",
        ET_schedule_bins_entries: Optional[List[Tuple[str, str]]] = None,
        ET_current_thorn_CodeParams_used: Optional[List[str]] = None,
        ET_other_thorn_CodeParams_used: Optional[List[str]] = None,
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    ) -> None:
        for attribute in [(name, "name"), (desc, "desc"), (body, "body")]:
            if not attribute[0]:
                raise ValueError(
                    f"Error in CFunction: '{attribute[1]}' attribute must be set."
                )
        if includes and not isinstance(includes, list):
            raise ValueError("includes must take the form of a list")
        if CoordSystem_for_wrapper_func and "griddata_struct" in params:
            raise ValueError(
                "griddata_struct cannot be passed into a CoordSystem-specific function"
            )

        self.subdirectory = subdirectory
        self.enable_simd = enable_simd
        self.includes = includes
        self.prefunc = prefunc
        self.desc = desc
        self.cfunc_type = cfunc_type
        self.name = name
        self.params = params
        self.include_CodeParameters_h = include_CodeParameters_h
        self.body = body
        self.postfunc = postfunc
        self.CoordSystem_for_wrapper_func = CoordSystem_for_wrapper_func
        self.ET_thorn_name = ET_thorn_name
        self.ET_schedule_bins_entries = ET_schedule_bins_entries
        self.ET_current_thorn_CodeParams_used = ET_current_thorn_CodeParams_used
        self.ET_other_thorn_CodeParams_used = ET_other_thorn_CodeParams_used
        self.cfunc_decorators = (
            f"{cfunc_decorators} " if cfunc_decorators != "" else cfunc_decorators
        )
        self.clang_format_options = clang_format_options

        self.function_prototype, self.raw_function, self.full_function = (
            self.generate_full_function()
        )

    @staticmethod
    def subdirectory_depth(subdirectory: str) -> int:
        """
        Calculate the depth of a given subdirectory by counting the number of folders in the provided path.
        It handles leading "./", trailing slashes, and consecutive slashes.

        :param subdirectory: The subdirectory path as a string.
        :return: The depth of the subdirectory.

        Example:
            >>> CFunction.subdirectory_depth('./folder1/folder2/folder3/')
            3
            >>> CFunction.subdirectory_depth('./folder1/folder2/folder3')
            3
            >>> CFunction.subdirectory_depth('folder1/folder2/folder3')
            3
            >>> CFunction.subdirectory_depth('folder1//')
            1
            >>> CFunction.subdirectory_depth('')
            0
            >>> CFunction.subdirectory_depth('.')
            0
            >>> CFunction.subdirectory_depth('.//')
            0
            >>> CFunction.subdirectory_depth('./folder1//folder2///')
            2
        """
        subdirectory = os.path.normpath(subdirectory)

        if subdirectory == ".":
            return 0

        folders = [folder for folder in subdirectory.split(os.sep) if folder]

        return len(folders)

    def generate_full_function(self) -> Tuple[str, str, str]:
        """
        Construct a full C function from a class instance.

        This method combines various components of a C function, including includes,
        pre-function definitions, function description, function prototype, and body,
        into a single, formatted C function string. It also optionally applies clang-format
        to the generated C function string based on the instance's clang format options.

        :return: A tuple containing two strings: the raw C function string and the clang-formatted C function string.

        :raises TypeError: If any item in the `includes` list is not a string.
        """
        rel_path_to_root_directory = ""
        for _ in range(self.subdirectory_depth(self.subdirectory)):
            rel_path_to_root_directory = os.path.join(rel_path_to_root_directory, "..")

        include_Cparams_str = ""
        if self.include_CodeParameters_h:
            CodeParameters_file_name = (
                "set_CodeParameters-simd.h"
                if self.enable_simd or "simd_width" in self.body
                else "set_CodeParameters.h"
            )
            if par.parval_from_str("Infrastructure") == "BHaH":
                include_Cparams_str = f'#include "{os.path.join(rel_path_to_root_directory, CodeParameters_file_name)}"\n'
            else:
                include_Cparams_str = f'#include "{CodeParameters_file_name}"\n'

        complete_func = ""

        if self.includes:
            for inc in self.includes:
                if not isinstance(inc, str):
                    raise TypeError(
                        f"Error in Cfunction(name={self.name}): includes must be a list of strings. Found includes = {self.includes}"
                    )

                if "<" in inc:
                    complete_func += f"#include {inc}\n"
                else:
                    if par.parval_from_str("Infrastructure") == "BHaH":
                        if any(
                            x in inc
                            for x in [
                                "BHaH_defines.h",
                                "BHaH_function_prototypes.h",
                                "simd_intrinsics.h",
                            ]
                        ):
                            inc = os.path.join(rel_path_to_root_directory, inc)
                    complete_func += f'#include "{inc}"\n'

        if self.prefunc:
            complete_func += f"{self.prefunc}\n"

        if self.desc:
            complete_func += f"/**\n{prefix_with_star(self.desc)}\n*/\n"

        function_prototype = (
            f"{self.cfunc_decorators}{self.cfunc_type} {self.name}({self.params});"
        )
        complete_func += f"{function_prototype.replace(';', '')} {{\n{include_Cparams_str}{self.body}}} // END FUNCTION {self.name}\n"

        complete_func += f"{self.postfunc}\n"

        return (
            function_prototype,
            complete_func,
            clang_format(complete_func, clang_format_options=self.clang_format_options),
        )


CFunction_dict: Dict[str, CFunction] = {}


def function_name_and_subdir_with_CoordSystem(
    subdirectory: str, name: str, CoordSystem_for_wrapper_func: str
) -> Tuple[str, str]:
    """
    Append a CoordSystem_for_wrapper_func string with a specific format to the provided name.

    :param subdirectory: The subdirectory within which we place this function.
    :param name: The wrapper function name.
    :param CoordSystem_for_wrapper_func: The coordinate system subdirectory string.
    :return: The coordinate-specific subdirectory and function name.

    DocTests:
        >>> function_name_and_subdir_with_CoordSystem(os.path.join("."), "xx_to_Cart", "SinhSpherical")
        ('./SinhSpherical', 'xx_to_Cart__rfm__SinhSpherical')
    """
    if CoordSystem_for_wrapper_func:
        return (
            os.path.join(subdirectory, CoordSystem_for_wrapper_func),
            f"{name}__rfm__{CoordSystem_for_wrapper_func}",
        )
    return subdirectory, name


def register_CFunction(
    subdirectory: str = os.path.join("."),
    enable_simd: bool = False,
    includes: Optional[List[str]] = None,
    prefunc: str = "",
    desc: str = "",
    cfunc_decorators: str = "",
    cfunc_type: str = "void",
    name: str = "",
    params: str = "",
    include_CodeParameters_h: bool = False,
    body: str = "",
    postfunc: str = "",
    CoordSystem_for_wrapper_func: str = "",
    ET_thorn_name: str = "",
    ET_schedule_bins_entries: Optional[List[Tuple[str, str]]] = None,
    ET_current_thorn_CodeParams_used: Optional[List[str]] = None,
    ET_other_thorn_CodeParams_used: Optional[List[str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Add a C function to a dictionary called CFunction_dict, using the provided parameters.

    :param subdirectory: Path from the root source directory to this C function. Defaults to the current directory.
    :param enable_simd: Boolean to enable SIMD. Default is False.
    :param includes: A list of strings representing include files.
    :param prefunc: A string containing code above the core function declaration. Defaults to an empty string.
    :param desc: A description of the function.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param cfunc_type: The C/C++ type of the function (e.g., void, int). Default is "void".
    :param name: The name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param include_CodeParameters_h: Boolean to enable C parameters. Default is False.
    :param body: The body of the function.
    :param postfunc: A string containing code after the core function declaration. Defaults to an empty string.
    :param CoordSystem_for_wrapper_func: (BHaH only) Coordinate system for the wrapper function. E.g., if set to Cartesian -> create subdirectory/name() wrapper function and subdirectory/Cartesian/name__rfm__Cartesian(). Defaults to an empty string.
    :param ET_thorn_name: (ET only) Thorn home for this function.
    :param ET_schedule_bins_entries: (ET only) List of tuples for Einstein Toolkit schedule.
    :param ET_current_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *this thorn's* param.ccl.
    :param ET_other_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *other thorn's* param.ccl.
    :param clang_format_options: Options for the clang-format tool. Defaults to "-style={BasedOnStyle: LLVM, ColumnLimit: 150}".

    :raises ValueError: If the name is already registered in CFunction_dict.

    DocTests:
        >>> register_CFunction(name="test_func", desc="test", body="return;")
        >>> "test_func" in CFunction_dict
        True
    """
    actual_subdirectory, actual_name = function_name_and_subdir_with_CoordSystem(
        subdirectory, name, CoordSystem_for_wrapper_func
    )
    if actual_name in CFunction_dict:
        raise ValueError(f"Error: already registered {actual_name} in CFunction_dict.")
    CFunction_dict[actual_name] = CFunction(
        subdirectory=actual_subdirectory,
        enable_simd=enable_simd,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=actual_name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
        postfunc=postfunc,
        CoordSystem_for_wrapper_func=CoordSystem_for_wrapper_func,
        ET_thorn_name=ET_thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
        ET_current_thorn_CodeParams_used=ET_current_thorn_CodeParams_used,
        ET_other_thorn_CodeParams_used=ET_other_thorn_CodeParams_used,
        cfunc_decorators=cfunc_decorators,
        clang_format_options=clang_format_options,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

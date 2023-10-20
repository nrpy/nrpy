"""
C function management and registration classes/functions.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
"""
import os
from typing import Optional, List, Dict, Tuple

import nrpy.params as par
from nrpy.helpers.generic import prefix_with_star, clang_format


class CFunction:
    r"""
    Stores a C function and associated properties.

    :param subdirectory: Path from the root source directory to this C function. Defaults to the current directory.
    :param enable_simd: Boolean to enable SIMD. Default is False.
    :param includes: A list of strings representing include files.
    :param prefunc: A string representing a pre-function code. Defaults to an empty string.
    :param desc: A description of the function.
    :param c_type: The C type of the function (e.g., void, int). Default is "void".
    :param name: The name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param include_CodeParameters_h: Boolean to enable C parameters. Default is False.
    :param body: The body of the function.
    :param CoordSystem_for_wrapper_func: (BHaH only) Coordinate system for the wrapper function. E.g., if set to Cartesian -> create subdirectory/name() wrapper function and subdirectory/Cartesian/name__rfm__Cartesian(). Defaults to an empty string.
    :param ET_thorn_name: (ET only) Thorn home for this function.
    :param ET_schedule_bins_entries: (ET only) List of (schedule bin, schedule entry) tuples for Einstein Toolkit schedule.ccl.
    :param ET_current_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *this thorn's* param.ccl.
    :param ET_other_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *other thorn's* param.ccl.
    :param clang_format_options: Options for the clang-format tool. Defaults to "-style={BasedOnStyle: LLVM, ColumnLimit: 200}".

    DocTests:
    >>> func = CFunction(desc="just a test... testing 1,2,3", name="main", params="", body="return 0;")
    >>> print(func.full_function)
    /*
     * just a test... testing 1,2,3
     */
    void main() { return 0; }
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
        c_type: str = "void",
        name: str = "",
        params: str = "",
        include_CodeParameters_h: bool = False,
        body: str = "",
        CoordSystem_for_wrapper_func: str = "",
        ET_thorn_name: str = "",
        ET_schedule_bins_entries: Optional[List[Tuple[str, str]]] = None,
        ET_current_thorn_CodeParams_used: Optional[List[str]] = None,
        ET_other_thorn_CodeParams_used: Optional[List[str]] = None,
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 200}",
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
        self.c_type = c_type
        self.name = name
        self.params = params
        self.include_CodeParameters_h = include_CodeParameters_h
        self.body = body
        self.CoordSystem_for_wrapper_func = CoordSystem_for_wrapper_func
        self.ET_thorn_name = ET_thorn_name
        self.ET_schedule_bins_entries = ET_schedule_bins_entries
        self.ET_current_thorn_CodeParams_used = ET_current_thorn_CodeParams_used
        self.ET_other_thorn_CodeParams_used = ET_other_thorn_CodeParams_used
        self.clang_format_options = clang_format_options

        self.function_prototype = f"{self.c_type} {self.name}({self.params});"
        self.raw_function, self.full_function = self.generate_full_function()

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
        # Remove the leading "./" if present
        if subdirectory.startswith("./"):
            subdirectory = subdirectory[2:]

        # Remove the trailing slash if present
        if subdirectory.endswith("/"):
            subdirectory = subdirectory[:-1]

        # If subdirectory is a single period, return 0
        if subdirectory == ".":
            return 0

        # Split by slashes and filter out any empty strings
        folders = [folder for folder in subdirectory.split("/") if folder]

        return len(folders)

    def generate_full_function(self) -> Tuple[str, str]:
        """Construct a full C function from a class instance."""
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
                        # BHaH-specific:
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
            complete_func += f"/*\n{prefix_with_star(self.desc)}\n*/\n"

        complete_func += f"{self.function_prototype.replace(';', '')} {{\n{include_Cparams_str}{self.body}}}\n"

        return complete_func, clang_format(
            complete_func, clang_format_options=self.clang_format_options
        )


# Contains a dictionary of CFunction objects
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
    c_type: str = "void",
    name: str = "",
    params: str = "",
    include_CodeParameters_h: bool = False,
    body: str = "",
    CoordSystem_for_wrapper_func: str = "",
    ET_thorn_name: str = "",
    ET_schedule_bins_entries: Optional[List[Tuple[str, str]]] = None,
    ET_current_thorn_CodeParams_used: Optional[List[str]] = None,
    ET_other_thorn_CodeParams_used: Optional[List[str]] = None,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 200}",
) -> None:
    """
    Add a C function to a dictionary called CFunction_dict, using the provided parameters.

    :param subdirectory: Path from the root source directory to this C function. Defaults to the current directory.
    :param enable_simd: Boolean to enable SIMD. Default is False.
    :param includes: A list of strings representing include files.
    :param prefunc: A string representing a pre-function code. Defaults to an empty string.
    :param desc: A description of the function.
    :param c_type: The C type of the function (e.g., void, int). Default is "void".
    :param name: The name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param include_CodeParameters_h: Boolean to enable C parameters. Default is False.
    :param body: The body of the function.
    :param CoordSystem_for_wrapper_func: (BHaH only) Coordinate system for the wrapper function. E.g., if set to Cartesian -> create subdirectory/name() wrapper function and subdirectory/Cartesian/name__rfm__Cartesian(). Defaults to an empty string.
    :param ET_thorn_name: (ET only) Thorn home for this function.
    :param ET_schedule_bins_entries: (ET only) List of tuples for Einstein Toolkit schedule.
    :param ET_current_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *this thorn's* param.ccl.
    :param ET_other_thorn_CodeParams_used: (ET only) List of CodeParameter names this function uses, for *other thorn's* param.ccl.
    :param clang_format_options: Options for the clang-format tool. Defaults to "-style={BasedOnStyle: LLVM, ColumnLimit: 200}".

    :raises ValueError: If the name is already registered in CFunction_dict.
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
        c_type=c_type,
        name=actual_name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
        CoordSystem_for_wrapper_func=CoordSystem_for_wrapper_func,
        ET_thorn_name=ET_thorn_name,
        ET_schedule_bins_entries=ET_schedule_bins_entries,
        ET_current_thorn_CodeParams_used=ET_current_thorn_CodeParams_used,
        ET_other_thorn_CodeParams_used=ET_other_thorn_CodeParams_used,
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

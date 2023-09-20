"""
C function management and registration classes/functions

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
"""

import os
from typing import Optional, List, Dict
from nrpy.helpers.generic import prefix_with_star, clang_format


class CFunction:
    r"""
    Represents a C function with all its properties.

    :param includes: A list of strings representing include files.
    :param prefunc: An optional string representing a pre-function code.
    :param desc: A description of the function.
    :param c_type: The C type of the function (e.g., void, int).
    :param name: The name of the function.
    :param params: A string representing the function's input parameters.
    :param preloop: An optional string representing a pre-loop code.
    :param body: A string representing the body of the function.
    :param postloop: An optional string representing a post-loop code.
    :param subdirectory: An optional string representing the path from the root source directory to this C function.
    :param include_CodeParameters_h: An optional boolean to enable C parameters.
    :param clang_format_options: Options for the clang-format tool.
    :param loop_options_kwargs: Any additional loop options passed as keyword arguments.

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
    ValueError: Error in CFunction: 'name', 'desc', and 'body' must be set. ...
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
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 200}",
    ) -> None:
        if not name or not desc or not body:
            raise ValueError(
                "Error in CFunction: 'name', 'desc', and 'body' must be set. "
                "Please provide appropriate values for these attributes."
            )
        if includes and not isinstance(includes, list):
            raise ValueError("includes must take the form of a list")

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
        self.clang_format_options = clang_format_options

        self.function_prototype = f"{self.c_type} {self.name}({self.params});"
        self.full_function = self.generate_full_function()

    @staticmethod
    def subdirectory_depth(subdirectory: str) -> int:
        """
        Calculates the depth of a given subdirectory by counting the number of folders in the provided path.
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

    def generate_full_function(self) -> str:
        """
        Constructs a full C function from a class instance.
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
            include_Cparams_str = f'#include "{os.path.join(rel_path_to_root_directory, CodeParameters_file_name)}"\n'

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

        return clang_format(
            complete_func, clang_format_options=self.clang_format_options
        )


# Contains a dictionary of CFunction objects
CFunction_dict: Dict[str, CFunction] = {}


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
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 200}",
) -> None:
    """
    Adds a C function to a dictionary called CFunction_dict, using the provided parameters.

    :param subdirectory: An optional string representing the path from the root source directory to this C function.
    :param enable_simd: An optional boolean to enable simd.
    :param includes: A list of strings representing include files.
    :param prefunc: An optional string representing a pre-function code.
    :param desc: A description of the function.
    :param c_type: The C type of the function (e.g., void, int).
    :param name: The name of the function.
    :param params: A string representing the function's input parameters.
    :param include_CodeParameters_h: An optional boolean to enable C parameters.
    :param body: A string representing the body of the function.
    :param clang_format_options: Options for the clang-format tool.

    :raises ValueError: If the name is already registered in CFunction_dict.
    """
    if name in CFunction_dict:
        raise ValueError(f"Error: already registered {name} in CFunction_dict.")
    CFunction_dict[name] = CFunction(
        subdirectory=subdirectory,
        enable_simd=enable_simd,
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        c_type=c_type,
        name=name,
        params=params,
        include_CodeParameters_h=include_CodeParameters_h,
        body=body,
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

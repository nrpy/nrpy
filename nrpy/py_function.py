"""
Provide classes and functions for managing and registering JAX-compatible Python functions.

Authors: Zachariah B. Etienne; zachetie **at** gmail **dot* com
         Ken Sible; ksible **at** outlook **dot* com
         Siddharth Mahesh; sidmahesh **at** gmail **dot* com
"""

import os
import textwrap
from typing import Dict, List, Optional


class PyFunction:
    r'''
    Stores a JAX-compatible Python function and associated properties.

    :param subdirectory: Path from the root source directory to this Python function. Defaults to the current directory.
    :param imports: A list of strings representing import files.
    :param prefunc: A string containing code above the core function declaration. Defaults to an empty string.
    :param desc: A *required* description of the function.
    :param pyfunc_decorators: Optional decorators for PyFunctions, e.g. JAX decorators, Equinox decorators, etc.
    :param name: The *required* name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param body: The *required* body of the function.
    :param postfunc: A string containing code below the core function definition. Defaults to an empty string.

    DocTests:
    >>> func = PyFunction(desc="just a test... testing 1,2,3", name="main", params="", body="return 0")
    >>> print(func.full_function)
    <BLANKLINE>
    def main():
        """
        just a test... testing 1,2,3
        """
        return 0
    <BLANKLINE>
    >>> func2 = PyFunction()  # doctest: +ELLIPSIS
    Traceback (most recent call last):
        ...
    ValueError: Error in PyFunction: 'name' attribute must be set.
    >>> imports = ["import jax","import jax.numpy as jnp"]
    >>> pyfunc_decorators = "@jax.jit"
    >>> name = "main"
    >>> desc = """Add two arrays
    ... # this function adds two arrays
    ... """
    >>> params = "a,b"
    >>> body = """res = jnp.add(a,b)
    ... return res"""
    >>> func = PyFunction(imports=imports,pyfunc_decorators=pyfunc_decorators,name=name,params=params,desc=desc,body=body)
    >>> print(func.full_function)
    import jax
    import jax.numpy as jnp
    <BLANKLINE>
    @jax.jit
    def main(a,b):
        """
        Add two arrays
        this function adds two arrays
        """
        res = jnp.add(a,b)
        return res
    <BLANKLINE>
    '''

    def __init__(
        self,
        subdirectory: str = os.path.join("."),
        imports: Optional[List[str]] = None,
        prefunc: str = "",
        desc: str = "",
        pyfunc_decorators: str = "",
        name: str = "",
        params: str = "",
        body: str = "",
        postfunc: str = "",
    ) -> None:
        for attribute in [(name, "name"), (desc, "desc"), (body, "body")]:
            if not attribute[0]:
                raise ValueError(
                    f"Error in PyFunction: '{attribute[1]}' attribute must be set."
                )
        if imports and not isinstance(imports, list):
            raise ValueError("Imports must take the form of a list")

        self.subdirectory = subdirectory
        self.imports = imports
        self.prefunc = prefunc
        self.desc = desc
        self.pyfunc_decorators = pyfunc_decorators
        self.name = name
        self.params = params
        self.body = body
        self.postfunc = postfunc

        self.full_function = self.generate_full_function()

    # Used within py_function to remove hashes from multi-line comments.
    @staticmethod
    def remove_hashes(input_string: str) -> str:
        r"""
        Wrap a multi-line string with \"\"\".

        This function standardizes comment blocks by ensuring the content
        is prefixed and suffixed with \"\"\". Any hashing (e.g., #) is removed.
        Trailing whitespace on each line is always removed.

        :param input_string: The input multi-line string to be reformatted.
        :return: The modified string with \"\"\" prefixed and suffixed.
        """
        # remove leading/trailing whitespace
        dedented = textwrap.dedent(input_string)
        # remove any hashes from the string
        hashes_removed = "\n".join(
            line.lstrip("#").strip() for line in dedented.splitlines()
        )
        # add triple quotes
        quotes_added = f'"""\n{hashes_removed}\n"""'
        # indent the string
        indented = textwrap.indent(quotes_added, "    ")
        return indented

    @staticmethod
    def indent_body(body: str) -> str:
        """
        Indent the body of a function with standard 4-space indentation.
        Preserves relative indentation between lines.

        :param body: The body of the function to be indented.
        :return: The indented body of the function.
        """
        body_as_lines = body.split("\n")
        out_body = []
        for body_line in body_as_lines:
            if body_line.strip() == "":
                out_body.append("\n")
                continue
            # Shift all non-blank lines right by four spaces, but do not
            # normalize or strip their existing indentation.
            out_body.append(f"    {body_line.rstrip()}\n")
        return "".join(out_body)

    @staticmethod
    def subdirectory_depth(subdirectory: str) -> int:
        """
        Calculate the depth of a given subdirectory by counting the number of folders in the provided path.
        It handles leading "./", trailing slashes, and consecutive slashes.

        :param subdirectory: The subdirectory path as a string.
        :return: The depth of the subdirectory.

        Example:
            >>> PyFunction.subdirectory_depth('./folder1/folder2/folder3/')
            3
            >>> PyFunction.subdirectory_depth('./folder1/folder2/folder3')
            3
            >>> PyFunction.subdirectory_depth('folder1/folder2/folder3')
            3
            >>> PyFunction.subdirectory_depth('folder1//')
            1
            >>> PyFunction.subdirectory_depth('')
            0
            >>> PyFunction.subdirectory_depth('.')
            0
            >>> PyFunction.subdirectory_depth('.//')
            0
            >>> PyFunction.subdirectory_depth('./folder1//folder2///')
            2
        """
        subdirectory = os.path.normpath(subdirectory)

        if subdirectory == ".":
            return 0

        folders = [folder for folder in subdirectory.split(os.sep) if folder]

        return len(folders)

    def generate_full_function(self) -> str:
        """
        Construct a full JAX-compatible Python function from a class instance.

        This method combines various components of a JAX-compatible Python function, including imports,
        pre-function definitions, function description, function prototype, and body,
        into a single, formatted JAX-compatible Python function string.

        :return: the raw Python function string.

        :raises TypeError: If any item in the `imports` list is not a string.
        """
        complete_func = ""

        if self.imports:
            for imp in self.imports:
                if not isinstance(imp, str):
                    raise TypeError(
                        f"Error in PyFunction(name={self.name}): imports must be a list of strings. Found imports = {self.imports}"
                    )

                complete_func += f"{imp}\n"

            # Add a blank line after imports.
            complete_func += "\n"

        if self.prefunc:
            # self.prefunc: Strip leading & trailing newlines, then add one newlines at start and two at the end.
            complete_func += "\n" + self.prefunc.strip("\n") + "\n"

        complete_func += f"{self.pyfunc_decorators}\ndef {self.name}({self.params}):\n"
        if self.desc:
            complete_func += f"{self.remove_hashes(self.desc)}\n"

        # self.body: Strip leading & trailing newlines, then add a single newline at the end of string. --v
        complete_func += f"{self.indent_body(self.body)}"
        # self.postfunc: Strip leading & trailing newlines, then add newlines at start and end.
        if self.postfunc != "":
            complete_func += "\n" + self.postfunc.strip("\n") + "\n"
        return complete_func


PyFunction_dict: Dict[str, PyFunction] = {}


def register_PyFunction(
    subdirectory: str = os.path.join("."),
    imports: Optional[List[str]] = None,
    prefunc: str = "",
    desc: str = "",
    pyfunc_decorators: str = "",
    name: str = "",
    params: str = "",
    body: str = "",
    postfunc: str = "",
) -> None:
    """
    Add a Python function to a dictionary called PyFunction_dict, using the provided parameters.

    :param subdirectory: Path from the root source directory to this Python function. Defaults to the current directory.
    :param imports: A list of strings representing import files.
    :param prefunc: A string containing code above the core function declaration. Defaults to an empty string.
    :param desc: A description of the function.
    :param pyfunc_decorators: Optional decorators for PyFunctions, e.g. CUDA identifiers, templates
    :param name: The name of the function.
    :param params: A string representing the function's input parameters. Defaults to an empty string.
    :param body: The body of the function.
    :param postfunc: A string containing code after the core function declaration. Defaults to an empty string.

    :raises ValueError: If the name is already registered in PyFunction_dict.

    DocTests:
        >>> register_PyFunction(name="test_func", desc="test", body="return;")
        >>> "test_func" in PyFunction_dict
        True
    """
    if name in PyFunction_dict:
        raise ValueError(f"Error: already registered {name} in PyFunction_dict.")
    PyFunction_dict[name] = PyFunction(
        subdirectory=subdirectory,
        imports=imports,
        prefunc=prefunc,
        desc=desc,
        pyfunc_decorators=pyfunc_decorators,
        name=name,
        params=params,
        body=body,
        postfunc=postfunc,
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

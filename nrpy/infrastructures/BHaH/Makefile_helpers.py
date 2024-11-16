"""
Constructs a Makefile for a BHaH C-code project based on CFunctions registered in the c_function NRPy module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import multiprocessing
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from nrpy.c_function import CFunction_dict
from nrpy.helpers.generic import clang_format


def output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir: str,
    project_name: str,
    exec_or_library_name: str = "",
    compiler_opt_option: str = "default",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "autodetect",
    create_lib: bool = False,
    static_lib: bool = False,
    lib_function_prefix: str = "",
    include_dirs: Optional[List[str]] = None,
    src_code_file_ext: str = "c",
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    """
    Output C functions registered to CFunction_dict and construct a Makefile for compiling C code.

    :param project_dir: The root directory of the C project.
    :param project_name: Name of the C project.
    :param exec_or_library_name: The name of the executable. If empty, same as project_name.
    :param compiler_opt_option: Compiler optimization option. Defaults to "default".
    :param addl_CFLAGS: Additional compiler flags.
    :param addl_libraries: Additional libraries to link.
    :param CC: C compiler to use. Defaults to "autodetect".
    :param create_lib: Whether to create a library. Defaults to False.
    :param static_lib: Whether the library created should be a static library. Defaults to False.
    :param lib_function_prefix: Prefix to add to library function names.
    :param include_dirs: List of include directories.
    :param src_code_file_ext: set what the file extension is for each code file.
    :param clang_format_options: Options for the clang-format tool.

    :raises ValueError: If the main() function is not defined in CFunction_dict.
    :raises FileNotFoundError: If the specified C compiler is not found.
    :raises TypeError: If addl_CFLAGS or include_dirs are not lists.
    :raises ValueError: If addl_CFLAGS or addl_libraries are specified incorrectly or OS is unsupported.

    DocTests:
        >>> from nrpy.c_function import register_CFunction
        >>> CFunction_dict.clear()
        >>> register_CFunction(
        ...     subdirectory="",
        ...     desc='Main function',
        ...     name='main',
        ...     cfunc_type='int',
        ...     params='int argc, char *argv[]',
        ...     body='return 0;',
        ... )
        >>> output_CFunctions_function_prototypes_and_construct_Makefile('/tmp/nrpy_BHaH_Makefile_doctest1', 'project_name')
    """
    if not create_lib and "main" not in CFunction_dict:
        raise ValueError(
            "Error: C codes will not compile if main() function not defined!\n"
            'Make sure that the main() function registered to CFunction_dict has name "main".'
        )
    if static_lib and not create_lib:
        raise ValueError(
            "`static_lib` set to True, but `create_lib` set to False. Both must be set to True to create a static library."
        )

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    if not exec_or_library_name:
        exec_or_library_name = project_name

    os_name = platform.system()
    if create_lib:
        if static_lib:
            ext = ".a"
        else:
            ext = ".so" if os_name == "Linux" else ".dylib"
        if not exec_or_library_name.endswith(ext):
            exec_or_library_name += ext

        # For shared libraries, we need -fPIC and -shared
        if not static_lib:
            if addl_CFLAGS is None:
                addl_CFLAGS = []
            for flag in ["-fPIC", "-shared"]:
                if flag not in addl_CFLAGS:
                    addl_CFLAGS.append(flag)
    else:
        if os_name not in ["Linux", "Darwin"]:
            raise ValueError(f"Unsupported operating system: {os_name}")

    Makefile_list_of_files: List[str] = []

    for name, CFunction in CFunction_dict.items():
        if lib_function_prefix:
            CFunction.name = f"{lib_function_prefix}{name}"
            CFunction.function_prototype = (
                f"{CFunction.cfunc_type} {CFunction.name}({CFunction.params});"
            )
            CFunction.function_prototype, _, CFunction.full_function = (
                CFunction.generate_full_function()
            )
            name = CFunction.name

        subdir_Path = Path(CFunction.subdirectory or ".")
        (project_Path / subdir_Path).mkdir(parents=True, exist_ok=True)
        c_file_path = project_Path / subdir_Path / f"{name}.{src_code_file_ext}"
        with open(c_file_path, "w", encoding="utf-8") as file:
            file.write(CFunction.full_function)
        Makefile_list_of_files.append(str(subdir_Path / f"{name}.{src_code_file_ext}"))

    # Output BHaH_function_prototypes.h
    with open(
        project_Path / "BHaH_function_prototypes.h", "w", encoding="utf-8"
    ) as file:
        outstr = "\n".join(
            CFunction_dict[key].function_prototype
            for key in sorted(CFunction_dict, key=str.lower)
        )
        file.write(clang_format(outstr, clang_format_options=clang_format_options))

    # Set compiler
    if CC == "autodetect":
        CC = "clang" if os_name == "Darwin" else "gcc"
    if shutil.which(CC) is None:
        raise FileNotFoundError(f"{CC} C compiler is not found")

    # Compiler flags
    CFLAGS_dict = {
        "default": "-std=gnu99 -O2 -march=native -g -Wall",
        "fast": "-std=gnu99 -O3 -funroll-loops -march=native -g -Wall",
        "debug": "-std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas",
        "nvcc": "-Xcompiler -fopenmp -Xcompiler -g -O2 -arch=native -O2 -Xcompiler=-march=native -Xcompiler -Wall --forward-unknown-to-host-compiler --Werror cross-execution-space-call --relocatable-device-code=true -allow-unsupported-compiler",
    }
    if addl_CFLAGS:
        if not isinstance(addl_CFLAGS, list):
            raise TypeError("addl_CFLAGS must be a list!")
        for key in CFLAGS_dict:
            CFLAGS_dict[key] += " " + " ".join(addl_CFLAGS)

    # Object files
    OBJ_FILES_str = "OBJ_FILES = " + " ".join(
        c_file.replace(f".{src_code_file_ext}", ".o")
        for c_file in sorted(Makefile_list_of_files, key=str.lower)
    )

    # Linker flags
    LDFLAGS_str = "LDFLAGS ="
    if addl_libraries:
        if not isinstance(addl_libraries, list):
            raise TypeError("addl_libraries must be a list!")
        LDFLAGS_str += " " + " ".join(addl_libraries)
    LDFLAGS_str += " -lm"

    # Include directories
    INCLUDEDIRS_str = "INCLUDEDIRS ="
    if include_dirs:
        if not isinstance(include_dirs, list):
            raise TypeError("include_dirs must be a list!")
        INCLUDEDIRS_str += " " + " ".join(f"-I{dir}" for dir in include_dirs)

    # Prepare the target rule
    if create_lib:
        if static_lib:
            # Static library
            target_rule = f"""{exec_or_library_name}: $(OBJ_FILES)
\tar rcs $@ $^

"""
        else:
            # Shared library
            target_rule = f"""{exec_or_library_name}: $(OBJ_FILES)
\t$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

"""
    else:
        # Executable
        target_rule = f"""{exec_or_library_name}: $(OBJ_FILES)
\t$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

"""

    # Construct the Makefile string
    Makefile_str = (
        rf"CC = {CC}  # Locally overwrites CC to {CC}\n"
        if CC == "nvcc"
        else f"CC ?= {CC}  # assigns the value CC to {CC} only if environment variable CC is not already set\n"
    )
    Makefile_str += f"""
CFLAGS = {CFLAGS_dict[compiler_opt_option]}
VALGRIND_CFLAGS = {CFLAGS_dict["debug"]}
{INCLUDEDIRS_str}
{LDFLAGS_str}
"""
    if not CC == "nvcc":
        Makefile_str += f"""
# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)
endif
"""
    else:
        Makefile_str += """
OPENMP_FLAG = -fopenmp
CFLAGS += $(OPENMP_FLAG)
LDFLAGS += $(OPENMP_FLAG)
"""
    Makefile_str += f"""
{OBJ_FILES_str}

all: {exec_or_library_name}

%.o: %.{src_code_file_ext} $(COMMON_HEADERS)
\t$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

{target_rule}
"""

    # Add the valgrind target
    Makefile_str += """valgrind: clean
\t$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
"""
    if not create_lib:
        Makefile_str += f"""\tvalgrind --track-origins=yes --leak-check=full -s ./{exec_or_library_name}

"""

    # Clean target
    Makefile_str += f"""# Use $(RM) to be cross-platform compatible.
clean:
\t$(RM) *.o */*.o *~ */*~ ./#* *.txt *.gp *.dat *.avi *.png {exec_or_library_name}
"""

    # Write the Makefile
    makefile_path = project_Path / "Makefile"
    with makefile_path.open("w", encoding="utf-8") as Makefile:
        Makefile.write(Makefile_str)


def compile_Makefile(
    project_dir: str,
    project_name: str,
    exec_or_library_name: str,
    compiler_opt_option: str = "fast",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "autodetect",
    attempt: int = 1,
) -> None:
    """
    Compile the Makefile using the given parameters and options.

    :param project_dir: Root directory for C code.
    :param project_name: Name of the project.
    :param exec_or_library_name: Name of the executable.
    :param compiler_opt_option: Compiler optimization option.
    :param addl_CFLAGS: Additional CFLAGS.
    :param addl_libraries: Additional libraries.
    :param CC: C compiler.
    :param attempt: Compilation attempt number.

    :raises FileNotFoundError: If the C compiler or make is not found.
    :raises RuntimeError: If compilation fails after two attempts.

    DocTests:
        >>> from nrpy.c_function import register_CFunction
        >>> CFunction_dict.clear()
        >>> register_CFunction(
        ...     name='main',
        ...     desc='Main function',
        ...     cfunc_type='int',
        ...     params='int argc, char *argv[]',
        ...     body='return 0;',
        ... )
        >>> compile_Makefile('/tmp/nrpy_BHaH_Makefile_doctest2/', 'project_name', 'executable')
    """
    if CC == "autodetect":
        os_name = platform.system()

        if os_name == "Darwin":
            CC = "clang"
        else:
            CC = "gcc"

    if shutil.which(CC) is None:
        raise FileNotFoundError(f"{CC} C compiler is not found")
    if shutil.which("make") is None:
        raise FileNotFoundError("make is not found")

    output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_or_library_name,
        compiler_opt_option=compiler_opt_option,
        addl_CFLAGS=addl_CFLAGS,
        addl_libraries=addl_libraries,
        CC=CC,
    )

    orig_working_directory = os.getcwd()
    os.chdir(project_dir)

    subprocess.run(
        f"make -j{multiprocessing.cpu_count() + 2}",
        stdout=subprocess.DEVNULL,
        shell=True,
        check=True,
    )

    os.chdir(orig_working_directory)

    exec_path = Path(project_dir).joinpath(exec_or_library_name)
    if not exec_path.is_file() and attempt == 1:
        # First clean up object files.
        for obj_file in Path(project_dir).glob("*.o"):
            obj_file.unlink()

        # Then retry compilation (recursion)
        compile_Makefile(
            project_dir=project_dir,
            project_name=project_name,
            exec_or_library_name=exec_or_library_name,
            compiler_opt_option="debug",
            addl_CFLAGS=addl_CFLAGS,
            addl_libraries=addl_libraries,
            CC=CC,
            attempt=2,
        )
    if not exec_path.is_file() and attempt == 2:
        raise RuntimeError("Compilation failed")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

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
from typing import Dict, List, Optional, Tuple

from nrpy.c_function import CFunction_dict
from nrpy.helpers.generic import clang_format

# Define constants for filenames to avoid magic strings
_BHAH_PROTOTYPES_H = "BHaH_function_prototypes.h"
_MAKEFILE = "Makefile"


def _autodetect_cc() -> str:
    """
    Autodetect the C compiler based on the operating system.

    :return: The detected C compiler ('clang' or 'gcc').
    """
    return "clang" if platform.system() == "Darwin" else "gcc"


def _validate_inputs(
    create_lib: bool,
    static_lib: bool,
    addl_CFLAGS: Optional[List[str]],
    include_dirs: Optional[List[str]],
) -> None:
    """
    Validate the inputs for Makefile generation.

    :param create_lib: Whether a library is being created.
    :param static_lib: Whether the library should be static.
    :param addl_CFLAGS: A list of additional compiler flags.
    :param include_dirs: A list of include directories.
    :raises ValueError: If main() is missing for an executable, for invalid library flags,
                        or if the OS is unsupported for executable creation.
    :raises TypeError: If addl_CFLAGS or include_dirs are not lists.
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
    if addl_CFLAGS and not isinstance(addl_CFLAGS, list):
        raise TypeError("addl_CFLAGS must be a list!")
    if include_dirs and not isinstance(include_dirs, list):
        raise TypeError("include_dirs must be a list!")
    if platform.system() not in ["Linux", "Darwin"] and not create_lib:
        raise ValueError(
            f"Unsupported operating system for executables: {platform.system()}"
        )


def _generate_c_files_and_header(
    project_path: Path, lib_function_prefix: str, src_code_file_ext: str
) -> List[str]:
    """
    Generate C source files and the main header file.

    :param project_path: The path to the project directory.
    :param lib_function_prefix: Prefix to add to library function names.
    :param src_code_file_ext: File extension for C source files.
    :return: A list of the paths to the generated C source files.
    """
    c_files: List[str] = []
    # Create C code files and directory structure
    for name, cfunc in CFunction_dict.items():
        if lib_function_prefix:
            # Update CFunction with the library prefix
            cfunc.name = f"{lib_function_prefix}{name}"
            (
                cfunc.function_prototype,
                _,
                cfunc.full_function,
            ) = cfunc.generate_full_function()
            name = cfunc.name

        subdir_path = project_path / (cfunc.subdirectory or ".")
        subdir_path.mkdir(parents=True, exist_ok=True)
        c_file_path = subdir_path / f"{name}.{src_code_file_ext}"
        with open(c_file_path, "w", encoding="utf-8") as file:
            file.write(cfunc.full_function)
        # Store the relative path for the Makefile
        c_files.append(
            str(Path(cfunc.subdirectory or ".") / f"{name}.{src_code_file_ext}")
        )

    # Generate BHaH_function_prototypes.h
    header_path = project_path / _BHAH_PROTOTYPES_H
    with open(header_path, "w", encoding="utf-8") as file:
        prototypes = "\n".join(
            CFunction_dict[key].function_prototype
            for key in sorted(CFunction_dict, key=str.lower)
        )
        file.write(clang_format(prototypes))

    return c_files


def _configure_compiler_and_flags(
    cc: str,
    create_lib: bool,
    static_lib: bool,
    addl_CFLAGS: Optional[List[str]],
) -> Tuple[str, Dict[str, str]]:
    """
    Determine the C compiler and configure compiler flags.

    :param cc: The C compiler to use, or "autodetect".
    :param create_lib: Whether a library is being created.
    :param static_lib: Whether the library should be static.
    :param addl_CFLAGS: A list of additional compiler flags.
    :return: A tuple containing the compiler and the CFLAGS dictionary.
    :raises FileNotFoundError: If the specified C compiler is not found.
    """
    if cc == "autodetect":
        cc = _autodetect_cc()
    if not shutil.which(cc):
        raise FileNotFoundError(f"{cc} C compiler is not found")

    updated_addl_CFLAGS = addl_CFLAGS.copy() if addl_CFLAGS else []
    if create_lib and not static_lib:
        # For shared libraries, we need -fPIC and -shared
        for flag in ["-fPIC", "-shared"]:
            if flag not in updated_addl_CFLAGS:
                updated_addl_CFLAGS.append(flag)

    cflags_dict = {
        "default": "-std=gnu99 -O2 -march=native -g -Wall -I.",
        "fast": "-std=gnu99 -O3 -funroll-loops -march=native -g -Wall -I.",
        "debug": "-I. -std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas",
        "nvcc": "-Xcompiler -fopenmp -Xcompiler -g -O2 -arch=native -O2 -Xcompiler=-march=native -Xcompiler -Wall --forward-unknown-to-host-compiler --Werror cross-execution-space-call --relocatable-device-code=true -allow-unsupported-compiler -I.",
    }

    if updated_addl_CFLAGS:
        cflags_str = " ".join(updated_addl_CFLAGS)
        for key in cflags_dict:
            cflags_dict[key] += f" {cflags_str}"

    return cc, cflags_dict


def _construct_makefile_content(
    cc: str,
    cflags: str,
    valgrind_cflags: str,
    include_dirs_str: str,
    ldflags_str: str,
    obj_files_str: str,
    exec_or_library_name: str,
    src_code_file_ext: str,
    target_rule: str,
    valgrind_rule: str,
    clean_rule: str,
    use_openmp: bool = True,
) -> str:
    """
    Construct the entire Makefile content using a template.

    :param cc: The C compiler executable name.
    :param cflags: The string of primary compiler flags.
    :param valgrind_cflags: The string of compiler flags for Valgrind builds.
    :param include_dirs_str: The Makefile string for include directories.
    :param ldflags_str: The Makefile string for linker flags.
    :param obj_files_str: The Makefile string for object files.
    :param exec_or_library_name: The name of the final target executable or library.
    :param src_code_file_ext: The file extension for source code files.
    :param target_rule: The Makefile rule for building the main target.
    :param valgrind_rule: The Makefile rule for the 'valgrind' target.
    :param clean_rule: The Makefile rule for the 'clean' target.
    :param use_openmp: If True, add OpenMP flags; if False, omit them.
    :return: The complete string content of the Makefile.
    """
    # Set the compiler assignment line based on whether it's nvcc
    cc_line = (
        f"CC = {cc}  # Locally overwrites CC to {cc}"
        if cc == "nvcc"
        else f"CC ?= {cc}  # assigns the value CC to {cc} only if environment variable CC is not already set"
    )

    if use_openmp:
        # Define the OpenMP detection block, which differs for nvcc
        if cc != "nvcc":
            openmp_block = """
    # Check for OpenMP support
    OPENMP_FLAG = -fopenmp
    COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

    ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
        CFLAGS += $(OPENMP_FLAG)
        LDFLAGS += $(OPENMP_FLAG)
    endif"""
        else:
            openmp_block = """
    OPENMP_FLAG = -fopenmp
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)"""
    else:
        openmp_block = ""  # No OpenMP flags

    # Assemble the final Makefile string using an f-string template for clarity
    return f"""{cc_line}

CFLAGS = {cflags}
VALGRIND_CFLAGS = {valgrind_cflags}
{include_dirs_str}
{ldflags_str}
{openmp_block.strip()}

{obj_files_str}

all: {exec_or_library_name}

%.o: %.{src_code_file_ext} $(COMMON_HEADERS)
\t$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

{target_rule}

{valgrind_rule}

{clean_rule}
"""


def output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir: str,
    project_name: str,
    exec_or_library_name: str = "",
    compiler_opt_option: str = "default",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_dirs_to_make: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "autodetect",
    create_lib: bool = False,
    static_lib: bool = False,
    lib_function_prefix: str = "",
    include_dirs: Optional[List[str]] = None,
    src_code_file_ext: str = "c",
    use_openmp: bool = True,
) -> None:
    """
    Output C functions registered to CFunction_dict and construct a Makefile for compiling C code.

    :param project_dir: The root directory of the C project.
    :param project_name: Name of the C project.
    :param exec_or_library_name: The name of the executable or library. If empty, same as project_name.
    :param compiler_opt_option: Compiler optimization option. Defaults to "default".
    :param addl_CFLAGS: Additional compiler flags.
    :param addl_dirs_to_make: Additional directories in which `make` should be run.
    :param addl_libraries: Additional libraries to link.
    :param CC: C compiler to use. Defaults to "autodetect".
    :param create_lib: Whether to create a library. Defaults to False.
    :param static_lib: Whether the library created should be a static library. Defaults to False.
    :param lib_function_prefix: Prefix to add to library function names.
    :param include_dirs: List of include directories.
    :param src_code_file_ext: Extension for C source files.
    :param use_openmp: If True, add OpenMP flags; if False, omit them.
    :raises TypeError: If 'addl_libraries' is not a list.

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
        >>> try:
        ...     # Test this function
        ...     output_CFunctions_function_prototypes_and_construct_Makefile('/tmp/nrpy_BHaH_Makefile_doctest1', 'project_name')
        ...     # Verify the content of the generated Makefile
        ...     with open('/tmp/nrpy_BHaH_Makefile_doctest1/Makefile', 'r') as f:
        ...         content = f.read()
        ...     assert 'CC ?= gcc' in content
        ...     assert 'OBJ_FILES = main.o' in content
        ...     assert 'all: project_name' in content
        ... finally:
        ...     # Clean up any created files
        ...     if Path('/tmp/nrpy_BHaH_Makefile_doctest1').exists():
        ...         shutil.rmtree('/tmp/nrpy_BHaH_Makefile_doctest1')
    """
    # Part 1: Validation and Initialization
    _validate_inputs(create_lib, static_lib, addl_CFLAGS, include_dirs)
    local_addl_dirs_to_make = addl_dirs_to_make or []

    final_exec_or_library_name = exec_or_library_name or project_name

    # Adjust exec_or_library_name for libraries
    if create_lib:
        os_name = platform.system()
        if static_lib:
            ext = ".a"
        else:
            ext = ".so" if os_name == "Linux" else ".dylib"
        if not final_exec_or_library_name.endswith(ext):
            final_exec_or_library_name += ext

    # Part 2: File Generation
    project_path = Path(project_dir)
    project_path.mkdir(parents=True, exist_ok=True)
    c_files = _generate_c_files_and_header(
        project_path, lib_function_prefix, src_code_file_ext
    )

    # Part 3: Compiler and Flag Configuration
    cc, cflags_dict = _configure_compiler_and_flags(
        CC, create_lib, static_lib, addl_CFLAGS
    )

    # Part 4: Construct Makefile Components
    # Object files string
    obj_files_str = "OBJ_FILES = " + " ".join(
        f.replace(f".{src_code_file_ext}", ".o") for f in sorted(c_files, key=str.lower)
    )

    # Linker flags string
    if addl_libraries and not isinstance(addl_libraries, list):
        raise TypeError("addl_libraries must be a list!")
    ldflags_list = addl_libraries or []
    ldflags_list.append("-lm")
    ldflags_str = "LDFLAGS = " + " ".join(ldflags_list)

    # Include directories string
    include_dirs_str = (
        "INCLUDEDIRS ="
        if not include_dirs
        else "INCLUDEDIRS = " + " ".join(f"-I{d}" for d in include_dirs)
    )

    # Main target rule
    if create_lib:
        if static_lib:
            target_rule = f"{final_exec_or_library_name}: $(OBJ_FILES)\n\tar rcs $@ $^"
        else:  # Shared library
            target_rule = f"{final_exec_or_library_name}: $(OBJ_FILES)\n\t$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)"
    else:  # Executable
        make_commands = "".join(f"\t$(MAKE) -C {d}\n" for d in local_addl_dirs_to_make)
        target_rule = f"{final_exec_or_library_name}: $(OBJ_FILES)\n{make_commands}\t$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)"

    # Valgrind target rule
    valgrind_make_cmds = "".join(
        f'\t$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" -C {d}\n'
        for d in local_addl_dirs_to_make
    )
    valgrind_rule = f'valgrind: clean\n{valgrind_make_cmds}\t$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all'
    if not create_lib:
        valgrind_rule += f"\n\tvalgrind --track-origins=yes --leak-check=full --show-leak-kinds=all -s ./{final_exec_or_library_name}"

    # Clean target rule
    clean_rule = f"""# Use $(RM) to be cross-platform compatible.
clean:
\t$(RM) *.o */*.o *~ */*~ */lib*.a ./#* *.txt *.gp *.dat *.avi *.png {final_exec_or_library_name}"""

    # Part 5: Assemble and Write Makefile
    makefile_content = _construct_makefile_content(
        cc=cc,
        cflags=cflags_dict[compiler_opt_option],
        valgrind_cflags=cflags_dict["debug"],
        include_dirs_str=include_dirs_str,
        ldflags_str=ldflags_str,
        obj_files_str=obj_files_str,
        exec_or_library_name=final_exec_or_library_name,
        src_code_file_ext=src_code_file_ext,
        target_rule=target_rule,
        valgrind_rule=valgrind_rule,
        clean_rule=clean_rule,
        use_openmp=use_openmp,
    )

    (project_path / _MAKEFILE).write_text(makefile_content, encoding="utf-8")


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
        >>> test_dir = '/tmp/nrpy_BHaH_Makefile_doctest2/'
        >>> try:
        ...     # Test compilation
        ...     compile_Makefile(test_dir, 'project_name', 'executable')
        ...     # Check if the executable exists
        ...     assert Path(test_dir).joinpath('executable').is_file()
        ... finally:
        ...     # Clean up
        ...     if Path(test_dir).exists():
        ...         shutil.rmtree(test_dir)
    """
    final_CC = CC
    if final_CC == "autodetect":
        final_CC = _autodetect_cc()

    if not shutil.which(final_CC):
        raise FileNotFoundError(f"{final_CC} C compiler is not found")
    if not shutil.which("make"):
        raise FileNotFoundError("make is not found")

    output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_or_library_name,
        compiler_opt_option=compiler_opt_option,
        addl_CFLAGS=addl_CFLAGS,
        addl_libraries=addl_libraries,
        CC=final_CC,
    )

    project_path = Path(project_dir)
    orig_working_directory = Path.cwd()
    os.chdir(project_path)

    try:
        subprocess.run(
            f"make -j{multiprocessing.cpu_count() + 2}",
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
            shell=True,
            check=True,
        )
    except subprocess.CalledProcessError as e:
        # If compilation fails, we might retry or just raise a more informative error
        if attempt == 2:
            raise RuntimeError(
                f"Compilation failed on second attempt. Stderr:\n{e.stderr.decode()}"
            ) from e
    finally:
        os.chdir(orig_working_directory)

    exec_path = project_path / exec_or_library_name
    if not exec_path.is_file():
        if attempt == 1:
            # First attempt failed. Clean up object files and retry with debug flags.
            for obj_file in project_path.glob("*.o"):
                obj_file.unlink()

            compile_Makefile(
                project_dir=project_dir,
                project_name=project_name,
                exec_or_library_name=exec_or_library_name,
                compiler_opt_option="debug",
                addl_CFLAGS=addl_CFLAGS,
                addl_libraries=addl_libraries,
                CC=final_CC,
                attempt=2,
            )
        else:
            # Second attempt also failed.
            raise RuntimeError("Compilation failed after two attempts.")


if __name__ == "__main__":
    import doctest
    import sys

    # Manually clear CFunction_dict and temp directories before running tests
    CFunction_dict.clear()
    for temp_dir in [
        "/tmp/nrpy_BHaH_Makefile_doctest1",
        "/tmp/nrpy_BHaH_Makefile_doctest2/",
    ]:
        if Path(temp_dir).exists():
            shutil.rmtree(temp_dir)

    results = doctest.testmod()

    # Clean up directories after tests
    for temp_dir in [
        "/tmp/nrpy_BHaH_Makefile_doctest1",
        "/tmp/nrpy_BHaH_Makefile_doctest2/",
    ]:
        if Path(temp_dir).exists():
            shutil.rmtree(temp_dir)

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

# nrpy/infrastructures/BHaH/Makefile_helpers.py
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
import warnings
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
) -> List[Tuple[str, List[str]]]:
    """
    Generate C source files and the main header file.

    :param project_path: The path to the project directory.
    :param lib_function_prefix: Prefix to add to library function names.
    :param src_code_file_ext: File extension for C source files.
    :return: Generated C source paths paired with their registered includes.
    """
    c_files_and_includes: List[Tuple[str, List[str]]] = []
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

        # Change the file extension from .cu to .cc if we're using SIMD intrinsics.
        this_src_code_file_ext = src_code_file_ext
        if src_code_file_ext == "cu" and (
            "simd_intrinsics.h" in cfunc.full_function
            or "interpolation_lagrange_uniform.h" in cfunc.full_function
        ):
            this_src_code_file_ext = "cc"

        subdir_path = project_path / (cfunc.subdirectory or ".")
        subdir_path.mkdir(parents=True, exist_ok=True)
        c_file_path = subdir_path / f"{name}.{this_src_code_file_ext}"
        with open(c_file_path, "w", encoding="utf-8") as file:
            file.write(cfunc.full_function)
        # Store the relative path for the Makefile
        c_files_and_includes.append(
            (
                str(
                    Path(cfunc.subdirectory or ".") / f"{name}.{this_src_code_file_ext}"
                ),
                list(cfunc.includes or []),
            )
        )

    # Generate BHaH_function_prototypes.h
    header_path = project_path / _BHAH_PROTOTYPES_H
    with open(header_path, "w", encoding="utf-8") as file:
        prototypes = "\n".join(
            CFunction_dict[key].function_prototype
            for key in sorted(CFunction_dict, key=str.lower)
        )
        file.write(clang_format(prototypes))

    return c_files_and_includes


def _configure_compiler_and_flags(
    cc: str,
    addl_CFLAGS: Optional[List[str]],
) -> Tuple[str, Dict[str, str], str]:
    """
    Determine the C compiler and configure compiler flags.

    :param cc: The C compiler to use, or "autodetect".
    :param addl_CFLAGS: A list of additional compiler flags.
    :return: The compiler, CFLAGS dictionary, and CXXFLAGS string.
    """
    if cc == "autodetect":
        cc = _autodetect_cc()
    if not shutil.which(cc):
        warnings.warn(f"{cc} C compiler is not found", UserWarning)

    updated_addl_CFLAGS = addl_CFLAGS.copy() if addl_CFLAGS else []

    cflags_dict = {
        "default": "-std=gnu99 -O2 -march=native -g -Wall",
        "fast": "-std=gnu99 -O3 -funroll-loops -march=native -g -Wall",
        "debug": "-std=gnu99 -O2 -g -Wall -Wno-unknown-pragmas",
        "nvcc": "-Xcompiler -g -O2 -arch=native -Xcompiler=-march=native -Xcompiler -Wall --forward-unknown-to-host-compiler --Werror cross-execution-space-call --relocatable-device-code=true -allow-unsupported-compiler",
    }

    if updated_addl_CFLAGS:
        cflags_str = " ".join(updated_addl_CFLAGS)
        for key in cflags_dict:
            cflags_dict[key] += f" {cflags_str}"

    cxxflags = "-O2 -g -Wall -Wno-unknown-pragmas -march=native"
    if cc != "nvcc" and updated_addl_CFLAGS:
        cxxflags += f" {' '.join(updated_addl_CFLAGS)}"

    return cc, cflags_dict, cxxflags


def _construct_makefile_content(
    cc: str,
    cflags: str,
    cxxflags: str,
    valgrind_cflags: str,
    cppflags: str,
    ldlibs: str,
    sources: List[str],
    direct_header_rules: List[str],
    exec_or_library_name: str,
    addl_dirs_to_make: List[str],
    create_lib: bool,
    static_lib: bool,
    use_openmp: bool = True,
) -> str:
    """
    Construct the entire Makefile content using a template.

    :param cc: The C compiler executable name.
    :param cflags: The string of primary compiler flags.
    :param cxxflags: The string of C++ compiler flags.
    :param valgrind_cflags: The string of compiler flags for Valgrind builds.
    :param cppflags: Project preprocessor flags.
    :param ldlibs: Libraries and opaque caller-provided link flags.
    :param sources: Registered generated source paths.
    :param direct_header_rules: Registry-derived direct-header dependency rules.
    :param exec_or_library_name: The name of the final target executable or library.
    :param addl_dirs_to_make: Additional generated projects to build.
    :param create_lib: Whether the target is a library.
    :param static_lib: Whether the target is a static library.
    :param use_openmp: If True, add OpenMP flags; if False, omit them.
    :return: The complete string content of the Makefile.
    """
    if cc == "nvcc":
        compiler_block = """# CUDA projects require nvcc unless CC is overridden on the command line.
CC := nvcc"""
    else:
        compiler_block = f"""# Replace GNU Make's built-in cc, while honoring environment and command-line overrides.
ifeq ($(origin CC),default)
CC := {cc}
endif"""

    openmp_default = "1" if use_openmp else "0"
    pic_cflags = "-fPIC" if create_lib and not static_lib else ""
    pic_nvccflags = "-Xcompiler -fPIC" if pic_cflags else ""
    if cc == "nvcc":
        openmp_block = f"""OPENMP ?= {openmp_default}

ifeq ($(OPENMP),1)
OMP_NVCCFLAGS := -Xcompiler -fopenmp
OMP_CXXFLAGS := -fopenmp
OMP_LDFLAGS := -Xcompiler -fopenmp
else
OMP_NVCCFLAGS :=
OMP_CXXFLAGS := -Wno-unknown-pragmas
OMP_LDFLAGS :=
endif

ALL_NVCCFLAGS = $(NVCCFLAGS) $(OMP_NVCCFLAGS) $(PIC_NVCCFLAGS)
ALL_CXXFLAGS = $(CXXFLAGS) $(OMP_CXXFLAGS) $(PIC_CFLAGS)
ALL_LDFLAGS = $(LDFLAGS) $(OMP_LDFLAGS)"""
        flags_block = f"""NVCCFLAGS = {cflags}
CXXFLAGS = {cxxflags}
VALGRIND_CFLAGS = {valgrind_cflags}
PIC_CFLAGS := {pic_cflags}
PIC_NVCCFLAGS := {pic_nvccflags}"""
    else:
        openmp_block = f"""OPENMP ?= {openmp_default}
OPENMP_FLAG := -fopenmp

ifeq ($(OPENMP),1)
OPENMP_SUPPORTED := $(shell printf '%s\\n' '#include <omp.h>' 'int main(void) {{ return omp_get_max_threads() < 1; }}' | $(CC) $(ALL_CPPFLAGS) $(OPENMP_FLAG) -x c - -o /dev/null >/dev/null 2>&1 && echo YES || echo NO)
endif

ifeq ($(OPENMP_SUPPORTED),YES)
OMP_CFLAGS := $(OPENMP_FLAG)
OMP_CXXFLAGS := $(OPENMP_FLAG)
OMP_LDFLAGS := $(OPENMP_FLAG)
else
OMP_CFLAGS := -Wno-unknown-pragmas
OMP_CXXFLAGS := -Wno-unknown-pragmas
OMP_LDFLAGS :=
endif

ALL_CFLAGS = $(CFLAGS) $(OMP_CFLAGS) $(PIC_CFLAGS)
ALL_CXXFLAGS = $(CXXFLAGS) $(OMP_CXXFLAGS) $(PIC_CFLAGS)
ALL_LDFLAGS = $(LDFLAGS) $(OMP_LDFLAGS)"""
        flags_block = f"""CFLAGS = {cflags}
CXXFLAGS = {cxxflags}
VALGRIND_CFLAGS = {valgrind_cflags}
PIC_CFLAGS := {pic_cflags}"""

    if sources:
        sources_block = "SOURCES := \\\n" + " \\\n".join(
            f"    {source}" for source in sources
        )
    else:
        sources_block = "SOURCES :="

    compile_rules_list: List[str] = []
    source_extensions = sorted({Path(source).suffix.lstrip(".") for source in sources})
    for source_extension in source_extensions:
        compiler_command = (
            "$(CXX) $(ALL_CPPFLAGS) $(ALL_CXXFLAGS)"
            if source_extension in ["cc", "cpp", "cxx"]
            else (
                "$(CC) $(ALL_CPPFLAGS) $(ALL_NVCCFLAGS)"
                if source_extension == "cu"
                else "$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS)"
            )
        )
        compile_rules_list.append(f"""%.o: %.{source_extension}
\t{compiler_command} $(DEPFLAGS) -c $< -o $@""")
    compile_rules = "\n\n".join(compile_rules_list)

    additional_projects_prerequisite = ""
    additional_projects_rule = ""
    object_order_only_rule = ""
    recursive_clean = ""
    phony_targets = ["all", "valgrind", "clean"]
    if addl_dirs_to_make:
        additional_projects_prerequisite = " additional-projects"
        phony_targets.append("additional-projects")
        make_commands = "\n".join(
            f"\t+$(MAKE) -C {directory}" for directory in addl_dirs_to_make
        )
        additional_projects_rule = f"""additional-projects:
{make_commands}
"""
        object_order_only_rule = "$(OBJECTS): | additional-projects\n"
        recursive_clean = "\n" + "\n".join(
            f"\t+$(MAKE) -C {directory} clean" for directory in addl_dirs_to_make
        )

    target_prerequisites = f"$(OBJECTS){additional_projects_prerequisite}"
    if create_lib and static_lib:
        target_recipe = "\t$(RM) $@\n\t$(AR) rcs $@ $(OBJECTS)"
    else:
        shared_flag = " -shared" if create_lib else ""
        nvcc_link_flags = " $(NVCCFLAGS)" if cc == "nvcc" else ""
        has_cxx_sources = bool({"cc", "cpp", "cxx"}.intersection(source_extensions))
        linker = "$(CC)" if cc == "nvcc" or not has_cxx_sources else "$(CXX)"
        target_recipe = (
            f"\t{linker}{nvcc_link_flags} $(ALL_LDFLAGS){shared_flag} "
            "-o $@ $(OBJECTS) $(ALL_LDLIBS)"
        )

    if cc == "nvcc":
        valgrind_rule = """valgrind:
\t@echo "Valgrind target is unavailable for CUDA builds."
\t@false"""
    else:
        run_valgrind = ""
        if not create_lib:
            run_valgrind = (
                f"\n\tOMP_NUM_THREADS=1 valgrind --track-origins=yes "
                "--leak-check=full --show-leak-kinds=all "
                f"--error-exitcode=1 -s ./{exec_or_library_name}"
            )
        valgrind_rule = f"""valgrind:
\t+$(MAKE) clean
\t+$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" OPENMP=0 all{run_valgrind}"""

    header_rules = "\n".join(direct_header_rules)
    if header_rules:
        header_rules += "\n"

    return f""".DEFAULT_GOAL := all
.DELETE_ON_ERROR:

{compiler_block}

PROJECT_CPPFLAGS := {cppflags}
ALL_CPPFLAGS = $(CPPFLAGS) $(PROJECT_CPPFLAGS)
{flags_block}
LDFLAGS ?=
PROJECT_LDLIBS := {ldlibs}
ALL_LDLIBS = $(LDLIBS) $(PROJECT_LDLIBS)
DEPFLAGS = -MMD -MP -MF $(@:.o=.d) -MT $@

{openmp_block}

{sources_block}
OBJECTS := $(addsuffix .o,$(basename $(SOURCES)))
DEPFILES := $(OBJECTS:.o=.d)

.PHONY: {" ".join(phony_targets)}

all: {exec_or_library_name}

{header_rules}{compile_rules}

$(OBJECTS): Makefile
{additional_projects_rule}{object_order_only_rule}{exec_or_library_name}: {target_prerequisites}
{target_recipe}

{valgrind_rule}

# Remove only files owned by this generated build.
clean:
\t$(RM) $(OBJECTS) $(DEPFILES) {exec_or_library_name}{recursive_clean}

-include $(DEPFILES)
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
    r"""
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

    Doctests:
        >>> from nrpy.c_function import register_CFunction
        >>> CFunction_dict.clear()
        >>> register_CFunction(
        ...     subdirectory="",
        ...     desc='Main function',
        ...     name='main',
        ...     cfunc_type='int',
        ...     params='int argc, char *argv[]',
        ...     includes=['project.h', '<stdio.h>'],
        ...     body='return 0;',
        ... )
        >>> try:
        ...     test_path = Path('/tmp/nrpy_BHaH_Makefile_doctest1')
        ...     test_path.mkdir(parents=True, exist_ok=True)
        ...     _ = test_path.joinpath('project.h').write_text('', encoding='utf-8')
        ...     # Test this function
        ...     output_CFunctions_function_prototypes_and_construct_Makefile(
        ...         '/tmp/nrpy_BHaH_Makefile_doctest1',
        ...         'project_name',
        ...         addl_dirs_to_make=['support'],
        ...     )
        ...     # Verify the content of the generated Makefile
        ...     with open('/tmp/nrpy_BHaH_Makefile_doctest1/Makefile', 'r') as f:
        ...         content = f.read()
        ...     assert 'ifeq ($(origin CC),default)' in content
        ...     assert 'CC := gcc' in content
        ...     assert 'SOURCES := \\\n    main.c' in content
        ...     assert 'OBJECTS := $(addsuffix .o,$(basename $(SOURCES)))' in content
        ...     assert 'DEPFILES := $(OBJECTS:.o=.d)' in content
        ...     assert 'main.o: project.h' in content
        ...     assert 'main.o: main.c' not in content
        ...     assert 'main.o: project.h <stdio.h>' not in content
        ...     assert 'DEPFLAGS = -MMD -MP -MF $(@:.o=.d) -MT $@' in content
        ...     assert 'PROJECT_CPPFLAGS := -I.' in content
        ...     assert 'ALL_CPPFLAGS = $(CPPFLAGS) $(PROJECT_CPPFLAGS)' in content
        ...     assert '$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) $(DEPFLAGS) -c $< -o $@' in content
        ...     assert '$(CC) $(ALL_LDFLAGS) -o $@ $(OBJECTS) $(ALL_LDLIBS)' in content
        ...     assert '-include $(DEPFILES)' in content
        ...     assert 'all: project_name' in content
        ...     assert '$(OBJECTS): Makefile' in content
        ...     assert '$(OBJECTS): | additional-projects' in content
        ...     assert 'project_name: $(OBJECTS) additional-projects' in content
        ...     assert '\t+$(MAKE) -C support' in content
        ...     makefile_lines = content.splitlines()
        ...     clean_recipe = makefile_lines[makefile_lines.index('clean:') + 1]
        ...     assert clean_recipe == '\t$(RM) $(OBJECTS) $(DEPFILES) project_name'
        ...     assert '\t+$(MAKE) -C support clean' in content
        ...     assert '*.out' not in content
        ...     cuda_path = test_path / 'cuda'
        ...     output_CFunctions_function_prototypes_and_construct_Makefile(
        ...         str(cuda_path),
        ...         'cuda_project',
        ...         CC='nvcc',
        ...         compiler_opt_option='nvcc',
        ...         src_code_file_ext='cu',
        ...     )
        ...     cuda_content = cuda_path.joinpath('Makefile').read_text(encoding='utf-8')
        ...     assert 'CC := nvcc' in cuda_content
        ...     assert '    main.cu' in cuda_content
        ...     assert '$(CC) $(ALL_CPPFLAGS) $(ALL_NVCCFLAGS) $(DEPFLAGS) -c $< -o $@' in cuda_content
        ...     assert '$(CC) $(NVCCFLAGS) $(ALL_LDFLAGS) -o $@ $(OBJECTS) $(ALL_LDLIBS)' in cuda_content
        ...     assert 'Valgrind target is unavailable for CUDA builds.' in cuda_content
        ... finally:
        ...     # Clean up any created files
        ...     if Path('/tmp/nrpy_BHaH_Makefile_doctest1').exists():
        ...         shutil.rmtree('/tmp/nrpy_BHaH_Makefile_doctest1')
    """
    # Step 1: Validate inputs and initialize local state
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

    # Step 2: Generate source files and create the output directory
    project_path = Path(project_dir)
    project_path.mkdir(parents=True, exist_ok=True)
    c_files_and_includes = _generate_c_files_and_header(
        project_path, lib_function_prefix, src_code_file_ext
    )

    # Step 3: Configure compiler and flags
    cc, cflags_dict, cxxflags = _configure_compiler_and_flags(CC, addl_CFLAGS)

    # Step 4: Construct Makefile components
    # Linker libraries and opaque caller-provided link flags
    if addl_libraries and not isinstance(addl_libraries, list):
        raise TypeError("addl_libraries must be a list!")
    ldlibs_list = list(addl_libraries or [])
    ldlibs_list.append("-lm")
    ldlibs = " ".join(ldlibs_list)

    cppflags = " ".join(
        ["-I."] + [f"-I{directory}" for directory in include_dirs or []]
    )

    # Explicit source, object, and project-local header dependency rules
    project_root = project_path.resolve()
    project_include_dirs = []
    for include_dir in include_dirs or []:
        include_dir_path = Path(include_dir)
        if not include_dir_path.is_absolute():
            include_dir_path = project_root / include_dir_path
        resolved_include_dir = include_dir_path.resolve()
        try:
            resolved_include_dir.relative_to(project_root)
        except ValueError:
            continue
        project_include_dirs.append(resolved_include_dir)

    direct_header_rules = []
    sorted_c_files_and_includes = sorted(
        c_files_and_includes, key=lambda item: item[0].lower()
    )
    for c_file, registered_includes in sorted_c_files_and_includes:
        c_file_path = Path(c_file)
        object_file = str(c_file_path.with_suffix(".o"))

        local_headers = []
        for registered_include in registered_includes:
            if registered_include.startswith("<") and registered_include.endswith(">"):
                include_name = registered_include[1:-1]
                include_search_dirs = [project_root] + project_include_dirs
            elif "<" in registered_include or ">" in registered_include:
                continue
            else:
                include_name = registered_include
                include_search_dirs = [
                    project_root / c_file_path.parent,
                    project_root,
                ] + project_include_dirs

            for include_search_dir in include_search_dirs:
                header_candidate = Path(
                    os.path.normpath(str(include_search_dir / include_name))
                )
                header_path = header_candidate.resolve()
                try:
                    relative_header_path = header_candidate.relative_to(project_root)
                    header_path.relative_to(project_root)
                except ValueError:
                    continue
                relative_header = str(relative_header_path)
                if header_path.is_file():
                    if relative_header not in local_headers:
                        local_headers.append(relative_header)
                    break

        if local_headers:
            direct_header_rules.append(f"{object_file}: {' '.join(local_headers)}")

    # Step 5: Assemble and write the Makefile
    makefile_content = _construct_makefile_content(
        cc=cc,
        cflags=cflags_dict[compiler_opt_option],
        cxxflags=cxxflags,
        valgrind_cflags=cflags_dict["debug"],
        cppflags=cppflags,
        ldlibs=ldlibs,
        sources=[source for source, _ in sorted_c_files_and_includes],
        direct_header_rules=direct_header_rules,
        exec_or_library_name=final_exec_or_library_name,
        addl_dirs_to_make=local_addl_dirs_to_make,
        create_lib=create_lib,
        static_lib=static_lib,
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
        warnings.warn(f"{final_CC} C compiler is not found", UserWarning)
    if not shutil.which("make"):
        warnings.warn("make is not found")

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

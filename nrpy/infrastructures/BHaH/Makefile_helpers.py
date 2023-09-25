"""
Construct Makefile for BHaH C-code project, based
 upon C functions registered through the c_function
 NRPy+ module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path
import os
import shutil
import subprocess
import multiprocessing
from typing import List, Optional
from nrpy.c_function import CFunction_dict


def output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir: str,
    project_name: str,
    exec_name: str = "",
    compiler_opt_option: str = "default",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "gcc",
    create_lib: bool = False,
    include_dirs: Optional[List[str]] = None,
) -> None:
    """
    Output C functions registered to CFunction_dict and construct a Makefile for compiling C code.

    :param project_dir: The root directory of the C project.
    :param project_name: Name of the C project.
    :param exec_name: The name of the executable. If set to empty string, same as project_name.
    :param compiler_opt_option: Compiler optimization option. Defaults to "default". Other options: "fast" and "debug"
    :param addl_CFLAGS: Additional compiler flags. Must be a list.
    :param addl_libraries: Additional libraries to link. Must be a list.
    :param CC: C compiler to use. Defaults to "gcc".
    :param create_lib: Whether to create a library. Defaults to False.
    :param include_dirs: List of include directories. Must be a list.

    :raises SystemExit: Exits if errors are encountered.
    """

    if not create_lib and "main" not in CFunction_dict:
        raise SystemExit(
            "output_CFunctions_function_prototypes_and_construct_Makefile() error: C codes will not compile if main() function not defined!\n"
            '    Make sure that the main() function registered to CFunction_dict has name "main".'
        )

    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    if exec_name == "":
        exec_name = project_name

    Makefile_list_of_files: List[str] = []

    def add_to_Makefile(project_Path: Path, path_and_file: str) -> Path:
        """
        Add the given path and file to the Makefile list of files and return the joined path.

        :param project_Path: The Path of the project
        :param path_and_file: The path and file, which can start with "./"
        :return: The joined path of project_dir and path_and_file without "./" prefix
        """
        # Removing "./" prefix if it exists
        if path_and_file.startswith("./"):
            path_and_file = path_and_file[2:]

        Makefile_list_of_files.append(path_and_file)
        return project_Path / path_and_file

    # Output all CFunctions to file. Keep track of all functions being output for the Makefile.
    list_of_uniq_functions: List[str] = []
    for name, CFunction in CFunction_dict.items():
        # Convention: Output all C functions starting in [CoordSystem] into the CoordSystem/ subdirectory.
        if "__rfm__" in name:
            subdir = name.split("__rfm__")[-1]
            (project_Path / subdir).mkdir(exist_ok=True)
            with open(
                add_to_Makefile(project_Path, f"{subdir}/{name}.c"),
                "w",
                encoding="utf-8",
            ) as file:
                file.write(CFunction.full_function)
        elif CFunction.subdirectory != "":
            subdir = CFunction.subdirectory
            with open(
                add_to_Makefile(project_Path, f"{subdir}/{name}.c"),
                "w",
                encoding="utf-8",
            ) as file:
                file.write(CFunction.full_function)
        else:
            with open(
                add_to_Makefile(project_Path, f"{name}.c"),
                "w",
                encoding="utf-8",
            ) as file:
                file.write(CFunction.full_function)
        list_of_uniq_functions += [name]

    # Output BHaH_function_prototypes.h
    with open(
        project_Path / "BHaH_function_prototypes.h", "w", encoding="utf-8"
    ) as file:
        for key in sorted(CFunction_dict.keys()):
            file.write(f"{CFunction_dict[key].function_prototype}\n")

    CFLAGS_dict = {
        "default": "-O2 -march=native -g -Wall -Wno-unused-variable",
        # FASTCFLAGS: -O3 causes AVX-2+ SIMD optimizations to be used on MoL update loops. -O2 drops to SSE2
        "fast": "-O3 -funroll-loops -march=native -g -Wall -Wno-unused-variable -std=gnu99",
        # DEBUGCFLAGS: OpenMP requires -fopenmp, and when disabling -fopenmp, unknown pragma warnings appear. -Wunknown-pragmas silences these warnings
        "debug": "-O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas",
    }

    if CC == "gcc":
        for key in CFLAGS_dict:
            CFLAGS_dict[key] += " -std=gnu99"

    if addl_CFLAGS is not None:
        if not isinstance(addl_CFLAGS, list):
            raise ValueError(
                "Error: output_CFunctions_function_prototypes_and_construct_Makefile(): addl_CFLAGS must be a list!"
            )
        for FLAG in addl_CFLAGS:
            for key in CFLAGS_dict:
                CFLAGS_dict[key] += f" {FLAG}"

    OBJ_FILES_str = "OBJ_FILES ="
    for c_file in sorted(Makefile_list_of_files):
        OBJ_FILES_str += " " + c_file.replace(".c", ".o")

    LDFLAGS_str = "LDFLAGS ="
    if addl_libraries is not None:
        if not isinstance(addl_libraries, list):
            raise ValueError(
                "Error: output_CFunctions_function_prototypes_and_construct_Makefile(): addl_libraries must be a list!"
            )
        for lib in addl_libraries:
            LDFLAGS_str += f" {lib}"
    LDFLAGS_str += " -lm"

    CFLAGS_str = f"CFLAGS = {CFLAGS_dict[compiler_opt_option]}\n"
    del CFLAGS_dict[compiler_opt_option]
    for value in CFLAGS_dict.values():
        CFLAGS_str += f"#CFLAGS = {value}\n"

    INCLUDEDIRS_str = "INCLUDEDIRS ="
    if include_dirs is not None:
        if not isinstance(include_dirs, list):
            raise TypeError(
                "Error: output_CFunctions_function_prototypes_and_construct_Makefile(): include_dirs must be a list!"
            )
        INCLUDEDIRS_str = " ".join(f"-I{include_dir}" for include_dir in include_dirs)

    # Below code is responsible for either writing a Makefile or a backup shell script depending on the conditions
    Makefile_str = f"""
CC = {CC}
{CFLAGS_str}
{INCLUDEDIRS_str}
{LDFLAGS_str}

# Check for OpenMP support
OPENMP_FLAG = -fopenmp
COMPILER_SUPPORTS_OPENMP := $(shell echo | $(CC) $(OPENMP_FLAG) -E - >/dev/null 2>&1 && echo YES || echo NO)

ifeq ($(COMPILER_SUPPORTS_OPENMP), YES)
    CFLAGS += $(OPENMP_FLAG)
    LDFLAGS += $(OPENMP_FLAG)  # -lgomp does not work with clang in termux
endif

{OBJ_FILES_str}

all: {exec_name}

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

{exec_name}: $(OBJ_FILES)
	$(CC) $^ -o $@ $(LDFLAGS)

clean:
	rm -f *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png {exec_name}
"""
    makefile_path = Path(project_Path) / "Makefile"
    with makefile_path.open("w", encoding="utf-8") as Makefile:
        Makefile.write(Makefile_str)


def compile_Makefile(
    project_dir: str,
    project_name: str,
    exec_name: str,
    compiler_opt_option: str = "fast",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "gcc",
    attempt: int = 1,
) -> None:
    """
    Compile the Makefile using the given parameters and options.

    :param project_dir: Root directory for C code.
    :param project_name: Name of the project.
    :param exec_name: Name of the executable.
    :param compiler_opt_option: Compiler optimization option (default: "fast").
    :param addl_CFLAGS: Additional CFLAGS.
    :param addl_libraries: Additional libraries.
    :param CC: C compiler (default: "gcc").
    :param attempt: Compilation attempt number (default: 1).
    """

    if shutil.which(CC) is None:
        raise FileNotFoundError(f"{CC} C compiler is not found")
    if shutil.which("make") is None:
        raise FileNotFoundError("make is not found")

    # Assuming output_CFunctions_function_prototypes_and_construct_Makefile is defined elsewhere
    output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_name=exec_name,
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

    exec_path = Path(project_dir).joinpath(exec_name)
    if not exec_path.is_file() and attempt == 1:
        print(
            "Optimized compilation FAILED. Removing optimizations (including OpenMP) and retrying with debug enabled..."
        )

        # First clean up object files.
        for obj_file in Path(project_dir).glob("*.o"):
            obj_file.unlink()

        # Then retry compilation (recursion)
        compile_Makefile(
            project_dir=project_dir,
            project_name=project_name,
            exec_name=exec_name,
            compiler_opt_option="debug",
            addl_CFLAGS=addl_CFLAGS,
            addl_libraries=addl_libraries,
            CC=CC,
            attempt=2,
        )
    if not exec_path.is_file() and attempt == 2:
        raise SystemExit("Compilation failed")
    print("Finished compilation.")

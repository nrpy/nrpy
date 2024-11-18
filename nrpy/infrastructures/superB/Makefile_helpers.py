"""
Construct Makefile for superB C-code project, based on CFunctions registered in the c_function NRPy+ module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
        Nishita Jadoo
        njadoo **at** uidaho **dot* edu

superB:
-removed "if not create_lib and "main" not in CFunction_dict: .."
-removed "if shutil.which(CC) is None: ..."
-changed file extensions from .c to .cpp
-changed "CC ?= {CC}" to "CC = {CC}"
-changed "	$(CC) $^ -o $@ $(LDFLAGS)" to  "$(CC) -language charm++ $^ -o $@ $(LDFLAGS)"
-added compilation and linking of timestepping.cpp. timestepping.ci, main.cpp and main.ci
-added *.decl.h *.def.h charmrun to "clean: .."
-removed openmp "# Check for OpenMP support..."
"""

import multiprocessing
import os
import platform
import shutil
import subprocess
from pathlib import Path
from typing import List, Optional

from nrpy.c_function import CFunction_dict


def output_CFunctions_function_prototypes_and_construct_Makefile(
    project_dir: str,
    project_name: str,
    exec_or_library_name: str = "",
    compiler_opt_option: str = "default",
    addl_CFLAGS: Optional[List[str]] = None,
    addl_libraries: Optional[List[str]] = None,
    CC: str = "autodetect",
    create_lib: bool = False,
    include_dirs: Optional[List[str]] = None,
) -> None:
    """
    Output C functions registered to CFunction_dict and construct a Makefile for compiling C code.

    :param project_dir: The root directory of the C project.
    :param project_name: Name of the C project.
    :param exec_or_library_name: The name of the executable. If set to empty string, same as project_name.
    :param compiler_opt_option: Compiler optimization option. Defaults to "default". Other options: "fast" and "debug"
    :param addl_CFLAGS: Additional compiler flags. Must be a list.
    :param addl_libraries: Additional libraries to link. Must be a list.
    :param CC: C compiler to use. Defaults to "autodetect" (clang if using Darwin, gcc otherwise)
    :param create_lib: Whether to create a library. Defaults to False.
    :param include_dirs: List of include directories. Must be a list.

    :raises TypeError: If addl_CFLAGS or include_dirs are not lists.
    :raises ValueError: If addl_CFLAGS or addl_libraries are specified incorrectly, or if OS unsupported.
    """
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    if exec_or_library_name == "":
        exec_or_library_name = project_name

    os_name = platform.system()
    if create_lib:
        if os_name == "Linux":
            ext = ".so"
        elif os_name == "Darwin":
            ext = ".dylib"
        else:
            raise ValueError(f"Sorry, {os_name} operating system not supported.")
        if not exec_or_library_name.endswith(ext):
            exec_or_library_name += ext

        def add_flag(flag_list: Optional[List[str]], flag: str) -> List[str]:
            """
            Check if a flag is in the list, add it if not.

            :param flag_list: The list to which the flag should be added.
            :param flag: The flag to add to the list.

            :return: The updated list with the flag added, if it was not already present.
            """
            if not flag_list:
                flag_list = []
            if flag not in flag_list:
                flag_list += [flag]
            return flag_list

        addl_CFLAGS = add_flag(addl_CFLAGS, "-fPIC")
        addl_libraries = add_flag(addl_libraries, "-fPIC")
        addl_libraries = add_flag(addl_libraries, "-shared")

    Makefile_list_of_files: List[str] = []

    def add_to_Makefile(project_Path: Path, subdir_Path: Path, file_name: str) -> Path:
        """
        Add the given path and file to the Makefile list of files and return the joined path.

        :param project_Path: The Path of the project.
        :param subdir_Path: The subdirectory path within the project.
        :param file_name: The filename.
        :return: Destination Path for the generated code, including project & subdirectory paths: e.g., path/to/code.c.
        """
        Makefile_list_of_files.append(str(subdir_Path / file_name))
        return project_Path / subdir_Path / file_name

    # Output all CFunctions to file. Keep track of all functions being output for the Makefile.
    list_of_uniq_functions: List[str] = []
    for name, CFunction in CFunction_dict.items():
        if CFunction.subdirectory:
            subdir_Path = Path(CFunction.subdirectory)
            (project_Path / subdir_Path).mkdir(parents=True, exist_ok=True)
            with open(
                add_to_Makefile(project_Path, subdir_Path, f"{name}.cpp"),
                "w",
                encoding="utf-8",
            ) as file:
                file.write(CFunction.full_function)
        else:
            with open(
                add_to_Makefile(project_Path, Path("."), f"{name}.cpp"),
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

    if CC == "autodetect":
        if os_name == "Darwin":
            CC = "clang"
        else:
            CC = "gcc"

    # ~ if shutil.which(CC) is None:
    # ~ raise FileNotFoundError(f"{CC} C compiler is not found")

    CFLAGS_dict = {
        "default": "-O2 -march=native -g",
        "fast": "-O3 -funroll-loops -march=native -g -Wall",
        "debug": "-O2 -g -Wall -Wno-unknown-pragmas",
    }

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
        OBJ_FILES_str += " " + c_file.replace(".cpp", ".o")

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
        INCLUDEDIRS_str += " " + " ".join(
            f"-I{include_dir}" for include_dir in include_dirs
        )

    # Below code is responsible for either writing a Makefile or a backup shell script depending on the conditions
    Makefile_str = f"""CC = {CC}
{CFLAGS_str}
VALGRIND_CFLAGS = {CFLAGS_dict["debug"]}
{INCLUDEDIRS_str}
{LDFLAGS_str}

{OBJ_FILES_str}

all: {exec_or_library_name}

%.o: %.cpp $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

timestepping.o: timestepping.cpp timestepping.h main.h timestepping.def.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

timestepping.h: timestepping.decl.h

timestepping.decl.h timestepping.def.h: timestepping.ci
	$(CC) $(CFLAGS) $(INCLUDEDIRS) timestepping.ci

main.o: main.cpp main.h main.decl.h main.def.h timestepping.decl.h
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

main.h: timestepping.decl.h main.decl.h

main.decl.h main.def.h: main.ci
	$(CC) $(CFLAGS) $(INCLUDEDIRS) main.ci

{exec_or_library_name}: $(OBJ_FILES) timestepping.o main.o
	$(CC) -language charm++ $^ -o $@ $(LDFLAGS)

# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png {exec_or_library_name} *.decl.h *.def.h charmrun
	$(RM) -r log
"""
    # Add the valgrind target
    Makefile_str += """valgrind: clean
\t$(MAKE) CFLAGS="$(VALGRIND_CFLAGS)" all
"""
    if not create_lib:
        Makefile_str += f"""\tvalgrind --track-origins=yes --leak-check=full -s ./{exec_or_library_name}

"""
    makefile_path = Path(project_Path) / "Makefile"
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
    :param compiler_opt_option: Compiler optimization option (default: "fast").
    :param addl_CFLAGS: Additional CFLAGS (default: None).
    :param addl_libraries: Additional libraries (default: None).
    :param CC: C compiler (default: "autodetect").
    :param attempt: Compilation attempt number (default: 1).

    :raises FileNotFoundError: If the C compiler or make is not found.
    :raises SystemExit: If compilation fails after two attempts.
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

    # Assuming output_CFunctions_function_prototypes_and_construct_Makefile is defined elsewhere
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
            exec_or_library_name=exec_or_library_name,
            compiler_opt_option="debug",
            addl_CFLAGS=addl_CFLAGS,
            addl_libraries=addl_libraries,
            CC=CC,
            attempt=2,
        )
    if not exec_path.is_file() and attempt == 2:
        raise SystemExit("Compilation failed")
    print("Finished compilation.")

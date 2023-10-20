"""
Construct Makefile for BHaH C-code project, based on CFunctions registered in the c_function NRPy+ module.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path
import os
import shutil
import subprocess
import platform
import multiprocessing
from typing import List, Optional
import cpuinfo  # type: ignore

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

    :raises SystemExit: Exits if errors are encountered.
    """
    if not create_lib and "main" not in CFunction_dict:
        raise SystemExit(
            "output_CFunctions_function_prototypes_and_construct_Makefile() error: C codes will not compile if main() function not defined!\n"
            '    Make sure that the main() function registered to CFunction_dict has name "main".'
        )

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
        if not exec_or_library_name.endswith(ext):
            exec_or_library_name += ext

        def add_flag(flag_list: Optional[List[str]], flag: str) -> List[str]:
            """Check if a flag is in the list, add it if not."""
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
                add_to_Makefile(project_Path, subdir_Path, f"{name}.c"),
                "w",
                encoding="utf-8",
            ) as file:
                file.write(CFunction.full_function)
        else:
            with open(
                add_to_Makefile(project_Path, Path("."), f"{name}.c"),
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

    if shutil.which(CC) is None:
        raise FileNotFoundError(f"{CC} C compiler is not found")

    CFLAGS_dict = {
        "default": "-O2 -march=native -g -Wall -Wno-unused-variable",
        # FASTCFLAGS: -O3 causes AVX-2+ SIMD optimizations to be used on MoL update loops. -O2 drops to SSE2
        "fast": "-O3 -funroll-loops -march=native -g -Wall -Wno-unused-variable -std=gnu99",
        # DEBUGCFLAGS: OpenMP requires -fopenmp, and when disabling -fopenmp, unknown pragma warnings appear. -Wunknown-pragmas silences these warnings
        "debug": "-O2 -g -Wall -Wno-unused-variable -Wno-unknown-pragmas",
    }

    if CC == "gcc":
        cpu_info = cpuinfo.get_cpu_info()
        if "flags" in cpu_info:
            if any("avx512" in flag for flag in cpu_info["flags"]):
                # Experiment:
                # if cpu_info.get("l2_cache_size"):
                #     try:
                #         l2_cache_size_KB = int(
                #             int(cpu_info.get("l2_cache_size")) / 1024
                #         )
                #         for key, value in CFLAGS_dict.items():
                #             CFLAGS_dict[
                #                 key
                #             ] += f" --param l2-cache-size={l2_cache_size_KB}"
                #     except ValueError:
                #         # Sometimes cpu_info.get("l2_cache_size") returns a value that has explicit units, like 2.5MiB,
                #         #  ignore these cases
                #         pass

                # -march=native hangs when using GCC on
                avx512_features = [
                    "avx512f",
                    "avx512pf",
                    "avx512er",
                    "avx512cd",
                    "avx512vl",
                    "avx512bw",
                    "avx512dq",
                    "avx512ifma",
                    "avx512vbmi",
                ]
                # manually add the avx-512 featureset
                avx512_compileflags = ""
                for feature in avx512_features:
                    if any(
                        flag.replace("_", "") == feature for flag in cpu_info["flags"]
                    ):
                        avx512_compileflags += f"-m{feature} "
                for key, value in CFLAGS_dict.items():
                    CFLAGS_dict[key] = CFLAGS_dict[key].replace(
                        "-march=native", avx512_compileflags
                    )

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
    Makefile_str = f"""CC ?= {CC}  # assigns the value CC to {CC} only if environment variable CC is not already set
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

all: {exec_or_library_name}

%.o: %.c $(COMMON_HEADERS)
	$(CC) $(CFLAGS) $(INCLUDEDIRS) -c $< -o $@

{exec_or_library_name}: $(OBJ_FILES)
	$(CC) $^ -o $@ $(LDFLAGS)

# Use $(RM) to be cross-platform compatible.
clean:
	$(RM) *.o */*.o *~ */*~ ./#* *.txt *.dat *.avi *.png {exec_or_library_name}
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

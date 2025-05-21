"""
Construct BHaH_gpu_defines from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

import sys
from pathlib import Path
from typing import Any, Dict, Union

import nrpy.grid as gri
import nrpy.params as par
from nrpy.helpers.generic import clang_format

if "DEVICE_THREAD_MACROS" not in par.glb_extras_dict:
    par.glb_extras_dict["DEVICE_THREAD_MACROS"] = {}
par.glb_extras_dict["DEVICE_THREAD_MACROS"].update(
    {
        "BHAH_THREADS_IN_X_DIR_DEFAULT": 32,
        "BHAH_THREADS_IN_Y_DIR_DEFAULT": 1,
        "BHAH_THREADS_IN_Z_DIR_DEFAULT": 1,
    }
)


def generate_declaration_str(
    decl_dict: Dict[str, Dict[str, str]], prefix: str = ""
) -> str:
    """
    Generate block string of header definitions.

    :param decl_dict: Dictionary of definitions and their properties
    :param prefix: optional prefix
    :returns: str
    """
    if prefix != "":
        prefix += " "
    decl_str: str = ""
    for var, sub_dict in decl_dict.items():
        suffix = f"[{sub_dict['array_size']}]" if sub_dict["array_size"] else ""
        # Standard declarations str
        decl_str += f"{sub_dict['comment']}"
        decl_str += f"{prefix} {sub_dict['type']} {var}{suffix};\n"
    return decl_str


class CUDA_BHaH_device_defines_h:
    r"""
    Generate and write to file the BHaH_CUDA_defines.h file.

    :param project_dir: Location to write file to
    :param additional_declarations_dict: Dictionary storing additional declaration dictionaries
    :param additional_macros_str: Block string of additional macro definitions
    :param num_streams: Number of CUDA streams to use
    :param nghosts: Number of ghost zones for the FD stencil

    >>> import nrpy.c_function as cfc
    >>> from nrpy.helpers.generic import validate_strings
    >>> cfc.CFunction_dict.clear()
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> d = CUDA_BHaH_device_defines_h(project_dir)
    >>> d.generate_output_str()
    >>> generated_str = d.file_output_str
    >>> validation_desc = "CUDA_BHaH_device_defines_h"
    >>> validate_strings(generated_str, validation_desc, file_ext="cu")
    """

    def __init__(
        self,
        project_dir: str,
        additional_declarations_dict: Union[Dict[str, Any], None] = None,
        additional_macros_str: Union[str, None] = None,
        num_streams: int = 3,
        nghosts: Union[int, None] = None,
        set_parity_on_aux: bool = False,
        set_parity_on_auxevol: bool = False,
    ) -> None:
        self.project_Path = Path(project_dir)
        self.project_Path.mkdir(parents=True, exist_ok=True)
        self.num_streams = num_streams
        self.additional_decl_dict = additional_declarations_dict
        self.additional_macros_str = additional_macros_str
        self.bhah_CUDA_defines_filename = "BHaH_device_defines.h"
        self.NGHOSTS = nghosts

        # Standard macros str
        self.macro_str = f"""
// Standard macro definitions
// We include the macro definition NUM_STREAMS since it is used for calculations in various
// algorithms in addition to defining the streams array
#define NUM_STREAMS {self.num_streams}
"""
        standard_decl_dict = {
            "d_params": {
                "type": "__constant__ params_struct",
                "array_size": "NUM_STREAMS",
                "comment": "// Device storage for grid parameters\n",
            },
            "d_commondata": {
                "type": "__constant__ commondata_struct",
                "array_size": "",
                "comment": "// Device storage for commondata\n",
            },
            "streams": {
                "type": "cudaStream_t",
                "array_size": "NUM_STREAMS",
                "comment": "",
            },
            "GPU_N_SMS": {
                "type": "size_t",
                "array_size": "",
                "comment": "",
            },
        }
        # First add human-readable gridfunction aliases (grid.py) to BHaH_defines dictionary.
        (
            evolved_variables_list,
            auxiliary_variables_list,
            auxevol_variables_list,
        ) = gri.BHaHGridFunction.gridfunction_lists()[0:3]

        if len(evolved_variables_list) > 0:
            standard_decl_dict["d_evol_gf_parity"] = {
                "type": "__constant__ int8_t",
                "array_size": f"{len(evolved_variables_list)}",
                "comment": "// Device storage for evolved gridfunction parity\n",
            }

        if set_parity_on_aux:
            if len(auxiliary_variables_list) > 0:
                standard_decl_dict["d_aux_gf_parity"] = {
                    "type": "__constant__ int8_t",
                    "array_size": f"{len(auxiliary_variables_list)}",
                    "comment": "// Device storage for evolved gridfunction parity\n",
                }

        if set_parity_on_auxevol:
            if len(auxevol_variables_list) > 0:
                standard_decl_dict["d_auxevol_gf_parity"] = {
                    "type": "__constant__ int8_t",
                    "array_size": f"{len(auxevol_variables_list)}",
                    "comment": "// Device storage for evolved gridfunction parity\n",
                }

        evolved_variables_list: list[str]
        (
            evolved_variables_list,
            _auxiliary_variables_list,
            _auxevol_variables_list,
        ) = gri.BHaHGridFunction.gridfunction_lists()
        # This device storage is only needed by some problems
        if evolved_variables_list:
            standard_decl_dict["d_gridfunctions_wavespeed"] = {
                "type": "__constant__ REAL",
                "array_size": "NUM_EVOL_GFS",
                "comment": "",
            }
            standard_decl_dict["d_gridfunctions_f_infinity"] = {
                "type": "__constant__ REAL",
                "array_size": "NUM_EVOL_GFS",
                "comment": "",
            }
        self.combined_decl_dict = standard_decl_dict
        self.decl_str: str = "// Standard declarations\n"
        self.decl_str += generate_declaration_str(
            self.combined_decl_dict, prefix="extern"
        )

        self.file_output_str = ""
        self.generate_output_str()
        self.write_to_file()

    def combine_declarations_dicts(self) -> None:
        """Add additional_decl_dict to combined_decl_dict."""
        if not self.additional_decl_dict is None:
            for k, v in self.additional_decl_dict.items():
                self.combined_decl_dict[k] = v

    def generate_output_str(self) -> None:
        """Generate block output str to prepare writing to file."""
        self.file_output_str = """// BHaH core header file, automatically generated from cuda.output_BHaH_defines_h,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""

        self.file_output_str += self.macro_str
        if self.additional_macros_str:
            self.file_output_str += (
                "\n\n// Additional Macros\n" + self.additional_macros_str
            )

        self.file_output_str += "\n\n" + self.decl_str

        if self.additional_decl_dict:
            self.file_output_str += "\n\n// Additional Declarations\n"
            self.file_output_str += generate_declaration_str(
                self.additional_decl_dict, prefix="extern"
            )
            self.combine_declarations_dicts()

        self.file_output_str += (
            "\n\n"
            + r"""
// CUDA Error checking macro only active if compiled with -DDEBUG
// Otherwise additional synchronization overhead will occur
#ifdef DEBUG
#define cudaCheckErrors(v, msg) \
    do { \
        cudaError_t __err = cudaGetLastError(); \
        if (__err != cudaSuccess) { \
            fprintf(stderr, "Fatal error: %s %s (%s at %s:%d)\n", \
                #v, msg, cudaGetErrorString(__err), \
                __FILE__, __LINE__); \
            fprintf(stderr, "*** FAILED - ABORTING\n"); \
            exit(1); \
        } \
    } while (0);
#else
#define cudaCheckErrors(v, msg)
#endif

#define BHAH_MEMCPY_HOST_TO_DEVICE(dest_ptr, src_ptr, sz) cudaMemcpy(dest_ptr, src_ptr, sz, cudaMemcpyHostToDevice);

"""
        )
        for k, thread_cnt in par.glb_extras_dict["DEVICE_THREAD_MACROS"].items():
            self.file_output_str += f"#define {k.upper()} {thread_cnt}\n"
        self.file_output_str = clang_format(self.file_output_str)

    def write_to_file(self) -> None:
        """Write file_output_str to header file."""
        bhah_gpu_defines_file = self.project_Path / self.bhah_CUDA_defines_filename
        with bhah_gpu_defines_file.open("w", encoding="utf-8") as file:
            file.write(self.file_output_str)


class BHaH_CUDA_global_init_h:
    r"""
    Generate and write to file the BHaH_CUDA_global_init.h file.

    :param project_dir: Location to write file to
    :param declarations_dict: Dictionary storing declaration dictionaries

    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> gpu_d = CUDA_BHaH_device_defines_h(project_dir)
    >>> gpu_init = BHaH_CUDA_global_init_h(project_dir,gpu_d.combined_decl_dict)
    >>> print(gpu_init.file_output_str)
    // BHaH core header file, automatically generated from cuda.output_BHaH_defines_h,
    //    DO NOT EDIT THIS FILE BY HAND.
    <BLANKLINE>
    // Initialize streams
    for (int i = 0; i < NUM_STREAMS; ++i) {
      cudaStreamCreate(&streams[i]);
    }
    // Copy parity array to device __constant__ memory
    cudaMemcpyToSymbol(d_evol_gf_parity, evol_gf_parity, 24 * sizeof(int8_t));
    cudaCheckErrors(copy, "Copy to d_evol_gf_parity failed");
    <BLANKLINE>
    """

    def __init__(
        self,
        project_dir: str,
        declarations_dict: Dict[str, Dict[str, str]],
        set_parity_on_aux: bool = False,
        set_parity_on_auxevol: bool = False,
    ) -> None:
        self.project_Path = Path(project_dir)
        self.project_Path.mkdir(parents=True, exist_ok=True)
        self.declarations_dict = declarations_dict
        self.filename = "BHaH_CUDA_global_init.h"

        self.file_output_str = """// BHaH core header file, automatically generated from cuda.output_BHaH_defines_h,
//    DO NOT EDIT THIS FILE BY HAND.\n\n"""

        self.file_output_str += """\n\n
// Initialize streams
for(int i = 0; i < NUM_STREAMS; ++i) {
    cudaStreamCreate(&streams[i]);
}"""

        # Ensure necessary global arrays, currently stored in BHaH_defines.h, are copied to
        # device __constant__ memory
        for constant_ary in [
            "d_evol_gf_parity",
            "d_aux_gf_parity",
            "d_auxevol_gf_parity",
            "d_gridfunctions_wavespeed",
            "d_gridfunctions_f_infinity"]:
            if constant_ary in declarations_dict.keys():
                host_ary_name = constant_ary.replace("d_", "")
                host_ary_type = declarations_dict[constant_ary]["type"].replace("__constant__", "").strip()
                self.file_output_str += f"""// Copy {host_ary_name} to device __constant__ memory
cudaMemcpyToSymbol({constant_ary}, {host_ary_name}, {declarations_dict[constant_ary]['array_size']} * sizeof({host_ary_type}));
cudaCheckErrors(copy, "Copy to {constant_ary} failed");
"""

        self.file_output_str = clang_format(self.file_output_str)
        self.write_to_file()

    def write_to_file(self) -> None:
        """Write file_output_str to header file."""
        output_file = self.project_Path / self.filename
        with output_file.open("w", encoding="utf-8") as file:
            file.write(self.file_output_str)


class BHaH_CUDA_global_defines_h:
    r"""
    Generate and write to file the BHaH_device_defines_h file.

    :param project_dir: Location to write file to
    :param declarations_dict: Dictionary storing declaration dictionaries

    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> gpu_d = CUDA_BHaH_device_defines_h(project_dir)
    >>> gpu_init = BHaH_CUDA_global_init_h(project_dir,gpu_d.combined_decl_dict)
    >>> print(gpu_init.file_output_str)
    // BHaH core header file, automatically generated from cuda.output_BHaH_defines_h,
    //    DO NOT EDIT THIS FILE BY HAND.
    <BLANKLINE>
    // Initialize streams
    for (int i = 0; i < NUM_STREAMS; ++i) {
      cudaStreamCreate(&streams[i]);
    }
    // Copy parity array to device __constant__ memory
    cudaMemcpyToSymbol(d_evol_gf_parity, evol_gf_parity, 24 * sizeof(int8_t));
    cudaCheckErrors(copy, "Copy to d_evol_gf_parity failed");
    <BLANKLINE>
    """

    def __init__(
        self,
        project_dir: str,
        declarations_dict: Dict[str, Dict[str, str]],
        **_: Any,
    ) -> None:
        self.project_Path = Path(project_dir)
        self.project_Path.mkdir(parents=True, exist_ok=True)
        self.declarations_dict = declarations_dict
        self.filename = "BHaH_global_device_defines.h"

        self.file_output_str = """// BHaH core header file, automatically generated from cuda.output_BHaH_defines_h,
//    DO NOT EDIT THIS FILE BY HAND.\n\n
#ifndef __BHAH_GLOBAL_DEVICE_DEFINES_H__
#include "BHaH_defines.h"
"""

        self.file_output_str += generate_declaration_str(self.declarations_dict)

        self.file_output_str += "#endif // __BHAH_GLOBAL_DEVICE_DEFINES_H__\n\n"
        self.file_output_str = clang_format(self.file_output_str)
        self.write_to_file()

    def write_to_file(self) -> None:
        """Write file_output_str to header file."""
        output_file = self.project_Path / self.filename
        with output_file.open("w", encoding="utf-8") as file:
            file.write(self.file_output_str)


def output_device_headers(
    project_dir: str,
    additional_declarations_dict: Union[Dict[str, Any], None] = None,
    additional_macros_str: Union[str, None] = None,
    num_streams: int = 3,
    nghosts: Union[int, None] = None,
) -> str:
    """
    Generate device specific header files.

    e.g. CUDA headers store global __constant__ variable declarations, initialize device storage,
    and include device specific macros, e.g. NUM_STREAMS.

    :param project_dir: Project directory
    :param additional_declarations_dict: Dictionary storing additional declaration dictionaries
    :param additional_macros_str: Block string of additional macro definitions
    :param num_streams: Number of CUDA streams to use
    :param nghosts: FD stencil radius
    :returns: header filename
    """
    parallelization = par.parval_from_str("parallelization")
    device_header_filename = ""
    if parallelization == "cuda":
        gpu_defines = CUDA_BHaH_device_defines_h(
            project_dir,
            additional_declarations_dict=additional_declarations_dict,
            additional_macros_str=additional_macros_str,
            num_streams=num_streams,
            nghosts=nghosts,
        )

        BHaH_CUDA_global_defines_h(
            project_dir,
            gpu_defines.combined_decl_dict,
        )

        BHaH_CUDA_global_init_h(
            project_dir,
            gpu_defines.combined_decl_dict,
        )
        device_header_filename = gpu_defines.bhah_CUDA_defines_filename

    return device_header_filename


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

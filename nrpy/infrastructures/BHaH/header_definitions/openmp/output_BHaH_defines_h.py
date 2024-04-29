"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from pathlib import Path
from typing import Optional, Dict, List

import nrpy.params as par
import nrpy.grid as gri
from nrpy.infrastructures.BHaH import griddata_commondata
from nrpy.helpers.generic import clang_format
from nrpy.infrastructures.BHaH.header_definitions.base_output_BHaH_defines_h import base_output_BHaH_defines_h

class output_BHaH_defines_h(base_output_BHaH_defines_h):
    r"""
    Output C code header file with macro definitions and other configurations for the project.

    :param project_dir: Directory where the project C code is output
    :param additional_includes: Additional header files to be included in the output
    :param REAL_means: The floating-point type to be used in the C code (default is "double")
    :param enable_simd: Flag to enable Single Instruction Multiple Data (SIMD) optimizations
    :param enable_rfm_precompute: A boolean value reflecting whether reference metric precomputation is enabled.
    :param fin_NGHOSTS_add_one_for_upwinding_or_KO: Option to add one extra ghost zone for upwinding
    :param supplemental_defines_dict: Additional key-value pairs to be included in the output file
    :param clang_format_options: Options for clang formatting.

    >>> from nrpy.infrastructures.BHaH.MoLtimestepping import MoL
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> MoL.register_CFunctions(register_MoL_step_forward_in_time=False)
    >>> project_dir = Path("/tmp", "tmp_BHaH_defines_h")
    >>> project_dir.mkdir(parents=True, exist_ok=True)
    >>> C=output_BHaH_defines_h(project_dir=str(project_dir))
    >>> expected_string = decompress_base64_to_string("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")
    >>> returned_string = (project_dir / "BHaH_defines.h").read_text()
    >>> if returned_string != expected_string:
    ...    compressed_str = compress_string_to_base64(returned_string)
    ...    error_message = "Trusted BHaH_defines.h string changed!\n"
    ...    error_message += "Here's the diff:\n" + diff_strings(expected_string, returned_string) + "\n"
    ...    raise ValueError(error_message + f"base64-encoded output: {compressed_str}")
    """    
    def __init__(
        self,
        project_dir: str,
        additional_includes: Optional[List[str]] = None,
        REAL_means: str = "double",
        enable_simd: bool = True,
        enable_rfm_precompute: bool = True,
        fin_NGHOSTS_add_one_for_upwinding_or_KO: bool = False,
        supplemental_defines_dict: Optional[Dict[str, str]] = None,
        clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
    ) -> None:

        super().__init__(
            project_dir,
            additional_includes=additional_includes,
            REAL_means=REAL_means,
            enable_simd=enable_simd,
            enable_rfm_precompute=enable_rfm_precompute,
            fin_NGHOSTS_add_one_for_upwinding_or_KO=fin_NGHOSTS_add_one_for_upwinding_or_KO,
            supplemental_defines_dict=supplemental_defines_dict,
            clang_format_options=clang_format_options,
        )
        self.register_define_blocks()
        self.generate_output_str()
        self.write_to_file()

if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

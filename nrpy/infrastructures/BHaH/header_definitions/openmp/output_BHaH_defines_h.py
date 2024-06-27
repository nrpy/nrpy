"""
Construct BHaH_defines.h from data registered to griddata_commondata, CodeParameters, and NRPyParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import sys
from typing import Dict, List, Optional

from nrpy.infrastructures.BHaH.header_definitions.base_output_BHaH_defines_h import (
    base_output_BHaH_defines_h,
)


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

    >>> from pathlib import Path
    >>> from nrpy.infrastructures.BHaH.MoLtimestepping.openmp import MoL
    >>> import nrpy.finite_difference as fin
    >>> from nrpy.helpers.generic import compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> _ = MoL.register_CFunctions(register_MoL_step_forward_in_time=False)
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

"""
Functions for setting params_struct and commondata_struct parameters.

Set to default values specified when registering them within NRPy+'s CodeParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path

import nrpy.params as par
from nrpy.helpers.generic import clang_format


def write_CodeParameters_simd_h_files(
    project_dir: str,
    thorn_name: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> None:
    r"""
    Generate C code to set C parameter constants and writes them to files.

    :param project_dir: The path of the project directory.
    :param clang_format_options: Options for clang_format.
    :return: None

    Doctests:
    >>> par.set_parval_from_str("Infrastructure", "ETLegacy")
    >>> import nrpy.c_function as cfc
    >>> _, __ = par.register_CodeParameters("CCTK_REAL", "CodeParameters_c_files", ["a", "pi_three_sigfigs"], [1, 3.14])
    >>> _int = par.register_CodeParameter("CCTK_INT", "CodeParameters_c_files", "j0", 1)
    >>> _leaveitbe = par.register_CodeParameter("CCTK_REAL", "CodeParameters_c_files", "leaveitbe", add_to_parfile=False, add_to_set_CodeParameters_h=False)
    >>> cfc.CFunction_dict.clear()
    >>> project_dir = Path("/tmp/et_project/")
    >>> write_CodeParameters_simd_h_files(project_dir=str(project_dir), thorn_name="thorn")
    >>> print((project_dir / "thorn" / "src"/ "set_CodeParameters-simd.h").read_text())
    const CCTK_REAL NOSIMDa = CCTK_ParameterGet("a", "CodeParameters_c_files", NULL);                               // CodeParameters_c_files::a
    const REAL_SIMD_ARRAY a = ConstSIMD(NOSIMDa);                                                                   // CodeParameters_c_files::a
    const CCTK_INT j0 = params->j0;                                                                                 // CodeParameters_c_files::j0
    const CCTK_REAL NOSIMDpi_three_sigfigs = CCTK_ParameterGet("pi_three_sigfigs", "CodeParameters_c_files", NULL); // CodeParameters_c_files::pi_three_sigfigs
    const REAL_SIMD_ARRAY pi_three_sigfigs = ConstSIMD(NOSIMDpi_three_sigfigs);                                     // CodeParameters_c_files::pi_three_sigfigs
    <BLANKLINE>
    """
    # Create output directory if it doesn't already exist
    src_Path = Path(project_dir) / thorn_name / "src"
    src_Path.mkdir(parents=True, exist_ok=True)

    # Next, output header file for setting C parameters to current values within functions.
    # Step 4: Output set_CodeParameters-simd.h
    set_CodeParameters_SIMD_str = ""
    for CPname, CodeParam in sorted(
        par.glb_code_params_dict.items(), key=lambda x: x[0].lower()
    ):
        # SIMD does not support char arrays.
        if (
            CodeParam.add_to_set_CodeParameters_h
            and "char" not in CodeParam.c_type_alias
        ):
            struct = "commondata" if CodeParam.commondata else "params"
            CPtype = CodeParam.c_type_alias
            comment = f"  // {CodeParam.module}::{CPname}"
            if CPtype == "CCTK_REAL" or CPtype == "REAL":
                c_output = f"""const CCTK_REAL NOSIMD{CPname} = CCTK_ParameterGet("{CPname}", "{CodeParam.module}", NULL);{comment}\n"""
                c_output += f"const REAL_SIMD_ARRAY {CPname} = ConstSIMD(NOSIMD{CPname});{comment}\n"
                set_CodeParameters_SIMD_str += c_output
            else:
                c_output = f"const {CPtype} {CPname} = {struct}->{CPname};{comment}\n"
                set_CodeParameters_SIMD_str += c_output

    header_file_simd_path = src_Path / "set_CodeParameters-simd.h"
    with header_file_simd_path.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(
                set_CodeParameters_SIMD_str,
                clang_format_options=clang_format_options,
            )
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

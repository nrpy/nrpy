"""
Functions for setting params_struct and commondata_struct parameters
  to default values specified when registering them within NRPy+'s
  CodeParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path
from typing import List
import nrpy.params as par
import nrpy.c_function as cfc
from nrpy.helpers.generic import clang_format


def register_CFunctions_params_commondata_struct_set_to_default() -> None:
    """
    Register a C function to set code parameters to their default values.

    :param None: This function takes no parameters.

    >>> _, __ = par.register_CodeParameters("REAL", "CodeParameters_c_files", ["a", "pi_three_sigfigs"], [1, 3.14], commondata=True)
    >>> ___ = par.register_CodeParameter("#define", "CodeParameters_c_files", "b", 0)
    >>> _leaveitbe = par.register_CodeParameter("REAL", "CodeParameters_c_files", "leaveitbe", add_to_parfile=False, add_to_set_CodeParameters_h=False)
    >>> _str = par.register_CodeParameter("char[100]", "CodeParameters_c_files", "some_string", "cheese")
    >>> _bool = par.register_CodeParameter("bool", "CodeParameters_c_files", "BHaH_is_amazing", True, add_to_set_CodeParameters_h=True)
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions_params_commondata_struct_set_to_default()
    >>> print(cfc.CFunction_dict["params_struct_set_to_default"].full_function)
    #include "BHaH_defines.h"
    /*
     * Set params_struct to default values specified within NRPy+.
     */
    void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
      // Loop over params structs:
      for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
        params_struct *restrict params = &griddata[grid].params;
        // Set params_struct variables to default
        params->BHaH_is_amazing = true;               // CodeParameters_c_files::BHaH_is_amazing
        snprintf(params->some_string, 100, "cheese"); // CodeParameters_c_files::some_string
      }
    }
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["commondata_struct_set_to_default"].full_function)
    #include "BHaH_defines.h"
    /*
     * Set commondata_struct to default values specified within NRPy+.
     */
    void commondata_struct_set_to_default(commondata_struct *restrict commondata) {
    <BLANKLINE>
      // Set commondata_struct variables to default
      commondata->a = 1;                   // CodeParameters_c_files::a
      commondata->pi_three_sigfigs = 3.14; // CodeParameters_c_files::pi_three_sigfigs
    }
    <BLANKLINE>
    """
    for function_name in ["commondata_struct", "params_struct"]:
        includes = ["BHaH_defines.h"]
        desc = f"Set {function_name} to default values specified within NRPy+."
        c_type = "void"
        name = f"{function_name}_set_to_default"
        params = "commondata_struct *restrict commondata"
        if function_name == "params_struct":
            params += ", griddata_struct *restrict griddata"

        struct_list: List[str] = []  # List to store individual struct elements

        for parname, CodeParam in par.glb_code_params_dict.items():
            if (function_name == "commondata_struct" and CodeParam.commondata) or (
                function_name == "params_struct" and not CodeParam.commondata
            ):
                if CodeParam.add_to_parfile and CodeParam.add_to_set_CodeParameters_h:
                    struct = "commondata" if CodeParam.commondata else "params"
                    CPtype = CodeParam.c_type_alias
                    comment = f"  // {CodeParam.module}::{parname}"
                    defaultval = CodeParam.defaultvalue
                    if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                        chararray_size = CPtype.split("[")[1].replace("]", "")
                        c_output = f'snprintf({struct}->{parname}, {chararray_size}, "{defaultval}");{comment}\n'
                    elif isinstance(defaultval, (bool, int, float)):
                        c_output = f"{struct}->{parname} = {str(defaultval).lower()};{comment}\n"
                    else:
                        raise ValueError(
                            f"{CodeParam.defaultvalue} is not a valid default value. Note that Booleans must be capitalized."
                        )
                    struct_list.append(c_output)

        # Sort the lines alphabetically and join them with line breaks
        body = ""
        if function_name == "params_struct":
            body += r"""// Loop over params structs:
    for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
    params_struct *restrict params = &griddata[grid].params;
    """
            body += "// Set params_struct variables to default\n"
            body += "".join(sorted(struct_list))
            body += "}\n"
        else:
            body += "\n// Set commondata_struct variables to default\n"
            body += "".join(sorted(struct_list))

        cfc.register_CFunction(
            include_CodeParameters_h=False,
            includes=includes,
            desc=desc,
            c_type=c_type,
            name=name,
            params=params,
            body=body,
        )


def write_CodeParameters_h_files(
    project_dir: str,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 0}",
) -> None:
    """
    Generates C code to set C parameter constants and writes them to files.

    :param project_name: The name of the project directory.
    :param project_root_dir: The root directory of the project.
    :param clang_format_options: Options for clang_format.
    >>> project_dir = Path("/tmp/tmp_project/")
    >>> write_CodeParameters_h_files(str(project_dir))
    >>> print((project_dir / 'set_CodeParameters.h').read_text())
    const REAL a = commondata->a;                               // CodeParameters_c_files::a
    const bool BHaH_is_amazing = params->BHaH_is_amazing;       // CodeParameters_c_files::BHaH_is_amazing
    const REAL pi_three_sigfigs = commondata->pi_three_sigfigs; // CodeParameters_c_files::pi_three_sigfigs
    char some_string[100];
    strncpy(some_string, params->some_string, 100); // CodeParameters_c_files::some_string
    <BLANKLINE>
    >>> print((project_dir / 'set_CodeParameters-nopointer.h').read_text())
    const REAL a = commondata.a;                               // CodeParameters_c_files::a
    const bool BHaH_is_amazing = params.BHaH_is_amazing;       // CodeParameters_c_files::BHaH_is_amazing
    const REAL pi_three_sigfigs = commondata.pi_three_sigfigs; // CodeParameters_c_files::pi_three_sigfigs
    char some_string[100];
    strncpy(some_string, params.some_string, 100); // CodeParameters_c_files::some_string
    <BLANKLINE>
    >>> print((project_dir / 'set_CodeParameters-simd.h').read_text())
    const REAL NOSIMDa = commondata->a;                                         // CodeParameters_c_files::a
    const REAL_SIMD_ARRAY a = ConstSIMD(NOSIMDa);                               // CodeParameters_c_files::a
    const bool BHaH_is_amazing = params->BHaH_is_amazing;                       // CodeParameters_c_files::BHaH_is_amazing
    const REAL NOSIMDpi_three_sigfigs = commondata->pi_three_sigfigs;           // CodeParameters_c_files::pi_three_sigfigs
    const REAL_SIMD_ARRAY pi_three_sigfigs = ConstSIMD(NOSIMDpi_three_sigfigs); // CodeParameters_c_files::pi_three_sigfigs
    <BLANKLINE>
    """
    # Create output directory if it doesn't already exist
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)

    # Step 4: Generate C code to set C parameter constants
    #         output to filename "set_CodeParameters.h" if enable_simd==False
    #         or "set_CodeParameters-simd.h" if enable_simd==True

    # Step 4.a: Output non-SIMD version, set_CodeParameters.h
    def gen_set_CodeParameters(pointerEnable: bool = True) -> str:
        returnstring = ""
        for CPname, CodeParam in sorted(
            par.glb_code_params_dict.items(), key=lambda x: x[0].lower()
        ):
            if CodeParam.add_to_set_CodeParameters_h:
                if CodeParam.commondata:
                    struct = "commondata"
                else:
                    struct = "params"
                # C parameter type, parameter name
                CPtype = CodeParam.c_type_alias
                # For efficiency reasons, set_CodeParameters*.h does not set char arrays;
                #   access those from the params struct directly.
                pointer = "->" if pointerEnable else "."

                comment = f"  // {CodeParam.module}::{CPname}"
                if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                    # Handle char array C type
                    CPsize = CPtype.split("[")[1].split("]")[0]
                    Coutput = f"char {CPname}[{CPsize}];  strncpy({CPname}, {struct}{pointer}{CPname}, {CPsize});{comment}\n"
                else:
                    # Handle all other C types
                    Coutput = f"const {CPtype} {CPname} = {struct}{pointer}{CPname};{comment}\n"

                returnstring += Coutput

        return returnstring

    # Next, output header files for setting C parameters to current values within functions.
    header_file_path = project_Path / "set_CodeParameters.h"
    with header_file_path.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(
                gen_set_CodeParameters(pointerEnable=True),
                clang_format_options=clang_format_options,
            )
        )

    header_file_nopointer_path = project_Path / "set_CodeParameters-nopointer.h"
    with header_file_nopointer_path.open("w", encoding="utf-8") as file:
        file.write(
            clang_format(
                gen_set_CodeParameters(pointerEnable=False),
                clang_format_options=clang_format_options,
            )
        )

    # Step 4.b: Output SIMD version, set_CodeParameters-simd.h
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
            if CPtype == "REAL":
                c_output = f"const REAL NOSIMD{CPname} = {struct}->{CPname};{comment}\n"
                c_output += f"const REAL_SIMD_ARRAY {CPname} = ConstSIMD(NOSIMD{CPname});{comment}\n"
                set_CodeParameters_SIMD_str += c_output
            else:
                c_output = f"const {CPtype} {CPname} = {struct}->{CPname};{comment}\n"
                set_CodeParameters_SIMD_str += c_output

    header_file_simd_path = project_Path / "set_CodeParameters-simd.h"
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

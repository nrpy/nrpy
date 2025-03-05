"""
Provide functions to set params_struct and commondata_struct parameters to their default values specified when registering them within NRPy+'s CodeParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from pathlib import Path
from typing import List

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.helpers.generic import clang_format


def register_CFunctions_params_commondata_struct_set_to_default() -> None:
    """
    Register a C function to set code parameters to their default values.

    :raises ValueError: If an invalid default value is provided for any parameter. This
        ensures that all parameters can be correctly initialized in the
        generated C code.

    DocTests:
    >>> _real_array = par.register_CodeParameter("REAL[5]", "CodeParameters_c_files", "bad_real_array", 0.0, commondata=True, add_to_set_CodeParameters_h=True)
    Traceback (most recent call last):
    ...
    ValueError: Parameter 'bad_real_array' of type 'REAL[5]': For REAL or int array parameters, commondata must be True, and add_to_set_CodeParameters_h must be False.
    >>> _int_array = par.register_CodeParameter("int[3]", "CodeParameters_c_files", "bad_int_array", 42, commondata=False, add_to_set_CodeParameters_h=False)
    Traceback (most recent call last):
    ...
    ValueError: Parameter 'bad_int_array' of type 'int[3]': For REAL or int array parameters, commondata must be True, and add_to_set_CodeParameters_h must be False.
    >>> _, __ = par.register_CodeParameters("REAL", "CodeParameters_c_files", ["a", "pi_three_sigfigs"], [1, 3.14], commondata=True)
    >>> ___ = par.register_CodeParameter("#define", "CodeParameters_c_files", "b", 0)
    >>> _leaveitbe = par.register_CodeParameter("REAL", "CodeParameters_c_files", "leaveitbe", add_to_parfile=False, add_to_set_CodeParameters_h=False)
    >>> _int = par.register_CodeParameter("int", "CodeParameters_c_files", "blah_int", 1, commondata=True, add_to_parfile=True, add_to_set_CodeParameters_h=False)
    >>> _str = par.register_CodeParameter("char[100]", "CodeParameters_c_files", "some_string", "cheese")
    >>> _bool = par.register_CodeParameter("bool", "CodeParameters_c_files", "BHaH_is_amazing", True, add_to_set_CodeParameters_h=True)
    >>> _real_array = par.register_CodeParameter("REAL[5]", "CodeParameters_c_files", "real_array", 0.0, commondata=True, add_to_set_CodeParameters_h=False)
    >>> _int_array = par.register_CodeParameter("int[2]", "CodeParameters_c_files", "int_array", [4, 2], commondata=True, add_to_set_CodeParameters_h=False)
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunctions_params_commondata_struct_set_to_default()
    >>> print(cfc.CFunction_dict["params_struct_set_to_default"].full_function)
    #include "BHaH_defines.h"
    /**
     * Set params_struct to default values specified within NRPy+.
     */
    void params_struct_set_to_default(commondata_struct *restrict commondata, griddata_struct *restrict griddata) {
      // Loop over params structs:
      for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
        params_struct *restrict params = &griddata[grid].params;
        // Set params_struct variables to default
        params->BHaH_is_amazing = true;               // CodeParameters_c_files::BHaH_is_amazing
        snprintf(params->some_string, 100, "cheese"); // CodeParameters_c_files::some_string
      } // END LOOP over grids
    } // END FUNCTION params_struct_set_to_default
    <BLANKLINE>
    >>> print(cfc.CFunction_dict["commondata_struct_set_to_default"].full_function)
    #include "BHaH_defines.h"
    /**
     * Set commondata_struct to default values specified within NRPy+.
     */
    void commondata_struct_set_to_default(commondata_struct *restrict commondata) {
    <BLANKLINE>
      // Set commondata_struct variables to default
      commondata->a = 1;                   // CodeParameters_c_files::a
      commondata->blah_int = 1;            // CodeParameters_c_files::blah_int
      commondata->pi_three_sigfigs = 3.14; // CodeParameters_c_files::pi_three_sigfigs
      {
        REAL temp_val_array[] = {0.0, 0.0, 0.0, 0.0, 0.0};
        memcpy(commondata->real_array, temp_val_array, sizeof(temp_val_array));
      } // CodeParameters_c_files::real_array
      {
        int temp_val_array[] = {4, 2};
        memcpy(commondata->int_array, temp_val_array, sizeof(temp_val_array));
      } // CodeParameters_c_files::int_array
    } // END FUNCTION commondata_struct_set_to_default
    <BLANKLINE>
    """
    for function_name in ["commondata_struct", "params_struct"]:
        includes = ["BHaH_defines.h"]
        desc = f"Set {function_name} to default values specified within NRPy+."
        cfunc_type = "void"
        name = f"{function_name}_set_to_default"
        params = "commondata_struct *restrict commondata"
        if function_name == "params_struct":
            params += ", griddata_struct *restrict griddata"

        struct_list: List[str] = []  # List to store individual struct elements

        for parname, CodeParam in par.glb_code_params_dict.items():
            if (function_name == "commondata_struct" and CodeParam.commondata) or (
                function_name == "params_struct" and not CodeParam.commondata
            ):
                if (
                    CodeParam.defaultvalue != "unset"
                    and CodeParam.cparam_type != "#define"
                ):
                    struct = "commondata" if CodeParam.commondata else "params"
                    CPtype = CodeParam.cparam_type
                    comment = f"  // {CodeParam.module}::{parname}"
                    defaultval = CodeParam.defaultvalue
                    if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                        chararray_size = CPtype.split("[")[1].replace("]", "")
                        c_output = f'snprintf({struct}->{parname}, {chararray_size}, "{defaultval}");{comment}\n'
                    elif "[" in CPtype and "]" in CPtype:
                        # Handle REAL[N] and int[N] arrays
                        c_output = "{\n"
                        c_output += f"  {CPtype.split('[')[0]} temp_val_array[] = {{ {', '.join(str(x) for x in defaultval)} }};\n"
                        c_output += f"  memcpy({struct}->{parname}, temp_val_array, sizeof(temp_val_array));\n"
                        c_output += f"}} {comment}\n"
                    elif isinstance(defaultval, (bool, int, float)) or (
                        CPtype == "REAL" and isinstance(defaultval, str)
                    ):
                        c_output = f"{struct}->{parname} = {str(defaultval).lower()};{comment}\n"
                    else:
                        raise ValueError(
                            f"{CodeParam.defaultvalue} is not a valid default value type ({type(CodeParam.defaultvalue)}), for parameter {CodeParam.module}::{parname}, commondata = {CodeParam.commondata}"
                        )
                    struct_list.append(c_output)

        # Sort the lines alphabetically and join them with line breaks
        body = ""
        if function_name == "params_struct":
            body += r"""// Loop over params structs:
  for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
    params_struct *restrict params = &griddata[grid].params;
    // Set params_struct variables to default
"""
            body += "".join(sorted(struct_list))
            body += "} // END LOOP over grids\n"
        else:
            body += "\n// Set commondata_struct variables to default\n"
            body += "".join(sorted(struct_list))

        cfc.register_CFunction(
            include_CodeParameters_h=False,
            includes=includes,
            desc=desc,
            cfunc_type=cfunc_type,
            name=name,
            params=params,
            body=body,
        )


def write_CodeParameters_h_files(
    project_dir: str,
    set_commondata_only: bool = False,
    clang_format_options: str = "-style={BasedOnStyle: LLVM, ColumnLimit: 150}",
) -> None:
    r"""
    Generate C code to set C parameter constants and write them to files.

    :param project_dir: The path of the project directory.
    :param set_commondata_only: If True, generate code parameters only if `commondata=True`.
        Useful for BHaH projects without grids, like SEOBNR.
    :param clang_format_options: Options for clang_format.

    DocTests:
    >>> project_dir = Path("/tmp/tmp_project/")
    >>> write_CodeParameters_h_files(str(project_dir))
    >>> print((project_dir / 'set_CodeParameters.h').read_text())
    MAYBE_UNUSED const REAL a = commondata->a;                               // CodeParameters_c_files::a
    MAYBE_UNUSED const bool BHaH_is_amazing = params->BHaH_is_amazing;       // CodeParameters_c_files::BHaH_is_amazing
    MAYBE_UNUSED const REAL pi_three_sigfigs = commondata->pi_three_sigfigs; // CodeParameters_c_files::pi_three_sigfigs
    char some_string[100];                                                   // CodeParameters_c_files::some_string
    {
      // Safely copy string with snprintf, which guarantees null termination
      snprintf(some_string, sizeof(some_string), "%s", params->some_string);
    }
    >>> print((project_dir / 'set_CodeParameters-nopointer.h').read_text())
    MAYBE_UNUSED const REAL a = commondata.a;                               // CodeParameters_c_files::a
    MAYBE_UNUSED const bool BHaH_is_amazing = params.BHaH_is_amazing;       // CodeParameters_c_files::BHaH_is_amazing
    MAYBE_UNUSED const REAL pi_three_sigfigs = commondata.pi_three_sigfigs; // CodeParameters_c_files::pi_three_sigfigs
    char some_string[100];                                                  // CodeParameters_c_files::some_string
    {
      // Safely copy string with snprintf, which guarantees null termination
      snprintf(some_string, sizeof(some_string), "%s", params.some_string);
    }
    >>> print((project_dir / 'set_CodeParameters-simd.h').read_text())
    const REAL NOSIMDa = commondata->a;                                                      // CodeParameters_c_files::a
    MAYBE_UNUSED const REAL_SIMD_ARRAY a = ConstSIMD(NOSIMDa);                               // CodeParameters_c_files::a
    MAYBE_UNUSED const bool BHaH_is_amazing = params->BHaH_is_amazing;                       // CodeParameters_c_files::BHaH_is_amazing
    const REAL NOSIMDpi_three_sigfigs = commondata->pi_three_sigfigs;                        // CodeParameters_c_files::pi_three_sigfigs
    MAYBE_UNUSED const REAL_SIMD_ARRAY pi_three_sigfigs = ConstSIMD(NOSIMDpi_three_sigfigs); // CodeParameters_c_files::pi_three_sigfigs
    <BLANKLINE>
    """
    # Create output directory if it doesn't already exist
    project_Path = Path(project_dir)
    project_Path.mkdir(parents=True, exist_ok=True)
    # Generate C code to set C parameter constants
    # output to filename "set_CodeParameters.h" if enable_simd==False
    # or "set_CodeParameters-simd.h" if enable_simd==True

    # Output non-SIMD version, set_CodeParameters.h
    def gen_set_CodeParameters(pointerEnable: bool = True) -> str:
        """
        Generate content for set_CodeParameters*.h based on the pointerEnable flag.

        :param pointerEnable: A boolean flag indicating whether to access parameters through pointers.
            If True, parameters are accessed through pointers (struct->param).
            If False, direct access is assumed (struct.param).

        :return: A string containing the C code to be included in set_CodeParameters*.h, setting simulation parameters according to their specification in NRPy+.
        """
        returnstring = ""
        for CPname, CodeParam in sorted(
            par.glb_code_params_dict.items(), key=lambda x: x[0].lower()
        ):
            if CodeParam.add_to_set_CodeParameters_h:
                struct = "commondata" if CodeParam.commondata else "params"
                if not (set_commondata_only and struct == "params"):
                    CPtype = CodeParam.cparam_type
                    pointer = "->" if pointerEnable else "."
                    comment = f"  // {CodeParam.module}::{CPname}"
                    if "char" in CPtype and "[" in CPtype and "]" in CPtype:
                        # Handle char array C type
                        CPsize = int(CPtype.split("[")[1].split("]")[0])
                        # Char arrays are never unused; we use them below.
                        Coutput = rf"""char {CPname}[{CPsize}]; {comment}
{{
  // Safely copy string with snprintf, which guarantees null termination
  snprintf({CPname}, sizeof({CPname}), "%s", {struct}{pointer}{CPname});
}}"""
                    elif "[" in CPtype and "]" in CPtype:
                        # Handle REAL[N] and int[N] arrays
                        base_type = CPtype.split("[")[0]
                        array_size = CPtype.split("[")[1].split("]")[0]
                        Coutput = rf"""{base_type} {CPname}[{array_size}];{comment}
{{
  // Copy {array_size} elements from {struct}{pointer}{CPname} to {CPname}
  for(int i=0; i<{array_size}; i++) {{
    {CPname}[i] = {struct}{pointer}{CPname}[i];
  }}
}}"""
                    else:
                        # Handle all other C types.
                        # MAYBE_UNUSED is used to avoid compiler warnings about unused variables from including set_CodeParameters_simd.h.
                        Coutput = f"MAYBE_UNUSED const {CPtype} {CPname} = {struct}{pointer}{CPname};{comment}\n"

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
        # SIMD does not support arrays (char, REAL, int, etc.)
        if (
            CodeParam.add_to_set_CodeParameters_h
            and "char" not in CodeParam.cparam_type
            and "[" not in CodeParam.cparam_type
        ):
            struct = "commondata" if CodeParam.commondata else "params"
            if not (set_commondata_only and struct == "params"):
                CPtype = CodeParam.cparam_type
                comment = f"  // {CodeParam.module}::{CPname}"
                if CPtype == "REAL":
                    c_output = (
                        f"const REAL NOSIMD{CPname} = {struct}->{CPname};{comment}\n"
                    )
                    c_output += f"MAYBE_UNUSED const REAL_SIMD_ARRAY {CPname} = ConstSIMD(NOSIMD{CPname});{comment}\n"
                    set_CodeParameters_SIMD_str += c_output
                else:
                    # MAYBE_UNUSED is used to avoid compiler warnings about unused variables from including set_CodeParameters_simd.h.
                    c_output = f"MAYBE_UNUSED const {CPtype} {CPname} = {struct}->{CPname};{comment}\n"
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

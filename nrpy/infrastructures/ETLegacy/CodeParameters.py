"""
Functions for setting params_struct and commondata_struct parameters.

Set to default values specified when registering them within NRPy+'s CodeParameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import Optional, List, Tuple

import nrpy.params as par


def read_CodeParameters(
    list_of_tuples__thorn_CodeParameter: Optional[List[Tuple[str, str]]] = None,
    enable_simd: bool = False,
    declare_invdxxs: bool = True,
) -> str:
    r"""
    Generate C code to set C parameter constants and writes them to files.

    :param list_of_tuples__thorn_CodeParameter: A list of tuples containing thorn and code parameter names.
    :param enable_simd: Flag to enable SIMD for the code parameter.
    :param declare_invdxxs: Flag to declare inverse delta x.

    :return: A string that contains the generated C code.

    Doctests:
    >>> _, __ = par.register_CodeParameters(c_type_alias="CCTK_REAL", module="ignore", names=["a", "b"], defaultvalues=0.125)
    >>> outstr = read_CodeParameters(list_of_tuples__thorn_CodeParameter=[("thorna", "a"), ("thornb", "b")],
    ...                              enable_simd=True, declare_invdxxs=True)
    >>> print(outstr)
    const CCTK_REAL *restrict NOSIMDa = CCTK_ParameterGet("a", "thorna", NULL);  // thorna::a
    const REAL_SIMD_ARRAY a = ConstSIMD(*NOSIMDa);  // thorna::a
    const CCTK_REAL *restrict NOSIMDb = CCTK_ParameterGet("b", "thornb", NULL);  // thornb::b
    const REAL_SIMD_ARRAY b = ConstSIMD(*NOSIMDb);  // thornb::b
    const CCTK_REAL NOSIMDinvdxx0 = 1.0/CCTK_DELTA_SPACE(0);
    const REAL_SIMD_ARRAY invdxx0 = ConstSIMD(NOSIMDinvdxx0);
    const CCTK_REAL NOSIMDinvdxx1 = 1.0/CCTK_DELTA_SPACE(1);
    const REAL_SIMD_ARRAY invdxx1 = ConstSIMD(NOSIMDinvdxx1);
    const CCTK_REAL NOSIMDinvdxx2 = 1.0/CCTK_DELTA_SPACE(2);
    const REAL_SIMD_ARRAY invdxx2 = ConstSIMD(NOSIMDinvdxx2);
    <BLANKLINE>
    >>> outstr = read_CodeParameters(list_of_tuples__thorn_CodeParameter=[("thorna", "a"), ("thornb", "b")],
    ...                              enable_simd=False, declare_invdxxs=True)
    >>> print(outstr)
    const CCTK_REAL a = CCTK_ParameterGet("a", "thorna", NULL);  // thorna::a
    const CCTK_REAL b = CCTK_ParameterGet("b", "thornb", NULL);  // thornb::b
    const CCTK_REAL invdxx0 = 1.0/CCTK_DELTA_SPACE(0);
    const CCTK_REAL invdxx1 = 1.0/CCTK_DELTA_SPACE(1);
    const CCTK_REAL invdxx2 = 1.0/CCTK_DELTA_SPACE(2);
    <BLANKLINE>
    >>> outstr = read_CodeParameters(list_of_tuples__thorn_CodeParameter=[("thorna", "a"), ("thornb", "b")],
    ...                              enable_simd=False, declare_invdxxs=False)
    >>> print(outstr)
    const CCTK_REAL a = CCTK_ParameterGet("a", "thorna", NULL);  // thorna::a
    const CCTK_REAL b = CCTK_ParameterGet("b", "thornb", NULL);  // thornb::b
    <BLANKLINE>
    """

    def read_CCTK_REAL_CodeParameter(CPname: str, CPthorn: str) -> str:
        """
        Reads the CCTK_REAL code parameter based on the given CPname and CPthorn.

        :param CPname: The name of the code parameter.
        :param CPthorn: The thorn to which the code parameter belongs.

        :return: A string that represents the code for reading the CCTK_REAL parameter.

        >>> read_CCTK_REAL_CodeParameter("invdxx1", "thorn_name")
        'const CCTK_REAL invdxx1 = 1.0/CCTK_DELTA_SPACE(1);\n'
        """
        read_str = ""
        # fmt: off
        if CPname.startswith("invdxx"):
            dirn = CPname[-1]
            if enable_simd:
                read_str += f"const CCTK_REAL NOSIMD{CPname} = 1.0/CCTK_DELTA_SPACE({dirn});\n"
                read_str += f"const REAL_SIMD_ARRAY {CPname} = ConstSIMD(NOSIMD{CPname});\n"
            else:
                read_str += f"const CCTK_REAL {CPname} = 1.0/CCTK_DELTA_SPACE({dirn});\n"
        else:
            CPcomment = f"  // {CPthorn}::{CPname}"
            if enable_simd:
                read_str += f"""const CCTK_REAL *restrict NOSIMD{CPname} = CCTK_ParameterGet("{CPname}", "{CPthorn}", NULL);{CPcomment}\n"""
                read_str += f"const REAL_SIMD_ARRAY {CPname} = ConstSIMD(*NOSIMD{CPname});{CPcomment}\n"
            else:
                read_str += f"""const CCTK_REAL {CPname} = CCTK_ParameterGet("{CPname}", "{CPthorn}", NULL);{CPcomment}\n"""
        # fmt: on

        return read_str

    read_CodeParameters = ""
    if list_of_tuples__thorn_CodeParameter:
        for CPthorn, CPname in sorted(list_of_tuples__thorn_CodeParameter):
            CodeParam = par.glb_code_params_dict[CPname]
            if "char" in CodeParam.c_type_alias:
                raise ValueError("Cannot declare a char array in SIMD.")
            CPtype = CodeParam.c_type_alias
            if CPtype == "CCTK_REAL" or CPtype == "REAL":
                read_CodeParameters += read_CCTK_REAL_CodeParameter(CPname, CPthorn)
            else:
                CPcomment = f"  // {CPthorn}::{CPname}"
                c_output = f'const {CPtype} {CPname} = CCTK_ParameterGet("{CPname}", "{CodeParam.module}", NULL);{CPcomment}\n'
                read_CodeParameters += c_output
    if declare_invdxxs:
        for dirn in range(3):
            read_CodeParameters += read_CCTK_REAL_CodeParameter(f"invdxx{dirn}", "")

    return read_CodeParameters


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

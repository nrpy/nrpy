"""
Register the helper that converts Cartesian-basis BSSN fields to ADM fields.

This file emits only ``BSSN_to_ADM_Cartesian``. The generated function reads
Cartesian-basis BSSN gridfunctions from memory, computes Cartesian-basis ADM
quantities through the NRPy equation modules, and writes them in a stable
requested slot order.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from __future__ import annotations

from typing import List, Tuple

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.infrastructures import BHaH


def register_CFunction_BSSN_to_ADM_Cartesian() -> None:
    r"""
    Register ``BSSN_to_ADM_Cartesian``.

    The output array uses lexicographically sorted public names:
    ``ADMALPHA``, ``ADMBETA*``, ``ADMG**``, and ``ADMK**``.

    Doctests:
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.helpers.generic import validate_strings
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> register_CFunction_BSSN_to_ADM_Cartesian()  # doctest: +ELLIPSIS
    Setting up BSSN_Quantities[Cartesian]...
    >>> generated_str = cfc.CFunction_dict["BSSN_to_ADM_Cartesian"].full_function
    >>> validation_desc = "BSSN_to_ADM_Cartesian__openmp"
    >>> validate_strings(generated_str, validation_desc, file_ext="c")
    """
    bssn_to_adm = BSSN_to_ADM(CoordSystem="Cartesian", enable_rfm_precompute=False)
    bssn_quantities = BSSN_quantities["Cartesian"]

    labels = ["X", "Y", "Z"]
    output_pairs: List[Tuple[str, sp.Expr]] = [("ADMALPHA", bssn_quantities.alpha)]
    for i in range(3):
        output_pairs.append((f"ADMBETA{labels[i]}", bssn_to_adm.betaU[i]))
    for i in range(3):
        for j in range(i, 3):
            output_pairs.append(
                (f"ADMG{labels[i]}{labels[j]}", bssn_to_adm.gammaDD[i][j])
            )
            output_pairs.append((f"ADMK{labels[i]}{labels[j]}", bssn_to_adm.KDD[i][j]))
    output_pairs = sorted(output_pairs, key=lambda pair: pair[0])

    prefunc = "\n".join(
        [f"#define {name} {idx}" for idx, (name, _) in enumerate(output_pairs)]
    )
    desc = r"""
@brief Convert Cartesian-basis BSSN gridfunctions to Cartesian-basis ADM gridfunctions.

The input is read from ``in_gfs`` and the ADM output is written to ``out_gfs``
at matching grid points using the requested ADM slot names.

@param[in] params Grid parameters used for loop extents.
@param[in] in_gfs Cartesian-basis BSSN input gridfunctions.
@param[out] out_gfs Cartesian-basis ADM output gridfunctions.
"""
    cfunc_type = "void"
    name = "BSSN_to_ADM_Cartesian"
    params = (
        "const params_struct *restrict params, "
        "const REAL *restrict in_gfs, REAL *restrict out_gfs"
    )

    loop_body = ""
    for scalar in ("alpha", "cf", "trK"):
        loop_body += (
            f"const REAL {scalar} = in_gfs[IDX4({scalar.upper()}GF, i0, i1, i2)];\n"
        )
    for i in range(3):
        loop_body += f"const REAL vetU{i} = in_gfs[IDX4(VETU{i}GF, i0, i1, i2)];\n"
    for gf in ("hDD", "aDD"):
        for i in range(3):
            for j in range(i, 3):
                loop_body += f"const REAL {gf}{i}{j} = in_gfs[IDX4({gf.upper()}{i}{j}GF, i0, i1, i2)];\n"

    loop_body += ccg.c_codegen(
        [expr for _, expr in output_pairs],
        [f"out_gfs[IDX4({name}, i0, i1, i2)]" for name, _ in output_pairs],
        include_braces=False,
        verbose=False,
    )

    body = """const int Nxx_plus_2NGHOSTS0 = params->Nxx_plus_2NGHOSTS0;
const int Nxx_plus_2NGHOSTS1 = params->Nxx_plus_2NGHOSTS1;
const int Nxx_plus_2NGHOSTS2 = params->Nxx_plus_2NGHOSTS2;

"""
    body += BHaH.simple_loop.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        enable_intrinsics=False,
        read_xxs=False,
    )

    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
        include_CodeParameters_h=False,
    )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")

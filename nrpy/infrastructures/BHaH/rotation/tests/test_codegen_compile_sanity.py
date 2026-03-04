"""Compile-sanity check integration test for SO(3) rotation codegen helpers."""

from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.params as par
from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis_from_axis_angle import (
    register_CFunction_rotate_BSSN_Cartesian_basis_from_axis_angle,
)
from nrpy.infrastructures.BHaH.general_relativity.rotation.rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src import (
    register_CFunction_rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src,
)
from nrpy.infrastructures.BHaH.Makefile_helpers import (
    output_CFunctions_function_prototypes_and_construct_Makefile,
)
from nrpy.infrastructures.BHaH.rotation.register_all import register_CFunctions


def run_codegen_compile_sanity() -> bool:
    r"""
    Generate and compile SO(3) rotation C helpers in a temporary project.

    Doctests:
    >>> run_codegen_compile_sanity()
    True
    """
    cc = shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        # Toolchain unavailable in this environment; treat as skipped-success.
        return True

    outdir = Path("/tmp/nrpy_SO3_codegen_compile_sanity")
    if outdir.exists():
        shutil.rmtree(outdir)

    cfc.CFunction_dict.clear()
    par.set_parval_from_str("parallelization", "openmp")

    register_CFunctions()
    register_CFunction_rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src()
    register_CFunction_rotate_BSSN_Cartesian_basis_from_axis_angle()

    output_CFunctions_function_prototypes_and_construct_Makefile(
        str(outdir),
        "SO3_codegen_sanity",
        create_lib=True,
        use_openmp=False,
        addl_CFLAGS=["-Wall", "-Wextra", "-Werror"],
        compiler_opt_option="debug",
    )

    (outdir / "BHaH_defines.h").write_text(
        "\n".join(
            [
                "#ifndef BHAH_DEFINES_H",
                "#define BHAH_DEFINES_H",
                "#include <math.h>",
                "#include <stdio.h>",
                "#include <stdlib.h>",
                "typedef double REAL;",
                "typedef struct {",
                "  REAL cumulative_regrid_xhatU[3];",
                "  REAL cumulative_regrid_yhatU[3];",
                "  REAL cumulative_regrid_zhatU[3];",
                "} commondata_struct;",
                "#endif",
                "",
            ]
        ),
        encoding="utf-8",
    )

    prototypes = (outdir / "BHaH_function_prototypes.h").read_text(encoding="utf-8")
    required_prototypes = (
        "unrotate_xCart_to_fixed_frame",
        "unrotate_find_two_nUs_and_dphis_to_return_to_fixed_frame",
        "rotate_BSSN_Cartesian_basis_from_DeltaR_dst_from_src",
        "rotate_BSSN_Cartesian_basis_from_axis_angle",
    )
    for symbol in required_prototypes:
        if symbol not in prototypes:
            raise AssertionError(
                f"Missing expected prototype in sanity check build: {symbol}"
            )

    for c_file in sorted(outdir.rglob("*.c")):
        cmd = [
            cc,
            "-std=c11",
            "-I",
            str(outdir),
            "-c",
            str(c_file),
            "-o",
            str(c_file.with_suffix(".o")),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=False)
        if result.returncode != 0:
            raise AssertionError(
                "Compile sanity check failed for "
                + str(c_file)
                + "\nSTDOUT:\n"
                + result.stdout
                + "\nSTDERR:\n"
                + result.stderr
            )

    return True


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

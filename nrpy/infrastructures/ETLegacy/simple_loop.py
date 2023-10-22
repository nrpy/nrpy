"""
Simple loop generation for use within the ETLegacy infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

import nrpy.helpers.loop as lp


def simple_loop(
    loop_body: str,
    enable_simd: bool = False,
    loop_region: str = "",
    enable_OpenMP: bool = True,
    OMP_custom_pragma: str = "",
    OMP_collapse: int = 1,
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param enable_simd: Enable SIMD support
    :param loop_region: Loop over all points on a numerical grid or just the interior
    :param enable_OpenMP: Enable loop parallelization using OpenMP
    :param OMP_custom_pragma: Enable loop parallelization using OpenMP with custom pragma
    :param OMP_collapse: Specifies the number of nested loops to collapse
    :return: Complete loop code, output as a string.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="all points")))
    #pragma omp parallel for
    for (int i2 = 0; i2 < cctk_lsh[0]; i2++) {
      for (int i1 = 0; i1 < cctk_lsh[1]; i1++) {
        for (int i0 = 0; i0 < cctk_lsh[0]; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = 0; i0 < cctk_lsh[0]; i0++)
      }   // END LOOP: for (int i1 = 0; i1 < cctk_lsh[1]; i1++)
    } // END LOOP: for (int i2 = 0; i2 < cctk_lsh[0]; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="interior", OMP_custom_pragma="#CUSTOM_OMP")))
    #CUSTOM_OMP
    for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2] - cctk_nghostzones[2]; i2++) {
      for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1] - cctk_nghostzones[1]; i1++) {
        for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0] - cctk_nghostzones[0]; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = cctk_nghostzones[0]; i0 < cctk_lsh[0]-cctk_nghostzones[0]; i0++)
      }   // END LOOP: for (int i1 = cctk_nghostzones[1]; i1 < cctk_lsh[1]-cctk_nghostzones[1]; i1++)
    } // END LOOP: for (int i2 = cctk_nghostzones[2]; i2 < cctk_lsh[2]-cctk_nghostzones[2]; i2++)
    <BLANKLINE>
    """
    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "":
        return loop_body
    if loop_region == "all points":
        i2i1i0_mins = ["0", "0", "0"]
        i2i1i0_maxs = ["cctk_lsh[0]", "cctk_lsh[1]", "cctk_lsh[0]"]
    # 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
    elif loop_region == "interior":
        i2i1i0_mins = [
            "cctk_nghostzones[2]",
            "cctk_nghostzones[1]",
            "cctk_nghostzones[0]",
        ]
        i2i1i0_maxs = [
            "cctk_lsh[2]-cctk_nghostzones[2]",
            "cctk_lsh[1]-cctk_nghostzones[1]",
            "cctk_lsh[0]-cctk_nghostzones[0]",
        ]
    else:
        raise ValueError(
            f'loop_region = {loop_region} unsupported. Choose "", "all points", or "interior"'
        )

    # 'DisableOpenMP': disable loop parallelization using OpenMP
    if enable_OpenMP or OMP_custom_pragma != "":
        if OMP_custom_pragma == "":
            pragma = "#pragma omp parallel for"
            if OMP_collapse > 1:
                pragma = f"#pragma omp parallel for collapse({OMP_collapse})"
        # 'OMP_custom_pragma': enable loop parallelization using OpenMP with custom pragma
        else:
            pragma = OMP_custom_pragma
    else:
        pragma = ""
    increment = ["1", "1", "simd_width"] if enable_simd else ["1", "1", "1"]

    return str(
        lp.loop(
            ["i2", "i1", "i0"],
            i2i1i0_mins,
            i2i1i0_maxs,
            increment,
            pragma=[pragma, "", ""],
            loop_body=loop_body,
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

"""
Simple loop generation for use within the CarpetX infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""
from typing import List


def simple_loop(
    loop_body: str,
    enable_simd: bool = False,
    loop_region: str = "",
    loop_centering: List[int] = [0, 0, 0],
    run_on_device: bool = True,
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param enable_simd: Enable SIMD support
    :param loop_region: Loop over all points on a numerical grid or just the interior
    :param enable_OpenMP: Enable loop parallelization using OpenMP
    :return: Complete loop code, output as a string.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="all points")))
    grid.loop_all_device<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      // <INTERIOR>
    }); // END LOOP: grid.loop_all_device<0, 0, 0>
    <BLANKLINE>
      >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="interior", run_on_device=False)))
    grid.loop_int<0, 0, 0>(
        grid.nghostzones,
        [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      // <INTERIOR>
    }); // END LOOP: grid.loop_int<0, 0, 0>
    <BLANKLINE>
    """
    loop_macro = "grid.loop_"
    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "all points":
        loop_macro += "all"
    # 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
    elif loop_region == "interior":
        loop_macro += "int"
    else:
        raise ValueError(
            f'loop_region = {loop_region} unsupported. Choose "all points" or "interior"'
        )

    if run_on_device:
        loop_macro += "_device"

    if loop_centering[0] not in [0, 1]:
        raise ValueError(
            f"loop_centering[0] = {loop_centering[0]} unsupported. Choose 0 or 1 for vertex- or cell-centering, respectively."
        )
    if loop_centering[1] not in [0, 1]:
        raise ValueError(
            f"loop_centering[1] = {loop_centering[1]} unsupported. Choose 0 or 1 for vertex- or cell-centering, respectively."
        )
    if loop_centering[2] not in [0, 1]:
        raise ValueError(
            f"loop_centering[2] = {loop_centering[2]} unsupported. Choose 0 or 1 for vertex- or cell-centering, respectively."
        )
    # loop_centering: set loop to properly choose stencil based on the given centering
    loop_macro += f"<{loop_centering[0]}, {loop_centering[1]}, {loop_centering[2]}>(\n"

    return str(
        loop_macro
        + "grid.nghostzones,\n"
        + "[=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {\n"
        + loop_body
        + "}); // END LOOP: "
        + loop_macro
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

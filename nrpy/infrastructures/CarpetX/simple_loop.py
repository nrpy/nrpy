"""
Simple loop generation for use within the CarpetX infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

from typing import Union, List


def simple_loop(
    loop_body: str,
    enable_simd: bool = False,
    loop_region: str = "",
    loop_centering: Union[List[int], None] = None,
    run_on_device: bool = True,
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: The body of the loop to be generated.
    :param enable_simd: Flag to enable SIMD support. Currently unsupported.
    :param loop_region: Specifies the loop region (e.g., 'all points', 'interior').
    :param loop_centering: Specifies the centering of the loop (0 = 'V', 1 = 'C').
    :param run_on_device: Flag to run loop kernel on device (GPU) instead of host (CPU).
    :return: The complete loop code as a string.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="all points")))
    grid.loop_all_device<0, 0, 0>(grid.nghostzones, [=] CCTK_DEVICE(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      // <INTERIOR>
    }); // END LOOP: grid.loop_all_device<0, 0, 0>(
    <BLANKLINE>
    >>> print(clang_format(simple_loop(loop_body='// <INTERIOR>', loop_region="interior", run_on_device=False)))
    grid.loop_int<0, 0, 0>(grid.nghostzones, [=] CCTK_HOST(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {
      // <INTERIOR>
    }); // END LOOP: grid.loop_int<0, 0, 0>(
    <BLANKLINE>
    """
    if loop_centering is None:
        loop_centering = [0, 0, 0]
    elif (
        not isinstance(loop_centering, list)
        or len(loop_centering) != 3
        or not all(isinstance(x, int) for x in loop_centering)
    ):
        raise ValueError(
            "loop_centering must be set to None (defaults to [0,0,0]) or be a list of three integers."
        )

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

    if enable_simd:
        raise ValueError(
            "CarpetX SIMD generation is not yet supported. Please generate with enable_simd=False."
        )

    if run_on_device:
        loop_macro += "_device"
        run_on = "DEVICE"
    else:
        run_on = "HOST"

    for i, center in enumerate(loop_centering):
        if center not in [0, 1]:
            raise ValueError(
                f"loop_centering[{i}] = {center} unsupported. Choose 0 or 1 for vertex- or cell-centering, respectively."
            )

    loop_macro += f"<{loop_centering[0]}, {loop_centering[1]}, {loop_centering[2]}>(\n"

    return f"""{loop_macro}grid.nghostzones,
[=] CCTK_{run_on}(const Loop::PointDesc &p) CCTK_ATTRIBUTE_ALWAYS_INLINE {{
{loop_body}
}}); // END LOOP: {loop_macro}
"""


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

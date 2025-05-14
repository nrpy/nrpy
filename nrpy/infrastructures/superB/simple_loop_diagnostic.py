"""
Simple loop generation for use within the superB infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
Contributor: Nishita Jadoo
Email: njadoo *at* uidaho *dot* edu

"""

from typing import List, Tuple, Union

import sympy as sp

import nrpy.indexedexp as ixp
from nrpy.infrastructures.BHaH.simple_loop import (
    append_1d_loop_body,
    compute_1d_loop_ranges,
    generate_1d_loop_header,
    generate_qsort_compare_string,
    max_numpts__i012_pts__numpts_2D,
)


def simple_loop_1D(
    CoordSystem: str,
    axis: str = "z",
) -> Tuple[str, str]:
    r"""
    Generate a simple 1D loop in C (for use inside of a function).

    :param CoordSystem: Coordinate system, e.g., "Cartesian".
    :param axis: Specifies the axis of output: either "x" or "z" axis.

    :return: Complete loop code, output as a string.

    :raises ValueError: If axis is not 'y' or 'z'.
    """

    NGHOSTS, Nxx_plus_2NGHOSTS, Nxx, i012_pts, numpts = compute_1d_loop_ranges(
        CoordSystem, axis
    )

    pragma = "#pragma omp parallel for\n"

    out_string = generate_1d_loop_header(axis, CoordSystem, numpts)

    # Loop body for storing results.
    loop_body_store_results = f"""{{
data_point_1d_struct dp1d;
dp1d.xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
dp1d.i0 = i0;
dp1d.i1 = i1;
dp1d.i2 = i2;
"""

    loop_body_store_results += "data_points[data_index] = dp1d; data_index++;\n}\n"

    out_string = append_1d_loop_body(
        out_string, loop_body_store_results, axis, Nxx, numpts, i012_pts, pragma
    )

    # Post-loop: qsort() along xCart_axis .
    prefunc_content = """
typedef struct {
REAL xCart_axis;
int i0;
int i1;
int i2;
} data_point_1d_struct;
"""
    prefunc_content += generate_qsort_compare_string()

    qsort = r"""
qsort(data_points, data_index, sizeof(data_point_1d_struct), compare);
"""
    out_string += qsort

    return prefunc_content, out_string


def simple_loop_2D(
    CoordSystem: str,
    plane: str = "yz",
) -> str:
    r"""
    Generate a simple 2D loop in the xy or yz plane in C (for use inside of a function).

    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param plane: Specifies the plane of output: either "xy" or "yz" (default)
    :return: Complete loop code, output as a string.
    :raises ValueError: If plane is not 'xy' or 'yz'.
    """
    pragma = "#pragma omp parallel for\n"
    max_numpts, i012_pts, numpts = max_numpts__i012_pts__numpts_2D(CoordSystem, plane)

    out_string = f"""// Define points for output along the {plane}-plane in {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
"""
    for i in range(3):
        if numpts[i] == max_numpts[i]:
            out_string += f"{pragma}for(int i{i}=NGHOSTS; i{i}<Nxx{i} + NGHOSTS; i{i}++) i{i}_pts[i{i}-NGHOSTS] = i{i};\n"
        else:
            for j, pt in enumerate(i012_pts[i]):
                out_string += f"i{i}_pts[{j}] = (int)({sp.ccode(pt)});\n"

    return out_string


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

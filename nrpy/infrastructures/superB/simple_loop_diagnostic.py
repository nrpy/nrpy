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
    if axis not in ["y", "z"]:
        raise ValueError(
            f"1D loop output only supports y or z axes. axis = {axis} not supported."
        )

    NGHOSTS = sp.Symbol("NGHOSTS", real=True)
    Nxx_plus_2NGHOSTS = ixp.declarerank1("Nxx_plus_2NGHOSTS")
    Nxx = ixp.declarerank1("Nxx")

    i0_pts: List[Union[int, sp.Expr]] = []
    i1_pts: List[Union[int, sp.Expr]] = []
    i2_pts: List[Union[int, sp.Expr]] = []

    if axis == "y":
        if "Cartesian" in CoordSystem:
            # x-axis == { x_mid, z_mid }
            i0_pts += [Nxx_plus_2NGHOSTS[0] / 2]
            i2_pts += [Nxx_plus_2NGHOSTS[2] / 2]
        elif "Spherical" in CoordSystem:
            # y-axis == { theta_mid, see yz-plane discussion for Spherical below for explanation of phi points }
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
            i2_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[2] - sp.Rational(1, 2)]
            i2_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[2] - sp.Rational(1, 2)]
        elif "Cylindrical" in CoordSystem:
            # Cylindrical: rho,phi,z
            # y-axis == { see yz-plane discussion for Spherical below for explanation of phi points, z_mid }
            i1_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[1] - sp.Rational(1, 2)]
            i1_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[1] - sp.Rational(1, 2)]
            i2_pts += [Nxx_plus_2NGHOSTS[2] / 2]
        elif "SymTP" in CoordSystem:
            # y-axis == { x1_mid, see yz-plane discussion for Spherical below for explanation of phi points }
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
            i2_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[2] - sp.Rational(1, 2)]
            i2_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[2] - sp.Rational(1, 2)]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    if axis == "z":
        if "Cartesian" in CoordSystem:
            # z-axis == { x_mid, y_mid }
            i0_pts += [Nxx_plus_2NGHOSTS[0] / 2]
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
        elif "Spherical" in CoordSystem:
            # z-axis == { th_min & th_max, phi_min  }
            i1_pts += [NGHOSTS]
            i1_pts += [Nxx_plus_2NGHOSTS[1] - NGHOSTS - 1]
            i2_pts += [NGHOSTS]
        elif "Cylindrical" in CoordSystem:
            # Cylindrical: rho,phi,z
            # z-axis == { rho_min & phi_min }
            i0_pts += [NGHOSTS]
            i1_pts += [NGHOSTS]
        elif "SymTP" in CoordSystem:
            # SymTP:
            # self.xx_to_Cart[2] = f(xx0) * sp.cos(self.xx[1])
            #  -> Aim for cos(xx1) = 1 -> xx1 = 0 & pi
            # z_axis == { xx1_min & xx1_max, xx2_min }
            # FIXME: Probably missing points between foci.
            i1_pts += [NGHOSTS]
            i1_pts += [Nxx_plus_2NGHOSTS[1] - NGHOSTS - 1]
            i2_pts += [NGHOSTS]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    i012_pts = [i0_pts, i1_pts, i2_pts]
    numpts: List[Union[int, sp.Expr]] = [0] * 3
    if (
        (i0_pts and i0_pts[0] == -1)
        or (i1_pts and i1_pts[0] == -1)
        or (i2_pts and i2_pts[0] == -1)
    ):
        numpts = [0, 0, 0]
    else:
        numpts[0] = len(i0_pts) if i0_pts else Nxx[0]
        numpts[1] = len(i1_pts) if i1_pts else Nxx[1]
        numpts[2] = len(i2_pts) if i2_pts else Nxx[2]

    pragma = "#pragma omp parallel for\n"
    out_string = f"""// Output along {axis}-axis in {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];

data_point_1d_struct data_points[numpts_i0 * numpts_i1 * numpts_i2];
int data_index = 0;
"""

    # Loop body for storing results.
    loop_body_store_results = f"""{{
data_point_1d_struct dp1d;
dp1d.xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
dp1d.i0 = i0;
dp1d.i1 = i1;
dp1d.i2 = i2;
"""

    # ~ for key, value in out_quantities_dict.items():
    # ~ loop_body_store_results += f"dp1d.{key[1]} = {value};\n"
    loop_body_store_results += "data_points[data_index] = dp1d; data_index++;\n}\n"

    # Main loop body.
    for i in range(3):
        if numpts[i] == Nxx[i]:
            out_string += f"{pragma}for(int i{i}=NGHOSTS; i{i}<Nxx{i} + NGHOSTS; i{i}++) i{i}_pts[i{i}-NGHOSTS] = i{i};\n"
        elif numpts == [0, 0, 0] and i == 0:
            out_string += (
                f"// CoordSystem = {CoordSystem} has no points on the {axis} axis!\n"
            )
        else:
            for j, pt in enumerate(i012_pts[i]):
                out_string += f"i{i}_pts[{j}] = (int)({sp.ccode(pt)});\n"
    out_string += f"""// Main loop:
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

  {loop_body_store_results}
}}
"""

    # Post-loop: qsort() along xCart_axis and output to file.
    prefunc_content = """
typedef struct {
REAL xCart_axis;
int i0;
int i1;
int i2;
"""
    prefunc_content += """} data_point_1d_struct;

// qsort() comparison function for 1D output.
static int compare(const void *a, const void *b) {
  REAL l = ((data_point_1d_struct *)a)->xCart_axis;
  REAL r = ((data_point_1d_struct *)b)->xCart_axis;
  return (l > r) - (l < r);
}
"""

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
    if plane not in ["xy", "yz"]:
        raise ValueError(
            f"2D loop output only supports xy or yz planes. plane = {plane} not supported."
        )

    NGHOSTS = sp.Symbol("NGHOSTS", real=True)
    Nxx_plus_2NGHOSTS = ixp.declarerank1("Nxx_plus_2NGHOSTS")
    Nxx = ixp.declarerank1("Nxx")

    i0_pts: List[sp.Expr] = []
    i1_pts: List[sp.Expr] = []
    i2_pts: List[sp.Expr] = []

    if plane == "xy":
        if "Cartesian" in CoordSystem or "Cylindrical" in CoordSystem:
            # xy-plane == { z_mid }, where z index is i2
            i2_pts += [Nxx_plus_2NGHOSTS[2] / 2]
        elif "Spherical" in CoordSystem or "SymTP" in CoordSystem:
            # xy-plane == { theta_mid }, where theta index is i1
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    if plane == "yz":
        if "Cartesian" in CoordSystem:
            i0_pts += [Nxx_plus_2NGHOSTS[0] / 2]
        elif "Cylindrical" in CoordSystem:
            # See documentation for Spherical below; Cylindrical-like coordinates choose xx1 = phi instead of xx2.
            i1_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[1] - sp.Rational(1, 2)]
            i1_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[1] - sp.Rational(1, 2)]
        elif "Spherical" in CoordSystem or "SymTP" in CoordSystem:
            # In Spherical/SymTP coordinates, phi_min=-PI and phi_max=+PI.
            # The yz plane is specifically at -PI/2 and +PI/2.
            # When Nphi=2, NGHOSTS and NGHOSTS+1 correspond to -PI/2 and +PI/2 *exactly*.
            # WARNING: For Nphi=4, the planes don't sample -PI/2 and +PI/2 exactly. Instead, they sample at:
            #          {-3/4, -1/4, +1/4, +3/4}*PI. The closest planes are at indices 1 and 3 for these values.
            #                   ^           ^  <--- closest planes; at 1 and 3
            #             ^           ^        <--- closest planes; at 0 and 2
            #          The same applies for Nphi=4, 8, 12, etc. It's best to choose Nphi as a multiple of 2 but not 4.
            # General formula for cell-centered grid in phi is:
            # xx_i = -PI + [(i-NGHOSTS) + 0.5] * (2*PI)/Nxx2
            # In symbolic math (using sympy), the closest planes for -PI/2 and +PI/2 can be calculated as follows:
            # index_at_PIo2, PI, Nxx2, NGHOSTS = symbols('index_at_PIo2 PI Nxx2 NGHOSTS', real=True)
            # expr = -PI + ((index_at_PIo2-NGHOSTS) + Rational(1,2)) * (2*PI)/Nxx2
            # solve(expr - PI/2, index_at_PIo2)  # solves for expr = PI/2
            # # NGHOSTS + 3*Nphi/4 - 1/2 -> Closest plane for +PI/2: NGHOSTS + 3*Nxx2/4 - 1/2
            # solve(expr - (-PI/2), index_at_PIo2)  # solves for expr = -PI/2
            # # NGHOSTS + Nphi/4 - 1/2 -> Closest plane for -PI/2: NGHOSTS + Nxx2/4 - 1/2
            i2_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[2] - sp.Rational(1, 2)]
            i2_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[2] - sp.Rational(1, 2)]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    i012_pts = [i0_pts, i1_pts, i2_pts]
    numpts: List[Union[int, sp.Expr]] = [0] * 3
    numpts[0] = len(i0_pts) if i0_pts else Nxx[0]
    numpts[1] = len(i1_pts) if i1_pts else Nxx[1]
    numpts[2] = len(i2_pts) if i2_pts else Nxx[2]
    pragma = "#pragma omp parallel for\n"
    out_string = f"""// Output data in {plane}-plane in {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
"""
    for i in range(3):
        if numpts[i] == Nxx[i]:
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

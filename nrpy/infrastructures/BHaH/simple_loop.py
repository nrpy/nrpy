"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

from typing import List, Union
import sympy as sp
import nrpy.helpers.loop as lp
import nrpy.indexedexp as ixp


def simple_loop(
    loop_body: str,
    enable_simd: bool = False,
    loop_region: str = "",
    read_xxs: bool = False,
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
    enable_OpenMP: bool = True,
    OMP_custom_pragma: str = "",
    OMP_collapse: int = 1,
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param enable_simd: Enable SIMD support
    :param loop_region: Loop over all points on a numerical grid or just the interior
    :param read_xxs: Read the xx[3][:] 1D coordinate arrays if interior dependency exists
    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param enable_rfm_precompute: Enable pre-computation of reference metric
    :param enable_OpenMP: Enable loop parallelization using OpenMP
    :param OMP_custom_pragma: Enable loop parallelization using OpenMP with custom pragma
    :param OMP_collapse: Specifies the number of nested loops to collapse
    :return: Complete loop code, output as a string.

    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="all points")))
    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
      }   // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior", OMP_custom_pragma="#CUSTOM_OMP")))
    #CUSTOM_OMP
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=3)))
    Setting up reference metric for CoordSystem = SinhSymTP.
    #pragma omp parallel for collapse(3)
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          const REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
          const REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
          const REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
          const REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
          const REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
          const REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
          const REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          const REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          const REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          const REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          const REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          const REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          const REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=2)))
    #pragma omp parallel for collapse(2)
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
        const REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
        const REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
        const REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
        const REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
        const REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
    <BLANKLINE>
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          const REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          const REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          const REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          const REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          const REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          const REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          const REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      }   // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    <BLANKLINE>
    """
    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "":
        return loop_body
    if loop_region == "all points":
        i2i1i0_mins = ["0", "0", "0"]
        i2i1i0_maxs = ["Nxx_plus_2NGHOSTS2", "Nxx_plus_2NGHOSTS1", "Nxx_plus_2NGHOSTS0"]
    # 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
    elif loop_region == "interior":
        i2i1i0_mins = ["NGHOSTS", "NGHOSTS", "NGHOSTS"]
        i2i1i0_maxs = ["NGHOSTS+Nxx2", "NGHOSTS+Nxx1", "NGHOSTS+Nxx0"]
    else:
        raise ValueError(
            f'loop_region = {loop_region} unsupported. Choose "", "all points", or "interior"'
        )

    read_rfm_xx_arrays = ["", "", ""]
    # 'Read_xxs': read the xx[3][:] 1D coordinate arrays, as some interior dependency exists
    if read_xxs:
        if not enable_simd:
            read_rfm_xx_arrays = [
                "const REAL xx0 = xx[0][i0];",
                "const REAL xx1 = xx[1][i1];",
                "const REAL xx2 = xx[2][i2];",
            ]
        else:
            raise ValueError("no innerSIMD support for Read_xxs (currently).")
    # 'enable_rfm_precompute': enable pre-computation of reference metric
    if enable_rfm_precompute:
        if read_xxs:
            raise ValueError(
                "enable_rfm_precompute and Read_xxs cannot both be enabled."
            )
        # pylint: disable=C0415
        from nrpy.infrastructures.BHaH import rfm_precompute

        rfmp = rfm_precompute.ReferenceMetricPrecompute(CoordSystem)
        if enable_simd:
            read_rfm_xx_arrays = [
                rfmp.readvr_SIMD_inner_str[0],
                rfmp.readvr_SIMD_outer_str[1],
                rfmp.readvr_SIMD_outer_str[2],
            ]
        else:
            read_rfm_xx_arrays = [
                rfmp.readvr_str[0],
                rfmp.readvr_str[1],
                rfmp.readvr_str[2],
            ]
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

    loop_body = read_rfm_xx_arrays[0] + loop_body
    prefix_loop_with = [pragma, read_rfm_xx_arrays[2], read_rfm_xx_arrays[1]]
    if OMP_collapse == 2:
        prefix_loop_with = [
            pragma,
            "",
            read_rfm_xx_arrays[2] + read_rfm_xx_arrays[1],
        ]
    elif OMP_collapse == 3:
        prefix_loop_with = [
            pragma,
            "",
            "",
        ]
        # above: loop_body = read_rfm_xx_arrays[0] + loop_body -----v
        loop_body = read_rfm_xx_arrays[2] + read_rfm_xx_arrays[1] + loop_body

    return str(
        lp.loop(
            ["i2", "i1", "i0"],
            i2i1i0_mins,
            i2i1i0_maxs,
            increment,
            prefix_loop_with,
            loop_body=loop_body,
        )
    )


def simple_loop_2D(
    loop_body: str,
    CoordSystem: str,
    plane: str = "yz",
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param plane: Specifies the plane of output: either "xy" or "yz" (default)
    :return: Complete loop code, output as a string.
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop_2D("// loop body", "Cartesian", plane="xy")))
    // Set 2D loops over xy-plane for Cartesian coordinates.
    const int numpts_i0 = Nxx0, numpts_i1 = Nxx1, numpts_i2 = 1;
    int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
    #pragma omp parallel for
    for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
      i0_pts[i0 - NGHOSTS] = i0;
    #pragma omp parallel for
    for (int i1 = NGHOSTS; i1 < Nxx1 + NGHOSTS; i1++)
      i1_pts[i1 - NGHOSTS] = i1;
    i2_pts[0] = (int)((1.0 / 2.0) * Nxx_plus_2NGHOSTS2);
    // Main loop:
    LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
      const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
      // loop body
    }
    >>> print(clang_format(simple_loop_2D("// loop body", "Spherical", plane="yz")))
    // Set 2D loops over yz-plane for Spherical coordinates.
    const int numpts_i0 = Nxx0, numpts_i1 = Nxx1, numpts_i2 = 2;
    int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
    #pragma omp parallel for
    for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
      i0_pts[i0 - NGHOSTS] = i0;
    #pragma omp parallel for
    for (int i1 = NGHOSTS; i1 < Nxx1 + NGHOSTS; i1++)
      i1_pts[i1 - NGHOSTS] = i1;
    i2_pts[0] = (int)(NGHOSTS + (1.0 / 4.0) * Nxx2 - 1.0 / 2.0);
    i2_pts[1] = (int)(NGHOSTS + (3.0 / 4.0) * Nxx2 - 1.0 / 2.0);
    // Main loop:
    LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
      const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
      // loop body
    }
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

    if plane == "yz":
        if CoordSystem == "Cartesian":
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

    i012_pts = [i0_pts, i1_pts, i2_pts]
    numpts: List[Union[int, sp.Expr]] = [
        0,
        0,
        0,
    ]
    numpts[0] = len(i0_pts) if i0_pts else Nxx[0]
    numpts[1] = len(i1_pts) if i1_pts else Nxx[1]
    numpts[2] = len(i2_pts) if i2_pts else Nxx[2]
    pragma = "#pragma omp parallel for\n"
    out_string = f"""// Set 2D loops over {plane}-plane for {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];
"""
    for i in range(3):
        if numpts[i] == Nxx[i]:
            out_string += f"{pragma}for(int i{i}=NGHOSTS; i{i}<Nxx{i} + NGHOSTS; i{i}++) i{i}_pts[i{i}-NGHOSTS] = i{i};\n"
        else:
            for j, pt in enumerate(i012_pts[i]):
                out_string += f"i{i}_pts[{j}] = (int)({sp.ccode(pt)});\n"
    out_string += f"""// Main loop:
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  {loop_body}
}}"""
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

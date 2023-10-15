"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

from typing import List, Union, Tuple, Dict
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

    Doctests:
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


def simple_loop_1D(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    axis: str = "z",
) -> Tuple[str, str]:
    r"""
    Generate a simple 1D loop in C (for use inside of a function).

    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param out_quantities_dict: Output quantities dictionary
    :param axis: Specifies the axis of output: either "x" or "z" axis
    :return: Complete loop code, output as a string.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4ASAAjddABfgfJ51OEu5G9+EUlvYp1dk7+K8eO40VXFOrC5KgQdAGZkA+nmwLa3aIodzb3cJ641+dM68R4YgftRQdXGatJLFn6FOce+L5rhgnsnMvkljHaochh2KqpK9sEbv1xn+MV3kddMhWMLGRwssQKqNsAeNMBHVziHlJECGFaTWWDnj0k0dlKF927gVSP7thdDf/zDIZ0bv0zRbyjdVGTg/KxnnF1Ad0wXSILsLCGvv2Q+mD6VTwU3yheYlrnYnZZeUsAUiWiGwXJ4en7FVF98GkBlyQ9/4BOYE6ITi3TA9cUIYyYJzXzosV9ABszAcYWklb30SwrPnEsGp4iwYTffIeX3iSf3V9y1LJhVEziPRR62Iw+Fg/YzIPke50zqQVLSO/YxHpgfRc9gHuEhDwE3b5LLUxZnRipI8NgCk7RY+1sxoabKNwrRjm4j32cSyRIDHHqwt9drJThsr3ik6Y2ljmkuyOT7sooYVVL8+npF7r4MPjP9p4XbQ4Wqd+AM/sqgnlfGw2MxoBrw8FN3L5WNZA4PYYbMwk/PLkNcIGMHvJavyOCqL1sEiRsmiBeTaFeG8084j7oJTtH1WVTeq5+bYNs1Denili0e6Hs6wv4CfL+bKj0v6/l4XI2qqxCOStvyh11EpUgL8H8tqGuKrFpWVvy3uYjHpBr6P9CWgWgR3M4720kHldCRgV/RVFahYlAfKOQAo5JGXTuzww4jQFBkJhhFBQtQNuItDBWKh+TFxUDTmOXMukPhkAAAAKGvTIYQgyTYAAdMEgQkAAE32S/ixxGf7AgAAAAAEWVo=")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4ATJAk5dABfgfJ51OEu5G9+EUlvYp1dk7+K8eO40VXFOwE0lywhL6Na6gbyqpgT1R8e50ML5tgBveOVyZl6U8ClAM+5RMuYtVE+1XueVLF0J9QDYPB/lrsfhIONJAKP/VtgkerMr8ejwimCGmGbFM8KQHLH2Qr6iPaXrqbVso37SpGwID3OAKugmkzy5qdleSgsULA4KPPQaDiyXM21X+kfxEmRnVFIAwwMH4R7KDHeTGEWgJ+8ZkEcu634G9HNinr10zZsG1ZyhJjqJXmI6SKwqSFkoUVkvmrB6n14wevAN//zJa8TdpUIgIMKkGlO10TjtsvWFOvDH4arjW3/2EgKRR1/QeN2YIRqaufz/z111Sul6q2z1mysFx7vNN1CiEZjIsnKQVd6bgB6D4gZcFnjyMl/vVGUoNssGHgbPA3zUJikPxrlY98Qb5aCclWep6Gp7hOwB4wUfvtQZ/TXEjF1Ue9XRAQ9vEt41eHUKxtjSewIvPf1U8lSS48132PyTeNPBOOI8+afNu5m7VzEnPYUI0oO+p0D5QDSXdcJVSgE0R7lWQb62t2peU0r7wVdcDLBOsQa7M7PhZzqJpTflY4+MnrTaQhBykeZD3ADQ9xTT8LiX3NvFcTmTQd/StqDRFfwW8qHrlwP3vM8a99SRy4JJ4v8/WFvpmwXQP++RNnK91WgUX8JJixfzivAh4IRtR+Y5c3a6EfhKxYoRg5mCnujk/dJTSs/rhsPzNpNx7Eb66gh3Z6MUNHhmhH1DWIHt+lDi4ZsVLrdu6HMakW7slMueKO9vAAAAnoB6kHdcIEoAAeoEygkAABR/K5exxGf7AgAAAAAEWVo=")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Spherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="z")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4ASEAjtdABfgfJ51OEu5G9+EUlvYp1dlAyuYeO40VXFOv/eTS7T39N1d1WAuEsg/wIVqdPkL2NmyE0Z0wib/s9J2rhoPSkyFAtdMiER/4Zf//bctEJLtW0auvD5Twse/NuGMT0C0Yy9gFVoSJfhxIi1c1w/cttPIWplpZd69Uybtrx1SYR1p31EyLpIeJVZqitJRNFJ0S5Rbxtz0WQ4aOJXV2WvgsCigjb0GdPZpdunqATvpCVL8BSkP2YXPuYuy55NnIOvuIfbGSew5g23Nsqf9SQK+UFJNx9K/fmkd9GBouOwKegr/AGP/yffxC23YC6zGlqHIT7r4CixsWKhYSK855u0xleA2Cspg12x8gVHCsuSKbxHoGH43W9OXYSQYN8Tg1yTgAvnDPVB4HLiOtKIAN+xvd+wApMIR6isVwqn24TozmgkCRTFTfM6XS0P8RPQwajr05OoyUYkdYTbhJOppY15vpA5wqrXuGCLy+GRdOR11OpXc1W7gg+tSWD562XzUXXBt9ryInnC4bqBTpYnZMfe8P1i8Iv6yQkHYXMKTCLeJVP208+aH1IamdDoxy++JttEjkpz9uDHlxwIRGlaY8aBft2XxxPR6a0VJ9GmIflyBf1ZfWoJRXOSHFTvgV48RrAkiWuDa03L4CdDW2JsMl+RqefYskbjfDVEtRdzR//bdEOJMso3AH+P3dOT3h8eceZqniq+PfjcM4gFiPwJ97yVU88Dn9TwZRE1HZ7wJ9Z780q+w2yBKEOqC5s343AsAAA2NhflX1595AAHXBIUJAAAMI7jsscRn+wIAAAAABFla")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"Diff:\n {diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
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
"""

    for key, value in out_quantities_dict.items():
        loop_body_store_results += f"dp1d.{key[1]} = {value};\n"
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
"""
    for key in out_quantities_dict.keys():
        prefunc_content += f"  {key[0]} {key[1]};\n"
    prefunc_content += """} data_point_1d_struct;

// qsort() comparison function for 1D output.
static int compare(const void *a, const void *b) {
  REAL l = ((data_point_1d_struct *)a)->xCart_axis;
  REAL r = ((data_point_1d_struct *)b)->xCart_axis;
  return (l > r) - (l < r);
}
"""

    qsort_and_output_to_file = r"""
qsort(data_points, data_index, sizeof(data_point_1d_struct), compare);

for (int i = 0; i < data_index; i++) {
  fprintf(outfile, "%.15e """

    for key in out_quantities_dict.keys():
        printf_c_type = "%.15e" if key[0] != "int" else "%d"
        qsort_and_output_to_file += f"{printf_c_type} "

    qsort_and_output_to_file = (
        f'{qsort_and_output_to_file[:-1]}\\n", data_points[i].xCart_axis, '
    )

    for key in out_quantities_dict.keys():
        qsort_and_output_to_file += f"data_points[i].{key[1]}, "

    qsort_and_output_to_file = f"{qsort_and_output_to_file[:-2]});\n}}\n"

    out_string += qsort_and_output_to_file

    return prefunc_content, out_string


def simple_loop_2D(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str = "yz",
) -> str:
    r"""
    Generate a simple 2D loop in the xy or yz plane in C (for use inside of a function).

    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param out_quantities_dict: Dictionary of quantities to be output
    :param plane: Specifies the plane of output: either "xy" or "yz" (default)
    :return: Complete loop code, output as a string.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> diag2d = clang_format(simple_loop_2D("Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4ANwAcNdABfgfJ51OEu5G9+HknfNr147qqTmBC3cI97UzHG3dLgKhgkSaNA4b0gwik0t5tWJAsdm3VWwiKtWc3tY77eSymZwSFkcARy74iuD1V89PH3ZOsNJwCBtX1mytk8Gmin3gskZ4wSKZZgPSc5pFGLs7OVOubWy/d0Lrc/ncXImAq/TSAcJ47NLlzM/3CwIU2eAHMiuxZJwmpmaSClGnE2qvwVqQNi8nvkxi5CqlqeYLqPKA9f8vDt7yZflaQEHrkw1RqU0L0IfFMpz5010pQ/9nMpdfkDimwxpuWsHu2cVv8TIZ/ESnVBWP7mE4GsECRgilpCvFmQhm9W/wvjaWVMWlOTVEtmUHg2EERMWLVHUHcz/v5b/ASCFs1Gs0IlD4VjAyulNbOdNIgGbbbeaOB7/0WxhmAEDncryySD4nFSRQwktzAgAOqBlU6bEt+1aRggHXko0pyEpSsTTK+THaTkiLDdMs5bxBwdTa0o5wY1ACoHcTKrbjjcKjadqzhdt7xs6nOC+yq53HGFcExgh/EDergbb5pjDHC2lvCQIws4s/Igg8sOg9mU9Igtwl3hxXagF/dcUpkcMLb21DPQtDcwGcHppgwAAAG8zxFlBx4kiAAHfA/EGAACHDcpdscRn+wIAAAAABFla")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="yz"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AO5Ac5dABfgfJ51OEu5G9+HknfNr147qqWllBPHI97UzHG9N+QT+URR78q1HU1SAj5n0bvig/XkRzY6UtW1SM0KnjsTwz4bVpsXCoU1zeR3qAriRz/uk4HZriJBGPIOm0e8WoWBi8RwETKvt6gklgr2I2dbftr5I0l340lYiBldQBKs3j3wyZDYfCdcY4AJk/JvA/IIROisi++rX10PvyW0tsJkzkeFPub0bpYMtqoyX79CNMEnDvGMujcXFb+hBHRLHJXIFWM6KQuFUnr7bfHZ7q7eq8TmBkj0+mitMIIjgMOgpXdTNv9uwGAEUTd89/jUMvMT6pn+i3T4WtOklcUpfXCEIjXorwIfFGBYmBXjyPppon9hbkHynYraJeBjWAayRs0cRtu4VRyiXKLqwPwtBV0YoCjy4vd8B09bc+/zt1L6/dOr+iea/oYl+cTHoa0kcKKoy4Yaj7L5sOLpwbfh3h3iUMD8yZ9ca0ZR4gSJjv/7WTjawCS+U9UTbVMtd8OSF+p1umNEt4iMOdcwv8nnb8wZO5eZ2xy9ck2gxC9AmJC8k0OauGJPtmdw9Uiv+WaH85zeAAVgfxeigrPqiupqOGE6uC0vE3GP4mnTLVj/0Ry1xQAAANYZRU9EmcEJAAHqA7oHAACSLqpEscRn+wIAAAAABFla")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4ANwAchdABfgfJ51OEu5G9+HknfNr147qqTmBC3cI97UzHG9N+QT+ocOVsGYU10j3aThYalczu1S8QtgITDuqHTXU1eapZEQw75OCzSgIofGTPSb4MqSNFlMjWOEYG0I710AGvwqxJmZ+I6L4zF7+PP9FYQ/d9u05L6IMNbDRYe9ph3SaiPkChicl8u/Vm+W3iuoTLsc2OFKVtP8t2ctj+5Ism3MCRN1ayJWME1y2FggLDxEjq7DGxfRhfgm8sS8boaDnAZxBrlKwL4hODHwQgdeNF/okej2WLE7ImUBGfXSrGKTSEBpluwwEffZSDJYHjP/La92PP+Np29gNrKzfUxSKtl/1yJ22RbkRDxffiWDt1z1k8gWNAMjYsGzCgrQFBXg8ffUILaRJs1H9PtWJlOfkeL6KCqW8pJ3YKan/UPfadOCiSuEk4Hp+X47FY+GKXdWkg4abjliDPVrE9gzCZC3lpoVz+XGuh3g1ir2KqB/luykcUtCSwgEjj8lIDYgzlwSIdvDCJ/lvXiSq5sTpza+EtIASZ67N+iqq5B98+3g7sjk2PzD6YGAwJ12YuT6LzMcj09oPouFn3b3QYUJLmhGCtjA8v7xV86PNiMsGQD+TryjJfwT/QAB5APxBgAA6fp3M7HEZ/sCAAAAAARZWg==")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
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
        if "Cartesian" in CoordSystem:
            # xy-plane == { z_mid }, where z index is i2
            i2_pts += [Nxx_plus_2NGHOSTS[2] / 2]
        elif (
            "Spherical" in CoordSystem
            or "SymTP" in CoordSystem
            or "Cylindrical" in CoordSystem
        ):
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

    out_string += """// Main loop:
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
  {
"""
    for key, value in out_quantities_dict.items():
        out_string += f"const {key[0]} {key[1]} = {value};\n"
    out_string += 'fprintf(outfile, "%.15e %.15e '

    for key in out_quantities_dict.keys():
        printf_c_type = "%.15e" if key[0] != "int" else "%d"
        out_string += f"{printf_c_type} "

    out_string = f'{out_string[:-1]}\\n", '

    if plane == "xy":
        out_string += "xCart[0], xCart[1], "
    elif plane == "yz":
        out_string += "xCart[1], xCart[2], "

    for key in out_quantities_dict.keys():
        out_string += f"{key[1]}, "

    out_string = f"{out_string[:-2]});\n}}\n}}\n"

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

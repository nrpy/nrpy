"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
"""

from typing import Dict, List, Sequence, Tuple, Union

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
    :param loop_region: Loop over "all points" or "interior" of a numerical grid.
    :param read_xxs: Read the xx[3][:] 1D coordinate arrays if interior dependency exists
    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param enable_rfm_precompute: Enable pre-computation of reference metric
    :param enable_OpenMP: Enable loop parallelization using OpenMP
    :param OMP_custom_pragma: Enable loop parallelization using OpenMP with custom pragma
    :param OMP_collapse: Specifies the number of nested loops to collapse
    :return: The complete loop code as a string.
    :raises ValueError: If `loop_region` is unsupported or if `read_xxs` and `enable_rfm_precompute` are both enabled.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="all points")))
    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
      } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior", OMP_custom_pragma="#CUSTOM_OMP")))
    #CUSTOM_OMP
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=3)))
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    #pragma omp parallel for collapse(3)
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          const MAYBE_UNUSED REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
          const MAYBE_UNUSED REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
          const MAYBE_UNUSED REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
          const MAYBE_UNUSED REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
          const MAYBE_UNUSED REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
          const MAYBE_UNUSED REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
          const MAYBE_UNUSED REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          const MAYBE_UNUSED REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          const MAYBE_UNUSED REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          const MAYBE_UNUSED REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < NGHOSTS+Nxx2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=2)))
    #pragma omp parallel for collapse(2)
    for (int i2 = NGHOSTS; i2 < NGHOSTS + Nxx2; i2++) {
      for (int i1 = NGHOSTS; i1 < NGHOSTS + Nxx1; i1++) {
        const MAYBE_UNUSED REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
        const MAYBE_UNUSED REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
        const MAYBE_UNUSED REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
        const MAYBE_UNUSED REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
        const MAYBE_UNUSED REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
        const MAYBE_UNUSED REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
    <BLANKLINE>
        for (int i0 = NGHOSTS; i0 < NGHOSTS + Nxx0; i0++) {
          const MAYBE_UNUSED REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          const MAYBE_UNUSED REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          const MAYBE_UNUSED REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          const MAYBE_UNUSED REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          const MAYBE_UNUSED REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < NGHOSTS+Nxx0; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < NGHOSTS+Nxx1; i1++)
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
    # offset upper bounds by 1
    elif loop_region == "interior plus one upper":
        i2i1i0_mins = ["NGHOSTS", "NGHOSTS", "NGHOSTS"]
        i2i1i0_maxs = ["NGHOSTS+Nxx2+1", "NGHOSTS+Nxx1+1", "NGHOSTS+Nxx0+1"]
    else:
        raise ValueError(
            f'loop_region = {loop_region} unsupported. Choose "", "all points", or "interior"'
        )

    read_rfm_xx_arrays = ["", "", ""]
    # 'Read_xxs': read the xx[3][:] 1D coordinate arrays, as some interior dependency exists
    if read_xxs:
        if not enable_simd:
            read_rfm_xx_arrays = [
                "const MAYBE_UNUSED REAL xx0 = xx[0][i0];",
                "const MAYBE_UNUSED REAL xx1 = xx[1][i1];",
                "const MAYBE_UNUSED REAL xx2 = xx[2][i2];",
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
    Generate a C code snippet to output data along a specified axis in a given coordinate system.
    The generated code includes a loop that considers the points closest to the specified axis
    in the provided coordinate system.

    :param CoordSystem: Specifies the coordinate system (e.g., "Cartesian" or "Spherical").
    :param out_quantities_dict: Dictionary containing quantities tooutput. Keys are tuples of (type, name)
                                and values are the C code expressions to compute the quantities.
    :param axis: The axis along which the output is generated; accepts either "y" or "z" (default is "z").
    :return: A tuple containing two strings: the first string represents the struct and comparison function definitions,
             and the second string represents the loop code to output data along the specified axis.
    :raises ValueError: If the provided axis is not "y" or "z",
                        or if the CoordSystem is not supported by this function.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AT1AmxdABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsG+ldk9FmdPbsp0g6RWC5bdGf2YzJ2K/2yuxKCa6CL2nIyRvDVjakbAbX5BMu291HujDeC3jpHflLFFuxS9CUE77JWgfezCN4uqB0IBXg01sUOwGQ2i+BDvk6FyLLzsDHFpcsvComWnDXUETQzVEOZKdTNYgFKLGGNyU9TQjwoS1YkGaTTGweGbY5H/3rJYLSgRQdHenZdtNFoOp+9J1Uj180hax8hJuSI7z+fVa5SfdfdFGFy45TsivVnzwazf0qr7jwpHvr+8zLE4IvP7awONSsmSBZCByG03t+3GAMhtl/LdR8f/VgBeK6njZ7VKg7psNtaNXDSERWJ7JgfXsIQqY0+nK+ed87Mo0iet7Z0JdEQgnUklobffK4Hi7BkVcG9SmDOpW6lGUCvMHl4AGmhuYyTkhv50rVTlGWvqELC2KmmoQXbX30HlPx8tYu15p1fjHJOA6lh8GB24lRndmsh/nNhreAggDe9WD18G0Q80enLlxpmXGPWu75tc765q75DO+DztI6FnMvHAeaikx4Jzslw7+H/Em7JHXk5mRdYCLbu85k8pLorsFlQ7XXzoJz4MvAeGZzV3SWi3a4gpl3dqVOJ/vCzDthCi4igjAoaItvxbg6f1EmqmcVOTBln2sk4ZWI2s7IQfUePg/sCERQt6D1Krs5gLZZqP31hCsLQGmgWHhGbl5YznRW6fWO1LEasTyrEnIeOj16OHPzi6v+yubnui20rSduFqb1tDoFjp2PMQUxSpNP0eZbTSVGGeTtaFEE5+00ABRou1R79IkKAAGIBfYJAADsv4pVscRn+wIAAAAABFla")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AU+AoNdABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsG+ldk9FmdPbsrjX64N9akGaoGUFEU5ORWsN5geYreRK86LYVVk5VrW5b1up8AoOcyENkotWRT5RPjzooP+Q6k4ipyEtZJrJd/N3xnTKblsyW9tSweJzngYjTMuNewwn9QJy+Alef+94DLS7fJxXb5xiCY9iYCgUlrO0XH5QCq60B9/9jsrkmmj0vldh4TbEalc9GJgwMvkshPkZqxB9DNrSg4AUJAb8wwag8Q8vaQ2K6Z6dvdO7SjrNpi40eG+S/fvFXS6ufn0z3TWEan+g1cI7PZfnPayTShUzVqO/RUQqU/+1+scsPlyX5pxUBrwMSzzv1SX1ylejzal7vPc8W0KOl1DpzBYMuDnlSuNMpKUk94v062nOik6IKAtnloKUOxXQckuXqqJq0DwsxEto10Fg1J8flB4iGQkxSzr6JsBaSjVqZAMs4rECWTIiro13/p2suO6IeVRgHwPSZM8l0YL85mRwYJwfwmGBfrnVcN9KysQ/h4n7Ps/JgewrHIeeJcEevEu+Bq9cj0iYxUgMGTowG+lybAc7+mtEAjiFkFGMZZ4zvfv75gkBj4D0bva2hff5VklOFJ0hcEBfDWmcX+1NlJsyzdL2v7Ti0NT2T8cs+CXBi5KW9wXMnTy54StxXVes5nQ3TvM1/V0fGWxPA237DakDRbL7axkxVf/jZmj+dOBEOJLJaGHidkbVq+V5tdLqbor7vN81FGhR+wN6afUJdp8Nh343uEQ+aKGgnY3TJH8A8qMoPvnLsdE19jpO35V5RaMmv9+RxwXBet4SzUqBIhOiPiKLB+oOhNN0AAPIhA50NSPdIAAGfBb8KAAAhJAOvscRn+wIAAAAABFla")
    >>> if diag1d != expected_string:
    ...     raise ValueError(f"{diff_strings(expected_string, diag1d)}\n base64-encoded output: {compress_string_to_base64(diag1d)}")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Spherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="z")[1])
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AT5AnJdABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsHGY1k9FmdPbsrhibF9Nv4BiDsCPObVIY7ZqQYVYp+DDPK6WTG4MxWRmkapDc/p1H6FPzwc2DnX8W5QSTRqYhnZDVNw5Mel0AuikwKDrF+7C+gE58Hvs1mL3Ikd8RS45EiAyBBheadz676ml8hYqhnxQ3hX7HWNx4aPUbn1MGOamPEvPuxux3UtB6G/YYMMoIkZAcbfTLThEP50XMhkvcCm7zW1fQ1tPHCBvwkckib+cMDWmyCB3b4erz76XsbGOwBNKWOzGHeeYfFeYZ5ob/qpsgF+QhP7Ab/47Noe+WXfzJ6bFJ+gk4ARqfvykWkqeBXrEkQpMAeK3+ra2t/dqUcFvibnZcsc+YOvwdadDwMyZyUXZA0i4qjm8AlR0x7P1bb3z4DsEcJqfOXsXpaV6eHtHloqzY19rfs5e+HhhBeHguaKGdqzetDlIRlpzVJyf2T8MFeEw+btvoSO32O3kd2XWiTMGgok6FDbJq7LlKJDEAZIEfYk0nDV+45cXBPI9IS9BEfHXPHKtHxTTY4gSpNo8PdC+GyqTGj0VW92qU/mjigkmjahkNAmjxV1aCDJKRYSpwIpHvUWazqCNqAuBf5UiFWb6YnP5E2KJ+sK1n0z/YCuUDSpprzoVnLeFbYbkvWih73GGjWCOAnNw2zfE4CsGsrEQHf+M0D7YjBmBqbgxHjG6vNSh2tfKiKw8Is4SwgIJNRILR2oumL51Hq6lHaJMY6f8QTObyDe9O2vES90TtANOMLZxnDQZKBrEebiGvXNhjo+ErTUSicnUAAAAAFyWFPQZVMbkAAY4F+gkAAEnjBcmxxGf7AgAAAAAEWVo=")
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
        elif "Wedge" in CoordSystem:
            # NO POINTS ON Y AXIS
            i0_pts += [-1]
            i1_pts += [-1]
            i2_pts += [-1]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    if axis == "z":
        if "Cartesian" in CoordSystem:
            # z-axis == { x_mid, y_mid }
            i0_pts += [Nxx_plus_2NGHOSTS[0] / 2]
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
        elif "Spherical" in CoordSystem:
            if "Ring" in CoordSystem:
                # NO POINTS ON Z AXIS
                i0_pts += [-1]
                i1_pts += [-1]
                i2_pts += [-1]
            else:
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
            # FIXME: Missing points between foci (not an easy fix).
            i1_pts += [NGHOSTS]
            i1_pts += [Nxx_plus_2NGHOSTS[1] - NGHOSTS - 1]
            i2_pts += [NGHOSTS]
        elif "Wedge" in CoordSystem:
            # Wedge-like: same as Spherical except x_new = +/-z_old, y_new = y_old, z_new = x_old
            # Thus the z-axis here is the same as the +x-axis in Spherical-like.
            # +x-axis == { theta_mid, phi={phi_mid} (since phi goes from -pi to pi) }
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
            i2_pts += [Nxx_plus_2NGHOSTS[2] / 2]
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
    out_string = f"""// Define points for output along the {axis}-axis in {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];

data_point_1d_struct data_points[numpts_i0 * numpts_i1 * numpts_i2];
int data_index = 0;
"""

    # Loop body for storing results.
    loop_body_store_results = ""
    # Continue to append to loop_body_store_results here...
    loop_body_store_results += f"""{{
// Store the data in the data_point_1d_struct
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
    out_string += f"""// Main loop to store data points along the specified axis
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);

  {loop_body_store_results}
}}
"""

    # Post-loop: qsort() along xCart_axis and output to file.
    prefunc_content = """// Struct to hold 1D data points
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


def max_numpts__i012_pts__numpts_2D(
    CoordSystem: str, plane: str, include_GHOSTS: bool = False
) -> Tuple[Sequence[sp.Expr], List[List[sp.Expr]], List[Union[int, sp.Expr]]]:
    """
    Determine the points and number of points for a 2D loop output in the specified coordinate system and plane.

    :param CoordSystem: Specifies the coordinate system (e.g., "Cartesian", "Spherical").
    :param plane: The plane to consider for the output; either "xy" or "yz".
    :param include_GHOSTS: Whether to include ghost zones in the computation of max_numpts. Defaults to False.
    :return: A tuple containing:
        - max_numpts: A list of the maximum number of points in each direction, either Nxx[3] or Nxx_plus_2NGHOSTS[3] depending on include_GHOSTS.
        - i012_pts: A list of points for each coordinate, structured as [[i0_pts], [i1_pts], [i2_pts]].
        - numpts: The number of interior points along each coordinate, structured as [numpts_x, numpts_y, numpts_z].
    :raises ValueError: If an unsupported `plane` or `CoordSystem` is provided.

    Doctests:
    >>> for Coord in ["SinhCartesian", "HoleySinhSpherical", "Wedge", "SinhCylindrical"]:
    ...     print(Coord, max_numpts__i012_pts__numpts_2D(Coord, "yz"))
    SinhCartesian ([Nxx0, Nxx1, Nxx2], [[Nxx_plus_2NGHOSTS0/2], [], []], [1, Nxx1, Nxx2])
    HoleySinhSpherical ([Nxx0, Nxx1, Nxx2], [[], [], [NGHOSTS + Nxx2/4 - 1/2, NGHOSTS + 3*Nxx2/4 - 1/2]], [Nxx0, Nxx1, 2])
    Wedge ([Nxx0, Nxx1, Nxx2], [[], [Nxx_plus_2NGHOSTS1/2], []], [Nxx0, 1, Nxx2])
    SinhCylindrical ([Nxx0, Nxx1, Nxx2], [[], [NGHOSTS + Nxx1/4 - 1/2, NGHOSTS + 3*Nxx1/4 - 1/2], []], [Nxx0, 2, Nxx2])
    >>> for Coord in ["SinhCartesian", "HoleySinhSpherical", "Wedge", "SinhCylindrical"]:
    ...     print(Coord, max_numpts__i012_pts__numpts_2D(Coord, "xy"))
    SinhCartesian ([Nxx0, Nxx1, Nxx2], [[], [], [Nxx_plus_2NGHOSTS2/2]], [Nxx0, Nxx1, 1])
    HoleySinhSpherical ([Nxx0, Nxx1, Nxx2], [[], [Nxx_plus_2NGHOSTS1/2], []], [Nxx0, 1, Nxx2])
    Wedge ([Nxx0, Nxx1, Nxx2], [[], [], [NGHOSTS + Nxx2/4 - 1/2, NGHOSTS + 3*Nxx2/4 - 1/2]], [Nxx0, Nxx1, 2])
    SinhCylindrical ([Nxx0, Nxx1, Nxx2], [[], [], [Nxx_plus_2NGHOSTS2/2]], [Nxx0, Nxx1, 1])
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
        elif "Wedge" in CoordSystem:
            # UWedgeHSinhSph: same as Spherical except x_new = -z_old, y_new = y_old, z_new = x_old
            # Thus the xy plane here is the same as the -z,y plane in Spherical-like.
            # LWedgeHSinhSph: same as Spherical except x_new = z_old, y_new = y_old, z_new = -x_old
            # Thus the xy plane here is the same as the z,y plane in Spherical-like
            # (for both see discussion below for yz plane in Spherical)
            i2_pts += [NGHOSTS + sp.Rational(1, 4) * Nxx[2] - sp.Rational(1, 2)]
            i2_pts += [NGHOSTS + sp.Rational(3, 4) * Nxx[2] - sp.Rational(1, 2)]
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
        elif "Wedge" in CoordSystem:
            # UWedgeHSinhSph: same as Spherical except x_new = -z_old, y_new = y_old, z_new = x_old
            # Thus the yz plane here is the same as the y,x plane in Spherical-like.
            # LWedgeHSinhSph: same as Spherical except x_new = z_old, y_new = y_old, z_new = -x_old
            # Thus the yz plane here is the same as the y,-x plane in Spherical-like
            # xy-plane == { theta_mid }, where theta index is i1
            i1_pts += [Nxx_plus_2NGHOSTS[1] / 2]
        else:
            raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    i012_pts = [i0_pts, i1_pts, i2_pts]
    max_numpts = Nxx[:]
    if include_GHOSTS:
        max_numpts = Nxx_plus_2NGHOSTS[:]

    numpts: List[Union[int, sp.Expr]] = [0] * 3
    numpts[0] = len(i0_pts) if i0_pts else max_numpts[0]
    numpts[1] = len(i1_pts) if i1_pts else max_numpts[1]
    numpts[2] = len(i2_pts) if i2_pts else max_numpts[2]
    return max_numpts, i012_pts, numpts


def simple_loop_2D(
    CoordSystem: str,
    out_quantities_dict: Dict[Tuple[str, str], str],
    plane: str = "yz",
) -> str:
    r"""
    Generate a C code snippet to output data in a specified 2D plane within a given coordinate system.
    The generated code includes a loop that considers the points closest to the specified plane
    in the provided coordinate system.

    :param CoordSystem: Specifies the coordinate system (e.g., "Cartesian" or "Spherical").
    :param out_quantities_dict: Dictionary containing output quantities. Keys are tuples of (type, name)
                                and values are the C code expressions to compute the quantities.
    :param plane: The plane along which the output is generated; accepts either "xy" or "yz" (default is "yz").
    :return: A string containing the complete loop code to output data in the specified plane.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, compress_string_to_base64, decompress_base64_to_string, diff_strings
    >>> diag2d = clang_format(simple_loop_2D("Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4APOAftdABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsG4xvVIEaWbxkX3h8/EJehlGF5zbTysnCr1xUlJMEf1xS4Esouq+kqM3dFCeiSE9KLTnRNjP8xwHRV5EjkdF9Ji2Hy2enb9GoPlneKNCL2RRr/TqFRCUVfoKiOVuN+geLMLo0B85GbnSI9zBYCvpsLs3c5XO+zsrh/XucMLc4uRaNfj7P8ZDnLJsK+wRcx4uJFx+rkna25NtKRiL3QtCBq1xsggBsRTNwfVbTh9eW4lV9bLeo/FJBk8/nToG2Aodz5di7d7CTjJbWu6eGxKZeXBYxkdnH45MU1JjC57hsAUK+N8iBcIkCac0OUtwGapBAhoSwR3p/T4hdD4VANJw+YWx2WOqQcnjyTIgR0T5BWxC4lKglJgUsKX6RhIcpePnhQXn/rt1oyB6aZGu3TPUNQ9mRmX80+CswFwQ59J5Hd+gj+iOSvfLFuqFqgUxSIaH7yO3YFqiqQYW6hOiWfhvItttDK0CQ8aml50KKYIp5dVrDlKkcnTGFQCEnEauBEZgjkYqhhSYyD9ZBKqCPAlkyVeOw8P3lNblD4gTIj8Ix1jEBnT6CtUG/8NB0T+gyZm/hYfccJ0fiLTWV1H3ffpFDVwA2a7fwt0C3yQT4cSdSF60koEVPTgAAAOF94MkNNdgAAZcEzwcAADO81x2xxGf7AgAAAAAEWVo=")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="yz"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4AQXAgldABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsHAm9tIEaWbxkX3iJfpkDNY8C3AcF4owUwxhZ1N50NCItxBlKVFj80aYMhsZ4PkyoOwyyvmAAMTbf01+rgVRy41HN2nLVCbkW0hxsedloR6IXLH+fGmJi764mFKuVGZLhfX+Tel+UBrQdO1IBsp/IcXVHn5TEFnapERVo0OfRF7JB8vARs3EpWRjs8Fcf4i/aEqkPmLxknRjub9f8xbgKCffuhVAjTpRP8V9pbnhxelBKeLmESENqTKb60kB5htC4++puzsjfuQ7reEf1hJKp0nZu0kbkuv2kQshySkF/o6+7Mi2bPxxb3KG81om068KUADAB8CieYaF32WO6DWMX6GI1tqfEELfMu9LyHaGeXy4HduPCE1GxeGZVZ0H/ITMHWiWtNujJsn+ZlVm4jxmSfGwBu07LVh/rRM/d20jjo6x5cNbJ8eujc0zJUmEEi1RVsJmez2frRDeBXbxe0hh4BhD1fRvr0AFlTdtJOU0JHJ2JGfSm1Yr9fj8XKD+grTBW971jCsydOjI8srbCqFcasnhwboxpZlcKCQb6CFmxZeKOFEo0DMP9RMAFl/8xOLas46NveAzGKmJlbtcNAi49kVOO6SV8PmB5DjYguHq3Czgwdx4kmUEzfugEySAZW4jk9YgAAAAAAJUIHs5wIG6JAAGlBJgIAACzauMJscRn+wIAAAAABFla")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> expected_string = decompress_base64_to_string("/Td6WFoAAATm1rRGAgAhARwAAAAQz1jM4APOAf9dABfgfIhlMYIN4cUPfyespQRzjUH+1RgEkebsTINj5AD0Ew8BQihsG4xvVIEaWbxkX3iJfpkDOlQAENwZ5pYH87gVaph5MAUkAZEfp+JtAaWbFdmwHpyodWegADbAuSnNidKn5pe3IwsI+doFc8n4j9JL1477t0UHac/0hpu7A5dnrcNg2ueX4Caw+7+MInN/Nab5zUoCpxLj3OK8Hep9zbPnrVQGtMmXlClwnI3oakoIqQNPo5SvIaamYnJfy5Anr2vJtlGDBTYgk70dNqPYMNE2QXL+iReM80EnlTcze+fkO89qlnpN8cJODamsANsDSSq+UGeokPPRtEcN7UHwLAuE394uh7hyfdcwE0WtOH50jDjodOlWrxPAIQLiF3airOTIheQbqJkaRXZcRT3amRSElh5eV5W7QpAh1Wf2teEemYDycl7k9hZvViYf+9INCI1Mqc3TGcMLBoT8F+ANs3wFPRGxEwUf/yBSE8K85+uQYTerPLBy4xJ0upJtm8gdVcjTX7gIWRnxohyfu2/bQE2Yk74pZx8BlK2BdeVmnVp1h6Yt2GQLNi8pDHk3Leyf/QMowCKbMDvnemvNAk9fwTKEygPUKgXgY2/WAszwW7PB0wKHZjUtY77YgCg89lfpEBOMG/kRqdVkTa/B+ajuaVw9PCXEyxUbuq5bONYqxwAAAAB+W+yG0xvoAAGbBM8HAABIfBVqscRn+wIAAAAABFla")
    >>> if diag2d != expected_string:
    ...     raise ValueError(f"Diff:\n{diff_strings(expected_string, diag2d)}\n base64-encoded output: {compress_string_to_base64(diag2d)}")
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

    out_string += """// Main loop to store data points in the specified plane
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  xx_to_Cart(commondata, params, xx, i0, i1, i2, xCart);
"""
    out_string += "{\n"
    out_string += "// Collect diagnostic data\n"
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

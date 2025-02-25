"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
Contributor: Samuel D. Tootle
Email: sdtootle **at** gmail **dot** com
"""

from typing import Dict, List, Sequence, Tuple, Union

import sympy as sp

import nrpy.helpers.loop as lp
import nrpy.indexedexp as ixp

implemented_loop_regions = ["", "all points", "interior"]


def implemented_loop_regions_err(loop_region: str) -> str:
    """
    Generate the string that is printed when a loop_region is not defined.

    :param loop_region: string denoting the intended loop region
    :returns: the error message describing the loop_region was not found
    """
    regions = [f'"{region}"' for region in implemented_loop_regions]
    regions[-1] = f" or {regions[-1]}"
    tmp_str = ",".join(regions)
    return f"loop_region = {loop_region} unsupported. Choose {tmp_str}."


def get_loop_region_ranges(
    loop_region: str, min_idx_prefix: Union[str, None] = None
) -> Tuple[List[str], List[str]]:
    """
    Return Loop region index ranges.

    :param loop_region: Loop region
    :param min_idx_prefix: String that specifies the starting value prefix (GPU specific)
    :returns: region indicies
    """
    i2i1i0_mins = ["", "", ""]
    i2i1i0_maxs = ["", "", ""]

    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "all points":
        if min_idx_prefix is not None:
            i2i1i0_mins = [f"{min_idx_prefix}{i}" for i in reversed(range(3))]
        else:
            i2i1i0_mins = ["0", "0", "0"]
        i2i1i0_maxs = ["Nxx_plus_2NGHOSTS2", "Nxx_plus_2NGHOSTS1", "Nxx_plus_2NGHOSTS0"]
    # 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
    elif loop_region == "interior":
        if not min_idx_prefix is None:
            i2i1i0_mins = [f"{min_idx_prefix}{i}+NGHOSTS" for i in reversed(range(3))]
        else:
            i2i1i0_mins = ["NGHOSTS", "NGHOSTS", "NGHOSTS"]
        i2i1i0_maxs = [
            "Nxx_plus_2NGHOSTS2 - NGHOSTS",
            "Nxx_plus_2NGHOSTS1 - NGHOSTS",
            "Nxx_plus_2NGHOSTS0 - NGHOSTS",
        ]

    return i2i1i0_mins, i2i1i0_maxs


def simple_loop(
    loop_body: str,
    enable_intrinsics: bool = False,
    loop_region: str = "",
    read_xxs: bool = False,
    CoordSystem: str = "Cartesian",
    enable_rfm_precompute: bool = False,
    enable_OpenMP: bool = True,
    OMP_custom_pragma: str = "",
    OMP_collapse: int = 1,
    parallelization: str = "openmp",
) -> str:
    """
    Generate a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param enable_intrinsics: Enable SIMD support
    :param loop_region: Loop over "all points" or "interior" of a numerical grid.
    :param read_xxs: Read the xx[3][:] 1D coordinate arrays if interior dependency exists
    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param enable_rfm_precompute: Enable pre-computation of reference metric
    :param enable_OpenMP: Enable loop parallelization using OpenMP
    :param OMP_custom_pragma: Enable loop parallelization using OpenMP with custom pragma
    :param OMP_collapse: Specifies the number of nested loops to collapse
    :param parallelization: Parallelization method to use. Default is "openmp".
    :return: The complete loop code as a string.
    :raises ValueError: If `loop_region` is unsupported or if `read_xxs` and `enable_rfm_precompute` are both enabled.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="all points")))
    #pragma omp parallel for
    for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++) {
      for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++) {
        for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    <BLANKLINE>
          // <INTERIOR>
        } // END LOOP: for (int i0 = 0; i0 < Nxx_plus_2NGHOSTS0; i0++)
      } // END LOOP: for (int i1 = 0; i1 < Nxx_plus_2NGHOSTS1; i1++)
    } // END LOOP: for (int i2 = 0; i2 < Nxx_plus_2NGHOSTS2; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior", OMP_custom_pragma="#CUSTOM_OMP")))
    #CUSTOM_OMP
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
    <BLANKLINE>
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=3)))
    Setting up reference_metric[SinhSymTP_rfm_precompute]...
    #pragma omp parallel for collapse(3)
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
          MAYBE_UNUSED const REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
          MAYBE_UNUSED const REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
          MAYBE_UNUSED const REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
          MAYBE_UNUSED const REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
          MAYBE_UNUSED const REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
          MAYBE_UNUSED const REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
          MAYBE_UNUSED const REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          MAYBE_UNUSED const REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          MAYBE_UNUSED const REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          MAYBE_UNUSED const REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
    <BLANKLINE>
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
    <BLANKLINE>
    >>> print(clang_format(simple_loop('// <INTERIOR>', loop_region="interior",
    ...       CoordSystem="SinhSymTP", enable_rfm_precompute=True, OMP_collapse=2)))
    #pragma omp parallel for collapse(2)
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        MAYBE_UNUSED const REAL f1_of_xx1 = rfmstruct->f1_of_xx1[i1];
        MAYBE_UNUSED const REAL f1_of_xx1__D1 = rfmstruct->f1_of_xx1__D1[i1];
        MAYBE_UNUSED const REAL f1_of_xx1__DD11 = rfmstruct->f1_of_xx1__DD11[i1];
        MAYBE_UNUSED const REAL f4_of_xx1 = rfmstruct->f4_of_xx1[i1];
        MAYBE_UNUSED const REAL f4_of_xx1__D1 = rfmstruct->f4_of_xx1__D1[i1];
        MAYBE_UNUSED const REAL f4_of_xx1__DD11 = rfmstruct->f4_of_xx1__DD11[i1];
    <BLANKLINE>
        for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++) {
          MAYBE_UNUSED const REAL f0_of_xx0 = rfmstruct->f0_of_xx0[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__D0 = rfmstruct->f0_of_xx0__D0[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__DD00 = rfmstruct->f0_of_xx0__DD00[i0];
          MAYBE_UNUSED const REAL f0_of_xx0__DDD000 = rfmstruct->f0_of_xx0__DDD000[i0];
          MAYBE_UNUSED const REAL f2_of_xx0 = rfmstruct->f2_of_xx0[i0];
          MAYBE_UNUSED const REAL f2_of_xx0__D0 = rfmstruct->f2_of_xx0__D0[i0];
          MAYBE_UNUSED const REAL f2_of_xx0__DD00 = rfmstruct->f2_of_xx0__DD00[i0];
    <BLANKLINE>
          // <INTERIOR>
        } // END LOOP: for (int i0 = NGHOSTS; i0 < Nxx_plus_2NGHOSTS0 - NGHOSTS; i0++)
      } // END LOOP: for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++)
    } // END LOOP: for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++)
    <BLANKLINE>
    """
    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "":
        return loop_body
    elif loop_region not in implemented_loop_regions:
        raise ValueError(implemented_loop_regions_err(loop_region))
    i2i1i0_mins, i2i1i0_maxs = get_loop_region_ranges(
        loop_region, min_idx_prefix="tid" if parallelization == "cuda" else None
    )

    read_rfm_xx_arrays = ["", "", ""]

    # SIMD is used only if intrinsics are enabled and we're not, currently, using CUDA.
    use_simd = enable_intrinsics and (parallelization != "cuda")

    # Determine if a reset is needed since it is only relevant to openmp
    OMP_collapse = OMP_collapse if parallelization == "openmp" else 1

    # 'Read_xxs': read the xx[3][:] 1D coordinate arrays, as some interior dependency exists
    if not use_simd and read_xxs:
        if parallelization == "cuda":
            read_rfm_xx_arrays = [
                "MAYBE_UNUSED const REAL xx0 = x0[i0];",
                "MAYBE_UNUSED const REAL xx1 = x1[i1];",
                "MAYBE_UNUSED const REAL xx2 = x2[i2];",
            ]
        else:
            read_rfm_xx_arrays = [
                "MAYBE_UNUSED const REAL xx0 = xx[0][i0];",
                "MAYBE_UNUSED const REAL xx1 = xx[1][i1];",
                "MAYBE_UNUSED const REAL xx2 = xx[2][i2];",
            ]
    elif read_xxs and use_simd:
        raise ValueError("no innerSIMD support for Read_xxs (currently).")

    # 'enable_rfm_precompute': enable pre-computation of reference metric
    if enable_rfm_precompute:
        if read_xxs:
            raise ValueError(
                "enable_rfm_precompute and Read_xxs cannot both be enabled."
            )
        # pylint: disable=C0415
        from nrpy.infrastructures.BHaH import rfm_precompute

        rfmp = rfm_precompute.ReferenceMetricPrecompute(
            CoordSystem, parallelization=parallelization
        )
        if enable_intrinsics:
            read_rfm_xx_arrays = [
                rfmp.readvr_intrinsics_inner_str[0],
                rfmp.readvr_intrinsics_outer_str[1],
                rfmp.readvr_intrinsics_outer_str[2],
            ]
        else:
            read_rfm_xx_arrays = [
                rfmp.readvr_str[0],
                rfmp.readvr_str[1],
                rfmp.readvr_str[2],
            ]
    # 'DisableOpenMP': disable loop parallelization using OpenMP
    if enable_OpenMP or OMP_custom_pragma != "":
        if OMP_custom_pragma == "" and parallelization == "openmp":
            pragma = "#pragma omp parallel for"
            if OMP_collapse > 1:
                pragma = f"#pragma omp parallel for collapse({OMP_collapse})"
        # 'OMP_custom_pragma': enable loop parallelization using OpenMP with custom pragma
        else:
            pragma = OMP_custom_pragma
    else:
        pragma = ""
    increment = (
        ["stride2", "stride1", "stride0"]
        if parallelization == "cuda"
        else (["1", "1", "simd_width"] if enable_intrinsics else ["1", "1", "1"])
    )

    loop_body = read_rfm_xx_arrays[0] + f"\n\n{loop_body}"
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

    full_loop_body = str(
        lp.loop(
            ["i2", "i1", "i0"],
            i2i1i0_mins,
            i2i1i0_maxs,
            increment,
            prefix_loop_with,
            loop_body=loop_body,
        )
    )

    return full_loop_body


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
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> validate_strings(diag1d, "Cartesian_y_axis")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> validate_strings(diag1d, "SinhSpherical_y_axis")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Spherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="z")[1])
    >>> validate_strings(diag1d, "Spherical_z_axis")
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
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> diag2d = clang_format(simple_loop_2D("Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> validate_strings(diag2d, "Cartesian_xy_plane")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="yz"))
    >>> validate_strings(diag2d, "SinhSpherical_yz_plane")
    >>> diag2d = clang_format(simple_loop_2D(CoordSystem="SinhSymTP", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, plane="xy"))
    >>> validate_strings(diag2d, "SinhSymTP_xy_plane")
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

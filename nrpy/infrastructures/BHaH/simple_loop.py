"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
Contributor: Samuel D. Tootle
Email: sdtootle **at** gmail **dot* com
"""

from typing import List, Tuple, Union

import nrpy.helpers.loop as lp
import nrpy.params as par

implemented_loop_regions = ["", "all points", "interior", "interior plus one upper"]


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
    Return Loop region index ranges using a data-driven approach for clarity.

    :param loop_region: Loop region.
    :param min_idx_prefix: String that specifies the starting value prefix (GPU specific).
    :returns: A tuple of two lists: loop minimums and loop maximums.
    """
    # Define base components for ranges
    if min_idx_prefix:
        all_points_min = [f"{min_idx_prefix}{i}" for i in reversed(range(3))]
        interior_min = [f"{min_idx_prefix}{i}+NGHOSTS" for i in reversed(range(3))]
    else:
        all_points_min = ["0", "0", "0"]
        interior_min = ["NGHOSTS", "NGHOSTS", "NGHOSTS"]

    all_points_max = [f"Nxx_plus_2NGHOSTS{i}" for i in reversed(range(3))]
    interior_max = [f"{m} - NGHOSTS" for m in all_points_max]
    interior_plus_one_max = [f"{m} + 1" for m in interior_max]

    # Configuration dictionary maps loop_region to its min/max ranges
    region_map = {
        "all points": (all_points_min, all_points_max),
        "interior": (interior_min, interior_max),
        "interior plus one upper": (interior_min, interior_plus_one_max),
    }

    # Default to empty ranges if loop_region is not found (e.g., "").
    return region_map.get(loop_region, (["", "", ""], ["", "", ""]))


def _simple_loop_pragma(
    enable_OpenMP: bool,
    OMP_custom_pragma: str,
    OMP_collapse: int,
    parallelization: str,
) -> str:
    """
    Generate the OpenMP pragma string.

    :param enable_OpenMP: Flag to enable OpenMP parallelization.
    :param OMP_custom_pragma: A custom OpenMP pragma string to use instead of the default.
    :param OMP_collapse: The number of loops to collapse in the pragma.
    :param parallelization: The parallelization strategy (e.g., "openmp", "cuda").
    :return: The generated OpenMP pragma string, or an empty string if not applicable.
    """
    if OMP_custom_pragma:
        return OMP_custom_pragma
    if enable_OpenMP and parallelization == "openmp":
        pragma = "#pragma omp parallel for"
        if OMP_collapse > 1:
            pragma += f" collapse({OMP_collapse})"
        return pragma
    return ""


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
    if not loop_region:
        return loop_body
    if loop_region not in implemented_loop_regions:
        raise ValueError(implemented_loop_regions_err(loop_region))

    if enable_rfm_precompute and read_xxs:
        raise ValueError("enable_rfm_precompute and Read_xxs cannot both be enabled.")

    parallelization = par.parval_from_str("parallelization")
    is_cuda = parallelization == "cuda"
    use_simd = enable_intrinsics and not is_cuda
    if read_xxs and use_simd:
        raise ValueError("no innerSIMD support for Read_xxs (currently).")

    increment = (
        ["stride2", "stride1", "stride0"]
        if is_cuda
        else ["1", "1", "simd_width" if use_simd else "1"]
    )

    min_idx_prefix = "tid" if is_cuda else None
    i2i1i0_mins, i2i1i0_maxs = get_loop_region_ranges(loop_region, min_idx_prefix)

    rfm_reads = ["", "", ""]
    if enable_rfm_precompute:
        # pylint: disable=C0415
        from nrpy.infrastructures import BHaH

        rfmp = BHaH.rfm_precompute.ReferenceMetricPrecompute(CoordSystem)
        rfm_reads = (
            [
                rfmp.readvr_intrinsics_inner_str[0],
                rfmp.readvr_intrinsics_outer_str[1],
                rfmp.readvr_intrinsics_outer_str[2],
            ]
            if enable_intrinsics
            else list(rfmp.readvr_str)
        )
    elif read_xxs:
        rfm_reads = (
            [
                "MAYBE_UNUSED const REAL xx0 = x0[i0];",
                "MAYBE_UNUSED const REAL xx1 = x1[i1];",
                "MAYBE_UNUSED const REAL xx2 = x2[i2];",
            ]
            if is_cuda
            else [
                "MAYBE_UNUSED const REAL xx0 = xx[0][i0];",
                "MAYBE_UNUSED const REAL xx1 = xx[1][i1];",
                "MAYBE_UNUSED const REAL xx2 = xx[2][i2];",
            ]
        )
    rfm_i0_code, rfm_i1_code, rfm_i2_code = rfm_reads

    effective_collapse = OMP_collapse if parallelization == "openmp" else 1
    pragma = _simple_loop_pragma(
        enable_OpenMP, OMP_custom_pragma, effective_collapse, parallelization
    )

    prefix_i2, prefix_i1, prefix_i0 = pragma, "", ""
    body_prefix_parts = []

    if effective_collapse == 3:
        body_prefix_parts.extend([rfm_i2_code, rfm_i1_code, rfm_i0_code])
    elif effective_collapse == 2:
        prefix_i0 = "".join(filter(None, [rfm_i2_code, rfm_i1_code]))
        body_prefix_parts.append(rfm_i0_code)
    else:  # Default case: collapse = 1
        prefix_i1 = rfm_i2_code
        prefix_i0 = rfm_i1_code
        body_prefix_parts.append(rfm_i0_code)

    prefix_loop_with = [prefix_i2, prefix_i1, prefix_i0]
    body_prefix = "".join(filter(None, body_prefix_parts))

    final_loop_body = f"{body_prefix}\n\n{loop_body}"

    return str(
        lp.loop(
            ["i2", "i1", "i0"],
            i2i1i0_mins,
            i2i1i0_maxs,
            increment,
            prefix_loop_with,
            loop_body=final_loop_body,
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

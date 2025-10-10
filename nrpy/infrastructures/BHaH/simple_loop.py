"""
Simple loop generation for use within the BHaH infrastructure.

Author: Zachariah B. Etienne
Email: zachetie **at** gmail **dot* com
Contributor: Ken Sible
Email: ksible *at* outlook *dot* com
Contributor: Samuel D. Tootle
Email: sdtootle **at** gmail **dot* com
"""

from typing import Dict, List, Sequence, Tuple, Union

import sympy as sp

import nrpy.helpers.loop as lp
import nrpy.indexedexp as ixp
import nrpy.params as par

implemented_loop_regions = ["", "all points", "interior", "interior plus one upper"]

# NRPy-level symbolic expressions for grid properties
_NGHOSTS = sp.Symbol("NGHOSTS", real=True)
_Nxx = ixp.declarerank1("Nxx")
_Nxx_plus_2NGHOSTS = ixp.declarerank1("Nxx_plus_2NGHOSTS")

# Data-driven definitions for 1D loop ranges.
_COORD_DEFINITIONS_1D = {
    "y": {
        # y-axis == { x_mid, z_mid }
        "Cartesian": {
            "i0": [_Nxx_plus_2NGHOSTS[0] / 2],
            "i2": [_Nxx_plus_2NGHOSTS[2] / 2],
        },
        # y-axis == { x_mid, z_mid }
        "Fisheye": {
            "i0": [_Nxx_plus_2NGHOSTS[0] / 2],
            "i2": [_Nxx_plus_2NGHOSTS[2] / 2],
        },
        # y-axis == { theta_mid, see yz-plane discussion for Spherical in 2D for explanation of phi points }
        "Spherical": {
            "i1": [_Nxx_plus_2NGHOSTS[1] / 2],
            "i2": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[2] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[2] - sp.Rational(1, 2),
            ],
        },
        # Cylindrical: rho,phi,z
        # y-axis == { see yz-plane discussion for Spherical in 2D for explanation of phi points, z_mid }
        "Cylindrical": {
            "i1": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[1] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[1] - sp.Rational(1, 2),
            ],
            "i2": [_Nxx_plus_2NGHOSTS[2] / 2],
        },
        # y-axis == { x1_mid, see yz-plane discussion for Spherical in 2D for explanation of phi points }
        "SymTP": {
            "i1": [_Nxx_plus_2NGHOSTS[1] / 2],
            "i2": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[2] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[2] - sp.Rational(1, 2),
            ],
        },
        # NO POINTS ON Y AXIS
        "Wedge": {"i0": [-1], "i1": [-1], "i2": [-1]},
    },
    "z": {
        # z-axis == { x_mid, y_mid }
        "Cartesian": {
            "i0": [_Nxx_plus_2NGHOSTS[0] / 2],
            "i1": [_Nxx_plus_2NGHOSTS[1] / 2],
        },
        # z-axis == { x_mid, y_mid }
        "Fisheye": {
            "i0": [_Nxx_plus_2NGHOSTS[0] / 2],
            "i1": [_Nxx_plus_2NGHOSTS[1] / 2],
        },
        # z-axis == { th_min & th_max, phi_min  }
        "Spherical": {
            "i1": [_NGHOSTS, _Nxx_plus_2NGHOSTS[1] - _NGHOSTS - 1],
            "i2": [_NGHOSTS],
        },
        # NO POINTS ON Z AXIS
        "Spherical_Ring": {"i0": [-1], "i1": [-1], "i2": [-1]},
        # Cylindrical: rho,phi,z
        # z-axis == { rho_min & phi_min }
        "Cylindrical": {"i0": [_NGHOSTS], "i1": [_NGHOSTS]},
        # SymTP:
        # self.xx_to_Cart[2] = f(xx0) * sp.cos(self.xx[1])
        #  -> Aim for cos(xx1) = 1 -> xx1 = 0 & pi
        # z_axis == { xx1_min & xx1_max, xx2_min }
        # FIXME: Missing points between foci (not an easy fix).
        "SymTP": {
            "i1": [_NGHOSTS, _Nxx_plus_2NGHOSTS[1] - _NGHOSTS - 1],
            "i2": [_NGHOSTS],
        },
        # Wedge-like: same as Spherical except x_new = +/-z_old, y_new = y_old, z_new = x_old
        # Thus the z-axis here is the same as the +x-axis in Spherical-like.
        # +x-axis == { theta_mid, phi={phi_mid} (since phi goes from -pi to pi) }
        "Wedge": {"i1": [_Nxx_plus_2NGHOSTS[1] / 2], "i2": [_Nxx_plus_2NGHOSTS[2] / 2]},
    },
}

# Data-driven definitions for 2D loop ranges.
_COORD_DEFINITIONS_2D = {
    "xy": {
        # xy-plane == { z_mid }, where z index is i2
        "Cartesian": {"i2": [_Nxx_plus_2NGHOSTS[2] / 2]},
        "Fisheye": {"i2": [_Nxx_plus_2NGHOSTS[2] / 2]},
        "Cylindrical": {"i2": [_Nxx_plus_2NGHOSTS[2] / 2]},
        # xy-plane == { theta_mid }, where theta index is i1
        "Spherical": {"i1": [_Nxx_plus_2NGHOSTS[1] / 2]},
        "SymTP": {"i1": [_Nxx_plus_2NGHOSTS[1] / 2]},
        # UWedgeHSinhSph: same as Spherical except x_new = -z_old, y_new = y_old, z_new = x_old
        # Thus the xy plane here is the same as the -z,y plane in Spherical-like.
        # LWedgeHSinhSph: same as Spherical except x_new = z_old, y_new = y_old, z_new = -x_old
        # Thus the xy plane here is the same as the z,y plane in Spherical-like
        # (for both see discussion below for yz plane in Spherical)
        "Wedge": {
            "i2": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[2] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[2] - sp.Rational(1, 2),
            ]
        },
    },
    "yz": {
        "Cartesian": {"i0": [_Nxx_plus_2NGHOSTS[0] / 2]},
        "Fisheye": {"i0": [_Nxx_plus_2NGHOSTS[0] / 2]},
        # See documentation for Spherical below; Cylindrical-like coordinates choose xx1 = phi instead of xx2.
        "Cylindrical": {
            "i1": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[1] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[1] - sp.Rational(1, 2),
            ]
        },
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
        "Spherical": {
            "i2": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[2] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[2] - sp.Rational(1, 2),
            ]
        },
        "SymTP": {
            "i2": [
                _NGHOSTS + sp.Rational(1, 4) * _Nxx[2] - sp.Rational(1, 2),
                _NGHOSTS + sp.Rational(3, 4) * _Nxx[2] - sp.Rational(1, 2),
            ]
        },
        # UWedgeHSinhSph: same as Spherical except x_new = -z_old, y_new = y_old, z_new = x_old
        # Thus the yz plane here is the same as the y,x plane in Spherical-like.
        # LWedgeHSinhSph: same as Spherical except x_new = z_old, y_new = y_old, z_new = -x_old
        # Thus the yz plane here is the same as the y,-x plane in Spherical-like
        # xy-plane == { theta_mid }, where theta index is i1
        "Wedge": {"i1": [_Nxx_plus_2NGHOSTS[1] / 2]},
    },
}


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


def compute_1d_loop_ranges(
    CoordSystem: str,
    axis: str,
) -> Tuple[
    Sequence[Union[int, sp.Expr]],
    List[List[Union[int, sp.Expr]]],
    List[Union[int, sp.Expr]],
]:
    r"""
    Compute the grid extents, per‐axis index lists, and point counts needed to extract a one‐dimensional diagnostic slice along the specified axis.

    :param CoordSystem: Specifies the coordinate system (e.g., "Cartesian", "Spherical").
    :param axis: Specifies the axis of output; accepts either "y" or "z".
    :return: Tuple containing:
             - Nxx array,
             - i012_pts list of [i0_pts, i1_pts, i2_pts],
             - numpts list of counts per direction.
    :raises ValueError: If the provided axis is not "y" or "z",
                        or if the CoordSystem is not supported by this function.
    """
    if axis not in ["y", "z"]:
        raise ValueError(
            f"1D loop output only supports y or z axes. axis = {axis} not supported."
        )

    # Select the correct configuration based on coordinate system and axis
    config = None
    # Handle special case of Spherical + Ring on z-axis first
    if axis == "z" and "Spherical" in CoordSystem and "Ring" in CoordSystem:
        config = _COORD_DEFINITIONS_1D["z"]["Spherical_Ring"]
    else:
        for family, family_config in _COORD_DEFINITIONS_1D[axis].items():
            if family in CoordSystem:
                config = family_config
                break

    if config is None:
        raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    i0_pts: List[Union[int, sp.Expr]] = config.get("i0", [])
    i1_pts: List[Union[int, sp.Expr]] = config.get("i1", [])
    i2_pts: List[Union[int, sp.Expr]] = config.get("i2", [])
    i012_pts = [i0_pts, i1_pts, i2_pts]

    # Calculate number of points
    if i0_pts and i0_pts[0] == -1:
        numpts: List[Union[int, sp.Expr]] = [0, 0, 0]
    else:
        numpts = [
            len(i0_pts) if i0_pts else _Nxx[0],
            len(i1_pts) if i1_pts else _Nxx[1],
            len(i2_pts) if i2_pts else _Nxx[2],
        ]

    return _Nxx, i012_pts, numpts


def generate_1d_loop_header(
    axis: str,
    CoordSystem: str,
    numpts: Sequence[Union[int, sp.Expr]],
) -> str:
    """
    Return the C code that defines the i*_pts arrays and data_points buffer for a 1D diagnostic loop.

    :param axis: The axis along which output is generated (e.g., "y" or "z").
    :param CoordSystem: The coordinate system name (e.g., "Cartesian").
    :param numpts: A sequence of three counts [numpts_i0, numpts_i1, numpts_i2].
    :return: The C code header as a string.
    """
    return f"""// Define points for output along the {axis}-axis in {CoordSystem} coordinates.
const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};
int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];

data_point_1d_struct data_points[numpts_i0 * numpts_i1 * numpts_i2];
int data_index = 0;
"""


def append_1d_loop_body(
    out_string: str,
    loop_body_store_results: str,
    axis: str,
    Nxx: Sequence[Union[int, sp.Expr]],
    numpts: Sequence[Union[int, sp.Expr]],
    i012_pts: List[List[Union[int, sp.Expr]]],
    pragma: str,
) -> str:
    r"""
    Append the index-array setup and main LOOP_NOOMP block to out_string.

    :param out_string: Prefix C code to which the loop body is added.
    :param loop_body_store_results: Code to store results in the innermost loop.
    :param axis: The axis name for diagnostic output ("y" or "z").
    :param Nxx: Sequence of full grid sizes per dimension including ghosts.
    :param numpts: Sequence of numbers of points per dimension.
    :param i012_pts: Three lists of explicit index points for each dimension.
    :param pragma: The loop prefix, e.g., OpenMP pragma.
    :return: Combined C code with index setup and main loop appended.

    Examples
    >>> import sympy as sp
    >>> from nrpy.helpers.generic import clang_format
    >>> Nxx0, Nxx1, Nxx2 = sp.symbols("Nxx0 Nxx1 Nxx2")

    >>> # Test Case 1: Full grid loop for the first dimension and explicit points for others.
    >>> # This also tests the conditional comment for 'data_points'.
    >>> params = {
    ...     "out_string": "",
    ...     "loop_body_store_results": "data_points[i0_pt] = xCart[0];",
    ...     "axis": "y",
    ...     "Nxx": [Nxx0, 128, 128],
    ...     "numpts": [Nxx0, 2, 2],
    ...     "i012_pts": [[], [63, 64], [63, 64]],
    ...     "pragma": "#pragma omp parallel for\n"
    ... }
    >>> result = clang_format(append_1d_loop_body(**params))
    >>> print(result)
    #pragma omp parallel for
    for (int i0 = NGHOSTS; i0 < Nxx0 + NGHOSTS; i0++)
      i0_pts[i0 - NGHOSTS] = i0;
    i1_pts[0] = (int)(63);
    i1_pts[1] = (int)(64);
    i2_pts[0] = (int)(63);
    i2_pts[1] = (int)(64);
    // Main loop to store data points along the specified axis
    LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
      const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
      const int idx3 = IDX3(i0, i1, i2);
      REAL xCart[3];
      REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
      xx_to_Cart(params, xOrig, xCart);
    <BLANKLINE>
      data_points[i0_pt] = xCart[0];
    }
    <BLANKLINE>

    >>> # Test Case 2: All dimensions use explicit points, including sympy expressions.
    >>> # Also tests appending to a non-empty out_string.
    >>> params = {
    ...     "out_string": "double time = 1.0;",
    ...     "loop_body_store_results": "results[idx3] = sin(xCart[0]);",
    ...     "axis": "z",
    ...     "Nxx": [64, 64, 64],
    ...     "numpts": [1, 2, 3],
    ...     "i012_pts": [[32], [sp.sympify(31), Nxx1 / 2], [10, 20, 30]],
    ...     "pragma": ""
    ... }
    >>> result = clang_format(append_1d_loop_body(**params))
    >>> print(result)
    double time = 1.0;
    i0_pts[0] = (int)(32);
    i1_pts[0] = (int)(31);
    i1_pts[1] = (int)((1.0 / 2.0) * Nxx1);
    i2_pts[0] = (int)(10);
    i2_pts[1] = (int)(20);
    i2_pts[2] = (int)(30);
    // Main loop
    LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
      const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
      const int idx3 = IDX3(i0, i1, i2);
      REAL xCart[3];
      REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
      xx_to_Cart(params, xOrig, xCart);
    <BLANKLINE>
      results[idx3] = sin(xCart[0]);
    }
    <BLANKLINE>

    >>> # Test Case 3: No points on the specified axis.
    >>> params = {
    ...     "out_string": "",
    ...     "loop_body_store_results": "return;",
    ...     "axis": "z",
    ...     "Nxx": [Nxx0, Nxx1, Nxx2],
    ...     "numpts": [0, 0, 0],
    ...     "i012_pts": [[], [], []],
    ...     "pragma": ""
    ... }
    >>> result = clang_format(append_1d_loop_body(**params))
    >>> print(result)
    // CoordSystem = {CoordSystem} has no points on the z axis!
    // Main loop
    LOOP_NOOMP(i0_pt, 0, numpts_i0, i1_pt, 0, numpts_i1, i2_pt, 0, numpts_i2) {
      const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
      const int idx3 = IDX3(i0, i1, i2);
      REAL xCart[3];
      REAL xOrig[3] = {xx[0][i0], xx[1][i1], xx[2][i2]};
      xx_to_Cart(params, xOrig, xCart);
    <BLANKLINE>
      return;
    }
    <BLANKLINE>
    """
    # Build the i*_pts arrays
    for i in range(3):
        if numpts[i] == Nxx[i]:
            out_string += (
                f"{pragma}"
                f"for(int i{i}=NGHOSTS; i{i}<Nxx{i} + NGHOSTS; i{i}++) "
                f"i{i}_pts[i{i}-NGHOSTS] = i{i};\n"
            )
        elif numpts == [0, 0, 0] and i == 0:
            out_string += (
                f"// CoordSystem = {{CoordSystem}} has no points on the {axis} axis!\n"
            )
        else:
            for j, pt in enumerate(i012_pts[i]):
                out_string += f"i{i}_pts[{j}] = (int)({sp.ccode(pt)});\n"

    # Append the main loop
    out_string += f"""// Main loop{(' to store data points along the specified axis' if 'data_points' in loop_body_store_results else '')}
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  REAL xOrig[3] = {{xx[0][i0], xx[1][i1], xx[2][i2]}};
  xx_to_Cart(params, xOrig, xCart);

  {loop_body_store_results}
}}
"""
    return out_string


def generate_qsort_compare_string() -> str:
    """
    Generate the qsort() comparison function for 1D data_point_1d_struct.

    :return: The C code for the comparison function that sorts by xCart_axis.

    DocTests:
    >>> s = generate_qsort_compare_string()
    >>> "// qsort() comparison function for 1D output." in s
    True
    >>> "static int compare" in s
    True
    """
    return """// qsort() comparison function for 1D output.
static int compare(const void *a, const void *b) {
  REAL l = ((data_point_1d_struct *)a)->xCart_axis;
  REAL r = ((data_point_1d_struct *)b)->xCart_axis;
  return (l > r) - (l < r);
}
"""


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

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Cartesian", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> validate_strings(diag1d, "Cartesian_y_axis")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="SinhSpherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="y")[1])
    >>> validate_strings(diag1d, "SinhSpherical_y_axis")
    >>> diag1d = clang_format(simple_loop_1D(CoordSystem="Spherical", out_quantities_dict = {("REAL", "log10HL"): "log10(fabs(diagnostic_output_gfs[IDX4pt(HGF, idx3)] + 1e-16))"}, axis="z")[1])
    >>> validate_strings(diag1d, "Spherical_z_axis")
    """
    Nxx, i012_pts, numpts = compute_1d_loop_ranges(CoordSystem, axis)

    struct_fields = "\n".join(
        f"  {ctype} {cname};" for ctype, cname in out_quantities_dict
    )
    prefunc_content = f"""// Struct to hold 1D data points
typedef struct {{
  REAL xCart_axis;
{struct_fields}
}} data_point_1d_struct;

{generate_qsort_compare_string()}"""

    out_string = generate_1d_loop_header(axis, CoordSystem, numpts)

    dp1d_assignments = "\n".join(
        f"dp1d.{name} = {expr};" for (_, name), expr in out_quantities_dict.items()
    )
    loop_body_store_results = f"""{{
// Store the data in the data_point_1d_struct
data_point_1d_struct dp1d;
dp1d.xCart_axis = {'xCart[1];' if axis == "y" else 'xCart[2];'}
{dp1d_assignments}
data_points[data_index] = dp1d; data_index++;
}}"""

    out_string = append_1d_loop_body(
        out_string,
        loop_body_store_results,
        axis,
        Nxx,
        numpts,
        i012_pts,
        "#pragma omp parallel for\n",
    )
    out_string = out_string.replace("{CoordSystem}", CoordSystem)

    printf_formats = " ".join(
        "%.15e" if ctype != "int" else "%d" for ctype, _ in out_quantities_dict
    )
    printf_vars = ", ".join(f"data_points[i].{name}" for _, name in out_quantities_dict)

    # Note the \\n to produce a literal \n in the C code. This fixes the doctest failure.
    qsort_and_output_to_file = f"""
qsort(data_points, data_index, sizeof(data_point_1d_struct), compare);

for (int i = 0; i < data_index; i++) {{
  fprintf(outfile, "%.15e {printf_formats}\\n", data_points[i].xCart_axis, {printf_vars});
}}
"""
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

    config = None
    plane_config = _COORD_DEFINITIONS_2D[plane]
    for family, family_config in plane_config.items():
        if family in CoordSystem:
            config = family_config
            break

    if config is None:
        raise ValueError(f"CoordSystem = {CoordSystem} not supported.")

    i0_pts: List[sp.Expr] = config.get("i0", [])
    i1_pts: List[sp.Expr] = config.get("i1", [])
    i2_pts: List[sp.Expr] = config.get("i2", [])
    i012_pts = [i0_pts, i1_pts, i2_pts]

    max_numpts = _Nxx_plus_2NGHOSTS[:] if include_GHOSTS else _Nxx[:]
    numpts: List[Union[int, sp.Expr]] = [
        len(i0_pts) if i0_pts else max_numpts[0],
        len(i1_pts) if i1_pts else max_numpts[1],
        len(i2_pts) if i2_pts else max_numpts[2],
    ]
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

    code_lines = [
        f"// Define points for output along the {plane}-plane in {CoordSystem} coordinates.",
        f"const int numpts_i0={numpts[0]}, numpts_i1={numpts[1]}, numpts_i2={numpts[2]};",
        "int i0_pts[numpts_i0], i1_pts[numpts_i1], i2_pts[numpts_i2];",
    ]
    for i in range(3):
        if numpts[i] == max_numpts[i]:
            code_lines.append(
                f"{pragma}for(int i{i}=NGHOSTS; i{i}<Nxx{i} + NGHOSTS; i{i}++) i{i}_pts[i{i}-NGHOSTS] = i{i};"
            )
        else:
            for j, pt in enumerate(i012_pts[i]):
                code_lines.append(f"i{i}_pts[{j}] = (int)({sp.ccode(pt)});")
    init_code = "\n".join(code_lines)

    declarations = "\n".join(
        [
            f"const {key[0]} {key[1]} = {value};"
            for key, value in out_quantities_dict.items()
        ]
    )
    printf_formats = " ".join(
        ["%.15e" if key[0] != "int" else "%d" for key in out_quantities_dict.keys()]
    )
    printf_vars = ", ".join(key[1] for key in out_quantities_dict.keys())
    coords_to_print = "xCart[0], xCart[1]" if plane == "xy" else "xCart[1], xCart[2]"

    # Note the \\n to produce a literal \n in the C code.
    main_loop = f"""// Main loop to store data points in the specified plane
LOOP_NOOMP(i0_pt,0,numpts_i0, i1_pt,0,numpts_i1, i2_pt,0,numpts_i2) {{
  const int i0 = i0_pts[i0_pt], i1 = i1_pts[i1_pt], i2 = i2_pts[i2_pt];
  const int idx3 = IDX3(i0, i1, i2);
  REAL xCart[3];
  REAL xOrig[3] = {{xx[0][i0], xx[1][i1], xx[2][i2]}};
  xx_to_Cart(params, xOrig, xCart);
  {{
    // Collect diagnostic data
    {declarations}
    fprintf(outfile, "%.15e %.15e {printf_formats}\\n", {coords_to_print}, {printf_vars});
  }}
}}
"""

    return f"{init_code}\n{main_loop}"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

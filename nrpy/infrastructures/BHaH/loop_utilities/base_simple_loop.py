"""
Base loop utilities for generating loops for use within the BHaH infrastructure.

Authors: 
Samuel Tootle
    Email: sdtootle **at** outlook **dot** com
Zachariah B. Etienne
    Email: zachetie **at** gmail **dot** com
Ken Sible
    Email: ksible **at** outlook **dot** com
"""

from typing import List, Tuple
import nrpy.helpers.loop as lp

implemented_loop_regions = ["", "all points", "interior"]


def implemented_loop_regions_err(loop_region: str) -> str:
    """
    Generate the string that is printed when a loop_region is not defined.
    
    :param loop_region: string denoting the intended loop region
    :returns: the error message describing the loop_region was not found
    """
    regions = [f'"{region}"' for region in implemented_loop_regions]
    regions[-1] = f"or {regions[-1]}"
    tmp_str = ",".join(regions)
    return f"loop_region = {loop_region} unsupported. Choose {tmp_str}."


def get_loop_region_ranges(loop_region: str, cuda: bool = False) -> Tuple[List[str], List[str]]:
    """
    Return Loop region index ranges.

    :param loop_region: Loop region
    :param cuda: Toggle whether the loop indices are modified for CUDA compatibility
    :returns: region indicies
    """
    i2i1i0_mins = ["", "", ""]
    i2i1i0_maxs = ["", "", ""]
    prefix = None if not cuda else "tid"

    # 'AllPoints': loop over all points on a numerical grid, including ghost zones
    if loop_region == "all points":
        if not prefix is None:
            i2i1i0_mins = [f"{prefix}{i}" for i in reversed(range(3))]
        else:
            i2i1i0_mins = ["0", "0", "0"]
        i2i1i0_maxs = ["Nxx_plus_2NGHOSTS2", "Nxx_plus_2NGHOSTS1", "Nxx_plus_2NGHOSTS0"]
    # 'InteriorPoints': loop over the interior of a numerical grid, i.e. exclude ghost zones
    elif loop_region == "interior":
        if not prefix is None:
            i2i1i0_mins = [f"{prefix}{i}+NGHOSTS" for i in reversed(range(3))]
        else:
            i2i1i0_mins = ["NGHOSTS", "NGHOSTS", "NGHOSTS"]
        i2i1i0_maxs = [
            "Nxx_plus_2NGHOSTS2 - NGHOSTS",
            "Nxx_plus_2NGHOSTS1 - NGHOSTS",
            "Nxx_plus_2NGHOSTS0 - NGHOSTS",
        ]

    return i2i1i0_mins, i2i1i0_maxs


class base_simple_loop:
    """
    Base class to facilitate generating a simple loop in C (for use inside of a function).

    :param loop_body: Loop body
    :param read_xxs: Read the xx[3][:] 1D coordinate arrays if interior dependency exists
    :param CoordSystem: Coordinate system, e.g., "Cartesian"
    :param enable_rfm_precompute: Enable pre-computation of reference metric
    :param fp_type: Floating point type, e.g., "double".
    :raises ValueError: If `loop_region` is unsupported or if `read_xxs` and `enable_rfm_precompute` are both enabled.
    """

    def __init__(
        self,
        loop_body: str,
        read_xxs: bool = False,
        CoordSystem: str = "Cartesian",
        enable_rfm_precompute: bool = False,
        fp_type: str = "double",
        loop_region: str = "",
        cuda: bool = False,
    ) -> None:
        self.loop_body = loop_body
        self.loop_region = loop_region
        self.read_xxs = read_xxs
        self.CoordSystem = CoordSystem
        self.enable_rfm_precompute = enable_rfm_precompute
        self.fp_type = fp_type

        self.full_loop_body = ""
        if self.loop_region == "":
            self.full_loop_body = self.loop_region
        if self.loop_region not in implemented_loop_regions:
            raise ValueError(implemented_loop_regions_err(self.loop_region))
        self.i2i1i0_mins, self.i2i1i0_maxs = get_loop_region_ranges(
            loop_region, cuda=cuda
        )
        self.prefix_loop_with = [""]

        self.read_rfm_xx_arrays = ["", "", ""]
        # 'Read_xxs': read the xx[3][:] 1D coordinate arrays, as some interior dependency exists
        if self.read_xxs:
            self.read_rfm_xx_arrays = [
                "const REAL xx0 = xx[0][i0];",
                "const REAL xx1 = xx[1][i1];",
                "const REAL xx2 = xx[2][i2];",
            ]

        # 'enable_rfm_precompute': enable pre-computation of reference metric
        if enable_rfm_precompute:
            if read_xxs:
                raise ValueError(
                    "enable_rfm_precompute and Read_xxs cannot both be enabled."
                )
            # pylint: disable=C0415
            from nrpy.infrastructures.BHaH import rfm_precompute

            self.rfmp = rfm_precompute.ReferenceMetricPrecompute(
                CoordSystem, fp_type=fp_type
            )

            self.read_rfm_xx_arrays = [
                self.rfmp.readvr_str[0],
                self.rfmp.readvr_str[1],
                self.rfmp.readvr_str[2],
            ]

        self.increment = ["1", "1", "1"]

    def initialize_based_on__read_rfm_xx_arrays(self) -> None:
        """Set prefix_loop_with based on initialized rfm object."""
        self.prefix_loop_with = [
            "",
            self.read_rfm_xx_arrays[2],
            self.read_rfm_xx_arrays[1],
        ]

        self.loop_body = self.read_rfm_xx_arrays[0] + self.loop_body

    def gen_loop_body(self) -> None:
        """Generate final loop body."""
        self.full_loop_body = str(
            lp.loop(
                ["i2", "i1", "i0"],
                self.i2i1i0_mins,
                self.i2i1i0_maxs,
                self.increment,
                self.prefix_loop_with,
                loop_body=self.loop_body,
            )
        )

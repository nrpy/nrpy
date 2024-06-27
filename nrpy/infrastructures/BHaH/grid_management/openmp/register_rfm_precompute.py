"""
CUDA implementation to register CFunctions for precomputed reference metric infrastructure.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import List

from nrpy.infrastructures.BHaH.grid_management.base_register_rfm_precompute import (
    base_register_CFunctions_rfm_precompute,
)


class register_CFunctions_rfm_precompute(base_register_CFunctions_rfm_precompute):
    """
    Cuda implementation to register C functions for reference metric precomputed lookup arrays.

    :param list_of_CoordSystems: List of coordinate systems to register the C functions.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(
        self, list_of_CoordSystems: List[str], fp_type: str = "double"
    ) -> None:
        super().__init__(list_of_CoordSystems, fp_type=fp_type)

        self.register_CFunction()

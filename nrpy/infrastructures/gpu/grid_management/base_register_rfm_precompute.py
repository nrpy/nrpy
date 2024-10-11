"""
Base to facilitate generating CFunctions for precomputed reference metric infrastructure.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
        Samuel D. Tootle
        sdtootle **at** gmail **dot** com
"""

from typing import Any, Dict, List

import nrpy.c_function as cfc
from nrpy.helpers.generic import superfast_uniq
from nrpy.infrastructures.BHaH.BHaH_defines_h import register_BHaH_defines
from nrpy.infrastructures.BHaH.rfm_precompute import ReferenceMetricPrecompute


class base_register_CFunctions_rfm_precompute:
    """
    Base class to facilitate registering C functions for reference metric precomputed lookup arrays.

    Here we allow for different parallelization schemes to override the default setup before
    registering the Cfunction.

    :param list_of_CoordSystems: List of coordinate systems to register the C functions.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(
        self, list_of_CoordSystems: List[str], fp_type: str = "double"
    ) -> None:

        self.combined_BHaH_defines_list = []
        self.list_of_CoordSystems = list_of_CoordSystems
        self.fp_type = fp_type
        self.includes = ["BHaH_defines.h"]
        self.function_dict: Dict[str, Any] = {}
        self.include_CodeParameters_h = True

        for CoordSystem in list_of_CoordSystems:
            rfm_precompute = ReferenceMetricPrecompute(CoordSystem, fp_type=fp_type)

            for func in [
                ("malloc", rfm_precompute.rfm_struct__malloc),
                ("defines", rfm_precompute.rfm_struct__define),
                ("free", rfm_precompute.rfm_struct__freemem),
            ]:

                desc = f"rfm_precompute_{func[0]}: reference metric precomputed lookup arrays: {func[0]}"
                cfunc_type = "void"
                name = "rfm_precompute_" + func[0]
                params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"
                if func[0] == "defines":
                    params += ", REAL *restrict xx[3]"

                body = " "
                body += func[1]

                self.function_dict[name] = {
                    "desc": desc,
                    "cfunc_type": cfunc_type,
                    "params": params,
                    "body": body,
                    "CoordSystem": CoordSystem,
                    "include_CodeParameters_h": self.include_CodeParameters_h,
                    "prefunc": "",
                }
                defines_list = [
                    s.replace("restrict", "") for s in rfm_precompute.BHaH_defines_list
                ]
                self.combined_BHaH_defines_list.extend(defines_list)

    def populate_BHaH_defines(self) -> None:
        """Add rfm_precompute struct to BHaH Defines."""
        BHaH_defines = "typedef struct __rfmstruct__ {\n"
        BHaH_defines += "\n".join(
            sorted(superfast_uniq(self.combined_BHaH_defines_list))
        )
        BHaH_defines += "\n} rfm_struct;\n"
        register_BHaH_defines("reference_metric", BHaH_defines)

    def register(self) -> None:
        """Register CFunctions."""
        self.populate_BHaH_defines()
        for name, val_dict in self.function_dict.items():
            include_CP = val_dict["include_CodeParameters_h"]
            cfc.register_CFunction(
                prefunc=val_dict["prefunc"],
                includes=self.includes,
                desc=val_dict["desc"],
                cfunc_type=val_dict["cfunc_type"],
                CoordSystem_for_wrapper_func=val_dict["CoordSystem"],
                name=name,
                params=val_dict["params"],
                include_CodeParameters_h=include_CP,
                body=val_dict["body"],
            )

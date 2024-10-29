"""
CUDA implementation to register CFunctions for precomputed reference metric infrastructure.

Authors: Samuel D. Tootle
        sdtootle **at** gmail **dot** com
        Zachariah B. Etienne
        zachetie **at** gmail **dot** com
"""

from typing import List

import nrpy.helpers.gpu_kernels.kernel_base as gputils
from nrpy.helpers.expr_tree import get_unique_expression_symbols
from nrpy.infrastructures.gpu.grid_management.base_register_rfm_precompute import (
    base_register_CFunctions_rfm_precompute,
)
from nrpy.infrastructures.gpu.grid_management.cuda.rfm_precompute import (
    ReferenceMetricPrecompute,
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

        for CoordSystem in list_of_CoordSystems:
            rfm_precompute = ReferenceMetricPrecompute(
                CoordSystem, fp_type=self.fp_type
            )

            for funcs in [
                ("malloc", rfm_precompute.rfm_struct__malloc),
                ("defines", rfm_precompute.rfm_struct__define),
                ("free", rfm_precompute.rfm_struct__freemem),
            ]:

                desc = f"rfm_precompute_{funcs[0]}: reference metric precomputed lookup arrays: {funcs[0]}"
                cfunc_type = "void"
                name = "rfm_precompute_" + funcs[0]
                params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"

                body = " "
                body += funcs[1]

                self.function_dict[name] = {
                    "desc": desc,
                    "cfunc_type": cfunc_type,
                    "params": params,
                    "body": body,
                    "CoordSystem": CoordSystem,
                    "prefunc": "",
                    "include_CodeParameters_h": self.include_CodeParameters_h,
                }

            for func, kernel_dicts in [
                ("defines", rfm_precompute.rfm_struct__define_kernel_dict),
            ]:

                desc = f"rfm_precompute_{func}: reference metric precomputed lookup arrays: {func}"
                cfunc_type = "void"
                name = "rfm_precompute_" + func
                params = "const commondata_struct *restrict commondata, const params_struct *restrict params, rfm_struct *restrict rfmstruct"
                params += ", REAL *restrict xx[3]"

                body = " "
                for i in range(3):
                    body += f"[[maybe_unused]] const REAL *restrict x{i} = xx[{i}];\n"
                prefunc = ""
                prefunc_defs = ""
                for i, (key_sym, kernel_dict) in enumerate(kernel_dicts.items()):
                    # prefunc_defs += f"REAL *restrict {key_sym} = rfmstruct->{key_sym};\n"
                    # These should all be in paramstruct?
                    unique_symbols = get_unique_expression_symbols(
                        kernel_dict["expr"], exclude=[f"xx{j}" for j in range(3)]
                    )
                    kernel_body = ""
                    kernel_body += "// Temporary parameters\n"
                    for sym in unique_symbols:
                        kernel_body += f"const REAL {sym} = d_params[streamid].{sym};\n"
                    kernel_body += kernel_dict["body"]
                    device_kernel = gputils.GPU_Kernel(
                        kernel_body,
                        {
                            f"{key_sym}": "REAL *restrict",
                            f'{kernel_dict["coord"]}': "const REAL *restrict",
                        },
                        f"{name}__{key_sym}_gpu",
                        launch_dict={
                            "blocks_per_grid": [],
                            "threads_per_block": ["32"],
                            "stream": f"(param_streamid + {i}) % nstreams",
                        },
                        fp_type=self.fp_type,
                        comments=f"GPU Kernel to precompute metric quantity {key_sym}.",
                    )
                    prefunc += device_kernel.CFunction.full_function
                    body += "{\n"
                    body += f"REAL *restrict {key_sym} = rfmstruct->{key_sym};\n"
                    body += "const size_t param_streamid = params->grid_idx % nstreams;\n"
                    body += device_kernel.launch_block
                    body += device_kernel.c_function_call().replace("(streamid", "(param_streamid")
                    body += "}\n"
                body = prefunc_defs + body

                self.function_dict[name] = {
                    "desc": desc,
                    "cfunc_type": cfunc_type,
                    "params": params,
                    "body": body,
                    "CoordSystem": CoordSystem,
                    "prefunc": prefunc,
                    "include_CodeParameters_h": True,
                }
        self.register()

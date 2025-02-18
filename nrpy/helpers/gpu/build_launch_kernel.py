"""
Module that abstracts the process of building and launching an host/device execution kernel.

Authors: Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

from typing import Any, Dict, Optional, Tuple

import nrpy.helpers.gpu.cuda_utilities as cuda_utils
from nrpy.helpers.gpu.gpu_kernel import GPU_Kernel, get_params_access


def build_and_launch_kernel(
    kernel_name: str,
    kernel_body: str,
    # Argument dictionaries for GPU vs host code
    arg_dict_cuda: Dict[str, str],
    arg_dict_host: Dict[str, str],
    parallelization: str = "openmp",
    # Optional function signature or other parameters
    cfunc_type: str = "static void",
    comments: str = "",
    # If you need different block/thread dims, pass them in:
    launch_dict: Optional[Dict[str, Any]] = None,
) -> Tuple[str, str]:
    """
    Generate kernels as prefuncs and the necessary launch body.
    Here we build the compute kernel using GPU_Kernel, and then
    append the function definition to `prefunc` and the kernel
    launch code to `body`.

    :param kernel_name: Name of the kernel function.
    :param kernel_body: Actual code inside the kernel.
    :param parallelization: "cuda" or "openmp" or other.
    :param arg_dict_cuda: Dictionary {arg_name: c_type} for CUDA version.
    :param arg_dict_host: Dictionary {arg_name: c_type} for host version.
    :param cfunc_type: e.g. "static void"
    :param comments: Kernel docstring or extra comments
    :param launch_dict: Dictionary to overload CUDA launch settings.
    :return: (prefunc, body) code strings.
    """
    # Prepare return strings
    prefunc = ""
    body = ""

    if parallelization == "cuda":
        params_access = get_params_access(parallelization)
        launch_dict = (
            cuda_utils.default_launch_dictionary if launch_dict is None else launch_dict
        )
        device_kernel = GPU_Kernel(
            kernel_body.replace("params->", params_access),
            arg_dict_cuda,
            kernel_name,
            launch_dict=launch_dict,
            comments=comments,
        )
        # Build the function definition:
        prefunc += device_kernel.CFunction.full_function

        # Build the launch call in `body`:
        body += "{\n"
        body += device_kernel.launch_block
        body += device_kernel.c_function_call()
        body += "}\n"

    else:
        # The "host" path
        device_kernel = GPU_Kernel(
            kernel_body,
            arg_dict_host,
            kernel_name,
            launch_dict=None,
            comments=comments,
            decorators="",
            cuda_check_error=False,
            streamid_param=False,
            cfunc_type=cfunc_type,
        )
        # Build the function definition:
        prefunc += device_kernel.CFunction.full_function

        # Host call is more straightforward:
        body += device_kernel.launch_block
        body += device_kernel.c_function_call()

    return prefunc, body

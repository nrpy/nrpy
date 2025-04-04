"""
Module that abstracts the process of building and launching an host/device execution kernel.

Authors: Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

from typing import Any, Dict, Optional, Tuple

import nrpy.helpers.parallelization.cuda_utilities as cuda_utils
from nrpy.helpers.parallelization.gpu_kernel import GPU_Kernel


def get_params_access(parallelization: str) -> str:
    """
    Return the appropriate params_struct-access prefix for CUDA vs. non-CUDA.
    E.g. 'd_params[streamid].' vs. 'params->' where 'd_params' is
    allocated in __constant__ memory rather a pointer passed as a function argument.

    :param parallelization: The parallelization method to use.
    :returns: The appropriate prefix for accessing the params struct.
    """
    if parallelization == "cuda":
        params_access = "d_params[streamid]."
    else:
        params_access = "params->"
    return params_access


def get_memory_malloc_function(parallelization: str) -> str:
    """
    Return the appropriate function to allocate memory.

    :param parallelization: The parallelization method to use.
    :returns: The appropriate function to allocate memory.
    """
    if parallelization == "cuda":
        malloc_func = "cudaMalloc"
    else:
        malloc_func = "malloc"
    return malloc_func


def get_memory_free_function(parallelization: str) -> str:
    """
    Return the appropriate function to free memory.

    :param parallelization: The parallelization method to use.
    :returns: The appropriate function to free memory.
    """
    if parallelization == "cuda":
        free_func = "cudaFree"
    else:
        free_func = "free"
    return free_func


def get_check_errors_str(
    parallelization: str, kernel_name: str, opt_msg: str = ""
) -> str:
    """
    Return the appropriate function to check for kernel errors.

    :param parallelization: The parallelization method to use.
    :param kernel_name: The name of the kernel function.
    :param opt_msg: Optional custom message to throw.
    :returns: The error checking string.
    """
    opt_msg = f"{kernel_name} failed." if opt_msg == "" else opt_msg

    if parallelization == "cuda":
        check_errors_str = f'cudaCheckErrors({kernel_name}, "{opt_msg}");'
    else:
        check_errors_str = ""
    return check_errors_str


def generate_kernel_and_launch_code(
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
    launchblock_with_braces: bool = False,
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
    :param launchblock_with_braces: If True, wrap the launch block in braces.

    :return: (prefunc, body) code strings.
    """
    # Prepare return strings
    prefunc = ""

    if parallelization == "cuda":
        params_access = get_params_access(parallelization)
        launch_dict = (
            cuda_utils.default_launch_dictionary if launch_dict is None else launch_dict
        )
        device_kernel = GPU_Kernel(
            kernel_body.replace("params->", params_access),
            arg_dict_cuda,
            f"{kernel_name}_gpu",
            launch_dict=launch_dict,
            comments=comments,
            streamid_param="stream" in launch_dict,
        )
        # Build the function definition:
        prefunc += device_kernel.CFunction.full_function

    else:
        # The "host" path
        device_kernel = GPU_Kernel(
            kernel_body,
            arg_dict_host,
            f"{kernel_name}_host",
            launch_dict=None,
            comments=comments,
            decorators="",
            cuda_check_error=False,
            streamid_param=False,
            cfunc_type=cfunc_type,
        )
        # Build the function definition:
        prefunc += device_kernel.CFunction.full_function

    # Build the launch call in `body`:
    body = device_kernel.launch_block + device_kernel.c_function_call()
    if launchblock_with_braces and len(body.splitlines()) > 1:
        body = f"""{{{body}
        }}"""

    return prefunc, body


def get_loop_parameters(
    parallelization: str, dim: int = 3, enable_intrinsics: bool = False
) -> str:
    """
    Return the appropriate loop parameters for CUDA vs. non-CUDA.

    :param parallelization: The parallelization method to use.
    :param dim: The number of dimensions to loop over.
    :param enable_intrinsics: Whether to modify str based on hardware intrinsics.
    :returns: The appropriate loop parameters.
    """
    loop_params = ""
    param_access = get_params_access(parallelization)
    for i in range(dim):
        loop_params += f"MAYBE_UNUSED const int Nxx_plus_2NGHOSTS{i} = {param_access}Nxx_plus_2NGHOSTS{i};\n"
    loop_params += "\n"

    for i in range(dim):
        loop_params += (
            f"MAYBE_UNUSED const REAL invdxx{i} = {param_access}invdxx{i};\n"
            if not enable_intrinsics
            else f"const REAL NOSIMDinvdxx{i} = {param_access}invdxx{i};\n"
            f"MAYBE_UNUSED const REAL_SIMD_ARRAY invdxx{i} = ConstSIMD(NOSIMDinvdxx{i});\n"
        )
    loop_params += "\n"

    if parallelization == "cuda":
        for i, coord in zip(range(dim), ["x", "y", "z"]):
            loop_params += f"MAYBE_UNUSED const int tid{i}  = blockIdx.{coord} * blockDim.{coord} + threadIdx.{coord};\n"
        loop_params += "\n"

        for i, coord in zip(range(dim), ["x", "y", "z"]):
            loop_params += f"MAYBE_UNUSED const int stride{i}  = blockDim.{coord} * gridDim.{coord};\n"
        loop_params += "\n"
        loop_params = loop_params.replace("SIMD", "CUDA")

    return loop_params

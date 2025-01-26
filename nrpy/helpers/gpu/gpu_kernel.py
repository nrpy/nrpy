"""
Module that provides the base functionality for generating GPU Kernels.

Authors: Samuel D. Tootle; sdtootle **at** gmail **dot** com
"""

from typing import Any, Dict, Union

import nrpy.c_function as cfc


class GPU_Kernel:
    """
    Class to Generate GPU Kernel code.

    While the current implementation is tuned with CUDA in mind, flexibility can easily
    be added or obtained by inheriting and overloading the implementation specific
    syntax.

    :param body: Kernel body
    :param params_dict: Dictionary storing function arguments as keys and types as Dictionary entry
    :param c_function_name: Kernel function name
    :param cfunc_type: C Function return type
    :param decorators: Function decorators i.e. Kernel type, templates, etc
    :param fp_type: Floating point type, i.e. double, float, long double
    :param comments: Additional comments to add to Function description
    :param launch_dict: Dictionary that stores kernel launch settings
    :param streamid_param: Toggle whether streamid is a kernel argument parameter
    :param cuda_check_error: Add CUDA error checking after kernel call. Default is True

    >>> kernel = GPU_Kernel(
    ... "*x = in;",
    ... {'x' : 'REAL *restrict', 'in' : 'const REAL'},
    ... 'basic_assignment_gpu',
    ... launch_dict = {
    ... 'blocks_per_grid' : [32],
    ... 'threads_per_block' : [128,28,1],
    ... },
    ... streamid_param = False,
    ... )
    >>> print(kernel.c_function_call())
    basic_assignment_gpu<<<blocks_per_grid,threads_per_block>>>(x, in);
    cudaCheckErrors(cudaKernel, "basic_assignment_gpu failure");
    <BLANKLINE>
    >>> print(kernel.CFunction.full_function)
    /**
     * GPU Kernel: basic_assignment_gpu.
     *
     */
    __global__ static void basic_assignment_gpu(REAL *restrict x, const REAL in) { *x = in; } // END FUNCTION basic_assignment_gpu
    <BLANKLINE>
    >>> print(kernel.launch_block)
    <BLANKLINE>
    const size_t threads_in_x_dir = 128;
    const size_t threads_in_y_dir = 28;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid(32,1,1);
    <BLANKLINE>
    >>> kernel = GPU_Kernel(
    ... "*x = in;",
    ... {'x' : 'REAL *restrict', 'in' : 'const REAL'},
    ... 'basic_assignment_gpu',
    ... launch_dict = {
    ... 'blocks_per_grid' : [],
    ... 'threads_per_block' : [128,28,1],
    ... },
    ... )
    >>> print(kernel.launch_block)
    <BLANKLINE>
    const size_t threads_in_x_dir = 128;
    const size_t threads_in_y_dir = 28;
    const size_t threads_in_z_dir = 1;
    dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);
    dim3 blocks_per_grid(
        (Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
        (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
        (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir
    );
    <BLANKLINE>
    """

    def __init__(
        self,
        body: str,
        params_dict: Dict[str, Any],
        c_function_name: str,
        cfunc_type: str = "static void",
        decorators: str = "__global__",
        fp_type: str = "double",
        comments: str = "",
        launch_dict: Union[Dict[str, Any], None] = None,
        streamid_param: bool = True,
        cuda_check_error: bool = True,
    ) -> None:
        self.body = body
        self.decorators = decorators
        self.params_dict = params_dict
        if self.decorators != "__host__" and streamid_param:
            self.params_dict = {"streamid": "const size_t", **params_dict}
        self.name = c_function_name
        self.cfunc_type = f"{decorators} {cfunc_type}"
        self.fp_type = fp_type
        self.cuda_check_error = cuda_check_error

        self.CFunction: cfc.CFunction
        self.desc: str = f"GPU Kernel: {self.name}.\n" + comments
        self.launch_dict = launch_dict
        self.launch_block: str = ""
        self.launch_settings: str = "("

        if self.decorators == "__global__" and launch_dict is None:
            raise ValueError(f"Error: {self.decorators} requires a launch_dict")
        self.generate_launch_block()

        self.param_list = [f"{v} {k}" for k, v in self.params_dict.items()]
        # Store CFunction
        self.CFunction = cfc.CFunction(
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=",".join(self.param_list),
            body=self.body,
        )

    def generate_launch_block(self) -> None:
        """Generate preceding launch block definitions for kernel function call."""
        if self.launch_dict is None:
            return
        threads_per_block = self.launch_dict["threads_per_block"]
        for _ in range(3 - len(threads_per_block)):
            threads_per_block += ["1"]
        block_def_str = f"""
const size_t threads_in_x_dir = {threads_per_block[0]};
const size_t threads_in_y_dir = {threads_per_block[1]};
const size_t threads_in_z_dir = {threads_per_block[2]};
dim3 threads_per_block(threads_in_x_dir, threads_in_y_dir, threads_in_z_dir);"""

        blocks_per_grid = self.launch_dict["blocks_per_grid"]
        if len(blocks_per_grid) > 0:
            for _ in range(3 - len(blocks_per_grid)):
                blocks_per_grid += [1]
            blocks_per_grid_str = ",".join(map(str, blocks_per_grid))
            grid_def_str = f"dim3 blocks_per_grid({blocks_per_grid_str});"
        else:
            grid_def_str = """dim3 blocks_per_grid(
    (Nxx_plus_2NGHOSTS0 + threads_in_x_dir - 1) / threads_in_x_dir,
    (Nxx_plus_2NGHOSTS1 + threads_in_y_dir - 1) / threads_in_y_dir,
    (Nxx_plus_2NGHOSTS2 + threads_in_z_dir - 1) / threads_in_z_dir
);"""

        # Determine if the stream needs to be added to launch
        stream_def_str = None
        if "stream" in self.launch_dict:
            if (
                self.launch_dict["stream"] == ""
                or self.launch_dict["stream"] == "default"
            ):
                stream_def_str = "size_t streamid = params->grid_idx % nstreams;"
            else:
                stream_def_str = f"size_t streamid = {self.launch_dict['stream']};"

        # Determine if the shared memory size needs to be added to launch
        # If a stream is specified, we need to at least set SM to 0
        sm_def_str = None
        if "sm" in self.launch_dict or not stream_def_str is None:
            if (
                not "sm" in self.launch_dict
                or self.launch_dict["sm"] == ""
                or self.launch_dict["sm"] == "default"
            ):
                sm_def_str = "size_t sm = 0;"
                self.launch_dict["sm"] = 0
            else:
                sm_def_str = f"size_t sm = {self.launch_dict['sm']};"

        self.launch_block = f"""{block_def_str}
{grid_def_str}
"""
        if not sm_def_str is None:
            self.launch_block += f"{sm_def_str}"
        if not stream_def_str is None:
            self.launch_block += f"{stream_def_str}"

        self.launch_settings = "<<<blocks_per_grid,threads_per_block"
        if not sm_def_str is None:
            self.launch_settings += ",sm"
        if not stream_def_str is None:
            self.launch_settings += ",streams[streamid]"
        self.launch_settings += ">>>("

    def c_function_call(self) -> str:
        """
        Generate the C function call for a given Kernel.

        :return: The C function call as a string.
        """
        c_function_call: str = self.name
        c_function_call += self.launch_settings

        for p in self.params_dict:
            c_function_call += f"{p}, "
        c_function_call = c_function_call[:-2] + ");\n"

        if self.cuda_check_error:
            msg = f"{self.name} failure"
            msg = f'cudaCheckErrors(cudaKernel, "{msg}")'
            c_function_call += f"{msg};\n"

        return c_function_call


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

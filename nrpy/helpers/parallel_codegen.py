"""
This NRPy+ helper module provides functionality for
 parallel code generation by registering and
 calling functions concurrently.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Callable, Dict, Optional, Tuple, Union, cast
from importlib import import_module
from concurrent.futures import ProcessPoolExecutor, as_completed
import nrpy.grid as gri
import nrpy.c_function as cfc
import nrpy.params as par


par.register_param(bool, __name__, "parallel_codegen_enable", False)
par.register_param(str, __name__, "parallel_codegen_stage", "register")


class ParallelCodeGen:
    def __init__(self, path: str, args: Dict[str, Any]):
        """
        Initializes a ParallelCodeGen object.

        :param path: Dot-separated path to the function.
        :param args: Dictionary containing arguments to pass to the function.
        """
        self.module_path = path.rsplit(".", 1)[0]
        self.function_name = path.rsplit(".", 1)[1]
        self.function_args = args


# Contains a dictionary of ParallelCodeGen objects
ParallelCodeGen_dict: Dict[str, ParallelCodeGen] = {}


NRPyEnv_type = Tuple[
    Dict[str, par.NRPyParameter],
    Dict[str, par.CodeParameter],
    Dict[str, cfc.CFunction],
    Dict[
        str,
        Union[
            gri.GridFunction,
            gri.BHaHGridFunction,
            gri.BaseETGridFunction,
            gri.CarpetXGridFunction,
        ],
    ],
]


def NRPyEnv() -> NRPyEnv_type:
    """
    :return: Tuple containing various global dictionaries
    """
    return (
        par.glb_params_dict,
        par.glb_code_params_dict,
        cfc.CFunction_dict,
        gri.glb_gridfcs_dict,
    )


def pcg_registration_phase() -> bool:
    """
    :return: Boolean indicating if the parallel code generation registration phase is active
    """
    return (
        cast(bool, par.parval_from_str("parallel_codegen_enable"))
        and par.parval_from_str("parallel_codegen_stage") == "register"
    )


def register_func_call(name: str, args: Dict[str, Any]) -> None:
    """
    Registers a function call if the parallel code generation registration phase is active.

    :param name: Name of the function
    :param args: Arguments to pass to the function
    :return: Boolean indicating if the function call was registered
    """
    if name + str(args) in ParallelCodeGen_dict:
        raise ValueError(f"Already registered {name + str(args)}.")
    ParallelCodeGen_dict[name + str(args)] = ParallelCodeGen(name, args)
    return None


def unpack_NRPy_environment_dict(NRPy_environment_dict: Dict[str, NRPyEnv_type]):
    for env in NRPy_environment_dict.values():
        par.glb_params_dict.update(env[0])
        par.glb_code_params_dict.update(env[1])
        cfc.CFunction_dict.update(env[2])
        gri.glb_gridfcs_dict.update(env[3])


def get_nested_function(module_path: str, function_name: str) -> Optional[Callable]:
    """
    Retrieves a nested function from a specified Python module.

    :param module_path: The dot-separated path to the Python module.
    :param function_name: The dot-separated path to the nested function within the module.
    :return: The nested function if found, otherwise None.
    """
    try:
        module = import_module(module_path)
    except ImportError as e:
        raise ImportError(f"Module could not be imported: {e}")

    function_parts = function_name.split(".")
    try:
        nested_obj = module
        for part in function_parts:
            nested_obj = getattr(nested_obj, part)
    except AttributeError as e:
        raise AttributeError(f"Error accessing nested object: {e}")

    if callable(nested_obj):
        return nested_obj
    else:
        raise TypeError(
            f"The specified path {function_name} did not lead to a callable function."
        )


def parallel_function_call(PCG: Any) -> NRPyEnv_type:
    """
    Calls the registered function specified by the given key and ParallelCodeGen object.

    :param key: The key corresponding to the function to be called.
    :param PCG: The ParallelCodeGen object containing function details.
    :param NRPy_environment_to_unpack: Dictionary for storing results of the function calls.
    """
    try:
        module_path = PCG.module_path
        function_name = PCG.function_name
        function_args = PCG.function_args

        function_to_call = get_nested_function(module_path, function_name)

        if function_to_call is None:
            raise KeyError(f"Function '{function_name}' not found.")

        return function_to_call(**function_args)

    except Exception as ex:
        import traceback

        tb_str = traceback.format_exception(
            etype=type(ex), value=ex, tb=ex.__traceback__
        )
        raise Exception(
            f"An error occurred while calling the function: {ex}\n Traceback: {''.join(tb_str)}"
        )


def do_parallel_codegen() -> None:
    """
    Performs parallel code generation by calling registered functions concurrently.
    """
    par.set_parval_from_str("parallel_codegen_stage", "codegen")

    NRPy_environment_to_unpack: Dict[str, Any] = {}

    with ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(parallel_function_call, value): key
            for key, value in ParallelCodeGen_dict.items()
        }

        for future in as_completed(futures):
            key = futures[future]
            try:
                NRPy_environment_to_unpack[key] = future.result()
                funcname_args = ParallelCodeGen_dict[key].function_name
                print(f"Worker associated with function '{funcname_args}' is done.")
            except Exception as e:
                raise Exception(
                    f"An error occurred in the process associated with key '{key}':\n {e}"
                )

    unpack_NRPy_environment_dict(NRPy_environment_to_unpack)

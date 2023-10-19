"""
Core functions that enable registering and calling functions in parallel.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Callable, Dict, Tuple, Union, cast
from importlib import import_module
import time

from multiprocess import Pool, Manager  # type: ignore # pylint: disable=E0611
import nrpy.grid as gri
import nrpy.c_function as cfc
import nrpy.params as par

par.register_param(bool, __name__, "parallel_codegen_enable", False)
par.register_param(str, __name__, "parallel_codegen_stage", "register")


class ParallelCodeGen:
    """Stores necessary information to call a function in parallel."""

    def __init__(self, path: str, args: Dict[str, Any]) -> None:
        """
        Initialize a ParallelCodeGen object.

        :param path: Path to the function.
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
            gri.ETLegacyGridFunction,
            gri.CarpetXGridFunction,
        ],
    ],
    Dict[str, Dict[str, Any]],
]


def NRPyEnv() -> NRPyEnv_type:
    """
    Retrieve a tuple containing various global dictionaries.

    :return: Tuple containing various global dictionaries.
    """
    return (
        par.glb_params_dict,
        par.glb_code_params_dict,
        cfc.CFunction_dict,
        gri.glb_gridfcs_dict,
        par.glb_extras_dict,
    )


def deep_update(d: Dict[Any, Any], u: Dict[Any, Any]) -> None:
    """
    Perform a deep update on a dictionary.

    :param d: The original dictionary to update.
    :param u: The dictionary containing new keys and values.

    Doctest:
    >>> original = {'a': 1, 'b': {'c': 2}}
    >>> new = {'a': 'new_a', 'b': {'d': 'new_d'}}
    >>> deep_update(original, new)
    >>> original == {'a': 'new_a', 'b': {'c': 2, 'd': 'new_d'}}
    True
    """
    for k, v in u.items():
        if isinstance(v, dict):
            d[k] = d.get(k, {})
            deep_update(d[k], v)
        else:
            d[k] = v


def unpack_NRPy_environment_dict(
    NRPy_environment_dict: Dict[str, NRPyEnv_type]
) -> None:
    """
    Unpack the NRPy environment dictionaries.

    :param NRPy_environment_dict: Dictionary containing NRPy environment types.
    """
    for env in NRPy_environment_dict.values():
        par.glb_params_dict.update(env[0])
        par.glb_code_params_dict.update(env[1])
        cfc.CFunction_dict.update(env[2])
        gri.glb_gridfcs_dict.update(env[3])
        deep_update(par.glb_extras_dict, env[4])


def pcg_registration_phase() -> bool:
    """
    Determine if the parallel code generation registration phase is active.

    :return: Boolean indicating if the registration phase is active.
    """
    return (
        cast(bool, par.parval_from_str("parallel_codegen_enable"))
        and par.parval_from_str("parallel_codegen_stage") == "register"
    )


def register_func_call(name: str, args: Dict[str, Any]) -> None:
    """
    Register a function call if the registration phase is active.

    :param name: Name of the function.
    :param args: Arguments to pass to the function.
    """
    if name + str(args) in ParallelCodeGen_dict:
        raise ValueError(f"Already registered {name + str(args)}.")

    ParallelCodeGen_dict[name + str(args)] = ParallelCodeGen(name, args)


def get_nested_function(
    module_path: str, function_name: str
) -> Callable[..., NRPyEnv_type]:
    """
    Retrieve a nested function from a specified Python module.

    :param module_path: The dot-separated path to the Python module.
    :param function_name: The dot-separated path to the nested function within the module.
    :return: The nested function if found.
    """
    try:
        module = import_module(module_path)
    except ImportError as e:
        raise ImportError(f"Module could not be imported: {e}") from e

    function_parts = function_name.split(".")
    try:
        nested_obj = module
        for part in function_parts:
            nested_obj = getattr(nested_obj, part)
    except AttributeError as e:
        raise AttributeError(f"Error accessing nested object: {e}") from e

    if not callable(nested_obj):
        raise TypeError(
            f"The specified path {function_name} did not lead to a callable function."
        )
    return nested_obj


def parallel_function_call(PCG: Any) -> NRPyEnv_type:
    """
    Call the registered function specified by the given ParallelCodeGen object.

    :param PCG: The ParallelCodeGen object containing function details.
    :return: The result of the function call, packed as NRPyEnv_type.
    """
    try:
        module_path = PCG.module_path
        function_name = PCG.function_name
        function_args = PCG.function_args

        function_to_call = get_nested_function(module_path, function_name)

        return function_to_call(**function_args)

    except (ImportError, AttributeError, TypeError) as ex:
        raise RuntimeError(
            f"An error occurred while calling the function: {ex}"
        ) from ex


def wrapper_func(args: Tuple[Dict[str, Any], str, Any]) -> Any:
    """
    Execute a given function in parallel, wrapping its call for error-handling and performance logging.

    This function serves as a bridge for parallel processing, ensuring that each task
    reports its completion time and handles any potential exceptions that may arise during execution.

    :param args: A tuple containing the shared dictionary, key, and value for each task.
    :return: The key and the result of the parallel_function_call.
    :raises RuntimeError: If any exception occurs during the task's execution.
    """
    shared_dict, key, value = args
    start_time = time.time()
    try:
        result = parallel_function_call(value)
        shared_dict[key] = result
        funcname_args = value.function_name
        print(
            f"In {(time.time()-start_time):.3f}s, worker completed task '{funcname_args}'"
        )
        return key, result
    except Exception as e:
        raise RuntimeError(
            f"An error occurred in the process associated with key '{key}':\n {e}"
        ) from e


def do_parallel_codegen() -> None:
    """Perform parallel code generation by calling registered functions concurrently."""
    if not par.parval_from_str("parallel_codegen_enable"):
        return

    par.set_parval_from_str("parallel_codegen_stage", "codegen")

    manager = Manager()
    NRPy_environment_to_unpack: Dict[str, Any] = manager.dict()

    with Pool() as pool:
        pool.map(
            wrapper_func,
            [
                (NRPy_environment_to_unpack, key, value)
                for key, value in ParallelCodeGen_dict.items()
            ],
        )

    unpack_NRPy_environment_dict(dict(NRPy_environment_to_unpack))


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

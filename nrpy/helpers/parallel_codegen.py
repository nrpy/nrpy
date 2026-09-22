# nrpy/helpers/parallel_codegen.py
"""
Core functions that enable registering and calling functions in parallel.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import logging
import time
from copy import deepcopy
from importlib import import_module
from multiprocessing import Pool, set_start_method
from typing import Any, Callable, Dict, List, Tuple, Union, cast

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
import nrpy.py_function as pyfc

logging.basicConfig(level=logging.DEBUG, format="%(message)s")
par.register_param(bool, __name__, "enable_parallel_codegen", False)
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
    Dict[str, pyfc.PyFunction],
    Dict[
        str,
        Union[
            gri.GridFunction,
            gri.BHaHGridFunction,
            gri.ETLegacyGridFunction,
            gri.CarpetXGridFunction,
            gri.DendroGridFunction,
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
        pyfc.PyFunction_dict,
        gri.glb_gridfcs_dict,
        par.glb_extras_dict,
    )


def deep_update(d: Dict[Any, Any], u: Dict[Any, Any]) -> None:
    """
    Perform a deep update on a dictionary.

    :param d: The original dictionary to update.
    :param u: The dictionary containing new keys and values.

    Doctests:
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


def _merge_registered_definitions(
    registry_name: str,
    destination: Dict[str, Any],
    additions: Dict[str, Any],
) -> None:
    """
    Merge one typed NRPy registry, rejecting conflicting definitions.

    :param registry_name: Human-readable registry type for error reporting.
    :param destination: Registry copy receiving validated definitions.
    :param additions: Definitions returned by one worker.
    :raises ValueError: If a name has two different definitions.
    """
    for name in sorted(additions):
        if name in destination and not (
            type(destination[name]) is type(additions[name])
            and vars(destination[name]) == vars(additions[name])
        ):
            raise ValueError(
                f"Parallel code generation produced conflicting {registry_name} "
                f"definitions for '{name}'."
            )
        destination[name] = additions[name]


def _merge_extras(destination: Dict[str, Any], additions: Dict[str, Any]) -> None:
    """
    Merge nested NRPy extras while rejecting conflicting leaf values.

    :param destination: Extras copy receiving validated values.
    :param additions: Extras returned by one worker.
    :raises ValueError: If a leaf name has two different values.
    """
    for name in sorted(additions):
        value = additions[name]
        if name not in destination:
            destination[name] = value
        elif isinstance(destination[name], dict) and isinstance(value, dict):
            _merge_extras(destination[name], value)
        elif destination[name] != value:
            raise ValueError(
                "Parallel code generation produced conflicting extras "
                f"definitions for '{name}'."
            )


def unpack_NRPy_environment_dict(
    NRPy_environment_dict: Dict[str, NRPyEnv_type],
) -> None:
    """
    Unpack the NRPy environment dictionaries.

    :param NRPy_environment_dict: Dictionary containing NRPy environment types.
    """
    merged_params = deepcopy(par.glb_params_dict)
    merged_code_params = deepcopy(par.glb_code_params_dict)
    merged_cfunctions = deepcopy(cfc.CFunction_dict)
    merged_pyfunctions = deepcopy(pyfc.PyFunction_dict)
    merged_gridfunctions = deepcopy(gri.glb_gridfcs_dict)
    merged_extras = deepcopy(par.glb_extras_dict)

    for task_key in sorted(NRPy_environment_dict):
        env = NRPy_environment_dict[task_key]
        _merge_registered_definitions("NRPyParameter", merged_params, env[0])
        _merge_registered_definitions("CodeParameter", merged_code_params, env[1])
        _merge_registered_definitions("CFunction", merged_cfunctions, env[2])
        _merge_registered_definitions("PyFunction", merged_pyfunctions, env[3])
        _merge_registered_definitions("gridfunction", merged_gridfunctions, env[4])
        _merge_extras(merged_extras, env[5])

    par.glb_params_dict.clear()
    par.glb_params_dict.update(merged_params)
    par.glb_code_params_dict.clear()
    par.glb_code_params_dict.update(merged_code_params)
    cfc.CFunction_dict.clear()
    cfc.CFunction_dict.update(merged_cfunctions)
    pyfc.PyFunction_dict.clear()
    pyfc.PyFunction_dict.update(merged_pyfunctions)
    gri.glb_gridfcs_dict.clear()
    gri.glb_gridfcs_dict.update(merged_gridfunctions)
    par.glb_extras_dict.clear()
    par.glb_extras_dict.update(merged_extras)


def pcg_registration_phase() -> bool:
    """
    Determine if the parallel code generation registration phase is active.

    :return: Boolean indicating if the registration phase is active.
    """
    return (
        cast(bool, par.parval_from_str("enable_parallel_codegen"))
        and par.parval_from_str("parallel_codegen_stage") == "register"
    )


def register_func_call(name: str, args: Dict[str, Any]) -> None:
    """
    Register a function call if the registration phase is active.

    :param name: Name of the function.
    :param args: Arguments to pass to the function.

    """
    task_key = f"{len(ParallelCodeGen_dict):08d}:{name}"
    ParallelCodeGen_dict[task_key] = ParallelCodeGen(name, args)


def get_nested_function(
    module_path: str, function_name: str
) -> Callable[..., NRPyEnv_type]:
    """
    Retrieve a nested function from a specified Python module.

    This function dynamically imports a Python module using its dot-separated path and
    then navigates through the module's attributes to find a nested function or callable
    object based on a dot-separated function name. If the module cannot be imported, an
    AttributeError occurs during navigation, or the final object is not callable, appropriate
    exceptions are raised.

    :param module_path: The dot-separated path to the Python module.
    :param function_name: The dot-separated path to the nested function within the module.
    :return: The nested function if found.
    :raises ImportError: If the module cannot be imported.
    :raises AttributeError: If an attribute error occurs during navigation to the nested function.
    :raises TypeError: If the specified path does not lead to a callable function.
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
    :raises RuntimeError: If an error occurs during the dynamic retrieval or the call of the function.
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


def wrapper_func(args: Tuple[str, Any]) -> Tuple[str, NRPyEnv_type]:
    """
    Execute a given function in parallel, wrapping its call for error-handling and performance logging.

    This wrapper calls one generated-code task in a worker process, records its
    completion time, and reports exceptions from that task.

    :param args: A tuple containing the shared dictionary, key, and value for each task.
    :return: The key and the result of the parallel_function_call.
    :raises RuntimeError: If any exception occurs during the task's execution.
    """
    key, value = args
    start_time = time.time()
    try:
        # logging.debug(f"Starting task with key: {key}")
        result = parallel_function_call(value)
        # logging.debug(f"Task {key} completed: {result}")
        funcname_args = value.function_name
        elapsed_time = time.time() - start_time
        logging.info(
            "In %.3fs, worker completed task '%s'", elapsed_time, funcname_args
        )
        return key, result
    except Exception as e:
        logging.exception(
            "An error occurred in the process associated with key '%s':", key
        )
        raise RuntimeError(
            f"An error occurred in the process associated with key '{key}':\n {e}"
        ) from e


def do_parallel_codegen() -> None:
    """
    Perform parallel code generation by calling registered functions concurrently.

    Worker and registry-merge failures propagate after restoring the
    code-generation stage. Parent registries remain unchanged on either failure.
    """
    if not par.parval_from_str("enable_parallel_codegen"):
        return

    tasks = sorted(ParallelCodeGen_dict.items())
    if not tasks:
        return
    previous_stage = par.parval_from_str("parallel_codegen_stage")
    par.set_parval_from_str("parallel_codegen_stage", "codegen")

    # By default, MacOS adopts the "spawn" method, which breaks the global environment.
    #   Linux, on the other hand, uses the "fork" method, which preserves the environment.
    #   Luckily MacOS does support "fork", though Apple advises against its use because fork
    #   without an exec step can cause weird breakage on macOS with certain libraries.
    #   Empirically, for NRPy, this seems to work fine.
    # Note that prior to using multiprocessing, we used the `multiprocess` library,
    #   which was very janky -- causing all sorts of race conditions.
    set_start_method("fork", force=True)
    # Starting more processes than independent tasks only adds fork, import,
    # and shutdown cost.  Large symbolic kernels can consume substantial memory,
    # so one worker per disjoint task is also the safe upper bound.
    try:
        with Pool(processes=len(tasks)) as pool:
            worker_results: List[Tuple[str, NRPyEnv_type]] = pool.map(
                wrapper_func,
                tasks,
            )
        unpack_NRPy_environment_dict(dict(worker_results))
    finally:
        par.set_parval_from_str("parallel_codegen_stage", previous_stage)
    ParallelCodeGen_dict.clear()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

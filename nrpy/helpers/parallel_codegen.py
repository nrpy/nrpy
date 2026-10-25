# nrpy/helpers/parallel_codegen.py
"""
Core functions that enable registering and calling functions in parallel.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import logging
import pickle
import time
from importlib import import_module
from multiprocessing import Pool, set_start_method
from typing import Any, Callable, Dict, Tuple, Union, cast

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


# NRPyEnv_type fixes a value type per dictionary, but the two helpers below treat all six
# dictionaries alike, so they use Dict[Any, Any].
def pickle_environment(env: NRPyEnv_type) -> Tuple[Dict[Any, Any], ...]:
    """
    Record the current contents of the NRPy environment dictionaries as pickled bytes.

    :param env: NRPy environment dictionaries.
    :return: One dictionary per input dictionary, holding the pickled value of each entry (nested dictionaries stay dictionaries).
    """

    def pickled(d: Dict[Any, Any]) -> Dict[Any, Any]:
        return {
            k: pickled(v) if isinstance(v, dict) else pickle.dumps(v)
            for k, v in d.items()
        }

    return tuple(pickled(d) for d in env)


def changed_environment(
    env: NRPyEnv_type, before: Tuple[Dict[Any, Any], ...]
) -> NRPyEnv_type:
    """
    Return the entries of the NRPy environment dictionaries that differ from an earlier record.

    :param env: NRPy environment dictionaries.
    :param before: Record of the same dictionaries from pickle_environment().
    :return: NRPy environment dictionaries holding only new or modified entries (only the extras dictionary is reduced recursively).
    """
    current: Tuple[Dict[Any, Any], ...] = env

    def changed(
        d: Dict[Any, Any], old: Dict[Any, Any], recurse: bool
    ) -> Dict[Any, Any]:
        result: Dict[Any, Any] = {}
        for k, v in d.items():
            old_v = old.get(k)
            if recurse and isinstance(v, dict):
                nested = changed(v, old_v if isinstance(old_v, dict) else {}, True)
                if nested or not isinstance(old_v, dict):
                    result[k] = nested
            elif old_v is None or isinstance(old_v, dict) or pickle.dumps(v) != old_v:
                result[k] = v
        return result

    return cast(
        NRPyEnv_type,
        tuple(changed(d, b, n == 5) for n, (d, b) in enumerate(zip(current, before))),
    )


def unpack_NRPy_environment_dict(
    NRPy_environment_dict: Dict[str, NRPyEnv_type],
) -> None:
    """
    Unpack the NRPy environment dictionaries, applying them in the order of the dictionary.

    :param NRPy_environment_dict: Dictionary containing NRPy environment types.
    """
    for env in NRPy_environment_dict.values():
        par.glb_params_dict.update(env[0])
        par.glb_code_params_dict.update(env[1])
        cfc.CFunction_dict.update(env[2])
        pyfc.PyFunction_dict.update(env[3])
        gri.glb_gridfcs_dict.update(env[4])
        deep_update(par.glb_extras_dict, env[5])


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

    :raises ValueError: If a function call with the same name and arguments has already been registered.
    """
    if name + str(args) in ParallelCodeGen_dict:
        raise ValueError(f"Already registered {name + str(args)}.")

    ParallelCodeGen_dict[name + str(args)] = ParallelCodeGen(name, args)


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


def wrapper_func(args: Tuple[str, ParallelCodeGen]) -> Tuple[str, NRPyEnv_type]:
    """
    Execute a given function in parallel, wrapping its call for error-handling and performance logging.

    This wrapper calls one generated-code task in a worker process, records its
    completion time, and reports exceptions from that task. A worker process runs
    several tasks in turn, so the result holds only the entries that this task
    added or modified.

    :param args: A tuple containing the key and value for each task.
    :return: The key and the entries of the NRPy environment changed by the parallel_function_call.
    :raises RuntimeError: If any exception occurs during the task's execution.
    """
    key, value = args
    start_time = time.time()
    try:
        before = pickle_environment(NRPyEnv())
        # logging.debug(f"Starting task with key: {key}")
        result = changed_environment(parallel_function_call(value), before)
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
    """Perform parallel code generation by calling registered functions concurrently."""
    if not par.parval_from_str("enable_parallel_codegen"):
        return

    par.set_parval_from_str("parallel_codegen_stage", "codegen")

    # By default, MacOS adopts the "spawn" method, which breaks the global environment.
    #   Linux, on the other hand, uses the "fork" method, which preserves the environment.
    #   Luckily MacOS does support "fork", though Apple advises against its use because fork
    #   without an exec step can cause weird breakage on macOS with certain libraries.
    #   Empirically, for NRPy, this seems to work fine.
    # Note that prior to using multiprocessing, we used the `multiprocess` library,
    #   which was very janky -- causing all sorts of race conditions.
    set_start_method("fork", force=True)

    with Pool() as pool:  # Set Pool(processes=1) to disable parallel codegen.
        # pool.map returns the results in registration order, whichever worker finishes first.
        worker_results = pool.map(wrapper_func, list(ParallelCodeGen_dict.items()))

    unpack_NRPy_environment_dict(dict(worker_results))


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

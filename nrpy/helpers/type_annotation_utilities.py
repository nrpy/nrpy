"""
Provide utility functions for use by multiple classes.

Author: Steven R. Brandt
        sbrandt **at** cct **dot** lsu **dot** edu
"""

from typing import Any, Tuple

try:
    # Ideally, use get_args from typing...
    from typing import get_args
except ImportError as ae:
    # But if you can't get get_args,
    # create our own. Works for Python 3.6
    def get_args(tp: Any) -> Tuple[Any, ...]:
        """
        Provide the functionality of get_args for earlier versions of Python.

        :param tp: The type to extract arguments from.
        :return: A tuple of arguments extracted from the type.

        >>> get_args(Literal["x", "y", "z"])
        ('x', 'y', 'z')
        >>> get_args(int)
        ()
        """
        if is_type_literal(tp):
            ret = tuple(tp.__values__)
        else:
            ret = tuple()
        return ret


import inspect
from json import dumps
from typing_extensions import (
    Literal,
)  # , get_args, get_origin <- Not in typing_extensions

_literal = Literal["x"]
_tuple = Tuple[str, int]


def is_type_literal(arg: Any) -> bool:
    """
    Determine whether a type matches Literal.

    :param arg: The argument to check.
    :return: True if arg is a Literal type, False otherwise.

    >>> is_type_literal(Literal["a", "b"])
    True
    >>> is_type_literal(int)
    False
    >>> is_type_literal(Tuple[str,int])
    False
    """
    # Checking that the type is the same works
    # on most platforms. In a few rare cases
    # we get 'typing._GenericAlias' for both.
    if type(_literal) is not type(_tuple):
        return type(_literal) is type(arg)

    # This, also, doesn't work everywhere.
    origin = getattr(_literal, "__origin__", None)

    assert origin is not None, "Cannot whether an object is a Literal"

    arg_origin = getattr(arg, "__origin__", None)
    if arg_origin is None:
        return False
    return origin is arg_origin


def validate_literal_arguments() -> None:
    """
    Check that the Literal type annotations of the calling function match the arguments the function was actually called with.

    :raises ValueError: If the value of any parameter doesn't match its allowed Literal values.

    >>> def foo(a:Literal["fish","bird"]):
    ...    validate_literal_arguments()
    ...    print("ok")
    >>> foo("fish")
    ok
    >>> try: foo("dog")
    ... except: "not ok"
    'not ok'
    """
    current_frame = inspect.currentframe()
    assert current_frame is not None
    calling_frame = current_frame.f_back
    assert calling_frame is not None
    calling_func_name = calling_frame.f_code.co_name
    # Hunt for the calling function in the current frame
    if "self" in calling_frame.f_locals:
        # Is it a member function?
        calling_self_ref = calling_frame.f_locals["self"]
        calling_func = getattr(calling_self_ref, calling_func_name)
        # full_name = calling_self_ref.__class__.__name__ + "." + calling_func_name
    else:
        # Is it a regular function? If so, is it locally defined?
        assert calling_frame.f_back is not None
        assert calling_frame.f_back.f_locals is not None
        calling_func = calling_frame.f_back.f_locals.get(calling_func_name, None)
        # If it is a regular function and not locally defined, it must be globally defined.
        if calling_func is None:
            calling_func = calling_frame.f_globals[calling_func_name]
        # full_name = calling_func_name
    signature = inspect.signature(calling_func)
    checked_pars = []
    for parameter_name in signature.parameters:
        if parameter_name == "self":
            continue
        parameter_annotation = signature.parameters[parameter_name].annotation
        if parameter_annotation is None:
            continue
        if not is_type_literal(parameter_annotation):
            continue
        parameter_value = calling_frame.f_locals[parameter_name]
        checked_pars += [(parameter_name, parameter_value)]
        allowed_values = get_args(parameter_annotation)
        if parameter_value not in allowed_values:
            raise ValueError(
                f"In function '{calling_func_name}': parameter '{parameter_name}' has value: '{parameter_value}', which is not in the allowed_values set: {allowed_values}"
            )


def generate_class_representation() -> str:
    """
    Generate a useful value for returning in a __repr__() function.

    :return: A string representation of the calling object.

    # Example usage inside a class:
    >>> class MyClass:
    ...     def __init__(self, name: str, age: int):
    ...         self.name = name
    ...         self.age = age
    ...     def __repr__(self):
    ...         return generate_class_representation()
    >>> repr(MyClass("John", 30))
    'MyClass(age=30, name="John")'
    """
    current_frame = inspect.currentframe()
    assert current_frame is not None
    calling_frame = current_frame.f_back
    assert calling_frame is not None
    # calling_func_name = calling_frame.f_code.co_name
    calling_self_ref = calling_frame.f_locals["self"]
    calling_class_name: str = calling_self_ref.__class__.__name__
    args = []
    for d in dir(calling_self_ref):
        if d.startswith("_"):
            continue
        v = getattr(calling_self_ref, d)
        if type(v) in [str, int]:
            args += [d + "=" + dumps(v)]
    sorted(args)
    return calling_class_name + "(" + ", ".join(args) + ")"


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

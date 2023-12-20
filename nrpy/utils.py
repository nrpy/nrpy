from typing_extensions import Literal #, get_args, get_origin <- Not in typing_extensions
import inspect
from json import dumps

_literal = Literal["x"]

def is_literal(arg):
    return type(_literal) == type(arg)

try:
    # Ideally, use get_args from typing...
    from typing import get_args
except ImportError as ae:
    # But if you can't get get_args,
    # create our own. Works for Python 3.6
    def get_args(lit):
        if is_literal(lit):
            return lit.__values__
        else:
            return []

def check_literals()->None:
    """
    Check that the Literal type annotations of the calling function match
    the arguments the function was actually called with.
    """
    calling_frame = inspect.currentframe().f_back
    calling_func_name = calling_frame.f_code.co_name
    # Hunt for the calling function in the current frame
    if 'self' in calling_frame.f_locals:
        # Is it a member function?
        calling_self_ref = calling_frame.f_locals['self']
        calling_func = getattr(calling_self_ref, calling_func_name)
        full_name = calling_self_ref.__class__.__name__ + "." + calling_func_name
    else:
        # Is it a regular function? If so, is it locally defined?
        calling_func = calling_frame.f_back.f_locals.get(calling_func_name, None)
        # If it is a regular function and not locally defined, it must be globally defined.
        if calling_func is None:
            calling_func = calling_frame.f_globals[calling_func_name]
        full_name = calling_func_name
    signature = inspect.signature(calling_func)
    checked_pars = []
    for parameter_name in signature.parameters:
        if parameter_name == "self":
            continue
        parameter_annotation = signature.parameters[parameter_name].annotation
        if parameter_annotation is None:
            continue
        if not is_literal(parameter_annotation):
            continue
        parameter_value = calling_frame.f_locals[parameter_name]
        checked_pars += [(parameter_name, parameter_value)]
        allowed_values = get_args(parameter_annotation)
        if parameter_value not in allowed_values:
            raise ValueError(f"In function '{calling_func_name}': parameter '{parameter_name}' has value: '{parameter_value}', which is not in the allowed_values set: {allowed_values}")

def get_repr()->str:
    """
    Generate a useful value for returning in a __repr__() function.
    Usage:
        class Foo:
           ...
        def __repr__(self):
           return get_repr()
    """
    calling_frame = inspect.currentframe().f_back
    calling_func_name = calling_frame.f_code.co_name
    calling_self_ref = calling_frame.f_locals['self']
    calling_class_name = calling_self_ref.__class__.__name__
    args = []
    for d in dir(calling_self_ref):
        if d.startswith("_"):
            continue
        v = getattr(calling_self_ref,d)
        if type(v) in [str,int]:
            args += [d+"="+dumps(v)]
    sorted(args)
    return calling_class_name + "(" + ", ".join(args) + ")"

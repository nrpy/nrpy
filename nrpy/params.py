"""
This core NRPy+ module is used for
 initializing, storing, and recalling
 parameters.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import warnings
from typing import List, Union, Any, Dict
import sympy as sp


very_verbose = False


class NRPyParameter:
    """
    Doctest:
    >>> param = NRPyParameter(py_type=str, module='module1', name='param0', value='defaultval1')
    >>> param.name in glb_params_dict
    True
    >>> print(glb_params_dict[param.name].__dict__)
    {'py_type': <class 'str'>, 'module': 'module1', 'name': 'param0', 'value': 'defaultval1'}
    """

    def __init__(self, py_type: Any, module: str, name: str, value: Any):
        self.py_type = py_type
        self.module = module
        self.name = name
        self.value = value

        if self.py_type not in (bool, int, float, str):
            raise ValueError("NRPy parameters must have the type int, float, or str")

        if self.name not in glb_params_dict:
            glb_params_dict[self.name] = self
        elif very_verbose:
            warnings.warn(
                f"NRPyParameter minor warning: Did nothing; already initialized parameter {self.module}::{self.name}"
            )


class CodeParameter:
    """
    Class for handling code parameters.

    Doctest:
    >>> param = CodeParameter(c_type_alias='int', module='module2', name='param1', defaultvalue='defaultval2')
    >>> param.name in glb_code_params_dict
    True
    """

    def __init__(
        self,
        c_type_alias: str,
        module: str,
        name: str,
        defaultvalue: Any = "unset",
        assumption: str = "Real",
        commondata: bool = False,
        add_to_parfile: bool = True,
        add_to_set_CodeParameters_h: bool = True,
        add_to_glb_code_params_dict: bool = True,
    ) -> None:
        """
        Initializes the CodeParameter object.

        :param c_type_alias: C type alias for the parameter.
        :param module: The module where the parameter is defined.
        :param name: The name of the parameter.
        :param defaultvalue: The default value for the parameter (default: "unset").
        :param assumption: Assumption related to the symbol; can be "Real" or "RealPositive" (default: "Real").
        :param commondata: Parameter is common to all grids (True), or each grid has its own value for this parameter (default: False).
        :param add_to_parfile: Allow parameter to be set within a parameter file (default: True).
        :param add_to_set_CodeParameters_h: Add parameter to set_CodeParameters*.h (default: True). Only applies for BHaH.
        :param add_to_glb_code_params_dict: Whether to add the parameter to the global code parameters dictionary (default: True).
        """

        self.c_type_alias = c_type_alias
        self.module = module
        self.name = name
        self.defaultvalue = defaultvalue
        self.assumption = assumption
        self.add_to_parfile = add_to_parfile
        self.add_to_set_CodeParameters_h = add_to_set_CodeParameters_h
        self.commondata = commondata
        self.add_to_glb_code_params_dict = add_to_glb_code_params_dict

        if c_type_alias == "#define":
            self.add_to_parfile = False
            self.add_to_set_CodeParameters_h = False

        if self.add_to_parfile and defaultvalue == "unset":
            raise ValueError(
                f"Parameter {self.name}: Must set a default value for all parameters with add_to_parfile=True"
            )

        if assumption == "Real":
            self.symbol = sp.Symbol(self.name, real=True)
        elif assumption == "RealPositive":
            self.symbol = sp.Symbol(self.name, real=True, positive=True)
        else:
            raise ValueError(
                f"Unsupported assumption: {self.assumption}. Only 'Real' and 'RealPositive' are supported."
            )

        if add_to_glb_code_params_dict:
            if self.name not in glb_code_params_dict:
                glb_code_params_dict[self.name] = self
            elif very_verbose:
                warnings.warn(
                    f"initialize_code_param() minor warning: Did nothing; already initialized parameter {self.module}::{self.name}"
                )


# Where we store NRPy+ parameters and default values of parameters.
glb_params_dict: Dict[str, NRPyParameter] = {}
# Where we store C runtime parameters and default values of parameters.
glb_code_params_dict: Dict[str, CodeParameter] = {}


def register_param(py_type: Any, module: str, name: str, value: Any) -> None:
    """
    Initialize a parameter, adding it to the global parameters list if not already present.

    Args:
    param: Parameter object to initialize.
    is_code_param: Boolean flag indicating if the parameter is a C parameter.

    Returns:
    None

    Doctest:
    >>> parname = "param1"
    >>> register_param(py_type=str, module='module1', name=parname, value='defaultval1')
    >>> parname in glb_params_dict
    True
    >>> print(glb_params_dict[parname].__dict__)
    {'py_type': <class 'str'>, 'module': 'module1', 'name': 'param1', 'value': 'defaultval1'}
    """
    NRPyParameter(py_type=py_type, module=module, name=name, value=value)


def parval_from_str(parameter_name: str) -> Any:
    """
    Applies to NRPyParameter objects only!
    Given a string, find the value of the parameter with that name in the global parameters list.
    Ignores any prefix including and before '::' in the parameter name.

    Args:
    parameter_name (str): Name of the parameter to find. Ignores any prefix including and before '::'.

    Returns:
    Any: Value of the parameter in the global parameters list. Raises a ValueError if the parameter is not found.

    Doctest:
    >>> try: parval_from_str('non_existent_param')
    ... except ValueError: print("test passes!")
    test passes!
    >>> register_param(py_type=int, module='blah', name='pvalfromstrtest', value=4)
    >>> print(parval_from_str('pvalfromstrtest'))
    4
    >>> print(parval_from_str('blah::pvalfromstrtest'))
    4
    """
    # Splitting the parameter_name by '::' and using the last part if '::' is present
    actual_param_name = parameter_name.split("::")[-1]

    if actual_param_name in glb_params_dict:
        return glb_params_dict[actual_param_name].value

    raise ValueError(
        f"Error: parameter {actual_param_name} not found in {str(glb_params_dict)}"
    )


def set_parval_from_str(parameter_name: str, new_value: Any) -> None:
    """
    Applies to NRPyParameter objects only!
    Given a string and a value, update the parameter with that name in the global parameters list to the new value.

    Args:
    parameter_name (str): Name of the parameter to update.
    new_value (Any): New value to be assigned to the parameter.

    Raises:
    ValueError: If the parameter with the specified name is not found in the global parameters list.

    Doctest:
    >>> try: set_parval_from_str('non_existent_param', -100)
    ... except ValueError: print("Test passes!")
    Test passes!
    """
    if parameter_name in glb_params_dict:
        glb_params_dict[parameter_name].value = new_value
    else:
        raise ValueError(
            f"Error: parameter {parameter_name} not found in {str(glb_params_dict)}"
        )


def register_CodeParameters(
    c_type_alias: str,
    module: str,
    names: List[str],
    defaultvalues: Union[str, int, float, List[Union[str, int, float]]] = "unset",
    assumption: str = "Real",
    commondata: bool = False,
    add_to_parfile: bool = True,
    add_to_set_CodeParameters_h: bool = True,
    add_to_glb_code_params_dict: bool = True,
) -> List[sp.Symbol]:
    """
    This function initializes the parameters and returns the symbolic representation of them.

    :param c_type_alias: C type alias for the parameter.
    :param module: The module where the parameter is defined.
    :param name: The name of the parameter.
    :param defaultvalues: A list of the default values for the parameters. If it's a single value, it will be duplicated for all parameters.
    :param assumption: Assumption related to the symbol; can be "Real" or "RealPositive" (default: "Real").
    :param commondata: Parameter is common to all grids (True), or each grid has its own value for this parameter (default: False).
    :param add_to_parfile: Allow parameter to be set within a parameter file (default: True).
    :param add_to_set_CodeParameters_h: Add parameter to set_CodeParameters*.h (default: True). Only applies for BHaH.
    :param add_to_glb_code_params_dict: Whether to add the parameter to the global code parameters dictionary (default: True).

    :return: A list of the symbolic parameters. If there's only one parameter, it directly returns the symbol.

    :raises ValueError: If the assumption is not supported, or if the lengths of names and defaultvalues do not match.

    Doctest:
    >>> glb_code_params_dict.clear()
    >>> a0,a1,b,c = register_CodeParameters(c_type_alias="REAL", module=__name__, names=["a0", "a1", "b", "c"], defaultvalues=1)
    >>> print(a0,a1,b,c)
    a0 a1 b c
    >>> outstr = ""
    >>> outstr += " ".join(param.name for _, param in glb_code_params_dict.items())
    >>> print(outstr)
    a0 a1 b c
    >>> outstr = ""
    >>> outstr += " ".join(str(param.defaultvalue) for _, param in glb_code_params_dict.items())
    >>> print(outstr)
    1 1 1 1
    """

    if not isinstance(names, list) or len(names) == 1:
        raise ValueError(
            "register_CodeParameters() expects a list of MULTIPLE code parameter names."
        )

    if not isinstance(defaultvalues, list):
        default_val_list = [defaultvalues] * len(names)
    elif len(defaultvalues) != len(names):
        raise ValueError(
            f"The lengths of names ({len(names)}) and defaultvalues ({len(defaultvalues)}) do not match."
        )
    else:
        default_val_list = defaultvalues

    symbols = []
    for name, default_val in zip(names, default_val_list):
        CP = CodeParameter(
            c_type_alias=c_type_alias,
            module=module,
            name=name,
            defaultvalue=default_val,
            assumption=assumption,
            commondata=commondata,
            add_to_parfile=add_to_parfile,
            add_to_set_CodeParameters_h=add_to_set_CodeParameters_h,
            add_to_glb_code_params_dict=add_to_glb_code_params_dict,
        )
        symbols.append(CP.symbol)

    return symbols


def register_CodeParameter(
    c_type_alias: str,
    module: str,
    name: str,
    defaultvalue: Union[str, int, float] = "unset",
    assumption: str = "Real",
    commondata: bool = False,
    add_to_parfile: bool = True,
    add_to_set_CodeParameters_h: bool = True,
    add_to_glb_code_params_dict: bool = True,
) -> sp.Symbol:
    """
    This function initializes the parameters and returns the symbolic representation of them.

    :param c_type_alias: C type alias for the parameter.
    :param module: The module where the parameter is defined.
    :param name: The name of the parameter.
    :param defaultvalue: The default value for the parameter (default: "unset").
    :param assumption: Assumption related to the symbol; can be "Real" or "RealPositive" (default: "Real").
    :param commondata: Parameter is common to all grids (True), or each grid has its own value for this parameter (default: False).
    :param add_to_parfile: Allow parameter to be set within a parameter file (default: True).
    :param add_to_set_CodeParameters_h: Add parameter to set_CodeParameters*.h (default: True). Only applies for BHaH.
    :param add_to_glb_code_params_dict: Whether to add the parameter to the global code parameters dictionary (default: True).

    :return: A symbolic parameter.

    :raises TypeError: If the name is not a simple string.

    Doctest:
    >>> glb_code_params_dict.clear()
    >>> a0 = register_CodeParameter(c_type_alias="REAL", module=__name__, name="a0", add_to_parfile=False, assumption="Real")
    >>> print(a0)
    a0
    >>> outstr = ""
    >>> outstr += " ".join(param.name for _, param in glb_code_params_dict.items())
    >>> print(outstr)
    a0
    >>> outstr = ""
    >>> outstr += " ".join(str(param.defaultvalue) for _, param in glb_code_params_dict.items())
    >>> print(outstr)
    unset
    >>> a1 = register_CodeParameter(c_type_alias="REAL", module=__name__, name="a1", assumption="Real")
    Traceback (most recent call last):
      ...
    ValueError: Parameter a1: Must set a default value for all parameters with add_to_parfile=True
    """
    if not isinstance(name, str):
        raise TypeError(
            "register_Codeparameter() expects the input name to be a simple string. Use register_Codeparameters() to register multiple CodeParameters."
        )

    CP = CodeParameter(
        c_type_alias=c_type_alias,
        module=module,
        name=name,
        defaultvalue=defaultvalue,
        assumption=assumption,
        commondata=commondata,
        add_to_parfile=add_to_parfile,
        add_to_set_CodeParameters_h=add_to_set_CodeParameters_h,
        add_to_glb_code_params_dict=add_to_glb_code_params_dict,
    )
    return CP.symbol


def adjust_CodeParam_default(CodeParameter_name: str, new_default: Any) -> None:
    """
    Adjusts the default value of a given code parameter.

    :param CodeParameter_name: The name of the code parameter to be adjusted.
    :param new_default: The new default value for the code parameter.

    :raises ValueError: If the given CodeParameter_name does not exist in glb_code_params_dict.

    >>> _ = register_CodeParameter("REAL", __name__, "dummy", 1.0)
    >>> glb_code_params_dict["dummy"].defaultvalue
    1.0
    >>> adjust_CodeParam_default("dummy", 2.0)
    >>> glb_code_params_dict["dummy"].defaultvalue
    2.0
    >>> adjust_CodeParam_default('param3', 'blah')
    Traceback (most recent call last):
    ...
    ValueError: Cannot adjust "param3" default value, as it does not appear in glb_code_params_dict
    """
    if CodeParameter_name in glb_code_params_dict:
        glb_code_params_dict[CodeParameter_name].defaultvalue = new_default
    else:
        raise ValueError(
            f'Cannot adjust "{CodeParameter_name}" default value, as it does not appear in glb_code_params_dict'
        )


# Valid default parameter.
register_param(str, __name__, "Infrastructure", "BHaH")

if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

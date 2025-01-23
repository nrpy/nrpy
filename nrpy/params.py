"""
Initialize, store, and recall parameters for NRPy+.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import warnings
from typing import Any, Dict, List, Union

import sympy as sp

very_verbose = False


class NRPyParameter:
    """
    Represent an NRPyParameter object.

    DocTests:
    >>> param = NRPyParameter(py_type=str, module='module1', name='param0', value='defaultval1', description='Test parameter.')
    >>> param.name in glb_params_dict
    True
    >>> print(glb_params_dict[param.name].__dict__)
    {'py_type': <class 'str'>, 'module': 'module1', 'name': 'param0', 'value': 'defaultval1', 'description': 'Test parameter.'}
    >>> param_no_desc = NRPyParameter(py_type=int, module='module1', name='param_no_desc', value=42)
    >>> print(glb_params_dict[param_no_desc.name].__dict__)
    {'py_type': <class 'int'>, 'module': 'module1', 'name': 'param_no_desc', 'value': 42, 'description': ''}
    """

    def __init__(
        self, py_type: Any, module: str, name: str, value: Any, description: str = ""
    ):
        """
        Initialize an NRPyParameter object with given properties.

        :param py_type: The type of the parameter, e.g., str, int, float.
        :param module: Name of the module the parameter belongs to.
        :param name: Name of the parameter.
        :param value: The default value of the parameter.
        :param description: Description of the parameter.
        :raises ValueError: If the py_type is not one of the valid types.
        """
        self.py_type = py_type
        self.module = module
        self.name = name
        self.value = value
        self.description = description

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
    Represent a code parameter.

    DocTests:
    >>> param = CodeParameter(cparam_type='int', module='module2', name='param1', defaultvalue='defaultval2', description='A code parameter.')
    >>> param.name in glb_code_params_dict
    True
    >>> param_no_desc = CodeParameter(cparam_type='float', module='module2', name='param_no_desc', defaultvalue=3.14)
    >>> print(glb_code_params_dict[param_no_desc.name].__dict__)
    {'cparam_type': 'float', 'module': 'module2', 'name': 'param_no_desc', 'defaultvalue': 3.14, 'assumption': 'Real', 'add_to_parfile': True, 'add_to_set_CodeParameters_h': True, 'commondata': False, 'add_to_glb_code_params_dict': True, 'description': '', 'symbol': param_no_desc}
    >>> try:
    ...     param_array = CodeParameter(cparam_type='REAL[5]', module='module3', name='param_array', commondata=False, add_to_set_CodeParameters_h=True, add_to_parfile=False)
    ... except ValueError as e:
    ...     print(e)
    Parameter 'param_array' of type 'REAL[5]': For REAL or int array parameters, commondata must be True, and add_to_set_CodeParameters_h must be False.
    >>> # Correct usage with commondata=True and add_to_set_CodeParameters_h=False
    >>> param_array_correct = CodeParameter(cparam_type='REAL[5]', module='module3', name='param_array_correct', defaultvalue=-10.1, commondata=True, add_to_set_CodeParameters_h=False, add_to_parfile=True)
    >>> param_array_correct.name in glb_code_params_dict
    True
    """

    def __init__(
        self,
        cparam_type: str,
        module: str,
        name: str,
        defaultvalue: Any = "unset",
        assumption: str = "Real",
        commondata: bool = False,
        add_to_parfile: bool = True,
        add_to_set_CodeParameters_h: bool = True,
        add_to_glb_code_params_dict: bool = True,
        description: str = "",
    ) -> None:
        """
        Initialize a CodeParameter with various properties and assumptions.

        :param cparam_type: The data type in C/C++ for the parameter.
        :param module: The module where this parameter is defined.
        :param name: Name of the parameter.
        :param defaultvalue: Default value of the parameter, "unset" by default.
        :param assumption: Type of symbolic assumption, "Real" or "RealPositive".
        :param commondata: If True, parameter is common across all grids. Defaults to False.
        :param add_to_parfile: If True, include this parameter in parameter files. Default is True.
        :param add_to_set_CodeParameters_h: If True, add to set_CodeParameters*.h, applicable for BHaH. Default is True.
        :param add_to_glb_code_params_dict: If True, add to global code parameters dictionary. Default is True.
        :param description: Description of the parameter.
        :raises ValueError: If `defaultvalue` is "unset" for a parameter intended for a parfile, if an unsupported
                            assumption is specified, or if defaultvalues array size is inconsistent with array size.
        """
        self.cparam_type = cparam_type
        self.module = module
        self.name = name
        self.defaultvalue = defaultvalue
        self.assumption = assumption
        self.add_to_parfile = add_to_parfile
        self.add_to_set_CodeParameters_h = add_to_set_CodeParameters_h
        self.commondata = commondata
        self.add_to_glb_code_params_dict = add_to_glb_code_params_dict
        self.description = description

        # Check if cparam_type is a REAL or int array
        if (
            (self.cparam_type.startswith("REAL") or self.cparam_type.startswith("int"))
            and "[" in self.cparam_type
            and "]" in self.cparam_type
        ):
            if not (self.commondata and not self.add_to_set_CodeParameters_h):
                raise ValueError(
                    f"Parameter '{self.name}' of type '{self.cparam_type}': "
                    "For REAL or int array parameters, commondata must be True, "
                    "and add_to_set_CodeParameters_h must be False."
                )
            arr_size = int(self.cparam_type.split("[")[1].split("]")[0])
            if isinstance(self.defaultvalue, list):
                if arr_size != len(self.defaultvalue):
                    raise ValueError(
                        f"{module}::{name}: Length of default values list {len(self.defaultvalue)} != size = {arr_size}."
                    )
            else:
                self.defaultvalue = [self.defaultvalue] * arr_size

        if cparam_type == "#define":
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
# Where we store dictionaries of dictionaries
glb_extras_dict: Dict[str, Dict[str, Any]] = {}


def register_param(
    py_type: Any, module: str, name: str, value: Any, description: str = ""
) -> None:
    """
    Initialize a parameter and add it to the global parameters list if not already present.

    :param py_type: The type of the parameter, e.g., str, int, float.
    :param module: Name of the module the parameter belongs to.
    :param name: Name of the parameter.
    :param value: The default value of the parameter.
    :param description: Description of the parameter.

    DocTests:
    >>> parname = "param1"
    >>> register_param(py_type=str, module='module1', name=parname, value='defaultval1', description='Test parameter.')
    >>> parname in glb_params_dict
    True
    >>> print(glb_params_dict[parname].__dict__)
    {'py_type': <class 'str'>, 'module': 'module1', 'name': 'param1', 'value': 'defaultval1', 'description': 'Test parameter.'}
    >>> register_param(py_type=int, module='module1', name='param_no_desc', value=42)
    >>> print(glb_params_dict['param_no_desc'].__dict__)
    {'py_type': <class 'int'>, 'module': 'module1', 'name': 'param_no_desc', 'value': 42, 'description': ''}
    """
    NRPyParameter(
        py_type=py_type, module=module, name=name, value=value, description=description
    )


def parval_from_str(parameter_name: str) -> Any:
    """
    Retrieve the parameter value from its name as a string.

    :param parameter_name: The name of the parameter. If the parameter name includes a module prefix ('module::name'),
                           only the part after '::' is considered.
    :return: The value of the parameter.
    :raises ValueError: If the specified parameter name does not exist in the global parameter dictionary.

    DocTests:
    >>> try:
    ...     parval_from_str('non_existent_param')
    ... except ValueError:
    ...     print("test passes!")
    test passes!
    >>> register_param(py_type=int, module='blah', name='pvalfromstrtest', value=4)
    >>> print(parval_from_str('pvalfromstrtest'))
    4
    >>> print(parval_from_str('blah::pvalfromstrtest'))
    4
    """
    actual_param_name = parameter_name.split("::")[-1]

    if actual_param_name in glb_params_dict:
        return glb_params_dict[actual_param_name].value

    raise ValueError(
        f"Error: parameter {actual_param_name} not found in {str(glb_params_dict)}"
    )


def set_parval_from_str(parameter_name: str, new_value: Any) -> None:
    """
    Apply to NRPyParameter objects only.

    Given a string and a value, this function updates the parameter with that name in the global parameters list
    to the new value.

    :param parameter_name: Name of the parameter to update.
    :param new_value: New value to be assigned to the parameter.
    :raises ValueError: If the parameter name is not found in the global parameters list.

    DocTests:
    >>> try:
    ...     set_parval_from_str('non_existent_param', -100)
    ... except ValueError:
    ...     print("Test passes!")
    Test passes!
    """
    if parameter_name in glb_params_dict:
        glb_params_dict[parameter_name].value = new_value
    else:
        raise ValueError(
            f"Error: parameter {parameter_name} not found in {str(glb_params_dict)}"
        )


def register_CodeParameters(
    cparam_type: str,
    module: str,
    names: List[str],
    defaultvalues: Union[str, int, float, List[Union[str, int, float]]] = "unset",
    assumption: str = "Real",
    commondata: bool = False,
    add_to_parfile: bool = True,
    add_to_set_CodeParameters_h: bool = True,
    add_to_glb_code_params_dict: bool = True,
    descriptions: Union[str, List[str]] = "",
) -> List[sp.Symbol]:
    """
    Initialize CodeParameters and return their symbolic representation.

    :param cparam_type: C/C++ type for the parameter.
    :param module: The module where the parameter is defined.
    :param names: The names of the parameters.
    :param defaultvalues: A list of the default values for the parameters. If it's a single value, it will be duplicated for all parameters.
    :param assumption: Assumption related to the symbol; can be "Real" or "RealPositive" (default: "Real").
    :param commondata: Parameter is common to all grids (True), or each grid has its own value for this parameter (default: False).
    :param add_to_parfile: Allow parameter to be set within a parameter file (default: True).
    :param add_to_set_CodeParameters_h: Add parameter to set_CodeParameters*.h (default: True). Only applies for BHaH.
    :param add_to_glb_code_params_dict: Whether to add the parameter to the global code parameters dictionary (default: True).
    :param descriptions: Either an empty string or a list of descriptions corresponding to each parameter.
    :return: A list of the symbolic parameters.

    :raises ValueError: If the assumption is not supported, if the lengths of names and defaultvalues do not match, or if descriptions are improperly provided.
    :raises TypeError: If descriptions are not a string or a list.

    DocTests:
    >>> glb_code_params_dict.clear()
    >>> # Register parameters without descriptions
    >>> a0, a1, b, c = register_CodeParameters(
    ...     cparam_type="REAL",
    ...     module=__name__,
    ...     names=["a0", "a1", "b", "c"],
    ...     defaultvalues=1
    ... )
    >>> print(a0, a1, b, c)
    a0 a1 b c
    >>> outstr = " ".join(param.name for param in glb_code_params_dict.values())
    >>> print(outstr)
    a0 a1 b c
    >>> outstr = " ".join(str(param.defaultvalue) for param in glb_code_params_dict.values())
    >>> print(outstr)
    1 1 1 1

    >>> glb_code_params_dict.clear()
    >>> # Register parameters with descriptions
    >>> descs = ["First parameter", "Second parameter", "Third parameter", "Fourth parameter"]
    >>> a0, a1, b, c = register_CodeParameters(
    ...     cparam_type="int",
    ...     module=__name__,
    ...     names=["a0", "a1", "b", "c"],
    ...     defaultvalues=[10, 20, 30, 40],
    ...     descriptions=descs
    ... )
    >>> print(a0, a1, b, c)
    a0 a1 b c
    >>> outstr = " ".join(param.name for param in glb_code_params_dict.values())
    >>> print(outstr)
    a0 a1 b c
    >>> outstr = " ".join(str(param.defaultvalue) for param in glb_code_params_dict.values())
    >>> print(outstr)
    10 20 30 40
    >>> outstr = " ".join(param.description for param in glb_code_params_dict.values())
    >>> print(outstr)
    First parameter Second parameter Third parameter Fourth parameter

    >>> glb_code_params_dict.clear()
    >>> # Register parameters with mismatched descriptions
    >>> try:
    ...     register_CodeParameters(
    ...         cparam_type="bool",
    ...         module=__name__,
    ...         names=["flag1", "flag2"],
    ...         defaultvalues=[True, False],
    ...         descriptions=["Enable feature"]
    ... )
    ... except ValueError as e:
    ...     print(e)
    The length of descriptions (1) does not match the number of names (2).
    >>> glb_code_params_dict.clear()
    >>> # Register parameters with non-empty single string description (should raise error)
    >>> try:
    ...     register_CodeParameters(
    ...         cparam_type="char",
    ...         module=__name__,
    ...         names=["char1", "char2"],
    ...         defaultvalues=["a", "b"],
    ...         descriptions="Single description"
    ... )
    ... except ValueError as e:
    ...     print(e)
    Description must be an empty string or a list of descriptions matching the number of parameters.
    """
    if not isinstance(names, list) or len(names) < 1:
        raise ValueError(
            "register_CodeParameters() expects a list of one or more code parameter names."
        )

    # Handle defaultvalues
    if not isinstance(defaultvalues, list):
        default_val_list = [defaultvalues] * len(names)
    elif len(defaultvalues) != len(names):
        raise ValueError(
            f"The lengths of names ({len(names)}) and defaultvalues ({len(defaultvalues)}) do not match."
        )
    else:
        default_val_list = defaultvalues

    # Handle descriptions
    if isinstance(descriptions, list):
        if len(descriptions) != len(names):
            raise ValueError(
                f"The length of descriptions ({len(descriptions)}) does not match the number of names ({len(names)})."
            )
        description_list = descriptions
    elif isinstance(descriptions, str):
        if descriptions != "":
            raise ValueError(
                "Description must be an empty string or a list of descriptions matching the number of parameters."
            )
        description_list = [""] * len(names)
    else:
        raise TypeError(
            "Description must be either an empty string or a list of description strings."
        )

    symbols = []
    for name, default_val, desc in zip(names, default_val_list, description_list):
        CP = CodeParameter(
            cparam_type=cparam_type,
            module=module,
            name=name,
            defaultvalue=default_val,
            assumption=assumption,
            commondata=commondata,
            add_to_parfile=add_to_parfile,
            add_to_set_CodeParameters_h=add_to_set_CodeParameters_h,
            add_to_glb_code_params_dict=add_to_glb_code_params_dict,
            description=desc,
        )
        symbols.append(CP.symbol)

        if add_to_glb_code_params_dict:
            glb_code_params_dict[name] = CP

    return symbols


def register_CodeParameter(
    cparam_type: str,
    module: str,
    name: str,
    defaultvalue: Union[str, int, float, List[Any]] = "unset",
    assumption: str = "Real",
    commondata: bool = False,
    add_to_parfile: bool = True,
    add_to_set_CodeParameters_h: bool = True,
    add_to_glb_code_params_dict: bool = True,
    description: str = "",
) -> sp.Symbol:
    """
    Initialize a CodeParameter and return its symbolic representation.

    :param cparam_type: C/C++ type for the parameter.
    :param module: The module where the parameter is defined.
    :param name: The name of the parameter.
    :param defaultvalue: The default value for the parameter (default: "unset").
    :param assumption: Assumption related to the symbol; can be "Real" or "RealPositive" (default: "Real").
    :param commondata: Parameter is common to all grids (True), or each grid has its own value for this parameter (default: False).
    :param add_to_parfile: Allow parameter to be set within a parameter file (default: True).
    :param add_to_set_CodeParameters_h: Add parameter to set_CodeParameters*.h (default: True). Only applies for BHaH.
    :param add_to_glb_code_params_dict: Whether to add the parameter to the global code parameters dictionary (default: True).
    :param description: Description of the parameter.
    :return: A symbolic parameter.

    :raises TypeError: If the name is not a simple string.

    DocTests:
    >>> glb_code_params_dict.clear()
    >>> a0 = register_CodeParameter(cparam_type="REAL", module=__name__, name="a0", add_to_parfile=False, assumption="Real", description='Test parameter.')
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
    >>> a1 = register_CodeParameter(cparam_type="REAL", module=__name__, name="a1", defaultvalue=0.0)
    >>> print(a1)
    a1
    >>> a2 = register_CodeParameter(cparam_type="REAL", module=__name__, name="a2", assumption="Real")
    Traceback (most recent call last):
      ...
    ValueError: Parameter a2: Must set a default value for all parameters with add_to_parfile=True
    >>> a3 = register_CodeParameter(cparam_type="REAL", module=__name__, name="a2", add_to_parfile=False)
    """
    if not isinstance(name, str):
        raise TypeError(
            "register_Codeparameter() expects the input name to be a simple string. Use register_Codeparameters() to register multiple CodeParameters."
        )

    CP = CodeParameter(
        cparam_type=cparam_type,
        module=module,
        name=name,
        defaultvalue=defaultvalue,
        assumption=assumption,
        commondata=commondata,
        add_to_parfile=add_to_parfile,
        add_to_set_CodeParameters_h=add_to_set_CodeParameters_h,
        add_to_glb_code_params_dict=add_to_glb_code_params_dict,
        description=description,
    )
    return CP.symbol


def adjust_CodeParam_default(CodeParameter_name: str, new_default: Any) -> None:
    """
    Adjust the default value of a given code parameter.

    :param CodeParameter_name: The name of the code parameter to be adjusted.
    :param new_default: The new default value for the code parameter.
    :raises ValueError: If the given CodeParameter_name does not exist in glb_code_params_dict.

    DocTests:
    >>> _ = register_CodeParameter("REAL", __name__, "dummy", 1.0)
    >>> glb_code_params_dict["dummy"].defaultvalue
    1.0
    >>> adjust_CodeParam_default("dummy", 2.0)
    >>> glb_code_params_dict["dummy"].defaultvalue
    2.0
    >>> adjust_CodeParam_default('param3', 'blah')
    Traceback (most recent call last):
    ...
    ValueError: Cannot adjust "param3" default value, as it does not appear in glb_code_params_dict: ['dummy', 'param1', 'param_array_correct', 'param_no_desc']
    >>> _ = register_CodeParameter("REAL[2]", __name__, "dummy_arr", [1.0, 2.5], commondata=True, add_to_set_CodeParameters_h=False)
    >>> glb_code_params_dict["dummy_arr"].defaultvalue
    [1.0, 2.5]
    >>> adjust_CodeParam_default("dummy_arr[0]", -1.0)
    >>> glb_code_params_dict["dummy_arr"].defaultvalue
    [-1.0, 2.5]
    >>> _ = register_CodeParameter("REAL[3]", __name__, "dummy_arr2", 0.25, commondata=True, add_to_set_CodeParameters_h=False)
    >>> glb_code_params_dict["dummy_arr2"].defaultvalue
    [0.25, 0.25, 0.25]
    >>> adjust_CodeParam_default("dummy_arr2[2]", -12.5)
    >>> glb_code_params_dict["dummy_arr2"].defaultvalue
    [0.25, 0.25, -12.5]
    """
    found_parameter = False
    if "[" in CodeParameter_name:
        for CPname, CPvalue in glb_code_params_dict.items():
            if CPname.split("[")[0] == CodeParameter_name.split("[")[0]:
                index = int(CodeParameter_name.split("[")[1].split("]")[0])
                CPvalue.defaultvalue[index] = new_default
                found_parameter = True
                break
    elif CodeParameter_name in glb_code_params_dict:
        glb_code_params_dict[CodeParameter_name].defaultvalue = new_default
        found_parameter = True

    if not found_parameter:
        raise ValueError(
            f'Cannot adjust "{CodeParameter_name}" default value, as it does not appear in glb_code_params_dict: {sorted(glb_code_params_dict.keys())}'
        )


# Valid default parameters.
register_param(str, __name__, "Infrastructure", "BHaH")
register_param(str, __name__, "fp_type", "double")

if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

"""
This module provides registration functions for
 numerical grid parameters and gridfunctions across
 various infrastructures.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""
from typing import List, Any, Union, Dict, Tuple, Optional, Sequence
import sympy as sp
import nrpy.indexedexp as ixp
import nrpy.params as par


Cart_origin = par.register_CodeParameters(
    "REAL", __name__, ["Cart_originx", "Cart_originy", "Cart_originz"], 0.0
)
_ = par.register_CodeParameter("int", __name__, "NUMGRIDS", 1, commondata=True)


# Core GridFunction class.
class GridFunction:
    """
    Core class for grid functions.
    """

    def __init__(self, name: str, group: str = "EVOL") -> None:
        self.c_type_alias: str = "double"
        self.name: str = name
        self.group: str = group

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        Catch-all function for gridfunctions that haven't been set up correctly.
        """
        return f"Please define read_gf_from_memory_Ccode_onept() inside grid.py. Inputs: [class instance] {i0_offset} {i1_offset} {i2_offset} {kwargs}"


class BHaHGridFunction(GridFunction):
    """
    Subclass for BlackHoles@Home grid functions.
    """

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
        gf_array_name: str = "in_gfs",
    ) -> None:
        super().__init__(name, group)
        self.c_type_alias = "REAL"  # always use REAL
        self.rank = rank
        self.f_infinity = f_infinity
        self.wavespeed = wavespeed
        self.is_basename = is_basename
        self.gf_array_name = gf_array_name

        if is_basename:
            self.verify_gridfunction_basename_is_valid(self.name)
        self.verify_gridfunction_group_is_valid()

    def verify_gridfunction_basename_is_valid(self, name: str) -> None:
        """
        Validates the gridfunction's base name. Raises ValueError if base name is
        zero-length or ends with an integer. Raises TypeError if the base name is not a string.

        Gridfunction base names must not end with integers. For instance, a rank-1
        gridfunction 'vecU1' has a valid base name 'vecU'. Further, a scalar grid-
        function named 'u2' would have an invalid base name 'u2'. This strict
        requirement facilitates quick identification of a gridfunction by its name,
        which is essential for inferring what a gridfunction represents at a glance.

        Args:
            name (str): The name of the gridfunction.

        Raises:
            ValueError: If the name is of zero length or ends with an integer.
            TypeError: If the name is not of type str.

        Examples:
            >>> try: GF = BHaHGridFunction(name="h1")
            ... except: "Attempt to create a gridfunction ending in an integer failed. Good."
            'Attempt to create a gridfunction ending in an integer failed. Good.'
        """
        # First check for zero-length basenames:
        if len(name) == 0:
            raise ValueError("Tried to register gridfunction without a name!")

        if not isinstance(name, str):
            raise TypeError("Gridfunction names must be strings")

        if len(name) > 0 and name[-1].isdigit():
            raise ValueError(
                f"Tried to register gridfunction with base name: {name}\n"
                f"To ensure that finite difference code generations do not get confused,\n"
                f"gridfunctions with base names ending in an integer are forbidden; pick a new name."
            )

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Validates the gridfunction group. Raises ValueError if the group is not valid.

        The valid groups are 'EVOL', 'AUX' and 'AUXEVOL'.

        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints.

        Raises:
            ValueError: If the group is not valid.

        Example:
            >>> try: GF = BHaHGridFunction(name="blah", group="INVALID")
            ... except: print("Errored out. This is good.")
            Errored out. This is good.
        """
        valid_groups = (
            "EVOL",
            "AUX",
            "AUXEVOL",
        )
        if self.group not in valid_groups:
            msg = (
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
                '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
                '    "AUX": for all other quantities needed at all gridpoints.\n'
            )
            raise ValueError(msg)

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        This function takes a grid function name, three integer offsets, and an optional keyword arguments dictionary,
        then returns a formatted string. If 'gf_array_name' is not provided in kwargs, it raises a ValueError.

        Parameters:
        gf_name (str): The grid function name
        i0_offset (int): Offset for the first index
        i1_offset (int): Offset for the second index
        i2_offset (int): Offset for the third index
        kwargs (dict): Optional keyword arguments, expected to have 'gf_array_name'

        Returns:
        str: Formatted string

        Raises:
        ValueError: If 'gf_array_name' is not provided in kwargs

        Doctests:
        >>> glb_gridfcs_dict.clear()
        >>> par.set_parval_from_str("Infrastructure", "BHaH")
        >>> abc = register_gridfunctions("abc")
        >>> glb_gridfcs_dict["abc"].read_gf_from_memory_Ccode_onept(1, 2, 3)
        'in_gfs[IDX4(ABCGF, i0+1, i1+2, i2+3)]'
        >>> defg = register_gridfunctions("defg")
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0, enable_simd=True, gf_array_name="My_Array")
        'ReadSIMD(&My_Array[IDX4(DEFGGF, i0, i1-1, i2)])'
        """
        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

        gf_array_name = kwargs.get("gf_array_name", self.gf_array_name)

        ret_string = f"{gf_array_name}[IDX4({self.name.upper()}GF, {i0}, {i1}, {i2})]"
        if kwargs.get("enable_simd"):
            return f"ReadSIMD(&{ret_string})"
        return ret_string

    @staticmethod
    def access_gf(
        gf_name: str,
        i0_offset: int = 0,
        i1_offset: int = 0,
        i2_offset: int = 0,
        gf_array_name: str = "in_gfs",
    ) -> str:
        """
        This function takes a grid function name, three integer offsets, and an optional keyword arguments dictionary,
        then returns a formatted string. If 'gf_array_name' is not provided in kwargs, it raises a ValueError.

        Parameters:
        gf_name (str): The grid function name
        i0_offset (int): Offset for the first index
        i1_offset (int): Offset for the second index
        i2_offset (int): Offset for the third index

        Returns:
        str: Formatted string

        Raises:
        ValueError: If 'gf_array_name' is not provided in kwargs

        Doctests:
        >>> BHaHGridFunction.access_gf("aa", 1,2,3)
        'in_gfs[IDX4(AAGF, i0+1, i1+2, i2+3)]'
        >>> BHaHGridFunction.access_gf("defg", 0, -1, 0, "My_Array")
        'My_Array[IDX4(DEFGGF, i0, i1-1, i2)]'
        """
        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

        return f"{gf_array_name}[IDX4({gf_name.upper()}GF, {i0}, {i1}, {i2})]"

    @staticmethod
    def gridfunction_lists() -> Tuple[List[str], List[str], List[str]]:
        """
        Generates sorted lists of gridfunction names for each group type: 'EVOL', 'AUX', and 'AUXEVOL'.

        The function creates a dictionary with group types as keys and corresponding gridfunction names
        as values. It iterates through a global dictionary of gridfunctions, 'glb_gridfcs_dict', and
        appends the gridfunction names to the corresponding group lists in 'groups'. The function then
        sorts each list and returns the sorted lists as a tuple.

        Note: This function assumes the existence of a global gridfunction dictionary named
        'glb_gridfcs_dict' whose keys are gridfunction names and values are gridfunction objects. Each
        gridfunction object should have 'group' and 'name' attributes.

        Returns:
            Tuple[List[str], List[str], List[str]]: A tuple containing lists of gridfunction names
            for each group: 'EVOL', 'AUX', and 'AUXEVOL', respectively.

        Example:
            gridfunction_lists()
            (['evol_gf1', 'evol_gf2'], ['aux_gf1', 'aux_gf2'], ['auxevol_gf1', 'auxevol_gf2'])
        """
        # Initialize dictionary for holding lists of gridfunction names for each group.
        groups: Dict[str, List[str]] = {
            "EVOL": [],
            "AUX": [],
            "AUXEVOL": [],
        }

        # Iterate through the global dictionary of gridfunctions.
        for _key, gf in glb_gridfcs_dict.items():
            if gf.group in groups:
                groups[gf.group].append(gf.name)

        # Sort the lists. Iterating through a copy of the keys to avoid modifying the dictionary while iterating.
        for group in list(groups.keys()):
            groups[group] = sorted(groups[group])

        # Pack the sorted lists into a tuple and return.
        return groups["EVOL"], groups["AUX"], groups["AUXEVOL"]

    @staticmethod
    def gridfunction_defines() -> str:
        """
        Generates a string representation of all the GridFunction definitions.
        This function also includes definitions for the values of each gridfunction at infinity and the characteristic wave speed.
        The output is formatted to be used in C-style languages.

        This function assumes the existence of a global gridfunction dictionary named 'glb_gridfcs_dict'. The keys of this
        dictionary should be gridfunction names, and the values should be instances of classes derived from GridFunction.

        Returns:
            str: The resulting string representation of the GridFunction definitions.
        """

        def define_gfs(name: str, gfs: List[str]) -> str:
            """
            Helper function to generate string representation for a list of GridFunctions.

            Args:
                name (str): The name of the group.
                gfs (List[str]): The list of GridFunctions.

            Returns:
                str: The resulting string representation.
            """
            num_gfs = len(gfs)
            defines = "\n".join(
                f"#define {gf.upper()}GF\t{i}" for i, gf in enumerate(gfs)
            )
            return f"\n// {name.upper()} VARIABLES:\n#define NUM_{name.upper()}_GFS {num_gfs}\n{defines}\n"

        (
            evolved_variables_list,
            auxiliary_variables_list,
            auxevol_variables_list,
        ) = BHaHGridFunction.gridfunction_lists()

        outstr = f"{define_gfs('EVOL', evolved_variables_list)}"
        outstr += f"{define_gfs('AUX', auxiliary_variables_list)}"
        outstr += f"{define_gfs('AUXEVOL', auxevol_variables_list)}"

        if evolved_variables_list:
            outstr += "\n\n// SET gridfunctions_f_infinity[i] = value of gridfunction i in the limit r->infinity:\n"
            f_infinity_str = ""
            f_wavespeed_str = ""
            for var_string in evolved_variables_list:
                gf = glb_gridfcs_dict[var_string]
                # mypy: doesn't consider that this func is only called for BHaHGridFunctions.
                f_infinity_str += str(gf.f_infinity) + ","  # type: ignore
                f_wavespeed_str += str(gf.wavespeed) + ","  # type: ignore

            outstr += f"static const REAL gridfunctions_f_infinity[NUM_EVOL_GFS] = {{ {f_infinity_str[:-1] } }};\n"
            outstr += "\n\n// SET gridfunctions_wavespeed[i] = gridfunction i's characteristic wave speed:\n"
            outstr += f"static const REAL gridfunctions_wavespeed[NUM_EVOL_GFS] = {{ {f_wavespeed_str[:-1]} }};\n"
        return outstr


class BaseETGridFunction(GridFunction):
    """
    Subclass for basic (non-CarpetX) Einstein Toolkit grid functions.
    """

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
    ) -> None:
        super().__init__(name, group)
        self.rank = rank
        self.f_infinity = f_infinity
        self.wavespeed = wavespeed
        self.is_basename = is_basename
        self.c_type_alias = "CCTK_REAL"

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Verifies that the 'group' attribute of the GridFunction instance is one of the valid groups.

        The valid groups are:
        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints.

        If the 'group' attribute is not one of these, a ValueError is raised.

        Raises:
            ValueError: If the 'group' attribute is not one of the valid groups.
        """
        valid_groups = (
            "EVOL",
            "AUX",
            "AUXEVOL",
        )
        if self.group not in valid_groups:
            msg = (
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
                '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
                '    "AUX": for all other quantities needed at all gridpoints.\n'
            )
            raise ValueError(msg)

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        This function takes a grid function name, three integer offsets, and an optional keyword arguments dictionary,
        then returns a formatted string. If 'gf_array_name' is not provided in kwargs, it raises a ValueError.

        Parameters:
        gf_name (str): The grid function name
        i0_offset (int): Offset for the first index
        i1_offset (int): Offset for the second index
        i2_offset (int): Offset for the third index
        kwargs (dict): Optional keyword arguments.

        Returns:
        str: Formatted string

        Doctests:
        >>> glb_gridfcs_dict.clear()
        >>> par.set_parval_from_str("Infrastructure", "BaseET")
        >>> abc = register_gridfunctions("abc", group="EVOL")
        >>> glb_gridfcs_dict["abc"].read_gf_from_memory_Ccode_onept(1, 2, 3)
        'abc[CCTK_GFINDEX3D(cctkGH, i0+1, i1+2, i2+3)]'
        >>> defg = register_gridfunctions("defg", group="EVOL")
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0)
        'defg[CCTK_GFINDEX3D(cctkGH, i0, i1-1, i2)]'
        """

        if kwargs:
            raise ValueError(
                "BaseETGridFunction.read_gf_from_memory_Ccode_onept() does not accept kwargs!"
            )

        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

        return f"{self.name}[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"


class CarpetXGridFunction(GridFunction):
    """Class for CarpetX gridfunction operations."""

    def __init__(
        self,
        name: str,
        group: str = "",
        rank: int = 0,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
        centering: str = "C",
    ) -> None:
        super().__init__(name, group)
        self.rank = rank
        self.f_infinity = f_infinity
        self.wavespeed = wavespeed
        self.is_basename = is_basename
        self.centering = centering

        group_suffixes = {
            "EXTERNAL": "_ext",
            "CORE": "_core",
            "TILE_TMP": "_tile_tmp",
            "SCALAR_TMP": "",
        }
        if self.group in group_suffixes:
            self.name += group_suffixes[self.group]

        # Dictionaries for storing grid function groups.
        # self.index_group = {}
        # self.rev_index_group = {}

        # Validation step.
        self.verify_gridfunction_centering_is_valid()
        self.verify_gridfunction_group_is_valid()

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Verifies that the 'group' attribute of the GridFunction instance is one of the valid groups.

        The valid groups are:
        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints,
        'EXTERNAL': for all quantities defined in other modules,
        'CORE': for all quantities defined inside the CarpetX driver,
        'TILE_TMP': for all temporary quantities defined for CarpetX tiles,
        'SCALAR_TMP': for all temporary quantities defined for doubles.

        If the 'group' attribute is not one of these, a ValueError is raised.

        Raises:
            ValueError: If the 'group' attribute is not one of the valid groups.
        """
        valid_groups = (
            "EVOL",
            "AUX",
            "AUXEVOL",
            "EXTERNAL",
            "CORE",
            "TILE_TMP",
            "SCALAR_TMP",
        )
        if self.group not in valid_groups:
            msg = (
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
                '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
                '    "AUX": for all other quantities needed at all gridpoints.\n'
                '    "EXTERNAL": for all quantities defined in other modules.\n'
                '    "CORE": for all quantities defined inside the CarpetX driver.\n'
                '    "TILE_TMP": for all temporary quantities defined for CarpetX tiles.\n'
                '    "SCALAR_TMP": for all temporary quantities defined for doubles.'
            )
            raise ValueError(msg)

    def verify_gridfunction_centering_is_valid(self) -> None:
        """
        Verifies that the 'centering' attribute of the GridFunction instance is valid.

        A valid 'centering' is a 3-character string where each character is either 'C' or 'V'.
        'C' stands for cell-centered quantities, 'V' stands for vertex-centered quantities.

        If the 'centering' attribute is not valid, a ValueError is raised.

        Raises:
            ValueError: If the 'centering' attribute is not valid.
        """
        if len(self.centering) != 3 or any(c not in ("C", "V") for c in self.centering):
            msg = (
                f"Unsupported gridfunction centering {self.centering}. Supported centerings are 3-character strings, each character being either:\n"
                '    "C": for cell-centered quantities,\n'
                '    "V": for vertex-centered quantities.'
            )
            raise ValueError(msg)


# Contains a list of gridfunction objects.
glb_gridfcs_dict: Dict[
    str, Union[GridFunction, BHaHGridFunction, BaseETGridFunction, CarpetXGridFunction]
] = {}


#####################################
# Outside the classes above, only define functions in grid.py that are used by multiple infrastructures


def register_gridfunctions(
    names: Union[str, List[str]], **kwargs: Any
) -> List[sp.Symbol]:
    """
    Register grid functions with a specified name or list of names.

    :param names: A name or list of names of grid functions to be registered.
    :param kwargs: Additional keyword arguments.
    :return: A list of sympy symbols representing the registered grid functions.

    Example:
    >>> glb_gridfcs_dict.clear()
    >>> gridfunc = register_gridfunctions('gridfunc')[0]
    >>> print(gridfunc)
    gridfunc
    >>> # glb_gridfcs_dict["gridfunc"] is of type BHaHGridFunction. __dict__.items() access all variables in the class.
    >>> for key, value in glb_gridfcs_dict["gridfunc"].__dict__.items():
    ...     print(key, value)
    c_type_alias CCTK_REAL
    name gridfunc
    group EVOL
    rank 0
    f_infinity 0.0
    wavespeed 1.0
    is_basename True
    >>> gridfunc1, gridfunc2 = register_gridfunctions(['gridfunc1', 'gridfunc2'], f_infinity=[1.0, 4.0])
    >>> for key, value in glb_gridfcs_dict["gridfunc1"].__dict__.items():
    ...     print(key, value)
    c_type_alias CCTK_REAL
    name gridfunc1
    group EVOL
    rank 0
    f_infinity [1.0, 4.0]
    wavespeed 1.0
    is_basename True
    """
    # Step 1: Convert names to a list if it's not already a list
    if not isinstance(names, list):
        names = [names]

    # Step 2: Process each gridfunction
    sympy_symbol_list: List[sp.Symbol] = []
    Infrastructure = par.parval_from_str("Infrastructure")
    for i, name in enumerate(names):
        if name in glb_gridfcs_dict:
            print(f"Warning: Gridfunction {name} is already registered.")
        else:
            gf: Union[BHaHGridFunction, BaseETGridFunction, CarpetXGridFunction]
            if Infrastructure == "BHaH":
                kwargs_modify = kwargs.copy()
                for param in ["f_infinity", "wavespeed"]:
                    if kwargs.get(param) and isinstance(kwargs.get(param), list):
                        # mypy: Once again bonks out after I've CONFIRMED kwargs.get(param) is not None and is a list!
                        kwargs_modify[param] = kwargs.get(param)[i]  # type: ignore
                gf = BHaHGridFunction(name, **kwargs_modify)
            elif Infrastructure == "BaseET":
                gf = BaseETGridFunction(name, **kwargs)
            elif Infrastructure == "CarpetX":
                gf = CarpetXGridFunction(name, **kwargs)
            else:
                raise ValueError(f"Infrastructure = {Infrastructure} unknown")

            glb_gridfcs_dict[name] = gf
        sympy_symbol_list += [sp.symbols(name, real=True)]

    # Step 3: Return the list of sympy symbols, which can have any number of elements >= 1
    return sympy_symbol_list


def register_gridfunctions_for_single_rank1(
    basename: str, **kwargs: Any
) -> Sequence[sp.Expr]:
    """
    Register gridfunctions for a single rank 1 variable.

    Args:
        basename: Base name of the gridfunction to register.
        **kwargs: Additional keyword arguments for registration.

    Returns:
        List of SymPy symbols.

    Doctests:
    >>> glb_gridfcs_dict.clear()
    >>> register_gridfunctions_for_single_rank1("betU")
    [betU0, betU1, betU2]

    Now make sure this was stored to glb_gridfcs_dict correctly:
    >>> outstr = ""
    >>> for key, gfobj in glb_gridfcs_dict.items():
    ...    outstr += gfobj.name + " "
    >>> print(outstr[:-1])
    betU0 betU1 betU2
    """

    dimension = kwargs.get("dimension", 3)
    if "dimension" in kwargs:
        del kwargs["dimension"]

    # Step 1: Declare a list of SymPy variables,
    #         where IDX_OBJ_TMP[i] = gf_basename+str(i)
    IDX_OBJ_TMP = ixp.declarerank1(basename, dimension=dimension, **kwargs)

    # Step 2: Register each gridfunction
    gf_list = []
    for i in range(3):
        gf_list += [str(IDX_OBJ_TMP[i])]
    # Manually adjust kwargs so that basename is not checked (we're not passing the basename!),
    # and rank is set to 1.
    kwargs["is_basename"] = False
    kwargs["rank"] = 1
    register_gridfunctions(gf_list, **kwargs)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


def register_gridfunctions_for_single_rank2(
    basename: str, symmetry: Optional[str] = None, **kwargs: Any
) -> Sequence[Sequence[sp.Expr]]:
    """
    Register gridfunctions for a single rank 2 variable.

    Args:
        basename: Base name of the gridfunction to register.
        symmetry: Optional symmetry for rank 2 gridfunctions.
        **kwargs: Additional keyword arguments for registration.

    Returns:
        2D list of SymPy symbols.

    Doctests:
    >>> glb_gridfcs_dict.clear()
    >>> register_gridfunctions_for_single_rank2("gDD", symmetry="sym01")
    [[gDD00, gDD01, gDD02], [gDD01, gDD11, gDD12], [gDD02, gDD12, gDD22]]

    Now make sure this was stored to glb_gridfcs_dict correctly:
    >>> outstr = ""
    >>> for key, gfobj in glb_gridfcs_dict.items():
    ...    outstr += gfobj.name + " "
    >>> print(outstr[:-1])
    gDD00 gDD01 gDD02 gDD11 gDD12 gDD22
    """

    dimension = kwargs.get("dimension", 3)
    if "dimension" in kwargs:
        del kwargs["dimension"]

    # Step 1: Declare a list of lists of SymPy variables,
    #         where IDX_OBJ_TMP[i][j] = gf_basename+str(i)+str(j)
    IDX_OBJ_TMP = ixp.declarerank2(
        basename, symmetry=symmetry, dimension=dimension, **kwargs
    )

    # Step 2: register each gridfunction, being careful not
    #         not to store duplicates due to rank-2 symmetries.
    gf_list: List[str] = []
    for i in range(dimension):
        for j in range(dimension):
            save = True
            for gf in gf_list:
                if gf == str(IDX_OBJ_TMP[i][j]):
                    save = False
            if save:
                gf_list.append(str(IDX_OBJ_TMP[i][j]))

    # Manually adjust kwargs so that basename is not checked (we're not passing the basename!),
    # and rank is set to 2.
    kwargs["is_basename"] = False
    kwargs["rank"] = 2
    register_gridfunctions(gf_list, **kwargs)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


# Outside the classes above, only define functions in grid.py that are used by multiple infrastructures
#####################################


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

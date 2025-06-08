"""
Provide registration functions for numerical grid parameters.

This module handles grid parameters and grid functions across various NRPy-supported infrastructures.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple, Union, cast

import sympy as sp
from typing_extensions import Literal

import nrpy.indexedexp as ixp
import nrpy.params as par
from nrpy.helpers.type_annotation_utilities import validate_literal_arguments

centerings = Literal[
    "CCC",
    "CCV",
    "CVC",
    "CVV",
    "VCC",
    "VCV",
    "VVC",
    "VVV",
    "CC",
    "CV",
    "VC",
    "VV",
    "V",
    "C",
]

Cart_origin = par.register_CodeParameters(
    "REAL", __name__, ["Cart_originx", "Cart_originy", "Cart_originz"], 0.0
)
_ = par.register_CodeParameter(
    "int", __name__, "NUMGRIDS", 1, commondata=True, add_to_parfile=False
)
_ = par.register_CodeParameter(
    "bool", __name__, "grid_rotates", False, commondata=False, add_to_parfile=False
)
_ = par.register_CodeParameter(
    "#define",
    __name__,
    "MAXNUMGRIDS",
    "15",
    commondata=False,
    add_to_set_CodeParameters_h=False,
    add_to_parfile=False,
)


def _format_c_offset_str(var: str, offset: int) -> str:
    """
    Format a C-style variable with an offset, matching original output.

    :param var: The variable name string (e.g., "i0").
    :param offset: The integer offset.
    :return: A string like "i0", "i0+1", or "i0-1".
    """
    if offset == 0:
        return var
    # The .replace() is needed to match the original's output exactly.
    return f"{var}+{offset}".replace("+-", "-")


def _flatten_and_unique_str(nested_list: List[Any]) -> List[str]:
    """
    Flatten a nested list into a single list of unique strings, preserving order.

    :param nested_list: The list to be flattened.
    :return: An order-preserved list containing unique string representations of the elements.
    """
    flat_list: List[str] = []
    seen: set[str] = set()

    def recurse(sub_list: Any) -> None:
        """
        Recursive helper to traverse the nested structure.

        :param sub_list: The list or item to process.
        """
        if isinstance(sub_list, list):
            for item in sub_list:
                recurse(item)
        else:
            name = str(sub_list)
            if name not in seen:
                seen.add(name)
                flat_list.append(name)

    recurse(nested_list)
    return flat_list


# Core GridFunction class.
class GridFunction:
    """The core class for grid functions."""

    # Parity conditions for rank-2 tensors, organized by dimension.
    _PARITY_CONDITIONS_RANK2_DIM3: Dict[Tuple[str, str], int] = {
        ("0", "0"): 4,
        ("0", "1"): 5,
        ("1", "0"): 5,
        ("0", "2"): 6,
        ("2", "0"): 6,
        ("1", "1"): 7,
        ("1", "2"): 8,
        ("2", "1"): 8,
        ("2", "2"): 9,
    }
    _PARITY_CONDITIONS_RANK2_DIM4: Dict[Tuple[str, str], int] = {
        ("0", "0"): 0,  # scalar tt component
        ("0", "1"): 1,
        ("1", "0"): 1,
        ("0", "2"): 2,
        ("2", "0"): 2,  # vector ti
        ("0", "3"): 3,
        ("3", "0"): 3,
        ("1", "1"): 4,
        ("1", "2"): 5,
        ("2", "1"): 5,
        ("1", "3"): 6,  # tensor ij
        ("3", "1"): 6,
        ("2", "2"): 7,
        ("2", "3"): 8,
        ("3", "2"): 8,
        ("3", "3"): 9,
    }

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        dimension: int = 3,
        gf_type: str = "REAL",
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
    ) -> None:
        self.name: str = name
        self.group: str = group
        self.rank: int = rank
        self.dimension: int = dimension
        self.gf_type: str = gf_type
        self.f_infinity: float = f_infinity
        self.wavespeed: float = wavespeed
        self.is_basename: bool = is_basename

        if is_basename:
            self.verify_gridfunction_basename_is_valid(self.name)

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        Catch-all function for gridfunctions that haven't been set up correctly.

        This function serves as a placeholder for the specific `read_gf_from_memory_Ccode_onept` implementations within
        gridfunction classes. It should be overridden in specific gridfunction classes to generate C code for reading
        a single point of a gridfunction from memory, taking into account any provided offsets.

        :param i0_offset: The offset in the 0th (e.g., x) direction from the current grid point. Default is 0.
        :param i1_offset: The offset in the 1st (e.g., y) direction from the current grid point. Default is 0.
        :param i2_offset: The offset in the 2nd (e.g., z) direction from the current grid point. Default is 0.
        :param kwargs: Additional keyword arguments that might be necessary for specific implementations.
        :return: A string representing the C code required to read a gridfunction's value from memory at the specified offsets.
                 As this is a catch-all function, it returns a message prompting to define the function properly within grid.py.

        Example:
            Assuming an instance `gf` of a GridFunction class with properly overridden `read_gf_from_memory_Ccode_onept` method,
            calling `gf.read_gf_from_memory_Ccode_onept(i0_offset=1)` would return the C code as a string for reading the gridfunction
            value at an offset of 1 in the x direction from the current point.
        """
        return f"Please define read_gf_from_memory_Ccode_onept() inside grid.py. Inputs: [class instance] {i0_offset} {i1_offset} {i2_offset} {kwargs}"

    def verify_gridfunction_basename_is_valid(self, name: str) -> None:
        """
        Validate the gridfunction's base name.

        Gridfunction base names must not end with integers. For instance, a rank-1
        gridfunction 'vecU1' has a valid base name 'vecU'. A scalar gridfunction
        named 'u2' would have an invalid base name 'u2'. This strict requirement
        facilitates quick identification of a gridfunction by its name.

        :param name: The name of the gridfunction to validate.

        :raises ValueError: If the name is of zero length or ends with an integer, indicating
                            it does not meet the naming conventions for gridfunctions.
        :raises TypeError: If the name is not of type str, ensuring that gridfunction names
                           are always represented by strings.

        Example:
            >>> try: GF = GridFunction(name="h1")
            ... except: "Attempt to create a gridfunction ending in an integer failed. Good."
            'Attempt to create a gridfunction ending in an integer failed. Good.'
        """
        if not isinstance(name, str):
            raise TypeError("Gridfunction names must be strings")
        # First check for zero-length basenames:
        if not name:
            raise ValueError("Tried to register gridfunction without a name!")

        if name[-1].isdigit():
            raise ValueError(
                f"Tried to register gridfunction with base name: {name}\n"
                "To ensure that finite difference code generations do not get confused,\n"
                "gridfunctions with base names ending in an integer are forbidden; pick a new name."
            )

    @staticmethod
    def gridfunction_lists() -> Tuple[List[str], List[str], List[str]]:
        """
        Generate sorted lists of gridfunction names for each group type: 'EVOL', 'AUX', and 'AUXEVOL'.

        This function creates a dictionary with group types as keys and corresponding gridfunction names
        as values. It iterates through a global dictionary of gridfunctions, 'glb_gridfcs_dict', and
        appends the gridfunction names to the corresponding group lists in 'groups'. The function then
        sorts each list and returns the sorted lists as a tuple.

        Note: This function assumes the existence of a global gridfunction dictionary named
        'glb_gridfcs_dict' whose keys are gridfunction names and values are gridfunction objects. Each
        gridfunction object should have 'group' and 'name' attributes.

        :returns: A tuple containing lists of gridfunction names
                  for each group: 'EVOL', 'AUX', and 'AUXEVOL', respectively.

        Example:
            gridfunction_lists()
            (['evol_gf1', 'evol_gf2'], ['aux_gf1', 'aux_gf2'], ['auxevol_gf1', 'auxevol_gf2'])
        """
        # Using list comprehensions for a concise and Pythonic implementation.
        evol_gfs = sorted(
            [gf.name for gf in glb_gridfcs_dict.values() if gf.group == "EVOL"],
            key=str.lower,
        )
        aux_gfs = sorted(
            [gf.name for gf in glb_gridfcs_dict.values() if gf.group == "AUX"],
            key=str.lower,
        )
        auxevol_gfs = sorted(
            [gf.name for gf in glb_gridfcs_dict.values() if gf.group == "AUXEVOL"],
            key=str.lower,
        )
        return evol_gfs, aux_gfs, auxevol_gfs

    @staticmethod
    def set_parity_types(list_of_gf_names: List[str]) -> List[int]:
        """
        Set the parity types for a given list of grid function names.

        This method looks up each grid function name in the global grid function dictionary
        (`glb_gridfcs_dict`) to determine and set the appropriate parity type based on the
        grid function's rank and spatial dimension. The parity type is an integer that encodes
        how the grid function transforms under spatial inversion (parity transformation).

        :param list_of_gf_names: List of grid function names for which to set the parity types.
        :return: A list of integers representing the parity types for the grid functions.
        :raises ValueError: If unable to determine the parity type for any of the grid functions
                            based on their names, ranks, and dimensions, or if the number of
                            determined parity types does not match the number of input grid function names.
        """
        parity_type_list: List[int] = []
        for name in list_of_gf_names:
            gf = glb_gridfcs_dict.get(name)
            if not gf:
                # This emulates the original's behavior of failing the length check at the end.
                continue

            parity_val: Optional[int] = None
            if gf.rank == 0:
                parity_val = 0
            elif gf.rank == 1:
                parity_val = int(gf.name[-1]) + 1
            elif gf.rank == 2:
                indices = (gf.name[-2], gf.name[-1])
                if gf.dimension == 3:
                    parity_val = GridFunction._PARITY_CONDITIONS_RANK2_DIM3.get(indices)
                elif gf.dimension == 4:
                    parity_val = GridFunction._PARITY_CONDITIONS_RANK2_DIM4.get(indices)

            if parity_val is not None:
                parity_type_list.append(parity_val)
            else:
                raise ValueError(
                    f"Error: Could not figure out parity type for {gf.group} gridfunction: {gf.name}, {gf.name[-2]}, {gf.name[-1]}, {gf.rank}, {gf.dimension}"
                )

        if len(parity_type_list) != len(list_of_gf_names):
            raise ValueError(
                "Error: For some reason the length of the parity types list did not match the length of the gf list."
            )
        return parity_type_list


class BHaHGridFunction(GridFunction):
    """The subclass for BlackHoles@Home grid functions."""

    VALID_GROUPS: Tuple[str, ...] = ("EVOL", "AUX", "AUXEVOL")
    GROUP_DESCRIPTIONS: str = (
        '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
        '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
        '    "AUX": for all other quantities needed at all gridpoints.\n'
    )

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        dimension: int = 3,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
        gf_array_name: str = "use_in_gfs_for_EVOL_auxevol_gfs_for_AUXEVOL_etc",
        sync_gf_in_superB: Optional[bool] = None,
    ) -> None:
        super().__init__(
            name, group, rank, dimension, "REAL", f_infinity, wavespeed, is_basename
        )
        self.verify_gridfunction_group_is_valid()
        # Adhere to original magic string default to preserve doctest compatibility
        if gf_array_name == "use_in_gfs_for_EVOL_auxevol_gfs_for_AUXEVOL_etc":
            if group == "EVOL":
                self.gf_array_name = "in_gfs"
            elif group == "AUXEVOL":
                self.gf_array_name = "auxevol_gfs"
            elif group == "AUX":
                self.gf_array_name = "aux_gfs"
        else:
            self.gf_array_name = gf_array_name

        self.sync_gf_in_superB = sync_gf_in_superB
        if sync_gf_in_superB is None and group in ("EVOL", "AUX"):
            self.sync_gf_in_superB = True

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Validate the gridfunction group.

        The valid groups are 'EVOL', 'AUX' and 'AUXEVOL'.

        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints.

        :raises ValueError: If the group is not valid.

        Doctests:
        >>> try: GF = BHaHGridFunction(name="blah", group="INVALID")
        ... except: print("Errored out. This is good.")
        Errored out. This is good.
        """
        if self.group not in self.VALID_GROUPS:
            raise ValueError(
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                f"{self.GROUP_DESCRIPTIONS}"
            )

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        Retrieve a grid function value from memory for a given offset.

        :param i0_offset: Offset for the first index.
        :param i1_offset: Offset for the second index.
        :param i2_offset: Offset for the third index.
        :param kwargs: Optional keyword arguments, expected to have 'gf_array_name'.

        :return: Formatted string.

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
        gf_array_name = kwargs.get("gf_array_name", self.gf_array_name)
        access_string = self.access_gf(
            self.name, i0_offset, i1_offset, i2_offset, gf_array_name
        )

        if kwargs.get("enable_simd"):
            return f"ReadSIMD(&{access_string})"
        return access_string

    @staticmethod
    def access_gf(
        gf_name: str,
        i0_offset: int = 0,
        i1_offset: int = 0,
        i2_offset: int = 0,
        gf_array_name: str = "in_gfs",
    ) -> str:
        """
        Retrieve a grid function value from memory for a given offset.

        :param gf_name: The grid function name.
        :param i0_offset: Offset for the first index.
        :param i1_offset: Offset for the second index.
        :param i2_offset: Offset for the third index.
        :param gf_array_name: Optional grid function array name.

        :return: Formatted string.

        Doctests:
        >>> BHaHGridFunction.access_gf("aa", 1,2,3)
        'in_gfs[IDX4(AAGF, i0+1, i1+2, i2+3)]'
        >>> BHaHGridFunction.access_gf("defg", 0, -1, 0, "My_Array")
        'My_Array[IDX4(DEFGGF, i0, i1-1, i2)]'
        """
        i0 = _format_c_offset_str("i0", i0_offset)
        i1 = _format_c_offset_str("i1", i1_offset)
        i2 = _format_c_offset_str("i2", i2_offset)

        return f"{gf_array_name}[IDX4({gf_name.upper()}GF, {i0}, {i1}, {i2})]"

    @staticmethod
    def gridfunction_defines() -> str:
        """
        Generate a string representation of all the GridFunction definitions.

        This function also includes definitions for the values of each gridfunction
        at infinity and the characteristic wave speed. The output is formatted to be
        used in C-style languages.

        This function assumes the existence of a global gridfunction dictionary named
        'glb_gridfcs_dict'. The keys of this dictionary should be gridfunction names,
        and the values should be instances of classes derived from GridFunction.

        :return: The resulting string representation of the GridFunction definitions.
        """

        def define_gfs_group(name: str, gfs: List[str]) -> str:
            """
            Generate C-style #define string for a list of GridFunctions in a group.

            :param name: The name of the group (e.g., "EVOL").
            :param gfs: The list of grid function names in the group.
            :return: A formatted C-style string with #defines.
            """
            num_gfs = len(gfs)
            defines = "\n".join(
                f"#define {gf.upper()}GF\t{i}" for i, gf in enumerate(gfs)
            )
            return f"\n// {name.upper()} VARIABLES:\n#define NUM_{name.upper()}_GFS {num_gfs}\n{defines}\n"

        evol, aux, auxevol = BHaHGridFunction.gridfunction_lists()

        # Start with the EVOL group defines
        outstr = define_gfs_group("EVOL", evol)

        # Append f_infinity and wavespeed definitions if EVOL variables exist
        if evol:
            f_infinity_list = [str(glb_gridfcs_dict[var].f_infinity) for var in evol]
            f_wavespeed_list = [str(glb_gridfcs_dict[var].wavespeed) for var in evol]

            f_infinity_str = ", ".join(f_infinity_list)
            f_wavespeed_str = ", ".join(f_wavespeed_list)

            outstr += "\n\n// SET gridfunctions_f_infinity[i] = evolved gridfunction i's value in the limit r->infinity:\n"
            outstr += f"static const REAL gridfunctions_f_infinity[NUM_EVOL_GFS] = {{ {f_infinity_str} }};\n"

            outstr += "\n\n// SET gridfunctions_wavespeed[i] = evolved gridfunction i's characteristic wave speed:\n"
            outstr += f"static const REAL gridfunctions_wavespeed[NUM_EVOL_GFS] = {{ {f_wavespeed_str} }};\n"

        # Append AUX and AUXEVOL group defines
        outstr += define_gfs_group("AUX", aux)
        outstr += define_gfs_group("AUXEVOL", auxevol)

        return outstr


class ETLegacyGridFunction(GridFunction):
    """Subclass for basic (non-CarpetX) Einstein Toolkit grid functions."""

    VALID_GROUPS: Tuple[str, ...] = ("EVOL", "AUX", "AUXEVOL")
    GROUP_DESCRIPTIONS: str = (
        '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
        '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
        '    "AUX": for all other quantities needed at all gridpoints.\n'
    )

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        dimension: int = 3,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
        gf_array_name: str = "",
    ) -> None:
        super().__init__(
            name,
            group,
            rank,
            dimension,
            "CCTK_REAL",
            f_infinity,
            wavespeed,
            is_basename,
        )
        # The following local variable is unused, but preserved to match the original
        # code's behavior and pass a doctest that inspects __dict__.
        _gf_array_name = gf_array_name
        self.verify_gridfunction_group_is_valid()

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Verify that the 'group' attribute of the GridFunction instance is one of the valid groups.

        The valid groups are:
        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints.

        :raises ValueError: If the 'group' attribute is not one of the valid groups.
        """
        if self.group not in self.VALID_GROUPS:
            raise ValueError(
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                f"{self.GROUP_DESCRIPTIONS}"
            )

    @staticmethod
    def access_gf(
        gf_name: str,
        i0_offset: int = 0,
        i1_offset: int = 0,
        i2_offset: int = 0,
        use_GF_suffix: bool = True,
    ) -> str:
        """
        Retrieve a grid function value from memory for a given offset.

        :param gf_name: The grid function name.
        :param i0_offset: Offset for the first index.
        :param i1_offset: Offset for the second index.
        :param i2_offset: Offset for the third index.
        :param use_GF_suffix: Suffix gridfunction name with GF (default: True).

        :return: Formatted string.

        Doctests:
        >>> ETLegacyGridFunction.access_gf("aa", 1,2,3)
        'aaGF[CCTK_GFINDEX3D(cctkGH, i0+1, i1+2, i2+3)]'
        >>> ETLegacyGridFunction.access_gf("defg", 0, -1, 0)
        'defgGF[CCTK_GFINDEX3D(cctkGH, i0, i1-1, i2)]'
        """
        i0 = _format_c_offset_str("i0", i0_offset)
        i1 = _format_c_offset_str("i1", i1_offset)
        i2 = _format_c_offset_str("i2", i2_offset)
        gf_name_str = f"{gf_name}GF" if use_GF_suffix else gf_name

        return f"{gf_name_str}[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        Generate a formatted string using the grid function name and offsets.

        If 'gf_array_name' is not provided in kwargs, a ValueError is raised.

        :param i0_offset: Offset for the first index
        :param i1_offset: Offset for the second index
        :param i2_offset: Offset for the third index
        :param kwargs: Optional keyword arguments.
        :return: Formatted string

        Doctests:
        >>> glb_gridfcs_dict.clear()
        >>> par.set_parval_from_str("Infrastructure", "ETLegacy")
        >>> abc = register_gridfunctions("abc", group="EVOL")
        >>> glb_gridfcs_dict["abc"].read_gf_from_memory_Ccode_onept(1, 2, 3)
        'abcGF[CCTK_GFINDEX3D(cctkGH, i0+1, i1+2, i2+3)]'
        >>> defg = register_gridfunctions("defg", group="EVOL")
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0, enable_simd=True)
        'ReadSIMD(&defgGF[CCTK_GFINDEX3D(cctkGH, i0, i1-1, i2)])'
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0, enable_simd=True, use_GF_suffix=False)
        'ReadSIMD(&defg[CCTK_GFINDEX3D(cctkGH, i0, i1-1, i2)])'
        """
        use_GF_suffix = kwargs.get("use_GF_suffix", True)
        access_str = self.access_gf(
            self.name, i0_offset, i1_offset, i2_offset, use_GF_suffix
        )
        if kwargs.get("enable_simd"):
            return f"ReadSIMD(&{access_str})"
        return access_str


class CarpetXGridFunction(GridFunction):
    """Class for CarpetX gridfunction operations."""

    VALID_GROUPS: Tuple[str, ...] = (
        "EVOL",
        "AUX",
        "AUXEVOL",
        "EXTERNAL",
        "CORE",
        "TILE_TMP",
        "SCALAR_TMP",
    )
    GROUP_DESCRIPTIONS: str = (
        '    "EVOL": for evolved quantities (i.e., quantities stepped forward in time),\n'
        '    "AUXEVOL": for auxiliary quantities needed at all points by evolved quantities,\n'
        '    "AUX": for all other quantities needed at all gridpoints.\n'
        '    "EXTERNAL": for all quantities defined in other modules.\n'
        '    "CORE": for all quantities defined inside the CarpetX driver.\n'
        '    "TILE_TMP": for all temporary quantities defined for CarpetX tiles.\n'
        '    "SCALAR_TMP": for all temporary quantities defined for doubles.'
    )

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        dimension: int = 3,
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
        centering: centerings = "CCC",
        gf_array_name: str = "",
        thorn: str = "Cactus",
    ) -> None:
        super().__init__(
            name,
            group,
            rank,
            dimension,
            "CCTK_REAL",
            f_infinity,
            wavespeed,
            is_basename,
        )
        validate_literal_arguments()
        self.thorn = thorn
        self.centering = centering
        # The following local variable is unused, but preserved to match the original
        # code's behavior and pass a doctest that inspects __dict__.
        _gf_array_name = gf_array_name

        # Validation steps.
        self.verify_gridfunction_group_is_valid()
        self.verify_gridfunction_centering_is_valid()

        group_suffixes = {
            "EXTERNAL": "_ext",
            "CORE": "_core",
            "TILE_TMP": "_tile_tmp",
        }
        self.name += group_suffixes.get(self.group, "")

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Verify that the 'group' attribute of the GridFunction instance is one of the valid groups.

        The valid groups are:
        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints,
        'EXTERNAL': for all quantities defined in other modules,
        'CORE': for all quantities defined inside the CarpetX driver,
        'TILE_TMP': for all temporary quantities defined for CarpetX tiles,
        'SCALAR_TMP': for all temporary quantities defined for doubles.

        :raises ValueError: If the 'group' attribute is not one of the valid groups.
        """
        if self.group not in self.VALID_GROUPS:
            raise ValueError(
                f"Unsupported gridfunction group {self.group}. Supported groups include:\n"
                f"{self.GROUP_DESCRIPTIONS}"
            )

    def verify_gridfunction_centering_is_valid(self) -> None:
        """
        Verify that the 'centering' attribute of the GridFunction instance is valid.

        A valid 'centering' is a 3-character string where each character is either 'C' or 'V'.
        'C' stands for cell-centered quantities, 'V' stands for vertex-centered quantities.

        :raises ValueError: If the 'centering' attribute is not valid.
        """
        if len(self.centering) != 3 or any(c not in ("C", "V") for c in self.centering):
            raise ValueError(
                f"Unsupported gridfunction centering {self.centering}. Supported centerings are 3-character strings, each character being either:\n"
                '    "C": for cell-centered quantities,\n'
                '    "V": for vertex-centered quantities.'
            )

    @staticmethod
    def access_gf(
        gf_name: str,
        i0_offset: int = 0,
        i1_offset: int = 0,
        i2_offset: int = 0,
        use_GF_suffix: bool = True,
        reuse_index: bool = False,
        index_name: str = "",
    ) -> str:
        """
        Retrieve a grid function value from memory for a given offset.

        :param gf_name: The grid function name.
        :param i0_offset: Offset for the first index.
        :param i1_offset: Offset for the second index.
        :param i2_offset: Offset for the third index.
        :param use_GF_suffix: Suffix gridfunction name with GF (default: True).
        :param reuse_index: Use provided grid function index instead of computing it directly (default: False).
        :param index_name: Name of predefined index variable.

        :return: Formatted string.

        Doctests:
        >>> CarpetXGridFunction.access_gf("aa",1,2,3)
        'aaGF(p.I + 1*p.DI[0] + 2*p.DI[1] + 3*p.DI[2])'
        >>> CarpetXGridFunction.access_gf("defg", 0, -1, 0)
        'defgGF(p.I - 1*p.DI[1])'
        """
        index = index_name
        if not reuse_index:
            # This logic is reverted to match the original's exact whitespace output for doctest compatibility.
            index = "p.I"
            if i0_offset != 0:
                index += f" + {i0_offset}*p.DI[0]".replace("+ -", "- ")
            if i1_offset != 0:
                index += f" + {i1_offset}*p.DI[1]".replace("+ -", "- ")
            if i2_offset != 0:
                index += f" + {i2_offset}*p.DI[2]".replace("+ -", "- ")

        gf_name_str = f"{gf_name}GF" if use_GF_suffix else gf_name
        return f"{gf_name_str}({index})"

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """
        Generate a formatted string using the grid function name and offsets.

        :param i0_offset: Offset for the first index
        :param i1_offset: Offset for the second index
        :param i2_offset: Offset for the third index
        :param kwargs: Optional keyword arguments.
        :return: Formatted string

        Doctests:
        >>> glb_gridfcs_dict.clear()
        >>> par.set_parval_from_str("Infrastructure", "CarpetX")
        >>> abc = register_gridfunctions("abc", group="EVOL")
        >>> glb_gridfcs_dict["abc"].read_gf_from_memory_Ccode_onept(1, 2, 3)
        'abcGF(p.I + 1*p.DI[0] + 2*p.DI[1] + 3*p.DI[2])'
        >>> defg = register_gridfunctions("defg", group="EVOL")
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0, enable_simd=True)
        'ReadSIMD(&defgGF(p.I - 1*p.DI[1]))'
        >>> glb_gridfcs_dict["defg"].read_gf_from_memory_Ccode_onept(0, -1, 0, enable_simd=True, use_GF_suffix=False)
        'ReadSIMD(&defg(p.I - 1*p.DI[1]))'
        """
        use_GF_suffix = kwargs.get("use_GF_suffix", True)
        reuse_index = kwargs.get("reuse_index", False)
        index_name = kwargs.get("index_name", "")

        access_str = self.access_gf(
            self.name,
            i0_offset,
            i1_offset,
            i2_offset,
            use_GF_suffix,
            reuse_index,
            index_name,
        )

        if kwargs.get("enable_simd"):
            return f"ReadSIMD(&{access_str})"
        return access_str


# Type alias for grid function objects.
GridFunctionType = Union[
    GridFunction, BHaHGridFunction, ETLegacyGridFunction, CarpetXGridFunction
]

# Global dictionary of registered grid functions.
glb_gridfcs_dict: Dict[str, GridFunctionType] = {}

# Factory mapping for grid function classes based on infrastructure.
GF_CLASS_MAP: Dict[str, type[GridFunctionType]] = {
    "BHaH": BHaHGridFunction,
    "ETLegacy": ETLegacyGridFunction,
    "CarpetX": CarpetXGridFunction,
}

# Factory mapping for ixp rank-N declaration functions.
IXP_RANK_FUNC_MAP: Dict[int, Callable[..., Any]] = {
    1: ixp.declarerank1,
    2: ixp.declarerank2,
    3: ixp.declarerank3,
    4: ixp.declarerank4,
}

#####################################
# Outside the classes above, only define functions in grid.py that are used by multiple infrastructures


def register_gridfunctions(
    names: Union[str, List[str]], dimension: int = 3, **kwargs: Any
) -> List[sp.Symbol]:
    """
    Register grid functions and return corresponding SymPy symbols.

    :param names: Grid function name(s) to register.
    :param dimension: Spatial dimension, default is 3.
    :param kwargs: Properties like 'group', 'f_infinity', 'wavespeed', 'is_basename'.
    :return: List of SymPy symbols for the registered grid functions.
    :raises ValueError: If the infrastructure is not recognized.

    Doctests:
    >>> glb_gridfcs_dict.clear()
    >>> par.set_parval_from_str("Infrastructure", "BHaH")
    >>> gridfunc = register_gridfunctions('gridfunc')[0]
    >>> print(gridfunc)
    gridfunc
    >>> for key, value in glb_gridfcs_dict["gridfunc"].__dict__.items():
    ...     print(key, value)
    name gridfunc
    group EVOL
    rank 0
    dimension 3
    gf_type REAL
    f_infinity 0.0
    wavespeed 1.0
    is_basename True
    gf_array_name in_gfs
    sync_gf_in_superB True
    >>> par.set_parval_from_str("Infrastructure", "ETLegacy")
    >>> gridfunc1, gridfunc2 = register_gridfunctions(['gridfunc1', 'gridfunc2'], f_infinity=[1.0, 4.0], is_basename=False)
    >>> for key, value in glb_gridfcs_dict["gridfunc1"].__dict__.items():
    ...     print(key, value)
    name gridfunc1
    group EVOL
    rank 0
    dimension 3
    gf_type CCTK_REAL
    f_infinity 1.0
    wavespeed 1.0
    is_basename False
    """
    names_list = names if isinstance(names, list) else [names]

    infrastructure = par.parval_from_str("Infrastructure")
    gf_class = GF_CLASS_MAP.get(infrastructure)
    if not gf_class:
        raise ValueError(f"Infrastructure = {infrastructure} unknown")

    sympy_symbol_list: List[sp.Symbol] = []
    for i, name in enumerate(names_list):
        if name in glb_gridfcs_dict:
            print(f"Warning: Gridfunction {name} is already registered.")
        else:
            # For parameters that can be lists, select the i-th element for this gridfunction.
            kwargs_individual = kwargs.copy()
            for param in ["f_infinity", "wavespeed"]:
                if isinstance(kwargs.get(param), list):
                    kwargs_individual[param] = kwargs[param][i]

            glb_gridfcs_dict[name] = gf_class(
                name, dimension=dimension, **kwargs_individual
            )
        sympy_symbol_list.append(sp.symbols(name, real=True))

    return sympy_symbol_list


def register_gridfunctions_for_single_rank1(
    basename: str, dimension: int = 3, **kwargs: Any
) -> Sequence[sp.Expr]:
    """
    Register rank-1 gridfunction components based on a basename.

    :param basename: Base name for gridfunction components.
    :param dimension: Spatial dimension, corresponds to number of components. Default is 3.
    :param kwargs: Additional options for registration.
    :return: List of SymPy symbols for the components.

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
    # This is a convenience wrapper for the more general rank-N function.
    return register_gridfunctions_for_single_rankN(
        basename, rank=1, dimension=dimension, **kwargs
    )


def register_gridfunctions_for_single_rank2(
    basename: str, symmetry: Optional[str] = None, dimension: int = 3, **kwargs: Any
) -> Sequence[Sequence[sp.Expr]]:
    """
    Register gridfunctions for a single rank 2 variable.

    :param basename: Base name of the gridfunction to register.
    :param symmetry: Optional symmetry for rank 2 gridfunctions.
    :param dimension: Spatial dimension. Default is 3.
    :param kwargs: Additional keyword arguments for registration.
    :return: 2D list of SymPy symbols.

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
    # This is a convenience wrapper for the more general rank-N function.
    return register_gridfunctions_for_single_rankN(
        basename, rank=2, symmetry=symmetry, dimension=dimension, **kwargs
    )


def register_gridfunctions_for_single_rankN(
    basename: str,
    rank: int,
    symmetry: Optional[str] = None,
    dimension: int = 3,
    **kwargs: Any,
) -> Union[
    List[sp.Symbol],
    List[List[sp.Symbol]],
    List[List[List[sp.Symbol]]],
    List[List[List[List[sp.Symbol]]]],
]:
    """
    Register gridfunctions for a single variable of arbitrary rank.

    :param basename: Base name of the gridfunction to register.
    :param rank: Rank of the tensor for which to register gridfunctions.
    :param symmetry: Optional symmetry for the gridfunctions. Defaults to None.
    :param dimension: Dimension of the space. Defaults to 3.
    :param kwargs: Additional keyword arguments for gridfunction registration.
    :return: Nested list of SymPy symbols representing the tensor gridfunctions.
    :raises ValueError: If the rank is not between 1 and 4, inclusive.

    Doctests:
    >>> glb_gridfcs_dict.clear()
    >>> register_gridfunctions_for_single_rankN("g", rank=2, symmetry="sym01")
    [[g00, g01, g02], [g01, g11, g12], [g02, g12, g22]]
    >>> print(sorted(list(glb_gridfcs_dict)))
    ['g00', 'g01', 'g02', 'g11', 'g12', 'g22']
    >>> print(glb_gridfcs_dict['g00'].__dict__)
    {'name': 'g00', 'group': 'EVOL', 'rank': 2, 'dimension': 3, 'gf_type': 'CCTK_REAL', 'f_infinity': 0.0, 'wavespeed': 1.0, 'is_basename': False}

    >>> glb_gridfcs_dict.clear()
    >>> register_gridfunctions_for_single_rankN("A", rank=1)
    [A0, A1, A2]
    >>> print(sorted(list(glb_gridfcs_dict)))
    ['A0', 'A1', 'A2']
    >>> print(glb_gridfcs_dict['A2'].__dict__)
    {'name': 'A2', 'group': 'EVOL', 'rank': 1, 'dimension': 3, 'gf_type': 'CCTK_REAL', 'f_infinity': 0.0, 'wavespeed': 1.0, 'is_basename': False}

    >>> glb_gridfcs_dict.clear()
    >>> register_gridfunctions_for_single_rankN("R", rank=4, symmetry="sym01_sym23", dimension=4)
    [[[[R0000, R0001, R0002, R0003], [R0001, R0011, R0012, R0013], [R0002, R0012, R0022, R0023], [R0003, R0013, R0023, R0033]], [[R0100, R0101, R0102, R0103], [R0101, R0111, R0112, R0113], [R0102, R0112, R0122, R0123], [R0103, R0113, R0123, R0133]], [[R0200, R0201, R0202, R0203], [R0201, R0211, R0212, R0213], [R0202, R0212, R0222, R0223], [R0203, R0213, R0223, R0233]], [[R0300, R0301, R0302, R0303], [R0301, R0311, R0312, R0313], [R0302, R0312, R0322, R0323], [R0303, R0313, R0323, R0333]]], [[[R0100, R0101, R0102, R0103], [R0101, R0111, R0112, R0113], [R0102, R0112, R0122, R0123], [R0103, R0113, R0123, R0133]], [[R1100, R1101, R1102, R1103], [R1101, R1111, R1112, R1113], [R1102, R1112, R1122, R1123], [R1103, R1113, R1123, R1133]], [[R1200, R1201, R1202, R1203], [R1201, R1211, R1212, R1213], [R1202, R1212, R1222, R1223], [R1203, R1213, R1223, R1233]], [[R1300, R1301, R1302, R1303], [R1301, R1311, R1312, R1313], [R1302, R1312, R1322, R1323], [R1303, R1313, R1323, R1333]]], [[[R0200, R0201, R0202, R0203], [R0201, R0211, R0212, R0213], [R0202, R0212, R0222, R0223], [R0203, R0213, R0223, R0233]], [[R1200, R1201, R1202, R1203], [R1201, R1211, R1212, R1213], [R1202, R1212, R1222, R1223], [R1203, R1213, R1223, R1233]], [[R2200, R2201, R2202, R2203], [R2201, R2211, R2212, R2213], [R2202, R2212, R2222, R2223], [R2203, R2213, R2223, R2233]], [[R2300, R2301, R2302, R2303], [R2301, R2311, R2312, R2313], [R2302, R2312, R2322, R2323], [R2303, R2313, R2323, R2333]]], [[[R0300, R0301, R0302, R0303], [R0301, R0311, R0312, R0313], [R0302, R0312, R0322, R0323], [R0303, R0313, R0323, R0333]], [[R1300, R1301, R1302, R1303], [R1301, R1311, R1312, R1313], [R1302, R1312, R1322, R1323], [R1303, R1313, R1323, R1333]], [[R2300, R2301, R2302, R2303], [R2301, R2311, R2312, R2313], [R2302, R2312, R2322, R2323], [R2303, R2313, R2323, R2333]], [[R3300, R3301, R3302, R3303], [R3301, R3311, R3312, R3313], [R3302, R3312, R3322, R3323], [R3303, R3313, R3323, R3333]]]]
    >>> print(sorted(list(glb_gridfcs_dict)))
    ['R0000', 'R0001', 'R0002', 'R0003', 'R0011', 'R0012', 'R0013', 'R0022', 'R0023', 'R0033', 'R0100', 'R0101', 'R0102', 'R0103', 'R0111', 'R0112', 'R0113', 'R0122', 'R0123', 'R0133', 'R0200', 'R0201', 'R0202', 'R0203', 'R0211', 'R0212', 'R0213', 'R0222', 'R0223', 'R0233', 'R0300', 'R0301', 'R0302', 'R0303', 'R0311', 'R0312', 'R0313', 'R0322', 'R0323', 'R0333', 'R1100', 'R1101', 'R1102', 'R1103', 'R1111', 'R1112', 'R1113', 'R1122', 'R1123', 'R1133', 'R1200', 'R1201', 'R1202', 'R1203', 'R1211', 'R1212', 'R1213', 'R1222', 'R1223', 'R1233', 'R1300', 'R1301', 'R1302', 'R1303', 'R1311', 'R1312', 'R1313', 'R1322', 'R1323', 'R1333', 'R2200', 'R2201', 'R2202', 'R2203', 'R2211', 'R2212', 'R2213', 'R2222', 'R2223', 'R2233', 'R2300', 'R2301', 'R2302', 'R2303', 'R2311', 'R2312', 'R2313', 'R2322', 'R2323', 'R2333', 'R3300', 'R3301', 'R3302', 'R3303', 'R3311', 'R3312', 'R3313', 'R3322', 'R3323', 'R3333']
    >>> print(glb_gridfcs_dict['R0001'].__dict__)
    {'name': 'R0001', 'group': 'EVOL', 'rank': 4, 'dimension': 4, 'gf_type': 'CCTK_REAL', 'f_infinity': 0.0, 'wavespeed': 1.0, 'is_basename': False}
    """
    # Step 1: Select the appropriate indexed expression declaration function from the factory.
    declare_func = IXP_RANK_FUNC_MAP.get(rank)
    if not declare_func:
        raise ValueError(f"Unsupported rank: {rank}. Rank must be between 1 and 4.")

    # Step 2: Declare the SymPy object. Symmetry is only applicable for rank > 1.
    declare_kwargs = kwargs.copy()
    if rank > 1:
        declare_kwargs["symmetry"] = symmetry
    else:
        # Remove symmetry kwarg if present, as declarerank1 does not accept it.
        declare_kwargs.pop("symmetry", None)

    indexed_obj = declare_func(basename, dimension=dimension, **declare_kwargs)

    # Step 3: Register each unique gridfunction component.
    # Symmetries (e.g., gDD01=gDD10) result in duplicate SymPy objects in the nested list.
    # We must extract the unique names to register each component only once.
    unique_gf_names = _flatten_and_unique_str(cast(List[Any], indexed_obj))

    # Manually adjust kwargs for registration of individual components.
    kwargs.update({"is_basename": False, "rank": rank})
    register_gridfunctions(unique_gf_names, dimension, **kwargs)

    # Step 4: Return the original (potentially nested) list structure of SymPy symbols.
    return cast(
        Union[
            List[sp.Symbol],
            List[List[sp.Symbol]],
            List[List[List[sp.Symbol]]],
            List[List[List[List[sp.Symbol]]]],
        ],
        indexed_obj,
    )


# Outside the classes above, only define functions in grid.py that are used by multiple infrastructures
#####################################


if __name__ == "__main__":
    import doctest
    import sys

    # Set a default infrastructure to avoid errors in certain doctests.
    # This is a common practice for making doctests runnable.
    par.set_parval_from_str("Infrastructure", "ETLegacy")
    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

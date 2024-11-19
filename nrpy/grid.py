"""
Provide registration functions for numerical grid parameters.

This module handles grid parameters and grid functions across various NRPy-supported infrastructures.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import Any, Dict, List, Optional, Sequence, Tuple, Union, cast

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


# Core GridFunction class.
class GridFunction:
    """The core class for grid functions."""

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
            # Sort case-insensitively to ensure consistent order, e.g., "RbarDD" doesn't appear before "cf".
            groups[group] = sorted(groups[group], key=lambda x: x.lower())

        # Pack the sorted lists into a tuple and return.
        return groups["EVOL"], groups["AUX"], groups["AUXEVOL"]

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
        parity_type: List[int] = []
        for name in list_of_gf_names:
            for gf in glb_gridfcs_dict.values():
                if gf.name == name:
                    parity_type__orig_len = len(parity_type)
                    if gf.rank == 0:
                        parity_type.append(0)
                    elif gf.rank == 1:
                        parity_type.append(int(gf.name[-1]) + 1)
                    elif gf.rank == 2:
                        idx0 = gf.name[-2]
                        idx1 = gf.name[-1]
                        parity_conditions: Dict[Tuple[str, str], int] = {}
                        if gf.dimension == 3:
                            parity_conditions = {
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
                        elif gf.dimension == 4:
                            parity_conditions = {
                                # scalar tt component
                                ("0", "0"): 0,
                                # vector ti components
                                ("0", "1"): 1,
                                ("1", "0"): 1,
                                ("0", "2"): 2,
                                ("2", "0"): 2,
                                ("0", "3"): 3,
                                ("3", "0"): 3,
                                # tensor ij components
                                ("1", "1"): 4,
                                ("1", "2"): 5,
                                ("2", "1"): 5,
                                ("1", "3"): 6,
                                ("3", "1"): 6,
                                ("2", "2"): 7,
                                ("2", "3"): 8,
                                ("3", "2"): 8,
                                ("3", "3"): 9,
                            }

                        parity_value: Union[int, None] = parity_conditions.get(
                            (idx0, idx1)
                        )
                        if parity_value is not None:
                            parity_type.append(parity_value)
                    if len(parity_type) == parity_type__orig_len:
                        raise ValueError(
                            f"Error: Could not figure out parity type for {gf.group} gridfunction: {gf.name}, {gf.name[-2]}, {gf.name[-1]}, {gf.rank}, {gf.dimension}"
                        )

        if len(parity_type) != len(list_of_gf_names):
            raise ValueError(
                "Error: For some reason the length of the parity types list did not match the length of the gf list."
            )
        return parity_type


class BHaHGridFunction(GridFunction):
    """The subclass for BlackHoles@Home grid functions."""

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
        if sync_gf_in_superB == None:
            if group == "EVOL" or group == "AUX":
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
        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

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

        def define_gfs(name: str, gfs: List[str]) -> str:
            """
            Generate string representation for a list of GridFunctions.

            :param name: The name of the group.
            :param gfs: The list of GridFunctions.
            :return: The resulting string representation.
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

        if evolved_variables_list:
            f_infinity_list = [
                str(glb_gridfcs_dict[var].f_infinity) for var in evolved_variables_list
            ]
            f_wavespeed_list = [
                str(glb_gridfcs_dict[var].wavespeed) for var in evolved_variables_list
            ]

            f_infinity_str = ", ".join(f_infinity_list)
            f_wavespeed_str = ", ".join(f_wavespeed_list)

            outstr += "\n\n// SET gridfunctions_f_infinity[i] = evolved gridfunction i's value in the limit r->infinity:\n"
            outstr += f"static const REAL gridfunctions_f_infinity[NUM_EVOL_GFS] = {{ {f_infinity_str} }};\n"

            outstr += "\n\n// SET gridfunctions_wavespeed[i] = evolved gridfunction i's characteristic wave speed:\n"
            outstr += f"static const REAL gridfunctions_wavespeed[NUM_EVOL_GFS] = {{ {f_wavespeed_str} }};\n"

        outstr += f"{define_gfs('AUX', auxiliary_variables_list)}"
        outstr += f"{define_gfs('AUXEVOL', auxevol_variables_list)}"

        return outstr


class ETLegacyGridFunction(GridFunction):
    """Subclass for basic (non-CarpetX) Einstein Toolkit grid functions."""

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
        _gf_array_name = gf_array_name

    def verify_gridfunction_group_is_valid(self) -> None:
        """
        Verify that the 'group' attribute of the GridFunction instance is one of the valid groups.

        The valid groups are:
        'EVOL': for evolved quantities (i.e., quantities stepped forward in time),
        'AUXEVOL': for auxiliary quantities needed at all points by evolved quantities,
        'AUX': for all other quantities needed at all gridpoints.

        :raises ValueError: If the 'group' attribute is not one of the valid groups.
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
        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

        access_str = f"{gf_name}[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"
        if use_GF_suffix:
            access_str = f"{gf_name}GF[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"
        return access_str

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
        i0 = f"i0+{i0_offset}".replace("+-", "-") if i0_offset != 0 else "i0"
        i1 = f"i1+{i1_offset}".replace("+-", "-") if i1_offset != 0 else "i1"
        i2 = f"i2+{i2_offset}".replace("+-", "-") if i2_offset != 0 else "i2"

        ret_string = f"{self.name}GF[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"
        # if use_GF_suffix defined AND set to False, remove GF suffix
        if "use_GF_suffix" in kwargs and not kwargs["use_GF_suffix"]:
            ret_string = f"{self.name}[CCTK_GFINDEX3D(cctkGH, {i0}, {i1}, {i2})]"
        if kwargs.get("enable_simd"):
            ret_string = f"ReadSIMD(&{ret_string})"

        return ret_string


class CarpetXGridFunction(GridFunction):
    """Class for CarpetX gridfunction operations."""

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
        _gf_array_name = gf_array_name
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
        Verify that the 'centering' attribute of the GridFunction instance is valid.

        A valid 'centering' is a 3-character string where each character is either 'C' or 'V'.
        'C' stands for cell-centered quantities, 'V' stands for vertex-centered quantities.

        :raises ValueError: If the 'centering' attribute is not valid.
        """
        if len(self.centering) != 3 or any(c not in ("C", "V") for c in self.centering):
            msg = (
                f"Unsupported gridfunction centering {self.centering}. Supported centerings are 3-character strings, each character being either:\n"
                '    "C": for cell-centered quantities,\n'
                '    "V": for vertex-centered quantities.'
            )
            raise ValueError(msg)

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
            index = "p.I"
            if i0_offset != 0:
                index += f" + {i0_offset}*p.DI[0]".replace("+ -", "- ")
            if i1_offset != 0:
                index += f" + {i1_offset}*p.DI[1]".replace("+ -", "- ")
            if i2_offset != 0:
                index += f" + {i2_offset}*p.DI[2]".replace("+ -", "- ")

        access_str = f"{gf_name}({index})"
        if use_GF_suffix:
            access_str = f"{gf_name}GF({index})"
        return access_str

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
        # if "reuse_index" in kwargs and not kwargs["reuse_index"]:
        index = "p.I"
        if i0_offset != 0:
            index += f" + {i0_offset}*p.DI[0]".replace("+ -", "- ")
        if i1_offset != 0:
            index += f" + {i1_offset}*p.DI[1]".replace("+ -", "- ")
        if i2_offset != 0:
            index += f" + {i2_offset}*p.DI[2]".replace("+ -", "- ")
        # else:
        #     if "index_name" in kwargs:
        #         index = kwargs["index_name"]
        #     else:
        #         #Error out
        #         exit(1)

        ret_string = f"{self.name}GF"
        # if use_GF_suffix defined AND set to True, add GF suffix
        if "use_GF_suffix" in kwargs and not kwargs["use_GF_suffix"]:
            ret_string = f"{self.name}"
        ret_string += f"({index})"

        if kwargs.get("enable_simd"):
            ret_string = f"ReadSIMD(&{ret_string})"

        return ret_string


# Contains a list of gridfunction objects.
glb_gridfcs_dict: Dict[
    str,
    Union[GridFunction, BHaHGridFunction, ETLegacyGridFunction, CarpetXGridFunction],
] = {}


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
            gf: Union[BHaHGridFunction, ETLegacyGridFunction, CarpetXGridFunction]
            kwargs_modify = kwargs.copy()
            for param in ["f_infinity", "wavespeed"]:
                if kwargs.get(param) and isinstance(kwargs.get(param), list):
                    # mypy: Once again bonks out after I've CONFIRMED kwargs.get(param) is not None and is a list!
                    kwargs_modify[param] = kwargs.get(param)[i]  # type: ignore
            if Infrastructure == "BHaH":
                gf = BHaHGridFunction(name, dimension=dimension, **kwargs_modify)
            elif Infrastructure == "ETLegacy":
                gf = ETLegacyGridFunction(name, dimension=dimension, **kwargs_modify)
            elif Infrastructure == "CarpetX":
                gf = CarpetXGridFunction(name, dimension=dimension, **kwargs_modify)
            else:
                raise ValueError(f"Infrastructure = {Infrastructure} unknown")

            glb_gridfcs_dict[name] = gf
        sympy_symbol_list += [sp.symbols(name, real=True)]

    # Step 3: Return the list of sympy symbols, which can have any number of elements >= 1
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
    register_gridfunctions(gf_list, dimension, **kwargs)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


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
    register_gridfunctions(gf_list, dimension, **kwargs)

    # Step 3: Return array of SymPy variables
    return IDX_OBJ_TMP


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

    def flatten_list(nested_list: List[Any]) -> List[Any]:
        """
        Flattens a nested list into a single list.

        :param nested_list: The list to be flattened.
        :return: A single list containing all the elements of the nested list.

        >>> flatten_list([[1, 2], [3, 4]])
        [1, 2, 3, 4]
        >>> flatten_list([[[1], [2]], [[3], [4]]])
        [1, 2, 3, 4]
        """
        if not isinstance(nested_list, list):
            return [nested_list]
        return [item for sublist in nested_list for item in flatten_list(sublist)]

    # Add type hint to IDX_OBJ_TMP, or some versions of mypy will assume its type from its first usage (List[sp.Symbol])
    IDX_OBJ_TMP: Union[
        List[sp.Symbol],
        List[List[sp.Symbol]],
        List[List[List[sp.Symbol]]],
        List[List[List[List[sp.Symbol]]]],
    ]

    # Determine the correct function to call based on rank
    if rank == 1:
        IDX_OBJ_TMP = cast(
            List[sp.Symbol], ixp.declarerank1(basename, dimension=dimension, **kwargs)
        )
    elif rank == 2:
        IDX_OBJ_TMP = cast(
            List[List[sp.Symbol]],
            ixp.declarerank2(
                basename, symmetry=symmetry, dimension=dimension, **kwargs
            ),
        )
    elif rank == 3:
        IDX_OBJ_TMP = cast(
            List[List[List[sp.Symbol]]],
            ixp.declarerank3(
                basename, symmetry=symmetry, dimension=dimension, **kwargs
            ),
        )
    elif rank == 4:
        IDX_OBJ_TMP = cast(
            List[List[List[List[sp.Symbol]]]],
            ixp.declarerank4(
                basename, symmetry=symmetry, dimension=dimension, **kwargs
            ),
        )
    else:
        raise ValueError(f"Unsupported rank: {rank}. Rank must be between 1 and 4.")

    # Flatten the list to register each gridfunction only once
    flat_list = [
        item
        for sublist in flatten_list(IDX_OBJ_TMP)
        for item in (sublist if isinstance(sublist, list) else [sublist])
    ]
    unique_gf_list = list(set(map(str, flat_list)))

    # Manually adjust kwargs for registration
    kwargs["is_basename"] = False
    kwargs["rank"] = rank
    register_gridfunctions(unique_gf_list, dimension, **kwargs)

    # Return the original nested list structure of SymPy symbols
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

"""
Provide registration functions for numerical grid parameters.

This module handles grid parameters and grid functions across various NRPy-supported infrastructures.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

from typing import List, Any, Union, Dict, Tuple, Optional, Sequence
from typing_extensions import Literal
import sympy as sp
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
_ = par.register_CodeParameter("int", __name__, "NUMGRIDS", 1, commondata=True)


# Core GridFunction class.
class GridFunction:
    """The core class for grid functions."""

    def __init__(
        self,
        name: str,
        group: str = "EVOL",
        rank: int = 0,
        dimension: int = 3,
        c_type_alias: str = "REAL",
        f_infinity: float = 0.0,
        wavespeed: float = 1.0,
        is_basename: bool = True,
    ) -> None:
        self.name: str = name
        self.group: str = group
        self.rank: int = rank
        self.dimension: int = dimension
        self.c_type_alias: str = c_type_alias
        self.f_infinity: float = f_infinity
        self.wavespeed: float = wavespeed
        self.is_basename: bool = is_basename

        if is_basename:
            self.verify_gridfunction_basename_is_valid(self.name)

    def read_gf_from_memory_Ccode_onept(
        self, i0_offset: int = 0, i1_offset: int = 0, i2_offset: int = 0, **kwargs: Any
    ) -> str:
        """Catch-all function for gridfunctions that haven't been set up correctly."""
        return f"Please define read_gf_from_memory_Ccode_onept() inside grid.py. Inputs: [class instance] {i0_offset} {i1_offset} {i2_offset} {kwargs}"

    def verify_gridfunction_basename_is_valid(self, name: str) -> None:
        """
        Validate the gridfunction's base name.

        Raises a ValueError if base name is zero-length or ends with an integer.
        Raises TypeError if the base name is not a string.

        Gridfunction base names must not end with integers. For instance, a rank-1
        gridfunction 'vecU1' has a valid base name 'vecU'. A scalar gridfunction
        named 'u2' would have an invalid base name 'u2'. This strict requirement
        facilitates quick identification of a gridfunction by its name.

        :param name: The name of the gridfunction.

        Raises
        ------
        ValueError
            If the name is of zero length or ends with an integer.
        TypeError
            If the name is not of type str.

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
            groups[group] = sorted(groups[group])

        # Pack the sorted lists into a tuple and return.
        return groups["EVOL"], groups["AUX"], groups["AUXEVOL"]

    @staticmethod
    def set_parity_types(list_of_gf_names: List[str]) -> List[int]:
        """
        Set the parity types for a given list of grid function names.

        :param list_of_gf_names: List of grid function names for which to set the parity types.
        :return: A list of integers representing the parity types for the grid functions.
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

        :raises ValueError: If 'gf_array_name' is not provided in kwargs.

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

        :raises ValueError: If 'gf_array_name' is not provided.

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

        :raises ValueError: If 'gf_array_name' is not provided.

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
        :raises ValueError: If 'gf_array_name' is not provided in kwargs or if kwargs are passed when they shouldn't be.

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

        :raises ValueError: If 'gf_array_name' is not provided.

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
        :raises ValueError: If 'gf_array_name' is not provided in kwargs or if kwargs are passed when they shouldn't be.

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
    Register grid functions with a specified name or list of names.

    :param names: A name or list of names of grid functions to be registered.
    :param kwargs: Additional keyword arguments.
    :return: A list of sympy symbols representing the registered grid functions.

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
    c_type_alias REAL
    f_infinity 0.0
    wavespeed 1.0
    is_basename True
    gf_array_name in_gfs
    >>> par.set_parval_from_str("Infrastructure", "ETLegacy")
    >>> gridfunc1, gridfunc2 = register_gridfunctions(['gridfunc1', 'gridfunc2'], f_infinity=[1.0, 4.0], is_basename=False)
    >>> for key, value in glb_gridfcs_dict["gridfunc1"].__dict__.items():
    ...     print(key, value)
    name gridfunc1
    group EVOL
    rank 0
    dimension 3
    c_type_alias CCTK_REAL
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
    Register gridfunctions for a single rank 1 variable.

    :param basename: Base name of the gridfunction to register.
    :param kwargs: Additional keyword arguments for registration.
    :return: List of SymPy symbols.

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

# nrpy/infrastructures/BHaH/CurviBoundaryConditions/BHaH_defines.py
"""
Defines core C data structures and constants for curvilinear boundary conditions.

This module generates C code for several key components:
- The `bc_struct`, which consolidates all boundary condition data.
- Helper structs (`innerpt_bc_struct`, `outerpt_bc_struct`) for inner and outer
  boundary points.
- Constant arrays that store the parity type for each grid function. Parity types
  are automatically inferred from grid function names.

All generated definitions are registered for inclusion in the main `BHaH_defines.h`
header file, making them globally available to the C code.

This process is documented in the NRPy tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.BHaH import BHaH_defines_h, griddata_commondata


# For example, if the gridfunction name ends with "01", then (based on the table in the
# NRPy Jupyter notebook corresponding to this Python module) the set_parity_types()
# function below will set the parity_type of that gridfunction to 5. We can be assured
# this is a rather robust algorithm, because gri.register_gridfunctions() in grid.py
# will throw an error if a gridfunction's base name ends in an integer. This strict
# syntax was added with the express purpose of making it easier to set parity types
# based solely on the gridfunction name.
#
# After each parity type is found, we store the parity type of each gridfunction to
# const int8_t arrays evol_gf_parity and aux_gf_parity, appended to the end of
# BHaH_defines.h.
def BHaH_defines_set_gridfunction_defines_with_parity_types(
    set_parity_on_aux: bool = False,
    set_parity_on_auxevol: bool = False,
    verbose: bool = True,
) -> str:
    """
    Set the grid function definitions with parity types and append them to the end of BHaH_defines.h.

    :param set_parity_on_aux: Flag to set parity on auxiliary variables. Default is False.
    :param set_parity_on_auxevol: Flag to set parity on auxevol variables. Default is False.
    :param verbose: Flag to control printing of details. Default is True.
    :return: A string containing the definitions for all grid functions with their parity types.
    """
    # First add human-readable gridfunction aliases (grid.py) to BHaH_defines dictionary.
    (
        evolved_variables_list,
        auxiliary_variables_list,
        auxevol_variables_list,
    ) = gri.BHaHGridFunction.gridfunction_lists()[0:3]

    # Next register outer_bc_type code parameter
    _ = par.CodeParameter(
        "char[50]", __name__, "outer_bc_type", "radiation", commondata=True
    )

    outstr = """
/* PARITY TYPES FOR EVOLVED (plus optional) GRIDFUNCTIONS.
 * SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */
"""
    if len(evolved_variables_list) > 0:
        evol_parity_type = gri.BHaHGridFunction.set_parity_types(evolved_variables_list)

        outstr += f"static const int8_t evol_gf_parity[{len(evolved_variables_list)}] = {{ {', '.join(map(str, evol_parity_type))} }};\n"
    if set_parity_on_aux:
        if len(auxiliary_variables_list) > 0:
            aux_parity_type = gri.BHaHGridFunction.set_parity_types(
                auxiliary_variables_list
            )
            outstr += f"static const int8_t aux_gf_parity[{len(auxiliary_variables_list)}] = {{ {', '.join(map(str, aux_parity_type))} }};\n"
    if set_parity_on_auxevol:
        if len(auxevol_variables_list) > 0:
            auxevol_parity_type = gri.BHaHGridFunction.set_parity_types(
                auxevol_variables_list
            )
            outstr += f"static const int8_t auxevol_gf_parity[{len(auxevol_variables_list)}] = {{ {', '.join(map(str, auxevol_parity_type))} }};\n"

    if verbose:
        for i, evolved_variable in enumerate(evolved_variables_list):
            print(
                f'Evolved gridfunction "{evolved_variable}" has parity type {evol_parity_type[i]}.'
            )
        if set_parity_on_aux:
            for i, auxiliary_variable in enumerate(auxiliary_variables_list):
                print(
                    f'Auxiliary gridfunction "{auxiliary_variable}" has parity type {aux_parity_type[i]}.'
                )
        if set_parity_on_auxevol:
            for i, auxevol_variable in enumerate(auxevol_variables_list):
                print(
                    f'AuxEvol gridfunction "{auxevol_variable}" has parity type {auxevol_parity_type[i]}.'
                )
    return outstr


def register_griddata_commondata() -> None:
    """
    Register the bcstruct's contribution to the griddata_struct and commondata_struct.

    This function registers the bcstruct, which contains all the data needed to perform
    boundary conditions in curvilinear coordinates, to the commondata structure.
    """
    griddata_commondata.register_griddata_commondata(
        __name__,
        "bc_struct bcstruct",
        "all data needed to apply boundary conditions in curvilinear coordinates",
    )


def register_BHaH_defines_h(
    set_parity_on_aux: bool = False,
    set_parity_on_auxevol: bool = False,
) -> None:
    """
    Register the bcstruct's contribution to the BHaH_defines.h file.

    This function registers the data structures needed for NRPy+ curvilinear boundary
    conditions in the BHaH_defines.h file, including structures for inner and outer
    boundary conditions and boundary loop bounds. It also sets the parity types for
    evolved, auxiliary, and auxiliary evolution grid functions.

    :param set_parity_on_aux: Flag to determine if parity should be set for auxiliary grid functions.
                              Default is False.
    :param set_parity_on_auxevol: Flag to determine if parity should be set for auxiliary evolution grid functions.
                                  Default is False.
    """
    CBC_BHd_str = r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

typedef struct __innerpt_bc_struct__ {
  int dstpt;  // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;  // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
  int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
} innerpt_bc_struct;

typedef struct __outerpt_bc_struct__ {
  short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct __bc_info_struct__ {
  int num_inner_boundary_points;  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3];  // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6];  // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct __bc_struct__ {
  innerpt_bc_struct *restrict inner_bc_array;  // information needed for updating each inner boundary point
  outerpt_bc_struct *restrict pure_outer_bc_array[NGHOSTS*3]; // information needed for updating each outer
  //                                                             boundary point
  bc_info_struct bc_info;  // stores number of inner and outer boundary points, needed for setting loop
  //                          bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""
    CBC_BHd_str += BHaH_defines_set_gridfunction_defines_with_parity_types(
        set_parity_on_aux=set_parity_on_aux,
        set_parity_on_auxevol=set_parity_on_auxevol,
        verbose=True,
    )
    BHaH_defines_h.register_BHaH_defines(
        __name__,
        CBC_BHd_str.replace(
            "*restrict",
            "*" if par.parval_from_str("parallelization") == "cuda" else "*restrict",
        ),
    )

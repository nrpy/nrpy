"""
Module providing base classes to support generating functions for setting up Curvilinear boundary conditions.

This is documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot** com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot** com
"""

# Step P1: Import needed NRPy+ core modules:
from typing import List, Tuple
import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import sympy.codegen.ast as sp_ast
import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.grid as gri  # NRPy+: Functions having to do with numerical grids
import nrpy.indexedexp as ixp  # NRPy+: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
import nrpy.finite_difference as fin  # NRPy+: Finite-difference module

# from nrpy.infrastructures.BHaH import griddata_commondata
# from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.validate_expressions.validate_expressions import check_zero
from nrpy.helpers.expr_tree import get_unique_expression_symbols

# Set unit-vector dot products (=parity) for each of the 10 parity condition types
def parity_conditions_symbolic_dot_products(
    CoordSystem: str, fp_type: str = "double"
) -> str:
    """
    Set unit-vector dot products (parity) for each of the 10 parity condition types.

    :param CoordSystem: The coordinate system for which to calculate the parity conditions.
    :param fp_type: Floating point type, e.g., "double".
    :return: C code string representing the unit vector dot products for the ten parity conditions.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    parity = ixp.zerorank1(dimension=10)
    UnitVectors_inner = ixp.zerorank2()
    xx0_inbounds, xx1_inbounds, xx2_inbounds = sp.symbols(
        "xx0_inbounds xx1_inbounds xx2_inbounds", real=True
    )
    for i in range(3):
        for j in range(3):
            UnitVectors_inner[i][j] = (
                rfm.UnitVectors[i][j]
                .subs(rfm.xx[0], xx0_inbounds)
                .subs(rfm.xx[1], xx1_inbounds)
                .subs(rfm.xx[2], xx2_inbounds)
            )
    # Type 0: scalar
    parity[0] = sp.sympify(1)
    # Type 1: i0-direction vector or one-form
    # Type 2: i1-direction vector or one-form
    # Type 3: i2-direction vector or one-form
    for i in range(3):
        for Type in range(1, 4):
            parity[Type] += (
                rfm.UnitVectors[Type - 1][i] * UnitVectors_inner[Type - 1][i]
            )
    # Type 4: i0i0-direction rank-2 tensor
    # parity[4] = parity[1]*parity[1]
    # Type 5: i0i1-direction rank-2 tensor
    # Type 6: i0i2-direction rank-2 tensor
    # Type 7: i1i1-direction rank-2 tensor
    # Type 8: i1i2-direction rank-2 tensor
    # Type 9: i2i2-direction rank-2 tensor
    count = 4
    for i in range(3):
        for j in range(i, 3):
            parity[count] = parity[i + 1] * parity[j + 1]
            count = count + 1

    lhs_strings = []
    for i in range(10):
        lhs_strings.append("REAL_parity_array[" + str(i) + "]")
    outstr = """/*
NRPy+ Curvilinear Boundary Conditions: Unit vector dot products for all
     ten parity conditions, in given coordinate system.
     Needed for automatically determining sign of tensor across coordinate boundary.
Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
*/
"""
    return outstr + ccg.c_codegen(parity, lhs_strings, fp_type=fp_type)


# For example, if the gridfunction name ends with "01", then (based on the table in the
# NRPy+ Jupyter notebook corresponding to this Python module) the set_parity_types()
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

    outstr = """
/* PARITY TYPES FOR EVOLVED (plus optional) GRIDFUNCTIONS.
 * SEE \"Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb\" FOR DEFINITIONS. */
"""
    if len(evolved_variables_list) > 0:
        evol_parity_type = gri.BHaHGridFunction.set_parity_types(evolved_variables_list)

        outstr += f"static const int8_t evol_gf_parity[{len(evolved_variables_list)}] = {{ {', '.join(map(str, evol_parity_type)) } }};\n"
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


# EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
    CoordSystem: str,
    fp_type: str = "double",
) -> str:
    """
    Map points between different coordinate systems.

    This function performs mapping between eigencoordinate form (x0,x1,x2) and its
    Cartesian equivalent (Cartx, Carty, Cartz), with the assumption of the same grid boundaries
    for both original and eigencoordinate systems. After mapping to Cartesian, it converts
    these coordinates back to an 'interior' point in eigencoordinate form. For cell-centered
    grids, this point aligns with a point on the numerical grid within round-off error.
    Finally, a check is done to ensure the conversion back to Cartesian matches the original
    values; a runtime error is thrown if not.

    :param CoordSystem: The coordinate system for mapping.
    :param fp_type: Floating point type, e.g., "double".

    :return: Body of the C code.

    :raises RuntimeError: If the conversion back to Cartesian coordinates does not match the original coordinates, indicating an error in the mapping process.
    """
    desc = """EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
  A coordinate system's "eigencoordinate" is the simplest member
  of its family; all spherical-like coordinate systems have
  Spherical as their eigencoordinate. The same is true for
  cylindrical-like (Cylindrical is eigencoordinate),
  Cartesian-like (Cartesian is the eigencoordinate), and
  SymTP-like (SymTP is the eigencoordinate) coordinates.

  For a given gridpoint (i0,i1,i2) and corresponding coordinate
  (x0,x1,x2), this function performs the dual mapping
  (x0,x1,x2) -> (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  Note that (x0,x1,x2) IS NOT ALWAYS equal to (x0,x1,x2)';
  For example consider in Spherical coordinates
  (x0,x1,x2)=(r,theta,phi)=(-0.1,pi/4,pi/4).
  This point will map to (x0,x1,x2)', in which x0>0,
  because the inversion r=sqrt(Cartx^2+Carty^2+Cartz^2)
  is always positive. In this case, (x0,x1,x2) is considered
  an *inner* boundary point, and on a cell-centered grid
  is guaranteed to map to a grid point in the grid interior;
  filling in this point requires copying data, and possibly
  multiplying by a +/- 1 if the data is from a gridfunction
  storing tensors/vectors.
"""
    cfunc_type = "static void"
    name = "EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
const int i0, const int i1, const int i2,
REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]"""
    type_literal = "" if fp_type == "double" else "f"
    body = r"""
  // This is a 3-step algorithm:
  // Step 1: (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  // Step 2: (Cartx,Carty,Cartz) -> (x0,x1,x2)'
  //         Find the interior eigencoordinate point
  //         (x0,x1,x2)' to which (Cartx,Carty,Cartz)
  //         maps, as well as the corresponding
  //         gridpoint integer index (i0,i1,i2). For
  //         cell-centered grids, (x0,x1,x2) will always
  //         overlap exactly (to roundoff error) a point
  //         on the numerical grid.
  // Step 3: Sanity check
  //         Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!
"""
    # Load up the EigenCoordinate corresponding to reference_metric::CoordSystem
    rfm_orig = refmetric.reference_metric[CoordSystem]
    rfm = refmetric.reference_metric[rfm_orig.EigenCoord]

    # Step 1: Output C code for the Eigen-Coordinate mapping from xx->Cartesian':
    body += (
        f"""
// Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
REAL xCart[3];  // where (x,y,z) is output
{{
// xx_to_Cart for EigenCoordinate {rfm.CoordSystem} (orig coord = {rfm_orig.CoordSystem}):
REAL xx0 = xx[0][i0];
REAL xx1 = xx[1][i1];
REAL xx2 = xx[2][i2];
"""
        + ccg.c_codegen(
            [rfm.xx_to_Cart[0], rfm.xx_to_Cart[1], rfm.xx_to_Cart[2]],
            ["xCart[0]", "xCart[1]", "xCart[2]"],
            fp_type=fp_type,
        )
        + "}\n"
    )
    body += r"""
REAL Cartx = xCart[0];
REAL Carty = xCart[1];
REAL Cartz = xCart[2];
"""

    # Step 2: Output C code for the Eigen-Coordinate mapping from Cartesian->xx':
    body += r"""
  // Step 2: Find the (i0_inbounds,i1_inbounds,i2_inbounds) corresponding to the above Cartesian coordinate.
  //   If (i0_inbounds,i1_inbounds,i2_inbounds) is in a ghost zone, then it must equal (i0,i1,i2), and
  //      the point is an outer boundary point.
  //   Otherwise (i0_inbounds,i1_inbounds,i2_inbounds) is in the grid interior, and data at (i0,i1,i2)
  //      must be replaced with data at (i0_inbounds,i1_inbounds,i2_inbounds), but multiplied by the
  //      appropriate parity condition (+/- 1).
  REAL Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
"""
    # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        raise RuntimeError(
            f"ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for "
            f"reference_metric::CoordSystem = {CoordSystem}. "
            "Boundary conditions in curvilinear coordinates REQUiRE this be set."
        )
    # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    body += f"  // Cart_to_xx for EigenCoordinate {rfm.CoordSystem} (orig coord = {rfm_orig.CoordSystem});\n"
    body += ccg.c_codegen(
        [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
        ["Cart_to_xx0_inbounds", "Cart_to_xx1_inbounds", "Cart_to_xx2_inbounds"],
        fp_type=fp_type,
    )
    body += rf"""
  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0{type_literal}/2.0{type_literal}))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((REAL)(0-NGHOSTS) + (1.0{type_literal}/2.0{type_literal}))*dxxi
  const REAL xxmin[3] = {{
    xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0{type_literal}/2.0{type_literal}))*dxx0,
    xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0{type_literal}/2.0{type_literal}))*dxx1,
    xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0{type_literal}/2.0{type_literal}))*dxx2 }};

  // Finally compute i{{0,1,2}}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0{type_literal}/2.0{type_literal})*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5{type_literal} );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0{type_literal}/2.0{type_literal})*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5{type_literal} );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0{type_literal}/2.0{type_literal})*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5{type_literal} );
"""

    # Restore reference_metric::CoordSystem back to the original CoordSystem

    # Step 3:
    body += f"""
  // Step 3: Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Cartx(x1(i1)),Cartx(x2(i2)))
  //         If not, error out!

  // Step 3.a: Compute {{x,y,z}}Cart_from_xx, as a
  //           function of i0,i1,i2
  REAL xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {{
    // xx_to_Cart for Coordinate {rfm_orig.CoordSystem}):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""
    body += ccg.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ["xCart_from_xx", "yCart_from_xx", "zCart_from_xx"],
        include_braces=False,
        fp_type=fp_type,
    )

    body += f"""  }}

  // Step 3.b: Compute {{x,y,z}}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  REAL xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {{
    // xx_to_Cart_inbounds for Coordinate {rfm_orig.CoordSystem}):
    REAL xx0 = xx[0][i0_inbounds];
    REAL xx1 = xx[1][i1_inbounds];
    REAL xx2 = xx[2][i2_inbounds];
"""
    body += ccg.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ["xCart_from_xx_inbounds", "yCart_from_xx_inbounds", "zCart_from_xx_inbounds"],
        include_braces=False,
        fp_type=fp_type,
    )

    body += r"""  }

  // Step 3.c: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
"""
    if fp_type == "float":
        body += r"""
#define EPS_REL 1e-6
"""
    else:
        body += r"""
#define EPS_REL 1e-8
"""
    body += (
        r"""
  const REAL norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (double)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {
    fprintf(stderr,"Error in """
        + rfm_orig.CoordSystem
        + r""" coordinate system: Inner boundary point does not map to grid interior point: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (double)xCart_from_xx,(double)yCart_from_xx,(double)zCart_from_xx,
            (double)xCart_from_xx_inbounds,(double)yCart_from_xx_inbounds,(double)zCart_from_xx_inbounds,
            xx[0][i0],xx[1][i1],xx[2][i2],
            xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
    exit(1);
  }

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
#undef EPS_REL
"""
    )
    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cf.full_function


# set_parity_for_inner_boundary_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__set_parity_for_inner_boundary_single_pt(
    CoordSystem: str, fp_type: str = "double"
) -> str:
    """
    Generate C code for setting the parity for inner boundary single point in a given coordinate system.

    :param CoordSystem: Coordinate system in which to set the parity
    :param fp_type: Floating point type, e.g., "double".
    :return: Full function C code as a string

    Doctest: FIXME
    """
    desc = """set_parity_for_inner_boundary_single_pt():
Given (x0,x1,x2)=(xx0,xx1,xx2) and
(x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
(see description of
EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
above for more details), here we compute the parity conditions
for all 10 tensor types supported by NRPy+."""
    cfunc_type = "static void"
    name = "set_parity_for_inner_boundary_single_pt"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
                const REAL xx0,const REAL xx1,const REAL xx2,  const REAL x0x1x2_inbounds[3], const int idx,
                innerpt_bc_struct *restrict innerpt_bc_arr"""
    if fp_type == "float":
        body = r"""
#define EPS_REL 1e-6
"""
    else:
        body = r"""
#define EPS_REL 1e-8
"""
    body += f"""
const REAL xx0_inbounds = x0x1x2_inbounds[0];
const REAL xx1_inbounds = x0x1x2_inbounds[1];
const REAL xx2_inbounds = x0x1x2_inbounds[2];

REAL REAL_parity_array[10];
{{
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+
{parity_conditions_symbolic_dot_products(CoordSystem, fp_type=fp_type)}
}}

// Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
for(int whichparity=0;whichparity<10;whichparity++) {{
    if( fabs(REAL_parity_array[whichparity]) < 1 - EPS_REL || fabs(REAL_parity_array[whichparity]) > 1 + EPS_REL ) {{
        fprintf(stderr,"Error at point (%e %e %e), which maps to (%e %e %e).\\n",
                xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
        fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\\n",
                REAL_parity_array[whichparity]);
        exit(1);
    }}
    for(int parity=0;parity<10;parity++) {{
        innerpt_bc_arr[idx].parity[parity] = 1;
        if(REAL_parity_array[parity] < 0) innerpt_bc_arr[idx].parity[parity] = -1;
    }}
}} // END for(int whichparity=0;whichparity<10;whichparity++)
#undef EPS_REL
"""
    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    )
    return cf.full_function


# bcstruct_set_up():
#      This function is documented in desc= and body= fields below.
class base_register_CFunction_bcstruct_set_up:
    """
    Base class for generating the function for setting up bcstruct.

    This function prescribes how inner and outer boundary points on the
    computational grid are filled, based on the given coordinate system (CoordSystem).

    :param CoordSystem: The coordinate system for which to set up boundary conditions.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:
        self.CoordSystem = CoordSystem
        self.fp_type = fp_type

        self.includes = [
            "BHaH_defines.h",
            "BHaH_function_prototypes.h",
        ]
        self.prefunc = (
            Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
                self.CoordSystem,
                fp_type=self.fp_type,
            )
        )
        self.prefunc += Cfunction__set_parity_for_inner_boundary_single_pt(
            self.CoordSystem, fp_type=self.fp_type
        )
        self.desc = r"""At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
  Recall that at each inner boundary point we must set innerpt_bc_struct:
    typedef struct __innerpt_bc_struct__ {
      int dstpt;  // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
       int srcpt;  // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
      int8_t parity[10];  // parity[10] is a calculation of dot products for the 10 independent parity types
    } innerpt_bc_struct;
  At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
    Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
        This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
        Cartesian coordinate (x,y,z), then finds the grid point
        (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
        corresponding to this Cartesian coordinate (x,y,z).
    If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an inner boundary point. We must set
        Set bcstruct->inner_bc_array for this point, which requires we specify
        both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
        conditions for this gridpoint. The latter is found & specified within the
        function set_parity_for_inner_boundary_single_pt().
    If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
        then we are at an outer boundary point. Take care of outer BCs in Step 2.
Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
  Recall that at each inner boundary point we must set outerpt_bc_struct:
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
  Outer boundary points are filled from the inside out, two faces at a time.
    E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
    including the ghostzones, with NGHOSTS=2.
    We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
    points in first, since they will in general (at least in the case of extrapolation
    outer BCs) depend on e.g., i0=2 and i0=3 points.
    Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
    since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
    Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
    depend on i0 min and max faces being filled. The remaining pattern goes like this:
    Upper x1 face: (i0={1,12},i1=12,i2={2,11})
    Lower x2 face: (i0={1,12},i1={1,12},i2=1)
    Upper x2 face: (i0={1,12},i1={1,12},i2=12)
    Lower x0 face: (i0=0,i1={1,12},i2={1,12})
    Upper x0 face: (i0=13,i1={1,12},i2={1,12})
    Lower x1 face: (i0={0,13},i1=0,i2={2,11})
    Upper x1 face: (i0={0,13},i1=13,i2={2,11})
    Lower x2 face: (i0={0,13},i1={0,13},i2=0)
    Upper x2 face: (i0={0,13},i1={0,13},i2=13)
  Note that we allocate a outerpt_bc_struct at *all* boundary points,
    regardless of whether the point is an outer or inner point. However
    the struct is set only at outer boundary points. This is slightly
    wasteful, but only in memory, not in CPU."""
        self.cfunc_type = "void"
        self.name = "bcstruct_set_up"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3], bc_struct *restrict bcstruct"


###############################
## apply_bcs_inner_only(): Apply inner boundary conditions.
##  Function is documented below in desc= and body=.
class base_register_CFunction_apply_bcs_inner_only:
    """
    Register C function for filling inner boundary points on the computational grid.
    Filling is prescribed by bcstruct.
    """

    def __init__(self) -> None:
        self.includes = ["BHaH_defines.h"]
        self.desc = r"""
    Apply BCs to inner boundary points only,
    using data stored in bcstruct->inner_bc_array.
    These structs are set in bcstruct_set_up().
    Inner boundary points map to either the grid
    interior ("pure inner") or to pure outer
    boundary points ("inner maps to outer").
    """
        self.cfunc_type = "void"
        self.name = "apply_bcs_inner_only"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"
        self.body = ""


###############################
## apply_bcs_outerextrap_and_inner(): Apply extrapolation outer boundary conditions.
##  Function is documented below in desc= and body=.
class base_register_CFunction_apply_bcs_outerextrap_and_inner:
    """Register C function for filling boundary points with extrapolation and prescribed bcstruct."""

    def __init__(self) -> None:
        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.desc = r"""#Suppose the outer boundary point is at the i0=max(i0) face. Then we fit known data at i0-3, i0-2, and i0-1
#  to the unique quadratic polynomial that passes through those points, and fill the data at
#  i0 with the value implied from the polynomial.
#As derived in nrpytutorial's Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
#  the coefficients must be f_{i0} = f_{i0-3} - 3 f_{i0-2} + 3 f_{i0-1}.
#  To check these coefficients are correct, consider
#  * f(x0 = constant. Then f_{i0} = f_{i0-3} <- CHECK!
#  * f(x) = x. WOLOG suppose x0=0. Then f_{i0} = (-3dx) - 3(-2dx) + 3(-dx) = + dx(-3+6-3) = 0 <- CHECK!
#  * f(x) = x^2. WOLOG suppose x0=0. Then f_{i0} = (-3dx)^2 - 3(-2dx)^2 + 3(-dx)^2 = + dx^2(9-12+3) = 0 <- CHECK!"""
        self.cfunc_type = "void"
        self.name = "apply_bcs_outerextrap_and_inner"
        self.params = "const commondata_struct *restrict commondata, const params_struct *restrict params, const bc_struct *restrict bcstruct, REAL *restrict gfs"


###############################
## RADIATION (NewRad-like) BOUNDARY CONDITIONS.
##  Functions are fully documented in nrpytutorial's
##   Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
##   as well as below, in desc= and body=.
# r_and_partial_xi_partial_r_derivs(): Compute r(x0,x1,x2) and dx^i / dr
class setup_Cfunction_r_and_partial_xi_partial_r_derivs:
    """
    Generate C code to compute the radial coordinate r(x0, x1, x2) and its derivatives.

    Compute the radial coordinate r(x0, x1, x2) and its partial derivatives
    partial x^i / partial r for a given coordinate system.

    :param CoordSystem: The coordinate system for which to compute r and its derivatives.
    :param fp_type: Floating point type, e.g., "double".
    :return: A string containing the generated C code for the function.
    """
    
    def __init__(
        self,
        CoordSystem: str, fp_type: str = "double"
    ) -> str:
        self.CoordSystem=CoordSystem
        self.fp_type=fp_type
        self.CFunction = None
        self.include_CodeParameters_h = True
        self.cfunc_decorators=""
        
        self.desc = "Compute r(xx0,xx1,xx2) and partial_r x^i."
        self.cfunc_type = "static inline void"
        self.name = "r_and_partial_xi_partial_r_derivs"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
        const REAL xx0,const REAL xx1,const REAL xx2,    REAL *r,
        REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r"""
        rfm = refmetric.reference_metric[CoordSystem]
        # sp.simplify(expr) is too slow here for SinhCylindrical
        self.expr_list = [
                rfm.xxSph[0],
                rfm.Jac_dUrfm_dDSphUD[0][0],
                rfm.Jac_dUrfm_dDSphUD[1][0],
                rfm.Jac_dUrfm_dDSphUD[2][0],
            ]
        self.unique_symbols = []
        for expr in self.expr_list:
            sub_list = get_unique_expression_symbols(expr, exclude=[f'xx{i}' for i in range(3)])
            self.unique_symbols += sub_list
        self.unique_symbols = sorted(list(set(self.unique_symbols)))
        self.body = ccg.c_codegen(
            self.expr_list,
            [
                "*r",
                "*partial_x0_partial_r",
                "*partial_x1_partial_r",
                "*partial_x2_partial_r",
            ],
            verbose=False,
            include_braces=False,
            fp_type=fp_type,
        )
        self.generate_CFunction()

    def generate_CFunction(self):
        self.CFunction = cfc.CFunction(
            subdirectory=self.CoordSystem,
            includes=[],
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=self.include_CodeParameters_h,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators
    )


# partial_r f term: generate finite-difference coefficients
#   for \partial_i f with arbitrary upwinding:
def get_arb_offset_FD_coeffs_indices(
    FDORDER: int, offset: int, deriv: int
) -> Tuple[List[float], List[int]]:
    """
    Generate finite-difference coefficients for partial derivatives with arbitrary upwinding.

    :param FDORDER: Order of the finite difference
    :param offset: Offset for upwinding
    :param deriv: Order of the derivative (e.g., 1 for 1st derivative)

    :return: A tuple containing the list of coefficients and the list of indices

    Doctest:
    >>> get_arb_offset_FD_coeffs_indices(3, 0, 1)
    ([1/24, -9/8, 9/8, -1/24], [-1, 0, 1, 2])
    """
    # deriv = 1 <-- 1st derivative
    Minv = fin.setup_FD_matrix__return_inverse(FDORDER + 1, offset)
    indices = []
    coeffs = []
    for i in range(FDORDER + 1):
        indices.append(i - int(FDORDER / 2) + offset)
        coeffs.append(Minv[i, deriv])
    return coeffs, indices


# partial_r f term: FD1_arbitrary_upwind(): C function to evaluate
#   partial_i f with arbitrary upwinding
class setup_Cfunction_FD1_arbitrary_upwind:
    """
    Set up the C function for computing the 1st derivative finite-difference.

    Supports arbitrary upwind for a given direction and order.

    :param dirn: Direction in which to compute the derivative.
    :param radiation_BC_fd_order: Finite difference order for radiation boundary condition.
                                  If -1, will use default finite difference order.
    :param fp_type: Floating point type, e.g., "double".
    :return: The full C function as a string.
    """
    def __init__(
        self,
        dirn: int,
        radiation_BC_fd_order: int = -1,
        fp_type: str = "double",
    ) -> None:
        self.dirn=dirn
        self.radiation_BC_fd_order=radiation_BC_fd_order
        self.fp_type=fp_type
        
        import sympy.codegen.ast as sp_ast
        self.default_FDORDER = par.parval_from_str("fd_order")
        if radiation_BC_fd_order == -1:
            radiation_BC_fd_order = self.default_FDORDER

        par.set_parval_from_str("fd_order", radiation_BC_fd_order)

        self.include_CodeParameters_h=True
        self.includes: List[str] = []
        self.desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
        self.cfunc_type = "static inline REAL"
        self.cfunc_decorators = ""
        self.name = f"FD1_arbitrary_upwind_x{dirn}_dirn"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const REAL *restrict gf,  const int i0,const int i1,const int i2, const int offset"""
        self.body = "switch(offset) {\n"

        tmp_list: List[int] = []
        fp_ccg_type = ccg.fp_type_to_sympy_type[fp_type]
        sp_type_alias = {sp_ast.real: fp_ccg_type}
        for offset in range(
            0, int(radiation_BC_fd_order // 2) + 1
        ):  # Use // for integer division
            tmp_list.append(offset)
            if offset > 0:
                tmp_list.append(-offset)

        for offset in tmp_list:
            self.body += f"case {offset}:\n"
            self.body += "  return ("
            coeffs, indices = get_arb_offset_FD_coeffs_indices(
                radiation_BC_fd_order, offset, 1
            )

            for i, coeff in enumerate(coeffs):
                if coeff == 0:
                    continue
                offset_str: str = str(indices[i])
                if i > 0:
                    self.body += "          "
                if offset_str == "0":
                    self.body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2)]\n"
                else:
                    if dirn == 0:
                        self.body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0+{offset_str},i1,i2)]\n"
                    elif dirn == 1:
                        self.body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1+{offset_str},i2)]\n"
                    elif dirn == 2:
                        self.body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2+{offset_str})]\n"

            self.body = self.body[:-1].replace("+-", "-") + f") * invdxx{dirn};\n"

        self.body += """}
    return 0.0 / 0.0;  // poison output if offset computed incorrectly
    """
        par.set_parval_from_str("fd_order", self.default_FDORDER)
        self.generate_CFunction()
        
    def generate_CFunction(self) -> None:
        "Generate CFunction from class parameters."        
        self.CFunction = cfc.CFunction(
            subdirectory="one_subdirectory_down",
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=self.include_CodeParameters_h,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators,
        )


# partial_r f term: Numerically evaluate partial_r f,
#   calling functions defined above.
class setup_Cfunction_compute_partial_r_f:
    """
    Set up a C function for computing the partial derivative of f with respect to r.

    :param CoordSystem: Coordinate system to be used for the computation
    :param radiation_BC_fd_order: Order of finite difference for radiation boundary conditions, default is -1
    :return: A C function for computing the partial derivative
    """    
    def __init__(
        self,
        CoordSystem: str, radiation_BC_fd_order: int = -1
    ) -> str:
        self.CoordSystem=CoordSystem
        self.radiation_BC_fd_order=radiation_BC_fd_order
        self.include_CodeParameters_h=True
        self.includes = []
        self.cfunc_decorators=""
        self.CFunction= None
        
        self.desc = "Compute \\partial_r f"
        self.cfunc_type = "static inline REAL"
        self.name = "compute_partial_r_f"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
REAL *restrict xx[3], const REAL *restrict gfs,
const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
const int FACEi0,const int FACEi1,const int FACEi2,
const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r"""

        self.default_FDORDER = par.parval_from_str("fd_order")
        if radiation_BC_fd_order == -1:
            radiation_BC_fd_order = self.default_FDORDER

        self.FD1_stencil_radius = int(radiation_BC_fd_order / 2)
        self.body = ""
        self.tmp_definitions = f"""///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = {self.FD1_stencil_radius}
  const int FD1_stencil_radius = {self.FD1_stencil_radius};

  const int ntot = Nxx_plus_2NGHOSTS0*Nxx_plus_2NGHOSTS1*Nxx_plus_2NGHOSTS2;
"""
        self.algorithm_header = f"""  
  ///////////////////////////////////////////////////////////
  // Next we'll compute partial_xi f, using a maximally-centered stencil.
  //   The {{i0,i1,i2}}_offset parameters set the offset of the maximally-centered
  //   stencil, such that an offset=0 implies a centered stencil.

  // CHECK: Nxx_plus_2NGHOSTS0=10; FD1_stencil_radius=2. Then Nxx_plus_2NGHOSTS0-FD1_stencil_radius-1 = 7
  //  if dest_i0 = 9, we get i0_offset=7-9=-2, so the (4th order) deriv
  //  stencil is: -4,-3,-2,-1,0

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 1, we get i0_offset = FD1_stencil_radius-dest_i0 = 1,
  //  so the (4th order) deriv stencil is: -1,0,1,2,3

  // CHECK: if FD1_stencil_radius=2 and dest_i0 = 0, we get i0_offset = FD1_stencil_radius-1 = 2,
  //  so the (4th order) deriv stencil is: 0,1,2,3,4
    """
        self.generate_CFunction()
        
    def regenerate_body(self):
        rfm = refmetric.reference_metric[self.CoordSystem]
        self.body = ""
        for i in range(3):
            si = str(i)
            if check_zero(rfm.Jac_dUrfm_dDSphUD[i][0]):
                self.body += f"  const REAL partial_x{si}_f=0.0;\n"
            else:
                self.body += (
                    f"  int i{si}_offset = FACEi{si};  // Shift stencil away from the face we're updating.\n"
                    f"  // Next adjust i{si}_offset so that FD stencil never goes out of bounds.\n"
                    f"  if(dest_i{si} < FD1_stencil_radius) i{si}_offset = FD1_stencil_radius-dest_i{si};\n"
                    f"  else if(dest_i{si} > (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1)) i{si}_offset = (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1) - dest_i{si};\n"
                    f"  const REAL partial_x{si}_f=FD1_arbitrary_upwind_x{si}_dirn(commondata, params,&gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i{si}_offset);\n"
                )
        self.body += "  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;\n"

    def generate_CFunction(self) -> None:
        "Generate CFunction from class parameters."
        self.regenerate_body()
        self.body = self.tmp_definitions + self.algorithm_header + self.body
        self.CFunction = cfc.CFunction(
            subdirectory="one_subdirectory_down",
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=self.include_CodeParameters_h,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators,
        )


# radiation_bcs(): Put it all together, for a single outer boundary point.
class setup_Cfunction_radiation_bcs:
    """
    Generate C code to apply radiation boundary conditions in a given coordinate system.

    :param CoordSystem: The coordinate system to use.
    :param radiation_BC_fd_order: Finite differencing order to use. Default is -1.
    :param fp_type: Floating point type, e.g., "double".
    :return: A string containing the generated C code for the function.
    """
    def __init__(
        self,
        CoordSystem: str,
        radiation_BC_fd_order: int = -1,
        fp_type: str = "double",
    ) -> str:
        self.CoordSystem=CoordSystem
        self.radiation_BC_fd_order=radiation_BC_fd_order
        self.fp_type=fp_type
        
        self.cfunc_decorators=""
        self.include_CodeParameters_h=True
        self.includes: List[str] = []
        self.prefunc = ""
        self.rfm = refmetric.reference_metric[self.CoordSystem]
        
        # These functions can be replaced to generate different prefunc strings
        self.upwind_setup_func = setup_Cfunction_FD1_arbitrary_upwind
        self.r_and_partial_xi_partial_r_derivs_prefunc_setup_func = setup_Cfunction_r_and_partial_xi_partial_r_derivs
        self.compute_partial_r_f_setup_func = setup_Cfunction_compute_partial_r_f
        
        # Initialize prefunc strings
        self.upwind_prefunc=""
        self.r_and_partial_xi_partial_r_derivs_prefunc=""
        self.compute_partial_r_f_prefunc=""
        
        self.desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
        self.cfunc_type = "static inline REAL"
        self.name = "radiation_bcs"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
        const bc_struct *restrict bcstruct,REAL *restrict xx[3],
        const REAL *restrict gfs, REAL *restrict gfs_rhss,
        const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
        const int dest_i0,const int dest_i1,const int dest_i2,
        const short FACEi0,const short FACEi1,const short FACEi2"""
        self.variable_defs = r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
REAL r, partial_x0_partial_r,partial_x1_partial_r,partial_x2_partial_r;
REAL r_int, partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int;
"""
        self.function_calls = """
r_and_partial_xi_partial_r_derivs(commondata, params,xx[0][dest_i0],xx[1][dest_i1],xx[2][dest_i2],
                                  &r, &partial_x0_partial_r, &partial_x1_partial_r,  &partial_x2_partial_r);
r_and_partial_xi_partial_r_derivs(commondata, params, xx[0][dest_i0_int], xx[1][dest_i1_int], xx[2][dest_i2_int],
                                  &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
const REAL partial_r_f     = compute_partial_r_f(commondata, params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
const REAL partial_r_f_int = compute_partial_r_f(commondata, params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);
"""
        self.algorithm_body = """
const int idx3 = IDX3(dest_i0,dest_i1,dest_i2);
const int idx3_int = IDX3(dest_i0_int,dest_i1_int,dest_i2_int);

const REAL partial_t_f_int = gfs_rhss[IDX4pt(which_gf, idx3_int)];

const REAL c = gf_wavespeed;
const REAL f_infinity = gf_f_infinity;
const REAL f     = gfs[IDX4pt(which_gf, idx3)];
const REAL f_int = gfs[IDX4pt(which_gf, idx3_int)];
const REAL partial_t_f_int_outgoing_wave = -c * (partial_r_f_int + (f_int - f_infinity) / r_int);

const REAL k = r_int*r_int*r_int * (partial_t_f_int - partial_t_f_int_outgoing_wave);

const REAL rinv = 1.0 / r;
const REAL partial_t_f_outgoing_wave = -c * (partial_r_f + (f - f_infinity) * rinv);

return partial_t_f_outgoing_wave + k * rinv*rinv*rinv;
"""
        self.body = ""
        self.generate_CFunction()
    def generate_upwind_prefunc(self) -> None:
        self.upwind_prefunc=""
        
        for i in range(3):
            # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
            if not check_zero(self.rfm.Jac_dUrfm_dDSphUD[i][0]):
                self.upwind_prefunc += self.upwind_setup_func(
                    dirn=i,
                    radiation_BC_fd_order=self.radiation_BC_fd_order,
                    fp_type=self.fp_type,
                ).CFunction.full_function
    
    def generate_r_and_partial_xi_partial_r_derivs_prefunc(self) -> None:
        self.r_and_partial_xi_partial_r_derivs_prefunc=""        
        self.r_and_partial_xi_partial_r_derivs_prefunc += self.r_and_partial_xi_partial_r_derivs_prefunc_setup_func(
            CoordSystem=self.CoordSystem,
            fp_type=self.fp_type,
        ).CFunction.full_function
        
    def generate_compute_partial_r_f_prefunc(self):
        self.compute_partial_r_f_prefunc=""
        self.compute_partial_r_f_prefunc += self.compute_partial_r_f_setup_func(
            CoordSystem=self.CoordSystem, radiation_BC_fd_order=self.radiation_BC_fd_order
        ).CFunction.full_function
        
    def generate_CFunction(self):
        self.generate_upwind_prefunc()
        self.generate_r_and_partial_xi_partial_r_derivs_prefunc()
        self.generate_compute_partial_r_f_prefunc()
        
        self.prefunc = self.upwind_prefunc
        self.prefunc += self.r_and_partial_xi_partial_r_derivs_prefunc
        self.prefunc += self.compute_partial_r_f_prefunc
        
        self.body = self.variable_defs
        self.body += self.function_calls
        self.body +=self.algorithm_body
        
        self.CFunction = cfc.CFunction(
            subdirectory=self.CoordSystem,
            includes=self.includes,
            prefunc=self.prefunc,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=self.include_CodeParameters_h,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators,
        )


# apply_bcs_outerradiation_and_inner():
#   Apply radiation BCs at outer boundary points, and
#   inner boundary conditions at inner boundary points.
class base_register_CFunction_apply_bcs_outerradiation_and_inner:
    """
    Register a C function to apply boundary conditions to both pure outer and inner boundary points.

    :param CoordSystem: The coordinate system to use.
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(
        self,
        CoordSystem: str,
        radiation_BC_fd_order: int = 2,
        fp_type: str = "double",
    ) -> None:
        self.CoordSystem = CoordSystem
        self.radiation_BC_fd_order = radiation_BC_fd_order
        self.fp_type = fp_type

        self.includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
        self.prefunc = setup_Cfunction_radiation_bcs(
            CoordSystem=CoordSystem,
            radiation_BC_fd_order=radiation_BC_fd_order,
            fp_type=fp_type,
        ).CFunction.full_function
        self.desc = """This function is responsible for applying boundary conditions (BCs) to both pure outer and inner
boundary points. In the first step, it parallelizes the task using OpenMP and starts by applying BCs to
the outer boundary points layer-by-layer, prioritizing the faces in the order x0, x1, x2. The second step
applies BCs to the inner boundary points, which may map either to the grid interior or to the outer boundary.
"""
        self.cfunc_type = "void"
        self.name = "apply_bcs_outerradiation_and_inner"
        self.params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const bc_struct *restrict bcstruct, REAL *restrict xx[3],
    const REAL custom_wavespeed[NUM_EVOL_GFS],
    const REAL custom_f_infinity[NUM_EVOL_GFS],
    REAL *restrict gfs, REAL *restrict rhs_gfs
"""


class base_CurviBoundaryConditions_register_C_functions:
    """
    Base class to generate functions responsible for handling boundary conditions.

    :param list_of_CoordSystems: List of coordinate systems to use.
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param set_parity_on_aux: If True, set parity on auxiliary grid functions.
    :param set_parity_on_auxevol: If True, set parity on auxiliary evolution grid functions.
    :param fp_type: Floating point type, e.g., "double".
    """

    def __init__(
        self,
        list_of_CoordSystems: List[str],
        radiation_BC_fd_order: int = 2,
        set_parity_on_aux: bool = False,
        set_parity_on_auxevol: bool = False,
        fp_type: str = "double",
    ) -> None:
        self.list_of_CoordSystems = list_of_CoordSystems
        self.radiation_BC_fd_order = radiation_BC_fd_order
        self.set_parity_on_aux = set_parity_on_aux
        self.set_parity_on_auxevol = set_parity_on_auxevol
        self.fp_type = fp_type

        # Register bcstruct's contribution to BHaH_defines.h:
        self.CBC_BHd_str = r"""
// NRPy+ Curvilinear Boundary Conditions: Core data structures
// Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

typedef struct __innerpt_bc_struct__ {
int dstpt;  // dstpt is the 3D grid index IDX3S(i0,i1,i2) of the inner boundary point (i0,i1,i2)
int srcpt;  // srcpt is the 3D grid index (a la IDX3S) to which the inner boundary point maps
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
        self.CBC_BHd_str += BHaH_defines_set_gridfunction_defines_with_parity_types(
            set_parity_on_aux=self.set_parity_on_aux,
            set_parity_on_auxevol=self.set_parity_on_auxevol,
            verbose=True,
        )


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

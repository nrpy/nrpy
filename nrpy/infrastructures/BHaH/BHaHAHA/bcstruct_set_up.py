"""
Module providing functions for setting up Curvilinear boundary conditions.

This is documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
"""

# Step P1: Import needed NRPy core modules:
from typing import Dict, List, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.indexedexp as ixp  # NRPy: Symbolic indexed expression (e.g., tensors, vectors, etc.) support
import nrpy.params as par
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.infrastructures import BHaH


# Set unit-vector dot products (=parity) for each of the 28 parity condition types
def parity_conditions_symbolic_dot_products(
    CoordSystem: str,
) -> str:
    """
    Set unit-vector dot products (parity) for each of the 28 parity condition types.

    :param CoordSystem: The coordinate system for which to calculate the parity conditions.
    :return: C code string representing the unit vector dot products for the 28 parity conditions.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    if rfm.CoordSystem.startswith("GeneralRFM"):
        outstr = """/*
NRPy Curvilinear Boundary Conditions: Unit vector dot products for all
     twenty-eight parity conditions, in given coordinate system.
     Needed for automatically determining sign of tensor across coordinate boundary.
Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
*/
"""
        return (
            outstr
            + "\n".join([f"REAL_parity_array[{i}] = 1.0;" for i in range(28)])
            + "\n"
        )
    parity = ixp.zerorank1(dimension=28)
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

    # Handling special case of partial_k gamma_{ij}: 6*3 parities = 18
    # Types 10-12: partial_{0-2} gamma_{00}
    # Types 13-15: partial_{0-2} gamma_{01}
    # Types 16-18: partial_{0-2} gamma_{02}
    # Types 19-21: partial_{0-2} gamma_{11}
    # Types 22-24: partial_{0-2} gamma_{12}
    # Types 25-27: partial_{0-2} gamma_{22}
    for k in range(3):
        for i in range(3):
            for j in range(i, 3):
                parity[count] = parity[i + 1] * parity[j + 1] * parity[k + 1]
                count = count + 1

    lhs_strings = []
    for i in range(28):
        lhs_strings.append("REAL_parity_array[" + str(i) + "]")
    outstr = """/*
NRPy Curvilinear Boundary Conditions: Unit vector dot products for all
     twenty-eight parity conditions, in given coordinate system.
     Needed for automatically determining the sign of tensors across coordinate boundaries.
Documented in: Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb
*/
"""
    return outstr + ccg.c_codegen(parity, lhs_strings)


# For example, if the gridfunction name ends with "01", then (based on the table in the
# NRPy Jupyter notebook corresponding to this Python module) the set_parity_types()
# function below will set the parity_type of that gridfunction to 5. We can be assured
# this is a rather robust algorithm because grid.register_gridfunctions() in grid.py
# will throw an error if a gridfunction's base name ends in an integer. This strict
# syntax was added with the express purpose of making it easier to set parity types
# based solely on the gridfunction name.


# After each parity type is found, we store the parity type of each gridfunction to
# const int8_t arrays evol_gf_parity and aux_gf_parity, appended to the end of
# BHaH_defines.h.
def BHaH_defines_set_gridfunction_defines_with_parity_types(
    grid_name: str,
    list_of_gf_names_ranks: List[Tuple[str, int]],
    verbose: bool = True,
) -> str:
    """
    Set the grid function definitions with parity types and append them to the end of BHaH_defines.h.

    :param grid_name: Name of grid, typically external_input, interp_src, or evol
    :param list_of_gf_names_ranks: List of tuples containing gridfunction names and their ranks.
    :param verbose: Flag to control printing of details. Default is True.
    :return: A string containing the definitions for all grid functions with their parity types.

    :raises ValueError: If the parity type for a given gridfunction name and rank cannot be determined.
    """
    outstr = """
/* PARITY TYPES FOR SRC GRID GRIDFUNCTIONS. */
"""
    parity_type: List[int] = []
    for name, rank in list_of_gf_names_ranks:
        parity_type__orig_len = len(parity_type)
        if rank == 0:
            parity_type.append(0)
        elif rank == 1:
            parity_type.append(int(name[-1]) + 1)
        elif rank == 2:
            idx0 = int(name[-2])
            idx1 = int(name[-1])
            parity_conditions2: Dict[Tuple[int, int], int] = {}
            count = 4
            for i in range(3):
                for j in range(i, 3):
                    parity_conditions2[(i, j)] = count
                    count = count + 1
            parity_value = parity_conditions2.get((idx0, idx1))
            if parity_value is not None:
                parity_type.append(parity_value)
        elif rank == 3:
            idx0 = int(name[-3])
            idx1 = int(name[-2])
            idx2 = int(name[-1])
            # Handling special case of partial_k gamma_{ij}: 6*3 parities = 18
            # Types 10-12: partial_{0-2} gamma_{00}
            # Types 13-15: partial_{0-2} gamma_{01}
            # Types 16-18: partial_{0-2} gamma_{02}
            # Types 19-21: partial_{0-2} gamma_{11}
            # Types 22-24: partial_{0-2} gamma_{12}
            # Types 25-27: partial_{0-2} gamma_{22}
            parity_conditions3: Dict[Tuple[int, int, int], int] = {}
            count = 10
            for k in range(3):
                for i in range(3):
                    for j in range(i, 3):
                        parity_conditions3[(k, i, j)] = count
                        count = count + 1
            parity_value = parity_conditions3.get((idx0, idx1, idx2))
            if parity_value is not None:
                parity_type.append(parity_value)

        if len(parity_type) == parity_type__orig_len:
            raise ValueError(
                f"Error: Could not figure out parity type for {grid_name} gridfunction: {name}, with rank = {rank}"
            )

    outstr += f"static const int8_t {grid_name}_gf_parity[{len(list_of_gf_names_ranks)}] = {{ {', '.join(map(str, parity_type))} }};\n"

    if verbose:
        for i, (name, rank) in enumerate(list_of_gf_names_ranks):
            print(
                f'Grid function "{name}" with rank {rank} has parity type {parity_type[i]}.'
            )
    return outstr


# EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
    CoordSystem: str,
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

    :return: Body of the C code.

    :raises RuntimeError: If the conversion back to Cartesian coordinates does not match the original coordinates, indicating an error in the mapping process.
    """
    desc = """EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
An "eigencoordinate" represents the simplest member of a given coordinate family,
often serving as the prototype for the entire family of coordinate systems.
For instance, all spherical-like coordinate systems share Spherical as their eigencoordinate.
Similarly, for cylindrical-like systems, Cylindrical is the eigencoordinate;
for Cartesian-like systems, Cartesian serves as the eigencoordinate;
and for SymTP-like systems (symmetric-TwoPunctures, a prolate-spheroidal coordinate system),
SymTP is the eigencoordinate.

This function performs a dual mapping for a given grid point (i0, i1, i2) and its corresponding
coordinate (x0, x1, x2). The transformation is carried out as:
(x0, x1, x2) -> (Cartx, Carty, Cartz) -> (x0', x1', x2').

Note that the resulting coordinates (x0', x1', x2') may not be identical to the original
(x0, x1, x2). For example, consider spherical coordinates (r, theta, phi) = (-0.1, pi/4, pi/4).
In this case, upon transformation, we obtain (r', theta', phi'), where r' is always positive,
as it is derived from the Euclidean distance, r = sqrt(Cartx^2 + Carty^2 + Cartz^2),
which is always non-negative. This makes the original (x0, x1, x2) an "inner boundary point."
On a cell-centered grid, such points will always map to a location within the interior of the grid.

The process of filling such boundary points requires copying data, and when the data represents
vectors or tensors, it might involve multiplying by either +1 or -1 to ensure proper orientation
and consistency in the transformation.
"""
    cfunc_type = "static int"
    name = "EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt"
    params = """const commondata_struct *restrict commondata, REAL *restrict xx[3],
const int i0, const int i1, const int i2,
REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]"""
    body = r"""
  // Step 0: Unpack grid spacings dxx0, dxx1, dxx2
  const REAL dxx0 = commondata->bcstruct_dxx0;
  const REAL dxx1 = commondata->bcstruct_dxx1;
  const REAL dxx2 = commondata->bcstruct_dxx2;

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
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Carty(x1(i1)),Cartz(x2(i2)))
  //         If not, error out!
"""
    # Load up the EigenCoordinate corresponding to reference_metric::CoordSystem
    rfm_orig = refmetric.reference_metric[CoordSystem]
    if rfm_orig.CoordSystem.startswith("GeneralRFM"):
        rfm = rfm_orig
    else:
        rfm = refmetric.reference_metric[rfm_orig.EigenCoord]

    xx_to_Cart_eigen = ccg.c_codegen(
        [rfm.xx_to_Cart[0], rfm.xx_to_Cart[1], rfm.xx_to_Cart[2]],
        ["xCart[0]", "xCart[1]", "xCart[2]"],
    )
    xx_to_Cart_orig = ccg.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ["xCart_from_xx", "yCart_from_xx", "zCart_from_xx"],
        include_braces=False,
    )
    xx_to_Cart_inbounds = ccg.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ["xCart_from_xx_inbounds", "yCart_from_xx_inbounds", "zCart_from_xx_inbounds"],
        include_braces=False,
    )
    cart_to_xx_eigen = ""
    if not rfm.CoordSystem.startswith("GeneralRFM"):
        cart_to_xx_eigen = ccg.c_codegen(
            [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
            ["Cart_to_xx0_inbounds", "Cart_to_xx1_inbounds", "Cart_to_xx2_inbounds"],
        )

    # Step 1: Output C code for the Eigen-Coordinate mapping from xx->Cartesian':
    body += f"""
// Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates
REAL xCart[3];  // where (x,y,z) is output
{{
    // xx_to_Cart for EigenCoordinate {rfm.CoordSystem} (original coord = {rfm_orig.CoordSystem}):
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
""" + xx_to_Cart_eigen + "}\n"
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
REAL Cart_to_xx0_inbounds, Cart_to_xx1_inbounds, Cart_to_xx2_inbounds;
"""
    # Step 2.a/2.b: Cartesian->xx map
    if rfm.CoordSystem.startswith("GeneralRFM"):
        body += f"""
  // BHaHAHA bcstruct_set_up() currently supports spherical eigencoordinates only.
  // Fail explicitly for GeneralRFM* instead of emitting code paths that require params.
  return BCSTRUCT_EIGENCOORD_FAILURE;
"""
    else:
        # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
        if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
            raise RuntimeError(
                f"ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for "
                f"reference_metric::CoordSystem = {CoordSystem}. "
                "Boundary conditions in curvilinear coordinates REQUIRE this to be set."
            )
        # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
        body += f"  // Cart_to_xx for EigenCoordinate {rfm.CoordSystem} (original coord = {rfm_orig.CoordSystem}):\n"
        body += cart_to_xx_eigen
    body += r"""
  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((REAL)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const REAL xxmin[3] = {
    xx[0][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx0,
    xx[1][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx1,
    xx[2][0] - ((REAL)(0-NGHOSTS) + (1.0/2.0))*dxx2 };

  // Finally compute i{0,1,2}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((REAL)NGHOSTS)*dxx0)/dxx0 + 0.5 );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((REAL)NGHOSTS)*dxx1)/dxx1 + 0.5 );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((REAL)NGHOSTS)*dxx2)/dxx2 + 0.5 );
"""

    # Restore reference_metric::CoordSystem back to the original CoordSystem

    # Step 3:
    body += f"""
  // Step 3: Convert x0(i0_inbounds),x1(i1_inbounds),x2(i2_inbounds) -> (Cartx,Carty,Cartz),
  //         and check that
  //         (Cartx,Carty,Cartz) == (Cartx(x0(i0)),Carty(x1(i1)),Cartz(x2(i2)))
  //         If not, error out!

  // Step 3.a: Compute {{x,y,z}}Cart_from_xx, as a
  //           function of i0,i1,i2
  REAL xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {{
    // xx_to_Cart for Coordinate {rfm_orig.CoordSystem}:
    REAL xx0 = xx[0][i0];
    REAL xx1 = xx[1][i1];
    REAL xx2 = xx[2][i2];
"""
    body += xx_to_Cart_orig

    body += f"""  }}

  // Step 3.b: Compute {{x,y,z}}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  REAL xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {{
    // xx_to_Cart_inbounds for Coordinate {rfm_orig.CoordSystem}:
    REAL xx0 = xx[0][i0_inbounds];
    REAL xx1 = xx[1][i1_inbounds];
    REAL xx2 = xx[2][i2_inbounds];
"""
    body += xx_to_Cart_inbounds

    body += r"""  }

  // Step 3.c: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
"""
    if par.parval_from_str("fp_type") == "float":
        body += r"""
#define EPS_REL 1e-6
"""
    else:
        body += r"""
#define EPS_REL 1e-8
"""
    body += r"""
  const REAL norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (double)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (double)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {
    // const int Nxx_plus_2NGHOSTS0 = commondata->bcstruct_Nxx_plus_2NGHOSTS0;
    // const int Nxx_plus_2NGHOSTS1 = commondata->bcstruct_Nxx_plus_2NGHOSTS1;
    // const int Nxx_plus_2NGHOSTS2 = commondata->bcstruct_Nxx_plus_2NGHOSTS2;
    // fprintf(stderr,"Error in Spherical coordinate system: Inner boundary point does not map to grid interior point: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
    //         (double)xCart_from_xx,(double)yCart_from_xx,(double)zCart_from_xx,
    //         (double)xCart_from_xx_inbounds,(double)yCart_from_xx_inbounds,(double)zCart_from_xx_inbounds,
    //         xx[0][i0],xx[1][i1],xx[2][i2],
    //         xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
    //         Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2);
    return BCSTRUCT_EIGENCOORD_FAILURE;
  }

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;

  return BHAHAHA_SUCCESS;
#undef EPS_REL
"""
    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cf.full_function


# set_parity_for_inner_boundary_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__set_parity_for_inner_boundary_single_pt(CoordSystem: str) -> str:
    """
    Generate C code for setting the parity for an inner boundary point in a given coordinate system.

    This function computes the parity conditions for all 10 tensor types supported by NRPy, plus 18 for h_{ij,k}.
    at a single grid point on the inner boundary.

    :param CoordSystem: Coordinate system in which to set the parity conditions.
    :return: Full C function code as a string.
    """
    desc = """set_parity_for_inner_boundary_single_pt():
Given (x0,x1,x2)=(xx0,xx1,xx2) and
(x0,x1,x2)'=(x0x1x2_inbounds[0],x0x1x2_inbounds[1],x0x1x2_inbounds[2])
(see description of
EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt()
above for more details), here we compute the parity conditions
for all 10 tensor types supported by NRPy, plus 18 for h_{ij,k}."""
    cfunc_type = "static int"
    name = "set_parity_for_inner_boundary_single_pt"
    params = """const commondata_struct *restrict commondata,
            const REAL xx0,const REAL xx1,const REAL xx2,  const REAL x0x1x2_inbounds[3], const int idx,
            innerpt_bc_struct *restrict innerpt_bc_arr"""
    if par.parval_from_str("fp_type") == "float":
        body = r"""
#define EPS_REL 1e-6
"""
    else:
        body = r"""
#define EPS_REL 1e-8
"""
    body += rf"""
const REAL xx0_inbounds = x0x1x2_inbounds[0];
const REAL xx1_inbounds = x0x1x2_inbounds[1];
const REAL xx2_inbounds = x0x1x2_inbounds[2];

REAL REAL_parity_array[28];
{{
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy
    {parity_conditions_symbolic_dot_products(CoordSystem)}
}}
    // Next perform sanity check on parity array output: should be +1 or -1 to within 8 significant digits:
    for(int whichparity=0; whichparity<28; whichparity++) {{
        if( fabs(REAL_parity_array[whichparity]) < 1 - EPS_REL || fabs(REAL_parity_array[whichparity]) > 1 + EPS_REL ) {{
            fprintf(stderr,"Error at point (%e %e %e), which maps to (%e %e %e).\n",
                xx0,xx1,xx2, xx0_inbounds,xx1_inbounds,xx2_inbounds);
            fprintf(stderr,"Parity evaluated to %e , which is not within 8 significant digits of +1 or -1.\n",
                REAL_parity_array[whichparity]);
            return BCSTRUCT_SET_PARITY_ERROR;
        }}
        innerpt_bc_arr[idx].parity[whichparity] = 1;
        if(REAL_parity_array[whichparity] < 0) innerpt_bc_arr[idx].parity[whichparity] = -1;
    }} // END for(int whichparity=0;whichparity<28;whichparity++)
  return BHAHAHA_SUCCESS;
#undef EPS_REL
"""
    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cf.full_function


# bcstruct_set_up():
#      This function is documented in desc= and body= fields below.
def register_CFunction_bcstruct_set_up(
    CoordSystem: str,
) -> None:
    """
    Register the C function responsible for setting up the source grid boundary condition structure (bcstruct).

    This function sets up the inner and outer boundary points on the computational grid for the given coordinate system.

    :param CoordSystem: The coordinate system for which boundary conditions are to be set up.
    """
    # Step 1: Register contributions of this function to BHaH_defines.h
    BHaH_defines = r"""
// NRPy Curvilinear Boundary Conditions: Core data structures
typedef struct {
  int dstpt;         // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
  int srcpt;         // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
  int8_t parity[28]; // parity[28] is a calculation of dot products for the 28 independent parity types
} innerpt_bc_struct;

typedef struct {
  short i0, i1, i2;              // the outer boundary point grid index (i0,i1,i2), on the 3D grid
  int8_t FACEX0, FACEX1, FACEX2; // 1-byte integers that store
  //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
  //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
  //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i2min face,
  //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
  //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
} outerpt_bc_struct;

typedef struct {
  int num_inner_boundary_points;                  // stores total number of inner boundary points
  int num_pure_outer_boundary_points[NGHOSTS][3]; // stores number of outer boundary points on each
  //                                                  ghostzone level and direction (update min and
  //                                                  max faces simultaneously on multiple cores)
  int bc_loop_bounds[NGHOSTS][6][6]; // stores outer boundary loop bounds. Unused after bcstruct_set_up()
} bc_info_struct;

typedef struct {
  innerpt_bc_struct *inner_bc_array;                   // information needed for updating each inner boundary point
  outerpt_bc_struct *pure_outer_bc_array[NGHOSTS * 3]; // information needed for updating each outer
  //                                                                        boundary point
  bc_info_struct bc_info; // stores number of inner and outer boundary points, needed for setting loop
  //                                  bounds and parallelizing over as many boundary points as possible.
} bc_struct;
"""
    BHaH.BHaH_defines_h.register_BHaH_defines(
        "nrpy.infrastructures.BHaH.CurviBoundaryConditions.BHaH_defines",
        # Puts these definitions near the top of BHaH_defines.h, before grid.
        BHaH_defines,
    )

    # Step 2: Register contributions of this function to the commondata struct.
    for i in range(3):
        BHaH.griddata_commondata.register_griddata_commondata(
            __name__,
            f"REAL bcstruct_Nxx_plus_2NGHOSTS{i}",
            f"The Nxx_plus_2NGHOSTS{i} used when setting up bcstruct",
            is_commondata=True,
        )
        BHaH.griddata_commondata.register_griddata_commondata(
            __name__,
            f"REAL bcstruct_dxx{i}",
            f"The dxx{i} used when setting up bcstruct",
            is_commondata=True,
        )

    # Step 3: Register bcstruct_set_up().
    includes = ["BHaH_defines.h"]
    prefunc = Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
        CoordSystem
    )
    prefunc += Cfunction__set_parity_for_inner_boundary_single_pt(CoordSystem)
    desc = r"""At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
*Step 1: Set up inner boundary structs bcstruct->inner_bc_array[].
*  Recall that at each inner boundary point we must set innerpt_bc_struct:
*    typedef struct {
*      int dstpt;  // dstpt is the 3D grid index IDX3(i0,i1,i2) of the inner boundary point (i0,i1,i2)
*      int srcpt;  // srcpt is the 3D grid index (a la IDX3) to which the inner boundary point maps
*      int8_t parity[28];  // parity[28] is a calculation of dot products for the 28 independent parity types
*    } innerpt_bc_struct;
*  At each ghostzone (i.e., each point within NGHOSTS points from grid boundary):
*    Call EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt().
*        This function converts the curvilinear coordinate (x0,x1,x2) to the corresponding
*        Cartesian coordinate (x,y,z), then finds the grid point
*        (i0_inbounds,i1_inbounds,i2_inbounds) in the grid interior or outer boundary
*        corresponding to this Cartesian coordinate (x,y,z).
*    If (i0,i1,i2) *is not* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
*        then we are at an inner boundary point. We must set
*        bcstruct->inner_bc_array for this point, which requires we specify
*        both (i0_inbounds,i1_inbounds,i2_inbounds) [just found!] and parity
*        conditions for this gridpoint. The latter is found & specified within the
*        function set_parity_for_inner_boundary_single_pt().
*    If (i0,i1,i2) *is* the same as (i0_inbounds,i1_inbounds,i2_inbounds),
*        then we are at an outer boundary point. Take care of outer BCs in Step 2.
*Step 2: Set up outer boundary structs bcstruct->outer_bc_array[which_gz][face][idx2d]:
*  Recall that at each outer boundary point we must set outerpt_bc_struct:
*    typedef struct {
*      short i0,i1,i2;  // the outer boundary point grid index (i0,i1,i2), on the 3D grid
*      int8_t FACEX0,FACEX1,FACEX2;  // 1-byte integers that store
*      //                               FACEX0,FACEX1,FACEX2 = +1, 0, 0 if on the i0=i0min face,
*      //                               FACEX0,FACEX1,FACEX2 = -1, 0, 0 if on the i0=i0max face,
*      //                               FACEX0,FACEX1,FACEX2 =  0,+1, 0 if on the i1=i1min face,
*      //                               FACEX0,FACEX1,FACEX2 =  0,-1, 0 if on the i1=i1max face,
*      //                               FACEX0,FACEX1,FACEX2 =  0, 0,+1 if on the i2=i2min face, or
*      //                               FACEX0,FACEX1,FACEX2 =  0, 0,-1 if on the i2=i2max face,
*    } outerpt_bc_struct;
*  Outer boundary points are filled from the inside out, two faces at a time.
*    E.g., consider a Cartesian coordinate grid that has 14 points in each direction,
*    including the ghostzones, with NGHOSTS=2.
*    We first fill in the lower x0 face with (i0=1,i1={2,11},i2={2,11}). We fill these
*    points in first, since they will in general (at least in the case of extrapolation
*    outer BCs) depend on e.g., i0=2 and i0=3 points.
*    Simultaneously we can fill in the upper x0 face with (i0=12,i1={2,11},i2={2,11}),
*    since these points depend only on e.g., i0=11 and i0=10 (again assuming extrap. BCs).
*    Next we can fill in the lower x1 face: (i0={1,12},i1=2,i2={2,11}). Notice these
*    depend on i0 min and max faces being filled. The remaining pattern goes like this:
*    Upper x1 face: (i0={1,12},i1=12,i2={2,11})
*    Lower x2 face: (i0={1,12},i1={1,12},i2=1)
*    Upper x2 face: (i0={1,12},i1={1,12},i2=12)
*    Lower x0 face: (i0=0,i1={1,12},i2={1,12})
*    Upper x0 face: (i0=13,i1={1,12},i2={1,12})
*    Lower x1 face: (i0={0,13},i1=0,i2={2,11})
*    Upper x1 face: (i0={0,13},i1=13,i2={2,11})
*    Lower x2 face: (i0={0,13},i1={0,13},i2=0)
*    Upper x2 face: (i0={0,13},i1={0,13},i2=13)
*  Note that we allocate an outerpt_bc_struct at *all* boundary points,
*    regardless of whether the point is an outer or inner point. However
*    the struct is set only at outer boundary points. This is slightly
*    wasteful, but only in memory, not in CPU."""
    cfunc_type = "int"
    name = "bcstruct_set_up"
    params = """const commondata_struct *restrict commondata, REAL *restrict xx[3],
                bc_struct *restrict bcstruct"""
    body = r"""
  const int Nxx_plus_2NGHOSTS0 = commondata->bcstruct_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->bcstruct_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->bcstruct_Nxx_plus_2NGHOSTS2;
  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  {
    // First count the number of inner points.
    bool error_flag = false;
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
      if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        if (EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds)) {
#pragma omp critical
          {
            error_flag = true;
          }
          continue; // Skip further processing.
        }
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        }
      }
    }
    if (error_flag)
      return BCSTRUCT_EIGENCOORD_FAILURE;

    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;
"""
    body += r"""
    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *restrict)malloc(sizeof(innerpt_bc_struct) * num_inner);
  }

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
      if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS)) {
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        if (EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds)) {
          return BCSTRUCT_EIGENCOORD_FAILURE;
        }
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3(i0, i1, i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3(i0i1i2_inbounds[0], i0i1i2_inbounds[1], i0i1i2_inbounds[2]);
          // printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          if (set_parity_for_inner_boundary_single_pt(commondata, xx[0][i0], xx[1][i1], xx[2][i2], x0x1x2_inbounds, which_inner,
                                                      bcstruct->inner_bc_array)) {
            return BCSTRUCT_SET_PARITY_ERROR;
          }

          which_inner++;
        }
      }
    }
  }

  ////////////////////////////////////////
  // STEP 2: SET UP OUTER BOUNDARY STRUCTS
  // First set up loop bounds for outer boundary condition updates,
  //   store to bc_info->bc_loop_bounds[which_gz][face][]. Also
  //   allocate memory for outer_bc_array[which_gz][face][]:
  int imin[3] = {NGHOSTS, NGHOSTS, NGHOSTS};
  int imax[3] = {Nxx_plus_2NGHOSTS0 - NGHOSTS, Nxx_plus_2NGHOSTS1 - NGHOSTS, Nxx_plus_2NGHOSTS2 - NGHOSTS};
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    const int x0min_face_range[6] = {imin[0] - 1, imin[0], imin[1], imax[1], imin[2], imax[2]};
    imin[0]--;
    const int x0max_face_range[6] = {imax[0], imax[0] + 1, imin[1], imax[1], imin[2], imax[2]};
    imax[0]++;
    const int x1min_face_range[6] = {imin[0], imax[0], imin[1] - 1, imin[1], imin[2], imax[2]};
    imin[1]--;
    const int x1max_face_range[6] = {imin[0], imax[0], imax[1], imax[1] + 1, imin[2], imax[2]};
    imax[1]++;
    const int x2min_face_range[6] = {imin[0], imax[0], imin[1], imax[1], imin[2] - 1, imin[2]};
    imin[2]--;
    const int x2max_face_range[6] = {imin[0], imax[0], imin[1], imax[1], imax[2], imax[2] + 1};
    imax[2]++;

    int face = 0;
    ////////////////////////
    // x0min and x0max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x0min and x0max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x0min_face_range[1] - x0min_face_range[0]) * (x0min_face_range[3] - x0min_face_range[2]) * (x0min_face_range[5] - x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i];
    }
    face++;
    // x0max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i];
    }
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x1min_face_range[1] - x1min_face_range[0]) * (x1min_face_range[3] - x1min_face_range[2]) * (x1min_face_range[5] - x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i];
    }
    face++;
    // x1max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i];
    }
    face++;
    ////////////////////////

    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    //                   Also, note that face/2 --v   offsets this factor of 2 ------------------------------------------v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *restrict)malloc(
        sizeof(outerpt_bc_struct) * 2 *
        ((x2min_face_range[1] - x2min_face_range[0]) * (x2min_face_range[3] - x2min_face_range[2]) * (x2min_face_range[5] - x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i];
    }
    face++;
    // x2max face: Set loop bounds & allocate memory for outer_bc_array:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i];
    }
    face++;
    ////////////////////////
  } // END LOOP over ghostzones

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++)
    for (int dirn = 0; dirn < 3; dirn++) {
      int idx2d = 0;
      // LOWER FACE: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn * 2;
#define IDX2D_BCS(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                                              \
  (((i0) - (i0min)) + ((i0max) - (i0min)) * (((i1) - (i1min)) + ((i1max) - (i1min)) * ((i2) - (i2min))))
        const int FACEX0 = (face == 0) - (face == 1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1 = (face == 2) - (face == 3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2 = (face == 4) - (face == 5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0, bcstruct->bc_info.bc_loop_bounds[which_gz][face][0], bcstruct->bc_info.bc_loop_bounds[which_gz][face][1], i1,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][2], bcstruct->bc_info.bc_loop_bounds[which_gz][face][3], i2,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][4], bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          if (EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds)) {
            return BCSTRUCT_EIGENCOORD_FAILURE;
          }
          if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      } // END LOOP over lower faces
      // UPPER FACE: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
      {
        const int face = dirn * 2 + 1;
        const int FACEX0 = (face == 0) - (face == 1); // +1 if face==0 ; -1 if face==1. Otherwise 0.
        const int FACEX1 = (face == 2) - (face == 3); // +1 if face==2 ; -1 if face==3. Otherwise 0.
        const int FACEX2 = (face == 4) - (face == 5); // +1 if face==4 ; -1 if face==5. Otherwise 0.
        LOOP_NOOMP(i0, bcstruct->bc_info.bc_loop_bounds[which_gz][face][0], bcstruct->bc_info.bc_loop_bounds[which_gz][face][1], i1,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][2], bcstruct->bc_info.bc_loop_bounds[which_gz][face][3], i2,
                   bcstruct->bc_info.bc_loop_bounds[which_gz][face][4], bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          if (EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds)) {
            return BCSTRUCT_EIGENCOORD_FAILURE;
          }
          if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          }
        }
      } // END LOOP over upper faces
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    } // END LOOPS over directions and ghost zone layers.

  return BHAHAHA_SUCCESS;
"""
    cfc.register_CFunction(
        subdirectory="",
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
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

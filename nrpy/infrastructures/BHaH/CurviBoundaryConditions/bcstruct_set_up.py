# nrpy/infrastructures/BHaH/CurviBoundaryConditions/bcstruct_set_up.py
"""
Generates C code to set up the boundary condition data structure (`bc_struct`).

This module creates the C function `bcstruct_set_up`, which populates the
`bc_struct` with all necessary information to handle inner and outer boundary
conditions on a curvilinear grid. The setup process involves:

1.  **Differentiating inner vs. outer boundaries:** It uses an "Eigencoordinate"
    mapping to robustly identify whether a ghost zone point maps to the grid
    interior (an inner boundary) or to itself (an outer boundary).
2.  **Storing inner boundary data:** For each inner boundary point, it stores the
    source grid point and the tensor-type-specific parity transformations.
3.  **Organizing outer boundary data:** It catalogs outer boundary points by face
    and layer to facilitate sequential application.

This pre-computation allows for efficient application of boundary conditions
at runtime by other functions.

This process is documented in the NRPy tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends

from nrpy import c_codegen, c_function, indexedexp
from nrpy import params as par
from nrpy import reference_metric


# Set unit-vector dot products (=parity) for each of the 10 parity condition types
def parity_conditions_symbolic_dot_products(
    CoordSystem: str,
) -> str:
    """
    Set unit-vector dot products (parity) for each of the 10 parity condition types.

    :param CoordSystem: The coordinate system for which to calculate the parity conditions.
    :return: C code string representing the unit vector dot products for the ten parity conditions.
    """
    rfm = reference_metric.reference_metric[CoordSystem]
    parity = indexedexp.zerorank1(dimension=10)
    UnitVectors_inner = indexedexp.zerorank2()
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

    # We enforce double or higher precision to avoid issues in finding
    # the inner/outer BC indicies
    fp_type = par.parval_from_str("fp_type")
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"
    calculation = c_codegen.c_codegen(parity, lhs_strings, fp_type_alias=fp_type_alias)

    return outstr + calculation


# set_parity_for_inner_boundary_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__set_parity_for_inner_boundary_single_pt(
    CoordSystem: str,
) -> str:
    """
    Generate C code for setting the parity for inner boundary single point in a given coordinate system.

    :param CoordSystem: Coordinate system in which to set the parity
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
    fp_type = par.parval_from_str("fp_type")
    if fp_type == "float":
        body = r"""
#define EPS_REL 1e-6
"""
    else:
        body = r"""
#define EPS_REL 1e-8
"""
    body += f"""
const DOUBLE xx0_inbounds = x0x1x2_inbounds[0];
const DOUBLE xx1_inbounds = x0x1x2_inbounds[1];
const DOUBLE xx2_inbounds = x0x1x2_inbounds[2];

DOUBLE REAL_parity_array[10];
{{
    // Evaluate dot products needed for setting parity
    //     conditions at a given point (xx0,xx1,xx2),
    //     using C code generated by NRPy+
{parity_conditions_symbolic_dot_products(CoordSystem)}
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
    cf = c_function.CFunction(
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


# EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt():
#      This function is documented in desc= and body= fields below.
def Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
    CoordSystem: str,
) -> str:
    """
    Identify and map ghost zone grid points to their corresponding inner boundary points.

    Map a reference metric grid point index (i0, i1, i2) in a ghost zone to an interior point index (i0, i1, i2)':
    (i0, i1, i2) -> (i0, i1, i2)',
    if it is an inner boundary point. If the grid point maps to itself; i.e.,
    (i0, i1, i2) -> (i0, i1, i2),
    it should have been marked as an outer boundary point. This process involves the following double-map:
    (x0, x1, x2) -> (Cartx, Carty, Cartz) -> (x0, x1, x2)'
    However, the second map from Cartesian to the reference metric does not always have a closed-form expression,
    and this simple algorithm will fail. To address this, we exploit the fact that an arbitrary reference metric and
    its "eigencoordinate" share the exact same index mapping:
    (i0, i1, i2) -> (i0, i1, i2)'
    Therefore, while the mapping
    (Cartx, Carty, Cartz) -> (x0, x1, x2)'
    may not be closed-form for the chosen CoordSystem, the eigencoordinate mapping
    (Cartx, Carty, Cartz) -> (x0, x1, x2)'
    will be, for all reference metrics in NRPy.

    :param CoordSystem: The coordinate system for mapping.

    :return: Body of the C code.

    :raises RuntimeError: If the conversion back to Cartesian coordinates does not match the original coordinates, indicating an error in the mapping process.
    """
    desc = """Map a reference metric grid point index (i0, i1, i2) in a ghost zone to an interior point index (i0, i1, i2)':
(i0, i1, i2) -> (i0, i1, i2)',
if it is an inner boundary point. If the grid point maps to itself; i.e.,
(i0, i1, i2) -> (i0, i1, i2),
it should have been marked as an outer boundary point. This process involves the following double-map:
(x0, x1, x2) -> (Cartx, Carty, Cartz) -> (x0, x1, x2)'
However, the second map from Cartesian to the reference metric does not always have a closed-form expression,
and this simple algorithm will fail. To address this, we exploit the fact that an arbitrary reference metric and
its "eigencoordinate" share the exact same index mapping:
(i0, i1, i2) -> (i0, i1, i2)'
Therefore, while the mapping
(Cartx, Carty, Cartz) -> (x0, x1, x2)'
may not be closed-form for the chosen CoordSystem, the eigencoordinate mapping
(Cartx, Carty, Cartz) -> (x0, x1, x2)'
will be, for all reference metrics in NRPy.

Definition of Eigencoordinate:

A coordinate system's "eigencoordinate" is the simplest member of its family:
- Spherical-like: Spherical
- Cylindrical-like: Cylindrical
- Cartesian-like: Cartesian
- SymTP-like: SymTP

Key Steps:
1. Convert to Cartesian Coordinates:
   - Transform (x0, x1, x2) from eigencoordinates to Cartesian (Cartx, Carty, Cartz).
2. Map Back to Eigencoordinates:
   - Convert (Cartx, Carty, Cartz) back to eigencoordinates (x0', x1', x2').
3. Sanity Check and Data Handling:
   - If (x0, x1, x2) != (x0', x1', x2'), the point is an inner boundary point.
     - For example, in Spherical coordinates, a negative radius becomes positive.
   - On a cell-centered grid, the mapped point lies within the grid interior.
     - Update the data at (i0, i1, i2) by copying from (i0_inbounds, i1_inbounds, i2_inbounds).
     - Apply a sign change (+1 or -1) if the data represents tensors or vectors.
"""
    cfunc_type = "static void"
    name = "EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params, REAL *restrict xx[3],
const int i0, const int i1, const int i2,
REAL x0x1x2_inbounds[3], int i0i1i2_inbounds[3]"""
    fp_type = par.parval_from_str("fp_type")
    fp_type_alias = "DOUBLE" if fp_type == "float" else "REAL"
    # Load up the EigenCoordinate corresponding to reference_metric::CoordSystem
    rfm_orig = reference_metric.reference_metric[CoordSystem]
    rfm = reference_metric.reference_metric[rfm_orig.EigenCoord]

    # Step 1: Output C code for the Eigen-Coordinate mapping from xx->Cartesian':
    body = rf"""
  // Step 1: Convert the (curvilinear) coordinate (x0,x1,x2) to Cartesian coordinates:
  //         (x0,x1,x2) -> (Cartx,Carty,Cartz)
  //         Find the Cartesian coordinate that (x0,x1,x2)
  //         maps to, assuming (x0,x1,x2) is the eigen-
  //         coordinate. Note that we assume (x0,x1,x2)
  //         has the same grid boundaries in both the
  //         original coordinate and the eigencoordinate.
  DOUBLE xCart[3];  // where (x,y,z) is output
  {{
    // xx_to_Cart for EigenCoordinate {rfm.CoordSystem} (orig coord = {rfm_orig.CoordSystem}):
    DOUBLE xx0 = xx[0][i0];
    DOUBLE xx1 = xx[1][i1];
    DOUBLE xx2 = xx[2][i2];
    {c_codegen.c_codegen([rfm.xx_to_Cart[0], rfm.xx_to_Cart[1], rfm.xx_to_Cart[2]],
                         ["xCart[0]", "xCart[1]", "xCart[2]"],
                         fp_type_alias=fp_type_alias)}
  }}
  DOUBLE Cartx = xCart[0];
  DOUBLE Carty = xCart[1];
  DOUBLE Cartz = xCart[2];

  // Step 2: Convert Cartesian coordinates (Cartx, Carty, Cartz) to eigencoordinates (x0, x1, x2)
  //         and determine the corresponding gridpoint indices (i0, i1, i2).
  //         - For cell-centered grids, (x0, x1, x2) will exactly (within roundoff error) align
  //           with a numerical grid point.
  //         - Identify the in-bounds indices (i0_inbounds, i1_inbounds, i2_inbounds) that
  //           correspond to the eigencoordinates.
  //         - Check the location of (i0_inbounds, i1_inbounds, i2_inbounds):
  //             - If it resides within a ghost zone:
  //                 - It must be equal to (i0, i1, i2), indicating that the point is an outer boundary point.
  //             - Otherwise:
  //                 - The indices are within the grid interior.
  //                 - Replace the data at (i0, i1, i2) with the data from (i0_inbounds, i1_inbounds, i2_inbounds),
  //                   multiplied by the appropriate parity condition (+1 or -1).
  DOUBLE Cart_to_xx0_inbounds,Cart_to_xx1_inbounds,Cart_to_xx2_inbounds;
"""
    # Step 2.a: Sanity check: First make sure that rfm.Cart_to_xx has been set. Error out if not!
    if rfm.Cart_to_xx[0] == 0 or rfm.Cart_to_xx[1] == 0 or rfm.Cart_to_xx[2] == 0:
        raise RuntimeError(
            f"ERROR: rfm.Cart_to_xx[], which maps Cartesian -> xx, has not been set for "
            f"reference_metric::CoordSystem = {CoordSystem}. "
            "Boundary conditions in curvilinear coordinates REQUiRE this be set."
        )
    # Step 2.b: Output C code for the Eigen-Coordinate mapping from Cartesian->xx:
    body += rf"""
  // Cart_to_xx for EigenCoordinate {rfm.CoordSystem} (orig coord = {rfm_orig.CoordSystem})
  {c_codegen.c_codegen(
        [rfm.Cart_to_xx[0], rfm.Cart_to_xx[1], rfm.Cart_to_xx[2]],
        ['Cart_to_xx0_inbounds', 'Cart_to_xx1_inbounds', 'Cart_to_xx2_inbounds'],
        fp_type_alias=fp_type_alias,
    )}
  // Next compute xxmin[i]. By definition,
  //    xx[i][j] = xxmin[i] + ((DOUBLE)(j-NGHOSTS) + (1.0/2.0))*dxxi;
  // -> xxmin[i] = xx[i][0] - ((DOUBLE)(0-NGHOSTS) + (1.0/2.0))*dxxi
  const DOUBLE xxmin[3] = {{
    xx[0][0] - ((DOUBLE)(0-NGHOSTS) + (1.0/2.0))*dxx0,
    xx[1][0] - ((DOUBLE)(0-NGHOSTS) + (1.0/2.0))*dxx1,
    xx[2][0] - ((DOUBLE)(0-NGHOSTS) + (1.0/2.0))*dxx2 }};

  // Finally compute i{{0,1,2}}_inbounds (add 0.5 to account for rounding down)
  const int i0_inbounds = (int)( (Cart_to_xx0_inbounds - xxmin[0] - (1.0/2.0)*dxx0 + ((DOUBLE)NGHOSTS)*dxx0)/dxx0 + 0.5 );
  const int i1_inbounds = (int)( (Cart_to_xx1_inbounds - xxmin[1] - (1.0/2.0)*dxx1 + ((DOUBLE)NGHOSTS)*dxx1)/dxx1 + 0.5 );
  const int i2_inbounds = (int)( (Cart_to_xx2_inbounds - xxmin[2] - (1.0/2.0)*dxx2 + ((DOUBLE)NGHOSTS)*dxx2)/dxx2 + 0.5 );

  // Step 3: Sanity Check
  //         - Convert eigencoordinates x0(i0_inbounds), x1(i1_inbounds), x2(i2_inbounds)
  //           to Cartesian coordinates (Cartx, Carty, Cartz).
  //         - Verify that the converted coordinates match the original mapping:
  //               (Cartx, Carty, Cartz) == (Cartx(x0(i0)), Carty(x1(i1)), Cartz(x2(i2))).
  //         - If the coordinates do not match, trigger an error.

  // Step 3.a: Compute {{x,y,z}}Cart_from_xx, as a function of i0,i1,i2
  DOUBLE xCart_from_xx, yCart_from_xx, zCart_from_xx;
  {{
    // xx_to_Cart for Coordinate {rfm_orig.CoordSystem}):
    DOUBLE xx0 = xx[0][i0];
    DOUBLE xx1 = xx[1][i1];
    DOUBLE xx2 = xx[2][i2];
    {c_codegen.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ['xCart_from_xx', 'yCart_from_xx', 'zCart_from_xx'],
        include_braces=False,
        fp_type_alias=fp_type_alias,
    )}
  }}

  // Step 3.b: Compute {{x,y,z}}Cart_from_xx_inbounds, as a
  //           function of i0_inbounds,i1_inbounds,i2_inbounds
  DOUBLE xCart_from_xx_inbounds, yCart_from_xx_inbounds, zCart_from_xx_inbounds;
  {{
    // xx_to_Cart_inbounds for Coordinate {rfm_orig.CoordSystem}):
    DOUBLE xx0 = xx[0][i0_inbounds];
    DOUBLE xx1 = xx[1][i1_inbounds];
    DOUBLE xx2 = xx[2][i2_inbounds];
    {c_codegen.c_codegen(
        [rfm_orig.xx_to_Cart[0], rfm_orig.xx_to_Cart[1], rfm_orig.xx_to_Cart[2]],
        ['xCart_from_xx_inbounds', 'yCart_from_xx_inbounds', 'zCart_from_xx_inbounds'],
        include_braces=False,
        fp_type_alias=fp_type_alias,
    )}
  }}

  // Step 3.c: Compare xCart_from_xx to xCart_from_xx_inbounds;
  //           they should be identical!!!
#define EPS_REL {"1e-6" if fp_type == "float" else "1e-8"}
  const DOUBLE norm_factor = sqrt(xCart_from_xx*xCart_from_xx + yCart_from_xx*yCart_from_xx + zCart_from_xx*zCart_from_xx) + 1e-15;
  if(fabs( (DOUBLE)(xCart_from_xx - xCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (DOUBLE)(yCart_from_xx - yCart_from_xx_inbounds) ) > EPS_REL * norm_factor ||
     fabs( (DOUBLE)(zCart_from_xx - zCart_from_xx_inbounds) ) > EPS_REL * norm_factor) {{
    fprintf(stderr,"Error in {rfm_orig.CoordSystem} coordinate system: Inner boundary point does not map to grid interior point: ( %.15e %.15e %.15e ) != ( %.15e %.15e %.15e ) | xx: %e %e %e -> %e %e %e | %d %d %d\n",
            (DOUBLE)xCart_from_xx,(DOUBLE)yCart_from_xx,(DOUBLE)zCart_from_xx,
            (DOUBLE)xCart_from_xx_inbounds,(DOUBLE)yCart_from_xx_inbounds,(DOUBLE)zCart_from_xx_inbounds,
            xx[0][i0],xx[1][i1],xx[2][i2],
            xx[0][i0_inbounds],xx[1][i1_inbounds],xx[2][i2_inbounds],
            Nxx_plus_2NGHOSTS0,Nxx_plus_2NGHOSTS1,Nxx_plus_2NGHOSTS2);
    exit(1);
  }}
#undef EPS_REL

  // Step 4: Set output arrays.
  x0x1x2_inbounds[0] = xx[0][i0_inbounds];
  x0x1x2_inbounds[1] = xx[1][i1_inbounds];
  x0x1x2_inbounds[2] = xx[2][i2_inbounds];
  i0i1i2_inbounds[0] = i0_inbounds;
  i0i1i2_inbounds[1] = i1_inbounds;
  i0i1i2_inbounds[2] = i2_inbounds;
"""
    return c_function.CFunction(
        subdirectory=CoordSystem,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=True,
        body=body,
    ).full_function


# bcstruct_set_up():
#      This function is documented in desc= and body= fields below.
def register_CFunction_bcstruct_set_up(
    CoordSystem: str, enable_masks: bool = False
) -> None:
    """
    Register C function for setting up bcstruct.

    This function prescribes how inner and outer boundary points on the
    computational grid are filled, based on the given coordinate system (CoordSystem).

    :param CoordSystem: The coordinate system for which to set up boundary conditions.
    :param enable_masks: If True, make bcstruct algorithm mask-aware.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "bcstruct_set_up__rfm"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       c_function.CFunction_dict.clear()
    ...       register_CFunction_bcstruct_set_up(CoordSystem)
    ...       generated_str = c_function.CFunction_dict[f'{name}__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[SymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Spherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    Setting up reference_metric[Cylindrical]...
    """
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
    ]
    prefunc = Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(
        CoordSystem=(CoordSystem if CoordSystem != "Fisheye" else "Cartesian")
    )
    prefunc += Cfunction__set_parity_for_inner_boundary_single_pt(
        CoordSystem=(CoordSystem if CoordSystem != "Fisheye" else "Cartesian")
    )
    parallelization = par.parval_from_str("parallelization")
    desc = r"""At each coordinate point (x0,x1,x2) situated at grid index (i0,i1,i2):
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
    cfunc_type = "void"
    name = "bcstruct_set_up"
    params = "const commondata_struct *restrict commondata, const params_struct *restrict params"
    params += ", const int8_t *restrict mask" if enable_masks else ""
    params += ", REAL *restrict xx[3], bc_struct *restrict bcstruct"
    params += (
        "" if parallelization != "cuda" else ", bc_struct *restrict bcstruct_device"
    )

    # Setup host-side struct to populate before copying to device
    body = """
  ////////////////////////////////////////
  // STEP 1: SET UP INNER BOUNDARY STRUCTS
  // First count the number of inner boundary points and allocate memory for inner_bc_array.
  {
    int num_inner = 0;
    LOOP_OMP("omp parallel for reduction(+:num_inner)", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
"""
    unset_mask_check = ""
    if enable_masks:
        unset_mask_check = "&& mask[IDX3(i0, i1, i2)] != UNSET"
    body += f"if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS) {unset_mask_check}) {{"
    body += """
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          // this is an inner boundary point, which maps either
          //  to the grid interior or to an outer boundary point
          num_inner++;
        } // END IF boundary point maps to itself, then outer; otherwise inner
      } // END IF point lies on grid boundary
    } // END LOOP over all points
    // Store num_inner to bc_info:
    bcstruct->bc_info.num_inner_boundary_points = num_inner;

    // Next allocate memory for inner_boundary_points:
    bcstruct->inner_bc_array = (innerpt_bc_struct *)malloc(sizeof(innerpt_bc_struct) * num_inner);
  } // END count number of inner boundary points and allocate memory for inner_bc_array.

  // Then set inner_bc_array:
  {
    int which_inner = 0;
    LOOP_NOOMP(i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int i0i1i2[3] = {i0, i1, i2};
"""
    body += f"if (!IS_IN_GRID_INTERIOR(i0i1i2, Nxx_plus_2NGHOSTS0, Nxx_plus_2NGHOSTS1, Nxx_plus_2NGHOSTS2, NGHOSTS) {unset_mask_check}) {{"
    body += r"""
        REAL x0x1x2_inbounds[3];
        int i0i1i2_inbounds[3];
        EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
        if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2]) {
          // this is a pure outer boundary point.
        } else {
          bcstruct->inner_bc_array[which_inner].dstpt = IDX3(i0, i1, i2);
          bcstruct->inner_bc_array[which_inner].srcpt = IDX3(i0i1i2_inbounds[0], i0i1i2_inbounds[1], i0i1i2_inbounds[2]);
          // printf("%d / %d\n",which_inner, bc_info->num_inner_boundary_points);
          set_parity_for_inner_boundary_single_pt(commondata, params, xx[0][i0], xx[1][i1], xx[2][i2], x0x1x2_inbounds, which_inner,
                                                  bcstruct->inner_bc_array);

          which_inner++;
        } // END IF boundary point maps to itself, then outer; otherwise inner
      } // END IF point lies on grid boundary
    } // END LOOP over all points
  } // END set inner_bc_array.

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
    //                   Also, note that    face/2 ---v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *)malloc(
    //   offsets this factor of 2 --v
        sizeof(outerpt_bc_struct) * 2 *
        ((x0min_face_range[1] - x0min_face_range[0]) * (x0min_face_range[3] - x0min_face_range[2]) * (x0min_face_range[5] - x0min_face_range[4])));
    // x0min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0min_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    // x0max face: Set loop bounds:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x0max_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    ////////////////////////

    ////////////////////////
    // x1min and x1max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x1min and x1max faces have exactly the same size.
    //                   Also, note that    face/2 ---v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *)malloc(
    //   offsets this factor of 2 --v
        sizeof(outerpt_bc_struct) * 2 *
        ((x1min_face_range[1] - x1min_face_range[0]) * (x1min_face_range[3] - x1min_face_range[2]) * (x1min_face_range[5] - x1min_face_range[4])));
    // x1min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1min_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    // x1max face: Set loop bounds:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x1max_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    ////////////////////////

    ////////////////////////
    // x2min and x2max faces: Allocate memory for outer_bc_array and set bc_loop_bounds:
    //                        Note that x2min and x2max faces have exactly the same size.
    //                   Also, note that    face/2 ---v
    bcstruct->pure_outer_bc_array[3 * which_gz + face / 2] = (outerpt_bc_struct *)malloc(
    //   offsets this factor of 2 --v
        sizeof(outerpt_bc_struct) * 2 *
        ((x2min_face_range[1] - x2min_face_range[0]) * (x2min_face_range[3] - x2min_face_range[2]) * (x2min_face_range[5] - x2min_face_range[4])));
    // x2min face: Can't set bc_info->bc_loop_bounds[which_gz][face] = { i0min,i0max, ... } since it's not const :(
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2min_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    // x2max face: Set loop bounds:
    for (int i = 0; i < 6; i++) {
      bcstruct->bc_info.bc_loop_bounds[which_gz][face][i] = x2max_face_range[i];
    } // END LOOP over all six faces of the grid
    face++;
    ////////////////////////
  } // END LOOP over ghost zones, from inner to outer.

  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {
    for (int dirn = 0; dirn < 3; dirn++) {
      int idx2d = 0;
      // LOWER FACES: dirn=0 -> x0min; dirn=1 -> x1min; dirn=2 -> x2min
      {
        const int face = dirn * 2;
#define IDX2D_BCS(i0, i0min, i0max, i1, i1min, i1max, i2, i2min, i2max)                                                                              \
  (((i0) - (i0min)) + ((i0max) - (i0min)) * (((i1) - (i1min)) + ((i1max) - (i1min)) * ((i2) - (i2min))))
        const int FACEX0 = (face == 0) - (face == 1); // +1 if face==0 (x0min) ; -1 if face==1 (x0max). Otherwise 0.
        const int FACEX1 = (face == 2) - (face == 3); // +1 if face==2 (x1min) ; -1 if face==3 (x1max). Otherwise 0.
        const int FACEX2 = (face == 4) - (face == 5); // +1 if face==4 (x2min) ; -1 if face==5 (x2max). Otherwise 0.
        LOOP_NOOMP(i0, bcstruct->bc_info.bc_loop_bounds[which_gz][face][0], bcstruct->bc_info.bc_loop_bounds[which_gz][face][1], //
                   i1, bcstruct->bc_info.bc_loop_bounds[which_gz][face][2], bcstruct->bc_info.bc_loop_bounds[which_gz][face][3], //
                   i2, bcstruct->bc_info.bc_loop_bounds[which_gz][face][4], bcstruct->bc_info.bc_loop_bounds[which_gz][face][5]) {
          REAL x0x1x2_inbounds[3];
          int i0i1i2_inbounds[3];
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
"""
    outer_boundary_mask_check = ""
    if enable_masks:
        outer_boundary_mask_check = "&& mask[IDX3(i0, i1, i2)] == OUTER_BOUNDARY"
    body += f"if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2] {outer_boundary_mask_check}) {{"
    body += """
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          } // END IF outer boundary point
        } // END LOOP over all boundary points on lower faces
      } // END BLOCK lower faces
      // UPPER FACES: dirn=0 -> x0max; dirn=1 -> x1max; dirn=2 -> x2max
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
          EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt(commondata, params, xx, i0, i1, i2, x0x1x2_inbounds, i0i1i2_inbounds);
"""
    body += f"if (i0 == i0i1i2_inbounds[0] && i1 == i0i1i2_inbounds[1] && i2 == i0i1i2_inbounds[2] {outer_boundary_mask_check}) {{"
    body += """
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i0 = i0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i1 = i1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].i2 = i2;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX0 = FACEX0;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX1 = FACEX1;
            bcstruct->pure_outer_bc_array[dirn + (3 * which_gz)][idx2d].FACEX2 = FACEX2;
            idx2d++;
          } // END IF outer boundary point
        } // END LOOP over all boundary points on upper faces
      } // END BLOCK upper faces
      bcstruct->bc_info.num_pure_outer_boundary_points[which_gz][dirn] = idx2d;
    } // END LOOP over three directions
  } // END LOOP over NGHOSTS ghost zones, from innermost to outermost ghost zones
"""
    if parallelization == "cuda":
        body += """
        int streamid = params->grid_idx % NUM_STREAMS;
        cpyHosttoDevice_bc_struct(bcstruct, bcstruct_device, streamid);
        cudaDeviceSynchronize();
"""
    c_function.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=True,
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

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
from typing import List

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy+ depends
import sympy.codegen.ast as sp_ast, os

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.params as par  # NRPy+: Parameter interface
import nrpy.reference_metric as refmetric  # NRPy+: Reference metric support
from nrpy.helpers.expr_tree import get_unique_expression_symbols
from nrpy.infrastructures.BHaH.CurviBoundaryConditions.CurviBoundaryConditions import (
    BHaH_defines_set_gridfunction_defines_with_parity_types,
    Cfunction__EigenCoord_set_x0x1x2_inbounds__i0i1i2_inbounds_single_pt,
    Cfunction__set_parity_for_inner_boundary_single_pt,
    get_arb_offset_FD_coeffs_indices,
)
from nrpy.validate_expressions.validate_expressions import check_zero


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
        self.body = ""

    def register(self) -> None:
        """Register CFunction."""
        _, actual_name = cfc.function_name_and_subdir_with_CoordSystem(
            os.path.join("."), self.name, self.CoordSystem
        )
        if not actual_name in cfc.CFunction_dict:
            cfc.register_CFunction(
                prefunc=self.prefunc,
                includes=self.includes,
                desc=self.desc,
                cfunc_type=self.cfunc_type,
                CoordSystem_for_wrapper_func=self.CoordSystem,
                name=self.name,
                params=self.params,
                include_CodeParameters_h=True,
                body=self.body,
            )


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
        self.prefunc = ""

    def register(self) -> None:
        """Register CFunction."""
        cfc.register_CFunction(
            prefunc=self.prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )


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
        self.body = ""
        self.prefunc = ""

    def register(self) -> None:
        """Register CFunction."""
        cfc.register_CFunction(
            prefunc=self.prefunc,
            includes=self.includes,
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=False,
            body=self.body,
        )


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

    def __init__(self, CoordSystem: str, fp_type: str = "double") -> None:
        self.CoordSystem = CoordSystem
        self.fp_type = fp_type
        self.CFunction: cfc.CFunction
        self.include_CodeParameters_h = True
        self.cfunc_decorators = ""

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
            sub_list = get_unique_expression_symbols(
                expr, exclude=[f"xx{i}" for i in range(3)]
            )
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

    def generate_CFunction(self) -> None:
        """Generate CFunction from class definitions."""
        self.CFunction = cfc.CFunction(
            subdirectory=self.CoordSystem,
            includes=[],
            desc=self.desc,
            cfunc_type=self.cfunc_type,
            name=self.name,
            params=self.params,
            include_CodeParameters_h=self.include_CodeParameters_h,
            body=self.body,
            cfunc_decorators=self.cfunc_decorators,
        )


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
    """

    def __init__(
        self,
        dirn: int,
        radiation_BC_fd_order: int = -1,
        fp_type: str = "double",
        rational_const_alias: str = "const",
        fp_type_alias: str = "REAL",
    ) -> None:
        self.dirn = dirn
        self.radiation_BC_fd_order = radiation_BC_fd_order
        self.fp_type = fp_type

        self.default_FDORDER = par.parval_from_str("fd_order")
        if radiation_BC_fd_order == -1:
            radiation_BC_fd_order = self.default_FDORDER

        par.set_parval_from_str("fd_order", radiation_BC_fd_order)

        self.include_CodeParameters_h = True
        self.includes: List[str] = []
        self.desc = (
            "Compute 1st derivative finite-difference derivative with arbitrary upwind"
        )
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
            self.body += f"case {offset}:\n {{\n"
            coeffs, indices = get_arb_offset_FD_coeffs_indices(
                radiation_BC_fd_order, offset, 1
            )

            # Build dictionary of coefficients to enable strong typing
            # and assignment to Rational declarations.
            rational_dict = {}
            for i, coeff in enumerate(coeffs):
                if coeff == 0:
                    continue

                # We store the absolute value in the dictionary to avoid
                # duplicate assignments
                sign = -1 if coeff <= 0 else 1
                decl_coeff = coeff * sign
                if decl_coeff in rational_dict:
                    continue
                v = f"Rational_decl__{decl_coeff.p}_{decl_coeff.q}"
                rational_dict[decl_coeff] = v
                RATIONAL_assignment = f"{rational_const_alias} {fp_type_alias} {str(v)}"
                self.body += sp.ccode(
                    decl_coeff,
                    assign_to=RATIONAL_assignment,
                    type_aliases=sp_type_alias,
                )
            self.body += "  return ("

            for i, coeff in enumerate(coeffs):
                if coeff == 0:
                    continue
                offset_str: str = str(indices[i])

                # Build key since it's abs(coeff)
                sign = -1 if coeff <= 0 else 1
                decl_coeff = coeff * sign

                # ensure the correct sign is applied
                # in the upwind algorithm
                sign_str = "-" if coeff < 0 else "+"

                if i > 0:
                    self.body += "          "
                if offset_str == "0":
                    self.body += (
                        f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1,i2)]\n"
                    )

                else:
                    if dirn == 0:
                        self.body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0+{offset_str},i1,i2)]\n"
                    elif dirn == 1:
                        self.body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1+{offset_str},i2)]\n"
                    elif dirn == 2:
                        self.body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1,i2+{offset_str})]\n"

            self.body = self.body[:-1].replace("+-", "-") + f") * invdxx{dirn};\n }}\n"

        self.body += """}
    return 0.0 / 0.0;  // poison output if offset computed incorrectly
    """
        par.set_parval_from_str("fd_order", self.default_FDORDER)
        self.generate_CFunction()

    def generate_CFunction(self) -> None:
        """Generate CFunction from class parameters."""
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

    def __init__(self, CoordSystem: str, radiation_BC_fd_order: int = -1) -> None:
        self.CoordSystem = CoordSystem
        self.radiation_BC_fd_order = radiation_BC_fd_order
        self.include_CodeParameters_h = True
        self.includes: List[str] = []
        self.cfunc_decorators = ""
        self.CFunction: cfc.CFunction

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
        self.algorithm_header = """
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

    def regenerate_body(self) -> None:
        """Regenerate self.body based on class parameters."""
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
        """Generate CFunction from class parameters."""
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
    ) -> None:
        self.CoordSystem = CoordSystem
        self.radiation_BC_fd_order = radiation_BC_fd_order
        self.fp_type = fp_type

        self.cfunc_decorators = ""
        self.include_CodeParameters_h = True
        self.includes: List[str] = []
        self.prefunc = ""
        self.rfm = refmetric.reference_metric[self.CoordSystem]

        # These functions can be replaced to generate different prefunc strings
        self.upwind_setup_func = setup_Cfunction_FD1_arbitrary_upwind
        self.r_and_partial_xi_partial_r_derivs_prefunc_setup_func = (
            setup_Cfunction_r_and_partial_xi_partial_r_derivs
        )
        self.compute_partial_r_f_setup_func = setup_Cfunction_compute_partial_r_f

        # Initialize prefunc strings
        self.upwind_prefunc = ""
        self.r_and_partial_xi_partial_r_derivs_prefunc = ""
        self.compute_partial_r_f_prefunc = ""

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
        """Generate Upwind prefunctions from expressions and class parameters."""
        self.upwind_prefunc = ""

        for i in range(3):
            # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
            if not check_zero(self.rfm.Jac_dUrfm_dDSphUD[i][0]):
                self.upwind_prefunc += self.upwind_setup_func(
                    dirn=i,
                    radiation_BC_fd_order=self.radiation_BC_fd_order,
                    fp_type=self.fp_type,
                    rational_const_alias="static constexpr",
                ).CFunction.full_function

    def generate_r_and_partial_xi_partial_r_derivs_prefunc(self) -> None:
        """Generate additional derivative prefunctions from expressions and class parameters."""
        self.r_and_partial_xi_partial_r_derivs_prefunc = ""
        self.r_and_partial_xi_partial_r_derivs_prefunc += (
            self.r_and_partial_xi_partial_r_derivs_prefunc_setup_func(
                CoordSystem=self.CoordSystem,
                fp_type=self.fp_type,
            ).CFunction.full_function
        )

    def generate_compute_partial_r_f_prefunc(self) -> None:
        """Generate additional derivative prefunctions from expressions and class parameters."""
        self.compute_partial_r_f_prefunc = ""
        self.compute_partial_r_f_prefunc += self.compute_partial_r_f_setup_func(
            CoordSystem=self.CoordSystem,
            radiation_BC_fd_order=self.radiation_BC_fd_order,
        ).CFunction.full_function

    def generate_CFunction(self) -> None:
        """Generate C Function from class parameters."""
        self.generate_upwind_prefunc()
        self.generate_r_and_partial_xi_partial_r_derivs_prefunc()
        self.generate_compute_partial_r_f_prefunc()

        self.prefunc = self.upwind_prefunc
        self.prefunc += self.r_and_partial_xi_partial_r_derivs_prefunc
        self.prefunc += self.compute_partial_r_f_prefunc

        self.body = self.variable_defs
        self.body += self.function_calls
        self.body += self.algorithm_body

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
        self.body = ""

    def register(self) -> None:
        """Register CFunction."""
        _, actual_name = cfc.function_name_and_subdir_with_CoordSystem(
            os.path.join("."), self.name, self.CoordSystem
        )
        if not actual_name in cfc.CFunction_dict:
            cfc.register_CFunction(
                includes=self.includes,
                prefunc=self.prefunc,
                desc=self.desc,
                cfunc_type=self.cfunc_type,
                CoordSystem_for_wrapper_func=self.CoordSystem,
                name=self.name,
                params=self.params,
                include_CodeParameters_h=False,
                body=self.body,
            )


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

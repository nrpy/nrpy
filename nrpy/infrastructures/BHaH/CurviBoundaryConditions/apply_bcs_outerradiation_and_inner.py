# nrpy/infrastructures/BHaH/CurviBoundaryConditions/apply_bcs_outerradiation_and_inner.py
"""
Generates C code for applying radiation boundary conditions on outer boundaries.

This module constructs the C functions for applying Sommerfeld-type (radiation)
boundary conditions on the outer boundaries of a curvilinear grid. It generates
code to compute radial derivatives (∂f/∂r) using the chain rule, with spatial
derivatives calculated via arbitrary-order, upwinded finite differencing.

The main function, `register_CFunction_apply_bcs_outerradiation_and_inner`,
registers a top-level C routine that applies the outer radiation conditions and
then the separate inner boundary conditions.

This process is documented in the tutorial:
Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
         Samuel D. Tootle
         sdtootle **at** gmail **dot* com
"""

from typing import List, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends
import sympy.codegen.ast as sp_ast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin  # NRPy: Finite-difference module
import nrpy.grid as gri
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.params as par  # NRPy: Parameter interface
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.helpers.expression_utils import get_unique_expression_symbols_as_strings
from nrpy.helpers.parallelization.gpu_kernel import GPUKernel
from nrpy.validate_expressions.validate_expressions import check_zero


def _radiation_r_and_partials_exprs(CoordSystem: str) -> List[sp.Expr]:
    """
    Build r and partial_r x^i expressions for radiation BCs.

    For spheroidal GeneralRFM fisheye, use provider-specific closed forms
    instead of the generic inverse spherical Jacobian column.

    :param CoordSystem: Coordinate system used by the radiation boundary.
    :return: Physical radius followed by its three Cartesian-direction
        derivative factors.
    """
    rfm = refmetric.reference_metric[CoordSystem]
    provider = getattr(rfm, "general_rfm_provider", None)
    provider_kind = getattr(rfm, "general_rfm_provider_name", None)

    if provider is not None and provider_kind == "spheroidal_fisheye":
        raw_x, raw_y, raw_z = rfm.xx
        lam_xy, lam_z, _, _, _ = provider.axis_maps_for_inverse(provider.r)
        phys_r = sp.sqrt(lam_xy**2 * (raw_x**2 + raw_y**2) + lam_z**2 * raw_z**2)
        radial_vector = provider.cartesian_radial_vector_in_raw_coordinates()
        return [phys_r] + [component / phys_r for component in radial_vector]

    return [
        rfm.xxSph[0],
        rfm.Jac_dUrfm_dDSphUD[0][0],
        rfm.Jac_dUrfm_dDSphUD[1][0],
        rfm.Jac_dUrfm_dDSphUD[2][0],
    ]


###############################
## RADIATION (NewRad-like) BOUNDARY CONDITIONS.
##  Functions are fully documented in nrpytutorial's
##   Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb,
##   as well as below, in desc= and body=.
# r_and_partial_xi_partial_r_derivs(): Compute r(x0,x1,x2) and dx^i / dr
def setup_Cfunction_r_and_partial_xi_partial_r_derivs(
    CoordSystem: str,
    cfunc_decorators: str = "",
    enable_xx_to_Cart_precompute: bool = False,
) -> str:
    """
    Generate C code to compute the radial coordinate r(x0, x1, x2) and its derivatives.

    Compute the radial coordinate r(x0, x1, x2) and its partial derivatives
    partial x^i / partial r for a given coordinate system. This helper is called
    only from bcstruct_set_up(), which is host code in all parallelizations, so
    it unconditionally takes a params_struct.

    :param CoordSystem: The coordinate system for which to compute r and its derivatives.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param enable_xx_to_Cart_precompute: Evaluate geometry from cached Cartesian
        positions instead of solving the coordinate map at each boundary call.
    :return: A string containing the generated C code for the function.
    """
    desc = "Compute r(xx0,xx1,xx2) and partial_r x^i."
    cfunc_type = "static inline void"
    name = "r_and_partial_xi_partial_r_derivs"
    params = (
        "const params_struct *restrict params, const REAL xx0,const REAL xx1,const REAL xx2, REAL *r,"
        "REAL *partial_x0_partial_r,REAL *partial_x1_partial_r,REAL *partial_x2_partial_r"
    )
    body = ""

    rfm = refmetric.reference_metric[CoordSystem]
    if enable_xx_to_Cart_precompute:
        n = rfm.general_rfm_provider_meta["num_transitions"]
        return cfc.CFunction(
            includes=[],
            desc=desc,
            cfunc_type=cfunc_type,
            name=name,
            params="const params_struct *restrict params, const REAL xCart[3], REAL *r, REAL *partial_x0_partial_r, REAL *partial_x1_partial_r, REAL *partial_x2_partial_r",
            include_CodeParameters_h=False,
            body=f"""
const REAL x[3]={{xCart[0]-params->Cart_originx, xCart[1]-params->Cart_originy, xCart[2]-params->Cart_originz}};
REAL J[3][3], K[3][3], g[3][3], dg[3][3][3], ddg[3][3][3][3];
inverse_fisheye_geometry__N{n}(params, x, false, J, K, g, dg, ddg);
*r = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
*partial_x0_partial_r = (J[0][0]*x[0]+J[0][1]*x[1]+J[0][2]*x[2])/(*r);
*partial_x1_partial_r = (J[1][0]*x[0]+J[1][1]*x[1]+J[1][2]*x[2])/(*r);
*partial_x2_partial_r = (J[2][0]*x[0]+J[2][1]*x[1]+J[2][2]*x[2])/(*r);
""",
        ).full_function
    # sp.simplify(expr) is too slow here for SinhCylindrical, and for
    # spheroidal GeneralRFM fisheye the generic inverse spherical Jacobian
    # column is far more complicated than needed.
    expr_list = _radiation_r_and_partials_exprs(CoordSystem)
    unique_symbols = []
    for expr in expr_list:
        sub_list = get_unique_expression_symbols_as_strings(
            expr, exclude=[f"xx{i}" for i in range(3)]
        )
        unique_symbols += sub_list
    unique_symbols = sorted(list(set(unique_symbols)))

    for param_sym in unique_symbols:
        body += f"const REAL {param_sym} = params->{param_sym};\n"
    body += "\n"
    body += ccg.c_codegen(
        expr_list,
        [
            "*r",
            "*partial_x0_partial_r",
            "*partial_x1_partial_r",
            "*partial_x2_partial_r",
        ],
        verbose=False,
        include_braces=False,
    )

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cf.full_function


# partial_r f term: generate finite-difference coefficients
#   for \partial_i f with arbitrary upwinding:
def get_arb_offset_FD_coeffs_indices(
    FDORDER: int, offset: int, deriv: int
) -> Tuple[List[sp.Number], List[int]]:
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
def setup_Cfunction_FD1_arbitrary_upwind(
    dirn: int,
    cfunc_decorators: str = "",
    radiation_BC_fd_order: int = -1,
    rational_const_alias: str = "static const",
) -> str:
    """
    Set up the C function for computing the 1st derivative finite-difference.

    Supports arbitrary upwind for a given direction and order.

    :param dirn: Direction in which to compute the derivative.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param radiation_BC_fd_order: Finite difference order for radiation boundary condition.
                                  If -1, will use default finite difference order.
    :param rational_const_alias: Set constant alias for rational numbers, default is "static const"
    :return: The full C function as a string.
    """
    default_FDORDER = par.parval_from_str("fd_order")
    if radiation_BC_fd_order == -1:
        radiation_BC_fd_order = default_FDORDER

    par.set_parval_from_str("fd_order", radiation_BC_fd_order)
    parallelization = par.parval_from_str("parallelization")

    includes: List[str] = []
    desc = "Compute 1st derivative finite-difference derivative with arbitrary upwind"
    cfunc_type = "static inline REAL"
    name = f"FD1_arbitrary_upwind_x{dirn}_dirn"
    params = (
        """const size_t streamid, const REAL *restrict gf, const int i0, const int i1, const int i2, const int offset"""
        if parallelization == "cuda"
        else """const params_struct *restrict params,
const REAL *restrict gf, const int i0,const int i1,const int i2, const int offset"""
    )
    body = ""
    cfunc_decorators = "__device__" if parallelization == "cuda" else ""

    body = f"{parallel_utils.get_loop_parameters(parallelization)}\n"
    body += "switch(offset) {\n"

    tmp_list: List[int] = []
    fp_ccg_type = ccg.fp_type_to_sympy_type[par.parval_from_str("fp_type")]
    sp_type_alias = {sp_ast.real: fp_ccg_type}
    for offset in range(
        0, int(radiation_BC_fd_order // 2) + 1
    ):  # Use // for integer division
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += f"case {offset}:\n {{\n"
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
            RATIONAL_assignment = f"{rational_const_alias} REAL {str(v)}"
            body += sp.ccode(
                decl_coeff,
                assign_to=RATIONAL_assignment,
                type_aliases=sp_type_alias,
            )
        body += "  return ("

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
                body += "          "
            if offset_str == "0":
                body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0+{offset_str},i1,i2)]\n"
                elif dirn == 1:
                    body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1+{offset_str},i2)]\n"
                elif dirn == 2:
                    body += f"{sign_str} {rational_dict[decl_coeff]}*gf[IDX3(i0,i1,i2+{offset_str})]\n"

        body = body[:-1].replace("+-", "-") + (
            f") * invdxx{dirn};\n }} // END BLOCK: stencil offset {offset}\n"
        )

    body += """} // END SWITCH: finite-difference stencil offset
return 0.0 / 0.0;  // poison output if offset computed incorrectly
"""

    cf = cfc.CFunction(
        subdirectory="one_subdirectory_down",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        cfunc_decorators=cfunc_decorators,
    )

    par.set_parval_from_str("fd_order", default_FDORDER)

    return cf.full_function


# partial_r f term: Numerically evaluate partial_r f,
#   calling functions defined above.
def setup_Cfunction_compute_partial_r_f(
    CoordSystem: str,
    cfunc_decorators: str = "",
    radiation_BC_fd_order: int = -1,
) -> str:
    """
    Set up a C function for computing the partial derivative of f with respect to r.

    :param CoordSystem: Coordinate system to be used for the computation
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param radiation_BC_fd_order: Order of finite difference for radiation boundary conditions, default is -1
    :return: A C function for computing the partial derivative
    """
    desc = "Compute \\partial_r f"
    cfunc_type = "static inline REAL"
    name = "compute_partial_r_f"
    parallelization = par.parval_from_str("parallelization")
    params = f"""{("const size_t streamid" if parallelization == "cuda" else "const params_struct *restrict params")},
REAL *restrict xx[3], const REAL *restrict gfs,
const int which_gf, const int dest_i0,const int dest_i1,const int dest_i2,
const int FACEi0,const int FACEi1,const int FACEi2,
const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r"""
    partials_exprs = _radiation_r_and_partials_exprs(CoordSystem)

    if parallelization == "cuda" and "device" not in cfunc_decorators:
        cfunc_decorators += " __device__"

    default_FDORDER = par.parval_from_str("fd_order")
    if radiation_BC_fd_order == -1:
        radiation_BC_fd_order = default_FDORDER

    FD1_stencil_radius = int(radiation_BC_fd_order / 2)

    body = f"""  ///////////////////////////////////////////////////////////

  // FD1_stencil_radius = radiation_BC_fd_order/2 = {FD1_stencil_radius}
  const int FD1_stencil_radius = {FD1_stencil_radius};
"""
    body += f"{parallel_utils.get_loop_parameters(parallelization)}\n"
    body += """const int ntot = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

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
    for i in range(3):
        si = str(i)
        # The implicit spheroidal map has radial dependence in all three directions;
        # its numerical component calls cannot be evaluated by check_zero.
        if getattr(
            refmetric.reference_metric[CoordSystem], "general_rfm_provider_name", ""
        ) != "inverse_spheroidal_fisheye" and check_zero(
            partials_exprs[i + 1], fixed_mpfs_for_free_symbols=True
        ):
            body += f"  const REAL partial_x{si}_f=0.0;\n"
        else:
            body += (
                f"  int i{si}_offset = FACEi{si};  // Shift stencil away from the face we're updating.\n"
                f"  // Next adjust i{si}_offset so that FD stencil never goes out of bounds.\n"
                f"  if(dest_i{si} < FD1_stencil_radius) i{si}_offset = FD1_stencil_radius-dest_i{si};\n"
                f"  else if(dest_i{si} > (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1)) i{si}_offset = (Nxx_plus_2NGHOSTS{si}-FD1_stencil_radius-1) - dest_i{si};\n"
                f"  const REAL partial_x{si}_f=FD1_arbitrary_upwind_x{si}_dirn(params, &gfs[which_gf*ntot],dest_i0,dest_i1,dest_i2,i{si}_offset);\n"
            ).replace(
                "params,", "streamid, " if parallelization == "cuda" else "params,"
            )
    body += "  return partial_x0_partial_r*partial_x0_f + partial_x1_partial_r*partial_x1_f + partial_x2_partial_r*partial_x2_f;\n"

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=[],
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
        cfunc_decorators=cfunc_decorators,
    )
    return cf.full_function


# radiation_bcs(): Put it all together, for a single outer boundary point.
def setup_Cfunction_radiation_bcs(
    CoordSystem: str,
    cfunc_decorators: str = "",
    radiation_BC_fd_order: int = -1,
    rational_const_alias: str = "static const",
    rhs_is_scalar: bool = False,
    enable_xx_to_Cart_precompute: bool = False,
) -> str:
    """
    Generate C code to apply radiation boundary conditions in a given coordinate system.

    :param CoordSystem: The coordinate system to use.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param radiation_BC_fd_order: Finite differencing order to use. Default is -1.
    :param rational_const_alias: Alias for rational constants. Default is "static const".
    :param rhs_is_scalar: Receive the neighboring RHS as a scalar instead of reading an RHS array.
    :param enable_xx_to_Cart_precompute: Receive geometry computed once per
        boundary point from cached Cartesian positions.
    :return: A string containing the generated C code for the function.
    """
    includes: List[str] = []
    prefunc = ""
    parallelization = par.parval_from_str("parallelization")
    partials_exprs = _radiation_r_and_partials_exprs(CoordSystem)
    if parallelization == "cuda" and "device" not in cfunc_decorators:
        cfunc_decorators += " __device__"
    for i in range(3):
        # Do not generate FD1_arbitrary_upwind_xj_dirn() if the symbolic expression for dxj/dr == 0!
        # Keep all three radial contributions for the implicit spheroidal map.
        if getattr(
            refmetric.reference_metric[CoordSystem], "general_rfm_provider_name", ""
        ) == "inverse_spheroidal_fisheye" or not check_zero(
            partials_exprs[i + 1], fixed_mpfs_for_free_symbols=True
        ):
            prefunc += setup_Cfunction_FD1_arbitrary_upwind(
                dirn=i,
                cfunc_decorators=cfunc_decorators,
                radiation_BC_fd_order=radiation_BC_fd_order,
                rational_const_alias=rational_const_alias,
            )
    if enable_xx_to_Cart_precompute:
        prefunc += setup_Cfunction_r_and_partial_xi_partial_r_derivs(
            CoordSystem=CoordSystem,
            cfunc_decorators=cfunc_decorators,
            enable_xx_to_Cart_precompute=enable_xx_to_Cart_precompute,
        )
    prefunc += setup_Cfunction_compute_partial_r_f(
        CoordSystem=CoordSystem,
        cfunc_decorators=cfunc_decorators,
        radiation_BC_fd_order=radiation_BC_fd_order,
    )
    desc = r"""*** Apply radiation BCs to all outer boundaries. ***
"""
    cfunc_type = "static inline REAL"
    name = "radiation_bcs"
    params = (
        """const size_t streamid, REAL *restrict xx[3],
        const REAL *restrict gfs, REAL *restrict gfs_rhss,
        const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
        const int dest_i0,const int dest_i1,const int dest_i2,
        const short FACEi0,const short FACEi1,const short FACEi2,
        const REAL r, const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r,
        const REAL r_int, const REAL partial_x0_partial_r_int, const REAL partial_x1_partial_r_int, const REAL partial_x2_partial_r_int"""
        if parallelization == "cuda"
        else """const params_struct *restrict params,
    REAL *restrict xx[3],
    const REAL *restrict gfs, REAL *restrict gfs_rhss,
    const int which_gf, const REAL gf_wavespeed, const REAL gf_f_infinity,
    const int dest_i0,const int dest_i1,const int dest_i2,
    const short FACEi0,const short FACEi1,const short FACEi2,
    const REAL r, const REAL partial_x0_partial_r, const REAL partial_x1_partial_r, const REAL partial_x2_partial_r,
    const REAL r_int, const REAL partial_x0_partial_r_int, const REAL partial_x1_partial_r_int, const REAL partial_x2_partial_r_int"""
    )

    param_access = parallel_utils.get_params_access(parallelization)
    body = ""
    body += f"{parallel_utils.get_loop_parameters(parallelization)}\n"
    body += r"""// Nearest "interior" neighbor of this gridpoint, based on current face
const int dest_i0_int=dest_i0+1*FACEi0, dest_i1_int=dest_i1+1*FACEi1, dest_i2_int=dest_i2+1*FACEi2;
const REAL partial_r_f     = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0,    dest_i1,    dest_i2,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r    ,partial_x1_partial_r    ,partial_x2_partial_r);
const REAL partial_r_f_int = compute_partial_r_f(params,xx,gfs, which_gf,dest_i0_int,dest_i1_int,dest_i2_int,
                                                 FACEi0,FACEi1,FACEi2,
                                                 partial_x0_partial_r_int,partial_x1_partial_r_int,partial_x2_partial_r_int);

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
    body = body.replace(
        "params,", "streamid," if parallelization == "cuda" else "params,"
    )

    if rhs_is_scalar:
        params = params.replace("REAL *restrict gfs_rhss", "const REAL partial_t_f_int")
        body = body.replace(
            "const REAL partial_t_f_int = gfs_rhss[IDX4pt(which_gf, idx3_int)];", ""
        )

    cf = cfc.CFunction(
        subdirectory=CoordSystem,
        includes=includes,
        prefunc=prefunc.replace(
            "params->",
            param_access,
        ),
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body.replace(
            "params->",
            param_access,
        ),
        cfunc_decorators=cfunc_decorators,
    )
    return cf.full_function


def setup_Cfunction_apply_bcs_pure_only(
    enable_xx_to_Cart_precompute: bool = False,
    xx_to_Cart_in_griddata: bool = False,
) -> Tuple[str, str]:
    """
    Generate the prefunction string for apply_bcs_pure_only.

    This requires a function that will launch the compute kernel as well
    as the compute kernel itself.

    :param xx_to_Cart_in_griddata: Store cached Cartesian positions in separate griddata arrays.
    :param enable_xx_to_Cart_precompute: Read cached Cartesian positions and
        compute boundary geometry before the evolved-gridfunction loop.
    :return: A tuple containing the prefunction and the compute kernel.
    """
    # Header details for function that will launch the GPU kernel
    desc = "Apply BCs to pure boundary points"

    params_dict = {
        "params": "const params_struct *restrict",
        "bcstruct": "const bc_struct *restrict",
        "xx": "REAL *restrict *",
        "gfs": "REAL *restrict",
        "rhs_gfs": "REAL *restrict",
        "custom_wavespeed": "const REAL *",
        "custom_f_infinity": "const REAL *",
    }
    if enable_xx_to_Cart_precompute:
        params_dict["xx_to_Cart" if xx_to_Cart_in_griddata else "auxevol_gfs"] = (
            "REAL *restrict *" if xx_to_Cart_in_griddata else "const REAL *restrict"
        )
    name = "apply_bcs_pure_only"
    cfunc_type = "static void"
    parallelization = par.parval_from_str("parallelization")

    # Specify compute kernel body
    kernel_body = f"{parallel_utils.get_loop_parameters(parallelization)}\n"

    if parallelization == "cuda":
        kernel_body += """
for (int idx2d = tid0; idx2d < num_pure_outer_boundary_points; idx2d+=stride0) {
"""
    else:
        kernel_body += """
#pragma omp parallel for
    for (int idx2d = 0; idx2d < num_pure_outer_boundary_points; idx2d++) {
"""
    kernel_body += """const short i0 = pure_outer_bc_array[idx2d].i0;
    const short i1 = pure_outer_bc_array[idx2d].i1;
    const short i2 = pure_outer_bc_array[idx2d].i2;
    const short FACEX0 = pure_outer_bc_array[idx2d].FACEX0;
    const short FACEX1 = pure_outer_bc_array[idx2d].FACEX1;
    const short FACEX2 = pure_outer_bc_array[idx2d].FACEX2;
    const int idx3 = IDX3(i0,i1,i2);
    REAL* xx[3] = {x0, x1, x2};
"""
    if enable_xx_to_Cart_precompute:
        cart = ", ".join(
            (
                f"xCart{i}[IDX3(i0,i1,i2)]"
                if xx_to_Cart_in_griddata
                else gri.BHaHGridFunction.access_gf(
                    f"xCart{i}", 0, 0, 0, gf_array_name="auxevol_gfs"
                )
            )
            for i in range(3)
        )
        cart_int = (
            cart.replace("i0", "i0_int").replace("i1", "i1_int").replace("i2", "i2_int")
        )
        kernel_body += f"""
    const int i0_int=i0+FACEX0, i1_int=i1+FACEX1, i2_int=i2+FACEX2;
    const REAL xCart[3]={{{cart}}}, xCart_int[3]={{{cart_int}}};
    REAL r, partial_x0_partial_r, partial_x1_partial_r, partial_x2_partial_r;
    REAL r_int, partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int;
    r_and_partial_xi_partial_r_derivs(params, xCart, &r, &partial_x0_partial_r, &partial_x1_partial_r, &partial_x2_partial_r);
    r_and_partial_xi_partial_r_derivs(params, xCart_int, &r_int, &partial_x0_partial_r_int, &partial_x1_partial_r_int, &partial_x2_partial_r_int);
"""
    else:
        kernel_body += """    const REAL r = pure_outer_bc_array[idx2d].r;
    const REAL partial_x0_partial_r = pure_outer_bc_array[idx2d].partial_x0_partial_r;
    const REAL partial_x1_partial_r = pure_outer_bc_array[idx2d].partial_x1_partial_r;
    const REAL partial_x2_partial_r = pure_outer_bc_array[idx2d].partial_x2_partial_r;
    const REAL r_int = pure_outer_bc_array[idx2d].r_int;
    const REAL partial_x0_partial_r_int = pure_outer_bc_array[idx2d].partial_x0_partial_r_int;
    const REAL partial_x1_partial_r_int = pure_outer_bc_array[idx2d].partial_x1_partial_r_int;
    const REAL partial_x2_partial_r_int = pure_outer_bc_array[idx2d].partial_x2_partial_r_int;
"""
    kernel_body += """    for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
        // *** Apply radiation BCs to all outer boundary points. ***
        rhs_gfs[IDX4pt(which_gf, idx3)] = radiation_bcs(params, xx, gfs, rhs_gfs, which_gf,
                                                        custom_wavespeed[which_gf], custom_f_infinity[which_gf],
                                                        i0,i1,i2, FACEX0,FACEX1,FACEX2,
                                                        r, partial_x0_partial_r, partial_x1_partial_r, partial_x2_partial_r,
                                                        r_int, partial_x0_partial_r_int, partial_x1_partial_r_int, partial_x2_partial_r_int);
    }
  }
"""
    kernel_body = kernel_body.replace(
        "params,", "streamid," if parallelization == "cuda" else "params,"
    ).replace("custom_", "d_gridfunctions_" if parallelization == "cuda" else "custom_")
    comments = "Apply BCs to pure points."
    # Prepare the argument dicts
    arg_dict_cuda = {
        "num_pure_outer_boundary_points": "const int",
        "which_gz": "const int",
        "dirn": "const int",
        "pure_outer_bc_array": "const outerpt_bc_struct *restrict",
        "gfs": "REAL *restrict",
        "rhs_gfs": "REAL *restrict",
        "x0": "REAL *restrict",
        "x1": "REAL *restrict",
        "x2": "REAL *restrict",
    }
    arg_dict_host = {
        "params": "const params_struct *restrict",
        **arg_dict_cuda,
        "custom_wavespeed": "const REAL *",
        "custom_f_infinity": "const REAL *",
    }
    if enable_xx_to_Cart_precompute:
        if xx_to_Cart_in_griddata:
            for i in range(3):
                arg_dict_host[f"xCart{i}"] = "const REAL *restrict"
        else:
            arg_dict_host["auxevol_gfs"] = "const REAL *restrict"
    prefunc, new_body = parallel_utils.generate_kernel_and_launch_code(
        name,
        kernel_body,
        arg_dict_cuda,
        arg_dict_host,
        parallelization=parallelization,
        comments=comments,
        launch_dict={
            "blocks_per_grid": [
                "NRPYMAX(1U, (num_pure_outer_boundary_points + threads_in_x_dir -1) / threads_in_x_dir)"
            ],
            "stream": "params->grid_idx % NUM_STREAMS",
        },
        cfunc_type=cfunc_type,
        thread_tiling_macro_suffix="CURVIBC_RAD_PURE",
    )

    # Specify the function body for the launch kernel
    kernel_launch_body = f"""
const bc_info_struct *bc_info = &bcstruct->bc_info;
REAL *restrict x0 = xx[0];
REAL *restrict x1 = xx[1];
REAL *restrict x2 = xx[2];
  for (int which_gz = 0; which_gz < NGHOSTS; which_gz++) {{
    for (int dirn = 0; dirn < 3; dirn++) {{
      if (bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {{
        size_t gz_idx = dirn + (3 * which_gz);
        const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
        int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
        {new_body}
      }}
    }}
  }}
"""

    if enable_xx_to_Cart_precompute and xx_to_Cart_in_griddata:
        kernel_launch_body = (
            "\n".join(
                f"const REAL *restrict xCart{i} = xx_to_Cart[{i}];" for i in range(3)
            )
            + "\n"
            + kernel_launch_body
        )
    launch_kernel = GPUKernel(
        kernel_launch_body,
        params_dict,
        f"{name}",
        decorators="",
        comments=desc,
        cuda_check_error=False,
        streamid_param=False,
    )

    prefunc += launch_kernel.CFunction.full_function

    return prefunc, launch_kernel.c_function_call()


# apply_bcs_outerradiation_and_inner():
#   Apply radiation BCs at outer boundary points, and
#   inner boundary conditions at inner boundary points.
def register_CFunction_apply_bcs_outerradiation_and_inner(
    CoordSystem: str,
    cfunc_decorators: str = "",
    radiation_BC_fd_order: int = 2,
    rational_const_alias: str = "static const",
    enable_xx_to_Cart_precompute: bool = False,
    xx_to_Cart_in_griddata: bool = False,
) -> None:
    """
    Register a C function to apply boundary conditions to both pure outer and inner boundary points.

    :param CoordSystem: The coordinate system to use.
    :param cfunc_decorators: Optional decorators for CFunctions, e.g. CUDA identifiers, templates
    :param radiation_BC_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    :param rational_const_alias: Alias for rational constants. Default is "static const".
    :param xx_to_Cart_in_griddata: Store cached Cartesian positions in separate griddata arrays.
    :param enable_xx_to_Cart_precompute: Read fixed Cartesian positions from
        the selected storage. Requires inverse spheroidal fisheye and matching precompute.

    :raises ValueError: If Cartesian precompute is requested for another provider
        or for CUDA.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> from nrpy.reference_metric import unittest_CoordSystems
    >>> supported_Parallelizations = ["openmp", "cuda"]
    >>> name = "apply_bcs_outerradiation_and_inner__rfm"
    >>> for parallelization in supported_Parallelizations:
    ...    par.set_parval_from_str("parallelization", parallelization)
    ...    for CoordSystem in unittest_CoordSystems:
    ...       cfc.CFunction_dict.clear()
    ...       register_CFunction_apply_bcs_outerradiation_and_inner(CoordSystem)
    ...       generated_str = cfc.CFunction_dict[f'{name}__{CoordSystem}'].full_function
    ...       validation_desc = f"{name}__{parallelization}__{CoordSystem}"
    ...       validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    Setting up reference_metric[SinhSymTP]...
    Setting up reference_metric[HoleySinhSpherical]...
    Setting up reference_metric[Cartesian]...
    Setting up reference_metric[SinhCylindricalv2n2]...
    """
    if enable_xx_to_Cart_precompute and (
        not CoordSystem.startswith("GeneralRFM_inverse_spheroidal_fisheyeN")
        or par.parval_from_str("parallelization") == "cuda"
    ):
        raise ValueError(
            "Cached radiation geometry requires inverse spheroidal fisheye on the CPU."
        )
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    prefunc = setup_Cfunction_radiation_bcs(
        CoordSystem=CoordSystem,
        cfunc_decorators=cfunc_decorators,
        radiation_BC_fd_order=radiation_BC_fd_order,
        rational_const_alias=rational_const_alias,
        enable_xx_to_Cart_precompute=enable_xx_to_Cart_precompute,
    )
    apply_bcs_pure_only_prefuncs, apply_bcs_pure_only_function_call = (
        setup_Cfunction_apply_bcs_pure_only(
            enable_xx_to_Cart_precompute, xx_to_Cart_in_griddata
        )
    )
    prefunc += apply_bcs_pure_only_prefuncs
    desc = """This function is responsible for applying boundary conditions (BCs) to both pure outer and inner
boundary points. In the first step, it parallelizes the task using OpenMP and starts by applying BCs to
the outer boundary points layer-by-layer, prioritizing the faces in the order x0, x1, x2. The second step
applies BCs to the inner boundary points, which may map either to the grid interior or to the outer boundary.
"""
    cfunc_type = "void"
    name = "apply_bcs_outerradiation_and_inner"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
    const bc_struct *restrict bcstruct, REAL *restrict xx[3],
    const REAL custom_wavespeed[NUM_EVOL_GFS],
    const REAL custom_f_infinity[NUM_EVOL_GFS],
    REAL *restrict gfs, REAL *restrict rhs_gfs"""
    if enable_xx_to_Cart_precompute:
        if xx_to_Cart_in_griddata:
            params = params.replace(
                "REAL *restrict xx[3],",
                "REAL *restrict xx[3], REAL *restrict xx_to_Cart[3],",
            )
        else:
            params += ", const REAL *restrict auxevol_gfs"

    body: str = ""
    if par.parval_from_str("parallelization") == "cuda":
        body = r"""
  // Update device constants
  cudaMemcpyToSymbol(d_gridfunctions_wavespeed, custom_wavespeed, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_wavespeed failed");
  cudaMemcpyToSymbol(d_gridfunctions_f_infinity, custom_f_infinity, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_f_infinity failed");
  """
    body += rf"""
    ////////////////////////////////////////////////////////
  // STEP 1 of 2: Apply BCs to pure outer boundary points.
  //              By "pure" we mean that these points are
  //              on the outer boundary and not also on
  //              an inner boundary.
  //              Here we fill in the innermost ghost zone
  //              layer first and move outward. At each
  //              layer, we fill in +/- x0 faces first,
  //              then +/- x1 faces, finally +/- x2 faces,
  //              filling in the edges as we go.
  {apply_bcs_pure_only_function_call}

  ///////////////////////////////////////////////////////
  // STEP 2 of 2: Apply BCs to inner boundary points.
  //              These map to either the grid interior
  //              ("pure inner") or to pure outer boundary
  //              points ("inner maps to outer"). Those
  //              that map to outer require that outer be
  //              populated first; hence this being
  //              STEP 2 OF 2.
  apply_bcs_inner_only(commondata, params, bcstruct, rhs_gfs); // <- apply inner BCs to RHS gfs only
"""
    cfc.register_CFunction(
        includes=includes,
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )


def register_CFunction_apply_bcs_outerradiation_and_inner_accumulate(
    CoordSystem: str,
    cfunc_decorators: str = "",
    radiation_BC_fd_order: int = 2,
    rational_const_alias: str = "static const",
) -> None:
    """
    Register radiation and inner BCs for an RHS accumulated into two existing arrays.

    Call with ``prepare=1`` before the interior RHS updates, then ``prepare=0``
    after them. Between the calls, the interior kernel adds ``stage_dt * RHS``
    and ``accum_dt * RHS`` to the two output arrays. The projected ``in_gfs``
    must remain read-only, and all three arrays must be distinct.

    Preparation subtracts the old boundary-source contributions in reverse
    dependency order: inner boundaries first, then outer layers and faces in
    reverse. Completion adds the updated sources and radiation terms in the
    ordinary forward order. This preserves arbitrary original output values,
    including values that do not already obey the boundary relations, without
    allocating a complete RHS array. An output with a zero coefficient is unchanged.

    :param CoordSystem: The coordinate system to use.
    :param cfunc_decorators: Optional C function decorators.
    :param radiation_BC_fd_order: Finite differencing order for radiation boundary conditions.
    :param rational_const_alias: Alias for rational constants.

    Doctests:
    >>> from nrpy.helpers.generic import validate_strings
    >>> for parallelization in ("openmp", "cuda"):
    ...     par.set_parval_from_str("parallelization", parallelization)
    ...     cfc.CFunction_dict.clear()
    ...     register_CFunction_apply_bcs_outerradiation_and_inner_accumulate("Cartesian")
    ...     name = "apply_bcs_outerradiation_and_inner_accumulate__rfm__Cartesian"
    ...     generated_str = cfc.CFunction_dict[name].full_function
    ...     validation_desc = f"apply_bcs_outerradiation_and_inner_accumulate__{parallelization}__Cartesian"
    ...     validate_strings(generated_str, validation_desc, file_ext="cu" if parallelization == "cuda" else "c")
    """
    # Step 1: Declare the accumulated-boundary function.
    desc = "Apply radiation and inner RHS boundary conditions while accumulating into two RK arrays."
    cfunc_type = "void"
    name = "apply_bcs_outerradiation_and_inner_accumulate"
    params = """const commondata_struct *restrict commondata, const params_struct *restrict params,
const bc_struct *restrict bcstruct, REAL *restrict xx[3],
const REAL custom_wavespeed[NUM_EVOL_GFS], const REAL custom_f_infinity[NUM_EVOL_GFS],
const REAL *restrict in_gfs, REAL *restrict stage_gfs, REAL *restrict accum_gfs,
const REAL stage_dt, const REAL accum_dt, const int prepare"""
    parallelization = par.parval_from_str("parallelization")
    # Step 2: Generate the shared radiation formula and the two boundary kernels.
    prefunc = setup_Cfunction_radiation_bcs(
        CoordSystem=CoordSystem,
        cfunc_decorators=cfunc_decorators,
        radiation_BC_fd_order=radiation_BC_fd_order,
        rational_const_alias=rational_const_alias,
        rhs_is_scalar=True,
    )
    # Step 2.a: Inner sources are interior or outer points, never inner destinations.
    kernel_body = parallel_utils.get_loop_parameters(parallelization)
    kernel_body += (
        "\nfor (int pt = tid0; pt < num_inner_boundary_points; pt += stride0) {\n"
        if parallelization == "cuda"
        else "\n#pragma omp parallel for\nfor (int pt = 0; pt < num_inner_boundary_points; pt++) {\n"
    )
    parity_array = "d_evol_gf_parity" if parallelization == "cuda" else "evol_gf_parity"
    kernel_body += f"""
  const innerpt_bc_struct *restrict bc = &inner_bc_array[pt];
  const REAL sign = prepare ? -1.0 : 1.0;
  for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {{
    const int dst = IDX4pt(which_gf, bc->dstpt);
    const int src = IDX4pt(which_gf, bc->srcpt);
    const REAL weight = sign * bc->parity[{parity_array}[which_gf]];
    if (stage_dt != 0.0) stage_gfs[dst] += weight * stage_gfs[src];
    if (accum_dt != 0.0) accum_gfs[dst] += weight * accum_gfs[src];
  }} // END LOOP: which_gf over evolved fields
}} // END LOOP: pt over inner boundaries
"""
    arg_dict_cuda = {
        "num_inner_boundary_points": "const int",
        "inner_bc_array": "const innerpt_bc_struct *restrict",
        "stage_gfs": "REAL *restrict",
        "accum_gfs": "REAL *restrict",
        "stage_dt": "const REAL",
        "accum_dt": "const REAL",
        "prepare": "const int",
    }
    inner_prefunc, inner_call = parallel_utils.generate_kernel_and_launch_code(
        "apply_bcs_inner_accumulate",
        kernel_body,
        arg_dict_cuda,
        {"params": "const params_struct *restrict", **arg_dict_cuda},
        parallelization=parallelization,
        comments="Prepare or complete inner-boundary RHS accumulation in two RK arrays.",
        launch_dict={
            "blocks_per_grid": [
                "NRPYMAX(1U, (num_inner_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir)"
            ],
            "stream": "params->grid_idx % NUM_STREAMS",
        },
        thread_tiling_macro_suffix="CURVIBC_INNER",
    )

    # Step 2.b: At a fixed input stage, q_ghost = weight * q_inner + beta.
    # Passing a zero neighboring RHS to radiation_bcs evaluates beta.
    kernel_body = parallel_utils.get_loop_parameters(parallelization)
    kernel_body += (
        "\nfor (int pt = tid0; pt < num_pure_outer_boundary_points; pt += stride0) {\n"
        if parallelization == "cuda"
        else "\n#pragma omp parallel for\nfor (int pt = 0; pt < num_pure_outer_boundary_points; pt++) {\n"
    )
    kernel_body += """
  const outerpt_bc_struct *restrict bc = &pure_outer_bc_array[pt];
  const int dstpt = IDX3(bc->i0, bc->i1, bc->i2);
  const int srcpt = IDX3(bc->i0 + bc->FACEX0, bc->i1 + bc->FACEX1, bc->i2 + bc->FACEX2);
  const REAL rinv = 1.0 / bc->r;
  const REAL weight = bc->r_int * bc->r_int * bc->r_int * rinv * rinv * rinv;
  const REAL sign = prepare ? -1.0 : 1.0;
  REAL *xx[3] = {x0, x1, x2};
  for (int which_gf = 0; which_gf < NUM_EVOL_GFS; which_gf++) {
    const int dst = IDX4pt(which_gf, dstpt);
    const int src = IDX4pt(which_gf, srcpt);
    if (stage_dt != 0.0) stage_gfs[dst] += sign * weight * stage_gfs[src];
    if (accum_dt != 0.0) accum_gfs[dst] += sign * weight * accum_gfs[src];
    if (!prepare) {
      const REAL beta = radiation_bcs(params, xx, in_gfs, 0.0, which_gf,
          custom_wavespeed[which_gf], custom_f_infinity[which_gf],
          bc->i0, bc->i1, bc->i2, bc->FACEX0, bc->FACEX1, bc->FACEX2,
          bc->r, bc->partial_x0_partial_r, bc->partial_x1_partial_r, bc->partial_x2_partial_r,
          bc->r_int, bc->partial_x0_partial_r_int, bc->partial_x1_partial_r_int, bc->partial_x2_partial_r_int);
      if (stage_dt != 0.0) stage_gfs[dst] += stage_dt * beta;
      if (accum_dt != 0.0) accum_gfs[dst] += accum_dt * beta;
    } // END IF: complete radiation increments
  } // END LOOP: which_gf over evolved fields
} // END LOOP: pt over outer boundaries
"""
    if parallelization == "cuda":
        kernel_body = kernel_body.replace("params,", "streamid,").replace(
            "custom_", "d_gridfunctions_"
        )
    arg_dict_cuda = {
        "num_pure_outer_boundary_points": "const int",
        "pure_outer_bc_array": "const outerpt_bc_struct *restrict",
        "in_gfs": "const REAL *restrict",
        "stage_gfs": "REAL *restrict",
        "accum_gfs": "REAL *restrict",
        "x0": "REAL *restrict",
        "x1": "REAL *restrict",
        "x2": "REAL *restrict",
        "stage_dt": "const REAL",
        "accum_dt": "const REAL",
        "prepare": "const int",
    }
    outer_prefunc, outer_call = parallel_utils.generate_kernel_and_launch_code(
        "apply_bcs_pure_accumulate",
        kernel_body,
        arg_dict_cuda,
        {
            "params": "const params_struct *restrict",
            **arg_dict_cuda,
            "custom_wavespeed": "const REAL *",
            "custom_f_infinity": "const REAL *",
        },
        parallelization=parallelization,
        comments="Prepare or complete radiation-boundary RHS accumulation in two RK arrays.",
        launch_dict={
            "blocks_per_grid": [
                "NRPYMAX(1U, (num_pure_outer_boundary_points + threads_in_x_dir - 1) / threads_in_x_dir)"
            ],
            "stream": "params->grid_idx % NUM_STREAMS",
        },
        thread_tiling_macro_suffix="CURVIBC_RAD_PURE",
    )
    # Step 3: Subtract sources in reverse order and restore them in forward order.
    prefunc += inner_prefunc + outer_prefunc
    body = ""
    if parallelization == "cuda":
        body += """
if (!prepare) {
  cudaMemcpyToSymbol(d_gridfunctions_wavespeed, custom_wavespeed, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_wavespeed failed");
  cudaMemcpyToSymbol(d_gridfunctions_f_infinity, custom_f_infinity, NUM_EVOL_GFS * sizeof(REAL));
  cudaCheckErrors(copy, "Copy to d_gridfunctions_f_infinity failed");
} // END IF: copy radiation parameters
"""
    body += f"""
const bc_info_struct *bc_info = &bcstruct->bc_info;
const int num_inner_boundary_points = bc_info->num_inner_boundary_points;
const innerpt_bc_struct *restrict inner_bc_array = bcstruct->inner_bc_array;
REAL *restrict x0 = xx[0];
REAL *restrict x1 = xx[1];
REAL *restrict x2 = xx[2];

// Inner sources may be pure outer points: subtract them before changing outer values.
if (prepare) {{
  {inner_call}
}} // END IF: subtract inner source values

// Reverse the complete layer/face dependency order only during preparation.
for (int layer = 0; layer < NGHOSTS; layer++) {{
  const int which_gz = prepare ? NGHOSTS - 1 - layer : layer;
  for (int face = 0; face < 3; face++) {{
    const int dirn = prepare ? 2 - face : face;
    const int num_pure_outer_boundary_points = bc_info->num_pure_outer_boundary_points[which_gz][dirn];
    if (num_pure_outer_boundary_points > 0) {{
      const int gz_idx = dirn + 3 * which_gz;
      const outerpt_bc_struct *restrict pure_outer_bc_array = bcstruct->pure_outer_bc_array[gz_idx];
      {outer_call}
    }} // END IF: nonempty outer boundary face
  }} // END LOOP: face over boundary directions
}} // END LOOP: layer over ghost zones

// All outer sources must have their final accumulated values before inner completion.
if (!prepare) {{
  {inner_call}
}} // END IF: restore inner source values
"""
    cfc.register_CFunction(
        includes=["BHaH_defines.h", "BHaH_function_prototypes.h"],
        prefunc=prefunc,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
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

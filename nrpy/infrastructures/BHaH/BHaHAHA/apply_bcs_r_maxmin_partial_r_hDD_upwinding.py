"""
Module providing functions for setting up Curvilinear boundary conditions.

This is documented in Tutorial-Start_to_Finish-Curvilinear_BCs.ipynb.

Authors: Zachariah B. Etienne
         zachetie **at** gmail **dot* com
         Terrence Pierre Jacques
"""

# Step P1: Import needed NRPy core modules:
from typing import List, Tuple

import sympy as sp  # SymPy: The Python computer algebra package upon which NRPy depends
import sympy.codegen.ast as sp_ast

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.finite_difference as fin  # NRPy: Finite-difference module
import nrpy.params as par  # NRPy: Parameter interface


###############################
# Compute partial derivatives of h_{ij} r ~ r_max face
###############################
# partial_r f term: generate finite-difference coefficients
#   for \partial_i f with arbitrary upwinding:
def get_arb_offset_FD_coeffs_indices(
    FDORDER: int, offset: int, deriv: int
) -> Tuple[List[float], List[int]]:
    """
    Generate finite-difference coefficients and their corresponding indices for an arbitrary upwinding.

    :param FDORDER: Order of the finite difference stencil.
    :param offset: Offset for upwinding. Positive for forward, negative for backward.
    :param deriv: Order of the derivative (1 for first derivative, etc.).
    :return: A tuple containing two lists:
             - List of finite-difference coefficients.
             - List of corresponding indices for each coefficient.

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
    upwinding_fd_order: int = -1,
) -> str:
    """
    Set up the C function for computing the first derivative using finite-difference with arbitrary upwinding.

    :param dirn: Direction (0, 1, or 2) in which the derivative is to be computed.
    :param upwinding_fd_order: Finite difference order for radiation boundary condition.
                               If -1, uses the default finite difference order.
    :return: Full C function code as a string.
    """
    default_FDORDER = par.parval_from_str("fd_order")
    if upwinding_fd_order == -1:
        upwinding_fd_order = default_FDORDER

    par.set_parval_from_str("fd_order", upwinding_fd_order)

    includes: List[str] = []
    desc = """Computes the first derivative of a grid function along the x0 direction using finite-difference
schemes with arbitrary upwind offsets. The offset determines which stencil to use:
  - Negative offsets correspond to backward (upwind) stencils.
  - Zero offset corresponds to a centered stencil.
  - Positive offsets correspond to forward (downwind) stencils.
Returns the computed first derivative at the given grid point."""
    cfunc_type = "static inline REAL"
    name = f"FD1_arbitrary_upwind_x{dirn}_dirn"
    params = """const commondata_struct *restrict commondata,
                const REAL *restrict gf, const int i0, const int i1, const int i2, const int offset"""
    body = """
  const MAYBE_UNUSED int Nxx_plus_2NGHOSTS0 = commondata->bcstruct_Nxx_plus_2NGHOSTS0;
  const MAYBE_UNUSED int Nxx_plus_2NGHOSTS1 = commondata->bcstruct_Nxx_plus_2NGHOSTS1;
  const MAYBE_UNUSED int Nxx_plus_2NGHOSTS2 = commondata->bcstruct_Nxx_plus_2NGHOSTS2;
  const REAL invdxx0 = 1.0 / (commondata->bcstruct_dxx0);
  switch(offset) {
"""

    tmp_list: List[int] = []
    fp_type = par.parval_from_str("fp_type")
    fp_ccg_type = ccg.fp_type_to_sympy_type[fp_type]
    sp_type_alias = {sp_ast.real: fp_ccg_type}
    for offset in range(
        0, int(upwinding_fd_order // 2) + 1
    ):  # Use // for integer division
        tmp_list.append(offset)
        if offset > 0:
            tmp_list.append(-offset)

    for offset in tmp_list:
        body += f"case {offset}:\n"
        body += "  return ("
        coeffs, indices = get_arb_offset_FD_coeffs_indices(
            upwinding_fd_order, offset, 1
        )

        for i, coeff in enumerate(coeffs):
            if coeff == 0:
                continue
            offset_str: str = str(indices[i])
            if i > 0:
                body += "          "
            if offset_str == "0":
                body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2)]\n"
            else:
                if dirn == 0:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0+{offset_str},i1,i2)]\n"
                elif dirn == 1:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1+{offset_str},i2)]\n"
                elif dirn == 2:
                    body += f"+{sp.ccode(coeff, type_aliases=sp_type_alias)}*gf[IDX3(i0,i1,i2+{offset_str})]\n"

        body = body[:-1].replace("+-", "-") + f") * invdxx{dirn};\n"

    body += """
  default:
    // Return NaN if offset is invalid
    return 0.0 / 0.0;
}
"""
    par.set_parval_from_str("fd_order", default_FDORDER)

    cf = cfc.CFunction(
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return cf.full_function


# apply_bcs_outer_partial_r_hDD_upwinding_and_inner():
#   Apply radiation BCs at outer boundary points, and
#   inner boundary conditions at inner boundary points.
def register_CFunction_apply_bcs_r_maxmin_partial_r_hDD_upwinding(
    upwinding_fd_order: int = 2,
) -> None:
    """
    Register a C function to apply boundary conditions at both pure outer and inner boundary points.

    :param upwinding_fd_order: Finite differencing order for the radiation boundary conditions. Default is 2.
    """
    includes = ["BHaH_defines.h"]
    prefunc = setup_Cfunction_FD1_arbitrary_upwind(
        dirn=0,
        upwinding_fd_order=upwinding_fd_order,
    )
    desc = """Applies boundary conditions to r_max and possibly r_min (outer) boundaries of the computational grid by computing
the radial derivative partial_r hDD using upwind finite-difference stencils. The function iterates over the boundary points
in the radial (x0) direction and computes the derivative for each component of hDD, ensuring that the chosen stencils do not
access grid points outside the computational domain. OpenMP parallelization is employed to optimize computations over
the angular directions (x1 and x2).

@param commondata - Pointer to common data structure containing boundary condition and grid information.
@param xx - Array of pointers to grid coordinate arrays.
@param gfs - Pointer to the grid functions array where derivatives are stored.
@param fill_r_min_ghosts - Boolean flag indicating if r_min boundary ghost zones should be filled.
@return - Void.
@note - Parallelizes angular computations to enhance performance and reduce computation time.
"""
    cfunc_type = "void"
    name = "apply_bcs_r_maxmin_partial_r_hDD_upwinding"
    params = """const commondata_struct *restrict commondata, REAL *restrict xx[3], REAL *restrict gfs,
                const bool fill_r_min_ghosts"""
    body = r"""
  const int Nxx_plus_2NGHOSTS0 = commondata->bcstruct_Nxx_plus_2NGHOSTS0;
  const int Nxx_plus_2NGHOSTS1 = commondata->bcstruct_Nxx_plus_2NGHOSTS1;
  const int Nxx_plus_2NGHOSTS2 = commondata->bcstruct_Nxx_plus_2NGHOSTS2;
  const int Nxxtot012 = Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2;

  //////////////////////////////////////////////////////////////
  // Evaluate partial_r hDD at r = r_max boundary

  // Iterate over the r_max boundary points in the radial (x0) direction.
  for (int i0 = Nxx_plus_2NGHOSTS0 - NGHOSTS; i0 < Nxx_plus_2NGHOSTS0; i0++) {
    // i0 = Nxx_plus_2NGHOSTS0 - 1 --> offset = -NGHOSTS
    // i0 = Nxx_plus_2NGHOSTS0 - 2 --> offset = -NGHOSTS + 1
    // -> i0 = Nxx_plus_2NGHOSTS0 - j --> offset = -NGHOSTS + (j-1)
    // -> j = Nxx_plus_2NGHOSTS0 - i0 --> offset = -NGHOSTS + ((Nxx_plus_2NGHOSTS0-i0) - 1)
    const int offset = -NGHOSTS + ((Nxx_plus_2NGHOSTS0 - i0) - 1);

    // Parallelize computations over the angular (x1 and x2) directions to leverage multi-core processing.
#pragma omp parallel for
    for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
      for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
        for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
          int base_gf = -1;

          // Map the interpolation source grid function to its corresponding base grid function.
          switch (which_gf) {
          case SRC_PARTIAL_D_HDD000GF:
            base_gf = SRC_HDD00GF;
            break;
          case SRC_PARTIAL_D_HDD001GF:
            base_gf = SRC_HDD01GF;
            break;
          case SRC_PARTIAL_D_HDD002GF:
            base_gf = SRC_HDD02GF;
            break;
          case SRC_PARTIAL_D_HDD011GF:
            base_gf = SRC_HDD11GF;
            break;
          case SRC_PARTIAL_D_HDD012GF:
            base_gf = SRC_HDD12GF;
            break;
          case SRC_PARTIAL_D_HDD022GF:
            base_gf = SRC_HDD22GF;
            break;
          case SRC_PARTIAL_D_WW0GF:
            base_gf = SRC_WWGF;
            break;
          default:
            // Skip processing for undefined grid function indices to maintain data integrity.
            break;
          } // END SWITCH to set base_gf

          if (base_gf != -1) {
            // Compute the radial derivative using the appropriate upwind stencil based on the offset.
            const REAL partial_x0_f = FD1_arbitrary_upwind_x0_dirn(commondata, &gfs[base_gf * Nxxtot012], i0, i1, i2, offset);
            // Store the computed derivative in the target grid function array.
            gfs[IDX4(which_gf, i0, i1, i2)] = partial_x0_f;
          } // END IF the derivative gridfunction needs to be set
        } // END LOOP over gridfunctions
      } // END LOOP over i1
    } // END LOOP over i2
  } // END LOOP over i0: iterating through r_max radial boundary points

  /////////////////////////////////////////////////////////////////
  // Evaluate partial_r hDD at r = r_min boundary if required
  if (fill_r_min_ghosts) {
    // Iterate over the r_min boundary points in the radial (x0) direction.
    for (int i0 = 0; i0 < NGHOSTS; i0++) {
      const int offset = NGHOSTS - i0;

      // Parallelize computations over the angular (x1 and x2) directions to leverage multi-core processing.
#pragma omp parallel for
      for (int i2 = NGHOSTS; i2 < Nxx_plus_2NGHOSTS2 - NGHOSTS; i2++) {
        for (int i1 = NGHOSTS; i1 < Nxx_plus_2NGHOSTS1 - NGHOSTS; i1++) {
          for (int which_gf = 0; which_gf < NUM_INTERP_SRC_GFS; which_gf++) {
            int base_gf = -1;

            // Map the interpolation source grid function to its corresponding base grid function.
            switch (which_gf) {
            case SRC_PARTIAL_D_HDD000GF:
              base_gf = SRC_HDD00GF;
              break;
            case SRC_PARTIAL_D_HDD001GF:
              base_gf = SRC_HDD01GF;
              break;
            case SRC_PARTIAL_D_HDD002GF:
              base_gf = SRC_HDD02GF;
              break;
            case SRC_PARTIAL_D_HDD011GF:
              base_gf = SRC_HDD11GF;
              break;
            case SRC_PARTIAL_D_HDD012GF:
              base_gf = SRC_HDD12GF;
              break;
            case SRC_PARTIAL_D_HDD022GF:
              base_gf = SRC_HDD22GF;
              break;
            case SRC_PARTIAL_D_WW0GF:
              base_gf = SRC_WWGF;
              break;
            default:
              // Skip processing for undefined grid function indices to maintain data integrity.
              break;
            } // END SWITCH to set base_gf

            if (base_gf != -1) {
              // Compute the radial derivative using the appropriate upwind stencil based on the offset.
              const REAL partial_x0_f = FD1_arbitrary_upwind_x0_dirn(commondata, &gfs[base_gf * Nxxtot012], i0, i1, i2, offset);
              // Store the computed derivative in the target grid function array.
              gfs[IDX4(which_gf, i0, i1, i2)] = partial_x0_f;
            } // END IF the derivative gridfunction needs to be set
          } // END LOOP over gridfunctions
        } // END LOOP over i1
      } // END LOOP over i2
    } // END LOOP: over i0, iterating through r_min radial boundary points
  } // END IF: fill_r_min_ghosts flag check
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

"""
C function to apply copy boundary conditions to scalar primitives and outflow boundary conditions to fluid velocities.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.basis_transforms.jacobians as bt
import nrpy.equations.grhd.Min_Max_and_Piecewise_Expressions as noif
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.BSSN_to_ADM import BSSN_to_ADM
from nrpy.validate_expressions.validate_expressions import check_zero


def _assert_zero(expr: sp.Expr) -> None:
    """
    Assert that the given symbolic expression evaluates to zero.

    :param expr: Expression to validate.
    :raises AssertionError: If the expression is not numerically zero.
    """
    if not check_zero(expr):
        raise AssertionError(f"Expression did not evaluate to zero: {expr}")


def _groovy_simplify(coord_system: str, expr: sp.Expr) -> sp.Expr:
    """
    Apply the symbolic simplification strategy used by GRoovy.

    :param coord_system: Coordinate system name.
    :param expr: Expression to simplify.
    :return: Simplified expression.
    """
    if "Sinh" in coord_system:
        return sp.together(expr)
    return sp.simplify(expr)


def register_CFunction_apply_copy_and_outflow_bcs(
    CoordSystem: str,
    enable_GoldenKernels: bool = False,
    evolving_temperature: bool = False,
    evolving_spacetime: bool = True,
    evolving_neutrinos: bool = False,
    evolving_entropy: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the outer-boundary update routine for primitive hydrodynamic data.

    Primitives are copied from the first interior point, and the
    velocity is projected to prevent inflow through the outer boundary. Inner
    boundary points are then populated from the precomputed parity map.

    :param CoordSystem: The coordinate system.
    :param enable_GoldenKernels: Boolean to enable Golden Kernels.
    :param evolving_temperature: Whether temperature is an evolved primitive.
    :param evolving_spacetime: Whether the spacetime metric is evolved in `in_gfs`.
    :param evolving_neutrinos: Whether NRPyLeakage variables are active.
    :param evolving_entropy: Whether entropy is an evolved primitive.
    :return: None if in registration phase, else the updated NRPy environment.
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Function skeleton
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    cfunc_type = "void"
    name = "apply_copy_and_outflow_bcs"
    desc = (
        "Apply copy boundary conditions to scalar primitives and outflow boundary "
        "conditions to the velocity."
    )
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "const ghl_parameters *restrict ghl_params, "
        "const bc_struct *restrict bcstruct, "
        "REAL *restrict xx[3], "
        "REAL *restrict in_gfs, "
        "REAL *restrict auxevol_gfs"
    )

    # Step 3: Initialize ADM quantities, basis transforms, and the reference metric.
    AitoB = BSSN_to_ADM(CoordSystem=CoordSystem)
    basis_transforms = bt.basis_transforms[CoordSystem]
    rfm = refmetric.reference_metric[CoordSystem]

    # Step 4: Construct the outward radial unit vector in the rfm basis.
    theta = rfm.xxSph[1]
    phi = rfm.xxSph[2]

    rhatU_Cart = ixp.zerorank1()
    rhatU_Cart[0] = sp.sin(theta) * sp.cos(phi)
    rhatU_Cart[1] = sp.sin(theta) * sp.sin(phi)
    rhatU_Cart[2] = sp.cos(theta)

    rhatU = basis_transforms.basis_transform_vectorU_from_Cartesian_to_rfmbasis(
        rhatU_Cart
    )

    normalization = cast(sp.Expr, sp.sympify(0))
    for i in range(3):
        for j in range(3):
            normalization += _groovy_simplify(
                CoordSystem, AitoB.gammaDD[i][j] * rhatU[i] * rhatU[j]
            )

    for i in range(3):
        rhatU[i] /= sp.sqrt(normalization)

    check_normalization = cast(sp.Expr, sp.sympify(0))
    for i in range(3):
        for j in range(3):
            check_normalization += _groovy_simplify(
                CoordSystem, AitoB.gammaDD[i][j] * rhatU[i] * rhatU[j]
            )
    _assert_zero(check_normalization - sp.sympify(1))

    # Step 5: Project the velocity along the radial direction and remove inflow.
    rescaledvU = ixp.declarerank1("rescaledvU")

    VU = ixp.zerorank1()
    for i in range(3):
        VU[i] = rescaledvU[i] * rfm.ReU[i]

    vr: sp.Expr = cast(sp.Expr, sp.sympify(0))
    for i in range(3):
        for j in range(3):
            vr += _groovy_simplify(CoordSystem, AitoB.gammaDD[i][j] * VU[i] * rhatU[j])

    vr = _groovy_simplify(CoordSystem, vr)
    outflow_check = noif.coord_less_bound(vr, sp.sympify(0))

    new_VU = ixp.zerorank1()
    new_rescaledvU = ixp.zerorank1()
    for i in range(3):
        new_VU[i] = VU[i] - outflow_check * vr * rhatU[i]
        if "Spherical" in CoordSystem:
            new_rescaledvU[i] = _groovy_simplify(CoordSystem, new_VU[i] / rfm.ReU[i])
        else:
            new_rescaledvU[i] = sp.together(new_VU[i] / rfm.ReU[i])

    # Step 6: Generate symbolic C code for the velocity update.
    output_vars = [new_rescaledvU[0], new_rescaledvU[1], new_rescaledvU[2]]
    vars_grid_access = [
        gri.BHaHGridFunction.access_gf("rescaledvU0", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU1", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU2", gf_array_name="auxevol_gfs"),
    ]
    expr_body = ccg.c_codegen(
        output_vars,
        vars_grid_access,
        enable_simd=False,
        enable_fd_codegen=True,
        enable_GoldenKernels=enable_GoldenKernels,
    )

    # Step 7: Build the boundary-condition C routine.
    metric_gf_array_name = "in_gfs" if evolving_spacetime else "auxevol_gfs"
    body = r"""

  // Unpack bc_info from bcstruct.
  const bc_info_struct *bc_info = &bcstruct->bc_info;
"""
    if evolving_temperature:
        if evolving_entropy:
            body += r"""
  const int NUM_PRIM_GFs = 8;
  const int prims_gfs[8] = {
      RHOBGF, PGF, RESCALEDVU0GF, RESCALEDVU1GF,
      RESCALEDVU2GF, YEGF, TEMPERATUREGF, SGF};
"""
        else:
            body += r"""
  const int NUM_PRIM_GFs = 7;
  const int prims_gfs[7] = {
      RHOBGF, PGF, RESCALEDVU0GF, RESCALEDVU1GF,
      RESCALEDVU2GF, YEGF, TEMPERATUREGF};
"""
    else:
        if evolving_entropy:
            body += r"""
  const int NUM_PRIM_GFs = 6;
  const int prims_gfs[6] = {
      RHOBGF, PGF, RESCALEDVU0GF, RESCALEDVU1GF, RESCALEDVU2GF, SGF};
"""
        else:
            body += r"""
  const int NUM_PRIM_GFs = 5;
  const int prims_gfs[5] = {
      RHOBGF, PGF, RESCALEDVU0GF, RESCALEDVU1GF, RESCALEDVU2GF};
"""
    if evolving_neutrinos:
        body += r"""

  const int NUM_LEAKAGE_GFs = 12;
  const int leakage_gfs[12] = {
      TAU_0_NUEGF,    TAU_1_NUEGF,
      TAU_0_ANUEGF,   TAU_1_ANUEGF,
      TAU_0_NUXGF,    TAU_1_NUXGF,
      KAPPA_0_NUEGF,  KAPPA_1_NUEGF,
      KAPPA_0_ANUEGF, KAPPA_1_ANUEGF,
      KAPPA_0_NUXGF,  KAPPA_1_NUXGF};
"""
    body += r"""

  ////////////////////////////////////////////////////////
  // Step 1 of 2: Apply outer-boundary conditions first.
  //              We sweep from the innermost ghost zone
  //              outward so edge and corner values are
  //              available as soon as they are needed.
  #pragma omp parallel
  {
    for(int which_gz=0;which_gz<NGHOSTS;which_gz++) for(int dirn=0;dirn<3;dirn++) {
      if(bc_info->num_pure_outer_boundary_points[which_gz][dirn] > 0) {
#pragma omp for
        for(int idx2d=0;
            idx2d<bc_info->num_pure_outer_boundary_points[which_gz][dirn];
            idx2d++) {
          const short i0 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i0;
          const short i1 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i1;
          const short i2 = bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].i2;
          const short FACEX0 =
              bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX0;
          const short FACEX1 =
              bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX1;
          const short FACEX2 =
              bcstruct->pure_outer_bc_array[dirn + (3*which_gz)][idx2d].FACEX2;
          const int idx = IDX3(i0,i1,i2);
          const int idx_offset1 = IDX3(i0+1*FACEX0,i1+1*FACEX1,i2+1*FACEX2);
          for(int which_gf=0;which_gf<NUM_PRIM_GFs;which_gf++) {
            auxevol_gfs[IDX4pt(prims_gfs[which_gf], idx)] =
                auxevol_gfs[IDX4pt(prims_gfs[which_gf], idx_offset1)];
          }
"""
    if evolving_neutrinos:
        body += r"""
          // Set leakage variables to zero on the outer boundary.
          for(int which_gf=0;which_gf<NUM_LEAKAGE_GFs;which_gf++) {
            auxevol_gfs[IDX4pt(leakage_gfs[which_gf], idx)] = 0.0;
          }
"""
    body += r"""

          // Update all velocity components together before recomputing u^t.
"""
    body += r"""
        const REAL xx0 = xx[0][i0];
        const REAL xx1 = xx[1][i1];
        const REAL xx2 = xx[2][i2];
"""
    body += expr_body + r"""

          compute_up_index_velocity_time_component_pointwise(
              commondata, params, &commondata->ghl_params,
              __METRIC_GFS__[IDX4pt(ALPHAGF, idx)],
              __METRIC_GFS__[IDX4pt(VETU0GF, idx)],
              __METRIC_GFS__[IDX4pt(VETU1GF, idx)],
              __METRIC_GFS__[IDX4pt(VETU2GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD00GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD01GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD02GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD11GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD12GF, idx)],
              __METRIC_GFS__[IDX4pt(HDD22GF, idx)],
              __METRIC_GFS__[IDX4pt(CFGF, idx)],
              &auxevol_gfs[IDX4pt(RESCALEDVU0GF, idx)],
              &auxevol_gfs[IDX4pt(RESCALEDVU1GF, idx)],
              &auxevol_gfs[IDX4pt(RESCALEDVU2GF, idx)],
              &auxevol_gfs[IDX4pt(U4UTGF, idx)]);
        }
      }
    }
  }

  ///////////////////////////////////////////////////////
  // Step 2 of 2: Populate inner-boundary points from
  //              their mapped source points once the
  //              outer boundary data are available.
""".replace("__METRIC_GFS__", metric_gf_array_name)

    if evolving_neutrinos:
        body += r"""
  // collapse(2) improves throughput here, especially in 2D.
#pragma omp parallel for collapse(2)
  for(int which_gf=0;which_gf<NUM_PRIM_GFs;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      auxevol_gfs[IDX4pt(prims_gfs[which_gf], dstpt)] =
          bcstruct->inner_bc_array[pt].parity[auxevol_gf_parity[prims_gfs[which_gf]]] *
          auxevol_gfs[IDX4pt(prims_gfs[which_gf], srcpt)];
    }
  }

#pragma omp parallel for collapse(2)
  for(int which_gf=0;which_gf<NUM_LEAKAGE_GFs;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      auxevol_gfs[IDX4pt(leakage_gfs[which_gf], dstpt)] =
          bcstruct->inner_bc_array[pt].parity[
              auxevol_gf_parity[leakage_gfs[which_gf]]] *
          auxevol_gfs[IDX4pt(leakage_gfs[which_gf], srcpt)];
    }
  }
"""
    else:
        body += r"""
  // collapse(2) improves throughput here, especially in 2D.
#pragma omp parallel for collapse(2)
  for(int which_gf=0;which_gf<NUM_PRIM_GFs;which_gf++) {
    for(int pt=0;pt<bc_info->num_inner_boundary_points;pt++) {
      const int dstpt = bcstruct->inner_bc_array[pt].dstpt;
      const int srcpt = bcstruct->inner_bc_array[pt].srcpt;
      auxevol_gfs[IDX4pt(prims_gfs[which_gf], dstpt)] =
          bcstruct->inner_bc_array[pt].parity[auxevol_gf_parity[prims_gfs[which_gf]]] *
          auxevol_gfs[IDX4pt(prims_gfs[which_gf], srcpt)];
    }
  }
"""

    # Step 8: Register the final C function.
    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
        name=name,
        params=params,
        body=body,
        enable_simd=False,
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

"""
C function registration for GRoovy shock-test initial data.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.basis_transforms.jacobians as bt
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.grhd.ShockTests import GRHD_ShockTests


def register_CFunction_shock_tests_initial_data(
    IDType: str,
    grhayl_setup_str: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    OMP_collapse: int,
    fp_type: str = "double",
    enable_GoldenKernels: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register a routine that fills GRoovy gridfunctions with shock-test initial data.

    The generated C code seeds the primitive variables symbolically, converts the
    flat Cartesian metric to the native GRoovy basis, and then applies GRHayL's
    primitive limiting and conservative-variable reconstruction.

    :param IDType: Shock-test initial data type.
    :param grhayl_setup_str: C code used to initialize the GRHayL EOS structs.
    :param CoordSystem: The coordinate system.
    :param enable_rfm_precompute: Whether to enable reference-metric
        precomputation.
    :param OMP_collapse: Level of OpenMP loop collapsing.
    :param fp_type: Floating-point type used in generated C code.
    :param enable_GoldenKernels: Whether to enable Golden Kernel generation.
    :return: None if in registration phase, else the updated NRPy environment.

    Doctests:
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> supported_Parallelizations = ["openmp"]
    >>> name = "initial_data"
    >>> for parallelization in supported_Parallelizations:
    ...     par.set_parval_from_str("parallelization", parallelization)
    ...     cfc.CFunction_dict.clear()
    ...     _ = register_CFunction_shock_tests_initial_data(
    ...         IDType="balsara0",
    ...         grhayl_setup_str="",
    ...         CoordSystem="Spherical",
    ...         enable_rfm_precompute=False,
    ...         OMP_collapse=1,
    ...     )
    ...     generated_str = clang_format(cfc.CFunction_dict[name].full_function)
    ...     validation_desc = (
    ...         f"shock_tests_initial_data__{parallelization}__Spherical__balsara0"
    ...     )
    ...     _ = validate_strings(generated_str, validation_desc, file_ext="c")
    Setting up reference_metric[Spherical]...
    Setting up basis_transforms[Spherical]...
    """
    # Step 1: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 2: Set up the C-function metadata.
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    desc = "Initial data for GRoovy shock tests"
    cfunc_type = "void"
    name = "initial_data"
    params = (
        "commondata_struct *restrict commondata, griddata_struct *restrict griddata"
    )

    # Step 3: Set up the background reference metric and shifted shock-test profile.
    rfm_key = CoordSystem + ("_rfm_precompute" if enable_rfm_precompute else "")
    rfm = refmetric.reference_metric[rfm_key]
    x = rfm.xx_to_Cart[0]
    y = rfm.xx_to_Cart[1]
    z = rfm.xx_to_Cart[2]

    offset = sp.Rational(11, 10)
    shifted_center = ixp.zerorank1()
    shifted_center[0] = x
    shifted_center[1] = y - offset
    shifted_center[2] = z - offset * sp.Rational(3, 10)

    rp = cast(
        sp.Expr,
        rfm.Cart_to_xx[0]
        .subs(rfm.Cartx, shifted_center[0])
        .subs(rfm.Carty, shifted_center[1])
        .subs(rfm.Cartz, shifted_center[2]),
    )

    shock_test_data = GRHD_ShockTests(
        IDType=IDType,
        x=z,
        r=rp,
        bound=sp.Integer(1),
        theta=rfm.xxSph[1],
    )

    # Step 4: Assemble the primitive and metric fields to seed on the grid.
    rescaledvU = ixp.zerorank1()
    output_vars: List[sp.Expr] = [
        shock_test_data.rho,
        shock_test_data.press,
        cast(sp.Expr, rescaledvU[0]),
        cast(sp.Expr, rescaledvU[1]),
        cast(sp.Expr, rescaledvU[2]),
    ]
    vars_grid_access: List[str] = [
        gri.BHaHGridFunction.access_gf("rhob", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("P", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU0", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU1", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("rescaledvU2", gf_array_name="auxevol_gfs"),
    ]

    gammaDD_Cart = ixp.zerorank2()
    gammaDD_Cart[0][0] = sp.Integer(1)
    gammaDD_Cart[1][1] = sp.Integer(1)
    gammaDD_Cart[2][2] = sp.Integer(1)

    KDD_Cart = ixp.zerorank2()
    betaU_Cart = ixp.zerorank1()
    alpha = sp.Integer(1)

    basis_transforms = bt.basis_transforms[CoordSystem]
    gammaDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        gammaDD_Cart
    )
    adm_to_bssn = ADM_to_BSSN(
        gammaDD=gammaDD,
        KDD=KDD_Cart,
        betaU=betaU_Cart,
        BU=betaU_Cart,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
    )

    metric_rhs: List[sp.Expr] = [
        adm_to_bssn.hDD[0][0],
        adm_to_bssn.hDD[0][1],
        adm_to_bssn.hDD[0][2],
        adm_to_bssn.hDD[1][1],
        adm_to_bssn.hDD[1][2],
        adm_to_bssn.hDD[2][2],
        adm_to_bssn.aDD[0][0],
        adm_to_bssn.aDD[0][1],
        adm_to_bssn.aDD[0][2],
        adm_to_bssn.aDD[1][1],
        adm_to_bssn.aDD[1][2],
        adm_to_bssn.aDD[2][2],
        cast(sp.Expr, betaU_Cart[0]),
        cast(sp.Expr, betaU_Cart[1]),
        cast(sp.Expr, betaU_Cart[2]),
        alpha,
        adm_to_bssn.cf,
        adm_to_bssn.trK,
    ]
    metric_lhs = [
        gri.BHaHGridFunction.access_gf("hDD00", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("hDD01", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("hDD02", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("hDD11", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("hDD12", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("hDD22", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD00", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD01", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD02", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD11", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD12", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("aDD22", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("vetU0", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("vetU1", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("vetU2", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("alpha", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("cf", gf_array_name="auxevol_gfs"),
        gri.BHaHGridFunction.access_gf("trK", gf_array_name="auxevol_gfs"),
    ]
    output_vars += metric_rhs
    vars_grid_access += metric_lhs

    # Step 5: Build the function body.
    body = grhayl_setup_str + r"""for(int grid=0; grid<commondata->NUMGRIDS; grid++) {
  // Unpack the current grid and its gridfunctions.
  params_struct *restrict params = &griddata[grid].params;
#include "set_CodeParameters.h"
  REAL *restrict xx[3];
  for (int ww = 0; ww < 3; ww++) {
    xx[ww] = griddata[grid].xx[ww];
  } // END LOOP: for ww over coordinate directions
  REAL *restrict in_gfs = griddata[grid].gridfuncs.y_n_gfs;
  REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;
"""
    if enable_rfm_precompute:
        body += r"""
  rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;
"""

    apply_grhayl_limits_code = r"""
      ghl_primitive_quantities prims;
      ghl_metric_quantities metric;
      ghl_conservative_quantities cons;

      // Rebuild GRHayL structs from the freshly seeded gridfunction data.
      basis_transform_rfm_basis_to_Cartesian(
          commondata, params, &prims, &cons, &metric, i0, i1, i2, xx,
          auxevol_gfs, in_gfs);

      bool speed_limited;
      ghl_enforce_primitive_limits_and_compute_u0(
          &commondata->ghl_params, &commondata->eos, &metric, &prims,
          &speed_limited);

      // Store the limited primitive variables and reconstructed conservatives.
      basis_transform_Cartesian_to_rfm_basis(
          commondata, params, &prims, &cons, i0, i1, i2, xx,
          auxevol_gfs, in_gfs);

      auxevol_gfs[IDX4(U4UTGF, i0, i1, i2)] = prims.u0;
"""

    # Step 6: Generate the pointwise loop that seeds the shock profile.
    loop_body = ccg.c_codegen(
        output_vars,
        vars_grid_access,
        enable_simd=False,
        fp_type=fp_type,
        enable_GoldenKernels=enable_GoldenKernels,
    )
    body += lp.simple_loop(
        loop_body=loop_body + apply_grhayl_limits_code,
        loop_region="all points",
        enable_intrinsics=False,
        CoordSystem=CoordSystem,
        enable_rfm_precompute=enable_rfm_precompute,
        read_xxs=not enable_rfm_precompute,
        OMP_collapse=OMP_collapse,
    )
    body += "} // END LOOP: for grid over grids\n"

    # Step 7: Register the generated C function.
    cfc.register_CFunction(
        include_CodeParameters_h=False,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        body=body,
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

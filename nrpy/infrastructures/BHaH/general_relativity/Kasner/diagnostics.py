"""
Kasner-specific diagnostics helpers for BHaH general relativity infrastructures.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.helpers.parallelization.utilities as parallel_utils
import nrpy.indexedexp as ixp
import nrpy.infrastructures.BHaH.simple_loop as lp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.basis_transforms.jacobians import BasisTransforms
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    kasner_adm_quantities,
)
from nrpy.helpers.expression_utils import (
    generate_definition_header,
    get_params_commondata_symbols_from_expr_list,
)


def register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics: bool = False,
    enable_psi4: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register ``diagnostic_gfs_set`` and inject Kasner-specific diagnostic channels.

    :param enable_interp_diagnostics: Whether interpolation diagnostics are enabled.
    :param enable_psi4: Whether Psi4 diagnostics are enabled.
    :return: Parallel-codegen registration token or ``None`` during registration phase.
    :raises ValueError: If ``EvolvedConformalFactor_cf`` is not one of ``phi``, ``chi``, or ``W``.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = [
        "BHaH_defines.h",
        "BHaH_function_prototypes.h",
        "diagnostics/diagnostic_gfs.h",
    ]
    desc = """
 * @file diagnostic_gfs_set.c
 * @brief Populate per-grid diagnostic arrays used by interpolation and integration routines.
 *
 * The function "diagnostic_gfs_set" loops over all grids and fills per-grid diagnostic arrays:
 *   1) Compute a residual-type diagnostic at all points using a helper that evaluates the
 *      finite-difference residual.
 *   2) If enabled at code-generation time, apply inner boundary conditions to that residual by
 *      copying from a source point to a destination point with a sign determined by the relevant
 *      parity, ensuring parity-consistent values near symmetry or excision boundaries.
 *   3) Copy selected evolved gridfunctions from the current time level (y_n_gfs) into designated
 *      diagnostic channels for downstream consumers.
 *   4) Set additional Kasner exact and recovered diagnostic channels.
 *   5) Set a per-point grid identifier channel to the grid index (converted to REAL).
 *
 * The routine assumes each per-grid output buffer is contiguous and large enough to store all
 * diagnostic channels:
 *     TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2
 * Loops over grid points may be parallelized with OpenMP if available.
 *
 * @param[in]  commondata
 *     Pointer to global simulation metadata (e.g., counters and configuration) accessed by helpers
 *     and used to determine the number of grids to process.
 * @param[in]  griddata
 *     Pointer to an array of per-grid data. For each grid, this provides parameters, coordinates,
 *     boundary condition metadata, and gridfunctions (including y_n_gfs and any auxiliary data)
 *     referenced by this routine and its helpers.
 * @param[out] diagnostic_gfs
 *     Array of per-grid output buffers. For each grid, diagnostic_gfs[grid] must point to a buffer
 *     of size TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2.
 *
 * @return void.
"""
    cfunc_type = "void"
    name = "diagnostic_gfs_set"
    params = "const commondata_struct *restrict commondata, const griddata_struct *restrict griddata, REAL *restrict diagnostic_gfs[MAXNUMGRIDS]"

    # Register every DIAG gridfunction this generator writes before building parity tables.
    kasner_diag_gridfunctions: Tuple[Tuple[str, str], ...] = (
        ("DIAG_ADDXX", "ADD00"),
        ("DIAG_ADDXY", "ADD01"),
        ("DIAG_ADDXZ", "ADD02"),
        ("DIAG_ADDYY", "ADD11"),
        ("DIAG_ADDYZ", "ADD12"),
        ("DIAG_ADDZZ", "ADD22"),
        ("DIAG_TRK", "trK"),
        ("DIAG_HDDXX", "HDD00"),
        ("DIAG_HDDXY", "HDD01"),
        ("DIAG_HDDXZ", "HDD02"),
        ("DIAG_HDDYY", "HDD11"),
        ("DIAG_HDDYZ", "HDD12"),
        ("DIAG_HDDZZ", "HDD22"),
        ("DIAG_EXACT_ADDXX", "EXACT_ADD00"),
        ("DIAG_EXACT_ADDXY", "EXACT_ADD01"),
        ("DIAG_EXACT_ADDXZ", "EXACT_ADD02"),
        ("DIAG_EXACT_ADDYY", "EXACT_ADD11"),
        ("DIAG_EXACT_ADDYZ", "EXACT_ADD12"),
        ("DIAG_EXACT_ADDZZ", "EXACT_ADD22"),
        ("DIAG_EXACT_W", "EXACT_W"),
        ("DIAG_EXACT_TRK", "EXACT_trK"),
        ("DIAG_EXACT_HDDXX", "EXACT_HDD00"),
        ("DIAG_EXACT_HDDXY", "EXACT_HDD01"),
        ("DIAG_EXACT_HDDXZ", "EXACT_HDD02"),
        ("DIAG_EXACT_HDDYY", "EXACT_HDD11"),
        ("DIAG_EXACT_HDDYZ", "EXACT_HDD12"),
        ("DIAG_EXACT_HDDZZ", "EXACT_HDD22"),
        ("DIAG_PX_RECOVERED", "p1_recovered"),
        ("DIAG_PY_RECOVERED", "p2_recovered"),
        ("DIAG_PZ_RECOVERED", "p3_recovered"),
        ("DIAG_EXACT_PX", "EXACT_p1"),
        ("DIAG_EXACT_PY", "EXACT_p2"),
        ("DIAG_EXACT_PZ", "EXACT_p3"),
    )
    sym01_components: Tuple[Tuple[int, int, str], ...] = (
        (0, 0, "XX"),
        (0, 1, "XY"),
        (0, 2, "XZ"),
        (1, 1, "YY"),
        (1, 2, "YZ"),
        (2, 2, "ZZ"),
    )
    for gf_name, gf_desc in kasner_diag_gridfunctions:
        gri.register_gridfunctions(names=[gf_name], desc_list=[gf_desc], group="DIAG")
    gri.register_gridfunctions(
        names="DIAG_HAMILTONIAN", desc="H_constraint", group="DIAG"
    )
    gri.register_gridfunctions(names="DIAG_MSQUARED", desc="M^2", group="DIAG")
    gri.register_gridfunctions(names="DIAG_LAPSE", desc="Lapse", group="DIAG")
    gri.register_gridfunctions(names="DIAG_W", desc="Conformal_factor_W", group="DIAG")
    gri.register_gridfunctions(names="DIAG_GRIDINDEX", desc="GridIndex", group="DIAG")
    if enable_psi4:
        gri.register_gridfunctions(names="DIAG_PSI4_RE", desc="Psi4_Re", group="DIAG")
        gri.register_gridfunctions(names="DIAG_PSI4_IM", desc="Psi4_Im", group="DIAG")
    gri.register_gridfunctions_for_single_rank2(
        "DIAG_RBARDD",
        desc="Ricci_tensor_component_RbarDD",
        symmetry="sym01",
        dimension=3,
        group="DIAG",
    )
    diag_gf_parity_types = gri.BHaHGridFunction.set_parity_types(
        sorted([v.name for v in gri.glb_gridfcs_dict.values() if v.group == "DIAG"])
    )

    body = f"MAYBE_UNUSED const int8_t diag_gf_parities[{len(diag_gf_parity_types)}] = {{ {', '.join(map(str, diag_gf_parity_types))} }};\n"
    body += """  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    const rfm_struct *restrict rfmstruct = griddata[grid].rfmstruct;
    SET_NXX_PLUS_2NGHOSTS_VARS(grid);
    const REAL *restrict y_n_gfs = griddata[grid].gridfuncs.y_n_gfs;
    MAYBE_UNUSED REAL *restrict auxevol_gfs = griddata[grid].gridfuncs.auxevol_gfs;

    // Poison diagnostic_gfs (for debugging purposes only; WARNING: this might make valgrind ineffective)
    // #pragma omp parallel for
    //     for (int ii = 0; ii < TOTAL_NUM_DIAG_GFS * Nxx_plus_2NGHOSTS0 * Nxx_plus_2NGHOSTS1 * Nxx_plus_2NGHOSTS2; ii++) {
    //       diagnostic_gfs[grid][ii] = NAN;
    //     } // END LOOP over all points & gridfunctions, poisoning diagnostic_gfs

    // Set Ricci and constraints gridfunctions
"""
    if par.parval_from_str("parallelization") == "cuda":
        body += """    Ricci_eval_host(params, rfmstruct, y_n_gfs, diagnostic_gfs[grid]);
"""
    else:
        body += """    Ricci_eval(params, rfmstruct, y_n_gfs, auxevol_gfs);
"""
    body += """    constraints_eval(commondata, params, rfmstruct, y_n_gfs, auxevol_gfs, diagnostic_gfs[grid]);
"""
    if enable_psi4:
        body += """
    // NOTE: Inner boundary conditions must be set before any interpolations are performed, whether for psi4 decomp. or interp diags.
    // Set psi4 gridfunctions
    psi4(commondata, params, (REAL * restrict*)griddata[grid].xx, y_n_gfs, auxevol_gfs, diagnostic_gfs[grid]);
    const int inner_bc_apply_gfs[] = {DIAG_PSI4_REGF, DIAG_PSI4_IMGF};
    const int num_inner_bc_apply_gfs = (int)(sizeof(inner_bc_apply_gfs) / sizeof(inner_bc_apply_gfs[0]));
    apply_bcs_inner_only_specific_gfs(commondata, params, &griddata[grid].bcstruct, diagnostic_gfs[grid], num_inner_bc_apply_gfs, diag_gf_parities,
                                      inner_bc_apply_gfs);
"""
    if enable_interp_diagnostics:
        body += """
    {
      // NOTE: Inner boundary conditions must be set before any interpolations are performed, whether for psi4 decomp. or interp diags.
      // Apply inner bcs to constraints needed to do interpolation correctly
      const int inner_bc_apply_gfs[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};
      const int num_inner_bc_apply_gfs = (int)(sizeof(inner_bc_apply_gfs) / sizeof(inner_bc_apply_gfs[0]));
      apply_bcs_inner_only_specific_gfs(commondata, params, &griddata[grid].bcstruct, diagnostic_gfs[grid], num_inner_bc_apply_gfs, diag_gf_parities,
                                        inner_bc_apply_gfs);
    } // END set inner BCs on desired GFs
"""
    # Build the Kasner symbolic kernel immediately before emitting it so the registration flow
    # reads in the same order as the generated C function body.
    CoordSystem = par.parval_from_str("CoordSystem_to_register_CodeParameters")
    _ = BSSN_quantities[CoordSystem]
    if "KASNER_t0" not in par.glb_code_params_dict:
        par.register_CodeParameters(
            "REAL",
            "nrpy.equations.general_relativity.InitialData_Cartesian",
            ["KASNER_t0", "KASNER_p1", "KASNER_p2", "KASNER_p3"],
            [1.0, -1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0],
            commondata=True,
        )
    hDD = ixp.declarerank2("hDD", symmetry="sym01")
    aDD = ixp.declarerank2("aDD", symmetry="sym01")
    cf, trK = sp.symbols("cf trK", real=True)
    KASNER_t_phys = sp.Symbol("KASNER_t_phys", real=True)
    KASNER_p1 = sp.Symbol("KASNER_p1", real=True)
    KASNER_p2 = sp.Symbol("KASNER_p2", real=True)
    KASNER_p3 = sp.Symbol("KASNER_p3", real=True)

    # First construct the exact Kasner BSSN fields in the target coordinate system.
    exact_gammaDD, exact_KDD, _alpha, exact_betaU, exact_BU = kasner_adm_quantities(
        KASNER_t_phys, KASNER_p1, KASNER_p2, KASNER_p3
    )
    exact_adm2bssn = ADM_to_BSSN(
        gammaDD=exact_gammaDD,
        KDD=exact_KDD,
        betaU=exact_betaU,
        BU=exact_BU,
        CoordSystem="Cartesian",
    )
    basis_transforms = BasisTransforms(CoordSystem)
    exact_gammabarDD = (
        basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
            exact_adm2bssn.gammabarDD
        )
    )
    exact_AbarDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        exact_adm2bssn.AbarDD
    )
    rfm = refmetric.reference_metric[CoordSystem]
    exact_hDD = ixp.zerorank2()
    exact_aDD = ixp.zerorank2()
    gammabarDD = ixp.zerorank2()
    AbarDD = ixp.zerorank2()
    gammaDD = ixp.zerorank2()
    KDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            exact_hDD[i][j] = (exact_gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][
                j
            ]
            exact_aDD[i][j] = exact_AbarDD[i][j] / rfm.ReDD[i][j]
    cf_evolution = par.parval_from_str("EvolvedConformalFactor_cf")
    if cf_evolution == "phi":
        exp4phi = sp.exp(4 * cf)
    elif cf_evolution == "chi":
        exp4phi = 1 / cf
    elif cf_evolution == "W":
        exp4phi = 1 / cf**2
    else:
        raise ValueError(f"Unsupported EvolvedConformalFactor_cf = {cf_evolution}")
    # Then rebuild the numerical BSSN tensors from the evolved gridfunctions so we can compare
    # numerical and exact quantities within one symbolic codegen pass.
    for i in range(3):
        for j in range(3):
            gammabarDD[i][j] = hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]
            gammaDD[i][j] = exp4phi * gammabarDD[i][j]
            KDD[i][j] = exp4phi * AbarDD[i][j] + sp.Rational(1, 3) * gammaDD[i][j] * trK
    gammaDD_cart = basis_transforms.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(
        gammaDD
    )
    KDD_cart = basis_transforms.basis_transform_tensorDD_from_rfmbasis_to_Cartesian(KDD)
    gammaUU_cart, _ = ixp.symm_matrix_inverter3x3(gammaDD_cart)
    KUD_cart = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            for k in range(3):
                KUD_cart[i][j] += gammaUU_cart[i][k] * KDD_cart[k][j]

    # Assemble one ordered expression/LHS list so c_codegen emits all Kasner diagnostics,
    # including copied numerical fields, exact fields, and recovered exponents, in one pass.
    expr_list: List[sp.Expr] = []
    lhs_list: List[str] = []
    for i, j, suffix in sym01_components:
        expr_list.append(aDD[i][j])
        lhs_list.append(f"diagnostic_gfs[grid][IDX4(DIAG_ADD{suffix}GF, i0, i1, i2)]")
    expr_list.append(trK)
    lhs_list.append("diagnostic_gfs[grid][IDX4(DIAG_TRKGF, i0, i1, i2)]")
    for i, j, suffix in sym01_components:
        expr_list.append(hDD[i][j])
        lhs_list.append(f"diagnostic_gfs[grid][IDX4(DIAG_HDD{suffix}GF, i0, i1, i2)]")
    for i, j, suffix in sym01_components:
        expr_list.append(exact_aDD[i][j])
        lhs_list.append(
            f"diagnostic_gfs[grid][IDX4(DIAG_EXACT_ADD{suffix}GF, i0, i1, i2)]"
        )
    expr_list.append(cast(sp.Expr, exact_adm2bssn.cf))
    lhs_list.append("diagnostic_gfs[grid][IDX4(DIAG_EXACT_WGF, i0, i1, i2)]")
    expr_list.append(cast(sp.Expr, exact_adm2bssn.trK))
    lhs_list.append("diagnostic_gfs[grid][IDX4(DIAG_EXACT_TRKGF, i0, i1, i2)]")
    for i, j, suffix in sym01_components:
        expr_list.append(exact_hDD[i][j])
        lhs_list.append(
            f"diagnostic_gfs[grid][IDX4(DIAG_EXACT_HDD{suffix}GF, i0, i1, i2)]"
        )
    for i, axis in enumerate(("X", "Y", "Z")):
        expr_list.append(-KASNER_t_phys * KUD_cart[i][i])
        lhs_list.append(
            f"diagnostic_gfs[grid][IDX4(DIAG_P{axis}_RECOVEREDGF, i0, i1, i2)]"
        )
    for exact_p, axis in zip((KASNER_p1, KASNER_p2, KASNER_p3), ("X", "Y", "Z")):
        expr_list.append(exact_p)
        lhs_list.append(f"diagnostic_gfs[grid][IDX4(DIAG_EXACT_P{axis}GF, i0, i1, i2)]")
    param_symbols, commondata_symbols = get_params_commondata_symbols_from_expr_list(
        expr_list, exclude=["xx0", "xx1", "xx2", "cf", "trK", "KASNER_t_phys"]
    )
    params_definitions = generate_definition_header(
        param_symbols,
        var_access=parallel_utils.get_params_access("openmp"),
    )
    commondata_definitions = generate_definition_header(
        commondata_symbols,
        var_access=parallel_utils.get_commondata_access("openmp"),
    )
    loop_body = f"""{params_definitions}
{commondata_definitions}
{ccg.c_codegen(
        expr_list,
        lhs_list,
        automatically_read_gf_data_from_memory=True,
        include_braces=False,
        verbose=False,
    )}"""
    body += f"""
    const REAL *restrict xx[3] = {{griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2]}};
    const REAL *restrict in_gfs = y_n_gfs;
    const REAL KASNER_t_phys = commondata->KASNER_t0 + commondata->time;
{lp.simple_loop(
        loop_body=loop_body,
        loop_region="all points",
        read_xxs=True,
    )
}
    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {{
      const int idx3 = IDX3(i0, i1, i2);
      diagnostic_gfs[grid][IDX4pt(DIAG_LAPSEGF, idx3)] = y_n_gfs[IDX4pt(ALPHAGF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_WGF, idx3)] = y_n_gfs[IDX4pt(CFGF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_GRIDINDEXGF, idx3)] = (REAL)grid;
    }} // END LOOP over all gridpoints to set lapse/W diagnostics
  }} // END LOOP over grids
"""

    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()


def register_CFunction_diagnostics_nearest() -> Union[None, pcg.NRPyEnv_type]:
    """
    Construct and register a Kasner-specific nearest-sampled diagnostics dispatcher.

    :return: Parallel-codegen registration token or ``None`` during registration phase.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h", "diagnostic_gfs.h"]
    desc = """
 * @brief Dispatch Kasner nearest-sampled diagnostics by invoking specialized helper routines for 0D, 1D, and 2D outputs.
 *
 * The diagnostics_nearest() dispatcher coordinates sampling of caller-provided diagnostic gridfunction data and
 * delegates output to three helper functions:
 *   - diagnostics_nearest_grid_center(): emits 0D diagnostics at the index triplet nearest the physical center.
 *   - diagnostics_nearest_1d_y_and_z_axes(): emits 1D diagnostics along lines nearest the y and z axes.
 *   - diagnostics_nearest_2d_xy_and_yz_planes(): emits 2D diagnostics across planes nearest the xy and yz planes.
 *
 * A single USER-EDIT block appears before the per-grid loop. In this Kasner version, the default selections include
 * the standard constraint diagnostics together with the Kasner exact and recovered channels. Users may still edit the
 * generated C arrays directly.
 *
 * @param[in] commondata
 *   Pointer to shared simulation metadata and runtime context, including NUMGRIDS and iteration/time information.
 *
 * @param[in] griddata
 *   Pointer to an array of per-grid data structures. For grid index "grid", griddata[grid] provides parameters,
 *   coordinates, and strides required by the diagnostics helper routines.
 *
 * @param[in] gridfuncs_diags
 *   Array of length MAXNUMGRIDS. For each grid index "grid", gridfuncs_diags[grid] must point to caller-owned
 *   REAL diagnostic gridfunction data that serve as the sampling source.
 *
 * @return void
 """
    cfunc_type = "void"
    name = "diagnostics_nearest"
    params = """commondata_struct *restrict commondata, griddata_struct *restrict griddata,
                const REAL *restrict gridfuncs_diags[MAXNUMGRIDS]"""

    body = r"""
  // --- USER-EDIT: Select diagnostic gridfunctions to sample (applies to all grids) ---

  // 0D diagnostics: nearest point to the grid center.
  const int which_gfs_0d[] = {
      DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF, DIAG_ADDXXGF, DIAG_ADDXYGF, DIAG_ADDXZGF, DIAG_ADDYYGF, DIAG_ADDYZGF,
      DIAG_ADDZZGF, DIAG_TRKGF, DIAG_HDDXXGF, DIAG_HDDXYGF, DIAG_HDDXZGF, DIAG_HDDYYGF, DIAG_HDDYZGF, DIAG_HDDZZGF,
      DIAG_EXACT_ADDXXGF, DIAG_EXACT_ADDXYGF, DIAG_EXACT_ADDXZGF, DIAG_EXACT_ADDYYGF, DIAG_EXACT_ADDYZGF,
      DIAG_EXACT_ADDZZGF, DIAG_EXACT_WGF, DIAG_EXACT_TRKGF, DIAG_EXACT_HDDXXGF, DIAG_EXACT_HDDXYGF, DIAG_EXACT_HDDXZGF,
      DIAG_EXACT_HDDYYGF, DIAG_EXACT_HDDYZGF, DIAG_EXACT_HDDZZGF, DIAG_PX_RECOVEREDGF, DIAG_PY_RECOVEREDGF,
      DIAG_PZ_RECOVEREDGF, DIAG_EXACT_PXGF, DIAG_EXACT_PYGF, DIAG_EXACT_PZGF };

  // 1D diagnostics: nearest lines to the y and z axes.
  const int which_gfs_1d[] = {
      DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF, DIAG_ADDXXGF, DIAG_ADDXYGF, DIAG_ADDXZGF, DIAG_ADDYYGF, DIAG_ADDYZGF,
      DIAG_ADDZZGF, DIAG_TRKGF, DIAG_HDDXXGF, DIAG_HDDXYGF, DIAG_HDDXZGF, DIAG_HDDYYGF, DIAG_HDDYZGF, DIAG_HDDZZGF,
      DIAG_EXACT_ADDXXGF, DIAG_EXACT_ADDXYGF, DIAG_EXACT_ADDXZGF, DIAG_EXACT_ADDYYGF, DIAG_EXACT_ADDYZGF,
      DIAG_EXACT_ADDZZGF, DIAG_EXACT_WGF, DIAG_EXACT_TRKGF, DIAG_EXACT_HDDXXGF, DIAG_EXACT_HDDXYGF, DIAG_EXACT_HDDXZGF,
      DIAG_EXACT_HDDYYGF, DIAG_EXACT_HDDYZGF, DIAG_EXACT_HDDZZGF, DIAG_PX_RECOVEREDGF, DIAG_PY_RECOVEREDGF,
      DIAG_PZ_RECOVEREDGF, DIAG_EXACT_PXGF, DIAG_EXACT_PYGF, DIAG_EXACT_PZGF };

  // 2D diagnostics: nearest planes to the xy and yz coordinate planes.
  const int which_gfs_2d[] = {
      DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF, DIAG_ADDXXGF, DIAG_ADDXYGF, DIAG_ADDXZGF, DIAG_ADDYYGF, DIAG_ADDYZGF,
      DIAG_ADDZZGF, DIAG_TRKGF, DIAG_HDDXXGF, DIAG_HDDXYGF, DIAG_HDDXZGF, DIAG_HDDYYGF, DIAG_HDDYZGF, DIAG_HDDZZGF,
      DIAG_EXACT_ADDXXGF, DIAG_EXACT_ADDXYGF, DIAG_EXACT_ADDXZGF, DIAG_EXACT_ADDYYGF, DIAG_EXACT_ADDYZGF,
      DIAG_EXACT_ADDZZGF, DIAG_EXACT_WGF, DIAG_EXACT_TRKGF, DIAG_EXACT_HDDXXGF, DIAG_EXACT_HDDXYGF, DIAG_EXACT_HDDXZGF,
      DIAG_EXACT_HDDYYGF, DIAG_EXACT_HDDYZGF, DIAG_EXACT_HDDZZGF, DIAG_PX_RECOVEREDGF, DIAG_PY_RECOVEREDGF,
      DIAG_PZ_RECOVEREDGF, DIAG_EXACT_PXGF, DIAG_EXACT_PYGF, DIAG_EXACT_PZGF };

  // --- END USER-EDIT ---

  for (int grid = 0; grid < commondata->NUMGRIDS; grid++) {
    const params_struct *restrict params = &griddata[grid].params;
    const REAL *restrict xx[3] = { griddata[grid].xx[0], griddata[grid].xx[1], griddata[grid].xx[2] };

    const int NUM_nearest_GFS_0d = (int)(sizeof which_gfs_0d / sizeof which_gfs_0d[0]);
    diagnostics_nearest_grid_center(commondata, grid, params, xx, NUM_nearest_GFS_0d, which_gfs_0d,
                                    diagnostic_gf_names, gridfuncs_diags);

    const int NUM_nearest_GFS_1d = (int)(sizeof which_gfs_1d / sizeof which_gfs_1d[0]);
    diagnostics_nearest_1d_y_and_z_axes(commondata, grid, params, xx, NUM_nearest_GFS_1d, which_gfs_1d,
                                        diagnostic_gf_names, gridfuncs_diags);

    const int NUM_nearest_GFS_2d = (int)(sizeof which_gfs_2d / sizeof which_gfs_2d[0]);
    diagnostics_nearest_2d_xy_and_yz_planes(commondata, grid, params, xx, NUM_nearest_GFS_2d, which_gfs_2d,
                                            diagnostic_gf_names, gridfuncs_diags);
  } // END loop over grids
"""
    cfc.register_CFunction(
        subdirectory="diagnostics",
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        name=name,
        params=params,
        include_CodeParameters_h=False,
        body=body,
    )
    return pcg.NRPyEnv()

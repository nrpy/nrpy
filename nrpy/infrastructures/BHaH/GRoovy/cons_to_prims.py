"""
C function registration for GRHayLHD conservative-to-primitive recovery in GRoovy.

Author: Terrence Pierre Jacques
        terrencepierrej **at** gmail **dot** com
"""

import textwrap
from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Union, cast

import sympy as sp
from typing_extensions import Literal

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.equations.basis_transforms.jacobians as bt
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric

# How often to output C2P diagnostics
_ = par.register_CodeParameter(
    "int",
    __name__,
    "C2P_diagnostics_every",
    2,
    commondata=True,
    add_to_parfile=True,
)


def _indent_block(block: str, spaces: int) -> str:
    """
    Indent a multi-line code block by a fixed number of spaces.

    :param block: The code block to indent.
    :param spaces: Number of spaces to prepend to each line.
    :return: The indented block, including a trailing newline when non-empty.
    """
    stripped = textwrap.dedent(block).strip("\n")
    if not stripped:
        return ""
    return textwrap.indent(stripped + "\n", " " * spaces)


def _select_recovery_mode(
    evolving_temperature: bool,
    evolving_entropy: bool,
) -> Literal["Hybrid", "HybridEntropy", "Tabulated", "TabulatedEntropy"]:
    """
    Determine the conservative-to-primitive recovery mode from the evolved variables.

    :param evolving_temperature: Whether temperature and electron fraction are evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :return: The recovery-mode label.
    """
    if evolving_temperature:
        if evolving_entropy:
            return "TabulatedEntropy"
        return "Tabulated"
    if evolving_entropy:
        return "HybridEntropy"
    return "Hybrid"


def _build_store_recovered_primitives_code(
    CoordSystem: str,
    evolving_temperature: bool,
    evolving_entropy: bool,
) -> str:
    """
    Build C code that stores recovered primitives back to the rfm-basis gridfunctions.

    Only primitive data are written here so the conservative averaging logic remains
    deterministic across threads. We intentionally do not call
    `basis_transform_Cartesian_to_rfm_basis()` in this first pass, since that
    helper also updates `evol_gfs`, and changing conservatives in-place would
    contaminate neighboring-point averaging fallbacks.

    :param CoordSystem: The coordinate system.
    :param evolving_temperature: Whether temperature and electron fraction are evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :return: The generated C code snippet.
    """
    # Step 1: Build the symbolic Cartesian-to-rfm velocity transform.
    rfm = refmetric.reference_metric[CoordSystem]
    basis_transforms = bt.basis_transforms[CoordSystem]

    vx, vy, vz = sp.symbols("vx vy vz", real=True)
    vU_Cart = [vx, vy, vz]
    vU = basis_transforms.basis_transform_vectorU_from_Cartesian_to_rfmbasis(vU_Cart)

    rescaledvU = ixp.zerorank1()
    for i in range(3):
        rescaledvU[i] = sp.together(vU[i] / rfm.ReU[i])

    # Step 2: Build the scalar stores and generated velocity assignments.
    pre_body = r"""
// Store the recovered primitive state without touching evol_gfs yet.
// Conserved variables are rewritten in a later pass once all points have
// finished their primitive recovery.
const REAL xx0 = xx[0][i];
const REAL xx1 = xx[1][j];
const REAL xx2 = xx[2][k];

auxevol_gfs[IDX4pt(RHOBGF, index)] = prims.rho;
auxevol_gfs[IDX4pt(PGF, index)] = prims.press;
"""
    if evolving_temperature:
        pre_body += r"""
auxevol_gfs[IDX4pt(YEGF, index)] = prims.Y_e;
auxevol_gfs[IDX4pt(TEMPERATUREGF, index)] = prims.temperature;
"""
    if evolving_entropy:
        pre_body += r"""
auxevol_gfs[IDX4pt(SGF, index)] = prims.entropy;
"""
    pre_body += r"""
const REAL vx = prims.vU[0];
const REAL vy = prims.vU[1];
const REAL vz = prims.vU[2];
"""

    lhs = [
        "auxevol_gfs[IDX4pt(RESCALEDVU0GF, index)]",
        "auxevol_gfs[IDX4pt(RESCALEDVU1GF, index)]",
        "auxevol_gfs[IDX4pt(RESCALEDVU2GF, index)]",
    ]
    rhs = [rescaledvU[0], rescaledvU[1], rescaledvU[2]]

    return pre_body + ccg.c_codegen(
        rhs,
        lhs,
        enable_simd=False,
        enable_fd_codegen=False,
        include_braces=True,
    )


def _build_standard_recovery_block(
    mode: Literal["Hybrid", "HybridEntropy", "Tabulated", "TabulatedEntropy"],
    nan_product: str,
    evolving_temperature: bool,
    evolving_entropy: bool,
) -> str:
    """
    Build the standard recovery block for the non-robust modes.

    :param mode: Recovery-mode label.
    :param nan_product: Product used for NaN detection after a recovery attempt.
    :param evolving_temperature: Whether temperature and electron fraction are evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :return: The generated C block.
    """
    # Step 1: Build shared averaging fragments for the optional channels.
    neighbor_init = r"""
ghl_conservative_quantities cons_neigh_avg, cons_avg;
cons_neigh_avg.rho = 0.0;
cons_neigh_avg.tau = 0.0;
cons_neigh_avg.SD[0] = 0.0;
cons_neigh_avg.SD[1] = 0.0;
cons_neigh_avg.SD[2] = 0.0;
"""
    neighbor_add = r"""
cons_neigh_avg.rho += cons_avg_loop.rho;
cons_neigh_avg.tau += cons_avg_loop.tau;
cons_neigh_avg.SD[0] += cons_avg_loop.SD[0];
cons_neigh_avg.SD[1] += cons_avg_loop.SD[1];
cons_neigh_avg.SD[2] += cons_avg_loop.SD[2];
"""
    avg_mix = r"""
cons_avg.rho = wfac * cons_neigh_avg.rho * inv_n_avg + cfac * cons_orig.rho;
cons_avg.tau = wfac * cons_neigh_avg.tau * inv_n_avg + cfac * cons_orig.tau;
cons_avg.SD[0] = wfac * cons_neigh_avg.SD[0] * inv_n_avg + cfac * cons_orig.SD[0];
cons_avg.SD[1] = wfac * cons_neigh_avg.SD[1] * inv_n_avg + cfac * cons_orig.SD[1];
cons_avg.SD[2] = wfac * cons_neigh_avg.SD[2] * inv_n_avg + cfac * cons_orig.SD[2];
"""

    if evolving_temperature:
        neighbor_init += "cons_neigh_avg.Y_e = 0.0;\n"
        neighbor_add += "cons_neigh_avg.Y_e += cons_avg_loop.Y_e;\n"
        avg_mix += "cons_avg.Y_e = wfac * cons_neigh_avg.Y_e * inv_n_avg + cfac * cons_orig.Y_e;\n"
    if evolving_entropy:
        neighbor_init += "cons_neigh_avg.entropy = 0.0;\n"
        neighbor_add += "cons_neigh_avg.entropy += cons_avg_loop.entropy;\n"
        avg_mix += (
            "cons_avg.entropy = wfac * cons_neigh_avg.entropy * inv_n_avg + "
            "cfac * cons_orig.entropy;\n"
        )

    # Step 2: Specialize the recovery call and hybrid-only fallbacks.
    if mode in ("Hybrid", "HybridEntropy"):
        recovery_call = "ghl_con2prim_hybrid_multi_method"
        pre_recovery = r"""
ghl_apply_conservative_limits(
    ghl_params, eos, &ADM_metric, &prims, &cons, &diagnostics);
"""
        hybrid_backup = r"""
    if (error != ghl_success) {
      pointcount_font++;
      cons = cons_orig;
      ghl_apply_conservative_limits(
          ghl_params, eos, &ADM_metric, &prims, &cons, &diagnostics);
      ghl_undensitize_conservatives(
          ADM_metric.sqrt_detgamma, &cons, &cons_undens);
      error = ghl_hybrid_Font1D(
          ghl_params, eos, &ADM_metric, &metric_aux,
          &cons_undens, &prims, &diagnostics);
      if (isnan("""
        hybrid_backup += nan_product
        hybrid_backup += r"""))
        error = ghl_error_c2p_singular;
    } // END IF: weighted-average recovery exhausted the hybrid multi-methods
"""
    else:
        recovery_call = "ghl_con2prim_tabulated_multi_method"
        pre_recovery = ""
        hybrid_backup = ""

    # Step 3: Build the full non-robust recovery logic.
    block = r"""
if (cons.rho > 0.0) {
"""
    block += _indent_block(pre_recovery, 2)
    block += r"""
  ghl_conservative_quantities cons_undens;
  ghl_undensitize_conservatives(
      ADM_metric.sqrt_detgamma, &cons, &cons_undens);
  error = """
    block += recovery_call
    block += r"""(
      ghl_params, eos, &ADM_metric, &metric_aux,
      &cons_undens, &prims, &diagnostics);
  if (isnan("""
    block += nan_product
    block += r"""))
    error = ghl_error_c2p_singular;
} else {
  ghl_set_prims_to_constant_atm(eos, &prims);
  rho_star_fix_applied++;
} // END IF: conservative density is positive

if (error != ghl_success) {
  pointcount_avg++;
  cons = cons_orig;
"""
    block += _indent_block(neighbor_init, 2)
    block += r"""
  const int iavg_min = MAX(imin, i - 1);
  const int javg_min = MAX(jmin, j - 1);
  const int kavg_min = MAX(kmin, k - 1);
  const int iavg_max = MIN(imax, i + 2);
  const int javg_max = MIN(jmax, j + 2);
  const int kavg_max = MIN(kmax, k + 2);

  int n_avg = 0;
  for (int kavg = kavg_min; kavg < kavg_max; kavg++) {
    for (int javg = javg_min; javg < javg_max; javg++) {
      for (int iavg = iavg_min; iavg < iavg_max; iavg++) {
        const int index_avg = IDX3(iavg, javg, kavg);
        if (index_avg == index)
          continue;

        ghl_conservative_quantities cons_avg_loop;
        ghl_metric_quantities ADM_metric_avg;
        basis_transform_rfm_basis_to_Cartesian__read_cons_only(
            commondata, params, &cons_avg_loop, &ADM_metric_avg,
            iavg, javg, kavg, xx, auxevol_gfs, evol_gfs);
"""
    block += _indent_block(neighbor_add, 8)
    block += r"""
        n_avg++;
      } // END LOOP: for iavg over neighboring x indices
    } // END LOOP: for javg over neighboring y indices
  } // END LOOP: for kavg over neighboring z indices

  int avg_weight = 1;
  while (error != ghl_success && avg_weight < 5) {
    const REAL wfac = ((REAL)avg_weight) / 4.0;
    const REAL cfac = 1.0 - wfac;
    const REAL inv_n_avg = 1.0 / (REAL)n_avg;
"""
    block += _indent_block(avg_mix, 4)
    block += r"""
    ghl_undensitize_conservatives(
        ADM_metric.sqrt_detgamma, &cons_avg, &cons_undens);
    error = """
    block += recovery_call
    block += r"""(
        ghl_params, eos, &ADM_metric, &metric_aux,
        &cons_undens, &prims, &diagnostics);
    avg_weight++;
    if (isnan("""
    block += nan_product
    block += r"""))
      error = ghl_error_c2p_singular;
  } // END WHILE: weighted averages remain available
"""
    block += _indent_block(hybrid_backup, 2)
    block += r"""
  if (error != ghl_success) {
    ghl_set_prims_to_constant_atm(eos, &prims);
    failures++;
    failures_inhoriz += in_horizon;
  } // END IF: all conservative-to-primitive fallbacks failed
} // END IF: initial conservative-to-primitive recovery failed
"""
    return block


def _build_tabulated_entropy_robust_block(nan_product: str) -> str:
    """
    Build the robust tabulated-entropy recovery block.

    :param nan_product: Product used for NaN detection after a recovery attempt.
    :return: The generated C block.
    """
    nan_product_1 = nan_product.replace("prims.", "prims1.")
    nan_product_2 = nan_product.replace("prims.", "prims2.")
    nan_product_3 = nan_product.replace("prims.", "prims3.")
    nan_product_4 = nan_product.replace("prims.", "prims4.")

    return (
        r"""
if (cons.rho > 0.0) {
  ghl_primitive_quantities prims1, prims2, prims3, prims4;
  ghl_conservative_quantities cons_undens1, cons_undens2, cons_undens3, cons_undens4;
  prims1 = prims;
  prims2 = prims;
  prims3 = prims;
  prims4 = prims;

  ghl_undensitize_conservatives(
      ADM_metric.sqrt_detgamma, &cons, &cons_undens1);
  cons_undens2 = cons_undens1;
  cons_undens3 = cons_undens1;
  cons_undens4 = cons_undens1;

  ghl_error_codes_t error1 = ghl_tabulated_Palenzuela1D_energy(
      ghl_params, eos, &ADM_metric, &metric_aux,
      &cons_undens1, &prims1, &diagnostics);
  ghl_error_codes_t error2 = ghl_tabulated_Newman1D_energy(
      ghl_params, eos, &ADM_metric, &metric_aux,
      &cons_undens2, &prims2, &diagnostics);
  ghl_error_codes_t error3 = ghl_tabulated_Palenzuela1D_entropy(
      ghl_params, eos, &ADM_metric, &metric_aux,
      &cons_undens3, &prims3, &diagnostics);
  ghl_error_codes_t error4 = ghl_tabulated_Newman1D_entropy(
      ghl_params, eos, &ADM_metric, &metric_aux,
      &cons_undens4, &prims4, &diagnostics);
  if (error1 == ghl_success && isnan("""
        + nan_product_1
        + r"""))
    error1 = ghl_error_c2p_singular;
  if (error2 == ghl_success && isnan("""
        + nan_product_2
        + r"""))
    error2 = ghl_error_c2p_singular;
  if (error3 == ghl_success && isnan("""
        + nan_product_3
        + r"""))
    error3 = ghl_error_c2p_singular;
  if (error4 == ghl_success && isnan("""
        + nan_product_4
        + r"""))
    error4 = ghl_error_c2p_singular;

  REAL err1 = 1e300;
  REAL err2 = 1e300;
  REAL err3 = 1e300;
  REAL err4 = 1e300;
  bool speed_limited_dummy = false;

  // Note: in the following lines we compare the conserved energy and electron
  // fraction between routines. We do not include entropy, since it is not
  // conserved at shocks. We ignore the diagnostic from
  // ghl_enforce_primitive_limits_and_compute_u0 here; a candidate that fails
  // primitive limiting is unlikely to be selected.

  if (error1 == ghl_success) {
    ghl_conservative_quantities cons_temp;
    ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, eos, &ADM_metric, &prims1, &speed_limited_dummy);
    ghl_compute_conservs(&ADM_metric, &metric_aux, &prims1, &cons_temp);
    err1 = fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100))
         + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
  } // END IF: Palenzuela energy recovery succeeded

  if (error2 == ghl_success) {
    ghl_conservative_quantities cons_temp;
    ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, eos, &ADM_metric, &prims2, &speed_limited_dummy);
    ghl_compute_conservs(&ADM_metric, &metric_aux, &prims2, &cons_temp);
    err2 = fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100))
         + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
  } // END IF: Newman energy recovery succeeded

  if (error3 == ghl_success) {
    ghl_conservative_quantities cons_temp;
    ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, eos, &ADM_metric, &prims3, &speed_limited_dummy);
    ghl_compute_conservs(&ADM_metric, &metric_aux, &prims3, &cons_temp);
    err3 = fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100))
         + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
  } // END IF: Palenzuela entropy recovery succeeded

  if (error4 == ghl_success) {
    ghl_conservative_quantities cons_temp;
    ghl_enforce_primitive_limits_and_compute_u0(
        ghl_params, eos, &ADM_metric, &prims4, &speed_limited_dummy);
    ghl_compute_conservs(&ADM_metric, &metric_aux, &prims4, &cons_temp);
    err4 = fabs((cons_temp.Y_e - cons_orig.Y_e) / (cons_orig.Y_e + 1e-100))
         + fabs((cons_temp.tau - cons_orig.tau) / (cons_orig.tau + 1e-100));
  } // END IF: Newman entropy recovery succeeded

  REAL min_err = err1;
  int best_method_index = (error1 == ghl_success) ? 1 : 0;
  if (err2 < min_err) {
    min_err = err2;
    best_method_index = 2;
  } // END IF: Newman energy gives a smaller mismatch
  if (err3 < min_err) {
    min_err = err3;
    best_method_index = 3;
  } // END IF: Palenzuela entropy gives a smaller mismatch
  if (err4 < min_err) {
    min_err = err4;
    best_method_index = 4;
  } // END IF: Newman entropy gives a smaller mismatch

  if (best_method_index == 0) {
    pointcount_avg++;
    ghl_conservative_quantities cons_neigh_avg, cons_avg;
    cons_neigh_avg.rho = 0.0;
    cons_neigh_avg.tau = 0.0;
    cons_neigh_avg.SD[0] = 0.0;
    cons_neigh_avg.SD[1] = 0.0;
    cons_neigh_avg.SD[2] = 0.0;
    cons_neigh_avg.entropy = 0.0;
    cons_neigh_avg.Y_e = 0.0;

    const int iavg_min = MAX(imin, i - 1);
    const int javg_min = MAX(jmin, j - 1);
    const int kavg_min = MAX(kmin, k - 1);
    const int iavg_max = MIN(imax, i + 2);
    const int javg_max = MIN(jmax, j + 2);
    const int kavg_max = MIN(kmax, k + 2);

    int n_avg = 0;
    for (int kavg = kavg_min; kavg < kavg_max; kavg++) {
      for (int javg = javg_min; javg < javg_max; javg++) {
        for (int iavg = iavg_min; iavg < iavg_max; iavg++) {
          const int index_avg = IDX3(iavg, javg, kavg);
          if (index_avg == index)
            continue;

          ghl_conservative_quantities cons_avg_loop;
          ghl_metric_quantities ADM_metric_avg;
          basis_transform_rfm_basis_to_Cartesian__read_cons_only(
              commondata, params, &cons_avg_loop, &ADM_metric_avg,
              iavg, javg, kavg, xx, auxevol_gfs, evol_gfs);
          cons_neigh_avg.rho += cons_avg_loop.rho;
          cons_neigh_avg.tau += cons_avg_loop.tau;
          cons_neigh_avg.SD[0] += cons_avg_loop.SD[0];
          cons_neigh_avg.SD[1] += cons_avg_loop.SD[1];
          cons_neigh_avg.SD[2] += cons_avg_loop.SD[2];
          cons_neigh_avg.entropy += cons_avg_loop.entropy;
          cons_neigh_avg.Y_e += cons_avg_loop.Y_e;
          n_avg++;
        } // END LOOP: for iavg over neighboring x indices
      } // END LOOP: for javg over neighboring y indices
    } // END LOOP: for kavg over neighboring z indices

    for (int avg_weight = 1; avg_weight <= 4 && best_method_index == 0; avg_weight++) {
      const REAL w_neigh = ((REAL)avg_weight) / 4.0;
      const REAL w_self = 1.0 - w_neigh;
      const REAL inv_n_avg = 1.0 / (REAL)n_avg;

      cons_avg.rho = (w_neigh * cons_neigh_avg.rho * inv_n_avg + w_self * cons_orig.rho);
      cons_avg.tau = (w_neigh * cons_neigh_avg.tau * inv_n_avg + w_self * cons_orig.tau);
      cons_avg.SD[0] = (w_neigh * cons_neigh_avg.SD[0] * inv_n_avg + w_self * cons_orig.SD[0]);
      cons_avg.SD[1] = (w_neigh * cons_neigh_avg.SD[1] * inv_n_avg + w_self * cons_orig.SD[1]);
      cons_avg.SD[2] = (w_neigh * cons_neigh_avg.SD[2] * inv_n_avg + w_self * cons_orig.SD[2]);
      cons_avg.entropy = (w_neigh * cons_neigh_avg.entropy * inv_n_avg + w_self * cons_orig.entropy);
      cons_avg.Y_e = (w_neigh * cons_neigh_avg.Y_e * inv_n_avg + w_self * cons_orig.Y_e);

      prims1 = prims;
      prims2 = prims;
      prims3 = prims;
      prims4 = prims;
      ghl_undensitize_conservatives(
          ADM_metric.sqrt_detgamma, &cons_avg, &cons_undens1);
      cons_undens2 = cons_undens1;
      cons_undens3 = cons_undens1;
      cons_undens4 = cons_undens1;

      error1 = ghl_tabulated_Palenzuela1D_energy(
          ghl_params, eos, &ADM_metric, &metric_aux,
          &cons_undens1, &prims1, &diagnostics);
      error2 = ghl_tabulated_Newman1D_energy(
          ghl_params, eos, &ADM_metric, &metric_aux,
          &cons_undens2, &prims2, &diagnostics);
      error3 = ghl_tabulated_Palenzuela1D_entropy(
          ghl_params, eos, &ADM_metric, &metric_aux,
          &cons_undens3, &prims3, &diagnostics);
      error4 = ghl_tabulated_Newman1D_entropy(
          ghl_params, eos, &ADM_metric, &metric_aux,
          &cons_undens4, &prims4, &diagnostics);
      if (error1 == ghl_success && isnan("""
        + nan_product_1
        + r"""))
        error1 = ghl_error_c2p_singular;
      if (error2 == ghl_success && isnan("""
        + nan_product_2
        + r"""))
        error2 = ghl_error_c2p_singular;
      if (error3 == ghl_success && isnan("""
        + nan_product_3
        + r"""))
        error3 = ghl_error_c2p_singular;
      if (error4 == ghl_success && isnan("""
        + nan_product_4
        + r"""))
        error4 = ghl_error_c2p_singular;

      err1 = 1e300;
      err2 = 1e300;
      err3 = 1e300;
      err4 = 1e300;

      if (error1 == ghl_success) {
        ghl_conservative_quantities cons_temp;
        ghl_enforce_primitive_limits_and_compute_u0(
            ghl_params, eos, &ADM_metric, &prims1, &speed_limited_dummy);
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims1, &cons_temp);
        err1 = fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100))
             + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
      } // END IF: averaged Palenzuela energy recovery succeeded

      if (error2 == ghl_success) {
        ghl_conservative_quantities cons_temp;
        ghl_enforce_primitive_limits_and_compute_u0(
            ghl_params, eos, &ADM_metric, &prims2, &speed_limited_dummy);
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims2, &cons_temp);
        err2 = fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100))
             + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
      } // END IF: averaged Newman energy recovery succeeded

      if (error3 == ghl_success) {
        ghl_conservative_quantities cons_temp;
        ghl_enforce_primitive_limits_and_compute_u0(
            ghl_params, eos, &ADM_metric, &prims3, &speed_limited_dummy);
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims3, &cons_temp);
        err3 = fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100))
             + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
      } // END IF: averaged Palenzuela entropy recovery succeeded

      if (error4 == ghl_success) {
        ghl_conservative_quantities cons_temp;
        ghl_enforce_primitive_limits_and_compute_u0(
            ghl_params, eos, &ADM_metric, &prims4, &speed_limited_dummy);
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims4, &cons_temp);
        err4 = fabs((cons_temp.Y_e - cons_avg.Y_e) / (cons_avg.Y_e + 1e-100))
             + fabs((cons_temp.tau - cons_avg.tau) / (cons_avg.tau + 1e-100));
      } // END IF: averaged Newman entropy recovery succeeded

      min_err = err1;
      best_method_index = (error1 == ghl_success) ? 1 : 0;
      if (err2 < min_err) {
        min_err = err2;
        best_method_index = 2;
      } // END IF: averaged Newman energy gives a smaller mismatch
      if (err3 < min_err) {
        min_err = err3;
        best_method_index = 3;
      } // END IF: averaged Palenzuela entropy gives a smaller mismatch
      if (err4 < min_err) {
        min_err = err4;
        best_method_index = 4;
      } // END IF: averaged Newman entropy gives a smaller mismatch
    } // END LOOP: for avg_weight over weighted averaging attempts
  } // END IF: all direct tabulated-entropy methods failed

  if (best_method_index == 1) {
    prims = prims1;
  } else if (best_method_index == 2) {
    prims = prims2;
  } else if (best_method_index == 3) {
    prims = prims3;
  } else if (best_method_index == 4) {
    prims = prims4;
  } else {
    ghl_set_prims_to_constant_atm(eos, &prims);
    failures++;
    failures_inhoriz += in_horizon;
  } // END IF: robust tabulated-entropy method selection completed
} else {
  ghl_set_prims_to_constant_atm(eos, &prims);
  rho_star_fix_applied++;
} // END IF: conservative density is positive
"""
    )


def _build_cons_to_prims_body(
    CoordSystem: str,
    mode: Literal["Hybrid", "HybridEntropy", "Tabulated", "TabulatedEntropy"],
    tabulated_entropy_robust: bool,
) -> str:
    """
    Build the C body for a conservative-to-primitive recovery mode.

    :param CoordSystem: The coordinate system.
    :param mode: Recovery-mode label.
    :param tabulated_entropy_robust: Whether to use the robust tabulated-entropy recovery.
    :return: The generated C function body.
    """
    # Step 1: Determine the optional channels and mode-specific symbolic snippets.
    evolving_temperature = mode in ("Tabulated", "TabulatedEntropy")
    evolving_entropy = mode in ("HybridEntropy", "TabulatedEntropy")
    store_primitives_code = _build_store_recovered_primitives_code(
        CoordSystem,
        evolving_temperature,
        evolving_entropy,
    )

    if mode == "Hybrid":
        guess_setup = r"""
if (prims.rho > 0.0) {
  REAL P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(
      eos, prims.rho, &P_cold, &eps_cold);
  const REAL eps_th = (prims.press - P_cold) / ((eos->Gamma_th - 1.0) * prims.rho);
  prims.eps = eps_cold + eps_th;
} // END IF: the primitive density is positive for the hybrid initial guess
"""
        nan_product = "prims.rho * prims.press * prims.eps * prims.vU[0] * prims.vU[1] * prims.vU[2]"
    elif mode == "HybridEntropy":
        guess_setup = r"""
if (prims.rho > 0.0) {
  REAL P_cold, eps_cold;
  ghl_hybrid_compute_P_cold_and_eps_cold(
      eos, prims.rho, &P_cold, &eps_cold);
  const REAL eps_th = (prims.press - P_cold) / ((eos->Gamma_th - 1.0) * prims.rho);
  prims.eps = eps_cold + eps_th;
} // END IF: the primitive density is positive for the hybrid-entropy initial guess
"""
        nan_product = (
            "prims.rho * prims.press * prims.eps * prims.vU[0] * prims.vU[1] * "
            "prims.vU[2] * prims.entropy"
        )
    elif mode == "Tabulated":
        guess_setup = r"""
ghl_tabulated_enforce_bounds_rho_Ye_P(
    eos, &prims.rho, &prims.Y_e, &prims.press);
ghl_tabulated_compute_eps_T_from_P(
    eos, prims.rho, prims.Y_e, prims.press, &prims.eps, &prims.temperature);
"""
        nan_product = (
            "prims.rho * prims.press * prims.eps * prims.vU[0] * prims.vU[1] * "
            "prims.vU[2] * prims.Y_e * prims.temperature"
        )
    else:
        guess_setup = r"""
ghl_tabulated_enforce_bounds_rho_Ye_P(
    eos, &prims.rho, &prims.Y_e, &prims.press);
ghl_tabulated_compute_eps_T_from_P(
    eos, prims.rho, prims.Y_e, prims.press, &prims.eps, &prims.temperature);
"""
        nan_product = (
            "prims.rho * prims.press * prims.eps * prims.vU[0] * prims.vU[1] * "
            "prims.vU[2] * prims.entropy * prims.Y_e * prims.temperature"
        )

    if tabulated_entropy_robust:
        recovery_block = _build_tabulated_entropy_robust_block(nan_product)
    else:
        recovery_block = _build_standard_recovery_block(
            mode,
            nan_product,
            evolving_temperature,
            evolving_entropy,
        )

    # Step 2: Build the optional conservative-error tracking fragments.
    error_decls = r"""
  REAL error_rho_numer = 0.0;
  REAL error_tau_numer = 0.0;
  REAL error_Sx_numer = 0.0;
  REAL error_Sy_numer = 0.0;
  REAL error_Sz_numer = 0.0;

  REAL error_rho_denom = 0.0;
  REAL error_tau_denom = 0.0;
  REAL error_Sx_denom = 0.0;
  REAL error_Sy_denom = 0.0;
  REAL error_Sz_denom = 0.0;
"""
    error_reduction = (
        "error_rho_numer, error_tau_numer, error_Sx_numer, error_Sy_numer, "
        "error_Sz_numer, error_rho_denom, error_tau_denom, error_Sx_denom, "
        "error_Sy_denom, error_Sz_denom"
    )
    error_update = r"""
      error_rho_numer += fabs(cons.rho - cons_orig.rho);
      error_tau_numer += fabs(cons.tau - cons_orig.tau);
      error_Sx_numer += fabs(cons.SD[0] - cons_orig.SD[0]);
      error_Sy_numer += fabs(cons.SD[1] - cons_orig.SD[1]);
      error_Sz_numer += fabs(cons.SD[2] - cons_orig.SD[2]);

      error_rho_denom += cons_orig.rho;
      error_tau_denom += cons_orig.tau;
      error_Sx_denom += fabs(cons_orig.SD[0]);
      error_Sy_denom += fabs(cons_orig.SD[1]);
      error_Sz_denom += fabs(cons_orig.SD[2]);
"""
    optional_error_calcs = ""
    optional_error_fmt = ""
    optional_error_args = ""

    if evolving_entropy:
        error_decls += r"""
  REAL error_ent_numer = 0.0;
  REAL error_ent_denom = 0.0;
"""
        error_reduction += ", error_ent_numer, error_ent_denom"
        error_update += r"""
      error_ent_numer += fabs(cons.entropy - cons_orig.entropy);
      error_ent_denom += cons_orig.entropy;
"""
        optional_error_calcs += r"""
    const REAL ent_error = (error_ent_denom == 0.0)
                               ? error_ent_numer
                               : error_ent_numer / error_ent_denom;
"""
        optional_error_fmt += " | entropy %.3e, %.3e"
        optional_error_args += r"""
        ent_error,
        error_ent_denom,
"""

    if evolving_temperature:
        error_decls += r"""
  REAL error_Ye_numer = 0.0;
  REAL error_Ye_denom = 0.0;
"""
        error_reduction += ", error_Ye_numer, error_Ye_denom"
        error_update += r"""
      error_Ye_numer += fabs(cons.Y_e - cons_orig.Y_e);
      error_Ye_denom += cons_orig.Y_e;
"""
        optional_error_calcs += r"""
    const REAL Ye_error = (error_Ye_denom == 0.0)
                              ? error_Ye_numer
                              : error_Ye_numer / error_Ye_denom;
"""
        optional_error_fmt += " | Y_e %.3e, %.3e"
        optional_error_args += r"""
        Ye_error,
        error_Ye_denom,
"""

    # Step 3: Specialize the recovery label and hybrid-only diagnostic counts.
    if tabulated_entropy_robust:
        recovery_name = "TabulatedEntropyRobust"
    else:
        recovery_name = mode

    if mode in ("Hybrid", "HybridEntropy") and not tabulated_entropy_robust:
        pointcount_font_decl = "  int pointcount_font = 0;\n"
        first_reduction = (
            "backup0, backup1, backup2, vel_limited_ptcount, "
            "rho_star_fix_applied, failures, failures_inhoriz, pointcount_inhoriz, "
            "n_iter, pointcount_avg, pointcount_font"
        )
        diagnostic_line2 = '        "         Averaged= %d Font1D= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\\n"\n'
        diagnostic_line2_args = r"""
        pointcount_avg,
        pointcount_font,
        failures,
        failures_inhoriz,
        pointcount_inhoriz,
        avg_iters,
"""
    else:
        pointcount_font_decl = ""
        first_reduction = (
            "backup0, backup1, backup2, vel_limited_ptcount, "
            "rho_star_fix_applied, failures, failures_inhoriz, pointcount_inhoriz, "
            "n_iter, pointcount_avg"
        )
        diagnostic_line2 = '        "         Averaged= %d | Failures: %d InHoriz= %d / %d | %.2f iters/gridpt\\n"\n'
        diagnostic_line2_args = r"""
        pointcount_avg,
        failures,
        failures_inhoriz,
        pointcount_inhoriz,
        avg_iters,
"""

    # Step 4: Assemble the final body from the mode-specific fragments.
    body = r"""
  // Step 1: Set up the interior-loop bounds and recovery diagnostics.
  const int imin = NGHOSTS;
  const int imax = Nxx_plus_2NGHOSTS0 - NGHOSTS;
  const int jmin = NGHOSTS;
  const int jmax = Nxx_plus_2NGHOSTS1 - NGHOSTS;
  const int kmin = NGHOSTS;
  const int kmax = Nxx_plus_2NGHOSTS2 - NGHOSTS;
  const int pointcount = (imax - imin) * (jmax - jmin) * (kmax - kmin);
  const char *recovery_name = """
    body += f'"{recovery_name}"'
    body += r""";

  int failures = 0;
  int failures_inhoriz = 0;
  int pointcount_inhoriz = 0;
  int vel_limited_ptcount = 0;
  int rho_star_fix_applied = 0;
  int backup0 = 0;
  int backup1 = 0;
  int backup2 = 0;
  int n_iter = 0;
  int pointcount_avg = 0;
"""
    body += pointcount_font_decl
    body += r"""
#pragma omp parallel for reduction(+: """
    body += first_reduction
    body += r""") schedule(static)
  for (int k = kmin; k < kmax; k++) {
    for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
        const int index = IDX3(i, j, k);

        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_orig;
        ghl_con2prim_diagnostics diagnostics;
        ghl_initialize_diagnostics(&diagnostics);

        ghl_metric_quantities ADM_metric;
        basis_transform_rfm_basis_to_Cartesian(
            commondata, params, &prims, &cons, &ADM_metric,
            i, j, k, xx, auxevol_gfs, evol_gfs);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        const int in_horizon = (ADM_metric.sqrt_detgamma > ghl_params->psi6threshold);
        pointcount_inhoriz += in_horizon;
        cons_orig = cons;
"""
    body += _indent_block(guess_setup, 8)
    body += r"""
        ghl_error_codes_t error = ghl_success;
"""
    body += _indent_block(recovery_block, 8)
    body += r"""
        ghl_enforce_primitive_limits_and_compute_u0(ghl_params, eos, &ADM_metric, &prims, &diagnostics.speed_limited);
        if (error != ghl_success) {
          ghl_set_prims_to_constant_atm(eos, &prims);
          failures++;
          failures_inhoriz += in_horizon;
          ghl_enforce_primitive_limits_and_compute_u0(
              ghl_params, eos, &ADM_metric, &prims, &diagnostics.speed_limited);
        } // END IF: primitive post-processing failed and atmosphere fallback was applied

        if (diagnostics.speed_limited)
          vel_limited_ptcount++;
        backup0 += diagnostics.backup[0];
        backup1 += diagnostics.backup[1];
        backup2 += diagnostics.backup[2];
        n_iter += diagnostics.n_iter;
"""
    body += _indent_block(store_primitives_code, 8)
    body += r"""
      } // END LOOP: for i over interior x indices
    } // END LOOP: for j over interior y indices
  } // END LOOP: for k over interior z indices

  // Step 2: Recompute conservatives from the recovered primitives in a separate pass.
"""
    body += error_decls
    body += r"""
#pragma omp parallel for reduction(+: """
    body += error_reduction
    body += r""") schedule(static)
  for (int k = kmin; k < kmax; k++) {
    for (int j = jmin; j < jmax; j++) {
      for (int i = imin; i < imax; i++) {
        const int index = IDX3(i, j, k);

        ghl_primitive_quantities prims;
        ghl_conservative_quantities cons, cons_orig;
        ghl_metric_quantities ADM_metric;
        basis_transform_rfm_basis_to_Cartesian(
            commondata, params, &prims, &cons, &ADM_metric,
            i, j, k, xx, auxevol_gfs, evol_gfs);

        ghl_ADM_aux_quantities metric_aux;
        ghl_compute_ADM_auxiliaries(&ADM_metric, &metric_aux);

        cons_orig = cons;
        ghl_compute_conservs(&ADM_metric, &metric_aux, &prims, &cons);
        basis_transform_Cartesian_to_rfm_basis(
            commondata, params, &prims, &cons, i, j, k, xx, auxevol_gfs, evol_gfs);
        auxevol_gfs[IDX4pt(U4UTGF, index)] = prims.u0;
"""
    body += error_update
    body += r"""
      } // END LOOP: for i over interior x indices
    } // END LOOP: for j over interior y indices
  } // END LOOP: for k over interior z indices

  // Step 3: Emit periodic conservative-to-primitive diagnostics.
  if (commondata->C2P_diagnostics_every > 0
      && (commondata->nn % commondata->C2P_diagnostics_every == 0)) {
    const REAL rho_error = (error_rho_denom == 0.0)
                               ? error_rho_numer
                               : error_rho_numer / error_rho_denom;
    const REAL tau_error = (error_tau_denom == 0.0)
                               ? error_tau_numer
                               : error_tau_numer / error_tau_denom;
    const REAL Sx_error = (error_Sx_denom == 0.0)
                              ? error_Sx_numer
                              : error_Sx_numer / error_Sx_denom;
    const REAL Sy_error = (error_Sy_denom == 0.0)
                              ? error_Sy_numer
                              : error_Sy_numer / error_Sy_denom;
    const REAL Sz_error = (error_Sz_denom == 0.0)
                              ? error_Sz_numer
                              : error_Sz_numer / error_Sz_denom;
"""
    body += optional_error_calcs
    body += r"""
    const REAL avg_iters = (pointcount == 0)
                               ? 0.0
                               : ((REAL)n_iter) / ((REAL)pointcount);

    printf(
        "C2P[%s]: NumPts= %d | Backups: %d %d %d | Fixes: VL= %d rho*= %d\n"
"""
    body += diagnostic_line2
    body += (
        '        "         Error, Sum: rho %.3e, %.3e | tau %.3e, %.3e'
        + optional_error_fmt
        + '\\n"\n'
    )
    body += r"""
        "                     Sx %.3e, %.3e | Sy %.3e, %.3e | Sz %.3e, %.3e\n",
        recovery_name,
        pointcount,
        backup0,
        backup1,
        backup2,
        vel_limited_ptcount,
        rho_star_fix_applied,
"""
    body += diagnostic_line2_args
    body += r"""
        rho_error,
        error_rho_denom,
        tau_error,
        error_tau_denom,
"""
    body += optional_error_args
    body += r"""
        Sx_error,
        error_Sx_denom,
        Sy_error,
        error_Sy_denom,
        Sz_error,
        error_Sz_denom);
  } // END IF: conservative-to-primitive diagnostics are enabled this timestep
"""
    return body


def register_CFunction_cons_to_prims(
    CoordSystem: str,
    evolving_temperature: bool = False,
    evolving_entropy: bool = False,
    tabulated_entropy_robust: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register the GRHayLHD conservative-to-primitive recovery routine for GRoovy.

    The supported output modes are:
    `Hybrid`, `Hybrid+Entropy`, `Tabulated`, `Tabulated+Entropy`, and
    `Tabulated+Entropy+Robust`.

    This implementation intentionally retains the current GRoovy behavior of
    always seeding the recovery from the existing primitive gridfunctions, then
    updating the surrounding recovery logic to follow the newer GRHayLHD_ET
    control flow as closely as possible in the BHaH setting.

    :param CoordSystem: The coordinate system.
    :param evolving_temperature: Whether temperature and electron fraction are evolved.
    :param evolving_entropy: Whether entropy is evolved.
    :param tabulated_entropy_robust: Whether to use the robust tabulated-entropy
        recovery that tries multiple tabulated methods explicitly.
    :return: None if in registration phase, else the updated NRPy environment.
    :raises ValueError: If `tabulated_entropy_robust` is requested outside the
        tabulated-entropy mode.

    Doctests:
    >>> import contextlib
    >>> import io
    >>> from nrpy.helpers.generic import clang_format, validate_strings
    >>> import nrpy.c_function as cfc
    >>> import nrpy.params as par
    >>> par.set_parval_from_str("parallelization", "openmp")
    >>> cfc.CFunction_dict.clear()
    >>> with contextlib.redirect_stdout(io.StringIO()):
    ...     _ = register_CFunction_cons_to_prims(
    ...         "Spherical",
    ...         evolving_temperature=True,
    ...         evolving_entropy=True,
    ...         tabulated_entropy_robust=True,
    ...     )
    >>> generated_str = clang_format(
    ...     cfc.CFunction_dict["cons_to_prims__rfm__Spherical"].full_function
    ... )
    >>> _ = validate_strings(
    ...     generated_str,
    ...     "cons_to_prims__openmp__Spherical__tabulated_entropy_robust",
    ...     file_ext="c",
    ... )
    """
    # Step 1: Validate the requested recovery mode.
    if tabulated_entropy_robust and not (evolving_temperature and evolving_entropy):
        raise ValueError(
            "tabulated_entropy_robust requires evolving_temperature=True "
            "and evolving_entropy=True."
        )

    # Step 2: Register the function with NRPy's parallel codegen infrastructure.
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    # Step 3: Build the mode-specific C function body.
    mode = _select_recovery_mode(evolving_temperature, evolving_entropy)
    strategy_desc = (
        "strategy with the robust tabulated-entropy fallback"
        if tabulated_entropy_robust
        else "strategy"
    )
    desc = f"""Recover primitive variables from conservative GRHD variables using the {mode} GRHayLHD {strategy_desc}.

@param[in] commondata Common simulation data and diagnostics settings.
@param[in] params Grid-local runtime parameters.
@param[in] ghl_params GRHayL hydrodynamics parameters.
@param[in] eos GRHayL equation-of-state parameters.
@param[in] xx Reference-metric coordinate arrays.
@param[in,out] evol_gfs Conservative gridfunctions.
@param[in,out] auxevol_gfs Primitive and auxiliary gridfunctions."""

    includes = ["BHaH_defines.h", "BHaH_function_prototypes.h"]
    body = _build_cons_to_prims_body(
        CoordSystem,
        mode,
        tabulated_entropy_robust,
    )

    # Step 4: Register the final C function.
    cfunc_type = "void"
    name = "cons_to_prims"
    params = (
        "const commondata_struct *restrict commondata, "
        "const params_struct *restrict params, "
        "const ghl_parameters *restrict ghl_params, "
        "const ghl_eos_parameters *restrict eos, "
        "REAL *restrict xx[3], "
        "REAL *restrict evol_gfs, "
        "REAL *restrict auxevol_gfs"
    )

    cfc.register_CFunction(
        include_CodeParameters_h=True,
        includes=includes,
        desc=desc,
        cfunc_type=cfunc_type,
        CoordSystem_for_wrapper_func=CoordSystem,
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

"""
Kasner-specific diagnostics helpers for BHaH infrastructures.

Author: Nishita Jadoo
        njadoo **at** uidaho **dot* com
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import Dict, List, Sequence, Tuple, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.basis_transforms.jacobians import BasisTransforms
from nrpy.equations.general_relativity.ADM_to_BSSN import ADM_to_BSSN
from nrpy.equations.general_relativity.InitialData_Cartesian import (
    kasner_adm_quantities,
)
from nrpy.infrastructures.BHaH.general_relativity import (
    diagnostic_gfs_set as gr_diagnostic_gfs_set,
)
from nrpy.infrastructures.BHaH.general_relativity import (
    diagnostics_nearest as gr_diagnostics_nearest,
)

# (gridfunction enum name, human-readable description)
KASNER_DIAG_GRIDFUNCTIONS: Tuple[Tuple[str, str], ...] = (
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

KASNER_NEAREST_DIAG_GFS_BHAH: Tuple[str, ...] = tuple(
    f"{name}GF" for name, _ in KASNER_DIAG_GRIDFUNCTIONS
)
KASNER_NEAREST_DIAG_GFS_SUPERB: Tuple[str, ...] = (
    "DIAG_ADDXXGF",
    "DIAG_ADDXYGF",
    "DIAG_ADDXZGF",
    "DIAG_ADDYYGF",
    "DIAG_ADDYZGF",
    "DIAG_ADDZZGF",
    "DIAG_TRKGF",
    "DIAG_HDDXXGF",
    "DIAG_HDDXYGF",
    "DIAG_HDDXZGF",
    "DIAG_HDDYYGF",
    "DIAG_HDDYZGF",
    "DIAG_HDDZZGF",
)


def register_kasner_diag_gridfunctions() -> None:
    """Register all Kasner-specific diagnostic gridfunctions in the ``DIAG`` group."""
    for name, desc in KASNER_DIAG_GRIDFUNCTIONS:
        gri.register_gridfunctions(names=[name], desc_list=[desc], group="DIAG")


def kasner_nearest_diag_gf_names_bhah() -> List[str]:
    """
    Return Kasner diagnostic GF enum names used by non-superB nearest diagnostics.

    :return: List of Kasner diagnostic gridfunction enum names.
    """
    return [
        "DIAG_HAMILTONIANGF",
        "DIAG_MSQUAREDGF",
        *list(KASNER_NEAREST_DIAG_GFS_BHAH),
    ]


def kasner_nearest_diag_gf_names_superb() -> List[str]:
    """
    Return Kasner diagnostic GF enum names used by superB nearest diagnostics.

    :return: List of Kasner diagnostic gridfunction enum names.
    """
    return list(KASNER_NEAREST_DIAG_GFS_SUPERB)


def _kasner_exact_codegen_targets() -> Sequence[str]:
    return (
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDXXGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDXYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDXZGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDYYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDYZGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_ADDZZGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_WGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_TRKGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDXXGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDXYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDXZGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDYYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDYZGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_HDDZZGF, idx3)]",
    )


def _kasner_exponent_recovery_codegen_targets() -> Sequence[str]:
    return (
        "diagnostic_gfs[grid][IDX4pt(DIAG_PX_RECOVEREDGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_PY_RECOVEREDGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_PZ_RECOVEREDGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_PXGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_PYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_PZGF, idx3)]",
    )


def _kasner_recovered_exponent_exprs(CoordSystem: str) -> Sequence[sp.Expr]:
    hDD = ixp.declarerank2("hDD", symmetry="sym01")
    aDD = ixp.declarerank2("aDD", symmetry="sym01")
    cf = sp.Symbol("cf", real=True)
    trK = sp.Symbol("trK", real=True)
    cf_evolution = par.parval_from_str("EvolvedConformalFactor_cf")
    if cf_evolution == "phi":
        exp4phi = sp.exp(4 * cf)
    elif cf_evolution == "chi":
        exp4phi = 1 / cf
    elif cf_evolution == "W":
        exp4phi = 1 / cf**2
    else:
        raise ValueError(f"Unsupported EvolvedConformalFactor_cf = {cf_evolution}")

    rfm = refmetric.reference_metric[CoordSystem]
    gammabarDD = ixp.zerorank2()
    AbarDD = ixp.zerorank2()
    gammaDD = ixp.zerorank2()
    KDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            gammabarDD[i][j] = hDD[i][j] * rfm.ReDD[i][j] + rfm.ghatDD[i][j]
            AbarDD[i][j] = aDD[i][j] * rfm.ReDD[i][j]
            gammaDD[i][j] = exp4phi * gammabarDD[i][j]
            KDD[i][j] = exp4phi * AbarDD[i][j] + sp.Rational(1, 3) * gammaDD[i][j] * trK

    basis_transforms = BasisTransforms(CoordSystem)
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

    t_phys = sp.Symbol("KASNER_t_phys", real=True)
    p1_exact = sp.Symbol("KASNER_p1", real=True)
    p2_exact = sp.Symbol("KASNER_p2", real=True)
    p3_exact = sp.Symbol("KASNER_p3", real=True)
    p_rec = [-t_phys * KUD_cart[i][i] for i in range(3)]
    return (
        p_rec[0],
        p_rec[1],
        p_rec[2],
        p1_exact,
        p2_exact,
        p3_exact,
    )


def _kasner_exact_bssn_exprs(CoordSystem: str) -> Dict[str, object]:
    t_phys = sp.Symbol("KASNER_t_phys", real=True)
    p1 = sp.Symbol("KASNER_p1", real=True)
    p2 = sp.Symbol("KASNER_p2", real=True)
    p3 = sp.Symbol("KASNER_p3", real=True)
    gammaDD, KDD, _alpha, betaU, BU = kasner_adm_quantities(t_phys, p1, p2, p3)
    adm2bssn = ADM_to_BSSN(
        gammaDD=gammaDD,
        KDD=KDD,
        betaU=betaU,
        BU=BU,
        CoordSystem="Cartesian",
    )

    basis_transforms = BasisTransforms(CoordSystem)
    gammabarDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        adm2bssn.gammabarDD
    )
    AbarDD = basis_transforms.basis_transform_tensorDD_from_Cartesian_to_rfmbasis(
        adm2bssn.AbarDD
    )

    rfm = refmetric.reference_metric[CoordSystem]
    hDD = ixp.zerorank2()
    for i in range(3):
        for j in range(3):
            hDD[i][j] = (gammabarDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[i][j]

    return {
        "cf": adm2bssn.cf,
        "trK": adm2bssn.trK,
        "hDD": hDD,
        "AbarDD": AbarDD,
    }


def build_kasner_diagnostic_gfs_set_body() -> str:
    """
    Build the Kasner-specific C body appended to ``diagnostic_gfs_set``.

    :return: C-code body string injected into ``diagnostic_gfs_set``.
    """
    CoordSystem = par.parval_from_str("CoordSystem_to_register_CodeParameters")
    exact_bssn = _kasner_exact_bssn_exprs(CoordSystem)
    exact_AbarDD = cast(List[List[sp.Expr]], exact_bssn["AbarDD"])
    rfm = refmetric.reference_metric[CoordSystem]
    exact_aDD = [
        [exact_AbarDD[i][j] / rfm.ReDD[i][j] for j in range(3)] for i in range(3)
    ]
    exact_hDD = cast(List[List[sp.Expr]], exact_bssn["hDD"])
    exact_cf = cast(sp.Expr, exact_bssn["cf"])
    exact_trK = cast(sp.Expr, exact_bssn["trK"])
    recovered_exponent_exprs = _kasner_recovered_exponent_exprs(CoordSystem)
    body = """
    const REAL KASNER_t_phys = commondata->KASNER_t0 + commondata->time;

    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int idx3 = IDX3(i0, i1, i2);
      const REAL xx0 = griddata[grid].xx[0][i0];
      const REAL xx1 = griddata[grid].xx[1][i1];
      const REAL xx2 = griddata[grid].xx[2][i2];
      const REAL hDD00 = y_n_gfs[IDX4pt(HDD00GF, idx3)];
      const REAL hDD01 = y_n_gfs[IDX4pt(HDD01GF, idx3)];
      const REAL hDD02 = y_n_gfs[IDX4pt(HDD02GF, idx3)];
      const REAL hDD11 = y_n_gfs[IDX4pt(HDD11GF, idx3)];
      const REAL hDD12 = y_n_gfs[IDX4pt(HDD12GF, idx3)];
      const REAL hDD22 = y_n_gfs[IDX4pt(HDD22GF, idx3)];
      const REAL aDD00 = y_n_gfs[IDX4pt(ADD00GF, idx3)];
      const REAL aDD01 = y_n_gfs[IDX4pt(ADD01GF, idx3)];
      const REAL aDD02 = y_n_gfs[IDX4pt(ADD02GF, idx3)];
      const REAL aDD11 = y_n_gfs[IDX4pt(ADD11GF, idx3)];
      const REAL aDD12 = y_n_gfs[IDX4pt(ADD12GF, idx3)];
      const REAL aDD22 = y_n_gfs[IDX4pt(ADD22GF, idx3)];
      const REAL cf = y_n_gfs[IDX4pt(CFGF, idx3)];
      const REAL trK = y_n_gfs[IDX4pt(TRKGF, idx3)];
#include "set_CodeParameters.h"
"""
    body += """
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDXXGF, idx3)] = y_n_gfs[IDX4pt(ADD00GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDXYGF, idx3)] = y_n_gfs[IDX4pt(ADD01GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDXZGF, idx3)] = y_n_gfs[IDX4pt(ADD02GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDYYGF, idx3)] = y_n_gfs[IDX4pt(ADD11GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDYZGF, idx3)] = y_n_gfs[IDX4pt(ADD12GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_ADDZZGF, idx3)] = y_n_gfs[IDX4pt(ADD22GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_TRKGF, idx3)] = y_n_gfs[IDX4pt(TRKGF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDXXGF, idx3)] = y_n_gfs[IDX4pt(HDD00GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDXYGF, idx3)] = y_n_gfs[IDX4pt(HDD01GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDXZGF, idx3)] = y_n_gfs[IDX4pt(HDD02GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDYYGF, idx3)] = y_n_gfs[IDX4pt(HDD11GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDYZGF, idx3)] = y_n_gfs[IDX4pt(HDD12GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_HDDZZGF, idx3)] = y_n_gfs[IDX4pt(HDD22GF, idx3)];
"""
    body += "      {\n      " + ccg.c_codegen(
        [
            exact_aDD[0][0],
            exact_aDD[0][1],
            exact_aDD[0][2],
            exact_aDD[1][1],
            exact_aDD[1][2],
            exact_aDD[2][2],
            exact_cf,
            exact_trK,
            exact_hDD[0][0],
            exact_hDD[0][1],
            exact_hDD[0][2],
            exact_hDD[1][1],
            exact_hDD[1][2],
            exact_hDD[2][2],
        ],
        list(_kasner_exact_codegen_targets()),
        verbose=False,
        include_braces=False,
    ).replace("\n", "\n      ")
    body += "\n      }\n"
    body += "      {\n      " + ccg.c_codegen(
        list(recovered_exponent_exprs),
        list(_kasner_exponent_recovery_codegen_targets()),
        verbose=False,
        include_braces=False,
    ).replace("\n", "\n      ")
    body += "\n      }\n"
    body += """
    } // END LOOP over all gridpoints to set Kasner diagnostics
"""
    return body


def _replace_cfunc_all(name: str, old: str, new: str) -> None:
    cfunc = cfc.CFunction_dict[name]
    if (
        old not in cfunc.body
        or old not in cfunc.raw_function
        or old not in cfunc.full_function
    ):
        raise ValueError(f"Expected snippet not found in CFunction '{name}': {old}")
    cfunc.body = cfunc.body.replace(old, new, 1)
    cfunc.raw_function = cfunc.raw_function.replace(old, new, 1)
    cfunc.full_function = cfunc.full_function.replace(old, new, 1)


def _insert_before_marker_all(name: str, marker: str, insert_text: str) -> None:
    cfunc = cfc.CFunction_dict[name]
    cfunc.body = cfunc.body.replace(marker, insert_text + marker)
    cfunc.raw_function = cfunc.raw_function.replace(marker, insert_text + marker)
    cfunc.full_function = cfunc.full_function.replace(marker, insert_text + marker)


def register_CFunction_diagnostic_gfs_set(
    enable_interp_diagnostics: bool = False,
    enable_psi4: bool = False,
    enable_T4munu: bool = False,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register ``diagnostic_gfs_set`` and inject Kasner-specific diagnostic channels.

    :param enable_interp_diagnostics: Whether interpolation diagnostics are enabled.
    :param enable_psi4: Whether Psi4 diagnostics are enabled.
    :param enable_T4munu: Whether stress-energy diagnostics are enabled.
    :return: Parallel-codegen registration token or ``None`` during registration phase.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    register_kasner_diag_gridfunctions()
    gr_diagnostic_gfs_set.register_CFunction_diagnostic_gfs_set(
        enable_interp_diagnostics=enable_interp_diagnostics,
        enable_psi4=enable_psi4,
        enable_T4munu=enable_T4munu,
    )
    cfunc = cfc.CFunction_dict["diagnostic_gfs_set"]
    kasner_body = build_kasner_diagnostic_gfs_set_body()
    grid_loop_end = "  } // END LOOP over grids\n"
    if kasner_body not in cfunc.body:
        _insert_before_marker_all("diagnostic_gfs_set", grid_loop_end, kasner_body)
    return pcg.NRPyEnv()


def register_CFunction_diagnostics_nearest() -> Union[None, pcg.NRPyEnv_type]:
    """
    Register ``diagnostics_nearest`` with Kasner diagnostic channels in 0D/1D/2D output lists.

    :return: Parallel-codegen registration token or ``None`` during registration phase.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None
    gr_diagnostics_nearest.register_CFunction_diagnostics_nearest()
    selected = ", ".join(kasner_nearest_diag_gf_names_bhah())
    _replace_cfunc_all(
        "diagnostics_nearest",
        "const int which_gfs_0d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};",
        f"const int which_gfs_0d[] = {{{selected}}};",
    )
    _replace_cfunc_all(
        "diagnostics_nearest",
        "const int which_gfs_1d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};",
        f"const int which_gfs_1d[] = {{{selected}}};",
    )
    _replace_cfunc_all(
        "diagnostics_nearest",
        "const int which_gfs_2d[] = {DIAG_HAMILTONIANGF, DIAG_MSQUAREDGF};",
        f"const int which_gfs_2d[] = {{{selected}}};",
    )
    return pcg.NRPyEnv()

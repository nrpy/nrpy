"""Kasner-specific diagnostics helpers for BHaH infrastructures."""

from __future__ import annotations

from typing import List, Sequence, Tuple

import nrpy.c_codegen as ccg
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures.BHaH.kasner_solution.kasner_exact import (
    kasner_exact_bssn_exprs,
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
    ("DIAG_LAMBDAUX", "lambdaU0"),
    ("DIAG_LAMBDAUY", "lambdaU1"),
    ("DIAG_LAMBDAUZ", "lambdaU2"),
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
    ("DIAG_EXACT_LAMBDAUX", "EXACT_lambdaU0"),
    ("DIAG_EXACT_LAMBDAUY", "EXACT_lambdaU1"),
    ("DIAG_EXACT_LAMBDAUZ", "EXACT_lambdaU2"),
)

KASNER_NEAREST_DIAG_GFS_BHAH: Tuple[str, ...] = tuple(f"{name}GF" for name, _ in KASNER_DIAG_GRIDFUNCTIONS)
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
        gri.register_gridfunctions(names=name, desc=desc, group="DIAG")


def kasner_nearest_diag_gf_names_bhah() -> List[str]:
    """Return Kasner diagnostic GF enum names used by non-superB nearest diagnostics."""
    return list(KASNER_NEAREST_DIAG_GFS_BHAH)


def kasner_nearest_diag_gf_names_superb() -> List[str]:
    """Return Kasner diagnostic GF enum names used by superB nearest diagnostics."""
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
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_LAMBDAUXGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_LAMBDAUYGF, idx3)]",
        "diagnostic_gfs[grid][IDX4pt(DIAG_EXACT_LAMBDAUZGF, idx3)]",
    )


def build_kasner_diagnostic_gfs_set_body() -> str:
    """Build the Kasner-specific C body appended to ``diagnostic_gfs_set``."""
    CoordSystem = par.parval_from_str("CoordSystem_to_register_CodeParameters")
    exact_bssn = kasner_exact_bssn_exprs(CoordSystem)
    exact_AbarDD = exact_bssn["AbarDD"]
    exact_hDD = exact_bssn["hDD"]
    exact_lambdaU = exact_bssn["lambdaU"]

    ghat_reads = ""
    if CoordSystem.startswith("GeneralRFM"):
        ghat_reads = """
      const REAL ghatDD00 = auxevol_gfs[IDX4pt(GHATDD00GF, idx3)];
      const REAL ghatDD01 = auxevol_gfs[IDX4pt(GHATDD01GF, idx3)];
      const REAL ghatDD02 = auxevol_gfs[IDX4pt(GHATDD02GF, idx3)];
      const REAL ghatDD11 = auxevol_gfs[IDX4pt(GHATDD11GF, idx3)];
      const REAL ghatDD12 = auxevol_gfs[IDX4pt(GHATDD12GF, idx3)];
      const REAL ghatDD22 = auxevol_gfs[IDX4pt(GHATDD22GF, idx3)];
      const REAL ghatDDdD000 = auxevol_gfs[IDX4pt(GHATDDDD000GF, idx3)];
      const REAL ghatDDdD001 = auxevol_gfs[IDX4pt(GHATDDDD001GF, idx3)];
      const REAL ghatDDdD002 = auxevol_gfs[IDX4pt(GHATDDDD002GF, idx3)];
      const REAL ghatDDdD010 = auxevol_gfs[IDX4pt(GHATDDDD010GF, idx3)];
      const REAL ghatDDdD011 = auxevol_gfs[IDX4pt(GHATDDDD011GF, idx3)];
      const REAL ghatDDdD012 = auxevol_gfs[IDX4pt(GHATDDDD012GF, idx3)];
      const REAL ghatDDdD020 = auxevol_gfs[IDX4pt(GHATDDDD020GF, idx3)];
      const REAL ghatDDdD021 = auxevol_gfs[IDX4pt(GHATDDDD021GF, idx3)];
      const REAL ghatDDdD022 = auxevol_gfs[IDX4pt(GHATDDDD022GF, idx3)];
      const REAL ghatDDdD110 = auxevol_gfs[IDX4pt(GHATDDDD110GF, idx3)];
      const REAL ghatDDdD111 = auxevol_gfs[IDX4pt(GHATDDDD111GF, idx3)];
      const REAL ghatDDdD112 = auxevol_gfs[IDX4pt(GHATDDDD112GF, idx3)];
      const REAL ghatDDdD120 = auxevol_gfs[IDX4pt(GHATDDDD120GF, idx3)];
      const REAL ghatDDdD121 = auxevol_gfs[IDX4pt(GHATDDDD121GF, idx3)];
      const REAL ghatDDdD122 = auxevol_gfs[IDX4pt(GHATDDDD122GF, idx3)];
      const REAL ghatDDdD220 = auxevol_gfs[IDX4pt(GHATDDDD220GF, idx3)];
      const REAL ghatDDdD221 = auxevol_gfs[IDX4pt(GHATDDDD221GF, idx3)];
      const REAL ghatDDdD222 = auxevol_gfs[IDX4pt(GHATDDDD222GF, idx3)];
"""

    body = """
    SET_NXX_PLUS_2NGHOSTS_VARS(grid);
    const REAL KASNER_t_phys = commondata->KASNER_t0 + commondata->time;
    const REAL KASNER_p1 = commondata->KASNER_p1;
    const REAL KASNER_p2 = commondata->KASNER_p2;
    const REAL KASNER_p3 = commondata->KASNER_p3;

    LOOP_OMP("omp parallel for", i0, 0, Nxx_plus_2NGHOSTS0, i1, 0, Nxx_plus_2NGHOSTS1, i2, 0, Nxx_plus_2NGHOSTS2) {
      const int idx3 = IDX3(i0, i1, i2);
      const REAL xx0 = griddata[grid].xx[0][i0];
      const REAL xx1 = griddata[grid].xx[1][i1];
      const REAL xx2 = griddata[grid].xx[2][i2];
#include "set_CodeParameters.h"
"""
    body += ghat_reads
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
      diagnostic_gfs[grid][IDX4pt(DIAG_LAMBDAUXGF, idx3)] = y_n_gfs[IDX4pt(LAMBDAU0GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_LAMBDAUYGF, idx3)] = y_n_gfs[IDX4pt(LAMBDAU1GF, idx3)];
      diagnostic_gfs[grid][IDX4pt(DIAG_LAMBDAUZGF, idx3)] = y_n_gfs[IDX4pt(LAMBDAU2GF, idx3)];
"""
    body += "      " + ccg.c_codegen(
        [
            exact_AbarDD[0][0],
            exact_AbarDD[0][1],
            exact_AbarDD[0][2],
            exact_AbarDD[1][1],
            exact_AbarDD[1][2],
            exact_AbarDD[2][2],
            exact_bssn["cf"],
            exact_bssn["trK"],
            exact_hDD[0][0],
            exact_hDD[0][1],
            exact_hDD[0][2],
            exact_hDD[1][1],
            exact_hDD[1][2],
            exact_hDD[2][2],
            exact_lambdaU[0],
            exact_lambdaU[1],
            exact_lambdaU[2],
        ],
        list(_kasner_exact_codegen_targets()),
        verbose=False,
        include_braces=False,
    ).replace("\n", "\n      ")
    body += """
    } // END LOOP over all gridpoints to set Kasner diagnostics
"""
    return body

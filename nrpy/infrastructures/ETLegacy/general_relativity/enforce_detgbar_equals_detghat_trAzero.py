"""
Generate function enforcing conformal metric determinant and A trace constraints.

Authors: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
         Samuel Cupp
"""

from inspect import currentframe as cfr
from types import FrameType as FT
from typing import List, Union, cast

import sympy as sp

import nrpy.c_codegen as ccg
import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.helpers.parallel_codegen as pcg
import nrpy.indexedexp as ixp
import nrpy.infrastructures.ETLegacy.simple_loop as lp
import nrpy.reference_metric as refmetric  # NRPy: Reference metric support
from nrpy.equations.general_relativity.BSSN_quantities import BSSN_quantities
from nrpy.infrastructures.ETLegacy.ETLegacy_include_header import (
    define_standard_includes,
)


def register_CFunction_enforce_detgbar_equals_detghat_trAzero(
    thorn_name: str,
    CoordSystem: str,
    enable_rfm_precompute: bool,
    fd_orders: List[int],
    OMP_collapse: int = 1,
) -> Union[None, pcg.NRPyEnv_type]:
    """
    Register combined determinant and trace-free conformal-tensor projection.

    :param thorn_name: The Einstein Toolkit thorn name.
    :param CoordSystem: The coordinate system to be used.
    :param enable_rfm_precompute: Whether to enable reference metric precomputation.
    :param fd_orders: Finite-difference orders of registered ADM-to-BSSN producers.
    :param OMP_collapse: Degree of OpenMP loop collapsing.

    :return: None if in registration phase, else the updated NRPy environment.
    """
    if pcg.pcg_registration_phase():
        pcg.register_func_call(f"{__name__}.{cast(FT, cfr()).f_code.co_name}", locals())
        return None

    Bq = BSSN_quantities[
        CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
    ]
    is_general_rfm = CoordSystem.startswith("GeneralRFM")

    desc = r"""Enforce det(gammabar) = det(gammahat) and tr(Abar) = 0."""
    name = f"{thorn_name}_enforce_detgbar_equals_detghat_trAzero"
    body = f"""  DECLARE_CCTK_ARGUMENTS_{name};
  DECLARE_CCTK_PARAMETERS;

"""

    hprimeDD = ixp.zerorank2()
    aprimeDD = ixp.zerorank2()
    if is_general_rfm:
        rfm = refmetric.reference_metric[
            CoordSystem + "_rfm_precompute" if enable_rfm_precompute else CoordSystem
        ]
        _, detgammabar = ixp.symm_matrix_inverter3x3(Bq.gammabarDD)
        nrpyAbs = sp.Function("nrpyAbs")
        q = (nrpyAbs(rfm.detgammahat) / detgammabar) ** sp.Rational(1, 3)
        gammabarprimeDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                gammabarprimeDD[i][j] = q * Bq.gammabarDD[i][j]
        gammabarprimeUU, _ = ixp.symm_matrix_inverter3x3(gammabarprimeDD)
        trAbarprime = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                trAbarprime += gammabarprimeUU[i][j] * Bq.AbarDD[i][j]
        for i in range(3):
            for j in range(3):
                hprimeDD[i][j] = (gammabarprimeDD[i][j] - rfm.ghatDD[i][j]) / rfm.ReDD[
                    i
                ][j]
                aprimeDD[i][j] = (
                    Bq.AbarDD[i][j]
                    - sp.Rational(1, 3) * gammabarprimeDD[i][j] * trAbarprime
                ) / rfm.ReDD[i][j]
    else:
        # Every standard RFM is orthogonal. With diagonal scale-factor matrix D,
        # gammabar = D H D, gammahat = D I D, and Abar = D a D, where H = I+h.
        # Thus det(gammabar)/det(gammahat) = det(H), so q = det(H)^(-1/3)
        # and h' = q H - I. Also gammabar'^(-1) = q^(-1) D^(-1) H^(-1)
        # D^(-1), hence tr(gammabar'^(-1) Abar) = q^(-1) (H^(-1):a).
        # The q factors cancel in gammabar' tr/3, leaving a' = a-H(H^(-1):a)/3.
        HDD = ixp.zerorank2()
        for i in range(3):
            for j in range(3):
                HDD[i][j] = Bq.hDD[i][j] + sp.KroneckerDelta(i, j)
        HUU, detH = ixp.symm_matrix_inverter3x3(HDD)
        q = detH ** sp.Rational(-1, 3)
        trA = sp.sympify(0)
        for i in range(3):
            for j in range(3):
                trA += HUU[i][j] * Bq.aDD[i][j]
        for i in range(3):
            for j in range(3):
                hprimeDD[i][j] = q * HDD[i][j] - sp.KroneckerDelta(i, j)
                aprimeDD[i][j] = Bq.aDD[i][j] - sp.Rational(1, 3) * HDD[i][j] * trA

    access_gfs: List[str] = []
    projected_exprs: List[sp.Expr] = []
    for basename, expressions in (("hDD", hprimeDD), ("aDD", aprimeDD)):
        for i in range(3):
            for j in range(i, 3):
                access_gfs.append(
                    gri.ETLegacyGridFunction.access_gf(gf_name=f"{basename}{i}{j}")
                )
                projected_exprs.append(expressions[i][j])

    # To evaluate the cube root, SIMD support requires e.g., SLEEF.
    #   Also need to be careful to not access memory out of bounds!
    #   After all this is a loop over ALL POINTS.
    #   Exercise to the reader: prove that for any reasonable grid,
    #   SIMD loops over grid interiors never write data out of bounds
    #   and are threadsafe for any reasonable number of threads.
    body += lp.simple_loop(
        loop_body=ccg.c_codegen(
            projected_exprs,
            access_gfs,
            enable_simd=False,
            automatically_read_gf_data_from_memory=True,
            enable_fd_codegen=True,
            enable_fd_functions=True,
        ),
        loop_region="all points",
        enable_simd=False,
        OMP_collapse=OMP_collapse,
    )

    reads = """hDD00GF, hDD01GF, hDD02GF, hDD11GF, hDD12GF, hDD22GF,
          aDD00GF, aDD01GF, aDD02GF, aDD11GF, aDD12GF, aDD22GF"""
    writes = """hDD00GF(everywhere), hDD01GF(everywhere), hDD02GF(everywhere),
          hDD11GF(everywhere), hDD12GF(everywhere), hDD22GF(everywhere),
          aDD00GF(everywhere), aDD01GF(everywhere), aDD02GF(everywhere),
          aDD11GF(everywhere), aDD12GF(everywhere), aDD22GF(everywhere)"""
    initial_producer_names = [
        f"{thorn_name}_ADM_to_BSSN_order_{fd_order}" for fd_order in fd_orders
    ]
    # Cactus requires a parenthesized whitespace-separated list when an AFTER
    # clause names multiple schedule items. A bare list fails CST parsing.
    initial_producers = " ".join(initial_producer_names)
    if len(initial_producer_names) > 1:
        initial_producers = f"({initial_producers})"
    schedule_initial = (
        "CCTK_INITIAL",
        f"""
schedule FUNC_NAME in CCTK_INITIAL after {initial_producers}
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project freshly initialized or restored BSSN/fCCZ4 state"
""",
    )
    schedule_poststep = (
        "MoL_PostStep",
        f"""
schedule FUNC_NAME in MoL_PostStep after {thorn_name}_evol_ApplyBCs
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project every updated BSSN/fCCZ4 state before its next RHS"
""",
    )
    schedule_pseudoevolution = (
        "MoL_PseudoEvolution",
        f"""
schedule FUNC_NAME in MoL_PseudoEvolution before {thorn_name}_BSSN_constraints
{{
  LANG: C
  READS:  {reads}
  WRITES: {writes}
}} "Project BSSN/fCCZ4 state before pseudo-evolution consumers"
""",
    )

    cfc.register_CFunction(
        subdirectory=thorn_name,
        includes=define_standard_includes(),
        desc=desc,
        cfunc_type="void",
        name=name,
        params="CCTK_ARGUMENTS",
        body=body,
        ET_thorn_name=thorn_name,
        ET_schedule_bins_entries=[
            schedule_initial,
            schedule_poststep,
            schedule_pseudoevolution,
        ],
    )
    return pcg.NRPyEnv()


if __name__ == "__main__":
    import nrpy.params as par
    import nrpy.validate_expressions.validate_expressions as ve

    validation_HDD = [
        [sp.Integer(4), sp.Integer(1), sp.Rational(1, 2)],
        [sp.Integer(1), sp.Integer(3), sp.Rational(1, 4)],
        [sp.Rational(1, 2), sp.Rational(1, 4), sp.Integer(2)],
    ]
    validation_aDD = [
        [sp.Integer(2), sp.Integer(3), sp.Integer(5)],
        [sp.Integer(3), sp.Integer(7), sp.Integer(11)],
        [sp.Integer(5), sp.Integer(11), sp.Integer(13)],
    ]
    validation_HUU, validation_detH = ixp.symm_matrix_inverter3x3(validation_HDD)
    validation_q = validation_detH ** sp.Rational(-1, 3)
    validation_HprimeDD = [
        [validation_q * validation_HDD[i][j] for j in range(3)] for i in range(3)
    ]
    validation_trA = sum(
        validation_HUU[i][j] * validation_aDD[i][j] for i in range(3) for j in range(3)
    )
    validation_HprimeUU, validation_detHprime = ixp.symm_matrix_inverter3x3(
        validation_HprimeDD
    )
    validation_aprimeDD = [
        [
            validation_aDD[i][j]
            - sp.Rational(1, 3) * validation_HDD[i][j] * validation_trA
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_trAprime = sum(
        validation_HprimeUU[i][j] * validation_aprimeDD[i][j]
        for i in range(3)
        for j in range(3)
    )
    ve.assert_equal(validation_detHprime, sp.sympify(1), suppress_message=True)
    ve.assert_equal(validation_trAprime, sp.sympify(0), suppress_message=True)

    validation_scaleD = [sp.Integer(2), sp.Integer(3), sp.Integer(5)]
    validation_gammahatDD = [
        [
            sp.KroneckerDelta(i, j) * validation_scaleD[i] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_gammabarDD = [
        [
            validation_scaleD[i] * validation_HDD[i][j] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_AbarDD = [
        [
            validation_scaleD[i] * validation_aDD[i][j] * validation_scaleD[j]
            for j in range(3)
        ]
        for i in range(3)
    ]
    _, validation_detgammahat = ixp.symm_matrix_inverter3x3(validation_gammahatDD)
    _, validation_detgammabar = ixp.symm_matrix_inverter3x3(validation_gammabarDD)
    validation_full_q = (
        validation_detgammahat / validation_detgammabar
    ) ** sp.Rational(1, 3)
    validation_gammabarprimeDD = [
        [validation_full_q * validation_gammabarDD[i][j] for j in range(3)]
        for i in range(3)
    ]
    validation_gammabarprimeUU = ixp.symm_matrix_inverter3x3(
        validation_gammabarprimeDD
    )[0]
    validation_full_trA = sum(
        validation_gammabarprimeUU[i][j] * validation_AbarDD[i][j]
        for i in range(3)
        for j in range(3)
    )
    validation_full_hprimeDD = [
        [
            (validation_gammabarprimeDD[i][j] - validation_gammahatDD[i][j])
            / (validation_scaleD[i] * validation_scaleD[j])
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_full_aprimeDD = [
        [
            (
                validation_AbarDD[i][j]
                - sp.Rational(1, 3)
                * validation_gammabarprimeDD[i][j]
                * validation_full_trA
            )
            / (validation_scaleD[i] * validation_scaleD[j])
            for j in range(3)
        ]
        for i in range(3)
    ]
    validation_reduced_hprimeDD = [
        [validation_HprimeDD[i][j] - sp.KroneckerDelta(i, j) for j in range(3)]
        for i in range(3)
    ]
    ve.assert_equal(
        {
            "hprimeDD": validation_full_hprimeDD,
            "aprimeDD": validation_full_aprimeDD,
        },
        {
            "hprimeDD": validation_reduced_hprimeDD,
            "aprimeDD": validation_aprimeDD,
        },
        suppress_message=True,
    )

    par.set_parval_from_str("Infrastructure", "ETLegacy")
    par.set_parval_from_str("enable_parallel_codegen", False)
    for validation_coord, validation_precompute in (
        ("Cartesian", False),
        ("SinhSpherical", True),
        ("GeneralRFM", True),
    ):
        cfc.CFunction_dict.clear()
        register_CFunction_enforce_detgbar_equals_detghat_trAzero(
            "Validation", validation_coord, validation_precompute, [2, 4]
        )
        generated = cfc.CFunction_dict[
            "Validation_enforce_detgbar_equals_detghat_trAzero"
        ].body
        if validation_coord.startswith("GeneralRFM"):
            for validation_i, validation_j in (
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 1),
                (1, 2),
                (2, 2),
            ):
                if (
                    f"const CCTK_REAL ghatDD{validation_i}{validation_j} ="
                    not in generated
                ):
                    raise AssertionError(
                        "GeneralRFM projection must read all six ghat components"
                    )
        elif "ghatDD" in generated or "_of_xx" in generated:
            raise AssertionError(
                "standard projection must not read reference-metric data"
            )
        validation_schedules = cfc.CFunction_dict[
            "Validation_enforce_detgbar_equals_detghat_trAzero"
        ].ET_schedule_bins_entries
        if validation_schedules is None:
            raise AssertionError("projection schedules must be registered")
        initial_schedule = dict(validation_schedules)["CCTK_INITIAL"]
        if (
            "after (Validation_ADM_to_BSSN_order_2 Validation_ADM_to_BSSN_order_4)"
            not in initial_schedule
            or "Validation_evol_ApplyBCs" in initial_schedule
        ):
            raise AssertionError(
                "initial projection must follow each exact ADM-to-BSSN producer"
            )
        poststep_schedule = dict(validation_schedules)["MoL_PostStep"]
        if "after Validation_evol_ApplyBCs" not in poststep_schedule:
            raise AssertionError(
                "poststep projection must follow evolved-variable boundary conditions"
            )
        pseudoevolution_schedule = dict(validation_schedules)["MoL_PseudoEvolution"]
        if (
            "before Validation_BSSN_constraints" not in pseudoevolution_schedule
            or "after Validation_aux_ApplyBCs" in pseudoevolution_schedule
        ):
            raise AssertionError(
                "projection must precede constraints without an auxiliary-boundary cycle"
            )
        if generated.count("Step 1 of 2") != 1 or generated.count("Step 2 of 2") != 1:
            raise AssertionError(
                "projection must contain one generated twelve-output block"
            )
        if generated.count("END LOOP: for i2 over") != 1:
            raise AssertionError("projection must contain one all-points loop nest")
        write_section = generated.split("Step 2 of 2", maxsplit=1)[1]
        if write_section.count("GF[CCTK_GFINDEX3D") != 12:
            raise AssertionError("projection must store twelve outputs in one block")
        first_store = generated.index(
            "GF[CCTK_GFINDEX3D", generated.index("Step 2 of 2")
        )
        for validation_basename in ("aDD", "hDD"):
            for validation_i, validation_j in (
                (0, 0),
                (0, 1),
                (0, 2),
                (1, 1),
                (1, 2),
                (2, 2),
            ):
                load = f"const CCTK_REAL {validation_basename}{validation_i}{validation_j} ="
                if load not in generated or generated.index(load) > first_store:
                    raise AssertionError(
                        "projection must load all raw inputs before stores"
                    )

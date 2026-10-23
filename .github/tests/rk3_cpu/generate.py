"""Generate the CPU BSSN, radiation, projection, and complete RK3 functions for CI."""

import argparse
from pathlib import Path

import nrpy.c_function as cfc
import nrpy.grid as gri
import nrpy.params as par
from nrpy.infrastructures import BHaH


def generate(output: Path, simd: bool) -> None:
    """
    Emit a standalone CPU integration program's functions and canonical headers.

    :param output: Owned directory for generated C and headers.
    :param simd: Whether RHS, Ricci, and RK updates use explicit SIMD.
    """
    output.mkdir(parents=True, exist_ok=True)
    par.set_parval_from_str("parallelization", "openmp")
    par.set_parval_from_str("fd_order", 4)
    par.set_parval_from_str("enable_parallel_codegen", False)
    par.set_parval_from_str("EvolvedConformalFactor_cf", "W")
    for axis in range(3):
        par.register_CodeParameter(
            "int",
            "rk3_ci",
            f"Nxx_plus_2NGHOSTS{axis}",
            add_to_parfile=False,
            add_to_set_CodeParameters_h=True,
        )
        par.register_CodeParameter(
            "REAL",
            "rk3_ci",
            f"invdxx{axis}",
            add_to_parfile=False,
            add_to_set_CodeParameters_h=True,
        )
    BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval(
        CoordSystem="Cartesian",
        enable_rfm_precompute=True,
        enable_RbarDD_gridfunctions=True,
        enable_T4munu=False,
        enable_intrinsics=simd,
        enable_fd_functions=True,
        LapseEvolutionOption="OnePlusLog",
        ShiftEvolutionOption="GammaDriving2ndOrder_Covariant",
        enable_KreissOliger_dissipation=True,
        enable_RK3_accumulation=True,
        enable_cpu_tiling=True,
    )
    BHaH.general_relativity.Ricci_eval.register_CFunction_Ricci_eval(
        CoordSystem="Cartesian",
        enable_intrinsics=simd,
        enable_fd_functions=True,
        OMP_collapse=1,
        enable_cpu_tiling=True,
    )
    BHaH.general_relativity.rhs_eval.register_CFunction_rhs_eval_with_Ricci("Cartesian")
    BHaH.general_relativity.enforce_detgbar_equals_detghat_trAzero.register_CFunction_enforce_detgbar_equals_detghat_trAzero(
        CoordSystem="Cartesian",
        enable_rfm_precompute=True,
        enable_fd_functions=True,
        OMP_collapse=1,
    )
    bc = BHaH.CurviBoundaryConditions
    bc.apply_bcs_inner_only.register_CFunction_apply_bcs_inner_only()
    bc.apply_bcs_outerextrap_and_inner.register_CFunction_apply_bcs_outerextrap_and_inner()
    bc.apply_bcs_outerradiation_and_inner.register_CFunction_apply_bcs_outerradiation_and_inner(
        "Cartesian", radiation_BC_fd_order=4
    )
    bc.apply_bcs_outerradiation_and_inner.register_CFunction_apply_bcs_outerradiation_and_inner_accumulate(
        "Cartesian", radiation_BC_fd_order=4
    )
    bc.BHaH_defines.register_BHaH_defines_h()
    BHaH.griddata_commondata.register_griddata_commondata(
        "rk3_ci", "bc_struct bcstruct", "inner and radiation boundary points"
    )
    par.register_CodeParameter(
        "char[50]", "rk3_ci", "outer_bc_type", "radiation", commondata=True
    )
    rhs, accumulation = BHaH.general_relativity.rhs_eval.MoL_rhs_strings(
        pre_rhs_string="Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);"
    )
    BHaH.MoLtimestepping.register_all.register_CFunctions(
        MoL_method="RK3 Ralston", register_MoL_step_forward_in_time=False
    )
    tables = BHaH.MoLtimestepping.rk_butcher_table_dictionary.generate_Butcher_tables()
    for method in ("Ralston", "Heun"):
        ordinary = rhs.replace(
            "Ricci_eval(params, rfmstruct, RK_INPUT_GFS, auxevol_gfs);", ""
        ).replace("rhs_eval(", "rhs_eval_with_Ricci(")
        post = """
if (strncmp(commondata->outer_bc_type, "extrapolation", 50) == 0)
  apply_bcs_outerextrap_and_inner(commondata, params, bcstruct, RK_OUTPUT_GFS);
else
  apply_bcs_inner_only(commondata, params, bcstruct, RK_OUTPUT_GFS);
enforce_detgbar_equals_detghat_trAzero(params, rfmstruct, RK_OUTPUT_GFS, auxevol_gfs);
"""
        for name in (
            "rhs_eval_with_Ricci",
            "rhs_eval_accumulate",
            "Ricci_eval",
            "apply_bcs_outerradiation_and_inner_accumulate",
            "apply_bcs_outerradiation_and_inner",
            "enforce_detgbar_equals_detghat_trAzero",
        ):
            ordinary = ordinary.replace(name + "(", name + "__rfm__Cartesian(")
            accumulation = accumulation.replace(name + "(", name + "__rfm__Cartesian(")
            post = post.replace(name + "(", name + "__rfm__Cartesian(")
        BHaH.MoLtimestepping.MoL_step_forward_in_time.register_CFunction_MoL_step_forward_in_time(
            Butcher_dict=tables,
            MoL_method="RK3 " + method,
            rhs_string=ordinary,
            rhs_accumulation_string=accumulation,
            post_rhs_string=post,
            enable_rfm_precompute=True,
            enable_curviBCs=True,
            enable_intrinsics=simd,
        )
        # Each complete driver has private substep functions in its own C file.
        driver = cfc.CFunction_dict.pop("MoL_step_forward_in_time")
        (output / ("MoL_" + method + ".c")).write_text(
            driver.full_function.replace(
                "MoL_step_forward_in_time(", "MoL_" + method + "("
            ),
            encoding="utf-8",
        )
        # Both RK3 tableaux have exactly these two intermediate evolved arrays.
        stages = BHaH.MoLtimestepping.rk_butcher_table_dictionary.intermediate_stage_gf_names_list(
            BHaH.MoLtimestepping.rk_butcher_table_dictionary.generate_Butcher_tables(),
            MoL_method="RK3 " + method,
        )
        assert stages == [
            "k1_or_y_nplus_a21_k1_or_y_nplus1_running_total_gfs",
            "k2_or_y_nplus_a32_k2_gfs",
        ], stages
    BHaH.BHaH_defines_h.register_BHaH_defines(
        "reference_metric", "typedef struct {int unused;} rfm_struct;\n"
    )
    BHaH.BHaH_defines_h.output_BHaH_defines_h(
        project_dir=str(output),
        additional_includes=["omp.h"],
        fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
    )
    BHaH.CodeParameters.write_CodeParameters_h_files(str(output))
    prototypes = '#include "BHaH_defines.h"\n'
    for name, function in cfc.CFunction_dict.items():
        (output / (name + ".c")).write_text(function.full_function, encoding="utf-8")
        prototypes += function.function_prototype + "\n"
    prototypes += "void MoL_Ralston(commondata_struct *, griddata_struct *);\n"
    prototypes += "void MoL_Heun(commondata_struct *, griddata_struct *);\n"
    (output / "BHaH_function_prototypes.h").write_text(prototypes, encoding="utf-8")
    # This is the same intrinsic header copied into generated BHaH projects.
    intrinsics = output / "intrinsics"
    intrinsics.mkdir(exist_ok=True)
    (intrinsics / "simd_intrinsics.h").write_bytes(
        (Path(gri.__file__).parent / "helpers" / "simd_intrinsics.h").read_bytes()
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--simd", action="store_true")
    args = parser.parse_args()
    generate(args.output, args.simd)

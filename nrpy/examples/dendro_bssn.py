"""
Generate the complete NRPy BSSN application beside Dendro-GR's BSSN_GR.

Run as ``python -m nrpy.examples.dendro_bssn``.

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com
"""

import argparse
import tempfile
from pathlib import Path
from typing import Dict

import nrpy.c_function as cfc
import nrpy.helpers.parallel_codegen as pcg
import nrpy.params as par
from nrpy.helpers.conditional_file_updater import ConditionalFileUpdater
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH.general_relativity import (
    ADM_Initial_Data_Reader__BSSN_Converter,
)
from nrpy.infrastructures.BHaH.general_relativity.TwoPunctures import (
    ID_persist_struct,
    TwoPunctures_lib,
)
from nrpy.infrastructures.Dendro import (
    CMakeLists,
    CodeParameters,
    Dendro_defines_h,
    checkpoint,
    constants_h,
    main_cpp,
    param_toml,
    solver_context,
    state_h,
    types_h,
)
from nrpy.infrastructures.Dendro.general_relativity import (
    ADM_to_BSSN,
    BSSN_constraints,
    BSSN_to_ADM,
    Ricci_eval,
    adm_quantities,
    adm_quantities_surface_data,
    apparent_horizon,
    diagnostics,
    enforce_detgbar_equals_detghat_trAzero,
    floor_the_lapse_and_conformal_factor,
    gravitational_waves,
    initial_data_lambdaU,
    physical_boundary,
    physical_boundary_ghosts,
    psi4_eval,
    rhs_eval,
    twopunctures,
)

SOLVER_NAME = "NRPy_BSSN_GR"
SOLVER_PREFIX = "BSSN"
SOLVER_STEM = "bssn"
SOLVER_NAMESPACE = "nrpy::bssn"
EXECUTABLE_NAME = "nrpyBssnSolver"
PROFILE_NAME = "bssn_twopunctures"
COORD_SYSTEM = "Cartesian"
LAPSE_EVOLUTION_OPTION = "OnePlusLog"
SHIFT_EVOLUTION_OPTION = "GammaDriving2ndOrder_Covariant__Hatted"
FD_ORDERS = (4, 6, 8)


def parse_args() -> argparse.Namespace:
    """
    Parse the generation options.

    :return: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Generate the complete NRPy BSSN application for Dendro-GR"
    )
    parser.add_argument(
        "--project-dir",
        type=Path,
        help="Dendro-GR root; by default, locate the nearest enclosing checkout",
    )
    parser.add_argument(
        "--fd-order",
        type=int,
        choices=FD_ORDERS,
        default=6,
        help="runtime finite-difference profile; all production orders are emitted",
    )
    parser.add_argument(
        "--conformal-factor",
        choices=("W", "chi"),
        default="W",
        help="evolved conformal factor (default: W)",
    )
    parser.add_argument("--ko", dest="ko", action="store_true")
    parser.add_argument("--no-ko", dest="ko", action="store_false")
    parser.set_defaults(ko=True)
    return parser.parse_args()


def main() -> None:
    """
    Generate the complete Dendro BSSN application.

    :raises ValueError: If no enclosing Dendro-GR root can be found.
    """
    args = parse_args()
    if args.project_dir is None:
        project_dir = next(
            (
                candidate
                for candidate in (Path.cwd(), *Path.cwd().parents)
                if (candidate / "CMakeLists.txt").is_file()
                and (candidate / "BSSN_GR").is_dir()
            ),
            None,
        )
        if project_dir is None:
            raise ValueError(
                "No enclosing Dendro-GR root contains CMakeLists.txt and BSSN_GR/. "
                "Pass --project-dir explicitly."
            )
    else:
        project_dir = args.project_dir.resolve()

    par.set_parval_from_str("Infrastructure", "Dendro")
    par.set_parval_from_str("fp_type", "double")
    par.set_parval_from_str("parallelization", "none")
    par.set_parval_from_str("fd_order", args.fd_order)
    par.set_parval_from_str("EvolvedConformalFactor_cf", args.conformal_factor)
    par.set_parval_from_str("detgbarOverdetghat_equals_one", True)
    par.set_parval_from_str("enable_parallel_codegen", True)
    par.register_CodeParameter(
        "REAL",
        "nrpy.equations.general_relativity.BSSN_gauge_RHSs",
        "eta",
        1.0,
        commondata=True,
    )
    state_h.register_canonical_gridfunctions(enable_fCCZ4=False)

    for fd_order in FD_ORDERS:
        Ricci_eval.register_CFunction_Ricci_eval(
            SOLVER_STEM, fd_order=fd_order, CoordSystem=COORD_SYSTEM
        )
        rhs_eval.register_CFunction_rhs_eval(
            SOLVER_STEM,
            fd_order=fd_order,
            enable_KreissOliger_dissipation=args.ko,
            CoordSystem=COORD_SYSTEM,
            LapseEvolutionOption=LAPSE_EVOLUTION_OPTION,
            ShiftEvolutionOption=SHIFT_EVOLUTION_OPTION,
            enable_SSL=True,
            enable_CAHD=True,
        )
        BSSN_constraints.register_CFunction_BSSN_constraints(
            SOLVER_STEM, fd_order=fd_order, CoordSystem=COORD_SYSTEM
        )
        initial_data_lambdaU.register_CFunction_initial_data_lambdaU(
            SOLVER_STEM, fd_order=fd_order, CoordSystem=COORD_SYSTEM
        )
        ADM_to_BSSN.register_CFunction_ADM_to_BSSN(
            SOLVER_STEM, fd_order=fd_order, CoordSystem=COORD_SYSTEM
        )
    BSSN_to_ADM.register_CFunction_BSSN_to_ADM(SOLVER_STEM, CoordSystem=COORD_SYSTEM)
    enforce_detgbar_equals_detghat_trAzero.register_CFunction_enforce_detgbar_equals_detghat_trAzero(
        SOLVER_STEM, CoordSystem=COORD_SYSTEM
    )
    floor_the_lapse_and_conformal_factor.register_CFunction_floor_the_lapse_and_conformal_factor(
        SOLVER_STEM
    )
    twopunctures.register_CFunction_twopunctures(SOLVER_STEM)
    ADM_Initial_Data_Reader__BSSN_Converter.register_BHaH_defines_h(
        ID_persist_struct.ID_persist_str()
    )
    physical_boundary.register_CFunction_physical_boundary(SOLVER_STEM)
    physical_boundary_ghosts.register_CFunction_physical_boundary_ghosts(SOLVER_STEM)
    diagnostics.register_CFunction_diagnostics(SOLVER_STEM)
    apparent_horizon.register_CFunction_apparent_horizon(SOLVER_STEM)
    psi4_eval.register_CFunction_psi4_eval(SOLVER_STEM)
    gravitational_waves.register_CFunction_gravitational_waves(SOLVER_STEM)
    adm_quantities_surface_data.register_CFunction_adm_quantities_surface_data(
        SOLVER_STEM
    )
    adm_quantities.register_CFunction_adm_quantities(SOLVER_STEM)
    pcg.do_parallel_codegen()

    state_h.validate_registered_state(enable_fCCZ4=False)
    CodeParameters.register_CFunctions_parameters(SOLVER_STEM, SOLVER_NAMESPACE)

    par.set_parval_from_str("parallelization", "openmp")
    with tempfile.TemporaryDirectory() as temporary_directory:
        BHaH_defines_h.output_BHaH_defines_h(
            project_dir=temporary_directory,
            enable_rfm_precompute=False,
        )
        bhah_defines = (Path(temporary_directory) / "BHaH_defines.h").read_text(
            encoding="utf-8"
        )
    bhah_defines = f"""#ifndef NRPY_PACKAGED_BHAH_DEFINES_H
#define NRPY_PACKAGED_BHAH_DEFINES_H
#include "TwoPunctures.h"
{bhah_defines}
#undef NUM_EVOL_GFS
#undef NUM_AUXEVOL_GFS
#undef NUM_AUX_GFS
#undef NUM_SCRATCH_GFS
#undef REAL
#undef DOUBLE
typedef double REAL;
typedef double DOUBLE;
#endif  // NRPY_PACKAGED_BHAH_DEFINES_H
"""
    par.set_parval_from_str("parallelization", "none")
    bhah_prototypes = [
        "#ifndef BHAH_FUNCTION_PROTOTYPES_H",
        "#define BHAH_FUNCTION_PROTOTYPES_H",
        "",
        '#include "BHaH_defines.h"',
        "",
    ]
    bhah_prototypes.extend(
        cfunction.function_prototype
        for name, cfunction in sorted(cfc.CFunction_dict.items())
        if cfunction.subdirectory == "TwoPunctures"
        or name
        in (
            "initialize_ID_persist_struct",
            "NRPyPN_quasicircular_momenta",
        )
    )
    bhah_prototypes.extend(("", "#endif  // BHAH_FUNCTION_PROTOTYPES_H", ""))

    ko_fd_order = args.fd_order - 2
    required_padding = args.fd_order // 2
    module_root = f"{SOLVER_NAME}/"
    application_sources = (
        f"src/{SOLVER_STEM}_main.cpp",
        f"src/{SOLVER_STEM}Ctx.cpp",
        "src/checkpoint.cpp",
    )
    artifacts: Dict[str, str] = {
        module_root + "twopunctures/include/BHaH_defines.h": bhah_defines,
        module_root
        + "twopunctures/include/BHaH_function_prototypes.h": "\n".join(bhah_prototypes),
        module_root
        + "twopunctures/include/TP_utilities.h": (
            Path(TwoPunctures_lib.__file__).parent / "TP_utilities.h"
        ).read_text(encoding="utf-8"),
        module_root
        + "twopunctures/include/TwoPunctures.h": (
            Path(TwoPunctures_lib.__file__).parent / "TwoPunctures.h"
        ).read_text(encoding="utf-8"),
        module_root
        + f"generated/include/{SOLVER_STEM}_types.h": types_h.output_types_h(
            SOLVER_STEM, SOLVER_NAMESPACE
        ),
        module_root
        + f"generated/include/{SOLVER_STEM}_constants.h": constants_h.output_constants_h(
            SOLVER_STEM,
            SOLVER_NAMESPACE,
            args.fd_order,
            ko_fd_order,
            args.fd_order,
            required_padding,
            args.ko,
            FD_ORDERS,
            {order: order // 2 for order in FD_ORDERS},
        ),
        module_root
        + f"generated/include/{SOLVER_STEM}_state.h": state_h.output_state_h(
            SOLVER_STEM, SOLVER_NAMESPACE
        ),
        module_root
        + f"generated/include/{SOLVER_STEM}_parameters.h": CodeParameters.output_parameters_h(
            SOLVER_STEM, SOLVER_NAMESPACE
        ),
        module_root
        + f"generated/include/{SOLVER_STEM}_defines.h": Dendro_defines_h.output_Dendro_defines_h(
            SOLVER_STEM
        ),
        module_root
        + f"include/{SOLVER_STEM}Ctx.h": solver_context.output_solver_context_h(
            SOLVER_STEM,
            SOLVER_NAMESPACE,
            args.fd_order,
            enable_fCCZ4=False,
        ),
        module_root
        + f"src/{SOLVER_STEM}Ctx.cpp": solver_context.output_solver_context_cpp(
            SOLVER_STEM,
            SOLVER_NAMESPACE,
            args.fd_order,
            enable_fCCZ4=False,
            enable_SSL=True,
            enable_CAHD=True,
        ),
        module_root
        + f"src/{SOLVER_STEM}_main.cpp": main_cpp.output_main_cpp(
            SOLVER_STEM, SOLVER_NAMESPACE, EXECUTABLE_NAME, PROFILE_NAME
        ),
        module_root
        + "src/checkpoint.cpp": checkpoint.output_checkpoint_cpp(
            SOLVER_STEM,
            SOLVER_NAMESPACE,
            "BSSN" if args.conformal_factor == "W" else "BSSN_chi",
            enable_fCCZ4=False,
        ),
        module_root
        + f"pars/{SOLVER_STEM}.toml": param_toml.generate_default_parfile(
            SOLVER_STEM,
            PROFILE_NAME,
            args.fd_order,
            ko_fd_order,
            args.fd_order,
            required_padding,
            args.ko,
        ),
    }
    artifacts.update(
        CMakeLists.output_CFunctions_function_prototypes_and_construct_CMakeLists(
            SOLVER_NAME,
            SOLVER_STEM,
            SOLVER_PREFIX,
            EXECUTABLE_NAME,
            application_sources,
        )
    )
    for relative_path, contents in sorted(artifacts.items()):
        target = project_dir / relative_path
        target.parent.mkdir(parents=True, exist_ok=True)
        with ConditionalFileUpdater(
            target,
            encoding="utf-8",
            do_format=target.suffix in (".h", ".hpp", ".c", ".cpp", ".cu"),
        ) as output_file:
            output_file.write(contents)

    solver_dir = project_dir / SOLVER_NAME
    print(f"Finished generating {SOLVER_NAME} in {solver_dir}.")
    print("Generated centered FD/KO profiles: 4/2, 6/4, 8/6.")
    print(f"Selected runtime profile: FD{args.fd_order}/KO{ko_fd_order}.")
    print(
        f"{args.conformal_factor} evolution, SSL, Dendro CAHD, "
        "TwoPunctures alpha=W, and eta=1 enabled."
    )


if __name__ == "__main__":
    main()

"""
Orchestrator for the CPU numerical-spacetime photon geodesic integrator.

This module constructs the C project that evolves batched photon trajectories
through a numerical spacetime sourced from the combined raytracing ``.bin``
generated from ``two_blackholes_collide.py --raytracing-time ...``.

The generated code keeps the existing Structure-of-Arrays photon pipeline and
adaptive RKF45 stepping, while metric and selected geometry evaluations use the
numerical interpolation helpers fed by the combined numerical-spacetime
dataset. The numerical spacetime is frozen to its final stored slice after
`t_numerical_end`.

This example currently generates a CPU/OpenMP project.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import os
import shutil
import sys

import sympy as sp

import nrpy.c_function as cfc
import nrpy.params as par
import nrpy.reference_metric as refmetric
from nrpy.equations.general_relativity.geodesics import geodesics as geo
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation import (
    azimuthal_symmetry_spatial_lagrange_interpolation,
    numerical_interpolation,
    temporal_lagrange_interpolation,
    time_window_manager_numerical,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_numerical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_non_terminal_plane_intersection,
    handle_terminal_plane_intersection,
    main_batch,
    normalization_constraint_photon_normalized,
    photon_momentum_to_normalized_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
)


def _require(condition: bool, message: str) -> None:
    """
    Raise ``ValueError`` if an argument validation condition fails.

    :param condition: Condition that must hold.
    :param message: Description of the violated constraint.
    :raises ValueError: If ``condition`` is false.
    """
    if not condition:
        raise ValueError(message)


SUPPORTED_NUMERICAL_COORDINATE_SYSTEMS = ("SinhCylindricalv2n2",)
SUPPORTED_INTERPOLATION_METHODS = ("g4DD", "g4DD_d0", "GammaUDD")


if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description=(
            "Generate a tiled numerical-spacetime photon raytracing project. "
            "The project emits v6 light-blueprint records with normalized "
            "image-sample coordinates."
        ),
        epilog="""To generate the desired numerical spacetime data, run:
python3 two_blackholes_collide.py --raytracing-time T_FINAL DIAGNOSTICS_OUTPUT_EVERY --raytracing-coord-system CoordSystem --raytracing-Nxx NXX0 NXX1 NXX2 --raytracing-domain ...

Then rerun this photon script using the .bin filename printed to the terminal by two_blackholes_collide.py, e.g.:
python3 photon_batch_geodesic_integrator_numerical.py --bin-name two_blackholes_collide_7p5_7p5_0p25_z1_0p5_z2_neg0p5_M1_0p5_M2_0p5_SinhCylindricalv2n2_sinhwrho_0p25_sinhwz_0p4_rho-slope_0p0625_z-slope_0p0625_72_2_12.bin --coord-system-numerical SinhCylindricalv2n2 --domain 7.5 0.25 0.4 0.0625 0.0625 --t-numerical-end 100.0 --dt-spacetime-data 0.5""",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--bin-name",
        type=str,
        required=True,
        help=(
            "Basename of the combined numerical-spacetime .bin file under "
            "project/raytracing_data."
        ),
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="Parent directory for the generated C project.",
    )
    parser.add_argument(
        "--coord-system-numerical",
        type=str,
        required=True,
        choices=SUPPORTED_NUMERICAL_COORDINATE_SYSTEMS,
        help="Coordinate system used by the numerical-spacetime dataset.",
    )
    parser.add_argument(
        "--domain",
        nargs="+",
        type=float,
        required=True,
        help="""SinhCylindricalv2n2 domain: GRID_PHYSICAL_SIZE SINHWRHO SINHWZ RHO_SLOPE Z_SLOPE.""",
    )
    parser.add_argument(
        "--t-numerical-end",
        type=float,
        required=True,
        help="""Coordinate time of the final stored numerical slice. After this time the numerical spacetime is frozen to that final slice.""",
    )
    parser.add_argument(
        "--dt-spacetime-data",
        type=float,
        required=True,
        help="""Nominal coordinate-time spacing; actual time gaps are approximately this size, while stored slice-table times remain authoritative.""",
    )
    parser.add_argument(
        "--t-start",
        type=float,
        required=True,
        help="Initial observer-event coordinate time for every ray in the batch.",
    )
    parser.add_argument(
        "--eom",
        type=str,
        choices=("geodesic", "normalized"),
        required=True,
        help="""Photon equations of motion: direct affine-parameter geodesics or normalized coordinate-time evolution.""",
    )
    parser.add_argument(
        "--interpolation-method",
        choices=SUPPORTED_INTERPOLATION_METHODS,
        required=True,
        help="""Numerical spacetime payload and interpolation method.""",
    )
    parser.add_argument(
        "--observer-position",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Observer coordinate position.",
    )
    parser.add_argument(
        "--observer-look-forward",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Observer look-forward direction seed.",
    )
    parser.add_argument(
        "--observer-up",
        nargs=3,
        type=float,
        required=True,
        metavar=("X", "Y", "Z"),
        help="Observer up/roll direction seed.",
    )
    parser.add_argument(
        "--observer-fov",
        nargs=2,
        type=float,
        required=True,
        metavar=("ALPHA_W", "ALPHA_H"),
        help="Observer width and height fields of view, in radians.",
    )
    parser.add_argument(
        "--scan-density",
        type=int,
        required=True,
        metavar="SAMPLES_PER_TILE_WIDTH",
        help="Width-side ray samples per tile; height-side count is derived.",
    )
    parser.add_argument(
        "--tile-counts",
        nargs=2,
        type=int,
        default=[1, 1],
        metavar=("TILES_WIDTH", "TILES_HEIGHT"),
        help="Number of angular image tiles; default: 1 1.",
    )
    parser.add_argument(
        "--escape-radius",
        type=float,
        required=True,
        metavar="R",
        help="Coordinate-radius termination threshold.",
    )
    parser.add_argument(
        "--terminal-plane-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="Terminal-plane center; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-normal",
        nargs=3,
        type=float,
        metavar=("NX", "NY", "NZ"),
        help="Terminal-plane normal; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-up",
        nargs=3,
        type=float,
        metavar=("UX", "UY", "UZ"),
        help="Terminal-plane up direction; requires all terminal-plane options.",
    )
    parser.add_argument(
        "--terminal-plane-radius",
        nargs=2,
        type=float,
        metavar=("MIN_RADIUS", "MAX_RADIUS"),
        help="Terminal-plane accepted coordinate-radius range.",
    )
    parser.add_argument(
        "--non-terminal-plane-center",
        nargs=3,
        type=float,
        metavar=("X", "Y", "Z"),
        help="Nonterminal-plane center; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--non-terminal-plane-normal",
        nargs=3,
        type=float,
        metavar=("NX", "NY", "NZ"),
        help="Nonterminal-plane normal; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--non-terminal-plane-up",
        nargs=3,
        type=float,
        metavar=("UX", "UY", "UZ"),
        help="Nonterminal-plane up direction; requires all nonterminal-plane options.",
    )
    parser.add_argument(
        "--initial-step",
        type=float,
        metavar="H",
        help="Initial RKF45 step magnitude; sign follows selected EOM.",
    )
    parser.add_argument(
        "--rkf45-tolerances",
        nargs=2,
        type=float,
        metavar=("ABSOLUTE", "RELATIVE"),
        help="RKF45 absolute and relative error tolerances.",
    )
    parser.add_argument(
        "--rkf45-step-range",
        nargs=2,
        type=float,
        metavar=("H_MIN", "H_MAX"),
        help="RKF45 minimum and maximum step magnitudes.",
    )
    parser.add_argument(
        "--rkf45-max-retries",
        type=int,
        metavar="N",
        help="Maximum RKF45 rejected-step retries; default: current value.",
    )
    parser.add_argument(
        "--evolution-measure-max",
        type=float,
        metavar="VALUE",
        help=(
            "Maximum evolution measure: abs(p^0) in direct mode, or "
            "abs(ln(abs(alpha*p^0))) in normalized mode."
        ),
    )
    parser.add_argument(
        "--rkf45-log-energy-tolerance",
        type=float,
        metavar="VALUE",
        help="Normalized-EOM RKF45 log-energy tolerance; invalid in direct mode.",
    )
    parser.add_argument(
        "--interpolation-half-widths",
        nargs=2,
        type=int,
        metavar=("SPATIAL", "TEMPORAL"),
        help="Numerical spatial and temporal interpolation half widths.",
    )
    parser.add_argument(
        "--time-window-slot-parameters",
        nargs=2,
        type=float,
        metavar=("T_MIN", "DELTA_T"),
        help="Numerical time-window slot-manager lower time and slot width.",
    )
    parser.add_argument(
        "--time-window-max-delta-t",
        type=float,
        metavar="DELTA_T",
        help="Maximum coordinate-time change per normalized RKF45 step.",
    )
    parser.add_argument(
        "--time-slice-stride",
        type=int,
        metavar="N",
        help="Retained numerical-spacetime time-slice stride.",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()
    normalized_eom = args.eom == "normalized"
    interpolation_method = args.interpolation_method
    rhs_uses_metric_derivatives = interpolation_method != "GammaUDD"

    _require(
        os.path.basename(args.bin_name) == args.bin_name,
        "--bin-name must be a filename only, not a path.",
    )
    _require(
        args.bin_name.endswith(".bin"),
        "--bin-name must name a numerical spacetime .bin file.",
    )
    _require(
        args.dt_spacetime_data > 0.0,
        "--dt-spacetime-data must be positive.",
    )
    _require(
        len(args.domain) > 0,
        "--domain requires at least one value.",
    )
    _require(
        args.t_start >= 0.0,
        "--t-start must be nonnegative.",
    )
    _require(
        all(value > 0.0 for value in args.observer_fov),
        "--observer-fov values must be positive.",
    )
    _require(args.scan_density >= 1, "--scan-density must be positive.")
    _require(
        all(value >= 1 for value in args.tile_counts),
        "--tile-counts values must be positive.",
    )
    _require(args.escape_radius > 0.0, "--escape-radius must be positive.")
    terminal_plane_values = (
        args.terminal_plane_center,
        args.terminal_plane_normal,
        args.terminal_plane_up,
        args.terminal_plane_radius,
    )
    _require(
        not any(value is not None for value in terminal_plane_values)
        or all(value is not None for value in terminal_plane_values),
        "Terminal-plane options must be supplied as one complete group.",
    )
    non_terminal_plane_values = (
        args.non_terminal_plane_center,
        args.non_terminal_plane_normal,
        args.non_terminal_plane_up,
    )
    _require(
        not any(value is not None for value in non_terminal_plane_values)
        or all(value is not None for value in non_terminal_plane_values),
        "Nonterminal-plane options must be supplied as one complete group.",
    )
    temporal_margin = 4.0 * args.dt_spacetime_data
    if args.t_start + temporal_margin > args.t_numerical_end:
        print(
            "WARNING: t-start is most likely too close to the numerical "
            "spacetime t_final for centered temporal interpolation using only "
            "the provided spacetime data. The last temporal slice of the .bin "
            "file will fill missing data. If the numerical evolution is not "
            "approximately static by its last slice, integration and photon "
            "trajectories may be inaccurate."
        )
        print(
            "         Set --t-start <= "
            f"{args.t_numerical_end - temporal_margin:.3f} "
            "(generated .par: t_start)."
        )
        if interpolation_method == "GammaUDD":
            print(
                "WARNING: GammaUDD users who do not lower --t-start should "
                "generate the numerical spacetime with "
                "'--raytracing-static-christoffels'."
            )
    domain = list(args.domain)
    grid_physical_size = None
    sinhw_numerical_rho = None
    sinhw_numerical_z = None
    rho_slope = None
    z_slope = None
    if args.coord_system_numerical == "SinhCylindricalv2n2":
        _require(
            len(domain) == 5,
            """--coord-system-numerical SinhCylindricalv2n2 requires --domain GRID_PHYSICAL_SIZE SINHWRHO SINHWZ RHO_SLOPE Z_SLOPE.""",
        )
        (
            grid_physical_size,
            sinhw_numerical_rho,
            sinhw_numerical_z,
            rho_slope,
            z_slope,
        ) = domain
        _require(sinhw_numerical_rho > 0.0, "SINHWRHO must be positive.")
        _require(sinhw_numerical_z > 0.0, "SINHWZ must be positive.")
        _require(rho_slope > 0.0, "RHO_SLOPE must be positive.")
        _require(z_slope > 0.0, "Z_SLOPE must be positive.")
    else:
        raise ValueError(
            f"Unsupported numerical coordinate system: {args.coord_system_numerical}."
        )
    assert grid_physical_size is not None
    _require(grid_physical_size > 0.0, "GRID_PHYSICAL_SIZE must be positive.")

    # Step 2: Define strict project constants and simulation targets
    project_name = f"photon_batch_geodesic_integrator_numerical_{interpolation_method}"
    exec_name = project_name
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))

    # Step 2.a: Define the target spacetime and particle type.
    SPACETIME = "Numerical"
    integrator_mode = "Numerical"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"
    coord_system_numerical = args.coord_system_numerical
    enable_simd = False

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Step 3.a: Select the BHaH OpenMP backend for this example.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")
    par.set_parval_from_str(
        "CoordSystem_to_register_CodeParameters", coord_system_numerical
    )

    # Step 4: Build the symbolic photon data used by numerical evolution.
    print(f" -> Assembling symbolic data for {GEO_KEY}...")
    # Step 4.a: Bypass __init__ because these symbolic helper methods build
    # placeholder-based expressions and do not require spacetime-specific state.
    generic_geodesic_equations = geo.GeodesicEquations.__new__(geo.GeodesicEquations)
    coordinate_symbols = list(sp.symbols("t x y z", real=True))
    if rhs_uses_metric_derivatives:
        geodesic_rhs = (
            generic_geodesic_equations.geodesic_eom_rhs_photon_normalized()
            if normalized_eom
            else generic_geodesic_equations.geodesic_eom_rhs_photon()
        )
    else:
        geodesic_rhs = (
            generic_geodesic_equations.geodesic_eom_rhs_photon_normalized_christoffel()
            if normalized_eom
            else generic_geodesic_equations.geodesic_eom_rhs_photon_christoffel()
        )
    normalization_constraint_expr = (
        generic_geodesic_equations.normalization_constraint_photon_normalized()
        if normalized_eom
        else generic_geodesic_equations.normalization_constraint()
    )

    # Step 5: Register C functions in split-pipeline order.
    print(" -> Registering C functions and local CodeParameters...")

    # Step 5.a: Register batch state definitions and initialization kernels.
    set_initial_conditions_kernel.register_photon_batch_structs()
    set_initial_conditions_kernel.set_initial_conditions_kernel(
        normalized_eom=normalized_eom
    )

    if normalized_eom:
        u_expr, PiD_exprs = (
            geo.GeodesicEquations.photon_momentum_to_normalized_quantities()
        )
        photon_momentum_to_normalized_kernel.photon_momentum_to_normalized_kernel(
            u_expr, PiD_exprs
        )

    # Step 5.b: Register normalization diagnostics.
    if normalized_eom:
        normalization_constraint_photon_normalized.normalization_constraint_photon_normalized(
            normalization_constraint_expr
        )
    else:
        normalization_constraint.normalization_constraint(
            normalization_constraint_expr, PARTICLE
        )

    # Step 5.c: Register numerical interpolation helpers.
    register_azimuthal_interp = (
        azimuthal_symmetry_spatial_lagrange_interpolation.register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation
    )
    register_azimuthal_interp(
        coord_system_numerical,
        interpolation_method=interpolation_method,
        enable_simd=enable_simd,
        project_dir=project_dir,
    )
    temporal_lagrange_interpolation.register_CFunction_temporal_lagrange_interpolation(
        interpolation_method=interpolation_method,
        enable_simd=enable_simd,
        project_dir=project_dir,
    )
    time_window_manager_numerical.time_window_manager_numerical(
        interpolation_method,
        include_batch_startup_validation=True,
    )
    numerical_interpolation.register_CFunction_numerical_interpolation(
        coord_system_numerical,
        interpolation_method=interpolation_method,
        enable_simd=enable_simd,
        project_dir=project_dir,
        normalized_eom=normalized_eom,
    )

    # Step 5.d: Register RKF45 evolution kernels.
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(
        geodesic_rhs,
        coordinate_symbols,
        rhs_uses_metric_derivatives=rhs_uses_metric_derivatives,
        normalized_eom=normalized_eom,
    )
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True,
        normalized_eom=normalized_eom,
    )

    # Step 5.e: Register event-detection and boundary-intersection kernels.
    find_event_time_and_state.find_event_time_and_state()
    handle_terminal_plane_intersection.handle_terminal_plane_intersection()
    handle_non_terminal_plane_intersection.handle_non_terminal_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel(
        normalized_eom=normalized_eom
    )
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal(
        normalized_eom=normalized_eom
    )

    # Step 5.f: Register project-level orchestration helpers.
    # The numerical interpolation and RKF45 finalization registrations above
    # already pull in the shared slot/time-window helpers on demand. Calling
    # them again here would append duplicate definitions into BHaH_defines.h.
    batch_integrator_numerical.batch_integrator_numerical(
        SPACETIME,
        coord_system_numerical,
        interpolation_method=interpolation_method,
        normalized_eom=normalized_eom,
    )
    main_batch.main(SPACETIME, integrator_mode, normalized_eom=normalized_eom)

    # Step 5.g: Remove helper registrations emitted only through other kernels.
    # The event manager emits its local event helpers through prefunc, so keeping
    # standalone registrations would add unused source files and prototypes.
    internal_funcs_to_remove = [
        "find_event_time_and_state",
        "handle_terminal_plane_intersection",
        "handle_non_terminal_plane_intersection",
    ]
    for internal_func in internal_funcs_to_remove:
        cfc.CFunction_dict.pop(internal_func, None)

    # Step 6: Override CodeParameter defaults before parfile generation.
    print(" -> Overriding desired photon CodeParameters before .par generation...")

    # Step 6.a: Set the numerical-spacetime data path.
    # This path intentionally stays under the repo-root raytracing cache.
    # ``--outdir`` controls only the generated photon project location.
    numerical_spacetime_bin_path = os.path.abspath(
        os.path.join(
            "project",
            "raytracing_data",
            args.bin_name,
        )
    )

    par.adjust_CodeParam_default(
        "numerical_spacetime_bin_path", numerical_spacetime_bin_path
    )

    # Step 6.b: Match the dataset reference-metric defaults used by the
    # combined raytracing file before metadata overrides Nxx/dxx/xxmin/xxmax.
    par.adjust_CodeParam_default("grid_physical_size", grid_physical_size)
    rfm = refmetric.reference_metric[coord_system_numerical]
    for param_name, grid_size_mapping in rfm.grid_physical_size_dict.items():
        if grid_size_mapping == "grid_physical_size":
            par.adjust_CodeParam_default(param_name, grid_physical_size)
        elif grid_size_mapping == "-grid_physical_size":
            par.adjust_CodeParam_default(param_name, -grid_physical_size)
        else:
            raise ValueError(
                f"""Unsupported grid_physical_size mapping '{grid_size_mapping}' for {coord_system_numerical}:{param_name}"""
            )
    assert sinhw_numerical_rho is not None
    assert sinhw_numerical_z is not None
    assert rho_slope is not None
    assert z_slope is not None
    par.adjust_CodeParam_default("SINHWRHO", sinhw_numerical_rho)
    par.adjust_CodeParam_default("SINHWZ", sinhw_numerical_z)
    par.adjust_CodeParam_default("rho_slope", rho_slope)
    par.adjust_CodeParam_default("z_slope", z_slope)

    # Step 6.c: Set interpolation half-width defaults.
    par.adjust_CodeParam_default("numerical_spacetime_spatial_interp_half_width", 3)
    par.adjust_CodeParam_default("numerical_spacetime_temporal_interp_half_width", 3)

    # Step 6.d: Set initial-condition defaults.
    par.adjust_CodeParam_default("t_start", args.t_start)
    par.adjust_CodeParam_default("scan_density", args.scan_density)

    # Step 6.e: Set batch-integrator and numerical-limit defaults.
    par.adjust_CodeParam_default(
        "evolution_measure_max", 3.0 if normalized_eom else 1000.0
    )
    par.adjust_CodeParam_default("perform_normalization_check", True)
    par.adjust_CodeParam_default("r_escape", args.escape_radius)

    # Step 6.f: Set numerical-spacetime time-range defaults.
    par.adjust_CodeParam_default("t_numerical_end", args.t_numerical_end)
    par.adjust_CodeParam_default("dt_numerical_spacetime_data", args.dt_spacetime_data)
    # Keep sparse retained-slice debugging disabled for normal runs.
    par.adjust_CodeParam_default("numerical_spacetime_time_slice_stride", 1)

    # Step 6.g: Cap coordinate-time growth per accepted RKF45 step.
    par.adjust_CodeParam_default("rkf45_max_delta_t", 0.5)

    # Step 6.h: Set time-window manager defaults.
    par.adjust_CodeParam_default("slot_manager_delta_t", 5.0)
    par.adjust_CodeParam_default("slot_manager_t_min", -150.0)

    # Step 6.i: Set observer and independent plane geometry from the CLI.
    for name, value in zip(
        ["observer_x", "observer_y", "observer_z"], args.observer_position
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        [
            "observer_look_forward_x",
            "observer_look_forward_y",
            "observer_look_forward_z",
        ],
        args.observer_look_forward,
    ):
        par.adjust_CodeParam_default(name, value)
    for name, value in zip(
        ["observer_up_x", "observer_up_y", "observer_up_z"], args.observer_up
    ):
        par.adjust_CodeParam_default(name, value)
    par.adjust_CodeParam_default("alpha_w", args.observer_fov[0])
    par.adjust_CodeParam_default("alpha_h", args.observer_fov[1])

    terminal_defaults = {
        "terminal_plane_center_x": -1.0e4,
        "terminal_plane_center_y": 0.0,
        "terminal_plane_center_z": 0.0,
        "terminal_plane_normal_x": 1.0,
        "terminal_plane_normal_y": 0.0,
        "terminal_plane_normal_z": 0.0,
        "terminal_plane_up_x": 0.0,
        "terminal_plane_up_y": 0.0,
        "terminal_plane_up_z": 1.0,
        "terminal_plane_min_coord_radius": 0.0,
        "terminal_plane_max_coord_radius": 1.0,
        "terminal_plane_enabled": False,
    }
    if args.terminal_plane_center is not None:
        terminal_defaults.update(
            {
                "terminal_plane_center_x": args.terminal_plane_center[0],
                "terminal_plane_center_y": args.terminal_plane_center[1],
                "terminal_plane_center_z": args.terminal_plane_center[2],
                "terminal_plane_normal_x": args.terminal_plane_normal[0],
                "terminal_plane_normal_y": args.terminal_plane_normal[1],
                "terminal_plane_normal_z": args.terminal_plane_normal[2],
                "terminal_plane_up_x": args.terminal_plane_up[0],
                "terminal_plane_up_y": args.terminal_plane_up[1],
                "terminal_plane_up_z": args.terminal_plane_up[2],
                "terminal_plane_min_coord_radius": args.terminal_plane_radius[0],
                "terminal_plane_max_coord_radius": args.terminal_plane_radius[1],
                "terminal_plane_enabled": True,
            }
        )
    non_terminal_defaults = {
        "non_terminal_plane_center_x": -1.0e4,
        "non_terminal_plane_center_y": 0.0,
        "non_terminal_plane_center_z": 0.0,
        "non_terminal_plane_normal_x": 1.0,
        "non_terminal_plane_normal_y": 0.0,
        "non_terminal_plane_normal_z": 0.0,
        "non_terminal_plane_up_x": 0.0,
        "non_terminal_plane_up_y": 0.0,
        "non_terminal_plane_up_z": 1.0,
        "non_terminal_plane_enabled": False,
    }
    if args.non_terminal_plane_center is not None:
        non_terminal_defaults.update(
            {
                "non_terminal_plane_center_x": args.non_terminal_plane_center[0],
                "non_terminal_plane_center_y": args.non_terminal_plane_center[1],
                "non_terminal_plane_center_z": args.non_terminal_plane_center[2],
                "non_terminal_plane_normal_x": args.non_terminal_plane_normal[0],
                "non_terminal_plane_normal_y": args.non_terminal_plane_normal[1],
                "non_terminal_plane_normal_z": args.non_terminal_plane_normal[2],
                "non_terminal_plane_up_x": args.non_terminal_plane_up[0],
                "non_terminal_plane_up_y": args.non_terminal_plane_up[1],
                "non_terminal_plane_up_z": args.non_terminal_plane_up[2],
                "non_terminal_plane_enabled": True,
            }
        )
    for name, value in {**terminal_defaults, **non_terminal_defaults}.items():
        par.adjust_CodeParam_default(name, value)

    # Step 6.k: Set nonterminal-plane tiling defaults.
    par.adjust_CodeParam_default("tiles_width", args.tile_counts[0])
    par.adjust_CodeParam_default("tiles_height", args.tile_counts[1])

    # Step 6.l: Set RKF45 controller defaults.
    par.adjust_CodeParam_default("initial_h", -0.05 if normalized_eom else 0.05)
    par.adjust_CodeParam_default("rkf45_absolute_error_tolerance", 1.0e-8)
    par.adjust_CodeParam_default("rkf45_error_tolerance", 1.0e-8)
    par.adjust_CodeParam_default("rkf45_h_max", 10.0)
    par.adjust_CodeParam_default("rkf45_h_min", 1.0e-4)
    if normalized_eom:
        par.adjust_CodeParam_default("rkf45_log_energy_tolerance", 1.0e0)

    if args.initial_step is not None:
        _require(args.initial_step > 0.0, "--initial-step must be positive.")
        par.adjust_CodeParam_default(
            "initial_h", -args.initial_step if normalized_eom else args.initial_step
        )
    if args.rkf45_tolerances is not None:
        absolute_tolerance, relative_tolerance = args.rkf45_tolerances
        _require(
            absolute_tolerance > 0.0 and relative_tolerance > 0.0,
            "RKF45 tolerances must be positive.",
        )
        par.adjust_CodeParam_default(
            "rkf45_absolute_error_tolerance", absolute_tolerance
        )
        par.adjust_CodeParam_default("rkf45_error_tolerance", relative_tolerance)
    if args.rkf45_step_range is not None:
        h_min, h_max = args.rkf45_step_range
        _require(0.0 < h_min <= h_max, "Invalid RKF45 step range.")
        par.adjust_CodeParam_default("rkf45_h_min", h_min)
        par.adjust_CodeParam_default("rkf45_h_max", h_max)
    if args.rkf45_max_retries is not None:
        _require(
            args.rkf45_max_retries >= 0, "--rkf45-max-retries must be nonnegative."
        )
        par.adjust_CodeParam_default("rkf45_max_retries", args.rkf45_max_retries)
    if args.evolution_measure_max is not None:
        _require(
            args.evolution_measure_max > 0.0,
            "--evolution-measure-max must be positive.",
        )
        par.adjust_CodeParam_default(
            "evolution_measure_max", args.evolution_measure_max
        )
    if args.rkf45_log_energy_tolerance is not None:
        _require(
            normalized_eom,
            "--rkf45-log-energy-tolerance requires --eom normalized.",
        )
        _require(
            args.rkf45_log_energy_tolerance > 0.0,
            "--rkf45-log-energy-tolerance must be positive.",
        )
        par.adjust_CodeParam_default(
            "rkf45_log_energy_tolerance", args.rkf45_log_energy_tolerance
        )
    if args.interpolation_half_widths is not None:
        spatial_width, temporal_width = args.interpolation_half_widths
        _require(
            spatial_width >= 1 and temporal_width >= 1,
            "Interpolation half widths must be positive.",
        )
        par.adjust_CodeParam_default(
            "numerical_spacetime_spatial_interp_half_width", spatial_width
        )
        par.adjust_CodeParam_default(
            "numerical_spacetime_temporal_interp_half_width", temporal_width
        )
    if args.time_window_slot_parameters is not None:
        slot_t_min, slot_delta_t = args.time_window_slot_parameters
        _require(slot_delta_t > 0.0, "Time-window slot delta_t must be positive.")
        par.adjust_CodeParam_default("slot_manager_t_min", slot_t_min)
        par.adjust_CodeParam_default("slot_manager_delta_t", slot_delta_t)
    if args.time_window_max_delta_t is not None:
        _require(
            args.time_window_max_delta_t > 0.0,
            "Time-window max delta_t must be positive.",
        )
        par.adjust_CodeParam_default("rkf45_max_delta_t", args.time_window_max_delta_t)
    if args.time_slice_stride is not None:
        _require(args.time_slice_stride >= 1, "Time-slice stride must be positive.")
        par.adjust_CodeParam_default(
            "numerical_spacetime_time_slice_stride", args.time_slice_stride
        )

    print(f" -> Numerical spacetime .bin path: {numerical_spacetime_bin_path}")
    print(f" -> Numerical coordinate system: {coord_system_numerical}")
    print(f" -> Interpolation method: {interpolation_method}")
    print(f" -> Photon equations of motion: {args.eom}")
    print(f" -> Numerical domain: {domain}")
    print(f" -> Numerical SINHWRHO: {sinhw_numerical_rho}")
    print(f" -> Numerical SINHWZ: {sinhw_numerical_z}")
    print(f" -> Numerical rho_slope: {rho_slope}")
    print(f" -> Numerical z_slope: {z_slope}")
    print(f" -> t_start: {args.t_start}")
    print(f" -> t_numerical_end: {args.t_numerical_end}")
    print(f" -> dt_spacetime_data: {args.dt_spacetime_data}")

    # Step 6.m: Generate C code for parameter handling.
    print(" -> Generating parameter handling code...")
    CPs.write_CodeParameters_h_files(project_dir=project_dir, set_commondata_only=True)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name,
        cmdline_inputs=[
            "terminal_plane_min_coord_radius",
            "terminal_plane_max_coord_radius",
            "tiles_width",
            "tiles_height",
            "scan_density",
        ],
    )

    # Step 7: Assemble the generated C project.
    print(" -> Assembling C project on disk...")

    # Step 7.a: Map the shared CUDA-style helpers onto OpenMP.
    cpu_macros = {
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))\n",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))\n",
        "BHAH_MALLOC_DEVICE(a, sz)": (
            "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)\n"
        ),
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)\n",
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))\n",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))\n",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))\n",
        "FusedMulAddCUDA(a, b, c)": (
            "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))\n"
        ),
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)\n",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)\n",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)\n",
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC\n",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline\n",
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": (
            "#define BHAH_WARP_ATOMIC_ADD(ptr, val) "
            '_Pragma("omp atomic") *(ptr) += (val)\n'
        ),
        "GLOBAL_COMMONDATA_EXTERN": (
            "// CPU passes commondata by reference, no global needed.\n"
        ),
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
    }

    # Step 7.b: Match the combined raytracing payload layout, which uses
    # NGHOSTS = 3.
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
        supplemental_defines_dict=cpu_macros,
    )

    compiler = "gcc"
    cflags = [
        "-fopenmp",
        "-O3",
        "-fno-omit-frame-pointer",
        "-DDEBUG",
        "-Wno-stringop-truncation",
    ]
    libs = ["-lm"]
    ext = "c"

    print(" -> Generating Makefile")

    # Step 7.c: Select the Makefile optimization profile.
    opt_option = "fast"

    Makefile.output_CFunctions_function_prototypes_and_construct_Makefile(
        project_dir=project_dir,
        project_name=project_name,
        exec_or_library_name=exec_name,
        compiler_opt_option=opt_option,
        addl_CFLAGS=cflags,
        addl_libraries=libs,
        CC=compiler,
        src_code_file_ext=ext,
    )
    # Step 8: Copy the v6 blueprint schema/reader/renderer/diagnostic helpers
    # and print usage instructions.
    vis_dir = os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "geodesic_visualizations"
    )
    vis_script_src = os.path.join(vis_dir, "visualize_lensed_image.py")
    config_src = os.path.join(vis_dir, "blueprint_config_and_schema.py")
    render_src = os.path.join(vis_dir, "render_lensed_image.py")
    blueprint_io_src = os.path.join(vis_dir, "blueprint_io.py")
    blueprint_analysis_src = os.path.join(vis_dir, "blueprint_analysis.py")

    for script_src in (
        vis_script_src,
        config_src,
        render_src,
        blueprint_io_src,
        blueprint_analysis_src,
    ):
        shutil.copy(script_src, project_dir)

    combiner_src = os.path.join(
        "nrpy",
        "infrastructures",
        "BHaH",
        "diagnostics",
        "combine_raytracing_time_slices.py",
    )
    shutil.copy(combiner_src, project_dir)

    # Step 8.a: Build helper commands from CodeParameter defaults.
    c_r_min = float(
        par.glb_code_params_dict["terminal_plane_min_coord_radius"].defaultvalue
    )
    c_r_max = float(
        par.glb_code_params_dict["terminal_plane_max_coord_radius"].defaultvalue
    )
    c_scan_density = int(par.glb_code_params_dict["scan_density"].defaultvalue)
    # Visualization output width; angular ray-sample density comes from
    # CodeParameters and the renderer maps normalized fractions to output pixels.
    c_pixel_width = 450
    data_prep_command = (
        "python3 combine_raytracing_time_slices.py "
        f'--input-dir "../two_blackholes_collide/raytracing_slices/{interpolation_method}" '
        '--pattern "raytracing_data_t????????.bin" '
        '--run-metadata "../two_blackholes_collide/raytracing_run_metadata.json" '
        f'--output "../raytracing_data/{args.bin_name}" --force'
    )
    parfile_path = os.path.join(project_dir, f"{project_name}.par")

    vis_command_parts = ["python3 visualize_lensed_image.py"]
    if args.terminal_plane_radius is not None:
        vis_command_parts.append(f"--terminal-plane-radius {c_r_min} {c_r_max}")
    vis_command_parts.append(f"--pixel-width {c_pixel_width}")
    vis_command = " ".join(vis_command_parts)

    print(f"Finished! Now go into {project_dir}.")
    print(f"    Parameter file can be found at {parfile_path}\n")
    print("    To prepare the required numerical spacetime data, run:")
    print(
        "    First generate the source evolution with "
        "two_blackholes_collide.py --raytracing-outputs."
    )
    print(f"    {data_prep_command}\n")
    print("    Then build and run the photon executable:")
    print("    make")
    print(f"    ./{exec_name}\n")
    print(
        """    To generate the lensed image after running the C executable, ensure you have the required Python packages:"""
    )
    print("    pip install matplotlib numpy numba Pillow\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print(f"    {vis_command}\n")
    print(
        "    Ray sampling uses "
        f"{c_scan_density} width-side samples per tile; height-side sampling "
        "is derived from the fields of view and tile-grid aspect ratio.\n"
    )
    print("    To use your own celestial sphere image, add:")
    print("    --sphere-image /path/to/celestial_sphere.png")
    print("    To use your own terminal-plane image, add:")
    print("    --terminal-plane-image /path/to/terminal_plane_image.png\n")

    blueprint_command = "python3 blueprint_analysis.py"

    print("    To run the blueprint diagnostic:")
    print(f"    {blueprint_command}\n")

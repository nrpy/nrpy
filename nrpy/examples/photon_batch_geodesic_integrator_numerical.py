r"""
Orchestrator for the CPU numerical-spacetime photon geodesic integrator.

This module constructs the C project that evolves batched photon trajectories
through a numerical spacetime sourced from the combined raytracing ``.bin``
generated from ``two_blackholes_collide.py --raytracing-outputs``.

The generated code keeps the existing Structure-of-Arrays photon pipeline and
adaptive RKF45 stepping, but all metric and Christoffel evaluations now come
from the numerical interpolation helpers instead of analytic spacetime kernels.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

# ##############################################################################
# PART 0: IMPORTS AND PATH SETUP
# ##############################################################################

import argparse
import math
import os
import shutil
import sys
from typing import Dict, Union

import sympy as sp

# NRPy core and helper modules for C code generation
import nrpy.c_function as cfc
import nrpy.params as par

# Physics/Math Generators (Symbolic definitions of geodesics)
from nrpy.equations.general_relativity.geodesics import geodesics as geo
from nrpy.examples.geodesic_helpers.combined_raytracing_bin_helper import (
    ensure_required_combined_bin,
)

# NRPy BlackHoles@Home (BHaH) infrastructure modules for C project management
from nrpy.infrastructures.BHaH import BHaH_defines_h
from nrpy.infrastructures.BHaH import CodeParameters as CPs
from nrpy.infrastructures.BHaH import Makefile_helpers as Makefile
from nrpy.infrastructures.BHaH import cmdline_input_and_parfiles

# C-Code Builder Functions (Registers specific physical and numerical routines)
from nrpy.infrastructures.BHaH.general_relativity.geodesics import (
    normalization_constraint,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.azimuthal_symmetry_spatial_lagrange_interpolation import (
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.numerical_interpolation import (
    register_CFunction_numerical_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.interpolation.temporal_lagrange_interpolation import (
    register_CFunction_temporal_lagrange_interpolation,
)
from nrpy.infrastructures.BHaH.general_relativity.geodesics.photon import (
    batch_integrator_numerical,
    calculate_and_fill_blueprint_data_universal,
    calculate_ode_rhs_kernel,
    event_detection_manager_kernel,
    find_event_time_and_state,
    handle_source_plane_intersection,
    handle_window_plane_intersection,
    main,
    p0_reverse_kernel,
    rkf45_finalize_and_control_kernel,
    rkf45_stage_update,
    set_initial_conditions_kernel,
)

# Establish the script directory for reliable relative path resolution
script_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(script_dir))


def _require(condition: bool, message: str) -> None:
    """
    Raise ``ValueError`` if a derived example constraint is violated.

    :param condition: Boolean condition that must hold.
    :param message: Error message describing the violated constraint.
    :raises ValueError: If ``condition`` is false.
    """
    if not condition:
        raise ValueError(message)


# ##############################################################################
# PART 1: MAIN CONFIGURATION
# ##############################################################################

if __name__ == "__main__":
    # Step 1: Configure command-line arguments for the generation pipeline
    parser = argparse.ArgumentParser(
        description="Generate the Split-Pipeline Photon Geodesic Integrator (Numerical)."
    )
    parser.add_argument(
        "--outdir",
        type=str,
        default="project",
        help="The parent directory where the generated C project will reside.",
    )
    args = parser.parse_args()

    # Step 2: Define strict project constants and simulation targets
    project_name = "photon_batch_geodesic_integrator_numerical"
    exec_name = "photon_batch_geodesic_integrator_numerical"
    project_dir = os.path.abspath(os.path.join(args.outdir, project_name))

    # Define the physical regime: Target spacetime and particle classification
    SPACETIME = "Numerical"
    integrator_mode = "Numerical"
    PARTICLE = "photon"
    GEO_KEY = f"{SPACETIME}_{PARTICLE}"
    dataset_coord_system = "Spherical"
    enable_simd = False

    # Step 3: Initialize the project directory and select the infrastructure backend
    print(f"Initializing project: {project_name}")
    shutil.rmtree(project_dir, ignore_errors=True)
    os.makedirs(project_dir, exist_ok=True)

    # Instruct NRPy to emit the CPU/OpenMP BHaH pipeline used by the numerical integrator.
    par.set_parval_from_str("Infrastructure", "BHaH")
    par.set_parval_from_str("parallelization", "openmp")

    # Step 4: Build the generic symbolic photon equations consumed by the
    # runtime numerical metric and Christoffel interpolation pipeline.
    print(f" -> Assembling symbolic data for {GEO_KEY}...")
    generic_geodesic_equations = geo.GeodesicEquations.__new__(geo.GeodesicEquations)
    coordinate_symbols = list(sp.symbols("t x y z", real=True))
    geodesic_rhs = generic_geodesic_equations.geodesic_eom_rhs_photon()
    p0_photon = generic_geodesic_equations.hamiltonian_constraint_photon()
    normalization_constraint_expr = (
        generic_geodesic_equations.normalization_constraint()
    )

    # ##########################################################################
    # Step 5: Execute Modules and Register C Functions (Split-Pipeline Order)
    # ##########################################################################
    print(" -> Registering C functions and local CodeParameters...")

    # --- Initialization Phase Kernels ---
    set_initial_conditions_kernel.set_initial_conditions_kernel(SPACETIME)

    p0_reverse_kernel.p0_reverse_kernel(p0_photon)

    # --- Fundamental Constraint Evaluation ---
    normalization_constraint.normalization_constraint(
        normalization_constraint_expr, PARTICLE
    )

    # --- Numerical Interpolation Pipeline ---
    register_CFunction_azimuthal_symmetry_spatial_lagrange_interpolation(
        dataset_coord_system, enable_simd=enable_simd, project_dir=project_dir
    )
    register_CFunction_temporal_lagrange_interpolation(
        enable_simd=enable_simd, project_dir=project_dir
    )
    register_CFunction_numerical_interpolation(
        dataset_coord_system, enable_simd=enable_simd, project_dir=project_dir
    )

    # --- Core RKF45 and Geodesic Update Kernels ---
    calculate_ode_rhs_kernel.calculate_ode_rhs_kernel(geodesic_rhs, coordinate_symbols)
    rkf45_stage_update.rkf45_stage_update()
    rkf45_finalize_and_control_kernel.rkf45_finalize_and_control_kernel(
        enable_numerical_time_window_step_cap=True
    )

    # --- Event Detection & Boundary Management ---
    find_event_time_and_state.find_event_time_and_state()
    handle_source_plane_intersection.handle_source_plane_intersection()
    handle_window_plane_intersection.handle_window_plane_intersection()
    event_detection_manager_kernel.event_detection_manager_kernel()
    calculate_and_fill_blueprint_data_universal.calculate_and_fill_blueprint_data_universal()

    # --- Infrastructure Helpers ---
    # The numerical interpolation and RKF45 finalization registrations above
    # already pull in the shared slot/time-window helpers on demand. Calling
    # them again here would append duplicate definitions into BHaH_defines.h.
    batch_integrator_numerical.batch_integrator_numerical(SPACETIME)
    main.main(SPACETIME, integrator_mode)

    # --- Native NRPy Cleanup ---
    # Remove helper registrations that the event manager inlines through its
    # prefunc so the project does not emit redundant standalone C files or
    # extra prototypes for functions that are no longer called directly.
    for internal_func in [
        "find_event_time_and_state",
        "handle_source_plane_intersection",
        "handle_window_plane_intersection",
    ]:
        cfc.CFunction_dict.pop(internal_func, None)

    # ##########################################################################
    # Step 5.5: OVERRIDE DEFAULT CODE PARAMETERS
    # ##########################################################################

    # Step 5.5.a: Define the user-facing physical and numerical controls.
    #
    # This is the only block ordinary users should edit. Most entries are photon
    # CodeParameter defaults. The remaining entries define the numerical
    # spacetime dataset contract needed by the derived relationships below.
    photon_code_param_defaults: Dict[str, Union[bool, float, int]] = {
        # Execution Initial Conditions
        "t_start": 69.0,
        "scan_density": 100,
        # Batch Integrator & Numerical Limits
        "p_t_max": 100.0,
        "perform_normalization_check": True,
        "r_escape": 15.0,
        "rkf45_max_delta_t": 1.0,
        "slot_manager_delta_t": 2.0,
        # Source Plane Geometric Mapping
        "source_plane_center_x": -10.0,
        "source_plane_center_y": 0.0,
        "source_plane_center_z": 0.0,
        "source_plane_normal_x": 1.0,
        "source_plane_normal_y": 0.0,
        "source_plane_normal_z": 0.0,
        "source_r_max": 30.0,
        "source_r_min": 0.0,
        "source_up_vec_x": 0.0,
        "source_up_vec_y": 1.0,
        "source_up_vec_z": 0.0,
        # Camera Window Geometric Mapping
        "camera_pos_x": 11.0,
        "camera_pos_y": 0.0,
        "camera_pos_z": 0.0,
        "original_window_center_x": 10.0,
        "original_window_center_y": 0.0,
        "original_window_center_z": 0.0,
        "window_height": 1.0,
        "window_up_vec_x": 0.0,
        "window_up_vec_y": 0.0,
        "window_up_vec_z": 1.0,
        "window_width": 1.0,
        "window_tiles_width": 1,
        "window_tiles_height": 1,
        # RKF45 Adaptive Control Tolerances
        "numerical_initial_h": 0.05,
        "rkf45_absolute_error_tolerance": 1.0e-5,
        "rkf45_error_tolerance": 1.0e-5,
        "rkf45_h_max": 10.0,
        "rkf45_h_min": 1.0e-15,
        # Lagrange Interpolation Order
        "numerical_spacetime_spatial_interp_order": 2,
        "numerical_spacetime_temporal_interp_order": 2,
    }
    # These affect only the trusted BBH data-generation helper path.
    # They are not photon CodeParameters and therefore must not be written
    # through par.glb_code_params_dict below.
    bbh_generation_defaults: Dict[str, int] = {
        "convergence_factor": 1,
    }

    # Numerical-spacetime dataset controls used by the radial and temporal
    # contract.
    #
    # The trusted BBH example uses a fixed spherical base split of
    # [72, 12, 2]. Users no longer choose the full Nxx list directly in this
    # photon script. Instead, they choose convergence_factor above, and the
    # actual generated BBH grid is derived as
    #
    #   Nxx = [72 * convergence_factor, 12 * convergence_factor, 2].
    #
    # This keeps the pre-tested BBH angular/radial shape intact while allowing
    # controlled uniform refinement of the non-symmetry directions.
    base_bbh_nxx = [72, 12, 2]
    dt_grids = 0.25

    # eta in the radial relationship. This safety factor scales the radial
    # RKF45/time-slot overshoot estimate, but it is not a temporal endpoint
    # padding. The temporal endpoint safety is already represented by the +2
    # in (N_t + 2) * dt_grids below.
    radial_overshoot_safety_factor = 1.0

    numerical_spacetime_spatial_interp_order = int(
        par.glb_code_params_dict[
            "numerical_spacetime_spatial_interp_order"
        ].defaultvalue
    )
    numerical_spacetime_temporal_interp_order = int(
        par.glb_code_params_dict[
            "numerical_spacetime_temporal_interp_order"
        ].defaultvalue
    )

    # Step 5.5.b: Derive the numerical-spacetime time-window relationships.
    #
    # Pull the relationship-driving photon parameters out of the user-facing
    # dictionary so that the equations below read like the derivation.
    r_escape = float(photon_code_param_defaults["r_escape"])
    rkf45_max_delta_t = float(photon_code_param_defaults["rkf45_max_delta_t"])
    slot_manager_delta_t = float(photon_code_param_defaults["slot_manager_delta_t"])
    t_start_photon = float(photon_code_param_defaults["t_start"])
    convergence_factor_raw = bbh_generation_defaults["convergence_factor"]
    if not isinstance(convergence_factor_raw, (int, float)):
        raise ValueError("convergence_factor must be numeric.")
    convergence_factor_float = float(convergence_factor_raw)
    if convergence_factor_float <= 0.0 or not convergence_factor_float.is_integer():
        raise ValueError("convergence_factor must be a positive integer.")
    convergence_factor = int(convergence_factor_float)
    Nxx = [
        base_bbh_nxx[0] * convergence_factor,
        base_bbh_nxx[1] * convergence_factor,
        base_bbh_nxx[2],
    ]

    n_r = int(Nxx[0])

    # The temporal interpolation stencil has order N_t and uses 2*N_t + 1 time
    # slices. In the worst case, interpolation near a time boundary can require
    # N_t + 1 slices on one side of the photon time. We add one extra diagnostic
    # interval for boundary alignment, floating-point roundoff, and rare endpoint
    # cases in large photon batches, giving the halo (N_t + 2) * dt_grids.
    temporal_halo = (numerical_spacetime_temporal_interp_order + 2) * dt_grids

    # Photons are reverse ray traced from t_start_photon toward smaller
    # coordinate time. Therefore the upper-time boundary of the BBH data does
    # not need rkf45_max_delta_t: photons move away from that boundary after
    # initialization. It still needs the temporal interpolation halo above
    # t_start_photon so the initial interpolation stencil is available.
    bbh_data_t_final = t_start_photon + temporal_halo

    # The BBH example is treated as always starting at coordinate time 0.0.
    # Therefore the photon-side lower-time contract cannot be specified
    # independently. Instead we solve the prior relationship
    #
    #   bbh_data_t_start = slot_manager_t_min - rkf45_max_delta_t - temporal_halo
    #
    # for slot_manager_t_min while fixing bbh_data_t_start = 0.0. This makes
    # the photon stopping threshold the derived quantity and guarantees that no
    # accepted photon step requires BBH data earlier than the first available
    # slice.
    bbh_data_t_start = 0.0
    slot_manager_t_min = rkf45_max_delta_t + temporal_halo

    # The BBH runtime target is one diagnostic interval beyond the requested
    # combined-data final time so that the raytracing-data output sequence can
    # include the final requested data slice.
    bbh_runtime_t_final = bbh_data_t_final + dt_grids

    # Step 5.5.c: Compute the ceiling-aware radial grid contract.
    #
    # The radial grid must extend beyond r_escape because photons are terminated
    # only after an accepted RKF45 step. Interpolation may also be needed during
    # the step or at the accepted endpoint, so merely setting
    # grid_physical_size = r_escape is not enough.
    #
    # Spatial interpolation of order N needs a far-side radial buffer of N + 1
    # cells in the worst case. With n_r = Nxx[0] and
    # N = numerical_spacetime_spatial_interp_order, this buffer is represented
    # by radial_interp_buffer_cells.
    radial_interp_buffer_cells = numerical_spacetime_spatial_interp_order + 1

    # For an outward null ray in the approximately flat far-field region, the
    # worst-case radial coordinate motion is estimated by Delta r <= Delta t.
    # The relevant coordinate-time lookahead is modeled as the slot-manager time
    # window plus the maximum RKF45 step-time cap. The eta factor accounts for
    # stage excursions, coordinate effects, and small mismatches between the
    # ideal estimate and the actual numerical evolution.
    radial_overshoot_buffer = radial_overshoot_safety_factor * (
        slot_manager_delta_t + rkf45_max_delta_t
    )

    # The radial extent condition is
    #
    #   R >= r_escape + (N + 1) * dr + T,
    #
    # where R is grid_physical_size, dr = R / n_r, and
    #
    #   n_r = 72 * convergence_factor,
    #
    # so equivalently dr = R / (72 * convergence_factor).
    #
    # T = eta * (slot_manager_delta_t + rkf45_max_delta_t).
    #
    # The formula is implicit because dr itself depends on R. In addition, the
    # time-overshoot buffer must be represented by an integer number of radial
    # cells. The ceiling-aware solve reserves
    #
    #   k_t = ceil(((n_r - (N + 1)) * T) / (r_escape + T))
    #
    # radial cells for the RKF45/time-slot overshoot.
    radial_buffer_denominator = r_escape + radial_overshoot_buffer
    radial_buffer_numerator = (
        n_r - radial_interp_buffer_cells
    ) * radial_overshoot_buffer
    radial_time_buffer_cells = math.ceil(
        radial_buffer_numerator / radial_buffer_denominator
    )

    # After reserving N + 1 interpolation cells and k_t overshoot cells, the
    # remaining radial cells determine the minimum physical grid size:
    #
    #   R_min = n_r * r_escape / (n_r - (N + 1) - k_t).
    #
    # Guards for invalid user choices are intentionally left out here; they can
    # be added later in the final validation block.
    radial_denominator = n_r - radial_interp_buffer_cells - radial_time_buffer_cells
    grid_physical_size = n_r * r_escape / radial_denominator

    print(" -> Numerical spacetime contract:")
    for name, value in [
        ("convergence_factor", float(convergence_factor)),
        ("n_r", float(n_r)),
        ("slot_manager_t_min", slot_manager_t_min),
        ("data_t_start", bbh_data_t_start),
        ("data_t_final", bbh_data_t_final),
        ("runtime_t_final", bbh_runtime_t_final),
        ("grid_physical_size", grid_physical_size),
        ("r_escape", r_escape),
    ]:
        print(f"      {name:<18} = {value:.5f}")

    # Step 5.5.d: Define fixed metadata for the generated BBH dataset.
    combined_format_magic = "NRPYRTSTACK4D"
    combined_format_version = 1
    source_format_version = 1
    serialized_real_bytes = 8

    num_grids = 1
    payload_includes_ghost_zones = 0
    target_basis = "Cartesian"

    payload_format_name = "Cartesian g4DD+Gamma4UDD"
    payload_layout = "time_major_stage1_aos"
    payload_loop_order = "i2maj_i0fast"
    point_record_real_count = 53
    serialized_point_record_bytes = point_record_real_count * serialized_real_bytes
    metric_component_count = 10
    christoffel_component_count = 40

    spatial_lookup_mode = "coordinate_table_only"
    enable_axisymmetry = True
    axisymmetry_axis = "z"
    requires_axisymmetry_rotation = True

    bbh_floating_point_precision = "double"
    bbh_parallelization = "openmp"
    enable_raytracing_outputs = True
    bh1_mass = 0.5
    bh2_mass = 0.5
    bh1_posn_z = 0.5
    bh2_posn_z = -0.5
    gamma_driving_eta = 1.0
    outer_bc_type = "radiation"

    # Step 5.5.e: Build the required combined-bin metadata contract.
    combined_bin_location = os.path.abspath(
        os.path.join(
            args.outdir,
            "raytracing_data",
            "combined_raytracing_data.bin",
        )
    )

    required_combined_bin_metadata = {
        "combined_file": {
            "format_magic": combined_format_magic,
            "combined_format_version": combined_format_version,
            "source_format_version": source_format_version,
            "endianness": "little",
            "serialized_real_bytes": serialized_real_bytes,
        },
        "grid": {
            "CoordSystem": dataset_coord_system,
            "Nxx": Nxx,
            "num_grids": num_grids,
            "payload_includes_ghost_zones": payload_includes_ghost_zones,
            "target_basis": target_basis,
            "grid_physical_size": grid_physical_size,
        },
        "payload": {
            "format_name": payload_format_name,
            "payload_layout": payload_layout,
            "loop_order": payload_loop_order,
            "point_record_real_count": point_record_real_count,
            "point_record_bytes": serialized_point_record_bytes,
            "record_component_count": point_record_real_count,
            "metric_component_count": metric_component_count,
            "christoffel_component_count": christoffel_component_count,
        },
        "time": {
            "t_start": bbh_data_t_start,
            "t_final": bbh_data_t_final,
            "dt": dt_grids,
            "absolute_tolerance": 1.0e-12,
        },
        "spatial_lookup": {
            "spatial_lookup_mode": spatial_lookup_mode,
            "axisymmetry_enabled": enable_axisymmetry,
            "axisymmetry_axis": axisymmetry_axis,
            "requires_axisymmetry_rotation": requires_axisymmetry_rotation,
        },
        "two_blackholes_run": {
            "floating_point_precision": bbh_floating_point_precision,
            "parallelization": bbh_parallelization,
            "raytracing_outputs_enabled": enable_raytracing_outputs,
            "convergence_factor": convergence_factor,
            "BH1_mass": bh1_mass,
            "BH2_mass": bh2_mass,
            "BH1_posn_z": bh1_posn_z,
            "BH2_posn_z": bh2_posn_z,
            "GammaDriving_eta": gamma_driving_eta,
            "outer_bc_type": outer_bc_type,
            "runtime_t_final": bbh_runtime_t_final,
            "diagnostics_output_every": dt_grids,
        },
        "generation": {
            "project_name": "two_blackholes_collide",
            "python_executable": sys.executable,
            "make_command": ["make"],
            "two_blackholes_example_script": os.path.join(
                script_dir,
                "two_blackholes_collide.py",
            ),
            "combine_raytracing_time_slices_script": os.path.join(
                repo_root,
                "nrpy",
                "infrastructures",
                "BHaH",
                "diagnostics",
                "combine_raytracing_time_slices.py",
            ),
            "stage1_raytracing_output_dir": os.path.abspath(
                os.path.join(
                    repo_root,
                    "project",
                    "two_blackholes_collide",
                )
            ),
            "stage1_raytracing_pattern": "raytracing_data_t*.bin",
            "executable_name": "two_blackholes_collide",
        },
    }

    # Step 5.5.f: Ensure the required numerical spacetime data file exists.
    print(
        " -> Checking for a compatible numerical spacetime .bin cache. "
        "If validation takes more than a few seconds, the cached data were "
        "missing or stale and the two_blackholes_collide pipeline is "
        "regenerating them now."
    )
    state_of_bin, combined_bin_location = ensure_required_combined_bin(
        required_metadata=required_combined_bin_metadata,
        combined_bin_location=combined_bin_location,
    )
    print(f" -> Numerical spacetime data: {state_of_bin}")
    print(f" -> Combined raytracing data path: {combined_bin_location}")

    # Step 5.5.g: Bind the validated dataset and photon CodeParameters.
    print(" -> Overriding desired CodeParameters before .par generation...")

    par.glb_code_params_dict["numerical_spacetime_bin_path"].defaultvalue = (
        combined_bin_location
    )
    par.glb_code_params_dict[
        "numerical_spacetime_spatial_interp_order"
    ].defaultvalue = numerical_spacetime_spatial_interp_order
    par.glb_code_params_dict[
        "numerical_spacetime_temporal_interp_order"
    ].defaultvalue = numerical_spacetime_temporal_interp_order

    for param_name, default_value in photon_code_param_defaults.items():
        par.glb_code_params_dict[param_name].defaultvalue = default_value
    par.glb_code_params_dict["slot_manager_t_min"].defaultvalue = slot_manager_t_min

    # Step 6: Generate C Code for Parameter Handling
    print(" -> Generating parameter handling code...")
    CPs.write_CodeParameters_h_files(project_dir=project_dir, set_commondata_only=True)
    CPs.register_CFunctions_params_commondata_struct_set_to_default()
    cmdline_input_and_parfiles.generate_default_parfile(
        project_dir=project_dir, project_name=project_name
    )
    cmdline_input_and_parfiles.register_CFunction_cmdline_input_and_parfile_parser(
        project_name=project_name,
        cmdline_inputs=[
            "source_r_min",
            "source_r_max",
            "window_width",
            "window_height",
            "window_tiles_width",
            "window_tiles_height",
            "scan_density",
        ],
    )

    # ##########################################################################
    # Step 7: Assemble Final C Project Infrastructure
    # ##########################################################################
    print(" -> Assembling C project on disk...")

    # Map the shared CUDA-flavored kernel abstraction macros onto CPU/OpenMP.
    cpu_macros = {
        # --- Memory Access Redirection ---
        "ReadCUDA(ptr)": "#define ReadCUDA(ptr) (*(ptr))\n",
        "WriteCUDA(ptr, val)": "#define WriteCUDA(ptr, val) (*(ptr) = (val))\n",
        # --- Device Memory Fallbacks ---
        "BHAH_MALLOC_DEVICE(a, sz)": "#define BHAH_MALLOC_DEVICE(a, sz) BHAH_MALLOC(a, sz)\n",
        "BHAH_FREE_DEVICE(a)": "#define BHAH_FREE_DEVICE(a) BHAH_FREE(a)\n",
        # --- Basic Arithmetic Intrinsics (Required by RKF45 Kernels) ---
        "MulCUDA(a, b)": "#define MulCUDA(a, b) ((a) * (b))\n",
        "DivCUDA(a, b)": "#define DivCUDA(a, b) ((a) / (b))\n",
        "AddCUDA(a, b)": "#define AddCUDA(a, b) ((a) + (b))\n",
        "FusedMulAddCUDA(a, b, c)": "#define FusedMulAddCUDA(a, b, c) ((a) * (b) + (c))\n",
        # --- Hardware-Specific Math Redirection ---
        "AbsCUDA(val)": "#define AbsCUDA(val) fabs(val)\n",
        "SqrtCUDA(val)": "#define SqrtCUDA(val) sqrt(val)\n",
        "PowCUDA(a, b)": "#define PowCUDA(a, b) pow(a, b)\n",
        # --- Function Decorators ---
        "BHAH_HD_FUNC": "#define BHAH_HD_FUNC\n",
        "BHAH_HD_INLINE": "#define BHAH_HD_INLINE static inline\n",
        # --- Parallelization & Scope Helpers ---
        "BHAH_WARP_ATOMIC_ADD(ptr, val)": '#define BHAH_WARP_ATOMIC_ADD(ptr, val) _Pragma("omp atomic") *(ptr) += (val)\n',
        "GLOBAL_COMMONDATA_EXTERN": "// CPU passes commondata by reference, no global needed.\n",
        "BHAH_DEVICE_SYNC()": "#define BHAH_DEVICE_SYNC() do {} while(0)\n",
    }

    # Match the combined numerical spacetime payload layout produced by the
    # trusted two-black-hole raytracing pipeline, which stores NGHOSTS = 3.
    BHaH_defines_h.output_BHaH_defines_h(
        project_dir=project_dir,
        enable_rfm_precompute=False,
        fin_NGHOSTS_add_one_for_upwinding_or_KO=True,
        supplemental_defines_dict=cpu_macros,
    )

    compiler = "gcc"
    cflags = ["-fopenmp", "-O3", "-fno-omit-frame-pointer", "-DDEBUG", "-Wno-stringop-truncation"]
    libs = ["-lm"]
    ext = "c"

    print(" -> Generating Makefile ")

    # Set optimization option string
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
    # ##########################################################################
    # PART 2: FINALIZE
    # ##########################################################################

    # Define the directory containing the visualization assets relative to the repository root
    vis_dir = os.path.join(
        "nrpy", "examples", "geodesic_helpers", "geodesic_visualizations"
    )

    # Locate the visualization script
    vis_script_src = os.path.join(vis_dir, "visualize_lensed_image.py")
    config_src = os.path.join(vis_dir, "blueprint_config_and_schema.py")
    render_src = os.path.join(vis_dir, "render_lensed_image.py")
    blueprint_analysis_src = os.path.join(vis_dir, "blueprint_analysis.py")

    for script_src in (
        vis_script_src,
        config_src,
        render_src,
        blueprint_analysis_src,
    ):
        shutil.copy(script_src, project_dir)

    # The inner disk radius ensures the texture mapping aligns with the computed initial conditions.
    c_r_min = float(par.glb_code_params_dict["source_r_min"].defaultvalue)

    # The outer disk radius ensures the mathematical bounds match the physical source geometry.
    c_r_max = float(par.glb_code_params_dict["source_r_max"].defaultvalue)

    # The camera window width scales the horizontal projection bounds.
    c_window_width = float(par.glb_code_params_dict["window_width"].defaultvalue)

    # The camera window height scales the vertical projection bounds.
    c_window_height = float(par.glb_code_params_dict["window_height"].defaultvalue)

    # The terminal execution string now references the local copy of the script.
    # Get tile dimensions to map the partitioned geometry blocks.
    c_tiles_width = int(par.glb_code_params_dict["window_tiles_width"].defaultvalue)
    c_tiles_height = int(par.glb_code_params_dict["window_tiles_height"].defaultvalue)

    # Set the baseline pixel resolution for the output render.
    c_pixel_width = 600

    vis_command = (
        f"python3 visualize_lensed_image.py "
        f"--source_r_min {c_r_min} "
        f"--source_r_max {c_r_max} "
        f"--window_width {c_window_width} "
        f"--window_height {c_window_height} "
        f"--window_tiles_width {c_tiles_width} "
        f"--window_tiles_height {c_tiles_height} "
        f"--pixel_width {c_pixel_width}"
    )

    print(
        f"Finished! Now go into {project_dir} and type `make` to build, then ./{exec_name} to run."
    )
    print(
        f"    Parameter file can be found at {os.path.join(project_dir, f'{project_name}.par')}\n"
    )
    print(
        "    To generate the lensed image after running the C executable, ensure you have the required Python packages:"
    )
    print("    pip install matplotlib numpy numba Pillow\n")
    print(
        "    Then, execute the visualization script directly from the project directory:"
    )
    print(f"    {vis_command}\n")

    blueprint_command = (
        f"python3 blueprint_analysis.py "
        f"--window_tiles_width {c_tiles_width} "
        f"--window_tiles_height {c_tiles_height} "
        f"--window_width {c_window_width} "
        f"--window_height {c_window_height}"
    )

    print("    To run the blueprint diagnostic:")
    print(f"    {blueprint_command}\n")

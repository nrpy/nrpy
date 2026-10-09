"""
Configure optional raytracing output for BHaH evolution examples.

This module stores example defaults, adds command-line options, validates the
selected grid and output times, and writes the files used to combine raytracing
time slices.

Author: Dalton J. Moone
"""

import argparse
import json
import math
import os
import shlex
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Mapping, Optional, Sequence, Tuple, Union

import nrpy.params as par
from nrpy.helpers.generic import copy_files

NameValue = Union[str, int, float]
RAYTRACING_MODES = ("g4DD", "g4DD_d0", "GammaUDD")
SLICE_PATTERN = "raytracing_data_t????????.bin"


@dataclass(frozen=True)
class RaytracingDefaults:
    """Default evolution and raytracing values for one example."""

    normal_coord_system: str
    raytracing_coord_system: str
    supported_coord_systems: Tuple[str, ...]
    normal_t_final: float
    normal_output_every: float
    raytracing_t_final: float
    raytracing_output_every: float
    normal_grid_physical_size: float
    raytracing_domains: Mapping[str, Tuple[float, ...]]
    raytracing_nxx: Mapping[str, Tuple[int, int, int]]


@dataclass(frozen=True)
class RaytracingOptions:
    """Validated values used for coordinate selection and diagnostics output."""

    enabled: bool
    data_mode: str
    static_christoffels: bool
    coord_system: str
    t_final: float
    output_every: float
    grid_physical_size: float
    nxx: Optional[Tuple[int, int, int]]
    domain: Optional[Tuple[float, ...]]
    example_name: str = ""
    initial_sep: float = 0.5
    initial_p_r: float = 0.0
    bh_positions: Tuple[float, float] = (0.5, -0.5)
    bh_masses: Tuple[float, float] = (0.5, 0.5)
    formulation: str = "BSSN"


@dataclass(frozen=True)
class RunFiles:
    """Paths written for one generated raytracing project."""

    combiner_path: Path
    metadata_path: Path
    shell_path: Path
    combined_paths: Mapping[str, Path]


EXAMPLE_DEFAULTS = {
    "blackhole_spectroscopy": RaytracingDefaults(
        normal_coord_system="SinhCylindrical",
        raytracing_coord_system="SinhCylindrical",
        supported_coord_systems=("SinhCylindrical", "SinhCylindricalv2n2"),
        normal_t_final=450.0,
        normal_output_every=0.5,
        raytracing_t_final=450.0,
        raytracing_output_every=0.5,
        normal_grid_physical_size=300.0,
        raytracing_domains={
            "SinhCylindrical": (300.0, 0.2, 0.2),
            "SinhCylindricalv2n2": (30.0, 0.075, 0.05, 1.0, 4.0),
        },
        raytracing_nxx={
            "SinhCylindrical": (400, 2, 1200),
            "SinhCylindricalv2n2": (162, 2, 256),
        },
    ),
    "two_blackholes_collide": RaytracingDefaults(
        normal_coord_system="Spherical",
        raytracing_coord_system="Spherical",
        supported_coord_systems=("Spherical", "SinhCylindricalv2n2"),
        normal_t_final=7.5,
        normal_output_every=0.25,
        raytracing_t_final=7.5,
        raytracing_output_every=0.25,
        normal_grid_physical_size=7.5,
        raytracing_domains={
            "Spherical": (7.5,),
            "SinhCylindricalv2n2": (30.0, 0.075, 0.05, 1.0, 4.0),
        },
        raytracing_nxx={
            "Spherical": (72, 12, 2),
            "SinhCylindricalv2n2": (162, 2, 256),
        },
    ),
}


def _defaults_for(example_name: str) -> RaytracingDefaults:
    """
    Select raytracing defaults for a supported example.

    :param example_name: Example module name without the ``.py`` suffix.
    :return: Defaults for the example.
    :raises ValueError: If the example has no raytracing setup.
    """
    if example_name not in EXAMPLE_DEFAULTS:
        raise ValueError(f"Unsupported raytracing example: {example_name}.")
    return EXAMPLE_DEFAULTS[example_name]


def add_raytracing_cli(parser: argparse.ArgumentParser, example_name: str) -> None:
    """
    Add common and example-specific raytracing options to a parser.

    :param parser: Argument parser owned by the example.
    :param example_name: Supported example module name without ``.py``.
    """
    _add_raytracing_arguments(parser, defaults=_defaults_for(example_name))
    if example_name == "blackhole_spectroscopy":
        parser.add_argument(
            "--initial-sep",
            type=float,
            default=None,
            help="Set the TwoPunctures separation in units of total mass.",
        )
        parser.add_argument(
            "--initial-p-r",
            type=float,
            default=None,
            help="Set radial puncture momentum; -1 is reserved for NRPyPN.",
        )
    else:
        parser.add_argument(
            "--raytracing-bhs",
            nargs=4,
            type=float,
            default=None,
            metavar=("Z_1", "Z_2", "M_1", "M_2"),
            help="Set Brill-Lindquist puncture positions and masses.",
        )


def check_raytracing_cli(
    args: argparse.Namespace, example_name: str
) -> RaytracingOptions:
    """
    Validate CLI values and select evolution parameters for one example.

    :param args: Parsed command-line arguments.
    :param example_name: Supported example module name without ``.py``.
    :return: Validated coordinate, grid, output, and initial-data values.
    :raises ValueError: If CLI values conflict or contain invalid physical values.
    """
    options = _resolve_raytracing_options(
        args,
        defaults=_defaults_for(example_name),
        parallelization="cuda" if args.cuda else "openmp",
        fp_type=args.floating_point_precision.lower(),
    )
    if example_name == "blackhole_spectroscopy":
        initial_sep = 0.5 if args.initial_sep is None else args.initial_sep
        initial_p_r = 0.0 if args.initial_p_r is None else args.initial_p_r
        if not math.isfinite(initial_sep) or initial_sep <= 0.0:
            raise ValueError("--initial-sep must be finite and positive.")
        if not math.isfinite(initial_p_r) or initial_p_r == -1.0:
            raise ValueError("--initial-p-r must be finite and different from -1.")
        return replace(
            options,
            example_name=example_name,
            initial_sep=initial_sep,
            initial_p_r=initial_p_r,
            bh_positions=(0.5 * initial_sep, -0.5 * initial_sep),
            formulation="fCCZ4" if args.fccz4 else "BSSN",
        )
    black_holes = (
        (0.5, -0.5, 0.5, 0.5)
        if args.raytracing_bhs is None
        else tuple(args.raytracing_bhs)
    )
    if not all(math.isfinite(value) for value in black_holes):
        raise ValueError("--raytracing-bhs values must be finite.")
    if black_holes[2] <= 0.0 or black_holes[3] <= 0.0:
        raise ValueError("--raytracing-bhs masses must be positive.")
    return replace(
        options,
        example_name=example_name,
        bh_positions=(black_holes[0], black_holes[1]),
        bh_masses=(black_holes[2], black_holes[3]),
    )


def _add_raytracing_arguments(
    parser: argparse.ArgumentParser, *, defaults: RaytracingDefaults
) -> None:
    """
    Add shared raytracing options to an example's parser.

    :param parser: Argument parser owned by the example.
    :param defaults: Supported coordinate systems for this example.
    """
    parser.add_argument(
        "--raytracing-time",
        nargs="*",
        type=float,
        default=None,
        metavar="VALUE",
        help="Write raytracing slices; optionally provide T_FINAL OUTPUT_EVERY.",
    )
    parser.add_argument(
        "--raytracing-data-mode",
        choices=(*RAYTRACING_MODES, "all"),
        default="g4DD",
        help="Select metric data written at each output time.",
    )
    parser.add_argument(
        "--raytracing-static-christoffels",
        action="store_true",
        help="Use static Christoffels for the final GammaUDD output.",
    )
    parser.add_argument(
        "--raytracing-coord-system",
        choices=defaults.supported_coord_systems,
        default=None,
        help="Select the coordinate system for raytracing output.",
    )
    parser.add_argument(
        "--raytracing-domain",
        nargs="+",
        type=float,
        default=None,
        metavar="VALUE",
        help="Set GRID_PHYSICAL_SIZE [SINHWRHO SINHWZ [RHO_SLOPE Z_SLOPE]].",
    )
    parser.add_argument(
        "--raytracing-Nxx",
        dest="raytracing_nxx",
        nargs=3,
        type=int,
        default=None,
        metavar=("NXX0", "NXX1", "NXX2"),
        help="Set three base grid dimensions; cylindrical grids need NXX1=2.",
    )


def _resolve_raytracing_options(
    args: argparse.Namespace,
    *,
    defaults: RaytracingDefaults,
    parallelization: str,
    fp_type: str,
) -> RaytracingOptions:
    """
    Resolve and validate raytracing values before C code generation.

    :param args: Parsed example command-line arguments.
    :param defaults: Normal and raytracing defaults for the example.
    :param parallelization: Selected BHaH host or CUDA build.
    :param fp_type: Selected floating-point type.
    :return: Values to use for the grid, evolution time, and diagnostics.
    :raises ValueError: If raytracing options or build settings are incompatible.
    """
    if fp_type not in ("float", "double"):
        raise ValueError("--floating_point_precision must be float or double.")
    if parallelization not in ("openmp", "cuda"):
        raise ValueError("parallelization must be openmp or cuda.")

    enabled = args.raytracing_time is not None or getattr(
        args, "raytracing_outputs", False
    )
    if not enabled:
        dependent_options = (
            ("--raytracing-data-mode", args.raytracing_data_mode != "g4DD"),
            ("--raytracing-static-christoffels", args.raytracing_static_christoffels),
            ("--raytracing-coord-system", args.raytracing_coord_system is not None),
            ("--raytracing-domain", args.raytracing_domain is not None),
            ("--raytracing-Nxx", args.raytracing_nxx is not None),
        )
        for option_name, supplied in dependent_options:
            if supplied:
                raise ValueError(f"{option_name} requires --raytracing-time.")
        return RaytracingOptions(
            enabled=False,
            data_mode="g4DD",
            static_christoffels=False,
            coord_system=defaults.normal_coord_system,
            t_final=defaults.normal_t_final,
            output_every=defaults.normal_output_every,
            grid_physical_size=defaults.normal_grid_physical_size,
            nxx=None,
            domain=None,
        )

    if parallelization != "openmp" or fp_type != "double":
        raise ValueError("--raytracing-time requires OpenMP and double precision.")
    if args.raytracing_static_christoffels and args.raytracing_data_mode not in (
        "GammaUDD",
        "all",
    ):
        raise ValueError("--raytracing-static-christoffels requires GammaUDD or all.")
    if args.raytracing_time is None or len(args.raytracing_time) == 0:
        t_final = defaults.raytracing_t_final
        output_every = defaults.raytracing_output_every
    elif len(args.raytracing_time) == 2:
        t_final, output_every = args.raytracing_time
    else:
        raise ValueError("--raytracing-time requires zero or two values.")
    if not all(
        math.isfinite(value) and value > 0.0 for value in (t_final, output_every)
    ):
        raise ValueError("--raytracing-time values must be finite and positive.")

    coord_system = args.raytracing_coord_system or defaults.raytracing_coord_system
    if coord_system not in defaults.supported_coord_systems:
        raise ValueError(f"Unsupported raytracing coordinate system: {coord_system}.")
    if coord_system not in defaults.raytracing_domains:
        raise ValueError(f"No --raytracing-domain default for {coord_system}.")
    if coord_system not in defaults.raytracing_nxx:
        raise ValueError(f"No --raytracing-Nxx default for {coord_system}.")
    domain = tuple(
        args.raytracing_domain
        if args.raytracing_domain is not None
        else defaults.raytracing_domains[coord_system]
    )
    expected_length = {
        "Spherical": 1,
        "SinhCylindrical": 3,
        "SinhCylindricalv2n2": 5,
    }.get(coord_system)
    if expected_length is None:
        raise ValueError(f"No raytracing domain layout for {coord_system}.")
    if len(domain) != expected_length:
        raise ValueError(
            f"--raytracing-domain needs {expected_length} values for {coord_system}."
        )
    if not all(math.isfinite(value) and value > 0.0 for value in domain):
        raise ValueError("--raytracing-domain values must be finite and positive.")
    nxx = tuple(
        args.raytracing_nxx
        if args.raytracing_nxx is not None
        else defaults.raytracing_nxx[coord_system]
    )
    if len(nxx) != 3 or any(value <= 0 for value in nxx):
        raise ValueError("--raytracing-Nxx needs three positive values.")
    if "Cylindrical" in coord_system and nxx[1] != 2:
        raise ValueError("--raytracing-Nxx needs NXX1=2 for cylindrical grids.")
    return RaytracingOptions(
        enabled=True,
        data_mode=args.raytracing_data_mode,
        static_christoffels=args.raytracing_static_christoffels,
        coord_system=coord_system,
        t_final=t_final,
        output_every=output_every,
        grid_physical_size=domain[0],
        nxx=nxx,
        domain=domain,
    )


def _name_token(value: NameValue) -> str:
    """
    Convert one physical parameter value to a file-name component.

    :param value: Parameter name or numerical value.
    :return: Component containing only letters, digits, underscores, or hyphens.
    :raises ValueError: If a component is empty or contains a path character.
    """
    token = (
        str(value).replace("-", "neg").replace(".", "p")
        if isinstance(value, (int, float))
        else value
    )
    if not token or not all(
        char.isascii() and (char.isalnum() or char in "_-") for char in token
    ):
        raise ValueError(f"Invalid combined-file name component: {value!r}.")
    return token


def _write_raytracing_run_files(
    *,
    project_dir: Union[str, Path],
    project_name: str,
    generator_script: str,
    executable_name: str,
    options: RaytracingOptions,
    initial_data: Mapping[str, object],
    name_values: Sequence[Tuple[str, NameValue]],
) -> Optional[RunFiles]:
    """
    Write raytracing metadata and a build, run, and combine shell script.

    :param project_dir: Generated BHaH project directory.
    :param project_name: Name recorded in the combined-file metadata.
    :param generator_script: Python example filename recorded in metadata.
    :param executable_name: Generated C executable launched by the shell script.
    :param options: Validated raytracing coordinate and output values.
    :param initial_data: Example-specific initial-data values for JSON metadata.
    :param name_values: Ordered initial-data values for combined-file names.
    :return: Written file paths, or None when raytracing is disabled.
    :raises ValueError: If a filename component or enabled option is invalid.
    :raises FileNotFoundError: If the generated project directory is missing.
    """
    if not options.enabled:
        return None
    if options.domain is None or options.nxx is None:
        raise ValueError("Enabled raytracing output needs a domain and Nxx.")
    if options.data_mode not in (*RAYTRACING_MODES, "all"):
        raise ValueError(f"Unsupported raytracing data mode: {options.data_mode}.")
    project_path = Path(project_dir)
    if not project_path.is_dir():
        raise FileNotFoundError(
            f"Generated project directory is missing: {project_path}."
        )
    for file_name in (project_name, executable_name):
        _name_token(file_name)
    if Path(generator_script).name != generator_script:
        raise ValueError("generator_script must be a Python filename, not a path.")
    if not generator_script.endswith(".py"):
        raise ValueError("generator_script must end in .py.")
    _name_token(generator_script[:-3])

    modes = RAYTRACING_MODES if options.data_mode == "all" else (options.data_mode,)
    stem_parts = [project_name]
    for name, value in name_values:
        stem_parts.extend((_name_token(name), _name_token(value)))
    stem_parts.extend(
        (
            "tf",
            _name_token(options.t_final),
            "dt",
            _name_token(options.output_every),
            _name_token(options.coord_system),
            "domain",
            *(_name_token(value) for value in options.domain),
            *(str(value) for value in options.nxx),
        )
    )
    if options.static_christoffels:
        stem_parts.append("staticGamma")
    stem = "_".join(stem_parts)
    combined_paths = {
        mode: project_path.parent / "raytracing_data" / f"{stem}_{mode}.bin"
        for mode in modes
    }
    metadata = {
        "generator_script": generator_script,
        "project_name": project_name,
        "initial_data": dict(initial_data),
        "grid_physical_size": options.grid_physical_size,
        "raytracing_coord_system": options.coord_system,
        "raytracing_domain": list(options.domain),
        "raytracing_Nxx": list(options.nxx),
        "raytracing_time": {
            "t_final": options.t_final,
            "diagnostics_output_every": options.output_every,
        },
        "raytracing_data_mode": options.data_mode,
        "raytracing_static_christoffels": options.static_christoffels,
        "slice_filename_pattern": SLICE_PATTERN,
        "datasets": {mode: path.name for mode, path in combined_paths.items()},
    }
    if options.data_mode == "all":
        metadata["generated_datasets"] = list(modes)
    else:
        metadata["combined_output_filename"] = combined_paths[options.data_mode].name

    copy_files(
        package="nrpy.infrastructures.BHaH.diagnostics",
        filenames_list=["combine_raytracing_time_slices.py"],
        project_dir=str(project_path),
        subdirectory="",
    )
    combiner_path = project_path / "combine_raytracing_time_slices.py"
    metadata_path = project_path / "raytracing_run_metadata.json"
    metadata_path.write_text(
        json.dumps(metadata, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )

    shell_lines = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        'cd -- "$(dirname -- "${BASH_SOURCE[0]}")"',
        "make",
        f"./{executable_name}",
        "mkdir -p ../raytracing_data",
    ]
    for mode, output_path in combined_paths.items():
        shell_lines.extend(
            (
                "python3 combine_raytracing_time_slices.py \\",
                f"  --input-dir {shlex.quote('raytracing_slices/' + mode)} \\",
                f"  --pattern {shlex.quote(SLICE_PATTERN)} \\",
                "  --run-metadata raytracing_run_metadata.json \\",
                f"  --output {shlex.quote(os.path.join('..', 'raytracing_data', output_path.name))} --force",
                f"echo {shlex.quote(os.path.join('..', 'raytracing_data', output_path.name))}",
            )
        )
    shell_path = project_path / "run_raytracing_data_pipeline.sh"
    shell_path.write_text("\n".join(shell_lines) + "\n", encoding="utf-8")
    shell_path.chmod(0o755)
    return RunFiles(combiner_path, metadata_path, shell_path, combined_paths)


def adjust_raytracing_code_parameters(options: RaytracingOptions) -> None:
    """
    Set reference-metric CodeParameters used by raytracing output.

    :param options: Validated raytracing coordinate and domain values.
    """
    if not options.enabled:
        return
    assert options.domain is not None
    if options.coord_system == "Spherical":
        return
    par.adjust_CodeParam_default("SINHWRHO", options.domain[1])
    par.adjust_CodeParam_default("SINHWZ", options.domain[2])
    if options.coord_system == "SinhCylindricalv2n2":
        par.adjust_CodeParam_default("rho_slope", options.domain[3])
        par.adjust_CodeParam_default("z_slope", options.domain[4])


def shell_generator(options: RaytracingOptions, project_dir: Union[str, Path]) -> str:
    """
    Write the build, evolution, and slice-combination commands.

    :param options: Validated values selected by ``check_raytracing_cli``.
    :param project_dir: Generated BHaH project directory.
    :return: Command to run from the generated project directory.
    """
    example_name = options.example_name
    _defaults_for(example_name)
    normal_command = (
        f"Finished! Now go into project/{example_name} and type `make` to build, "
        f"then ./{example_name} to run."
    )
    if not options.enabled:
        return normal_command
    name_values: Sequence[Tuple[str, NameValue]]
    if example_name == "blackhole_spectroscopy":
        initial_data = {
            "type": "TwoPunctures",
            "orientation": "legacy_swap_xz",
            "evolution_formulation": options.formulation,
            "initial_sep": options.initial_sep,
            "initial_p_r": options.initial_p_r,
            "initial_p_t": 0.0,
            "mass_ratio": 1.0,
            "target_adm_mass_M": 0.5,
            "target_adm_mass_m": 0.5,
            "TP_npoints_A": 48,
            "TP_npoints_B": 48,
            "TP_npoints_phi": 4,
            "bbhxy_BH_m_chix": 0.0,
            "bbhxy_BH_M_chix": 0.0,
        }
        name_values = (
            ("formulation", options.formulation),
            ("TP_sep", options.initial_sep),
            ("pr", options.initial_p_r),
            ("q", 1.0),
        )
    else:
        initial_data = {
            "type": "BrillLindquist",
            "BH1_mass": options.bh_masses[0],
            "BH1_posn_z": options.bh_positions[0],
            "BH2_mass": options.bh_masses[1],
            "BH2_posn_z": options.bh_positions[1],
        }
        name_values = (
            ("z1", options.bh_positions[0]),
            ("z2", options.bh_positions[1]),
            ("M1", options.bh_masses[0]),
            ("M2", options.bh_masses[1]),
        )
    run_files = _write_raytracing_run_files(
        project_dir=project_dir,
        project_name=example_name,
        generator_script=f"{example_name}.py",
        executable_name=example_name,
        options=options,
        initial_data=initial_data,
        name_values=name_values,
    )
    assert run_files is not None
    output_lines = [
        f"Finished! Go to {run_files.shell_path.parent} and run ./run_raytracing_data_pipeline.sh."
    ]
    output_lines.extend(
        f"    Combined {mode} data: {path}"
        for mode, path in run_files.combined_paths.items()
    )
    return "\n".join(output_lines)

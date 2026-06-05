"""
Helper functions for validating and generating combined raytracing bins.

This module provides a small helper layer for numerical geodesic examples that
need a combined spacetime ``.bin`` file before code generation proceeds. The
helper checks whether a cached combined file exists, validates only the
lightweight metadata needed for compatibility, regenerates the file through the
existing two-black-hole example when needed, and writes a small JSON sidecar
manifest next to the generated ``.bin``.

The helper intentionally contains no photon-specific logic. Its only job is to
ensure that the required combined numerical spacetime data file exists and
matches the expected metadata contract.

Author: OpenAI
        support **at** openai **dot** com
"""

import json
import math
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, cast

REPO_ROOT = Path(__file__).resolve().parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from nrpy.infrastructures.BHaH.diagnostics import combine_raytracing_time_slices as crt

STATE_ALREADY_CORRECT = ".bin already existed and was correct"
STATE_REBUILT_STALE = (
    ".bin already existed but meta data wasn't correct, so correct .bin was generated"
)
STATE_GENERATED_MISSING = ".bin file didn't exist and was generated"

MANIFEST_VERSION = 1
DEFAULT_POINT_RECORD_REAL_COUNT = 53
DEFAULT_POINT_RECORD_BYTES = 424
DEFAULT_RECORD_COMPONENT_COUNT = 53
DEFAULT_METRIC_COMPONENT_COUNT = 10
DEFAULT_CHRISTOFFEL_COMPONENT_COUNT = 40
DEFAULT_TWO_BLACKHOLES_PROJECT_NAME = "two_blackholes_collide"
DEFAULT_STAGE1_PATTERN = "raytracing_data_t*.bin"


def _as_float(value: object) -> float:
    if not isinstance(value, (int, float)):
        raise ValueError("Expected a numeric metadata value.")
    return float(value)


def _normalize_required_metadata(
    required_metadata: Dict[str, object],
) -> Dict[str, object]:
    combined_file = cast(Dict[str, object], required_metadata["combined_file"])
    grid = cast(Dict[str, object], required_metadata["grid"])
    payload = cast(Dict[str, object], required_metadata["payload"])
    time = cast(Dict[str, object], required_metadata["time"])
    spatial_lookup = cast(Dict[str, object], required_metadata["spatial_lookup"])
    two_blackholes_run = cast(
        Dict[str, object], required_metadata["two_blackholes_run"]
    )
    generation = cast(Dict[str, object], required_metadata["generation"])

    project_name = cast(
        str,
        generation.get("project_name", DEFAULT_TWO_BLACKHOLES_PROJECT_NAME),
    )
    project_dir = (REPO_ROOT / "project" / project_name).resolve()
    stage1_dir = Path(
        cast(str, generation.get("stage1_raytracing_output_dir", str(project_dir)))
    ).expanduser()
    if not stage1_dir.is_absolute():
        stage1_dir = (REPO_ROOT / stage1_dir).resolve()
    else:
        stage1_dir = stage1_dir.resolve()

    normalized = {
        "combined_file": {
            "format_magic": cast(
                str,
                combined_file.get(
                    "format_magic",
                    crt.COMBINED_MAGIC.split(b"\0", 1)[0].decode("ascii"),
                ),
            ),
            "combined_format_version": cast(
                int,
                combined_file.get(
                    "combined_format_version",
                    crt.COMBINED_FORMAT_VERSION,
                ),
            ),
            "source_format_version": cast(
                int,
                combined_file.get("source_format_version", crt.STAGE1_FORMAT_VERSION),
            ),
            "endianness": cast(str, combined_file.get("endianness", "little")),
            "serialized_real_bytes": cast(
                int,
                combined_file.get("serialized_real_bytes", crt.SERIALIZED_REAL_BYTES),
            ),
        },
        "grid": {
            "CoordSystem": cast(str, grid.get("CoordSystem", "Spherical")),
            "Nxx": cast(List[int], grid.get("Nxx", [72, 12, 2])),
            "num_grids": cast(int, grid.get("num_grids", 1)),
            "payload_includes_ghost_zones": cast(
                int,
                grid.get("payload_includes_ghost_zones", 0),
            ),
            "target_basis": cast(str, grid.get("target_basis", "Cartesian")),
            "grid_physical_size": _as_float(grid.get("grid_physical_size", 7.5)),
        },
        "payload": {
            "format_name": cast(
                str,
                payload.get("format_name", "Cartesian g4DD+Gamma4UDD"),
            ),
            "payload_layout": cast(
                str,
                payload.get("payload_layout", "time_major_stage1_aos"),
            ),
            "loop_order": cast(str, payload.get("loop_order", "i2maj_i0fast")),
            "point_record_real_count": cast(
                int,
                payload.get(
                    "point_record_real_count",
                    DEFAULT_POINT_RECORD_REAL_COUNT,
                ),
            ),
            "point_record_bytes": cast(
                int,
                payload.get("point_record_bytes", DEFAULT_POINT_RECORD_BYTES),
            ),
            "record_component_count": cast(
                int,
                payload.get("record_component_count", DEFAULT_RECORD_COMPONENT_COUNT),
            ),
            "metric_component_count": cast(
                int,
                payload.get("metric_component_count", DEFAULT_METRIC_COMPONENT_COUNT),
            ),
            "christoffel_component_count": cast(
                int,
                payload.get(
                    "christoffel_component_count",
                    DEFAULT_CHRISTOFFEL_COMPONENT_COUNT,
                ),
            ),
        },
        "time": {
            "t_start": _as_float(time.get("t_start", 0.0)),
            "t_final": _as_float(time.get("t_final", 7.5)),
            "absolute_tolerance": _as_float(time.get("absolute_tolerance", 1.0e-12)),
            "t_final_tolerance": _as_float(
                time.get(
                    "t_final_tolerance",
                    (
                        0.6 * _as_float(time["dt"])
                        if "dt" in time
                        else time.get("absolute_tolerance", 1.0e-12)
                    ),
                )
            ),
        },
        "spatial_lookup": {
            "spatial_lookup_mode": cast(
                str,
                spatial_lookup.get("spatial_lookup_mode", "coordinate_table_only"),
            ),
            "axisymmetry_enabled": cast(
                bool,
                spatial_lookup.get("axisymmetry_enabled", True),
            ),
            "axisymmetry_axis": cast(
                str,
                spatial_lookup.get("axisymmetry_axis", "z"),
            ),
            "requires_axisymmetry_rotation": cast(
                bool,
                spatial_lookup.get("requires_axisymmetry_rotation", True),
            ),
        },
        "two_blackholes_run": {
            "floating_point_precision": cast(
                str,
                two_blackholes_run.get("floating_point_precision", "double"),
            ),
            "parallelization": cast(
                str,
                two_blackholes_run.get("parallelization", "openmp"),
            ),
            "raytracing_outputs_enabled": cast(
                bool,
                two_blackholes_run.get("raytracing_outputs_enabled", True),
            ),
            "BH1_mass": _as_float(two_blackholes_run.get("BH1_mass", 0.5)),
            "BH2_mass": _as_float(two_blackholes_run.get("BH2_mass", 0.5)),
            "BH1_posn_z": _as_float(two_blackholes_run.get("BH1_posn_z", 0.5)),
            "BH2_posn_z": _as_float(two_blackholes_run.get("BH2_posn_z", -0.5)),
            "GammaDriving_eta": _as_float(
                two_blackholes_run.get("GammaDriving_eta", 1.0)
            ),
            "outer_bc_type": cast(
                str,
                two_blackholes_run.get("outer_bc_type", "radiation"),
            ),
            "runtime_t_final": _as_float(
                two_blackholes_run.get(
                    "runtime_t_final",
                    _as_float(time.get("t_final", 7.5)),
                )
            ),
            "diagnostics_output_every": _as_float(
                two_blackholes_run.get("diagnostics_output_every", 0.25)
            ),
        },
        "generation": {
            "python_executable": cast(
                str, generation.get("python_executable", sys.executable)
            ),
            "make_command": cast(List[str], generation.get("make_command", ["make"])),
            "two_blackholes_example_script": str(
                Path(
                    cast(
                        str,
                        generation.get(
                            "two_blackholes_example_script",
                            str(
                                REPO_ROOT
                                / "nrpy"
                                / "examples"
                                / "two_blackholes_collide.py"
                            ),
                        ),
                    )
                )
                .expanduser()
                .resolve()
            ),
            "combine_raytracing_time_slices_script": str(
                Path(
                    cast(
                        str,
                        generation.get(
                            "combine_raytracing_time_slices_script",
                            str(
                                REPO_ROOT
                                / "nrpy"
                                / "infrastructures"
                                / "BHaH"
                                / "diagnostics"
                                / "combine_raytracing_time_slices.py"
                            ),
                        ),
                    )
                )
                .expanduser()
                .resolve()
            ),
            "project_dir": str(project_dir),
            "stage1_raytracing_output_dir": str(stage1_dir),
            "stage1_raytracing_pattern": cast(
                str,
                generation.get("stage1_raytracing_pattern", DEFAULT_STAGE1_PATTERN),
            ),
            "executable_name": cast(
                str, generation.get("executable_name", project_name)
            ),
        },
    }
    _validate_supported_generation(normalized)
    return normalized


def _validate_supported_generation(normalized: Dict[str, object]) -> None:
    grid = cast(Dict[str, object], normalized["grid"])
    run = cast(Dict[str, object], normalized["two_blackholes_run"])
    if grid["CoordSystem"] != "Spherical":
        raise ValueError("This helper currently supports only CoordSystem='Spherical'.")
    if run["parallelization"] != "openmp":
        raise ValueError(
            "This helper currently supports only parallelization='openmp'."
        )
    if run["floating_point_precision"] != "double":
        raise ValueError(
            "This helper currently supports only floating_point_precision='double'."
        )
    if run["raytracing_outputs_enabled"] is not True:
        raise ValueError("raytracing_outputs_enabled must be True.")
    if run["outer_bc_type"] != "radiation":
        raise ValueError(
            "This helper currently supports only outer_bc_type='radiation'."
        )


def _manifest_path(combined_bin_path: Path) -> Path:
    return combined_bin_path.with_suffix(combined_bin_path.suffix + ".manifest.json")


def _manifest_metadata(normalized: Dict[str, object]) -> Dict[str, object]:
    return {
        "combined_file": normalized["combined_file"],
        "grid": normalized["grid"],
        "payload": normalized["payload"],
        "time": normalized["time"],
        "spatial_lookup": normalized["spatial_lookup"],
        "two_blackholes_run": normalized["two_blackholes_run"],
    }


def _manifest_matches(manifest_path: Path, normalized: Dict[str, object]) -> bool:
    if not manifest_path.is_file():
        return False
    with manifest_path.open("r", encoding="utf-8") as file_pointer:
        manifest = json.load(file_pointer)
    if not isinstance(manifest, dict):
        return False
    if manifest.get("manifest_version") != MANIFEST_VERSION:
        return False
    if manifest.get("cache_metadata") != _manifest_metadata(normalized):
        return False
    return True


def _write_manifest(manifest_path: Path, normalized: Dict[str, object]) -> None:
    manifest = {
        "manifest_version": MANIFEST_VERSION,
        "cache_metadata": _manifest_metadata(normalized),
    }
    temp_fd, temp_name = tempfile.mkstemp(
        prefix=f"{manifest_path.name}.tmp.",
        dir=str(manifest_path.parent),
    )
    temp_path = Path(temp_name)
    try:
        with os.fdopen(temp_fd, "w", encoding="utf-8") as file_pointer:
            json.dump(manifest, file_pointer, indent=2, sort_keys=True)
            file_pointer.write("\n")
            file_pointer.flush()
            os.fsync(file_pointer.fileno())
        os.replace(temp_path, manifest_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _parse_combined_bin_metadata(combined_bin_path: Path) -> Dict[str, object]:
    file_size = combined_bin_path.stat().st_size
    with combined_bin_path.open("rb") as file_pointer:
        header = crt._unpack_fixed_header(  # pylint: disable=protected-access
            file_pointer.read(crt.FIXED_HEADER_BYTES)
        )
        file_pointer.seek(header.metadata_offset)
        metadata = json.loads(file_pointer.read(header.metadata_bytes).decode("utf-8"))
        if not isinstance(metadata, dict):
            raise RuntimeError("Combined metadata JSON was not a dictionary.")
        file_pointer.seek(header.slice_table_offset)
        slice_times: List[float] = []
        for _ in range(header.num_time_slices):
            entry = crt._unpack_slice_entry(  # pylint: disable=protected-access
                file_pointer.read(crt.SLICE_TABLE_ENTRY_BYTES)
            )
            slice_times.append(entry.simulation_time)
    return {
        "header": header,
        "metadata": metadata,
        "slice_times": slice_times,
        "file_size": file_size,
    }


def _combined_bin_matches(
    normalized: Dict[str, object], parsed_metadata: Dict[str, object]
) -> bool:
    combined_file = cast(Dict[str, object], normalized["combined_file"])
    grid = cast(Dict[str, object], normalized["grid"])
    payload = cast(Dict[str, object], normalized["payload"])
    time = cast(Dict[str, object], normalized["time"])
    spatial_lookup = cast(Dict[str, object], normalized["spatial_lookup"])
    run = cast(Dict[str, object], normalized["two_blackholes_run"])

    header = cast(crt.CombinedHeaderInfo, parsed_metadata["header"])
    metadata = cast(Dict[str, object], parsed_metadata["metadata"])
    slice_times = cast(List[float], parsed_metadata["slice_times"])
    file_size = cast(int, parsed_metadata["file_size"])

    if header.magic != combined_file["format_magic"]:
        return False
    if header.combined_format_version != combined_file["combined_format_version"]:
        return False
    if header.source_format_version != combined_file["source_format_version"]:
        return False
    if header.endian_tag != crt.ENDIAN_TAG:
        return False
    if header.serialized_real_bytes != combined_file["serialized_real_bytes"]:
        return False
    if header.num_grids != grid["num_grids"]:
        return False
    if header.payload_includes_ghost_zones != grid["payload_includes_ghost_zones"]:
        return False
    if list(header.Nxx) != grid["Nxx"]:
        return False
    if header.point_record_real_count != payload["point_record_real_count"]:
        return False
    if header.point_record_bytes != payload["point_record_bytes"]:
        return False
    if header.record_component_count != payload["record_component_count"]:
        return False
    if header.metric_component_count != payload["metric_component_count"]:
        return False
    if header.christoffel_component_count != payload["christoffel_component_count"]:
        return False
    if header.total_file_bytes != file_size:
        return False
    if metadata.get("format_name") != payload["format_name"]:
        return False
    if metadata.get("payload_layout") != payload["payload_layout"]:
        return False
    if metadata.get("loop_order") != payload["loop_order"]:
        return False
    if metadata.get("source_coord_system") != grid["CoordSystem"]:
        return False
    if metadata.get("target_basis") != grid["target_basis"]:
        return False
    if metadata.get("spatial_lookup_mode") != spatial_lookup["spatial_lookup_mode"]:
        return False

    axisymmetry = metadata.get("axisymmetry")
    if not isinstance(axisymmetry, dict):
        return False
    if axisymmetry.get("enabled") != spatial_lookup["axisymmetry_enabled"]:
        return False
    if axisymmetry.get("symmetry_axis") != spatial_lookup["axisymmetry_axis"]:
        return False
    if (
        axisymmetry.get("requires_rotation_for_arbitrary_phi")
        != spatial_lookup["requires_axisymmetry_rotation"]
    ):
        return False

    if len(slice_times) == 0:
        return False
    previous_time: Optional[float] = None
    for time_value in slice_times:
        if previous_time is not None and time_value <= previous_time:
            return False
        previous_time = time_value
    if not math.isclose(
        slice_times[0],
        cast(float, time["t_start"]),
        rel_tol=0.0,
        abs_tol=cast(float, time["absolute_tolerance"]),
    ):
        return False
    required_t_final = cast(float, time["t_final"])
    t_final_tolerance = cast(float, time["t_final_tolerance"])
    runtime_t_final = cast(float, run["runtime_t_final"])
    # The photon integrator needs coverage through required_t_final, but the
    # BBH driver intentionally runs about one output interval longer because
    # diagnostics are emitted before each forward step.
    if slice_times[-1] < required_t_final - t_final_tolerance:
        return False
    if slice_times[-1] > runtime_t_final + t_final_tolerance:
        return False
    return True


def _format_parfile_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        return repr(value)
    if isinstance(value, str):
        escaped = value.replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(value, tuple):
        return "{ " + ", ".join(_format_parfile_value(entry) for entry in value) + " }"
    raise ValueError(f"Unsupported parfile value type: {type(value).__name__}")


def _rewrite_parfile(parfile_path: Path, replacements: Dict[str, object]) -> None:
    parameter_line_pattern = re.compile(
        r"^([A-Za-z0-9_]+(?:\[\d+\])?)\s*=\s*(.*?)\s*(#.*)?$"
    )
    seen_keys: Dict[str, bool] = {key: False for key in replacements}
    updated_lines: List[str] = []
    with parfile_path.open("r", encoding="utf-8") as file_pointer:
        for line in file_pointer:
            match = parameter_line_pattern.match(line.rstrip("\n"))
            if match is None:
                updated_lines.append(line)
                continue
            parameter_name = match.group(1)
            trailing_comment = match.group(3) or ""
            if parameter_name not in replacements:
                updated_lines.append(line)
                continue
            replacement = _format_parfile_value(replacements[parameter_name])
            updated_lines.append(
                f"{parameter_name} = {replacement} {trailing_comment}".rstrip() + "\n"
            )
            seen_keys[parameter_name] = True

    missing_keys = [key for key, found in seen_keys.items() if not found]
    if missing_keys:
        raise RuntimeError(
            f"Could not rewrite parfile {parfile_path}; missing keys: {missing_keys}"
        )

    temp_fd, temp_name = tempfile.mkstemp(
        prefix=f"{parfile_path.name}.tmp.",
        dir=str(parfile_path.parent),
    )
    temp_path = Path(temp_name)
    try:
        with os.fdopen(temp_fd, "w", encoding="utf-8") as file_pointer:
            file_pointer.writelines(updated_lines)
            file_pointer.flush()
            os.fsync(file_pointer.fileno())
        os.replace(temp_path, parfile_path)
    finally:
        if temp_path.exists():
            temp_path.unlink()


def _build_parfile_replacements(normalized: Dict[str, object]) -> Dict[str, object]:
    run = cast(Dict[str, object], normalized["two_blackholes_run"])
    bh1_mass = cast(float, run["BH1_mass"])
    bh2_mass = cast(float, run["BH2_mass"])
    bh1_posn_z = cast(float, run["BH1_posn_z"])
    bh2_posn_z = cast(float, run["BH2_posn_z"])
    mass_total = bh1_mass + bh2_mass
    return {
        "eta": cast(float, run["GammaDriving_eta"]),
        "BH1_mass": bh1_mass,
        "BH2_mass": bh2_mass,
        "BH1_posn_z": bh1_posn_z,
        "BH2_posn_z": bh2_posn_z,
        "outer_bc_type": cast(str, run["outer_bc_type"]),
        "t_final": cast(float, run["runtime_t_final"]),
        "diagnostics_output_every": cast(float, run["diagnostics_output_every"]),
        "bah_initial_grid_z_center[3]": (bh1_posn_z, bh2_posn_z, 0.0),
        "bah_M_scale[3]": (bh1_mass, bh2_mass, mass_total),
        "bah_max_search_radius[3]": (
            0.6 * bh1_mass,
            0.6 * bh2_mass,
            1.9 * mass_total,
        ),
    }


def _run_subprocess(command: Sequence[str], cwd: Path, label: str) -> None:
    try:
        subprocess.run(list(command), cwd=str(cwd), check=True)
    except subprocess.CalledProcessError as error:
        raise RuntimeError(
            f"{label} failed with exit code {error.returncode}: {list(command)!r}"
        ) from error


def _combine_stage1_outputs(
    normalized: Dict[str, object], combined_bin_path: Path
) -> None:
    generation = cast(Dict[str, object], normalized["generation"])
    spatial_lookup = cast(Dict[str, object], normalized["spatial_lookup"])

    command: List[str] = [
        cast(str, generation["python_executable"]),
        cast(str, generation["combine_raytracing_time_slices_script"]),
        "--input-dir",
        cast(str, generation["stage1_raytracing_output_dir"]),
        "--pattern",
        cast(str, generation["stage1_raytracing_pattern"]),
        "--output",
        str(combined_bin_path),
        "--force",
        "--spatial-lookup-mode",
        cast(str, spatial_lookup["spatial_lookup_mode"]),
        "--include-coordinate-table",
        "--validate-coordinate-table",
    ]
    if cast(bool, spatial_lookup["axisymmetry_enabled"]):
        command.extend(
            [
                "--axisymmetry-enabled",
                "--axisymmetry-axis",
                cast(str, spatial_lookup["axisymmetry_axis"]),
            ]
        )
        command.append(
            "--requires-axisymmetry-rotation"
            if cast(bool, spatial_lookup["requires_axisymmetry_rotation"])
            else "--no-requires-axisymmetry-rotation"
        )
    else:
        command.append("--no-axisymmetry")
    _run_subprocess(command, REPO_ROOT, "combine_raytracing_time_slices.py")


def check_combined_bin_metadata(
    required_metadata: Dict[str, object],
    combined_bin_location: str,
) -> Tuple[bool, bool]:
    """
    Check whether a combined `.bin` exists and matches the required metadata.

    The helper validates only lightweight compatibility metadata: fixed-header
    fields, JSON metadata, and slice-table times. Time gaps are required to be
    strictly positive; they are not required to be uniform. The first slice
    must match ``t_start``. The last slice must cover the required photon-side
    ``t_final`` and must not overshoot the BBH runtime stop time by more than
    the configured tolerance.

    :param required_metadata: Metadata describing the required combined bin.
    :param combined_bin_location: Expected combined bin path.
    :return: Tuple ``(bin_exists, metadata_is_correct)``.
    """
    normalized = _normalize_required_metadata(required_metadata)
    combined_bin_path = Path(combined_bin_location).expanduser().resolve()
    if not combined_bin_path.is_file():
        return False, False
    try:
        parsed_metadata = _parse_combined_bin_metadata(combined_bin_path)
    except Exception:  # pylint: disable=broad-except
        return True, False
    if not _combined_bin_matches(normalized, parsed_metadata):
        return True, False
    if not _manifest_matches(_manifest_path(combined_bin_path), normalized):
        return True, False
    return True, True


def generate_required_combined_bin(
    required_metadata: Dict[str, object],
    combined_bin_location: str,
) -> str:
    """
    Generate the required combined raytracing bin.

    The helper drives ``two_blackholes_collide.py`` with raytracing outputs
    enabled, rewrites the generated parfile for runtime values, builds the
    project, runs the executable, combines the stage-1 files, and writes a
    small sidecar manifest next to the final combined `.bin`.

    :param required_metadata: Metadata describing the required combined bin.
    :param combined_bin_location: Final combined bin path.
    :return: Final combined bin path.
    """
    normalized = _normalize_required_metadata(required_metadata)
    generation = cast(Dict[str, object], normalized["generation"])
    run = cast(Dict[str, object], normalized["two_blackholes_run"])

    combined_bin_path = Path(combined_bin_location).expanduser().resolve()
    combined_bin_path.parent.mkdir(parents=True, exist_ok=True)

    codegen_command = [
        cast(str, generation["python_executable"]),
        cast(str, generation["two_blackholes_example_script"]),
        "--raytracing-outputs",
        "--floating_point_precision",
        cast(str, run["floating_point_precision"]),
    ]
    _run_subprocess(
        codegen_command, REPO_ROOT, "two_blackholes_collide.py code generation"
    )

    project_dir = Path(cast(str, generation["project_dir"]))
    executable_name = cast(str, generation["executable_name"])
    parfile_path = project_dir / f"{executable_name}.par"
    executable_path = project_dir / executable_name

    _rewrite_parfile(parfile_path, _build_parfile_replacements(normalized))
    _run_subprocess(
        cast(List[str], generation["make_command"]),
        project_dir,
        "two_blackholes_collide make",
    )
    _run_subprocess(
        [str(executable_path), str(parfile_path)],
        project_dir,
        "two_blackholes_collide runtime evolution",
    )
    _combine_stage1_outputs(normalized, combined_bin_path)
    _write_manifest(_manifest_path(combined_bin_path), normalized)
    return str(combined_bin_path)


def ensure_required_combined_bin(
    required_metadata: Dict[str, object],
    combined_bin_location: str,
) -> Tuple[str, str]:
    """
    Ensure that a compatible combined raytracing bin exists.

    :param required_metadata: Metadata describing the required combined bin.
    :param combined_bin_location: Final combined bin path.
    :return: Tuple ``(state_of_bin, combined_bin_location)``.
    :raises RuntimeError: If generation completes but the final combined bin
        still fails validation.
    """
    combined_bin_path = str(Path(combined_bin_location).expanduser().resolve())
    bin_exists, metadata_is_correct = check_combined_bin_metadata(
        required_metadata,
        combined_bin_path,
    )
    if not bin_exists:
        final_path = generate_required_combined_bin(
            required_metadata,
            combined_bin_path,
        )
        _, metadata_is_correct = check_combined_bin_metadata(
            required_metadata,
            final_path,
        )
        if not metadata_is_correct:
            raise RuntimeError("Generated combined .bin did not pass final validation.")
        return STATE_GENERATED_MISSING, final_path

    if metadata_is_correct:
        return STATE_ALREADY_CORRECT, combined_bin_path

    final_path = generate_required_combined_bin(
        required_metadata,
        combined_bin_path,
    )
    _, metadata_is_correct = check_combined_bin_metadata(
        required_metadata,
        final_path,
    )
    if not metadata_is_correct:
        raise RuntimeError("Regenerated combined .bin did not pass final validation.")
    return STATE_REBUILT_STALE, final_path


if __name__ == "__main__":
    import doctest

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

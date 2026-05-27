"""
Combine stage-1 raytracing time-slice files into one immutable container.

This module parses the binary files written by
``nrpy.infrastructures.BHaH.diagnostics.output_raytracing_data``, validates
their headers strictly, sorts the inputs by physical simulation time, and
writes one read-only combined container for downstream raytracing.

The combined file stores only copied stage-1 point-record payloads plus the
metadata needed for stage 3 to locate time slices and recover the payload-local
logical-grid indexing convention. It does not recompute metric data,
Christoffels, or coordinate transforms, and it does not build a true spatial
index.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import json
import math
import os
import struct
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Sequence, Set, Tuple, cast

STAGE1_MAGIC = "NRPYRTDATA4D"
STAGE1_HEADER_LABEL = "NRPy RT spacetime data"
STAGE1_FORMAT_NAME = "Cartesian g4DD+Gamma4UDD"
STAGE1_TARGET_BASIS = "Cartesian"
STAGE1_LOOP_ORDER = "i2maj_i0fast"
STAGE1_FORMAT_VERSION = 1
STAGE1_HEADER_SIZE = 2568

COMBINED_MAGIC = b"NRPYRTSTACK4D\0\0\0"
GEOMETRY_MAGIC = b"NRPYRTGEOMV1\0\0\0\0"

COMBINED_FORMAT_VERSION = 1
GEOMETRY_FORMAT_VERSION = 1
ENDIAN_TAG = 0x01020304

FIXED_HEADER_BYTES = 4096
DEFAULT_ALIGNMENT = 4096
SLICE_TABLE_ENTRY_BYTES = 128
GEOMETRY_HEADER_BYTES = 512

SERIALIZED_REAL_BYTES = 8
POINT_RECORD_REAL_COUNT = 53
POINT_RECORD_BYTES = 424
RECORD_COMPONENT_COUNT = 53
METRIC_COMPONENT_COUNT = 10
CHRISTOFFEL_COMPONENT_COUNT = 40

RECORD_LAYOUT_TIME_MAJOR_STAGE1_AOS = 1

SPATIAL_LOOKUP_NATIVE_LOGICAL_GRID = 1
SPATIAL_LOOKUP_COORDINATE_TABLE_ONLY = 2
SPATIAL_LOOKUP_REQUIRES_READER_SPATIAL_INDEX = 3

GEOMETRY_COORD_TABLE_NONE = 0
GEOMETRY_COORD_TABLE_CARTESIAN_POINT_COORDS_AOS_F64 = 1

SPATIAL_LOOKUP_MODE_NAMES = {
    SPATIAL_LOOKUP_NATIVE_LOGICAL_GRID: "native_logical_grid",
    SPATIAL_LOOKUP_COORDINATE_TABLE_ONLY: "coordinate_table_only",
    SPATIAL_LOOKUP_REQUIRES_READER_SPATIAL_INDEX: "requires_reader_spatial_index",
}
SPATIAL_LOOKUP_MODE_IDS = {
    "native_logical_grid": SPATIAL_LOOKUP_NATIVE_LOGICAL_GRID,
    "coordinate_table_only": SPATIAL_LOOKUP_COORDINATE_TABLE_ONLY,
    "requires_reader_spatial_index": SPATIAL_LOOKUP_REQUIRES_READER_SPATIAL_INDEX,
}

GRID_POINT_POSITION_CONVENTION_ID = 1
HEADER_FLAGS_V1 = 0
GEOMETRY_FLAGS_V1 = 0

RECORD_COMPONENT_NAMES = (
    "x",
    "y",
    "z",
    "g4DD00",
    "g4DD01",
    "g4DD02",
    "g4DD03",
    "g4DD11",
    "g4DD12",
    "g4DD13",
    "g4DD22",
    "g4DD23",
    "g4DD33",
    "Gamma4UDD000",
    "Gamma4UDD001",
    "Gamma4UDD002",
    "Gamma4UDD003",
    "Gamma4UDD011",
    "Gamma4UDD012",
    "Gamma4UDD013",
    "Gamma4UDD022",
    "Gamma4UDD023",
    "Gamma4UDD033",
    "Gamma4UDD100",
    "Gamma4UDD101",
    "Gamma4UDD102",
    "Gamma4UDD103",
    "Gamma4UDD111",
    "Gamma4UDD112",
    "Gamma4UDD113",
    "Gamma4UDD122",
    "Gamma4UDD123",
    "Gamma4UDD133",
    "Gamma4UDD200",
    "Gamma4UDD201",
    "Gamma4UDD202",
    "Gamma4UDD203",
    "Gamma4UDD211",
    "Gamma4UDD212",
    "Gamma4UDD213",
    "Gamma4UDD222",
    "Gamma4UDD223",
    "Gamma4UDD233",
    "Gamma4UDD300",
    "Gamma4UDD301",
    "Gamma4UDD302",
    "Gamma4UDD303",
    "Gamma4UDD311",
    "Gamma4UDD312",
    "Gamma4UDD313",
    "Gamma4UDD322",
    "Gamma4UDD323",
    "Gamma4UDD333",
)
METRIC_COMPONENT_NAMES = (
    "g4DD00",
    "g4DD01",
    "g4DD02",
    "g4DD03",
    "g4DD11",
    "g4DD12",
    "g4DD13",
    "g4DD22",
    "g4DD23",
    "g4DD33",
)
CHRISTOFFEL_COMPONENT_NAMES = RECORD_COMPONENT_NAMES[13:]


@dataclass(frozen=True)
class Stage1Info:
    """
    Parsed stage-1 file metadata and offsets.

    This record mirrors the on-disk header plus the separately stored physical
    simulation time so later validation and combined-file layout code can work
    with one immutable schema object.
    """

    path: Path
    magic: str
    format_version: int
    header_size: int
    output_index: int
    num_grids: int
    serialized_real_bytes: int
    file_is_little_endian: int
    time_variable_is_f64: int
    record_component_count: int
    metric_component_count: int
    christoffel_component_count: int
    point_record_real_count: int
    point_record_bytes: int
    point_record_count: int
    payload_includes_ghost_zones: int
    Nxx: Tuple[int, int, int]
    Nxx_plus_2NGHOSTS: Tuple[int, int, int]
    NGHOSTS: int
    payload_i_count: Tuple[int, int, int]
    payload_i_start: Tuple[int, int, int]
    payload_i_end: Tuple[int, int, int]
    dxx: Tuple[float, float, float]
    invdxx: Tuple[float, float, float]
    xxmin: Tuple[float, float, float]
    xxmax: Tuple[float, float, float]
    cart_origin: Tuple[float, float, float]
    simulation_time_offset: int
    point_records_offset: int
    point_records_bytes: int
    total_file_bytes: int
    actual_file_size: int
    header_label: str
    format_name: str
    target_basis: str
    source_coord_system: str
    loop_order: str
    cf_convention: str
    record_component_names: Tuple[str, ...]
    metric_component_names: Tuple[str, ...]
    christoffel_component_names: Tuple[str, ...]
    simulation_time: float


@dataclass(frozen=True)
class Layout:
    """
    Concrete byte layout for the combined container.

    The combiner computes this repeatedly until metadata-dependent offsets
    stabilize, then uses the final values verbatim when writing the file.
    """

    alignment_bytes: int
    metadata_offset: int
    metadata_bytes: int
    slice_table_offset: int
    slice_table_bytes: int
    geometry_block_offset: int
    geometry_block_bytes: int
    first_payload_offset: int
    payload_bytes_per_slice: int
    payload_stride_bytes: int
    payload_bytes_total: int
    total_file_bytes: int


@dataclass(frozen=True)
class SliceEntry:
    """
    One fixed-width slice-table entry for the combined container.

    Each entry records the sorted simulation time together with the absolute
    payload offset and provenance offsets from the original stage-1 file.
    """

    output_index: int
    slice_flags: int
    simulation_time: float
    payload_offset: int
    payload_bytes: int
    point_record_count: int
    point_record_bytes: int
    source_file_size: int
    source_header_size: int
    source_simulation_time_offset: int
    source_point_records_offset: int
    source_point_records_bytes: int
    source_total_file_bytes: int


@dataclass(frozen=True)
class AxisymmetryMetadata:
    """
    Reader-facing axisymmetry contract stored in combined metadata.

    Stage 2 records only conventions here. Stage 3 remains responsible for any
    interpolation and arbitrary-phi tensor rotation logic.
    """

    enabled: bool
    symmetry_axis: str
    phi_sample_count: int
    stored_phi_samples: Tuple[float, ...]
    requires_rotation_for_arbitrary_phi: bool


@dataclass(frozen=True)
class GeometryBlockData:
    """
    GeometryBlockV1 bytes plus derived summary values.

    The summary fields let metadata generation and validation reuse the same
    geometry facts without reparsing the raw block bytes.
    """

    block_bytes: bytes
    coordinate_table_present: bool
    coordinate_table_bytes: int
    coordinate_table_kind: int
    coordinate_table_point_count: int
    cartesian_bbox_min: Optional[Tuple[float, float, float]]
    cartesian_bbox_max: Optional[Tuple[float, float, float]]


@dataclass(frozen=True)
class CombinedHeaderInfo:
    """Parsed FixedHeaderV1 fields from a combined container."""

    magic: str
    combined_format_version: int
    source_format_version: int
    endian_tag: int
    header_flags: int
    fixed_header_bytes: int
    alignment_bytes: int
    metadata_offset: int
    metadata_bytes: int
    slice_table_offset: int
    slice_table_entry_bytes: int
    num_time_slices: int
    geometry_block_offset: int
    geometry_block_bytes: int
    first_payload_offset: int
    payload_bytes_per_slice: int
    payload_stride_bytes: int
    payload_bytes_total: int
    total_file_bytes: int
    point_record_count: int
    serialized_real_bytes: int
    point_record_real_count: int
    point_record_bytes: int
    record_component_count: int
    metric_component_count: int
    christoffel_component_count: int
    file_is_little_endian: int
    time_variable_is_f64: int
    payload_includes_ghost_zones: int
    num_grids: int
    NGHOSTS: int
    record_layout_id: int
    spatial_lookup_mode_id: int
    axisymmetry_flag: int
    phi_sample_count: int
    geometry_format_version: int
    reserved_u32: int
    Nxx: Tuple[int, int, int]
    Nxx_plus_2NGHOSTS: Tuple[int, int, int]
    payload_i_count: Tuple[int, int, int]
    payload_i_start: Tuple[int, int, int]
    payload_i_end: Tuple[int, int, int]
    dxx: Tuple[float, float, float]
    invdxx: Tuple[float, float, float]
    xxmin: Tuple[float, float, float]
    xxmax: Tuple[float, float, float]
    cart_origin: Tuple[float, float, float]


@dataclass(frozen=True)
class GeometryHeaderInfo:
    """Parsed GeometryHeaderV1 fields from a combined container."""

    magic: str
    geometry_format_version: int
    geometry_flags: int
    spatial_lookup_mode_id: int
    coordinate_table_kind: int
    geometry_header_bytes: int
    coordinate_table_offset_from_block: int
    coordinate_table_point_count: int
    coordinate_table_components: int
    coordinate_table_real_bytes: int
    coordinate_table_bytes: int
    grid_point_position_convention_id: int
    reserved_u32: int
    reserved_u64: Tuple[int, int, int, int, int, int, int, int]


@dataclass(frozen=True)
class SliceEntryInfo:
    """Parsed fixed-width slice-table entry from a combined container."""

    output_index: int
    slice_flags: int
    simulation_time: float
    payload_offset: int
    payload_bytes: int
    point_record_count: int
    point_record_bytes: int
    reserved_u32: int
    source_file_size: int
    source_header_size: int
    source_simulation_time_offset: int
    source_point_records_offset: int
    source_point_records_bytes: int
    source_total_file_bytes: int
    reserved_u64: Tuple[int, int, int, int]


def align_up(n: int, alignment: int) -> int:
    """
    Round n up to the next alignment boundary.

    :param n: Integer byte count to round up.
    :param alignment: Power-of-two byte alignment.
    :return: n rounded up to a multiple of alignment.

    Doctests:
    >>> align_up(13, 8)
    16
    """
    return ((n + alignment - 1) // alignment) * alignment


def _validate_power_of_two(value: int, label: str) -> None:
    """
    Validate that an integer is a positive power of two.

    :param value: Value to validate.
    :param label: Human-readable label for error reporting.
    :raises ValueError: If the value is not a positive power of two.
    """
    if value <= 0 or (value & (value - 1)) != 0:
        raise ValueError(f"{label} must be a positive power of two, got {value}.")


def _read_exact(fp: BinaryIO, size: int, label: str) -> bytes:
    """
    Read an exact byte count from a file object.

    :param fp: Open binary file object.
    :param size: Number of bytes to read.
    :param label: Human-readable label for error reporting.
    :return: Bytes read from the file.
    :raises RuntimeError: If EOF is encountered early.
    """
    data = fp.read(size)
    if len(data) != size:
        raise RuntimeError(f"Unexpected EOF while reading {label}.")
    return data


def _read_fixed_string(fp: BinaryIO, size: int, label: str) -> str:
    """
    Read a null-padded fixed-width ASCII string.

    :param fp: Open binary file object.
    :param size: On-disk field width.
    :param label: Human-readable label for error reporting.
    :return: Decoded string with trailing null bytes removed.
    """
    raw = _read_exact(fp, size, label)
    return raw.split(b"\0", 1)[0].decode("ascii")


def _read_u32(fp: BinaryIO, label: str) -> int:
    """
    Read one little-endian uint32.

    :param fp: Open binary file object.
    :param label: Human-readable label for error reporting.
    :return: Parsed integer value.
    """
    return cast(int, struct.unpack("<I", _read_exact(fp, 4, label))[0])


def _read_u64(fp: BinaryIO, label: str) -> int:
    """
    Read one little-endian uint64.

    :param fp: Open binary file object.
    :param label: Human-readable label for error reporting.
    :return: Parsed integer value.
    """
    return cast(int, struct.unpack("<Q", _read_exact(fp, 8, label))[0])


def _read_f64(fp: BinaryIO, label: str) -> float:
    """
    Read one little-endian float64.

    :param fp: Open binary file object.
    :param label: Human-readable label for error reporting.
    :return: Parsed floating-point value.
    """
    return cast(float, struct.unpack("<d", _read_exact(fp, 8, label))[0])


def _pack_fixed_string(text: str, size: int) -> bytes:
    r"""
    Pack an ASCII string into a null-padded fixed-width field.

    :param text: ASCII string to encode.
    :param size: Exact output field width.
    :return: Null-padded bytes of length size.
    :raises ValueError: If the string does not fit.

    Doctests:
    >>> _pack_fixed_string("nrpy", 8)
    b'nrpy\x00\x00\x00\x00'
    """
    encoded = text.encode("ascii")
    if len(encoded) >= size:
        raise ValueError(f"String '{text}' does not fit in a {size}-byte field.")
    return encoded + b"\0" * (size - len(encoded))


def _require_u32(value: int, label: str) -> int:
    """
    Validate that an integer fits in uint32.

    :param value: Integer to validate.
    :param label: Human-readable label for error reporting.
    :return: The validated integer.
    :raises RuntimeError: If the integer does not fit in uint32.
    """
    if not 0 <= value <= 0xFFFFFFFF:
        raise RuntimeError(f"{label}={value} does not fit in uint32.")
    return value


def _require_u64(value: int, label: str) -> int:
    """
    Validate that an integer fits in uint64.

    :param value: Integer to validate.
    :param label: Human-readable label for error reporting.
    :return: The validated integer.
    :raises RuntimeError: If the integer does not fit in uint64.
    """
    if not 0 <= value <= 0xFFFFFFFFFFFFFFFF:
        raise RuntimeError(f"{label}={value} does not fit in uint64.")
    return value


def _require_i64(value: int, label: str) -> int:
    """
    Validate that an integer fits in int64.

    :param value: Integer to validate.
    :param label: Human-readable label for error reporting.
    :return: The validated integer.
    :raises RuntimeError: If the integer does not fit in int64.
    """
    if not -(1 << 63) <= value < (1 << 63):
        raise RuntimeError(f"{label}={value} does not fit in int64.")
    return value


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    :return: Parsed CLI namespace.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Combine stage-1 raytracing time-slice files into one immutable "
            "stacked container."
        )
    )
    parser.add_argument("input_files", nargs="*", help="Explicit input stage-1 files.")
    parser.add_argument(
        "--output",
        default="combined_raytracing_data.bin",
        help="Output combined file path.",
    )
    parser.add_argument(
        "--input-dir",
        default=None,
        help="Directory to search when using --pattern.",
    )
    parser.add_argument(
        "--pattern",
        default="raytracing_data_t*.bin",
        help="Input glob pattern used inside --input-dir.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Allow overwriting an existing output file.",
    )
    parser.add_argument(
        "--alignment",
        type=int,
        default=DEFAULT_ALIGNMENT,
        help="Alignment in bytes for major blocks and payload slices.",
    )
    parser.add_argument(
        "--include-coordinate-table",
        dest="include_coordinate_table",
        action="store_true",
        default=True,
        help="Include cartesian_point_coords[point][3] in GeometryBlockV1.",
    )
    parser.add_argument(
        "--no-coordinate-table",
        dest="include_coordinate_table",
        action="store_false",
        help="Do not include a GeometryBlock coordinate table.",
    )
    parser.add_argument(
        "--validate-coordinate-table",
        dest="validate_coordinate_table",
        action="store_true",
        default=True,
        help="Verify later slices reuse the first-slice x,y,z coordinates.",
    )
    parser.add_argument(
        "--no-validate-coordinate-table",
        dest="validate_coordinate_table",
        action="store_false",
        help="Skip coordinate table validation.",
    )
    parser.add_argument(
        "--spatial-lookup-mode",
        choices=sorted(SPATIAL_LOOKUP_MODE_IDS.keys()),
        help="Declared reader lookup mode for the combined container.",
    )
    parser.add_argument(
        "--native-inverse-map-name",
        default=None,
        help=(
            "Required inverse-map identifier when "
            "--spatial-lookup-mode=native_logical_grid."
        ),
    )
    parser.add_argument(
        "--axisymmetry-enabled",
        dest="axisymmetry_enabled",
        action="store_true",
        default=True,
        help=(
            "Record axisymmetry metadata in the combined file. Default matches "
            "the current axisymmetric two-phi-sample pipeline."
        ),
    )
    parser.add_argument(
        "--no-axisymmetry",
        dest="axisymmetry_enabled",
        action="store_false",
        help="Do not record axisymmetry metadata.",
    )
    parser.add_argument(
        "--axisymmetry-axis",
        choices=("x", "y", "z"),
        default="z",
        help="Axisymmetry axis name.",
    )
    parser.add_argument(
        "--phi-samples",
        default="0,3.141592653589793",
        help="Comma-separated stored phi sample values.",
    )
    parser.add_argument(
        "--requires-axisymmetry-rotation",
        dest="requires_axisymmetry_rotation",
        action="store_true",
        default=True,
        help="Record that arbitrary-phi queries require reader-side rotation.",
    )
    parser.add_argument(
        "--no-requires-axisymmetry-rotation",
        dest="requires_axisymmetry_rotation",
        action="store_false",
        help="Record that arbitrary-phi rotation is not required.",
    )
    parser.add_argument(
        "--inspect",
        help="Inspect an existing combined file instead of generating one.",
    )
    return parser.parse_args()


def discover_input_paths(args: argparse.Namespace) -> List[Path]:
    """
    Discover explicit and globbed input files.

    :param args: Parsed CLI namespace.
    :return: Sorted unique input file paths.
    :raises RuntimeError: If an explicit input path does not exist or is temporary.
    """
    explicit_paths = [Path(path_str) for path_str in args.input_files]
    resolved: Dict[Path, Path] = {}
    for path in explicit_paths:
        if not path.is_file():
            raise RuntimeError(f"Explicit input file not found: {path}")
        if ".tmp." in path.name:
            raise RuntimeError(f"Explicit input file appears temporary: {path}")
        resolved[path.resolve()] = path.resolve()
    if args.input_dir is not None or not explicit_paths:
        input_dir = Path(args.input_dir) if args.input_dir is not None else Path(".")
        for path in input_dir.glob(args.pattern):
            if ".tmp." in path.name:
                continue
            if path.is_file():
                resolved[path.resolve()] = path.resolve()
    return sorted(resolved.values())


def parse_stage1_file(path: Path) -> Stage1Info:
    """
    Parse one stage-1 raytracing output file exactly.

    :param path: Stage-1 binary file path.
    :return: Parsed header and simulation-time metadata.
    :raises RuntimeError: If the file is truncated or internally inconsistent.
    """
    # Step 1: Read the fixed-width header exactly in the stage-1 writer order.
    actual_file_size = path.stat().st_size
    with path.open("rb") as fp:
        magic = _read_fixed_string(fp, 16, "stage-1 magic")

        format_version = _read_u32(fp, "format_version")
        header_size = _read_u32(fp, "header_size")
        output_index = _read_u32(fp, "output_index")
        num_grids = _read_u32(fp, "num_grids")

        serialized_real_bytes = _read_u32(fp, "serialized_real_bytes")
        record_component_count = _read_u32(fp, "record_component_count")
        metric_component_count = _read_u32(fp, "metric_component_count")
        christoffel_component_count = _read_u32(fp, "christoffel_component_count")
        point_record_real_count = _read_u32(fp, "point_record_real_count")
        point_record_bytes = _read_u32(fp, "point_record_bytes")
        payload_includes_ghost_zones = _read_u32(fp, "payload_includes_ghost_zones")
        file_is_little_endian = _read_u32(fp, "file_is_little_endian")
        time_variable_is_f64 = _read_u32(fp, "time_variable_is_f64")
        _ = _read_u32(fp, "reserved_u32")

        # Step 1.a: Reject unexpected schema counts before reading any
        #           variable-length string tables from the file.
        if record_component_count != RECORD_COMPONENT_COUNT:
            raise RuntimeError(
                f"{path}: record_component_count must be {RECORD_COMPONENT_COUNT}, "
                f"found {record_component_count}."
            )
        if metric_component_count != METRIC_COMPONENT_COUNT:
            raise RuntimeError(
                f"{path}: metric_component_count must be {METRIC_COMPONENT_COUNT}, "
                f"found {metric_component_count}."
            )
        if christoffel_component_count != CHRISTOFFEL_COMPONENT_COUNT:
            raise RuntimeError(
                f"{path}: christoffel_component_count must be "
                f"{CHRISTOFFEL_COMPONENT_COUNT}, found "
                f"{christoffel_component_count}."
            )

        Nxx = (
            _read_u32(fp, "Nxx0"),
            _read_u32(fp, "Nxx1"),
            _read_u32(fp, "Nxx2"),
        )
        Nxx_plus_2NGHOSTS = (
            _read_u32(fp, "Nxx_plus_2NGHOSTS0"),
            _read_u32(fp, "Nxx_plus_2NGHOSTS1"),
            _read_u32(fp, "Nxx_plus_2NGHOSTS2"),
        )

        point_record_count = _read_u64(fp, "point_record_count")
        simulation_time_offset = _read_u64(fp, "simulation_time_offset")
        point_records_offset = _read_u64(fp, "point_records_offset")
        point_records_bytes = _read_u64(fp, "point_records_bytes")
        total_file_bytes = _read_u64(fp, "total_file_bytes")
        NGHOSTS = _read_u64(fp, "NGHOSTS")

        payload_i_count = (
            _read_u64(fp, "payload_i0_count"),
            _read_u64(fp, "payload_i1_count"),
            _read_u64(fp, "payload_i2_count"),
        )
        payload_i_start = (
            _read_u64(fp, "payload_i0_start"),
            _read_u64(fp, "payload_i1_start"),
            _read_u64(fp, "payload_i2_start"),
        )
        payload_i_end = (
            _read_u64(fp, "payload_i0_end"),
            _read_u64(fp, "payload_i1_end"),
            _read_u64(fp, "payload_i2_end"),
        )

        dxx = (
            _read_f64(fp, "dxx0"),
            _read_f64(fp, "dxx1"),
            _read_f64(fp, "dxx2"),
        )
        invdxx = (
            _read_f64(fp, "invdxx0"),
            _read_f64(fp, "invdxx1"),
            _read_f64(fp, "invdxx2"),
        )
        xxmin = (
            _read_f64(fp, "xxmin0"),
            _read_f64(fp, "xxmin1"),
            _read_f64(fp, "xxmin2"),
        )
        xxmax = (
            _read_f64(fp, "xxmax0"),
            _read_f64(fp, "xxmax1"),
            _read_f64(fp, "xxmax2"),
        )
        cart_origin = (
            _read_f64(fp, "cart_originx"),
            _read_f64(fp, "cart_originy"),
            _read_f64(fp, "cart_originz"),
        )

        header_label = _read_fixed_string(fp, 32, "header_label")
        format_name = _read_fixed_string(fp, 32, "format_name")
        target_basis = _read_fixed_string(fp, 16, "target_basis")
        source_coord_system = _read_fixed_string(fp, 32, "source_coord_system")
        loop_order = _read_fixed_string(fp, 16, "loop_order")
        cf_convention = _read_fixed_string(fp, 32, "cf_convention")

        record_component_names = tuple(
            _read_fixed_string(fp, 24, f"record_component_name[{i}]")
            for i in range(record_component_count)
        )
        metric_component_names = tuple(
            _read_fixed_string(fp, 16, f"metric_component_name[{i}]")
            for i in range(metric_component_count)
        )
        christoffel_component_names = tuple(
            _read_fixed_string(fp, 16, f"christoffel_component_name[{i}]")
            for i in range(christoffel_component_count)
        )

        # Step 2: Confirm the parsed header consumed the exact advertised size.
        if fp.tell() != header_size:
            raise RuntimeError(
                f"{path}: parsed stage-1 header length {fp.tell()} did not match "
                f"header_size={header_size}."
            )

        # Step 3: Recover the separately stored physical simulation time.
        fp.seek(simulation_time_offset)
        simulation_time = _read_f64(fp, "simulation_time")

    return Stage1Info(
        path=path,
        magic=magic,
        format_version=format_version,
        header_size=header_size,
        output_index=output_index,
        num_grids=num_grids,
        serialized_real_bytes=serialized_real_bytes,
        file_is_little_endian=file_is_little_endian,
        time_variable_is_f64=time_variable_is_f64,
        record_component_count=record_component_count,
        metric_component_count=metric_component_count,
        christoffel_component_count=christoffel_component_count,
        point_record_real_count=point_record_real_count,
        point_record_bytes=point_record_bytes,
        point_record_count=point_record_count,
        payload_includes_ghost_zones=payload_includes_ghost_zones,
        Nxx=Nxx,
        Nxx_plus_2NGHOSTS=Nxx_plus_2NGHOSTS,
        NGHOSTS=NGHOSTS,
        payload_i_count=payload_i_count,
        payload_i_start=payload_i_start,
        payload_i_end=payload_i_end,
        dxx=dxx,
        invdxx=invdxx,
        xxmin=xxmin,
        xxmax=xxmax,
        cart_origin=cart_origin,
        simulation_time_offset=simulation_time_offset,
        point_records_offset=point_records_offset,
        point_records_bytes=point_records_bytes,
        total_file_bytes=total_file_bytes,
        actual_file_size=actual_file_size,
        header_label=header_label,
        format_name=format_name,
        target_basis=target_basis,
        source_coord_system=source_coord_system,
        loop_order=loop_order,
        cf_convention=cf_convention,
        record_component_names=record_component_names,
        metric_component_names=metric_component_names,
        christoffel_component_names=christoffel_component_names,
        simulation_time=simulation_time,
    )


def validate_stage1_file_internal(info: Stage1Info) -> None:
    """
    Validate one parsed stage-1 file against the documented schema.

    :param info: Parsed stage-1 file metadata.
    :raises RuntimeError: If the file fails validation.
    """
    if info.magic != STAGE1_MAGIC:
        raise RuntimeError(f"{info.path}: unexpected magic '{info.magic}'.")
    if info.format_version != STAGE1_FORMAT_VERSION:
        raise RuntimeError(
            f"{info.path}: expected format_version={STAGE1_FORMAT_VERSION}, "
            f"found {info.format_version}."
        )
    if info.header_size != STAGE1_HEADER_SIZE:
        raise RuntimeError(
            f"{info.path}: expected header_size={STAGE1_HEADER_SIZE}, "
            f"found {info.header_size}."
        )
    if info.header_label != STAGE1_HEADER_LABEL:
        raise RuntimeError(
            f"{info.path}: expected header_label='{STAGE1_HEADER_LABEL}', "
            f"found '{info.header_label}'."
        )
    if info.format_name != STAGE1_FORMAT_NAME:
        raise RuntimeError(
            f"{info.path}: expected format_name='{STAGE1_FORMAT_NAME}', "
            f"found '{info.format_name}'."
        )
    if info.target_basis != STAGE1_TARGET_BASIS:
        raise RuntimeError(
            f"{info.path}: expected target_basis='{STAGE1_TARGET_BASIS}', "
            f"found '{info.target_basis}'."
        )
    if info.loop_order != STAGE1_LOOP_ORDER:
        raise RuntimeError(
            f"{info.path}: expected loop_order='{STAGE1_LOOP_ORDER}', "
            f"found '{info.loop_order}'."
        )
    if info.num_grids != 1:
        raise RuntimeError(f"{info.path}: num_grids must be 1, found {info.num_grids}.")
    if info.serialized_real_bytes != SERIALIZED_REAL_BYTES:
        raise RuntimeError(
            f"{info.path}: serialized_real_bytes must be {SERIALIZED_REAL_BYTES}, "
            f"found {info.serialized_real_bytes}."
        )
    if info.file_is_little_endian != 1:
        raise RuntimeError(f"{info.path}: expected little-endian payload.")
    if info.time_variable_is_f64 != 1:
        raise RuntimeError(f"{info.path}: expected time_variable_is_f64=1.")
    if info.record_component_count != RECORD_COMPONENT_COUNT:
        raise RuntimeError(
            f"{info.path}: record_component_count must be {RECORD_COMPONENT_COUNT}, "
            f"found {info.record_component_count}."
        )
    if info.metric_component_count != METRIC_COMPONENT_COUNT:
        raise RuntimeError(
            f"{info.path}: metric_component_count must be {METRIC_COMPONENT_COUNT}, "
            f"found {info.metric_component_count}."
        )
    if info.christoffel_component_count != CHRISTOFFEL_COMPONENT_COUNT:
        raise RuntimeError(
            f"{info.path}: christoffel_component_count must be "
            f"{CHRISTOFFEL_COMPONENT_COUNT}, found "
            f"{info.christoffel_component_count}."
        )
    if info.point_record_real_count != POINT_RECORD_REAL_COUNT:
        raise RuntimeError(
            f"{info.path}: point_record_real_count must be {POINT_RECORD_REAL_COUNT}, "
            f"found {info.point_record_real_count}."
        )
    if info.point_record_bytes != POINT_RECORD_BYTES:
        raise RuntimeError(
            f"{info.path}: point_record_bytes must be {POINT_RECORD_BYTES}, "
            f"found {info.point_record_bytes}."
        )
    if info.payload_includes_ghost_zones != 0:
        raise RuntimeError(f"{info.path}: payload_includes_ghost_zones must be 0.")
    if info.record_component_names != RECORD_COMPONENT_NAMES:
        raise RuntimeError(
            f"{info.path}: record component names do not match v1 schema."
        )
    if info.metric_component_names != METRIC_COMPONENT_NAMES:
        raise RuntimeError(
            f"{info.path}: metric component names do not match v1 schema."
        )
    if info.christoffel_component_names != CHRISTOFFEL_COMPONENT_NAMES:
        raise RuntimeError(
            f"{info.path}: Christoffel component names do not match v1 schema."
        )
    if info.output_index < 0:
        raise RuntimeError(f"{info.path}: output_index must be nonnegative.")
    if info.point_record_count <= 0:
        raise RuntimeError(f"{info.path}: point_record_count must be positive.")
    if any(count <= 0 for count in info.payload_i_count):
        raise RuntimeError(f"{info.path}: payload_i_count entries must be positive.")
    expected_point_record_count = (
        info.payload_i_count[0] * info.payload_i_count[1] * info.payload_i_count[2]
    )
    if info.point_record_count != expected_point_record_count:
        raise RuntimeError(
            f"{info.path}: point_record_count={info.point_record_count} did not equal "
            f"payload_i_count product={expected_point_record_count}."
        )
    expected_point_records_bytes = info.point_record_count * info.point_record_bytes
    if info.point_records_bytes != expected_point_records_bytes:
        raise RuntimeError(
            f"{info.path}: point_records_bytes={info.point_records_bytes} did not "
            f"equal point_record_count*point_record_bytes={expected_point_records_bytes}."
        )
    if info.simulation_time_offset != info.header_size:
        raise RuntimeError(
            f"{info.path}: simulation_time_offset={info.simulation_time_offset} did "
            f"not equal header_size={info.header_size}."
        )
    if info.point_records_offset != info.header_size + 8:
        raise RuntimeError(
            f"{info.path}: point_records_offset={info.point_records_offset} did not "
            f"equal header_size + 8 = {info.header_size + 8}."
        )
    if info.total_file_bytes != info.point_records_offset + info.point_records_bytes:
        raise RuntimeError(
            f"{info.path}: total_file_bytes={info.total_file_bytes} did not equal "
            "point_records_offset + point_records_bytes."
        )
    if info.actual_file_size != info.total_file_bytes:
        raise RuntimeError(
            f"{info.path}: actual size {info.actual_file_size} did not match "
            f"total_file_bytes={info.total_file_bytes}."
        )
    # Step 4: Enforce the stage-1 v1 interior-only payload convention. The
    #         combined format records payload_i_start/count/end for reader
    #         clarity, but this v1 combiner intentionally rejects future
    #         stage-1 variants with ghost zones, partial domains, or shifted
    #         payload windows.
    if info.payload_i_count != info.Nxx:
        raise RuntimeError(
            f"{info.path}: payload_i_count={info.payload_i_count} did not equal "
            f"Nxx={info.Nxx}."
        )
    expected_nxx_plus_2nghosts = tuple(count + 2 * info.NGHOSTS for count in info.Nxx)
    if info.Nxx_plus_2NGHOSTS != expected_nxx_plus_2nghosts:
        raise RuntimeError(
            f"{info.path}: Nxx_plus_2NGHOSTS={info.Nxx_plus_2NGHOSTS} did not equal "
            f"Nxx + 2 * NGHOSTS = {expected_nxx_plus_2nghosts}."
        )
    expected_payload_i_start = (info.NGHOSTS, info.NGHOSTS, info.NGHOSTS)
    if info.payload_i_start != expected_payload_i_start:
        raise RuntimeError(
            f"{info.path}: payload_i_start={info.payload_i_start} did not equal "
            f"(NGHOSTS, NGHOSTS, NGHOSTS)={expected_payload_i_start}."
        )
    expected_payload_i_end = tuple(
        start + count
        for start, count in zip(info.payload_i_start, info.payload_i_count)
    )
    if info.payload_i_end != expected_payload_i_end:
        raise RuntimeError(
            f"{info.path}: payload_i_end={info.payload_i_end} did not equal "
            f"payload_i_start + payload_i_count = {expected_payload_i_end}."
        )
    for field_name in ("dxx", "invdxx", "xxmin", "xxmax", "cart_origin"):
        values = getattr(info, field_name)
        if not all(math.isfinite(value) for value in values):
            raise RuntimeError(f"{info.path}: {field_name} entries must be finite.")
    if not all(value > 0.0 for value in info.dxx):
        raise RuntimeError(f"{info.path}: dxx entries must be positive.")
    if not all(value > 0.0 for value in info.invdxx):
        raise RuntimeError(f"{info.path}: invdxx entries must be positive.")
    if not math.isfinite(info.simulation_time):
        raise RuntimeError(f"{info.path}: simulation_time must be finite.")


def validate_stage1_compatible(base: Stage1Info, other: Stage1Info) -> None:
    """
    Validate that two stage-1 files are container-compatible.

    :param base: Reference stage-1 file metadata.
    :param other: Candidate stage-1 file metadata.
    :raises RuntimeError: If any required field differs.
    """
    comparable_fields = (
        "magic",
        "format_version",
        "num_grids",
        "serialized_real_bytes",
        "file_is_little_endian",
        "time_variable_is_f64",
        "record_component_count",
        "metric_component_count",
        "christoffel_component_count",
        "point_record_real_count",
        "point_record_bytes",
        "point_record_count",
        "payload_includes_ghost_zones",
        "Nxx",
        "Nxx_plus_2NGHOSTS",
        "NGHOSTS",
        "payload_i_count",
        "payload_i_start",
        "payload_i_end",
        "dxx",
        "invdxx",
        "xxmin",
        "xxmax",
        "cart_origin",
        "header_label",
        "format_name",
        "target_basis",
        "source_coord_system",
        "loop_order",
        "cf_convention",
        "record_component_names",
        "metric_component_names",
        "christoffel_component_names",
    )
    for field_name in comparable_fields:
        if getattr(base, field_name) != getattr(other, field_name):
            raise RuntimeError(
                f"{other.path}: field '{field_name}' differed from "
                f"{base.path.name}."
            )


def validate_sorted_slices(infos: Sequence[Stage1Info]) -> None:
    """
    Validate uniqueness and strict time ordering after sorting.

    :param infos: Stage-1 infos sorted by simulation_time.
    :raises RuntimeError: If times or output indices are duplicated.
    """
    seen_indices: Set[int] = set()
    seen_times: Set[float] = set()
    previous_time: Optional[float] = None
    for info in infos:
        if info.output_index in seen_indices:
            raise RuntimeError(
                f"Duplicate output_index {info.output_index} found in {info.path}."
            )
        if info.simulation_time in seen_times:
            raise RuntimeError(
                f"Duplicate simulation_time {info.simulation_time:.17g} found in "
                f"{info.path}."
            )
        if previous_time is not None and info.simulation_time <= previous_time:
            raise RuntimeError("Simulation times must be strictly increasing.")
        seen_indices.add(info.output_index)
        seen_times.add(info.simulation_time)
        previous_time = info.simulation_time


def determine_axisymmetry_metadata(args: argparse.Namespace) -> AxisymmetryMetadata:
    """
    Build axisymmetry metadata from CLI options.

    :param args: Parsed CLI namespace.
    :return: Axisymmetry metadata record.
    :raises ValueError: If axisymmetry metadata is internally inconsistent.
    """
    if not args.axisymmetry_enabled:
        return AxisymmetryMetadata(
            enabled=False,
            symmetry_axis=args.axisymmetry_axis,
            phi_sample_count=0,
            stored_phi_samples=tuple(),
            requires_rotation_for_arbitrary_phi=False,
        )

    phi_samples = tuple(
        float(token.strip()) for token in args.phi_samples.split(",") if token.strip()
    )
    if not phi_samples:
        raise ValueError("Axisymmetry is enabled but --phi-samples was empty.")
    if not all(math.isfinite(value) for value in phi_samples):
        raise ValueError("All phi samples must be finite.")
    if len(set(phi_samples)) != len(phi_samples):
        raise ValueError("Phi samples must be unique.")
    return AxisymmetryMetadata(
        enabled=True,
        symmetry_axis=args.axisymmetry_axis,
        phi_sample_count=len(phi_samples),
        stored_phi_samples=phi_samples,
        requires_rotation_for_arbitrary_phi=args.requires_axisymmetry_rotation,
    )


def determine_spatial_lookup_mode(
    args: argparse.Namespace, include_coordinate_table: bool
) -> Tuple[int, str]:
    """
    Determine the advertised reader lookup mode.

    :param args: Parsed CLI namespace.
    :param include_coordinate_table: Whether a geometry coordinate table is present.
    :return: Tuple of integer mode id and string mode name.
    :raises ValueError: If the requested lookup mode conflicts with available metadata.

    Doctests:
    >>> args = argparse.Namespace(
    ...     spatial_lookup_mode=None,
    ...     native_inverse_map_name=None,
    ... )
    >>> determine_spatial_lookup_mode(args, True)
    (2, 'coordinate_table_only')
    >>> args = argparse.Namespace(
    ...     spatial_lookup_mode="coordinate_table_only",
    ...     native_inverse_map_name=None,
    ... )
    >>> determine_spatial_lookup_mode(args, False)
    Traceback (most recent call last):
    ...
    ValueError: --spatial-lookup-mode=coordinate_table_only requires --include-coordinate-table.
    """
    if args.spatial_lookup_mode is not None:
        mode_name = args.spatial_lookup_mode
    elif include_coordinate_table:
        mode_name = "coordinate_table_only"
    else:
        mode_name = "requires_reader_spatial_index"
    if mode_name == "native_logical_grid" and args.native_inverse_map_name is None:
        raise ValueError(
            "--spatial-lookup-mode=native_logical_grid requires "
            "--native-inverse-map-name."
        )
    if mode_name == "coordinate_table_only" and not include_coordinate_table:
        raise ValueError(
            "--spatial-lookup-mode=coordinate_table_only requires "
            "--include-coordinate-table."
        )
    return SPATIAL_LOOKUP_MODE_IDS[mode_name], mode_name


def extract_coordinate_table_from_first_slice(base: Stage1Info) -> bytes:
    """
    Build the GeometryBlock coordinate table from the first slice payload.

    :param base: First time slice metadata.
    :return: Little-endian float64 coordinate table bytes.
    """
    coordinate_table = bytearray(base.point_record_count * 3 * SERIALIZED_REAL_BYTES)
    with base.path.open("rb") as fp:
        fp.seek(base.point_records_offset)
        for point_index in range(base.point_record_count):
            record = _read_exact(fp, 3 * SERIALIZED_REAL_BYTES, "point coordinates")
            fp.seek(base.point_record_bytes - 3 * SERIALIZED_REAL_BYTES, os.SEEK_CUR)
            dest_offset = point_index * 3 * SERIALIZED_REAL_BYTES
            coordinate_table[dest_offset : dest_offset + 24] = record
    return bytes(coordinate_table)


def _compute_cartesian_bbox_from_coordinate_table(
    coordinate_table: bytes,
) -> Tuple[Tuple[float, float, float], Tuple[float, float, float]]:
    """
    Compute Cartesian bounds from the coordinate table.

    :param coordinate_table: Packed coordinate table bytes.
    :return: Tuple of Cartesian bbox min and max triples.
    :raises RuntimeError: If the coordinate table is empty or contains non-finite data.
    """
    if not coordinate_table:
        raise RuntimeError("Coordinate table must not be empty.")
    unpacked = struct.iter_unpack("<ddd", coordinate_table)
    first_point = next(unpacked)
    xmin, ymin, zmin = first_point
    xmax, ymax, zmax = first_point
    if not all(math.isfinite(value) for value in first_point):
        raise RuntimeError("Coordinate table contained a non-finite Cartesian point.")
    for x_coord, y_coord, z_coord in unpacked:
        if not all(math.isfinite(value) for value in (x_coord, y_coord, z_coord)):
            raise RuntimeError(
                "Coordinate table contained a non-finite Cartesian point."
            )
        xmin = min(xmin, x_coord)
        ymin = min(ymin, y_coord)
        zmin = min(zmin, z_coord)
        xmax = max(xmax, x_coord)
        ymax = max(ymax, y_coord)
        zmax = max(zmax, z_coord)
    return (xmin, ymin, zmin), (xmax, ymax, zmax)


def validate_coordinate_table_against_slice(
    coordinate_table: bytes, info: Stage1Info
) -> None:
    """
    Validate that later slices reuse the same coordinates as the first slice.

    :param coordinate_table: Coordinate table bytes from the first slice.
    :param info: Later stage-1 slice to validate.
    :raises RuntimeError: If a coordinate record is truncated or any x,y,z triple differs.
    """
    with info.path.open("rb") as fp:
        fp.seek(info.point_records_offset)
        for point_index in range(info.point_record_count):
            record = _read_exact(fp, 3 * SERIALIZED_REAL_BYTES, "point coordinates")
            fp.seek(info.point_record_bytes - 3 * SERIALIZED_REAL_BYTES, os.SEEK_CUR)
            src_offset = point_index * 3 * SERIALIZED_REAL_BYTES
            if record != coordinate_table[src_offset : src_offset + 24]:
                raise RuntimeError(
                    f"{info.path}: coordinate table mismatch at point_index="
                    f"{point_index}."
                )


def build_geometry_block(
    base: Stage1Info,
    infos: Sequence[Stage1Info],
    args: argparse.Namespace,
    spatial_lookup_mode_id: int,
) -> GeometryBlockData:
    """
    Build GeometryBlockV1 bytes and derived geometry metadata.

    :param base: First stage-1 slice metadata.
    :param infos: Sorted stage-1 slice metadata.
    :param args: Parsed CLI namespace.
    :param spatial_lookup_mode_id: Declared reader lookup mode id.
    :return: Geometry block bytes plus derived summary values.
    :raises RuntimeError: If geometry-header packing or coordinate validation fails.
    """
    # Step 1: Optionally cache one Cartesian coordinate table. This is not a
    #         true spatial index; it exists so stage 3 can validate or build
    #         its own lookup structure for the nonuniform Cartesian embedding
    #         without scanning the time-major payload.
    if args.include_coordinate_table:
        coordinate_table = extract_coordinate_table_from_first_slice(base)
        if args.validate_coordinate_table:
            for info in infos[1:]:
                validate_coordinate_table_against_slice(coordinate_table, info)
        coordinate_table_kind = GEOMETRY_COORD_TABLE_CARTESIAN_POINT_COORDS_AOS_F64
        coordinate_table_offset_from_block = GEOMETRY_HEADER_BYTES
        coordinate_table_point_count = base.point_record_count
        coordinate_table_components = 3
        coordinate_table_real_bytes = SERIALIZED_REAL_BYTES
        coordinate_table_bytes = len(coordinate_table)
        cartesian_bbox_min, cartesian_bbox_max = (
            _compute_cartesian_bbox_from_coordinate_table(coordinate_table)
        )
    else:
        coordinate_table = b""
        coordinate_table_kind = GEOMETRY_COORD_TABLE_NONE
        coordinate_table_offset_from_block = 0
        coordinate_table_point_count = 0
        coordinate_table_components = 0
        coordinate_table_real_bytes = 0
        coordinate_table_bytes = 0
        cartesian_bbox_min = None
        cartesian_bbox_max = None

    # Step 2: Serialize the fixed GeometryHeaderV1 fields in little-endian order.
    geometry_header = bytearray()
    geometry_header.extend(GEOMETRY_MAGIC)
    geometry_header.extend(struct.pack("<I", GEOMETRY_FORMAT_VERSION))
    geometry_header.extend(struct.pack("<I", GEOMETRY_FLAGS_V1))
    geometry_header.extend(
        struct.pack(
            "<I",
            _require_u32(spatial_lookup_mode_id, "spatial_lookup_mode_id"),
        )
    )
    geometry_header.extend(
        struct.pack("<I", _require_u32(coordinate_table_kind, "coordinate_table_kind"))
    )
    geometry_header.extend(
        struct.pack("<Q", _require_u64(GEOMETRY_HEADER_BYTES, "geometry_header_bytes"))
    )
    geometry_header.extend(
        struct.pack(
            "<Q",
            _require_u64(
                coordinate_table_offset_from_block,
                "coordinate_table_offset_from_block",
            ),
        )
    )
    geometry_header.extend(
        struct.pack(
            "<Q",
            _require_u64(coordinate_table_point_count, "coordinate_table_point_count"),
        )
    )
    geometry_header.extend(
        struct.pack(
            "<I",
            _require_u32(coordinate_table_components, "coordinate_table_components"),
        )
    )
    geometry_header.extend(
        struct.pack(
            "<I",
            _require_u32(coordinate_table_real_bytes, "coordinate_table_real_bytes"),
        )
    )
    geometry_header.extend(
        struct.pack(
            "<Q", _require_u64(coordinate_table_bytes, "coordinate_table_bytes")
        )
    )
    geometry_header.extend(
        struct.pack(
            "<I",
            _require_u32(
                GRID_POINT_POSITION_CONVENTION_ID,
                "grid_point_position_convention_id",
            ),
        )
    )
    geometry_header.extend(struct.pack("<I", 0))
    geometry_header.extend(struct.pack("<8Q", *([0] * 8)))
    if len(geometry_header) > GEOMETRY_HEADER_BYTES:
        raise RuntimeError("GeometryHeaderV1 overflowed its fixed size.")
    geometry_header.extend(b"\0" * (GEOMETRY_HEADER_BYTES - len(geometry_header)))

    return GeometryBlockData(
        block_bytes=bytes(geometry_header) + coordinate_table,
        coordinate_table_present=args.include_coordinate_table,
        coordinate_table_bytes=coordinate_table_bytes,
        coordinate_table_kind=coordinate_table_kind,
        coordinate_table_point_count=coordinate_table_point_count,
        cartesian_bbox_min=cartesian_bbox_min,
        cartesian_bbox_max=cartesian_bbox_max,
    )


def compute_layout(
    base: Stage1Info,
    infos: Sequence[Stage1Info],
    metadata_bytes_len: int,
    geometry_block_bytes_len: int,
    args: argparse.Namespace,
) -> Layout:
    """
    Compute the combined container layout.

    :param base: Reference stage-1 file metadata.
    :param infos: Sorted stage-1 slice metadata.
    :param metadata_bytes_len: Metadata JSON byte count.
    :param geometry_block_bytes_len: Geometry block byte count.
    :param args: Parsed CLI namespace.
    :return: Concrete layout offsets and sizes.
    """
    # Step 1: Align each major block and every payload slice for direct mmap/pread access.
    _validate_power_of_two(args.alignment, "alignment")
    metadata_offset = FIXED_HEADER_BYTES
    slice_table_offset = align_up(metadata_offset + metadata_bytes_len, args.alignment)
    slice_table_bytes = len(infos) * SLICE_TABLE_ENTRY_BYTES
    geometry_block_offset = align_up(
        slice_table_offset + slice_table_bytes, args.alignment
    )
    first_payload_offset = align_up(
        geometry_block_offset + geometry_block_bytes_len, args.alignment
    )
    payload_bytes_per_slice = base.point_record_count * base.point_record_bytes
    payload_stride_bytes = align_up(payload_bytes_per_slice, args.alignment)
    payload_bytes_total = len(infos) * payload_stride_bytes
    total_file_bytes = first_payload_offset + payload_bytes_total
    return Layout(
        alignment_bytes=args.alignment,
        metadata_offset=metadata_offset,
        metadata_bytes=metadata_bytes_len,
        slice_table_offset=slice_table_offset,
        slice_table_bytes=slice_table_bytes,
        geometry_block_offset=geometry_block_offset,
        geometry_block_bytes=geometry_block_bytes_len,
        first_payload_offset=first_payload_offset,
        payload_bytes_per_slice=payload_bytes_per_slice,
        payload_stride_bytes=payload_stride_bytes,
        payload_bytes_total=payload_bytes_total,
        total_file_bytes=total_file_bytes,
    )


def build_slice_entries(
    base: Stage1Info, infos: Sequence[Stage1Info], layout: Layout
) -> List[SliceEntry]:
    """
    Build sorted SliceTableV1 entries.

    :param base: Reference stage-1 file metadata.
    :param infos: Sorted stage-1 slice metadata.
    :param layout: Final combined container layout.
    :return: Slice table entries.
    """
    # Step 1: Record payload offsets after the sort-by-time order is finalized.
    entries: List[SliceEntry] = []
    for slice_index, info in enumerate(infos):
        entries.append(
            SliceEntry(
                output_index=info.output_index,
                slice_flags=0,
                simulation_time=info.simulation_time,
                payload_offset=layout.first_payload_offset
                + slice_index * layout.payload_stride_bytes,
                payload_bytes=layout.payload_bytes_per_slice,
                point_record_count=base.point_record_count,
                point_record_bytes=base.point_record_bytes,
                source_file_size=info.actual_file_size,
                source_header_size=info.header_size,
                source_simulation_time_offset=info.simulation_time_offset,
                source_point_records_offset=info.point_records_offset,
                source_point_records_bytes=info.point_records_bytes,
                source_total_file_bytes=info.total_file_bytes,
            )
        )
    return entries


def build_metadata_json(
    base: Stage1Info,
    infos: Sequence[Stage1Info],
    layout: Layout,
    axisymmetry: AxisymmetryMetadata,
    geometry_data: GeometryBlockData,
    spatial_lookup_mode_name: str,
    native_inverse_map_name: Optional[str],
) -> bytes:
    """
    Build deterministic JSON metadata for the combined container.

    :param base: Reference stage-1 file metadata.
    :param infos: Sorted stage-1 slice metadata.
    :param layout: Final combined layout.
    :param axisymmetry: Axisymmetry metadata.
    :param geometry_data: Derived geometry metadata.
    :param spatial_lookup_mode_name: Reader lookup mode name.
    :param native_inverse_map_name: Reader-declared native inverse-map identifier.
    :return: UTF-8 encoded deterministic JSON metadata.
    """
    # Step 1: Emit a verbose, deterministic schema description for readers and debugging.
    metadata = {
        "combined_format_name": "NRPY raytracing stacked time-slice container",
        "combined_format_version": COMBINED_FORMAT_VERSION,
        "source_format_name": "NRPY stage-1 raytracing data",
        "source_format_version": base.format_version,
        "payload_layout": "time_major_stage1_aos",
        "record_dtype": "<f8",
        "record_component_count": base.record_component_count,
        "point_record_real_count": base.point_record_real_count,
        "point_record_bytes": base.point_record_bytes,
        "loop_order": base.loop_order,
        "format_name": base.format_name,
        "source_coord_system": base.source_coord_system,
        "target_basis": base.target_basis,
        "cf_convention": base.cf_convention,
        "spatial_lookup_mode": spatial_lookup_mode_name,
        "geometry_block": {
            "name": "GeometryBlockV1",
            "contains_true_spatial_index": False,
            "coordinate_table_present": geometry_data.coordinate_table_present,
            "coordinate_table_layout": (
                "cartesian_point_coords[point][3]"
                if geometry_data.coordinate_table_present
                else "none"
            ),
            "coordinate_table_dtype": (
                "<f8" if geometry_data.coordinate_table_present else "none"
            ),
        },
        "axisymmetry": {
            "enabled": axisymmetry.enabled,
            "symmetry_axis": axisymmetry.symmetry_axis,
            "phi_sample_count": axisymmetry.phi_sample_count,
            "stored_phi_samples": list(axisymmetry.stored_phi_samples),
            "requires_rotation_for_arbitrary_phi": (
                axisymmetry.requires_rotation_for_arbitrary_phi
            ),
        },
        "uniformity": {
            "uniform_in_native_coordinates": True,
            "uniform_in_cartesian_coordinates": False,
        },
        "geometry": {
            "grid_topology": "structured_logical_grid",
            "native_coordinate_system": base.source_coord_system,
            "stored_point_coordinates": "Cartesian",
            "stored_tensor_basis": "Cartesian",
            "lookup_warning": (
                "Do not use Cartesian x/y/z with dxx/invdxx directly for lookup."
            ),
        },
        "coordinate_map": {
            "native_coordinates": base.source_coord_system,
            "cartesian_coordinates": "x,y,z",
            "inverse_map_available_to_reader": (
                spatial_lookup_mode_name == "native_logical_grid"
            ),
            "inverse_map_name": native_inverse_map_name,
            "reader_must_not_use_cartesian_uniform_lookup": True,
        },
        "time_lookup_contract": {
            "slice_table_order": "strictly increasing simulation_time",
            "interpolation_owner": "stage3",
            "note": (
                "Stage 2 stores sorted times only. Stage 3 must choose the "
                "lower-inclusive or upper-inclusive bracketing convention used "
                "by the photon integrator."
            ),
        },
        "fixed_header_field_meanings": {
            "xxmin": "native/logical coordinate lower bounds, not Cartesian bounds",
            "xxmax": "native/logical coordinate upper bounds, not Cartesian bounds",
            "dxx": "native/logical uniform grid spacing",
            "invdxx": "inverse native/logical grid spacing",
        },
        "native_grid": {
            "Nxx": list(base.Nxx),
            "Nxx_plus_2NGHOSTS": list(base.Nxx_plus_2NGHOSTS),
            "NGHOSTS": base.NGHOSTS,
            "payload_i_count": list(base.payload_i_count),
            "payload_i_start": list(base.payload_i_start),
            "payload_i_end": list(base.payload_i_end),
            "native_xxmin": list(base.xxmin),
            "native_xxmax": list(base.xxmax),
            "native_dxx": list(base.dxx),
            "native_invdxx": list(base.invdxx),
            "cart_origin": list(base.cart_origin),
        },
        "grid_point_position_convention": {
            "id": GRID_POINT_POSITION_CONVENTION_ID,
            "meaning": (
                "Point records correspond to exported interior logical-grid "
                "points with payload-local indices j = i - payload_i_start."
            ),
        },
        "cartesian_bounds": (
            {
                "min": list(geometry_data.cartesian_bbox_min),
                "max": list(geometry_data.cartesian_bbox_max),
            }
            if geometry_data.cartesian_bbox_min is not None
            and geometry_data.cartesian_bbox_max is not None
            else None
        ),
        "reader_offset_contract": {
            "record_layout_id": "time_major_stage1_aos",
            "lookup_precondition": (
                "The i0,i1,i2 indices must already be native/logical grid "
                "indices. Do not compute them using Cartesian x,y,z with "
                "dxx/invdxx."
            ),
            "point_index_formula": (
                "j0 = i0 - payload_i_start[0]; "
                "j1 = i1 - payload_i_start[1]; "
                "j2 = i2 - payload_i_start[2]; "
                "point_index = j0 + payload_i_count[0] * "
                "(j1 + payload_i_count[1] * j2)"
            ),
            "record_offset_formula": (
                "record_offset = slice_table[s].payload_offset + "
                "point_index * point_record_bytes"
            ),
            "component_offset_formula": (
                "component_offset = record_offset + "
                "component_index * serialized_real_bytes"
            ),
            "metric_component_offset_base": 3,
            "christoffel_component_offset_base": 13,
        },
        "component_names": list(base.record_component_names),
        "metric_component_names": list(base.metric_component_names),
        "christoffel_component_names": list(base.christoffel_component_names),
        "layout_summary": {
            "alignment_bytes": layout.alignment_bytes,
            "num_time_slices": len(infos),
            "metadata_offset": layout.metadata_offset,
            "slice_table_offset": layout.slice_table_offset,
            "geometry_block_offset": layout.geometry_block_offset,
            "first_payload_offset": layout.first_payload_offset,
            "payload_bytes_per_slice": layout.payload_bytes_per_slice,
            "payload_stride_bytes": layout.payload_stride_bytes,
            "payload_bytes_total": layout.payload_bytes_total,
            "total_file_bytes": layout.total_file_bytes,
        },
        "sources": [
            {
                "path": os.path.relpath(info.path, Path.cwd()),
                "output_index": info.output_index,
                "simulation_time": info.simulation_time,
                "file_size": info.actual_file_size,
                "header_size": info.header_size,
                "simulation_time_offset": info.simulation_time_offset,
                "point_records_offset": info.point_records_offset,
                "point_records_bytes": info.point_records_bytes,
                "total_file_bytes": info.total_file_bytes,
            }
            for info in infos
        ],
    }
    return json.dumps(
        metadata,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def pack_slice_entry(entry: SliceEntry) -> bytes:
    """
    Pack one SliceEntryV1.

    :param entry: Slice table entry values.
    :return: Fixed-width packed entry bytes.
    :raises RuntimeError: If the packed entry size differs from the schema.
    """
    packed = struct.pack(
        "<IIdQQQIIQQQQQQ4Q",
        _require_u32(entry.output_index, "slice_entry.output_index"),
        _require_u32(entry.slice_flags, "slice_entry.slice_flags"),
        entry.simulation_time,
        _require_u64(entry.payload_offset, "slice_entry.payload_offset"),
        _require_u64(entry.payload_bytes, "slice_entry.payload_bytes"),
        _require_u64(entry.point_record_count, "slice_entry.point_record_count"),
        _require_u32(entry.point_record_bytes, "slice_entry.point_record_bytes"),
        0,
        _require_u64(entry.source_file_size, "slice_entry.source_file_size"),
        _require_u64(entry.source_header_size, "slice_entry.source_header_size"),
        _require_u64(
            entry.source_simulation_time_offset,
            "slice_entry.source_simulation_time_offset",
        ),
        _require_u64(
            entry.source_point_records_offset,
            "slice_entry.source_point_records_offset",
        ),
        _require_u64(
            entry.source_point_records_bytes,
            "slice_entry.source_point_records_bytes",
        ),
        _require_u64(
            entry.source_total_file_bytes, "slice_entry.source_total_file_bytes"
        ),
        0,
        0,
        0,
        0,
    )
    if len(packed) != SLICE_TABLE_ENTRY_BYTES:
        raise RuntimeError(
            f"SliceEntryV1 packed to {len(packed)} bytes, expected "
            f"{SLICE_TABLE_ENTRY_BYTES}."
        )
    return packed


def pack_fixed_header(
    base: Stage1Info,
    layout: Layout,
    axisymmetry: AxisymmetryMetadata,
    spatial_lookup_mode_id: int,
    num_time_slices: int,
) -> bytes:
    """
    Pack the fixed 4096-byte combined container header.

    :param base: Reference stage-1 file metadata.
    :param layout: Final combined layout.
    :param axisymmetry: Axisymmetry metadata.
    :param spatial_lookup_mode_id: Declared reader lookup mode id.
    :param num_time_slices: Number of sorted slices in the container.
    :return: 4096-byte fixed header.
    :raises RuntimeError: If the packed header exceeds its fixed allocation.
    """
    # Step 1: Pack only the hot-path numeric fields needed for direct offset math.
    header = bytearray()
    header.extend(COMBINED_MAGIC)
    header.extend(
        struct.pack(
            "<I", _require_u32(COMBINED_FORMAT_VERSION, "combined_format_version")
        )
    )
    header.extend(
        struct.pack("<I", _require_u32(base.format_version, "source_format_version"))
    )
    header.extend(struct.pack("<I", _require_u32(ENDIAN_TAG, "endian_tag")))
    header.extend(struct.pack("<I", _require_u32(HEADER_FLAGS_V1, "header_flags")))
    header.extend(
        struct.pack("<Q", _require_u64(FIXED_HEADER_BYTES, "fixed_header_bytes"))
    )
    header.extend(
        struct.pack("<Q", _require_u64(layout.alignment_bytes, "alignment_bytes"))
    )
    header.extend(
        struct.pack("<Q", _require_u64(layout.metadata_offset, "metadata_offset"))
    )
    header.extend(
        struct.pack("<Q", _require_u64(layout.metadata_bytes, "metadata_bytes"))
    )
    header.extend(
        struct.pack("<Q", _require_u64(layout.slice_table_offset, "slice_table_offset"))
    )
    header.extend(
        struct.pack(
            "<Q", _require_u64(SLICE_TABLE_ENTRY_BYTES, "slice_table_entry_bytes")
        )
    )
    header.extend(struct.pack("<Q", _require_u64(num_time_slices, "num_time_slices")))
    header.extend(
        struct.pack(
            "<Q",
            _require_u64(layout.geometry_block_offset, "geometry_block_offset"),
        )
    )
    header.extend(
        struct.pack(
            "<Q", _require_u64(layout.geometry_block_bytes, "geometry_block_bytes")
        )
    )
    header.extend(
        struct.pack(
            "<Q", _require_u64(layout.first_payload_offset, "first_payload_offset")
        )
    )
    header.extend(
        struct.pack(
            "<Q",
            _require_u64(layout.payload_bytes_per_slice, "payload_bytes_per_slice"),
        )
    )
    header.extend(
        struct.pack(
            "<Q", _require_u64(layout.payload_stride_bytes, "payload_stride_bytes")
        )
    )
    header.extend(
        struct.pack(
            "<Q", _require_u64(layout.payload_bytes_total, "payload_bytes_total")
        )
    )
    header.extend(
        struct.pack("<Q", _require_u64(layout.total_file_bytes, "total_file_bytes"))
    )
    header.extend(
        struct.pack("<Q", _require_u64(base.point_record_count, "point_record_count"))
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.serialized_real_bytes, "serialized_real_bytes")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.point_record_real_count, "point_record_real_count")
        )
    )
    header.extend(
        struct.pack("<I", _require_u32(base.point_record_bytes, "point_record_bytes"))
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.record_component_count, "record_component_count")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.metric_component_count, "metric_component_count")
        )
    )
    header.extend(
        struct.pack(
            "<I",
            _require_u32(
                base.christoffel_component_count, "christoffel_component_count"
            ),
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.file_is_little_endian, "file_is_little_endian")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(base.time_variable_is_f64, "time_variable_is_f64")
        )
    )
    header.extend(
        struct.pack(
            "<I",
            _require_u32(
                base.payload_includes_ghost_zones, "payload_includes_ghost_zones"
            ),
        )
    )
    header.extend(struct.pack("<I", _require_u32(base.num_grids, "num_grids")))
    header.extend(struct.pack("<I", _require_u32(base.NGHOSTS, "NGHOSTS")))
    header.extend(
        struct.pack(
            "<I", _require_u32(RECORD_LAYOUT_TIME_MAJOR_STAGE1_AOS, "record_layout_id")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(spatial_lookup_mode_id, "spatial_lookup_mode_id")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(1 if axisymmetry.enabled else 0, "axisymmetry_flag")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(axisymmetry.phi_sample_count, "phi_sample_count")
        )
    )
    header.extend(
        struct.pack(
            "<I", _require_u32(GEOMETRY_FORMAT_VERSION, "geometry_format_version")
        )
    )
    header.extend(struct.pack("<I", 0))
    header.extend(
        struct.pack("<3Q", *[_require_u64(value, "Nxx") for value in base.Nxx])
    )
    header.extend(
        struct.pack(
            "<3Q",
            *[
                _require_u64(value, "Nxx_plus_2NGHOSTS")
                for value in base.Nxx_plus_2NGHOSTS
            ],
        )
    )
    header.extend(
        struct.pack(
            "<3Q",
            *[_require_u64(value, "payload_i_count") for value in base.payload_i_count],
        )
    )
    header.extend(
        struct.pack(
            "<3q",
            *[_require_i64(value, "payload_i_start") for value in base.payload_i_start],
        )
    )
    header.extend(
        struct.pack(
            "<3q",
            *[_require_i64(value, "payload_i_end") for value in base.payload_i_end],
        )
    )
    header.extend(struct.pack("<3d", *base.dxx))
    header.extend(struct.pack("<3d", *base.invdxx))
    header.extend(struct.pack("<3d", *base.xxmin))
    header.extend(struct.pack("<3d", *base.xxmax))
    header.extend(struct.pack("<3d", *base.cart_origin))
    if len(header) > FIXED_HEADER_BYTES:
        raise RuntimeError("FixedHeaderV1 overflowed its 4096-byte allocation.")
    header.extend(b"\0" * (FIXED_HEADER_BYTES - len(header)))
    return bytes(header)


def copy_payload(
    src_path: Path,
    src_offset: int,
    num_bytes: int,
    dst_file: BinaryIO,
    dst_offset: int,
) -> None:
    """
    Copy one stage-1 payload region into the combined output.

    :param src_path: Source stage-1 path.
    :param src_offset: Source payload offset.
    :param num_bytes: Number of bytes to copy.
    :param dst_file: Open destination file.
    :param dst_offset: Destination payload offset.
    :raises RuntimeError: If EOF is encountered while copying.
    """
    chunk_size = 16 * 1024 * 1024
    with src_path.open("rb") as src:
        src.seek(src_offset)
        dst_file.seek(dst_offset)
        remaining = num_bytes
        while remaining > 0:
            chunk = src.read(min(chunk_size, remaining))
            if not chunk:
                raise RuntimeError(
                    f"Unexpected EOF while copying payload from {src_path}."
                )
            dst_file.write(chunk)
            remaining -= len(chunk)


def _fsync_parent_directory(path: Path) -> None:
    """
    Fsync the directory containing a file path.

    :param path: Final output path.
    """
    directory_fd = os.open(str(path.parent), os.O_RDONLY)
    try:
        os.fsync(directory_fd)
    finally:
        os.close(directory_fd)


def write_combined_file_atomically(
    output_path: Path,
    base: Stage1Info,
    infos: Sequence[Stage1Info],
    layout: Layout,
    metadata_bytes: bytes,
    geometry_data: GeometryBlockData,
    slice_entries: Sequence[SliceEntry],
    axisymmetry: AxisymmetryMetadata,
    spatial_lookup_mode_id: int,
    force: bool,
) -> None:
    """
    Write the combined container atomically.

    :param output_path: Final destination path.
    :param base: Reference stage-1 metadata.
    :param infos: Sorted stage-1 metadata.
    :param layout: Final layout.
    :param metadata_bytes: Metadata JSON bytes.
    :param geometry_data: Geometry block bytes and summary metadata.
    :param slice_entries: Slice table entries.
    :param axisymmetry: Axisymmetry metadata.
    :param spatial_lookup_mode_id: Declared reader lookup mode id.
    :param force: Whether overwriting an existing file is allowed.
    :raises RuntimeError: If overwrite rules fail during the final install step.
    :raises Exception: Propagates I/O failures after best-effort temporary cleanup.
    """
    # Step 1: Write into a same-directory temporary file so final installation
    #         stays atomic on the target filesystem.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_name = tempfile.mkstemp(
        prefix=f"{output_path.name}.tmp.", dir=str(output_path.parent)
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(tmp_fd, "wb") as out:
            # Step 2: Pre-size the file so aligned gaps are already zero-filled.
            _require_u64(layout.total_file_bytes, "layout.total_file_bytes")
            out.truncate(layout.total_file_bytes)
            header_bytes = pack_fixed_header(
                base=base,
                layout=layout,
                axisymmetry=axisymmetry,
                spatial_lookup_mode_id=spatial_lookup_mode_id,
                num_time_slices=len(infos),
            )
            out.seek(0)
            out.write(header_bytes)
            out.seek(layout.metadata_offset)
            out.write(metadata_bytes)
            out.seek(layout.slice_table_offset)
            for entry in slice_entries:
                out.write(pack_slice_entry(entry))
            out.seek(layout.geometry_block_offset)
            out.write(geometry_data.block_bytes)
            # Step 3: Copy only the stage-1 payload region for each sorted slice.
            for slice_index, info in enumerate(infos):
                copy_payload(
                    src_path=info.path,
                    src_offset=info.point_records_offset,
                    num_bytes=info.point_records_bytes,
                    dst_file=out,
                    dst_offset=slice_entries[slice_index].payload_offset,
                )
            out.flush()
            os.fsync(out.fileno())
        current_mode = tmp_path.stat().st_mode & 0o777
        # Step 3.a: mkstemp creates an owner-only file, so this normally
        #           installs the final container as owner-read-only. Use a
        #           wider chmod policy if readers run under another user or group.
        os.chmod(tmp_path, current_mode & ~0o222)

        # Step 4: Install the completed file only after data and metadata are durable.
        if force:
            os.replace(tmp_path, output_path)
        else:
            # Step 4.a: In non-force mode, use link+unlink instead of replace so
            #           an existing output cannot be overwritten by a race
            #           between the final existence check and installation.
            try:
                os.link(tmp_path, output_path)
            except FileExistsError as exc:
                raise RuntimeError(
                    f"Output file {output_path} exists. Use --force to overwrite."
                ) from exc
            tmp_path.unlink()
        _fsync_parent_directory(output_path)
    except Exception:
        if tmp_path.exists():
            tmp_path.unlink()
        raise


def _unpack_fixed_header(header_bytes: bytes) -> CombinedHeaderInfo:
    """
    Unpack FixedHeaderV1 into an immutable parsed-header record.

    :param header_bytes: First 4096 bytes of the combined file.
    :return: Parsed fixed header fields.
    """
    # Step 1: Consume the fixed header in the exact serialization order.
    offset = 0

    def unpack(fmt: str) -> Tuple[object, ...]:
        nonlocal offset
        size = struct.calcsize(fmt)
        values = struct.unpack_from(fmt, header_bytes, offset)
        offset += size
        return values

    magic_bytes = cast(bytes, unpack("<16s")[0])
    magic = magic_bytes.split(b"\0", 1)[0].decode("ascii")
    combined_u32 = cast(Tuple[int, int, int, int], unpack("<4I"))
    combined_u64 = cast(
        Tuple[
            int, int, int, int, int, int, int, int, int, int, int, int, int, int, int
        ],
        unpack("<15Q"),
    )
    trailing_u32 = cast(
        Tuple[
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
        ],
        unpack("<17I"),
    )
    return CombinedHeaderInfo(
        magic=magic,
        combined_format_version=combined_u32[0],
        source_format_version=combined_u32[1],
        endian_tag=combined_u32[2],
        header_flags=combined_u32[3],
        fixed_header_bytes=combined_u64[0],
        alignment_bytes=combined_u64[1],
        metadata_offset=combined_u64[2],
        metadata_bytes=combined_u64[3],
        slice_table_offset=combined_u64[4],
        slice_table_entry_bytes=combined_u64[5],
        num_time_slices=combined_u64[6],
        geometry_block_offset=combined_u64[7],
        geometry_block_bytes=combined_u64[8],
        first_payload_offset=combined_u64[9],
        payload_bytes_per_slice=combined_u64[10],
        payload_stride_bytes=combined_u64[11],
        payload_bytes_total=combined_u64[12],
        total_file_bytes=combined_u64[13],
        point_record_count=combined_u64[14],
        serialized_real_bytes=trailing_u32[0],
        point_record_real_count=trailing_u32[1],
        point_record_bytes=trailing_u32[2],
        record_component_count=trailing_u32[3],
        metric_component_count=trailing_u32[4],
        christoffel_component_count=trailing_u32[5],
        file_is_little_endian=trailing_u32[6],
        time_variable_is_f64=trailing_u32[7],
        payload_includes_ghost_zones=trailing_u32[8],
        num_grids=trailing_u32[9],
        NGHOSTS=trailing_u32[10],
        record_layout_id=trailing_u32[11],
        spatial_lookup_mode_id=trailing_u32[12],
        axisymmetry_flag=trailing_u32[13],
        phi_sample_count=trailing_u32[14],
        geometry_format_version=trailing_u32[15],
        reserved_u32=trailing_u32[16],
        Nxx=cast(Tuple[int, int, int], unpack("<3Q")),
        Nxx_plus_2NGHOSTS=cast(Tuple[int, int, int], unpack("<3Q")),
        payload_i_count=cast(Tuple[int, int, int], unpack("<3Q")),
        payload_i_start=cast(Tuple[int, int, int], unpack("<3q")),
        payload_i_end=cast(Tuple[int, int, int], unpack("<3q")),
        dxx=cast(Tuple[float, float, float], unpack("<3d")),
        invdxx=cast(Tuple[float, float, float], unpack("<3d")),
        xxmin=cast(Tuple[float, float, float], unpack("<3d")),
        xxmax=cast(Tuple[float, float, float], unpack("<3d")),
        cart_origin=cast(Tuple[float, float, float], unpack("<3d")),
    )


def _unpack_geometry_header(geometry_header_bytes: bytes) -> GeometryHeaderInfo:
    """
    Unpack GeometryHeaderV1 into an immutable parsed-header record.

    :param geometry_header_bytes: First 512 bytes of the geometry block.
    :return: Parsed geometry header fields.
    """
    # Step 1: Consume the fixed geometry header in the exact serialization order.
    offset = 0

    def unpack(fmt: str) -> Tuple[object, ...]:
        nonlocal offset
        size = struct.calcsize(fmt)
        values = struct.unpack_from(fmt, geometry_header_bytes, offset)
        offset += size
        return values

    magic_bytes = cast(bytes, unpack("<16s")[0])
    magic = magic_bytes.split(b"\0", 1)[0].decode("ascii")
    header_u32 = cast(Tuple[int, int, int, int], unpack("<4I"))
    header_u64 = cast(Tuple[int, int, int], unpack("<3Q"))
    coord_u32 = cast(Tuple[int, int], unpack("<2I"))
    coordinate_table_bytes = cast(int, unpack("<Q")[0])
    footer_u32 = cast(Tuple[int, int], unpack("<2I"))
    return GeometryHeaderInfo(
        magic=magic,
        geometry_format_version=header_u32[0],
        geometry_flags=header_u32[1],
        spatial_lookup_mode_id=header_u32[2],
        coordinate_table_kind=header_u32[3],
        geometry_header_bytes=header_u64[0],
        coordinate_table_offset_from_block=header_u64[1],
        coordinate_table_point_count=header_u64[2],
        coordinate_table_components=coord_u32[0],
        coordinate_table_real_bytes=coord_u32[1],
        coordinate_table_bytes=coordinate_table_bytes,
        grid_point_position_convention_id=footer_u32[0],
        reserved_u32=footer_u32[1],
        reserved_u64=cast(Tuple[int, int, int, int, int, int, int, int], unpack("<8Q")),
    )


def _unpack_slice_entry(raw_entry: bytes) -> SliceEntryInfo:
    """
    Unpack one fixed-width slice-table entry.

    :param raw_entry: One 128-byte slice-table record.
    :return: Parsed slice-table entry fields.
    """
    unpacked = cast(
        Tuple[
            int,
            int,
            float,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
            int,
        ],
        struct.unpack("<IIdQQQIIQQQQQQ4Q", raw_entry),
    )
    return SliceEntryInfo(
        output_index=unpacked[0],
        slice_flags=unpacked[1],
        simulation_time=unpacked[2],
        payload_offset=unpacked[3],
        payload_bytes=unpacked[4],
        point_record_count=unpacked[5],
        point_record_bytes=unpacked[6],
        reserved_u32=unpacked[7],
        source_file_size=unpacked[8],
        source_header_size=unpacked[9],
        source_simulation_time_offset=unpacked[10],
        source_point_records_offset=unpacked[11],
        source_point_records_bytes=unpacked[12],
        source_total_file_bytes=unpacked[13],
        reserved_u64=(unpacked[14], unpacked[15], unpacked[16], unpacked[17]),
    )


def inspect_combined_file(path: Path) -> None:
    """
    Inspect a combined file and print a concise summary.

    :param path: Combined container path.
    :raises RuntimeError: If validation fails.
    """
    # Step 1: Parse the fixed header and the JSON/geometry blocks they point to.
    file_size = path.stat().st_size
    with path.open("rb") as fp:
        header_bytes = _read_exact(fp, FIXED_HEADER_BYTES, "fixed header")
        header = _unpack_fixed_header(header_bytes)
        if header.magic != COMBINED_MAGIC.split(b"\0", 1)[0].decode("ascii"):
            raise RuntimeError(f"{path}: combined magic was invalid.")
        if header.combined_format_version != COMBINED_FORMAT_VERSION:
            raise RuntimeError(f"{path}: unsupported combined_format_version.")
        if header.source_format_version != STAGE1_FORMAT_VERSION:
            raise RuntimeError(f"{path}: unsupported source_format_version.")
        if header.endian_tag != ENDIAN_TAG:
            raise RuntimeError(f"{path}: endian_tag was invalid.")
        if header.header_flags != HEADER_FLAGS_V1:
            raise RuntimeError(f"{path}: header_flags was invalid.")
        if header.fixed_header_bytes != FIXED_HEADER_BYTES:
            raise RuntimeError(f"{path}: fixed_header_bytes was invalid.")
        if header.slice_table_entry_bytes != SLICE_TABLE_ENTRY_BYTES:
            raise RuntimeError(f"{path}: slice_table_entry_bytes was invalid.")
        if header.record_layout_id != RECORD_LAYOUT_TIME_MAJOR_STAGE1_AOS:
            raise RuntimeError(f"{path}: record_layout_id was invalid.")
        if header.serialized_real_bytes != SERIALIZED_REAL_BYTES:
            raise RuntimeError(f"{path}: serialized_real_bytes was invalid.")
        if header.point_record_real_count != POINT_RECORD_REAL_COUNT:
            raise RuntimeError(f"{path}: point_record_real_count was invalid.")
        if header.point_record_bytes != POINT_RECORD_BYTES:
            raise RuntimeError(f"{path}: point_record_bytes was invalid.")
        if header.record_component_count != RECORD_COMPONENT_COUNT:
            raise RuntimeError(f"{path}: record_component_count was invalid.")
        if header.metric_component_count != METRIC_COMPONENT_COUNT:
            raise RuntimeError(f"{path}: metric_component_count was invalid.")
        if header.christoffel_component_count != CHRISTOFFEL_COMPONENT_COUNT:
            raise RuntimeError(f"{path}: christoffel_component_count was invalid.")
        if header.file_is_little_endian != 1:
            raise RuntimeError(f"{path}: file_is_little_endian was invalid.")
        if header.time_variable_is_f64 != 1:
            raise RuntimeError(f"{path}: time_variable_is_f64 was invalid.")
        if header.payload_includes_ghost_zones != 0:
            raise RuntimeError(f"{path}: payload_includes_ghost_zones was invalid.")
        if header.num_grids != 1:
            raise RuntimeError(f"{path}: num_grids was invalid.")
        if header.reserved_u32 != 0:
            raise RuntimeError(f"{path}: reserved_u32 was invalid.")
        if header.axisymmetry_flag not in (0, 1):
            raise RuntimeError(f"{path}: axisymmetry_flag was invalid.")
        if header.geometry_format_version != GEOMETRY_FORMAT_VERSION:
            raise RuntimeError(f"{path}: geometry_format_version was invalid.")
        _validate_power_of_two(header.alignment_bytes, "header.alignment_bytes")
        if file_size != header.total_file_bytes:
            raise RuntimeError(
                f"{path}: file size {file_size} did not match total_file_bytes="
                f"{header.total_file_bytes}."
            )
        if header.spatial_lookup_mode_id not in SPATIAL_LOOKUP_MODE_NAMES:
            raise RuntimeError(
                f"{path}: unknown spatial_lookup_mode_id={header.spatial_lookup_mode_id}."
            )
        if header.metadata_offset != FIXED_HEADER_BYTES:
            raise RuntimeError(f"{path}: metadata_offset was invalid.")
        if header.slice_table_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: slice_table_offset was not aligned.")
        if header.geometry_block_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: geometry_block_offset was not aligned.")
        if header.metadata_offset + header.metadata_bytes > file_size:
            raise RuntimeError(f"{path}: metadata block exceeded file bounds.")
        if (
            header.slice_table_offset
            + header.num_time_slices * header.slice_table_entry_bytes
            > file_size
        ):
            raise RuntimeError(f"{path}: slice table exceeded file bounds.")
        if header.geometry_block_offset + header.geometry_block_bytes > file_size:
            raise RuntimeError(f"{path}: geometry block exceeded file bounds.")
        if header.first_payload_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: first_payload_offset was not aligned.")
        if header.payload_bytes_per_slice != (
            header.point_record_count * header.point_record_bytes
        ):
            raise RuntimeError(f"{path}: payload_bytes_per_slice was invalid.")
        if header.payload_stride_bytes < header.payload_bytes_per_slice:
            raise RuntimeError(f"{path}: payload_stride_bytes was invalid.")
        if header.payload_stride_bytes % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: payload_stride_bytes was not aligned.")
        if header.payload_bytes_total != (
            header.num_time_slices * header.payload_stride_bytes
        ):
            raise RuntimeError(f"{path}: payload_bytes_total was invalid.")
        if header.total_file_bytes != (
            header.first_payload_offset + header.payload_bytes_total
        ):
            raise RuntimeError(f"{path}: total_file_bytes was invalid.")
        metadata_end = header.metadata_offset + header.metadata_bytes
        slice_table_end = (
            header.slice_table_offset
            + header.num_time_slices * header.slice_table_entry_bytes
        )
        geometry_block_end = header.geometry_block_offset + header.geometry_block_bytes
        if metadata_end > header.slice_table_offset:
            raise RuntimeError(f"{path}: metadata block overlapped the slice table.")
        if slice_table_end > header.geometry_block_offset:
            raise RuntimeError(f"{path}: slice table overlapped the geometry block.")
        if geometry_block_end > header.first_payload_offset:
            raise RuntimeError(f"{path}: geometry block overlapped the payload region.")

        fp.seek(header.metadata_offset)
        metadata_bytes = _read_exact(fp, header.metadata_bytes, "metadata JSON")
        metadata = json.loads(metadata_bytes.decode("utf-8"))

        fp.seek(header.geometry_block_offset)
        geometry_header_bytes = _read_exact(
            fp, GEOMETRY_HEADER_BYTES, "geometry header"
        )
        geometry_header = _unpack_geometry_header(geometry_header_bytes)
        if geometry_header.magic != GEOMETRY_MAGIC.split(b"\0", 1)[0].decode("ascii"):
            raise RuntimeError(f"{path}: geometry magic was invalid.")
        if geometry_header.geometry_format_version != GEOMETRY_FORMAT_VERSION:
            raise RuntimeError(f"{path}: geometry header version was invalid.")
        if geometry_header.geometry_flags != GEOMETRY_FLAGS_V1:
            raise RuntimeError(f"{path}: geometry_flags was invalid.")
        if geometry_header.geometry_header_bytes != GEOMETRY_HEADER_BYTES:
            raise RuntimeError(f"{path}: geometry_header_bytes was invalid.")
        if geometry_header.spatial_lookup_mode_id != header.spatial_lookup_mode_id:
            raise RuntimeError(
                f"{path}: geometry spatial_lookup_mode_id differed from header."
            )
        if geometry_header.grid_point_position_convention_id != (
            GRID_POINT_POSITION_CONVENTION_ID
        ):
            raise RuntimeError(
                f"{path}: grid_point_position_convention_id was invalid."
            )
        if geometry_header.reserved_u32 != 0 or any(geometry_header.reserved_u64):
            raise RuntimeError(f"{path}: geometry reserved fields were nonzero.")
        if geometry_header.coordinate_table_kind == GEOMETRY_COORD_TABLE_NONE:
            if geometry_header.coordinate_table_offset_from_block != 0:
                raise RuntimeError(
                    f"{path}: coordinate_table_offset_from_block was invalid."
                )
            if geometry_header.coordinate_table_point_count != 0:
                raise RuntimeError(f"{path}: coordinate_table_point_count was invalid.")
            if geometry_header.coordinate_table_components != 0:
                raise RuntimeError(f"{path}: coordinate_table_components was invalid.")
            if geometry_header.coordinate_table_real_bytes != 0:
                raise RuntimeError(f"{path}: coordinate_table_real_bytes was invalid.")
            if geometry_header.coordinate_table_bytes != 0:
                raise RuntimeError(f"{path}: coordinate_table_bytes was invalid.")
        elif (
            geometry_header.coordinate_table_kind
            == GEOMETRY_COORD_TABLE_CARTESIAN_POINT_COORDS_AOS_F64
        ):
            if (
                geometry_header.coordinate_table_offset_from_block
                != GEOMETRY_HEADER_BYTES
            ):
                raise RuntimeError(
                    f"{path}: coordinate_table_offset_from_block was invalid."
                )
            if (
                geometry_header.coordinate_table_point_count
                != header.point_record_count
            ):
                raise RuntimeError(f"{path}: coordinate_table_point_count was invalid.")
            if geometry_header.coordinate_table_components != 3:
                raise RuntimeError(f"{path}: coordinate_table_components was invalid.")
            if geometry_header.coordinate_table_real_bytes != SERIALIZED_REAL_BYTES:
                raise RuntimeError(f"{path}: coordinate_table_real_bytes was invalid.")
            if geometry_header.coordinate_table_bytes != (
                header.point_record_count * 3 * SERIALIZED_REAL_BYTES
            ):
                raise RuntimeError(f"{path}: coordinate_table_bytes was invalid.")
        else:
            raise RuntimeError(
                f"{path}: unknown coordinate_table_kind="
                f"{geometry_header.coordinate_table_kind}."
            )
        if (
            geometry_header.coordinate_table_offset_from_block
            + geometry_header.coordinate_table_bytes
            > header.geometry_block_bytes
        ):
            raise RuntimeError(
                f"{path}: coordinate table exceeded geometry block bounds."
            )

        fp.seek(header.slice_table_offset)
        slice_entries: List[SliceEntryInfo] = []
        num_time_slices = header.num_time_slices
        for _ in range(num_time_slices):
            raw_entry = _read_exact(fp, SLICE_TABLE_ENTRY_BYTES, "slice table entry")
            slice_entries.append(_unpack_slice_entry(raw_entry))

    # Step 2: Validate that the offsets and times are consistent with file bounds.
    times = [entry.simulation_time for entry in slice_entries]
    previous_time: Optional[float] = None
    for time_value in times:
        if previous_time is not None and time_value <= previous_time:
            raise RuntimeError(f"{path}: slice times were not strictly increasing.")
        previous_time = time_value
    for entry in slice_entries:
        if entry.slice_flags != 0:
            raise RuntimeError(f"{path}: slice_flags was invalid.")
        if entry.reserved_u32 != 0 or any(entry.reserved_u64):
            raise RuntimeError(f"{path}: slice reserved fields were nonzero.")
        payload_offset = entry.payload_offset
        payload_bytes = entry.payload_bytes
        if payload_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: slice payload offset was not aligned.")
        if payload_bytes != header.payload_bytes_per_slice:
            raise RuntimeError(f"{path}: slice payload_bytes differed from header.")
        if entry.point_record_count != header.point_record_count:
            raise RuntimeError(
                f"{path}: slice point_record_count differed from header."
            )
        if entry.point_record_bytes != header.point_record_bytes:
            raise RuntimeError(
                f"{path}: slice point_record_bytes differed from header."
            )
        if payload_offset + payload_bytes > file_size:
            raise RuntimeError(f"{path}: slice payload exceeded file bounds.")
    for slice_index, entry in enumerate(slice_entries):
        expected_payload_offset = (
            header.first_payload_offset + slice_index * header.payload_stride_bytes
        )
        if entry.payload_offset != expected_payload_offset:
            raise RuntimeError(
                f"{path}: slice payload offset did not match stride formula."
            )

    # Step 3: Print a concise human-readable summary for quick inspection.
    print(f"magic: {header.magic}")
    print(f"combined_format_version: {header.combined_format_version}")
    print(f"source_format_version: {header.source_format_version}")
    print(f"total_file_bytes: {header.total_file_bytes}")
    print(f"num_time_slices: {header.num_time_slices}")
    if times:
        print("time range:")
        print(f"  first simulation_time: {times[0]:.17g}")
        print(f"  last simulation_time: {times[-1]:.17g}")
    print("payload:")
    print(f"  first_payload_offset: {header.first_payload_offset}")
    print(f"  payload_bytes_per_slice: {header.payload_bytes_per_slice}")
    print(f"  payload_stride_bytes: {header.payload_stride_bytes}")
    print(f"  point_record_count: {header.point_record_count}")
    print(f"  point_record_bytes: {header.point_record_bytes}")
    print("native/logical grid:")
    print(f"  Nxx: {header.Nxx}")
    print(f"  payload_i_start: {header.payload_i_start}")
    print(f"  payload_i_count: {header.payload_i_count}")
    print(f"  payload_i_end: {header.payload_i_end}")
    print(f"  native_dxx: {header.dxx}")
    print(f"  native_xxmin: {header.xxmin}")
    print(f"  native_xxmax: {header.xxmax}")
    print("  warning: dxx/invdxx are not Cartesian spacings.")
    print("geometry:")
    print(
        "  spatial_lookup_mode: "
        f"{SPATIAL_LOOKUP_MODE_NAMES[header.spatial_lookup_mode_id]}"
    )
    print(
        "  coordinate_table_present: "
        f"{geometry_header.coordinate_table_kind != GEOMETRY_COORD_TABLE_NONE}"
    )
    print(f"  coordinate_table_bytes: {geometry_header.coordinate_table_bytes}")
    print(f"  cartesian_bounds: {metadata.get('cartesian_bounds')}")
    axisymmetry = metadata.get("axisymmetry", {})
    print("axisymmetry:")
    print(f"  enabled: {axisymmetry.get('enabled')}")
    print(f"  symmetry_axis: {axisymmetry.get('symmetry_axis')}")
    print(f"  phi_sample_count: {axisymmetry.get('phi_sample_count')}")
    print(
        "  requires_rotation_for_arbitrary_phi: "
        f"{axisymmetry.get('requires_rotation_for_arbitrary_phi')}"
    )

    head_entries = slice_entries[: min(3, len(slice_entries))]
    tail_entries = slice_entries[-min(3, len(slice_entries)) :]
    print("first few slice table entries:")
    for entry in head_entries:
        print(
            "  "
            f"output_index={entry.output_index} "
            f"simulation_time={entry.simulation_time:.17g} "
            f"payload_offset={entry.payload_offset} payload_bytes={entry.payload_bytes}"
        )
    print("last few slice table entries:")
    for entry in tail_entries:
        print(
            "  "
            f"output_index={entry.output_index} "
            f"simulation_time={entry.simulation_time:.17g} "
            f"payload_offset={entry.payload_offset} payload_bytes={entry.payload_bytes}"
        )


def main() -> None:
    """
    Drive CLI execution.

    :return: None.
    :raises RuntimeError: If discovery, validation, layout, or output writing fails.
    """
    # Step 1: Dispatch inspect mode before any stage-1 discovery or validation.
    args = parse_args()
    if args.inspect:
        inspect_combined_file(Path(args.inspect))
        return

    # Step 2: Parse and validate every stage-1 source against one strict schema.
    input_paths = discover_input_paths(args)
    if not input_paths:
        raise RuntimeError("No input stage-1 raytracing files found.")

    output_path = Path(args.output).resolve()
    if output_path in input_paths:
        raise RuntimeError("Output path must not also be one of the input files.")
    if output_path.exists() and not args.force:
        raise RuntimeError(
            f"Output file {output_path} exists. Use --force to overwrite."
        )

    infos = [parse_stage1_file(path) for path in input_paths]
    for info in infos:
        validate_stage1_file_internal(info)
    infos.sort(key=lambda info: info.simulation_time)
    validate_sorted_slices(infos)
    base = infos[0]
    for info in infos[1:]:
        validate_stage1_compatible(base, info)

    # Step 3: Build the side metadata blocks before computing the final layout.
    axisymmetry = determine_axisymmetry_metadata(args)
    spatial_lookup_mode_id, spatial_lookup_mode_name = determine_spatial_lookup_mode(
        args, args.include_coordinate_table
    )
    geometry_data = build_geometry_block(base, infos, args, spatial_lookup_mode_id)
    layout = compute_layout(
        base,
        infos,
        metadata_bytes_len=0,
        geometry_block_bytes_len=len(geometry_data.block_bytes),
        args=args,
    )
    # Step 4: Iterate until the metadata size and aligned offsets stop changing.
    metadata_bytes = b""
    for _ in range(8):
        metadata_bytes = build_metadata_json(
            base=base,
            infos=infos,
            layout=layout,
            axisymmetry=axisymmetry,
            geometry_data=geometry_data,
            spatial_lookup_mode_name=spatial_lookup_mode_name,
            native_inverse_map_name=args.native_inverse_map_name,
        )
        next_layout = compute_layout(
            base,
            infos,
            metadata_bytes_len=len(metadata_bytes),
            geometry_block_bytes_len=len(geometry_data.block_bytes),
            args=args,
        )
        if next_layout == layout:
            break
        layout = next_layout
    else:
        raise RuntimeError("Combined layout did not stabilize after metadata rebuild.")

    final_layout = layout
    if final_layout.metadata_bytes != len(metadata_bytes):
        raise RuntimeError(
            "Final layout metadata_bytes did not match the serialized metadata size."
        )
    if final_layout.geometry_block_bytes != len(geometry_data.block_bytes):
        raise RuntimeError(
            "Final layout geometry_block_bytes did not match the geometry block size."
        )

    # Step 5: Materialize the slice table and atomically install the output file.
    slice_entries = build_slice_entries(base, infos, final_layout)
    if len(slice_entries) != len(infos):
        raise RuntimeError(
            "Slice table entry count did not match the input slice count."
        )
    _require_u64(final_layout.total_file_bytes, "final_layout.total_file_bytes")
    write_combined_file_atomically(
        output_path=output_path,
        base=base,
        infos=infos,
        layout=final_layout,
        metadata_bytes=metadata_bytes,
        geometry_data=geometry_data,
        slice_entries=slice_entries,
        axisymmetry=axisymmetry,
        spatial_lookup_mode_id=spatial_lookup_mode_id,
        force=args.force,
    )


if __name__ == "__main__":
    if len(sys.argv) == 1:
        import doctest

        results = doctest.testmod()

        if results.failed > 0:
            print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
            sys.exit(1)
        print(f"Doctest passed: All {results.attempted} test(s) passed")
    else:
        main()

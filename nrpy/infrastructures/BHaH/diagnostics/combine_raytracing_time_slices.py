"""
Combine raytracing time-slice files into one storage container.

This module parses the binary files written by
``nrpy.infrastructures.BHaH.diagnostics.output_raytracing_data``, reads enough
of each input slice to locate the physical simulation time and point-record
payload bytes, validates byte-size and offset integrity, sorts slices by
``simulation_time``, copies payload bytes verbatim, and embeds run metadata.

The combined file is a storage combiner. It does not infer geometry,
axisymmetry, coordinate maps, spatial lookup behavior, ghost-zone conventions,
or reader interpolation policy.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import json
import math
import os
import struct
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import BinaryIO, Dict, List, Optional, Sequence, Tuple, cast

INPUT_SLICE_MAGIC = "NRPYRTDATA4D"
COMBINED_MAGIC = "NRPYRTSTACK4D"

FIXED_HEADER_BYTES = 4096
DEFAULT_ALIGNMENT = 4096
FIXED_HEADER_FORMAT = "<16s10Q"
SLICE_TABLE_ENTRY_FORMAT = "<Qd10Q"
SLICE_TABLE_ENTRY_BYTES = struct.calcsize(SLICE_TABLE_ENTRY_FORMAT)


@dataclass(frozen=True)
class GridMetadata:
    """
    Grid metadata parsed from one source slice.

    The combiner does not validate that these values match across slices.
    When serialized into the combined metadata, they are recorded explicitly
    as first-slice provenance only.
    """

    Nxx: Tuple[int, int, int]
    Nxx_plus_2NGHOSTS: Tuple[int, int, int]
    NGHOSTS: int
    dxx: Tuple[float, float, float]
    invdxx: Tuple[float, float, float]
    xxmin: Tuple[float, float, float]
    xxmax: Tuple[float, float, float]
    cart_origin: Tuple[float, float, float]


@dataclass(frozen=True)
class InputSliceInfo:
    """
    Minimal parsed source-slice metadata.

    Only byte-location fields, simulation time, and the per-slice grid
    metadata retained for first-slice provenance are stored.
    """

    path: Path
    magic: str
    header_size: int
    output_index: int
    point_record_count: int
    point_record_bytes: int
    simulation_time_offset: int
    point_records_offset: int
    point_records_bytes: int
    total_file_bytes: int
    actual_file_size: int
    simulation_time: float
    grid_metadata: GridMetadata


@dataclass(frozen=True)
class Layout:
    """
    Concrete byte layout for the combined container.

    The layout depends on metadata size because the JSON metadata block stores
    combined payload offsets. The combiner therefore rebuilds layout until the
    metadata size stabilizes.
    """

    alignment_bytes: int
    metadata_offset: int
    metadata_bytes: int
    slice_table_offset: int
    slice_table_bytes: int
    first_payload_offset: int
    payload_bytes_total: int
    total_file_bytes: int
    payload_offsets: Tuple[int, ...]


@dataclass(frozen=True)
class SliceEntry:
    """
    One slice-table entry in the combined container.

    Each entry stores the slice time, the copied payload location, and compact
    provenance facts from the source file.
    """

    output_index: int
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
class CombinedHeaderInfo:
    """Parsed fixed-header fields from a combined container."""

    magic: str
    fixed_header_bytes: int
    alignment_bytes: int
    metadata_offset: int
    metadata_bytes: int
    slice_table_offset: int
    slice_table_entry_bytes: int
    num_time_slices: int
    first_payload_offset: int
    payload_bytes_total: int
    total_file_bytes: int


@dataclass(frozen=True)
class SliceEntryInfo:
    """Parsed slice-table entry from a combined container."""

    output_index: int
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


def parse_args() -> argparse.Namespace:
    """
    Parse command-line arguments.

    :return: Parsed CLI namespace.
    """
    parser = argparse.ArgumentParser(
        description=("Combine raytracing time-slice files into one storage container.")
    )
    parser.add_argument("input_files", nargs="*", help="Explicit input slice files.")
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
        help="Alignment in bytes for the slice table and payload regions.",
    )
    parser.add_argument(
        "--run-metadata",
        default=None,
        help="""JSON file containing run metadata to embed verbatim. Required unless --inspect is used.""",
    )
    parser.add_argument(
        "--check-finite-payload",
        action="store_true",
        help=(
            "While copying payload bytes, reject NaN or Inf values interpreted "
            "as little-endian float64."
        ),
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


def parse_input_slice_file(path: Path) -> InputSliceInfo:
    """
    Parse one source raytracing slice file.

    The parser reads only the source-header fields needed for byte integrity,
    simulation-time sorting, and first-slice grid provenance, then seeks to
    ``header_size`` instead of walking the source string tables.

    :param path: Input binary file path.
    :return: Parsed source-slice metadata.
    :raises RuntimeError: If the file is truncated.
    """
    # Step 1: Read the source header in writer order, retaining only the
    #         fields needed later.
    actual_file_size = path.stat().st_size
    with path.open("rb") as fp:
        magic = _read_exact(fp, 16, "input magic").split(b"\0", 1)[0].decode("ascii")
        header_size = _read_u32(fp, "header_size")
        if header_size <= 0:
            raise RuntimeError(f"{path}: header_size must be positive.")
        if header_size > actual_file_size:
            raise RuntimeError(
                f"{path}: header_size={header_size} exceeded the actual file size."
            )
        output_index = _read_u32(fp, "output_index")
        _ = _read_u32(fp, "num_grids")

        _ = _read_u32(fp, "serialized_real_bytes")
        _read_u32(fp, "record_component_count")
        _read_u32(fp, "metric_component_count")
        _read_u32(fp, "christoffel_component_count")
        _ = _read_u32(fp, "point_record_real_count")
        point_record_bytes = _read_u32(fp, "point_record_bytes")
        _ = _read_u32(fp, "payload_includes_ghost_zones")
        _ = _read_u32(fp, "file_is_little_endian")
        _ = _read_u32(fp, "time_variable_is_f64")
        _ = _read_u32(fp, "reserved_u32")

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

        fp.seek(72, os.SEEK_CUR)

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

        # Step 2: Skip the remainder of the source header and land exactly on
        #         the advertised header boundary.
        if fp.tell() > header_size:
            raise RuntimeError(
                f"{path}: parsed input header length {fp.tell()} did not match "
                f"header_size={header_size}."
            )
        fp.seek(header_size)

        # Step 3: Read the physical simulation time from its recorded offset.
        fp.seek(simulation_time_offset)
        simulation_time = _read_f64(fp, "simulation_time")

    grid_metadata = GridMetadata(
        Nxx=Nxx,
        Nxx_plus_2NGHOSTS=Nxx_plus_2NGHOSTS,
        NGHOSTS=NGHOSTS,
        dxx=dxx,
        invdxx=invdxx,
        xxmin=xxmin,
        xxmax=xxmax,
        cart_origin=cart_origin,
    )
    return InputSliceInfo(
        path=path,
        magic=magic,
        header_size=header_size,
        output_index=output_index,
        point_record_count=point_record_count,
        point_record_bytes=point_record_bytes,
        simulation_time_offset=simulation_time_offset,
        point_records_offset=point_records_offset,
        point_records_bytes=point_records_bytes,
        total_file_bytes=total_file_bytes,
        actual_file_size=actual_file_size,
        simulation_time=simulation_time,
        grid_metadata=grid_metadata,
    )


def validate_input_slice_binary_integrity(info: InputSliceInfo) -> None:
    """
    Validate source-slice byte integrity only.

    :param info: Parsed source-slice metadata.
    :raises RuntimeError: If the file fails binary-integrity checks.
    """
    if info.magic != INPUT_SLICE_MAGIC:
        raise RuntimeError(f"{info.path}: unexpected magic '{info.magic}'.")
    if info.header_size <= 0:
        raise RuntimeError(f"{info.path}: header_size must be positive.")
    if info.header_size > info.actual_file_size:
        raise RuntimeError(
            f"{info.path}: header_size={info.header_size} exceeded the actual file size."
        )
    if info.simulation_time_offset + 8 > info.actual_file_size:
        raise RuntimeError(
            f"{info.path}: simulation_time_offset exceeded the actual file size."
        )
    if info.point_records_offset > info.actual_file_size:
        raise RuntimeError(
            f"{info.path}: point_records_offset exceeded the actual file size."
        )
    if info.point_record_count <= 0:
        raise RuntimeError(f"{info.path}: point_record_count must be positive.")
    if info.point_record_bytes <= 0:
        raise RuntimeError(f"{info.path}: point_record_bytes must be positive.")
    if info.point_records_bytes <= 0:
        raise RuntimeError(f"{info.path}: point_records_bytes must be positive.")
    expected_point_records_bytes = info.point_record_count * info.point_record_bytes
    if info.point_records_bytes != expected_point_records_bytes:
        raise RuntimeError(
            f"{info.path}: point_records_bytes={info.point_records_bytes} did not "
            f"equal point_record_count*point_record_bytes="
            f"{expected_point_records_bytes}."
        )
    if info.point_records_offset + info.point_records_bytes != info.total_file_bytes:
        raise RuntimeError(
            f"{info.path}: point_records_offset + point_records_bytes did not "
            "equal total_file_bytes."
        )
    if info.actual_file_size != info.total_file_bytes:
        raise RuntimeError(
            f"{info.path}: actual size {info.actual_file_size} did not match "
            f"total_file_bytes={info.total_file_bytes}."
        )
    if not math.isfinite(info.simulation_time):
        raise RuntimeError(f"{info.path}: simulation_time must be finite.")


def compute_layout(
    infos: Sequence[InputSliceInfo],
    metadata_bytes_len: int,
    alignment_bytes: int,
) -> Layout:
    """
    Compute the combined container layout.

    :param infos: Sorted input-slice metadata.
    :param metadata_bytes_len: Metadata JSON byte count.
    :param alignment_bytes: Requested block alignment.
    :return: Concrete layout offsets and sizes.
    """
    # Step 1: Align the slice table and each payload start independently.
    _validate_power_of_two(alignment_bytes, "alignment")
    metadata_offset = FIXED_HEADER_BYTES
    slice_table_offset = align_up(metadata_offset + metadata_bytes_len, alignment_bytes)
    slice_table_bytes = len(infos) * SLICE_TABLE_ENTRY_BYTES
    first_payload_offset = align_up(
        slice_table_offset + slice_table_bytes, alignment_bytes
    )

    payload_offsets: List[int] = []
    cursor = first_payload_offset
    for info in infos:
        cursor = align_up(cursor, alignment_bytes)
        payload_offsets.append(cursor)
        cursor += info.point_records_bytes

    total_file_bytes = cursor
    payload_bytes_total = total_file_bytes - first_payload_offset
    return Layout(
        alignment_bytes=alignment_bytes,
        metadata_offset=metadata_offset,
        metadata_bytes=metadata_bytes_len,
        slice_table_offset=slice_table_offset,
        slice_table_bytes=slice_table_bytes,
        first_payload_offset=first_payload_offset,
        payload_bytes_total=payload_bytes_total,
        total_file_bytes=total_file_bytes,
        payload_offsets=tuple(payload_offsets),
    )


def build_slice_entries(
    infos: Sequence[InputSliceInfo], layout: Layout
) -> List[SliceEntry]:
    """
    Build the combined slice table.

    :param infos: Sorted input-slice metadata.
    :param layout: Final combined container layout.
    :return: Slice-table entries.
    """
    # Step 1: Record per-slice offsets and source provenance after sorting.
    entries: List[SliceEntry] = []
    for slice_index, info in enumerate(infos):
        entries.append(
            SliceEntry(
                output_index=info.output_index,
                simulation_time=info.simulation_time,
                payload_offset=layout.payload_offsets[slice_index],
                payload_bytes=info.point_records_bytes,
                point_record_count=info.point_record_count,
                point_record_bytes=info.point_record_bytes,
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
    infos: Sequence[InputSliceInfo],
    slice_entries: Sequence[SliceEntry],
    layout: Layout,
    run_metadata: Dict[str, object],
) -> bytes:
    """
    Build deterministic JSON metadata for the combined container.

    :param infos: Sorted input-slice metadata.
    :param slice_entries: Slice-table entries matching ``infos``.
    :param layout: Final combined container layout.
    :param run_metadata: Embedded run-metadata object.
    :return: UTF-8 encoded deterministic JSON metadata.
    """
    # Step 1: Emit only storage, first-slice grid, and provenance metadata.
    grid_metadata = infos[0].grid_metadata
    metadata = {
        "container_role": "storage_combiner",
        "payload_policy": "copied_verbatim_from_input_slices",
        "sort_key": "simulation_time",
        "first_slice_grid_metadata": {
            "source_file": os.path.relpath(infos[0].path, Path.cwd()),
            "Nxx": list(grid_metadata.Nxx),
            "Nxx_plus_2NGHOSTS": list(grid_metadata.Nxx_plus_2NGHOSTS),
            "NGHOSTS": grid_metadata.NGHOSTS,
            "dxx": list(grid_metadata.dxx),
            "invdxx": list(grid_metadata.invdxx),
            "xxmin": list(grid_metadata.xxmin),
            "xxmax": list(grid_metadata.xxmax),
            "cart_origin": list(grid_metadata.cart_origin),
        },
        "run_metadata": run_metadata,
        "layout_summary": {
            "alignment_bytes": layout.alignment_bytes,
            "metadata_offset": layout.metadata_offset,
            "slice_table_offset": layout.slice_table_offset,
            "slice_table_entry_bytes": SLICE_TABLE_ENTRY_BYTES,
            "num_time_slices": len(infos),
            "first_payload_offset": layout.first_payload_offset,
            "payload_bytes_total": layout.payload_bytes_total,
            "total_file_bytes": layout.total_file_bytes,
        },
        "slices": [
            {
                "source_file": os.path.relpath(info.path, Path.cwd()),
                "output_index": info.output_index,
                "simulation_time": info.simulation_time,
                "source_file_size": info.actual_file_size,
                "source_header_size": info.header_size,
                "source_simulation_time_offset": info.simulation_time_offset,
                "source_point_records_offset": info.point_records_offset,
                "source_point_records_bytes": info.point_records_bytes,
                "combined_payload_offset": entry.payload_offset,
                "combined_payload_bytes": entry.payload_bytes,
            }
            for info, entry in zip(infos, slice_entries)
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
    Pack one slice-table entry.

    :param entry: Slice-table entry values.
    :return: Fixed-width packed entry bytes.
    :raises RuntimeError: If the packed entry size differs from the schema.
    """
    packed = struct.pack(
        SLICE_TABLE_ENTRY_FORMAT,
        _require_u64(entry.output_index, "slice_entry.output_index"),
        entry.simulation_time,
        _require_u64(entry.payload_offset, "slice_entry.payload_offset"),
        _require_u64(entry.payload_bytes, "slice_entry.payload_bytes"),
        _require_u64(entry.point_record_count, "slice_entry.point_record_count"),
        _require_u64(entry.point_record_bytes, "slice_entry.point_record_bytes"),
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
            entry.source_total_file_bytes,
            "slice_entry.source_total_file_bytes",
        ),
    )
    if len(packed) != SLICE_TABLE_ENTRY_BYTES:
        raise RuntimeError(
            f"slice table entry packed to {len(packed)} bytes, expected "
            f"{SLICE_TABLE_ENTRY_BYTES}."
        )
    return packed


def pack_fixed_header(layout: Layout, num_time_slices: int) -> bytes:
    """
    Pack the fixed combined-container header.

    :param layout: Final combined layout.
    :param num_time_slices: Number of sorted slices in the container.
    :return: ``FIXED_HEADER_BYTES`` header bytes.
    :raises RuntimeError: If the packed header exceeds its fixed allocation.
    """
    # Step 1: Pack only the navigation fields needed to find later blocks.
    magic_bytes = COMBINED_MAGIC.encode("ascii")
    if len(magic_bytes) >= 16:
        raise RuntimeError("Combined magic did not fit in the fixed header.")
    header = struct.pack(
        FIXED_HEADER_FORMAT,
        magic_bytes + b"\0" * (16 - len(magic_bytes)),
        _require_u64(FIXED_HEADER_BYTES, "fixed_header_bytes"),
        _require_u64(layout.alignment_bytes, "alignment_bytes"),
        _require_u64(layout.metadata_offset, "metadata_offset"),
        _require_u64(layout.metadata_bytes, "metadata_bytes"),
        _require_u64(layout.slice_table_offset, "slice_table_offset"),
        _require_u64(SLICE_TABLE_ENTRY_BYTES, "slice_table_entry_bytes"),
        _require_u64(num_time_slices, "num_time_slices"),
        _require_u64(layout.first_payload_offset, "first_payload_offset"),
        _require_u64(layout.payload_bytes_total, "payload_bytes_total"),
        _require_u64(layout.total_file_bytes, "total_file_bytes"),
    )
    if len(header) > FIXED_HEADER_BYTES:
        raise RuntimeError("Fixed header overflowed its fixed allocation.")
    return header + b"\0" * (FIXED_HEADER_BYTES - len(header))


def copy_payload(
    src_path: Path,
    src_offset: int,
    num_bytes: int,
    dst_file: BinaryIO,
    dst_offset: int,
    check_finite_payload: bool,
) -> None:
    """
    Copy one source payload region into the combined output.

    :param src_path: Source input-slice path.
    :param src_offset: Source payload offset.
    :param num_bytes: Number of bytes to copy.
    :param dst_file: Open destination file.
    :param dst_offset: Destination payload offset.
    :param check_finite_payload: Whether to validate float64 finiteness while copying.
    :raises RuntimeError: If EOF is encountered while copying.
    """
    # Step 1: Stream the payload directly from source to destination.
    chunk_size = 16 * 1024 * 1024
    if check_finite_payload and num_bytes % 8 != 0:
        raise RuntimeError(
            f"{src_path}: payload byte count {num_bytes} was not divisible by 8."
        )

    with src_path.open("rb") as src:
        src.seek(src_offset)
        dst_file.seek(dst_offset)
        remaining = num_bytes
        copied_value_count = 0
        while remaining > 0:
            chunk = src.read(min(chunk_size, remaining))
            if len(chunk) == 0:
                raise RuntimeError(
                    f"Unexpected EOF while copying payload from {src_path}."
                )
            if check_finite_payload:
                if len(chunk) % 8 != 0:
                    raise RuntimeError(
                        f"{src_path}: copied payload chunk size {len(chunk)} was not "
                        "divisible by 8."
                    )
                for local_index, unpacked in enumerate(struct.iter_unpack("<d", chunk)):
                    if not math.isfinite(unpacked[0]):
                        raise RuntimeError(
                            f"{src_path}: payload contained a non-finite float64 "
                            f"value at index {copied_value_count + local_index}."
                        )
                copied_value_count += len(chunk) // 8
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
    infos: Sequence[InputSliceInfo],
    layout: Layout,
    metadata_bytes: bytes,
    slice_entries: Sequence[SliceEntry],
    force: bool,
    check_finite_payload: bool,
) -> None:
    """
    Write the combined container atomically.

    :param output_path: Final destination path.
    :param infos: Sorted input metadata.
    :param layout: Final layout.
    :param metadata_bytes: Metadata JSON bytes.
    :param slice_entries: Slice-table entries.
    :param force: Whether overwriting an existing file is allowed.
    :param check_finite_payload: Whether to validate payload finiteness while copying.
    :raises RuntimeError: If overwrite rules fail during the final install step.
    :raises Exception: Propagates I/O failures after best-effort temporary cleanup.
    """
    # Step 1: Write into a same-directory temporary file so the final install
    #         step stays atomic on the destination filesystem.
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_fd, tmp_name = tempfile.mkstemp(
        prefix=f"{output_path.name}.tmp.",
        dir=str(output_path.parent),
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(tmp_fd, "wb") as out:
            # Step 2: Pre-size the file so aligned gaps are already zero-filled.
            out.truncate(layout.total_file_bytes)
            out.seek(0)
            out.write(pack_fixed_header(layout, len(infos)))
            out.seek(layout.metadata_offset)
            out.write(metadata_bytes)
            out.seek(layout.slice_table_offset)
            for entry in slice_entries:
                out.write(pack_slice_entry(entry))

            # Step 3: Copy each sorted source payload byte-for-byte.
            for slice_index, info in enumerate(infos):
                copy_payload(
                    src_path=info.path,
                    src_offset=info.point_records_offset,
                    num_bytes=info.point_records_bytes,
                    dst_file=out,
                    dst_offset=slice_entries[slice_index].payload_offset,
                    check_finite_payload=check_finite_payload,
                )
            out.flush()
            os.fsync(out.fileno())

        current_mode = tmp_path.stat().st_mode & 0o777
        os.chmod(tmp_path, current_mode & ~0o222)

        # Step 4: Install the completed file only after data and metadata are durable.
        if force:
            os.replace(tmp_path, output_path)
        else:
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
    Unpack the combined fixed header.

    :param header_bytes: First ``FIXED_HEADER_BYTES`` bytes of the combined file.
    :return: Parsed fixed-header fields.
    """
    unpacked = cast(
        Tuple[bytes, int, int, int, int, int, int, int, int, int, int],
        struct.unpack_from(FIXED_HEADER_FORMAT, header_bytes, 0),
    )
    return CombinedHeaderInfo(
        magic=unpacked[0].split(b"\0", 1)[0].decode("ascii"),
        fixed_header_bytes=unpacked[1],
        alignment_bytes=unpacked[2],
        metadata_offset=unpacked[3],
        metadata_bytes=unpacked[4],
        slice_table_offset=unpacked[5],
        slice_table_entry_bytes=unpacked[6],
        num_time_slices=unpacked[7],
        first_payload_offset=unpacked[8],
        payload_bytes_total=unpacked[9],
        total_file_bytes=unpacked[10],
    )


def _unpack_slice_entry(raw_entry: bytes) -> SliceEntryInfo:
    """
    Unpack one fixed-width slice-table entry.

    :param raw_entry: One slice-table record.
    :return: Parsed slice-table entry fields.
    """
    unpacked = cast(
        Tuple[int, float, int, int, int, int, int, int, int, int, int, int],
        struct.unpack(SLICE_TABLE_ENTRY_FORMAT, raw_entry),
    )
    return SliceEntryInfo(
        output_index=unpacked[0],
        simulation_time=unpacked[1],
        payload_offset=unpacked[2],
        payload_bytes=unpacked[3],
        point_record_count=unpacked[4],
        point_record_bytes=unpacked[5],
        source_file_size=unpacked[6],
        source_header_size=unpacked[7],
        source_simulation_time_offset=unpacked[8],
        source_point_records_offset=unpacked[9],
        source_point_records_bytes=unpacked[10],
        source_total_file_bytes=unpacked[11],
    )


def inspect_combined_file(path: Path) -> None:
    """
    Inspect a combined file and print a concise summary.

    :param path: Combined container path.
    :raises RuntimeError: If validation fails.
    """
    # Step 1: Parse the fixed header, metadata block, and slice table.
    file_size = path.stat().st_size
    with path.open("rb") as fp:
        header_bytes = _read_exact(fp, FIXED_HEADER_BYTES, "fixed header")
        header = _unpack_fixed_header(header_bytes)
        if header.magic != COMBINED_MAGIC:
            raise RuntimeError(f"{path}: combined magic was invalid.")
        if header.fixed_header_bytes != FIXED_HEADER_BYTES:
            raise RuntimeError(f"{path}: fixed_header_bytes was invalid.")
        _validate_power_of_two(header.alignment_bytes, "header.alignment_bytes")
        if header.metadata_offset != FIXED_HEADER_BYTES:
            raise RuntimeError(f"{path}: metadata_offset was invalid.")
        if header.slice_table_entry_bytes != SLICE_TABLE_ENTRY_BYTES:
            raise RuntimeError(f"{path}: slice_table_entry_bytes was invalid.")
        if file_size != header.total_file_bytes:
            raise RuntimeError(
                f"{path}: file size {file_size} did not match total_file_bytes="
                f"{header.total_file_bytes}."
            )
        if header.metadata_offset + header.metadata_bytes > file_size:
            raise RuntimeError(f"{path}: metadata block exceeded file bounds.")
        slice_table_bytes = header.num_time_slices * header.slice_table_entry_bytes
        if header.slice_table_offset + slice_table_bytes > file_size:
            raise RuntimeError(f"{path}: slice table exceeded file bounds.")
        if header.slice_table_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: slice_table_offset was not aligned.")
        if header.first_payload_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: first_payload_offset was not aligned.")
        if header.total_file_bytes - header.first_payload_offset != (
            header.payload_bytes_total
        ):
            raise RuntimeError(f"{path}: payload_bytes_total was invalid.")
        if header.metadata_offset + header.metadata_bytes > header.slice_table_offset:
            raise RuntimeError(f"{path}: metadata block overlapped the slice table.")
        if header.slice_table_offset + slice_table_bytes > header.first_payload_offset:
            raise RuntimeError(f"{path}: slice table overlapped the payload region.")

        fp.seek(header.metadata_offset)
        metadata = json.loads(
            _read_exact(fp, header.metadata_bytes, "metadata JSON").decode("utf-8")
        )
        if not isinstance(metadata, dict):
            raise RuntimeError(f"{path}: metadata JSON was not an object.")

        fp.seek(header.slice_table_offset)
        slice_entries: List[SliceEntryInfo] = []
        for _ in range(header.num_time_slices):
            raw_entry = _read_exact(fp, SLICE_TABLE_ENTRY_BYTES, "slice table entry")
            slice_entries.append(_unpack_slice_entry(raw_entry))

    # Step 2: Validate per-slice bounds and strict time ordering.
    previous_time: Optional[float] = None
    previous_payload_end = header.first_payload_offset
    for entry in slice_entries:
        if not math.isfinite(entry.simulation_time):
            raise RuntimeError(f"{path}: slice simulation_time was not finite.")
        if previous_time is not None and entry.simulation_time <= previous_time:
            raise RuntimeError(f"{path}: slice times were not strictly increasing.")
        previous_time = entry.simulation_time

        if entry.payload_offset % header.alignment_bytes != 0:
            raise RuntimeError(f"{path}: slice payload offset was not aligned.")
        if entry.payload_bytes <= 0:
            raise RuntimeError(f"{path}: slice payload_bytes must be positive.")
        if entry.point_record_count <= 0:
            raise RuntimeError(f"{path}: slice point_record_count must be positive.")
        if entry.point_record_bytes <= 0:
            raise RuntimeError(f"{path}: slice point_record_bytes must be positive.")
        expected_payload_bytes = entry.point_record_count * entry.point_record_bytes
        if entry.payload_bytes != expected_payload_bytes:
            raise RuntimeError(
                f"{path}: slice payload_bytes did not match "
                "point_record_count*point_record_bytes."
            )
        if entry.payload_offset < header.first_payload_offset:
            raise RuntimeError(
                f"{path}: slice payload offset started before first_payload_offset."
            )
        if entry.payload_offset < previous_payload_end:
            raise RuntimeError(f"{path}: slice payloads overlapped.")
        if entry.payload_offset + entry.payload_bytes > file_size:
            raise RuntimeError(f"{path}: slice payload exceeded file bounds.")
        previous_payload_end = entry.payload_offset + entry.payload_bytes

    # Step 3: Print a concise human-readable summary for quick inspection.
    times = [entry.simulation_time for entry in slice_entries]
    run_metadata = metadata.get("run_metadata", {})
    first_slice_grid_metadata = metadata.get("first_slice_grid_metadata", {})

    print(f"magic: {header.magic}")
    print(f"total_file_bytes: {header.total_file_bytes}")
    print(f"num_time_slices: {header.num_time_slices}")
    print(f"metadata_bytes: {header.metadata_bytes}")
    print(f"slice_table_offset: {header.slice_table_offset}")
    print(f"first_payload_offset: {header.first_payload_offset}")
    if times:
        print("time range:")
        print(f"  first simulation_time: {times[0]:.17g}")
        print(f"  last simulation_time: {times[-1]:.17g}")
    if isinstance(first_slice_grid_metadata, dict) and first_slice_grid_metadata:
        print("first_slice_grid_metadata:")
        for key in (
            "Nxx",
            "Nxx_plus_2NGHOSTS",
            "NGHOSTS",
            "dxx",
            "xxmin",
            "xxmax",
        ):
            if key in first_slice_grid_metadata:
                print(f"  {key}: {first_slice_grid_metadata[key]}")
    print("run_metadata summary:")
    if isinstance(run_metadata, dict) and run_metadata:
        for key in (
            "generator_script",
            "project_name",
            "raytracing_coord_system",
            "raytracing_Nxx",
            "raytracing_time",
            "slice_filename_pattern",
        ):
            if key in run_metadata:
                print(f"  {key}: {run_metadata[key]}")
    else:
        print("  {}")

    head_entries = slice_entries[: min(3, len(slice_entries))]
    tail_entries = slice_entries[-min(3, len(slice_entries)) :]
    print("first few slice table entries:")
    for entry in head_entries:
        print(
            "  "
            f"output_index={entry.output_index} "
            f"simulation_time={entry.simulation_time:.17g} "
            f"payload_offset={entry.payload_offset} "
            f"payload_bytes={entry.payload_bytes}"
        )
    print("last few slice table entries:")
    for entry in tail_entries:
        print(
            "  "
            f"output_index={entry.output_index} "
            f"simulation_time={entry.simulation_time:.17g} "
            f"payload_offset={entry.payload_offset} "
            f"payload_bytes={entry.payload_bytes}"
        )


def main() -> None:
    """
    Drive CLI execution.

    :return: None.
    :raises RuntimeError: If discovery, validation, layout, or output writing fails.
    """
    # Step 1: Dispatch inspect mode before any input discovery.
    args = parse_args()
    if args.inspect:
        inspect_combined_file(Path(args.inspect))
        return

    # Step 2: Discover inputs and enforce output-path safety.
    input_paths = discover_input_paths(args)
    if not input_paths:
        raise RuntimeError("No input raytracing files found.")

    output_path = Path(args.output).expanduser().resolve()
    if output_path in input_paths:
        raise RuntimeError("Output path must not also be one of the input files.")
    if output_path.exists() and not args.force:
        raise RuntimeError(
            f"Output file {output_path} exists. Use --force to overwrite."
        )
    if args.run_metadata is None:
        raise RuntimeError("--run-metadata is required when combining files.")
    run_metadata_path = Path(args.run_metadata).expanduser().resolve()
    with run_metadata_path.open("r", encoding="utf-8") as fp:
        run_metadata_obj = json.load(fp)
    if not isinstance(run_metadata_obj, dict):
        raise RuntimeError(f"{run_metadata_path}: run metadata JSON must be an object.")
    run_metadata = cast(Dict[str, object], run_metadata_obj)

    # Step 3: Parse and validate each source slice independently.
    infos = [parse_input_slice_file(path) for path in input_paths]
    for info in infos:
        validate_input_slice_binary_integrity(info)
    infos.sort(key=lambda info: info.simulation_time)
    previous_time: Optional[float] = None
    for info in infos:
        if not math.isfinite(info.simulation_time):
            raise RuntimeError(f"{info.path}: simulation_time must be finite.")
        if previous_time is not None and info.simulation_time <= previous_time:
            raise RuntimeError("simulation_time values must be strictly increasing.")
        previous_time = info.simulation_time

    # Step 4: Iterate until metadata size and payload offsets stabilize.
    layout = compute_layout(infos, metadata_bytes_len=0, alignment_bytes=args.alignment)
    metadata_bytes = b""
    slice_entries: List[SliceEntry] = []
    for _ in range(8):
        slice_entries = build_slice_entries(infos, layout)
        metadata_bytes = build_metadata_json(
            infos=infos,
            slice_entries=slice_entries,
            layout=layout,
            run_metadata=run_metadata,
        )
        next_layout = compute_layout(
            infos,
            metadata_bytes_len=len(metadata_bytes),
            alignment_bytes=args.alignment,
        )
        if next_layout == layout:
            break
        layout = next_layout
    else:
        raise RuntimeError("Combined layout did not stabilize after metadata rebuild.")

    final_layout = layout
    slice_entries = build_slice_entries(infos, final_layout)
    metadata_bytes = build_metadata_json(
        infos=infos,
        slice_entries=slice_entries,
        layout=final_layout,
        run_metadata=run_metadata,
    )
    if final_layout.metadata_bytes != len(metadata_bytes):
        raise RuntimeError(
            "Final layout metadata_bytes did not match the serialized metadata size."
        )

    # Step 5: Atomically materialize the combined container.
    write_combined_file_atomically(
        output_path=output_path,
        infos=infos,
        layout=final_layout,
        metadata_bytes=metadata_bytes,
        slice_entries=slice_entries,
        force=args.force,
        check_finite_payload=args.check_finite_payload,
    )


if __name__ == "__main__":
    import doctest
    import sys

    if len(sys.argv) == 1:
        results = doctest.testmod()

        if results.failed > 0:
            print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
            sys.exit(1)
        else:
            print(f"Doctest passed: All {results.attempted} test(s) passed")
    else:
        main()

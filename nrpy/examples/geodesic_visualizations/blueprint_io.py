r"""
Read and validate native geodesic blueprint binary files.

The reader is shared by the renderer and diagnostics so the native record layout,
tile identity, record counts, and terminal-status rules have one owner.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com

Doctests:
>>> import tempfile
>>> from pathlib import Path
>>> with tempfile.TemporaryDirectory() as directory:
...     path = Path(directory) / "light_blueprint_01_00.bin"
...     raw = struct.pack(cfg.BLUEPRINT_HEADER_FORMAT, cfg.BLUEPRINT_MAGIC,
...                       cfg.BLUEPRINT_SCHEMA_VERSION, cfg.BLUEPRINT_HEADER_SIZE,
...                       cfg.BLUEPRINT_RECORD_SIZE, 1, 0, 2, 1, 4,
...                       1.0, 1.0, 0.0, 1.0, 0.0, 1.0)
...     _ = path.write_bytes(raw + bytes(4 * cfg.BLUEPRINT_RECORD_SIZE))
...     puncture_path = Path(directory) / "light_blueprint_puncture_01_00.bin"
...     nan = float("nan")
...     valid = b"".join(struct.pack("=BBd", *record) for record in
...                      [(0, 2, nan), (1, 1, -0.125), (2, 0, 0.25), (3, 2, nan)])
...     _ = puncture_path.write_bytes(valid)
...     assert read_puncture_labels(str(path)).tolist() == [0, 1, 2, 3]
...     records = read_puncture_records(str(path))
...     assert records["radial_gap"][1:3].tolist() == [-0.125, 0.25]
...     assert np.isnan(records["radial_gap"][[0, 3]]).all()
...     assert not records.flags.writeable
...     del records
...     for invalid in (bytes([0, 1, 2, 3]), bytes(8), valid[:-1], valid + b"\0",
...                     valid[:30] + struct.pack("=BBd", 4, 2, nan),
...                     valid[:30] + struct.pack("=BBd", 3, 3, nan),
...                     valid[:30] + struct.pack("=BBd", 3, 2, float("inf"))):
...         _ = puncture_path.write_bytes(invalid)
...         try:
...             read_puncture_labels(str(path))
...         except ValueError:
...             pass
...         else:
...             raise AssertionError("Invalid puncture labels were accepted")
...     puncture_path.unlink()
...     try:
...         read_puncture_labels(str(path))
...     except FileNotFoundError:
...         pass
...     else:
...         raise AssertionError("Missing puncture labels were accepted")
>>> with tempfile.TemporaryDirectory() as directory:
...     path = Path(directory) / "light_blueprint_00_00.bin"
...     record = np.zeros(1, dtype=cfg.BLUEPRINT_DTYPE)
...     record["termination_type"] = 12
...     record["y_t"], record["z_t"] = 3, np.nan
...     def check_tile(version, horizon_id):
...         record["y_t"] = horizon_id
...         raw = struct.pack("=8s7IQ6d", b"NRPYBP01", version, 92, 124,
...                           0, 0, 1, 1, 1, 1., 1., 0., 1., 0., 1.)
...         _ = path.write_bytes(raw + record.tobytes())
...         for header, _, chunk in iter_blueprint_chunks(str(path), 1):
...             assert header.schema_version == version
...             assert chunk["y_t"].tolist() == [horizon_id]
...     check_tile(10, 3)
...     for version, horizon_id in ((9, 3), (10, 1.5), (10, 4)):
...         try:
...             check_tile(version, horizon_id)
...         except ValueError:
...             pass
...         else:
...             raise AssertionError("Invalid horizon-stop record was accepted")
"""

import os
import struct
from typing import Iterator, NamedTuple, Optional, Tuple

import numpy as np
import numpy.typing as npt

try:
    import blueprint_config_and_schema as cfg  # type: ignore
except ImportError:
    from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg


class BlueprintHeader(NamedTuple):
    """Tile identity/counts, camera fields of view, and global image-plane bounds."""

    magic: bytes
    schema_version: int
    header_size: int
    record_size: int
    tile_x: int
    tile_y: int
    tiles_width: int
    tiles_height: int
    record_count: int
    alpha_w: float
    alpha_h: float
    a_min: float
    a_max: float
    b_min: float
    b_max: float


def read_blueprint_header(
    filename: str, expected_tile: Optional[Tuple[int, int]] = None
) -> BlueprintHeader:
    """
    Read and validate a blueprint header and exact payload length.

    :param filename: Native blueprint binary-file path.
    :param expected_tile: Optional expected ``(tile_x, tile_y)`` pair.
    :return: Validated blueprint metadata.
    :raises ValueError: If the binary file is malformed or has the wrong tile.

    Doctests:
    >>> import tempfile
    >>> from pathlib import Path
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     path = Path(directory) / "light_blueprint_00_00.bin"
    ...     raw = struct.pack("=8s7IQ6d", b"NRPYBP01", 9, 92, 124,
    ...                       0, 0, 1, 1, 1, 1.2, 0.8, 0.25, 0.75, 0.625, 0.875)
    ...     _ = path.write_bytes(raw + bytes(124))
    ...     header = read_blueprint_header(str(path), expected_tile=(0, 0))
    ...     assert (header.a_min, header.a_max, header.b_min, header.b_max) == (0.25, 0.75, 0.625, 0.875)
    ...     assert header.header_size == 92 and header.record_count == 1
    """
    with open(filename, "rb") as blueprint_file:
        raw_header = blueprint_file.read(cfg.BLUEPRINT_HEADER_SIZE)
    if len(raw_header) != cfg.BLUEPRINT_HEADER_SIZE:
        raise ValueError(f"Blueprint '{filename}' has a truncated header")
    values = struct.unpack(cfg.BLUEPRINT_HEADER_FORMAT, raw_header)
    header = BlueprintHeader(*values)
    if header.magic != cfg.BLUEPRINT_MAGIC:
        raise ValueError(f"Blueprint '{filename}' has an invalid magic value")
    if header.schema_version not in (
        cfg.BLUEPRINT_SCHEMA_VERSION,
        cfg.BLUEPRINT_HORIZON_VERSION,
    ):
        raise ValueError(
            f"Blueprint '{filename}' uses version {header.schema_version}; "
            f"expected {cfg.BLUEPRINT_SCHEMA_VERSION} or {cfg.BLUEPRINT_HORIZON_VERSION}."
        )
    if header.header_size != cfg.BLUEPRINT_HEADER_SIZE:
        raise ValueError(f"Blueprint '{filename}' has an invalid header size")
    if header.record_size != cfg.BLUEPRINT_DTYPE.itemsize:
        raise ValueError(f"Blueprint '{filename}' has an invalid record size")
    if header.tiles_width <= 0 or header.tiles_height <= 0:
        raise ValueError(f"Blueprint '{filename}' has invalid tile dimensions")
    if header.tile_x < 0 or header.tile_x >= header.tiles_width:
        raise ValueError(f"Blueprint '{filename}' has an invalid tile x index")
    if header.tile_y < 0 or header.tile_y >= header.tiles_height:
        raise ValueError(f"Blueprint '{filename}' has an invalid tile y index")
    if not np.isfinite(header.alpha_w) or not np.isfinite(header.alpha_h):
        raise ValueError(f"Blueprint '{filename}' has non-finite fields of view")
    if not (0.0 < header.alpha_w < np.pi and 0.0 < header.alpha_h < np.pi):
        raise ValueError(f"Blueprint '{filename}' has invalid fields of view")
    if header.record_count <= 0:
        raise ValueError(f"Blueprint '{filename}' has an invalid record count")
    if expected_tile is not None and (header.tile_x, header.tile_y) != expected_tile:
        raise ValueError(f"Blueprint '{filename}' has an unexpected tile identity")
    expected_size = cfg.BLUEPRINT_HEADER_SIZE + header.record_count * header.record_size
    actual_size = os.path.getsize(filename)
    if actual_size != expected_size:
        raise ValueError(
            f"Blueprint '{filename}' has size {actual_size}; expected {expected_size}"
        )
    return header


def read_puncture_records(
    filename: str, header: Optional[BlueprintHeader] = None
) -> npt.NDArray[np.void]:
    """
    Validate and map puncture IDs, horizon statuses, and signed radial gaps.

    :param filename: Main blueprint tile path, used to locate its label file.
    :param header: Previously validated main header, or None to read it here.
    :return: Read-only mapped records in main blueprint record order.
    :raises FileNotFoundError: If the matching puncture label file is missing.
    :raises ValueError: If the record count, ID, status, or radial gap is invalid.
    """
    if header is None:
        header = read_blueprint_header(filename)
    puncture_path = os.path.join(
        os.path.dirname(filename),
        cfg.PUNCTURE_FILENAME_TEMPLATE.format(
            tile_x=header.tile_x, tile_y=header.tile_y
        ),
    )
    expected_size = header.record_count * cfg.PUNCTURE_RECORD_DTYPE.itemsize
    try:
        actual_size = os.path.getsize(puncture_path)
    except FileNotFoundError as exc:
        raise FileNotFoundError(
            f"Missing puncture blueprint '{puncture_path}'. Generate photon labels "
            "using --puncture-diagnostics-dir before requesting puncture coloring."
        ) from exc
    if actual_size != expected_size:
        raise ValueError(
            f"Puncture blueprint '{puncture_path}' has size {actual_size}; "
            f"expected {expected_size} bytes, ten bytes per main blueprint record. "
            "Regenerate older puncture blueprints using --puncture-diagnostics-dir."
        )
    records: npt.NDArray[np.void] = np.memmap(
        puncture_path, dtype=cfg.PUNCTURE_RECORD_DTYPE, mode="r"
    )
    labels: npt.NDArray[np.uint8] = records["puncture_id"]
    for start in range(0, header.record_count, cfg.CHUNK_SIZE):
        if np.any(labels[start : start + cfg.CHUNK_SIZE] > cfg.MERGED_PUNCTURE):
            raise ValueError(
                f"Puncture blueprint '{puncture_path}' has an invalid label"
            )
        if np.any(
            records["horizon_status"][start : start + cfg.CHUNK_SIZE]
            > cfg.HORIZON_NOT_CHECKED
        ):
            raise ValueError(
                f"Puncture blueprint '{puncture_path}' has an invalid horizon status"
            )
        if np.any(np.isinf(records["radial_gap"][start : start + cfg.CHUNK_SIZE])):
            raise ValueError(
                f"Puncture blueprint '{puncture_path}' has an infinite radial gap"
            )
    return records


def read_puncture_labels(
    filename: str, header: Optional[BlueprintHeader] = None
) -> npt.NDArray[np.uint8]:
    """
    Read validated puncture IDs for the existing image-coloring functions.

    Missing files and invalid records raise the errors from read_puncture_records.

    :param filename: Main blueprint tile path, used to locate its puncture file.
    :param header: Previously validated main header, or None to read it here.
    :return: Read-only mapped labels in main blueprint record order.
    """
    labels: npt.NDArray[np.uint8] = read_puncture_records(filename, header)[
        "puncture_id"
    ]
    return labels


def iter_blueprint_chunks(
    filename: str,
    chunk_records: int,
    expected_tile: Optional[Tuple[int, int]] = None,
    allow_active_termination: bool = False,
) -> Iterator[Tuple[BlueprintHeader, int, npt.NDArray[np.void]]]:
    """
    Stream validated blueprint records with their tile-local ordinal.

    :param filename: Native blueprint binary-file path.
    :param chunk_records: Maximum records yielded per chunk.
    :param expected_tile: Optional expected ``(tile_x, tile_y)`` pair.
    :param allow_active_termination: Permit rays that still have active status.
    :yield: Header, starting ordinal, and structured record chunks.
    :raises ValueError: If the chunk size or binary file is invalid.
    """
    if chunk_records <= 0:
        raise ValueError("chunk_records must be positive")
    header = read_blueprint_header(filename, expected_tile)
    permitted_termination_types = cfg.FINAL_TERMINATION_TYPES
    if allow_active_termination:
        permitted_termination_types = (*permitted_termination_types, cfg.ACTIVE)
    with open(filename, "rb") as blueprint_file:
        blueprint_file.seek(cfg.BLUEPRINT_HEADER_SIZE)
        record_start = 0
        while record_start < header.record_count:
            records_to_read = min(chunk_records, header.record_count - record_start)
            expected_bytes = records_to_read * header.record_size
            raw_records = blueprint_file.read(expected_bytes)
            if len(raw_records) != expected_bytes:
                raise ValueError(f"Blueprint '{filename}' has a truncated payload")
            records = np.frombuffer(raw_records, dtype=cfg.BLUEPRINT_DTYPE)
            width_fraction = records["image_width_fraction"]
            height_fraction = records["image_height_fraction"]
            if np.any(
                ~np.isfinite(width_fraction)
                | ~np.isfinite(height_fraction)
                | (width_fraction < 0.0)
                | (width_fraction > 1.0)
                | (height_fraction < 0.0)
                | (height_fraction > 1.0)
            ):
                raise ValueError(
                    f"Blueprint '{filename}' has invalid normalized image coordinates"
                )
            termination_types = records["termination_type"]
            if np.any(~np.isin(termination_types, permitted_termination_types)):
                raise ValueError(
                    f"Blueprint '{filename}' has an invalid termination type"
                )
            horizon_stops = (
                termination_types == cfg.STOP_CONDITION_INSIDE_APPARENT_HORIZON
            )
            if np.any(horizon_stops):
                if header.schema_version != cfg.BLUEPRINT_HORIZON_VERSION:
                    raise ValueError(
                        f"Blueprint '{filename}' requires version 10 for horizon stops"
                    )
                if np.any(
                    ~np.isin(
                        records["y_t"][horizon_stops],
                        (cfg.PUNCTURE_1, cfg.PUNCTURE_2, cfg.MERGED_PUNCTURE),
                    )
                    | ~np.isnan(records["z_t"][horizon_stops])
                ):
                    raise ValueError(
                        f"Blueprint '{filename}' has an invalid horizon-stop record"
                    )
            yield (
                header,
                record_start,
                records,
            )
            record_start += records_to_read


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

"""
Read and validate native geodesic blueprint artifacts.

The reader is shared by the renderer and diagnostics so native artifact layout,
tile identity, record counts, and terminal-status rules have one owner.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
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
    """Validated tile identity/counts and fields of view from one header."""

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


def _unpack_header(raw_header: bytes, filename: str) -> BlueprintHeader:
    """
    Decode and validate one native blueprint header.

    :param raw_header: Raw header bytes read from the artifact.
    :param filename: Artifact name used in error messages.
    :return: Validated blueprint metadata.
    :raises ValueError: If header size, magic, version, or dimensions are invalid.
    """
    if len(raw_header) != cfg.BLUEPRINT_HEADER_SIZE:
        raise ValueError(f"Blueprint '{filename}' has a truncated header")
    values = struct.unpack(cfg.BLUEPRINT_HEADER_FORMAT, raw_header)
    header = BlueprintHeader(*values)
    if header.magic != cfg.BLUEPRINT_MAGIC:
        raise ValueError(f"Blueprint '{filename}' has an invalid magic value")
    if header.schema_version != cfg.BLUEPRINT_SCHEMA_VERSION:
        raise ValueError(f"Blueprint '{filename}' has an unsupported schema version")
    if header.header_size != cfg.BLUEPRINT_HEADER_SIZE:
        raise ValueError(f"Blueprint '{filename}' has an invalid header size")
    if header.record_size != cfg.BLUEPRINT_RECORD_SIZE:
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
    return header


def read_blueprint_header(
    filename: str, expected_tile: Optional[Tuple[int, int]] = None
) -> BlueprintHeader:
    """
    Read and validate a blueprint header and exact payload length.

    :param filename: Native blueprint artifact path.
    :param expected_tile: Optional expected ``(tile_x, tile_y)`` pair.
    :return: Validated blueprint metadata.
    :raises ValueError: If the artifact is malformed or has the wrong tile.
    """
    with open(filename, "rb") as blueprint_file:
        header = _unpack_header(
            blueprint_file.read(cfg.BLUEPRINT_HEADER_SIZE), filename
        )
    if expected_tile is not None and (header.tile_x, header.tile_y) != expected_tile:
        raise ValueError(f"Blueprint '{filename}' has an unexpected tile identity")
    expected_size = cfg.BLUEPRINT_HEADER_SIZE + (
        header.record_count * cfg.BLUEPRINT_RECORD_SIZE
    )
    actual_size = os.path.getsize(filename)
    if actual_size != expected_size:
        raise ValueError(
            f"Blueprint '{filename}' has size {actual_size}; expected {expected_size}"
        )
    return header


def iter_blueprint_chunks(
    filename: str,
    chunk_records: int,
    expected_tile: Optional[Tuple[int, int]] = None,
) -> Iterator[Tuple[BlueprintHeader, int, npt.NDArray[np.void]]]:
    """
    Stream validated blueprint records with their tile-local ordinal.

    :param filename: Native blueprint artifact path.
    :param chunk_records: Maximum records yielded per chunk.
    :param expected_tile: Optional expected ``(tile_x, tile_y)`` pair.
    :yield: Header, starting ordinal, and structured record chunks.
    :raises ValueError: If the chunk size or artifact is invalid.
    """
    if chunk_records <= 0:
        raise ValueError("chunk_records must be positive")
    header = read_blueprint_header(filename, expected_tile)
    with open(filename, "rb") as blueprint_file:
        blueprint_file.seek(cfg.BLUEPRINT_HEADER_SIZE)
        record_start = 0
        while record_start < header.record_count:
            records_to_read = min(chunk_records, header.record_count - record_start)
            raw_records = blueprint_file.read(
                records_to_read * cfg.BLUEPRINT_RECORD_SIZE
            )
            expected_bytes = records_to_read * cfg.BLUEPRINT_RECORD_SIZE
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
            if np.any(
                (termination_types < cfg.STOP_CONDITION_COORD_RADIUS_EXCEEDED)
                | (termination_types > cfg.FAILURE_GENERIC)
            ):
                raise ValueError(
                    f"Blueprint '{filename}' has an invalid termination type"
                )
            yield (
                header,
                record_start,
                records,
            )
            record_start += records_to_read

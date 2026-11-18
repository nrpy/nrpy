"""Store BHaHAHA angular radius grids at the saved spacetime output times."""

import argparse
import importlib
import json
import math
import os
import struct
import tempfile
from pathlib import Path
from typing import IO, Dict, List, Optional, Tuple

import numpy as np

# The generated evolution directory contains the existing spacetime combiner.
spacetime = importlib.import_module(
    "combine_raytracing_time_slices"
    if Path(__file__).with_name("combine_raytracing_time_slices.py").is_file()
    else "nrpy.infrastructures.BHaH.diagnostics.combine_raytracing_time_slices"
)
HEADER = struct.Struct("<8sQdQ")
SLICE = struct.Struct("<Qd3Q")
SURFACE = struct.Struct("<QII9d")
MAGIC = b"NRPYAH01"
HorizonRow = Tuple[int, Tuple[float, float, float]]
SliceTimes = List[Tuple[int, float]]


def _spacetime_times(path: Path, raw_directory: Path) -> SliceTimes:
    """
    Read output indices and times without loading metric arrays.

    :param path: Combined spacetime filename.
    :param raw_directory: Directory of the corresponding uncombined slices.
    :return: Output indices and times in strictly increasing time order.
    :raises ValueError: If the slice table or times are invalid.
    """
    result: SliceTimes = []
    if path.is_file():
        with path.open("rb") as stream:
            header = spacetime._unpack_fixed_header(
                stream.read(spacetime.FIXED_HEADER_BYTES)
            )
            if (
                header.magic not in spacetime.FORMAT_BY_COMBINED_MAGIC
                or header.slice_table_entry_bytes != spacetime.SLICE_TABLE_ENTRY_BYTES
                or header.slice_table_offset < spacetime.FIXED_HEADER_BYTES
                or header.slice_table_offset
                + header.num_time_slices * header.slice_table_entry_bytes
                > path.stat().st_size
            ):
                raise ValueError(f"{path}: invalid spacetime slice table")
            stream.seek(header.slice_table_offset)
            for _ in range(header.num_time_slices):
                entry = spacetime._unpack_slice_entry(
                    stream.read(header.slice_table_entry_bytes)
                )
                result.append((entry.output_index, entry.simulation_time))
    else:
        for source in raw_directory.glob("raytracing_data_t????????.bin"):
            info = spacetime.parse_input_slice_file(source)
            result.append((info.output_index, info.simulation_time))
        result.sort(key=lambda entry: entry[1])
    if (
        not result
        or len({index for index, _ in result}) != len(result)
        or any(not math.isfinite(time) for _, time in result)
        or any(a[1] >= b[1] for a, b in zip(result, result[1:]))
    ):
        raise ValueError(
            f"{path}: missing slices or invalid spacetime times or output indices"
        )
    return result


def _read_track(path: Path) -> Dict[str, HorizonRow]:
    """
    Read iterations and centroids keyed by the written diagnostic time.

    :param path: BHaHAHA diagnostic table.
    :return: Rows keyed by time at BHaHAHA's printed precision.
    :raises ValueError: If a diagnostic time or iteration is repeated.
    """
    rows: Dict[str, HorizonRow] = {}
    iterations = set()
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if not fields or fields[0].startswith("#"):
            continue
        iteration = int(fields[0])
        time = f"{float(fields[1]):.3f}"
        x, y, z = map(float, fields[2:5])
        if (
            time in rows
            or iteration in iterations
            or not all(map(math.isfinite, (float(time), x, y, z)))
        ):
            raise ValueError(f"{path}: invalid or repeated diagnostic row at t={time}")
        rows[time] = (iteration, (x, y, z))
        iterations.add(iteration)
    return rows


def _write_surface(stream: IO[bytes], path: Path, row: HorizonRow) -> None:
    r"""
    Recover the radial origin and write a cell-centered angular grid.

    :param stream: Output binary file at the start of the shape record.
    :param path: BHaHAHA Cartesian ring file.
    :param row: Simulation iteration and diagnostic centroid.
    :raises ValueError: If the rings do not describe a positive radial grid.

    Doctests:
    >>> from io import BytesIO
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     shape = Path(directory) / "h.t0000000.ah1.gp"
    ...     rings = []
    ...     for theta in (math.pi / 4, 3 * math.pi / 4):
    ...         ring = []
    ...         for j in (0, 1, 2, 3, 0):
    ...             phi = -math.pi + (j + 0.5) * math.pi / 2
    ...             x = 7 + 2 * math.sin(theta) * math.cos(phi)
    ...             y = 2 * math.sin(theta) * math.sin(phi)
    ...             z = 2 * math.cos(theta)
    ...             ring.append(f"{x:.17g} {y:.17g} {z:.17g}")
    ...         rings.append("\n".join(ring))
    ...     _ = shape.write_text("\n\n".join(rings), encoding="ascii")
    ...     stream = BytesIO()
    ...     _write_surface(stream, shape, (0, (6.0, 0.0, 0.0)))
    ...     data = stream.getvalue()
    >>> SURFACE.unpack_from(data)[:4]
    (0, 2, 4, 6.0)
    >>> np.allclose(SURFACE.unpack_from(data)[6:9], [7, 0, 0])
    True
    >>> np.allclose(struct.unpack_from("<8d", data, SURFACE.size), 2)
    True
    """
    rings: List[List[List[float]]] = []
    points: List[List[float]] = []
    for line in path.read_text(encoding="ascii").splitlines() + [""]:
        line = line.strip()
        if line.startswith("#"):
            continue
        if line:
            points.append([float(value) for value in line.split()])
        elif points:
            rings.append(points)
            points = []
    if (
        len(rings) < 2
        or len(rings[0]) < 3
        or (len(rings[0]) - 1) % 2
        or any(len(ring) != len(rings[0]) for ring in rings)
        or any(len(point) != 3 for ring in rings for point in ring)
    ):
        raise ValueError(f"{path}: invalid cell-centered angular grid")
    xyz = np.asarray(rings)
    if not np.isfinite(xyz).all():
        raise ValueError(f"{path}: nonfinite saved Cartesian point")
    roundoff = 1.0e-9 * float(np.max(np.abs(xyz)))
    if np.any(np.abs(xyz[:, 0] - xyz[:, -1]) > roundoff):
        raise ValueError(f"{path}: azimuthal ring does not close")
    ntheta, nphi = len(rings), len(rings[0]) - 1
    theta = (np.arange(ntheta)[:, None] + 0.5) * np.pi / ntheta
    phi = -np.pi + (np.arange(nphi)[None, :] + 0.5) * (2 * np.pi / nphi)
    direction = np.stack(
        (
            np.sin(theta) * np.cos(phi),
            np.sin(theta) * np.sin(phi),
            np.broadcast_to(np.cos(theta), (ntheta, nphi)),
        ),
        axis=-1,
    )
    # Project perpendicular to each radial direction to solve for the common origin.
    projection = np.eye(3) - direction[..., :, None] * direction[..., None, :]
    xyz = xyz[:, :-1]
    origin = np.linalg.solve(
        projection.sum(axis=(0, 1)), np.einsum("ijab,ijb->a", projection, xyz)
    )
    delta = xyz - origin
    radius = np.sum(delta * direction, axis=-1)
    if np.any(radius <= 0.0) or np.any(
        np.linalg.norm(delta - radius[..., None] * direction, axis=-1) > roundoff
    ):
        raise ValueError(f"{path}: saved points do not define a positive radial grid")
    iteration, centroid = row
    stream.write(
        SURFACE.pack(
            iteration,
            ntheta,
            nphi,
            *centroid,
            *origin,
            roundoff,
            radius[0].mean(),
            radius[-1].mean(),
        )
    )
    stream.write(radius.astype("<f8", copy=False).tobytes())


def combine(
    project_dir: Path, example_name: str, output_dir: Path, force: bool
) -> Path:
    """
    Write one apparent-horizon file for all spacetime modes of a run.

    :param project_dir: Evolution directory containing horizon files and metadata.
    :param example_name: Evolution example name.
    :param output_dir: Destination directory.
    :param force: Whether to replace an existing destination.
    :return: Written filename.
    :raises ValueError: If the run metadata or saved times disagree.
    :raises FileExistsError: If the destination exists and replacement is disabled.
    """
    metadata = json.loads(
        (project_dir / "raytracing_run_metadata.json").read_text(encoding="utf-8")
    )
    if metadata["project_name"] != example_name:
        raise ValueError(f"{project_dir}: metadata does not identify {example_name}")
    datasets: Dict[str, str] = metadata["datasets"]
    stems = {
        name[: -len(f"_{mode}.bin")]
        for mode, name in datasets.items()
        if Path(name).name == name and name.endswith(f"_{mode}.bin")
    }
    if len(stems) != 1 or any(
        name != f"{next(iter(stems))}_{mode}.bin" for mode, name in datasets.items()
    ):
        raise ValueError(f"{project_dir}: spacetime filenames do not share a run name")
    destination = output_dir / f"{stems.pop()}_apparent_horizons.bin"
    if destination.exists() and not force:
        raise FileExistsError(f"{destination} exists. Use --force to overwrite.")
    tables = [
        _spacetime_times(
            project_dir.parent / "raytracing_data" / filename,
            project_dir / "raytracing_slices" / mode,
        )
        for mode, filename in datasets.items()
    ]
    times = tables[0]
    if any(table != times for table in tables[1:]):
        raise ValueError(f"{project_dir}: spacetime modes have different slice tables")
    tracks = [
        _read_track(project_dir / f"BHaHAHA_diagnostics.ah{h}.gp") for h in (1, 2, 3)
    ]
    printed = {f"{time:.3f}": time for _, time in times}
    if len(printed) != len(times):
        raise ValueError(
            f"{project_dir}: spacetime times overlap at BHaHAHA's printed precision"
        )
    unmatched = set().union(*tracks).difference(printed)
    if unmatched:
        raise ValueError(
            f"{project_dir}: no matching spacetime output at t={min(unmatched, key=float)}"
        )
    common_time = min((printed[key] for key in tracks[2]), default=math.inf)
    coverage: List[List[float]] = [[], [], []]
    missing = [0, 0, 0]
    table = bytearray()
    output_dir.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=output_dir, prefix=f".{destination.name}.", delete=False
    ) as stream:
        temporary = Path(stream.name)
        try:
            stream.seek(HEADER.size + len(times) * SLICE.size)
            for output_index, time in times:
                offsets = [0, 0, 0]
                for horizon, track in enumerate(tracks):
                    row = track.get(f"{time:.3f}")
                    if row is not None:
                        path = project_dir / f"h.t{row[0]:07d}.ah{horizon + 1}.gp"
                        if path.is_file():
                            offsets[horizon] = stream.tell()
                            _write_surface(stream.file, path, row)
                            coverage[horizon].append(time)
                    if not offsets[horizon] and (
                        (time >= common_time) == (horizon == 2)
                    ):
                        missing[horizon] += 1
                table.extend(SLICE.pack(output_index, time, *offsets))
            metadata_offset = stream.tell()
            stream.write(
                json.dumps(
                    {
                        "project_name": example_name,
                        "datasets": datasets,
                        "common_horizon_detection_time": (
                            common_time if tracks[2] else None
                        ),
                        "horizon_time_ranges": [
                            [min(values), max(values)] if values else None
                            for values in coverage
                        ],
                    },
                    sort_keys=True,
                ).encode("utf-8")
            )
            stream.seek(0)
            stream.write(HEADER.pack(MAGIC, len(times), common_time, metadata_offset))
            stream.write(table)
            stream.flush()
            os.fsync(stream.fileno())
            if force:
                temporary.replace(destination)
            else:
                os.link(temporary, destination)
        finally:
            try:
                temporary.unlink()
            except FileNotFoundError:
                # A successful replacement already moved the temporary file.
                pass
    print(f"{destination} (missing required shapes for horizons 1, 2, 3: {missing})")
    return destination


def main(argv: Optional[List[str]] = None) -> None:
    """
    Combine the saved apparent horizon shapes for an evolution.

    :param argv: Command-line arguments, defaulting to the process arguments.
    """
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-dir", type=Path, required=True)
    parser.add_argument(
        "--example-name",
        choices=("blackhole_spectroscopy", "two_blackholes_collide"),
        required=True,
    )
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)
    project_dir = args.project_dir.resolve()
    combine(
        project_dir,
        args.example_name,
        args.output_dir or project_dir.parent / "raytracing_data",
        args.force,
    )


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        main()
    else:
        import doctest

        sys.exit(doctest.testmod().failed)

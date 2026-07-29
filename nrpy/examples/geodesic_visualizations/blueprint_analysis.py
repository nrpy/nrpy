"""
Validate and visualize native geodesic blueprint artifacts.

This utility checks the native tile headers and streams records for integrity,
termination-status, and coordinate diagnostics before plotting sampled data.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
from contextlib import ExitStack
from typing import BinaryIO, Dict, List, Optional, Tuple, cast

import numpy as np
import numpy.typing as npt

try:
    import blueprint_config_and_schema as cfg  # type: ignore
    import blueprint_io  # type: ignore
except ImportError:
    from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
    from nrpy.examples.geodesic_visualizations import blueprint_io

script_dir = os.path.dirname(os.path.abspath(__file__))

TERMINATION_TYPE_INFO: Dict[int, Tuple[str, str]] = {
    cfg.STOP_CONDITION_COORD_RADIUS_EXCEEDED: (
        "STOP_CONDITION_COORD_RADIUS_EXCEEDED",
        "Ray exceeded the coordinate-radius escape limit",
    ),
    cfg.STOP_CONDITION_TERMINAL_PLANE: (
        "STOP_CONDITION_TERMINAL_PLANE",
        "Ray hit the accretion disk / terminal plane",
    ),
    cfg.STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED: (
        "STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED",
        "Evolution measure exceeded its configured limit",
    ),
    cfg.FAILURE_RKF45_REJECTION_LIMIT: (
        "FAILURE_RKF45_REJECTION_LIMIT",
        "RKF45 rejected too many consecutive steps",
    ),
    cfg.STOP_CONDITION_T_MAX_EXCEEDED: (
        "STOP_CONDITION_T_MAX_EXCEEDED",
        "Ray exceeded the maximum allowed coordinate time",
    ),
    cfg.FAILURE_SLOT_MANAGER_ERROR: (
        "FAILURE_SLOT_MANAGER_ERROR",
        "Slot manager failed to handle the ray",
    ),
    cfg.FAILURE_GENERIC: ("FAILURE_GENERIC", "Unspecified integration failure"),
    cfg.ACTIVE: (
        "ACTIVE",
        "Ray is still being processed (should not appear in final blueprints)",
    ),
    cfg.REJECTED: (
        "REJECTED",
        "Ray is in a rejected RKF45 stage (not a final status)",
    ),
}


def _termination_type_info(enum_value: int) -> Tuple[str, str]:
    """
    Return the name and description for a termination enum.

    :param enum_value: Integer termination enum from a blueprint record.
    :return: Symbolic name and description, including an unknown-enum fallback.
    """
    return TERMINATION_TYPE_INFO.get(
        enum_value,
        (f"UNKNOWN_ENUM_{enum_value}", "Unknown termination status"),
    )


def _termination_type_label(enum_value: int) -> str:
    """
    Format a termination enum for plot legends.

    :param enum_value: Integer termination enum from a blueprint record.
    :return: Human-readable enum number and symbolic name.
    """
    name, _ = _termination_type_info(enum_value)
    return f"Type {enum_value}: {name}"


def _print_termination_diagnostics(
    enum_counts: Dict[int, int], total_rays: int
) -> None:
    """
    Print observed termination counts and every configured enum definition.

    :param enum_counts: Number of records observed for each termination enum.
    :param total_rays: Total number of records used for percentages.
    """
    print("--- Raw Termination Enums in Binary ---")
    for enum_value in sorted(enum_counts):
        count = enum_counts[enum_value]
        name, _ = _termination_type_info(enum_value)
        percentage = count / total_rays * 100.0 if total_rays else 0.0
        print(
            f"  Raw Enum {enum_value:2d} ({name}): "
            f"{count:12,} rays ({percentage:6.2f}%)"
        )

    print("\n--- Termination Enum Definitions ---")
    for enum_value in sorted(TERMINATION_TYPE_INFO):
        name, description = TERMINATION_TYPE_INFO[enum_value]
        print(f"  Enum {enum_value:2d} ({name}): {description}")

    configured_status = cfg.STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED
    configured_name, _ = _termination_type_info(configured_status)
    print(
        f"\n  Configured evolution-measure status = {configured_status} "
        f"({configured_name})"
    )


def _calculate_subsample_rate(total_expected_records: int, max_viz_points: int) -> int:
    """
    Calculate a visualization subsampling rate.

    :param total_expected_records: Number of records across all tiles.
    :param max_viz_points: Maximum number of points to retain.
    :return: Subsampling rate, always at least one.
    """
    return max(1, (total_expected_records + max_viz_points - 1) // max_viz_points)


def _has_nonterminal_diagnostic(data: npt.NDArray[np.void]) -> npt.NDArray[np.bool_]:
    """
    Identify records containing a nonterminal-plane crossing.

    The v6 blueprint record has no separate nonterminal event flag.  Its
    ``non_terminal_plane_lambda`` field is initialized to zero and is written
    only after a successful nonterminal crossing, so a nonzero finite value is
    the serialized crossing marker.  This is independent of
    ``termination_type``: a nonterminal crossing does not terminate a ray.

    :param data: Structured blueprint records.
    :return: Boolean mask for records with usable nonterminal diagnostics.
    """
    return cast(
        npt.NDArray[np.bool_],
        (
            np.isfinite(data["non_terminal_plane_lambda"])
            & (data["non_terminal_plane_lambda"] != 0.0)
            & np.isfinite(data["y_nt"])
            & np.isfinite(data["z_nt"])
        ),
    )


def plot_heatmaps(data: npt.NDArray[np.void]) -> None:
    """
    Generate diagnostic heatmaps for the three recorded event classes.

    Nonterminal coordinates are plotted only for records that crossed the
    nonterminal plane. Terminal coordinates are plotted only for records whose
    final status is ``STOP_CONDITION_TERMINAL_PLANE``. These plots never place pixels;
    image placement belongs to the renderer and uses normalized image-sample
    fractions.

    :param data: Structured blueprint records to plot.
    """
    # pylint: disable=import-outside-toplevel, import-error
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    nonterminal_mask = _has_nonterminal_diagnostic(data)
    if np.any(nonterminal_mask):
        hb0 = axes[0].hexbin(
            data["y_nt"][nonterminal_mask],
            data["z_nt"][nonterminal_mask],
            gridsize=50,
            cmap="viridis",
            mincnt=1,
        )
        fig.colorbar(hb0, ax=axes[0], label="Ray Count")
    axes[0].set_title("Nonterminal Plane Diagnostics\ny_nt vs z_nt")
    axes[0].set_xlabel("y_nt")
    axes[0].set_ylabel("z_nt")

    terminal_mask = (
        (data["termination_type"] == cfg.STOP_CONDITION_TERMINAL_PLANE)
        & np.isfinite(data["y_t"])
        & np.isfinite(data["z_t"])
    )
    if np.any(terminal_mask):
        hb1 = axes[1].hexbin(
            data["y_t"][terminal_mask],
            data["z_t"][terminal_mask],
            gridsize=50,
            cmap="inferno",
            mincnt=1,
        )
        fig.colorbar(hb1, ax=axes[1], label="Ray Count")
    axes[1].set_title("Terminal Plane\ny_t vs z_t")
    axes[1].set_xlabel("y_t")
    axes[1].set_ylabel("z_t")

    sphere_mask = (
        (data["termination_type"] == cfg.STOP_CONDITION_COORD_RADIUS_EXCEEDED)
        & np.isfinite(data["final_phi"])
        & np.isfinite(data["final_theta"])
    )
    if np.any(sphere_mask):
        hb2 = axes[2].hexbin(
            data["final_phi"][sphere_mask],
            data["final_theta"][sphere_mask],
            gridsize=50,
            cmap="magma",
            mincnt=1,
        )
        fig.colorbar(hb2, ax=axes[2], label="Ray Count")
    axes[2].set_title("Celestial Sphere\n$\\phi$ vs $\\theta$")
    axes[2].set_xlabel("$\\phi$")
    axes[2].set_ylabel("$\\theta$")
    plt.tight_layout()
    plt.show()


def plot_norm_abs_log_histogram(
    termination_types: npt.NDArray[np.int32],
    norm_abs_values: npt.NDArray[np.float64],
) -> None:
    """
    Plot normalization magnitudes grouped by terminal status.

    :param termination_types: Termination enum for each ray.
    :param norm_abs_values: Absolute normalization values for each ray.
    """
    # pylint: disable=import-outside-toplevel, import-error
    import matplotlib.pyplot as plt

    plt.figure(figsize=(10, 6))
    termination_colors = plt.get_cmap("tab10")
    for enum_value in np.unique(termination_types):
        values = norm_abs_values[termination_types == enum_value]
        positive_values = values[values > 0.0]
        if len(positive_values) > 0:
            enum_int = int(enum_value)
            plt.hist(
                np.log10(positive_values),
                bins=100,
                alpha=0.5,
                color=termination_colors(enum_int % 10),
                label=_termination_type_label(enum_int),
            )
    plt.title("Log-Scale Normalization Magnitude by Termination Type")
    plt.xlabel("$\\log_{10}(|g_{{\\mu\\nu}}p^\\mu p^\\nu|)$")
    plt.yscale("log", base=10)
    plt.ylabel("$\\log_{10}(N_{\\mathrm{photons}})$")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def diagnose_blueprint() -> None:
    """
    Read every metadata-described tile and print integrity/physics diagnostics.

    :raises FileNotFoundError: If an expected tile is missing.
    :raises ValueError: If an artifact or normalization sidecar is malformed.
    """
    first_tile = os.path.join(script_dir, "light_blueprint_00_00.bin")
    if not os.path.exists(first_tile):
        raise FileNotFoundError(f"Expected first blueprint tile: {first_tile}")
    first_header = blueprint_io.read_blueprint_header(first_tile, expected_tile=(0, 0))
    tiles_width = first_header.tiles_width
    tiles_height = first_header.tiles_height

    print("=================================================================")
    print(" BLUEPRINT DIAGNOSTICS & VISUALIZATION")
    print(f" Grid: {tiles_width}x{tiles_height} Tiles")
    print("=================================================================")

    blueprint_files: List[str] = []
    headers: List[blueprint_io.BlueprintHeader] = []
    missing_sidecars: List[str] = []
    total_expected_records = 0
    for tile_x in range(tiles_width):
        for tile_y in range(tiles_height):
            filename = f"light_blueprint_{tile_x:02d}_{tile_y:02d}.bin"
            filepath = os.path.join(script_dir, filename)
            if not os.path.exists(filepath):
                raise FileNotFoundError(f"Expected tile missing: {filepath}")
            header = blueprint_io.read_blueprint_header(
                filepath, expected_tile=(tile_x, tile_y)
            )
            if header.tiles_width != tiles_width or header.tiles_height != tiles_height:
                raise ValueError(f"Blueprint '{filename}' has an unexpected tile grid")
            blueprint_files.append(filepath)
            headers.append(header)
            total_expected_records += header.record_count

            sidecar_name = cfg.BLUEPRINT_NORM_ABS_FILENAME_TEMPLATE.format(
                tile_x=tile_x, tile_y=tile_y
            )
            sidecar_path = os.path.join(script_dir, sidecar_name)
            if not os.path.exists(sidecar_path):
                missing_sidecars.append(sidecar_name)
            elif (
                os.path.getsize(sidecar_path)
                != header.record_count * cfg.BLUEPRINT_NORM_ABS_DTYPE.itemsize
            ):
                raise ValueError(
                    f"Normalization sidecar '{sidecar_name}' has an invalid size"
                )

    max_viz_points = 2_000_000
    subsample_rate = _calculate_subsample_rate(total_expected_records, max_viz_points)
    have_sidecars = not missing_sidecars
    print(f"Total records detected: {total_expected_records:,}")
    print(
        f"Streaming data and sampling 1 in every {subsample_rate} rays for visualization..."
    )
    if not have_sidecars:
        print("[!] Missing normalization sidecars; skipping normalization histogram.")

    total_rays = 0
    in_view = 0
    enum_counts: Dict[int, int] = {}
    first_records: List[npt.NDArray[np.void]] = []
    sampled_records: List[npt.NDArray[np.void]] = []
    sampled_norm_abs: List[npt.NDArray[np.float64]] = []
    for filepath, header in zip(blueprint_files, headers):
        sidecar_path = os.path.join(
            script_dir,
            cfg.BLUEPRINT_NORM_ABS_FILENAME_TEMPLATE.format(
                tile_x=header.tile_x, tile_y=header.tile_y
            ),
        )
        with ExitStack() as sidecar_stack:
            sidecar_file: Optional[BinaryIO] = None
            if have_sidecars:
                sidecar_file = sidecar_stack.enter_context(open(sidecar_path, "rb"))
            for _, _, chunk_data in blueprint_io.iter_blueprint_chunks(
                filepath, cfg.CHUNK_SIZE
            ):
                norm_abs_chunk = None
                if sidecar_file is not None:
                    expected_bytes = (
                        len(chunk_data) * cfg.BLUEPRINT_NORM_ABS_DTYPE.itemsize
                    )
                    norm_abs_bytes = sidecar_file.read(expected_bytes)
                    if len(norm_abs_bytes) != expected_bytes:
                        raise ValueError(
                            f"Normalization sidecar for '{filepath}' is truncated"
                        )
                    norm_abs_chunk = np.frombuffer(
                        norm_abs_bytes, dtype=cfg.BLUEPRINT_NORM_ABS_DTYPE
                    )

                unique_enums, counts = np.unique(
                    chunk_data["termination_type"], return_counts=True
                )
                for enum_value, count in zip(unique_enums, counts):
                    enum_int = int(enum_value)
                    enum_counts[enum_int] = enum_counts.get(enum_int, 0) + int(count)
                total_rays += len(chunk_data)
                nonterminal_mask = _has_nonterminal_diagnostic(chunk_data)
                in_view += int(np.sum(nonterminal_mask))
                if len(first_records) < 9:
                    first_records.extend(chunk_data[: 9 - len(first_records)])
                sampled_records.append(chunk_data[::subsample_rate])
                if norm_abs_chunk is not None:
                    sampled_norm_abs.append(norm_abs_chunk[::subsample_rate])
            if sidecar_file is not None and sidecar_file.read(1):
                raise ValueError(
                    f"Normalization sidecar for '{filepath}' has extra data"
                )

    _print_termination_diagnostics(enum_counts, total_rays)
    print("\n--- Nonterminal-Plane Diagnostics ---")
    print(
        "  Rays with a recorded nonterminal-plane crossing: "
        f"{in_view:,} ({in_view/total_rays*100:.2f}%)"
    )
    print("\n--- First 5 Raw Records ---")
    print(f"{'Ray#':<6} | {'Term Type':<38} | " f"{'y_nt':>17} | {'z_nt':>17}")
    print("-" * 71)
    for index, record in enumerate(first_records[:5]):
        enum_value = int(record["termination_type"])
        enum_name, _ = _termination_type_info(enum_value)
        print(
            f"{index:<6} | {enum_value} ({enum_name}) | "
            f"{record['y_nt']:>17.5f} | "
            f"{record['z_nt']:>17.5f}"
        )

    if sampled_records:
        viz_data = np.concatenate(sampled_records)
        plot_heatmaps(viz_data)
        if have_sidecars:
            plot_norm_abs_log_histogram(
                viz_data["termination_type"], np.concatenate(sampled_norm_abs)
            )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=(
            "Diagnose blueprint artifacts using tile counts and fields of view "
            "stored in their headers."
        )
    )
    parser.parse_args()

    first_blueprint_tile = os.path.join(script_dir, "light_blueprint_00_00.bin")
    if not os.path.exists(first_blueprint_tile):
        print(
            f"No native blueprint artifacts found in {script_dir}; nothing to diagnose."
        )
    else:
        diagnose_blueprint()

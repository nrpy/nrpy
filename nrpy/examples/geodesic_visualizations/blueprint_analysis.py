# nrpy/examples/geodesic_visualizations/blueprint_analysis.py
"""
Defines the blueprint diagnostic and visualization suite.

This utility validates the tiled 'light_blueprint_XX_YY.zip' archives produced by
the C-based Geodesic Integrator.

It performs three functions:
1. Data Integrity: Parses the binary file using the structured 'BLUEPRINT_DTYPE' to
   ensure the data alignment matches the C 'blueprint_data_t' struct.
2. Statistical Validation: Calculates the distribution of termination types (hits on
   the source plane vs. escape to the celestial sphere) to verify physics enums.
3. Spatial Visualization: Generates hexbin heatmaps for window coordinates, source
   plane intercepts, and celestial sphere angles to check for coordinate mapping
   errors or integration 'blind spots'.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
from typing import Dict, List, Optional, Tuple, cast

import numpy as np
import numpy.typing as npt

try:
    import blueprint_config_and_schema as cfg  # type: ignore
except ImportError:
    import nrpy.examples.geodesic_visualizations.blueprint_config_and_schema as cfg

script_dir = os.path.dirname(os.path.abspath(__file__))


def _calculate_subsample_rate(total_expected_records: int, max_viz_points: int) -> int:
    """
    Calculate the subsampling rate to cap the total visualization points.

    >>> _calculate_subsample_rate(1_000_000, 2_000_000)
    1
    >>> _calculate_subsample_rate(2_000_000, 2_000_000)
    1
    >>> _calculate_subsample_rate(2_000_001, 2_000_000)
    2
    >>> _calculate_subsample_rate(4_000_000, 2_000_000)
    2
    >>> _calculate_subsample_rate(0, 2_000_000)
    1

    :param total_expected_records: Total number of records across all files.
    :param max_viz_points: The maximum number of points to visualize.
    :return: The subsampling rate (always >= 1).
    """
    return max(1, (total_expected_records + max_viz_points - 1) // max_viz_points)


def plot_heatmaps(data: "npt.NDArray[np.void]") -> None:
    """
    Generate hexbin heatmaps for window, source plane, and sphere coordinates.

    Hexbins are used instead of scatter plots to handle large datasets efficiently
    without overlapping individual points (overplotting).

    :param data: The structured NumPy array containing ray termination data.
    """
    # pylint: disable=import-outside-toplevel, import-error
    import matplotlib.pyplot as plt

    # Preamble: Descriptive Physical Variable Mapping.
    # Extract coordinates from the structured array for plotting.
    y_w, z_w = data["y_w"], data["z_w"]
    y_s, z_s = data["y_s"], data["z_s"]
    phi, theta = data["final_phi"], data["final_theta"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Step 1: Window Coordinates (y_w, z_w) - All rays.
    # Validates the initial ray distribution at the camera's field of view.
    # This should typically show a uniform grid or box representing the 'scan_density'.
    hb0 = axes[0].hexbin(y_w, z_w, gridsize=50, cmap="viridis", mincnt=1)
    axes[0].set_title("Window Plane (Camera)\n$y_w$ vs $z_w$")
    axes[0].set_xlabel("$y_w$")
    axes[0].set_ylabel("$z_w$")
    fig.colorbar(hb0, ax=axes[0], label="Ray Count")

    # Step 2: Source Plane (y_s, z_s) - Filtering for TERM_SOURCE_PLANE.
    # Visualizes the geometry of the accretion disk hits. Expect a 'ring' or 'annulus'
    # structure corresponding to [source_r_min, source_r_max].
    source_mask = data["termination_type"] == cfg.TERM_SOURCE_PLANE
    if np.any(source_mask):
        hb1 = axes[1].hexbin(
            y_s[source_mask],
            z_s[source_mask],
            gridsize=50,
            cmap="plasma",
            mincnt=1,
        )
        axes[1].set_title("Source Plane Hits\n$y_s$ vs $z_s$")
        fig.colorbar(hb1, ax=axes[1], label="Ray Count")
    else:
        # Diagnostic warning: if no hits appear, the C event detection might be failing.
        axes[1].text(
            0.5,
            0.5,
            "No Source Plane Hits\n(Check Enum Mapping)",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
        )
    axes[1].set_xlabel("$y_s$")
    axes[1].set_ylabel("$z_s$")

    # Step 3: Celestial Sphere (phi, theta) - Filtering for coordinate-radius exits.
    # Maps rays that escaped to 'infinity'.
    # theta is the polar angle [0, pi], phi is the azimuthal angle [0, 2pi].
    sphere_mask = data["termination_type"] == cfg.TERM_COORD_RADIUS_EXCEEDED
    if np.any(sphere_mask):
        hb2 = axes[2].hexbin(
            phi[sphere_mask],
            theta[sphere_mask],
            gridsize=50,
            cmap="magma",
            mincnt=1,
        )
        axes[2].set_title(r"Celestial Sphere Hits" + "\n" + r"$\phi$ vs $\theta$")
        fig.colorbar(hb2, ax=axes[2], label="Ray Count")
    else:
        axes[2].text(
            0.5,
            0.5,
            "No Sphere Hits\n(Check Enum Mapping)",
            ha="center",
            va="center",
            transform=axes[2].transAxes,
        )
    axes[2].set_xlabel(r"$\phi$")
    axes[2].set_ylabel(r"$\theta$")

    plt.tight_layout()
    print("[i] Displaying heatmaps. Close the window to continue...")
    plt.show()


def plot_norm_abs_log_histogram(
    termination_types: "npt.NDArray[np.int32]",
    norm_abs_values: "npt.NDArray[np.float64]",
    num_bins: int = 80,
) -> None:
    """
    Generate a log-scaled histogram of absolute norm values by termination enum.

    This helper treats both inputs as one-dimensional ray streams. The arrays must
    remain index-aligned so each normalization value corresponds to the same ray as
    its termination enum.

    :param termination_types: Per-ray termination enum values.
    :param norm_abs_values: Per-ray absolute norm values.
    :param num_bins: Number of logarithmic bins to use.
    :raises ValueError: If num_bins is non-positive or the inputs have mismatched sizes.
    """
    # pylint: disable=import-outside-toplevel, import-error
    import matplotlib.pyplot as plt

    # Step 1: Validate and flatten the per-ray input streams.
    if num_bins < 1:
        raise ValueError("num_bins must be positive.")

    termination_type_by_ray = np.asarray(termination_types, dtype=np.int32).ravel()
    norm_abs_by_ray = np.asarray(norm_abs_values, dtype=np.float64).ravel()

    if termination_type_by_ray.size != norm_abs_by_ray.size:
        raise ValueError(
            "termination_types and norm_abs_values must contain the same number "
            "of per-ray entries."
        )

    if norm_abs_by_ray.size == 0:
        print("[!] No normalization values available to plot.")
        return

    # Step 2: Keep only positive finite values that can appear on a log axis.
    finite_positive_mask = np.isfinite(norm_abs_by_ray) & (norm_abs_by_ray > 0.0)
    if not np.any(finite_positive_mask):
        print("[!] No finite positive normalization values available to plot.")
        return

    plotted_termination_types = termination_type_by_ray[finite_positive_mask]
    plotted_norm_abs = norm_abs_by_ray[finite_positive_mask]
    skipped_count = norm_abs_by_ray.size - int(np.count_nonzero(finite_positive_mask))

    # Step 3: Build shared logarithmic bins so enum distributions are comparable.
    min_norm_abs = float(np.min(plotted_norm_abs))
    max_norm_abs = float(np.max(plotted_norm_abs))
    if min_norm_abs == max_norm_abs:
        lower_edge = max(np.nextafter(0.0, 1.0), min_norm_abs / 2.0)
        upper_edge = max_norm_abs * 2.0
        if upper_edge <= lower_edge:
            upper_edge = np.nextafter(lower_edge, np.inf)
        bin_edges = np.array([lower_edge, upper_edge], dtype=np.float64)
    else:
        bin_edges = np.logspace(
            np.log10(min_norm_abs),
            np.log10(max_norm_abs),
            num_bins + 1,
        )
    bin_edges_list = cast(List[float], bin_edges.tolist())

    enum_label_by_value: Dict[int, str] = {
        cfg.TERM_COORD_RADIUS_EXCEEDED: "TERM_COORD_RADIUS_EXCEEDED",
        cfg.TERM_SOURCE_PLANE: "TERM_SOURCE_PLANE",
        cfg.TERM_ENERGY_LIMIT_EXCEEDED: "TERM_ENERGY_LIMIT_EXCEEDED",
        cfg.TERM_RKF45_REJECTION_LIMIT: "TERM_RKF45_REJECTION_LIMIT",
        cfg.TERM_T_MAX_EXCEEDED: "TERM_T_MAX_EXCEEDED",
        cfg.TERM_SLOT_MANAGER_ERROR: "TERM_SLOT_MANAGER_ERROR",
        cfg.TERM_FAILURE: "TERM_FAILURE",
        cfg.TERM_ACTIVE: "TERM_ACTIVE",
        cfg.TERM_REJECTED: "TERM_REJECTED",
    }

    # Step 4: Draw one histogram trace per enum, forcing energy-limit failures to black.
    fig, ax = plt.subplots(figsize=(10, 6))
    color_cycle = plt.rcParams["axes.prop_cycle"].by_key().get("color", ["C0"])
    color_index = 0

    for enum_value_raw in sorted(np.unique(plotted_termination_types)):
        enum_value = int(enum_value_raw)
        enum_mask = plotted_termination_types == enum_value
        enum_norm_abs = plotted_norm_abs[enum_mask]

        if enum_value == cfg.TERM_ENERGY_LIMIT_EXCEEDED:
            color = "black"
        else:
            color = color_cycle[color_index % len(color_cycle)]
            color_index += 1

        enum_label = enum_label_by_value.get(enum_value, f"Raw Enum {enum_value}")
        ax.hist(
            enum_norm_abs,
            bins=bin_edges_list,
            histtype="step",
            linewidth=1.8,
            color=color,
            label=f"{enum_label} ({enum_norm_abs.size:,})",
        )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_title(r"Absolute Norm Distribution by Termination Type")
    ax.set_xlabel(r"$|\mathrm{norm}|$")
    ax.set_ylabel("Ray Count")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

    if skipped_count > 0:
        print(
            f"[i] Skipped {skipped_count:,} non-positive or non-finite "
            "normalization values."
        )

    plt.tight_layout()
    print("[i] Displaying |norm| log histogram. Close the window to continue...")
    plt.show()


def diagnose_blueprint(
    window_tiles_width: int = 1,
    window_tiles_height: int = 1,
    window_width: float = 1.0,
    window_height: float = 1.0,
) -> None:
    """
    Read the binary blueprint files and print diagnostics to the console.

    This helps identify numerical instabilities or logic errors in the integrator.

    :param window_tiles_width: Number of horizontal tiles generated.
    :param window_tiles_height: Number of vertical tiles generated.
    :param window_width: The physical width of the camera's local window.
    :param window_height: The physical height of the camera's local window.
    :raises ValueError: If a required binary member is missing or sidecar sizes mismatch.
    """
    # pylint: disable=import-outside-toplevel
    import zipfile

    print("=================================================================")
    print(" BLUEPRINT DIAGNOSTICS & VISUALIZATION")
    print(f" Grid: {window_tiles_width}x{window_tiles_height} Tiles")
    print("=================================================================")

    blueprint_archives: List[Tuple[str, Optional[str]]] = []
    total_expected_records = 0
    missing_norm_abs_archives: List[str] = []

    # Step 1: Identify valid tiles and pre-calculate total records for safe subsampling.
    for i in range(window_tiles_width):
        for j in range(window_tiles_height):
            zip_filename = f"light_blueprint_{i:02d}_{j:02d}.zip"
            zip_filepath = os.path.join(script_dir, zip_filename)
            if os.path.exists(zip_filepath):
                norm_abs_zip_filename = cfg.BLUEPRINT_NORM_ABS_ARCHIVE_TEMPLATE.format(
                    tile_x=i, tile_y=j
                )
                norm_abs_zip_filepath = os.path.join(script_dir, norm_abs_zip_filename)
                if os.path.exists(norm_abs_zip_filepath):
                    blueprint_archives.append((zip_filepath, norm_abs_zip_filepath))
                else:
                    blueprint_archives.append((zip_filepath, None))
                    missing_norm_abs_archives.append(norm_abs_zip_filename)

                with zipfile.ZipFile(zip_filepath, "r") as zf:
                    blueprint_bin_filename = zip_filename.replace(".zip", ".bin")
                    blueprint_members = [
                        name
                        for name in zf.namelist()
                        if name.endswith(blueprint_bin_filename)
                    ]
                    if not blueprint_members:
                        raise ValueError(
                            f"Could not find '{blueprint_bin_filename}' in '{zip_filename}'."
                        )
                    itemsize = cfg.BLUEPRINT_DTYPE.itemsize
                    file_size = zf.getinfo(blueprint_members[0]).file_size
                    total_expected_records += file_size // itemsize

    if not blueprint_archives:
        print("[!] ERROR: No blueprint zip files found to analyze.")
        return

    # Step 2: To prevent OOM errors during plotting, cap visualization to ~2 million points.
    max_viz_points = 2_000_000
    subsample_rate = _calculate_subsample_rate(total_expected_records, max_viz_points)
    have_complete_norm_abs_sidecars = len(missing_norm_abs_archives) == 0

    total_rays = 0
    in_view = 0
    enum_counts: Dict[int, int] = {}
    first_records: List["npt.NDArray[np.void]"] = []
    sampled_data_list: List["npt.NDArray[np.void]"] = []
    sampled_norm_abs_list: List["npt.NDArray[np.float64]"] = []

    chunk_size = cfg.CHUNK_SIZE
    chunk_bytes = chunk_size * cfg.BLUEPRINT_DTYPE.itemsize
    norm_abs_chunk_itemsize = cfg.BLUEPRINT_NORM_ABS_DTYPE.itemsize

    print(f"Total records detected: {total_expected_records:,}")
    print(
        f"Streaming data and sampling 1 in every {subsample_rate} rays for "
        "visualization...\n"
    )
    if not have_complete_norm_abs_sidecars:
        print(
            "[!] Missing one or more normalization sidecar archives. "
            "Skipping |norm| histogram generation."
        )

    # Step 3: Stream the data chunk by chunk.
    for zip_filepath, paired_norm_abs_zip_filepath in blueprint_archives:
        zip_filename = os.path.basename(zip_filepath)
        blueprint_bin_filename = zip_filename.replace(".zip", ".bin")
        norm_abs_bin_filename = ""
        if paired_norm_abs_zip_filepath is not None:
            norm_abs_bin_filename = os.path.basename(
                paired_norm_abs_zip_filepath
            ).replace(".zip", ".bin")

        with zipfile.ZipFile(zip_filepath, "r") as blueprint_zf:
            blueprint_members = [
                name
                for name in blueprint_zf.namelist()
                if name.endswith(blueprint_bin_filename)
            ]
            if not blueprint_members:
                continue

            if (
                have_complete_norm_abs_sidecars
                and paired_norm_abs_zip_filepath is not None
            ):
                with zipfile.ZipFile(paired_norm_abs_zip_filepath, "r") as norm_abs_zf:
                    norm_abs_members = [
                        name
                        for name in norm_abs_zf.namelist()
                        if name.endswith(norm_abs_bin_filename)
                    ]
                    if not norm_abs_members:
                        raise ValueError(
                            f"Could not find '{norm_abs_bin_filename}' in "
                            f"'{os.path.basename(paired_norm_abs_zip_filepath)}'."
                        )

                    with blueprint_zf.open(
                        blueprint_members[0], "r"
                    ) as blueprint_file, norm_abs_zf.open(
                        norm_abs_members[0], "r"
                    ) as norm_abs_file:
                        while True:
                            raw_bytes = blueprint_file.read(chunk_bytes)
                            if not raw_bytes:
                                break

                            chunk_data = np.frombuffer(
                                raw_bytes, dtype=cfg.BLUEPRINT_DTYPE
                            )
                            chunk_norm_abs_bytes = norm_abs_file.read(
                                len(chunk_data) * norm_abs_chunk_itemsize
                            )
                            expected_norm_abs_bytes = (
                                len(chunk_data) * norm_abs_chunk_itemsize
                            )
                            if len(chunk_norm_abs_bytes) != expected_norm_abs_bytes:
                                raise ValueError(
                                    f"Normalization sidecar '{norm_abs_bin_filename}' "
                                    f"is misaligned with '{blueprint_bin_filename}'."
                                )
                            chunk_norm_abs = np.frombuffer(
                                chunk_norm_abs_bytes,
                                dtype=cfg.BLUEPRINT_NORM_ABS_DTYPE,
                            )

                            # Step 3.a: Accumulate Enum Mappings.
                            unique_enums, counts = np.unique(
                                chunk_data["termination_type"], return_counts=True
                            )
                            for e, c in zip(unique_enums, counts):
                                enum_counts[e] = enum_counts.get(e, 0) + c

                            # Step 3.b: Accumulate Window Coordinate Bounds.
                            total_rays += len(chunk_data)
                            half_w = window_width / 2.0
                            half_h = window_height / 2.0
                            in_view += int(
                                np.sum(
                                    (chunk_data["y_w"] >= -half_w)
                                    & (chunk_data["y_w"] < half_w)
                                    & (chunk_data["z_w"] >= -half_h)
                                    & (chunk_data["z_w"] < half_h)
                                )
                            )

                            # Step 3.c: Accumulate Raw Sample Data.
                            if len(first_records) < 9:
                                first_records.extend(
                                    chunk_data[: 9 - len(first_records)]
                                )

                            # Step 3.d: Subsample paired data for memory-safe plotting.
                            sampled_data_list.append(chunk_data[::subsample_rate])
                            sampled_norm_abs_list.append(
                                chunk_norm_abs[::subsample_rate]
                            )

                        if norm_abs_file.read(1):
                            raise ValueError(
                                f"Normalization sidecar '{norm_abs_bin_filename}' "
                                f"contains extra records beyond '{blueprint_bin_filename}'."
                            )
            else:
                with blueprint_zf.open(blueprint_members[0], "r") as blueprint_file:
                    while True:
                        raw_bytes = blueprint_file.read(chunk_bytes)
                        if not raw_bytes:
                            break

                        chunk_data = np.frombuffer(raw_bytes, dtype=cfg.BLUEPRINT_DTYPE)

                        # Step 3.a: Accumulate Enum Mappings.
                        unique_enums, counts = np.unique(
                            chunk_data["termination_type"], return_counts=True
                        )
                        for e, c in zip(unique_enums, counts):
                            enum_counts[e] = enum_counts.get(e, 0) + c

                        # Step 3.b: Accumulate Window Coordinate Bounds.
                        total_rays += len(chunk_data)
                        half_w = window_width / 2.0
                        half_h = window_height / 2.0
                        in_view += int(
                            np.sum(
                                (chunk_data["y_w"] >= -half_w)
                                & (chunk_data["y_w"] < half_w)
                                & (chunk_data["z_w"] >= -half_h)
                                & (chunk_data["z_w"] < half_h)
                            )
                        )

                        # Step 3.c: Accumulate Raw Sample Data.
                        if len(first_records) < 9:
                            first_records.extend(chunk_data[: 9 - len(first_records)])

                        # Step 3.d: Subsample the chunk for memory-safe plotting.
                        sampled_data_list.append(chunk_data[::subsample_rate])

    # Crucial for verifying that C enums and Python enums match.
    # If the C code terminates at the sphere but uses ID 9, this check will flag it.
    print("--- 1. Raw Termination Enums in Binary ---")
    for e in sorted(enum_counts.keys()):
        c = enum_counts[e]
        print(f"  Raw Enum {e:2d}: {c:12,} rays ({c/total_rays*100:6.2f}%)")

    print(
        f"\n  [Config Current]: COORD_RADIUS_EXCEEDED = {cfg.TERM_COORD_RADIUS_EXCEEDED}, "
        f"SOURCE_PLANE = {cfg.TERM_SOURCE_PLANE}, "
        f"ENERGY_LIMIT_EXCEEDED = {cfg.TERM_ENERGY_LIMIT_EXCEEDED}"
    )
    print(
        "  -> If your raw enums above do NOT match these, update blueprint_config_and_schema.py."
    )

    # Validates that the output rays actually align with the window dimensions
    # specified in the Python configuration.
    print("\n--- 2. Window Coordinate Bounds (y_w, z_w) ---")
    print(f"  Rays inside Renderer FOV: {in_view:,} ({in_view/total_rays*100:.2f}%)")

    # Direct sanity check of the first 5 records.
    # Useful for spotting immediate NANs or impossible values (e.g., negative theta).
    print("\n--- 3. First 5 Raw Records ---")
    header = (
        f"{'Ray#':<6} | {'Term Type':<9} | {'y_w':>8} | {'z_w':>8} | "
        f"{'y_s / theta':>11} | {'z_s / phi':>11}"
    )
    print(header)
    print("-" * len(header))

    for i, rec in enumerate(first_records):
        tt = rec["termination_type"]
        # Context-aware display: show plane coords for disk hits, angles for
        # sphere hits.
        val_1 = rec["y_s"] if tt == cfg.TERM_SOURCE_PLANE else rec["final_theta"]
        val_2 = rec["z_s"] if tt == cfg.TERM_SOURCE_PLANE else rec["final_phi"]
        print(
            f"{i:<6} | {tt:<9} | {rec['y_w']:>8.5} | {rec['z_w']:>8.5f} | "
            f"{val_1:>11.5f} | {val_2:>11.5f}"
        )

    # Step 4: Launch Visualization.
    print("\n[i] Generating heatmaps...")
    if sampled_data_list:
        viz_data = np.concatenate(sampled_data_list)
        plot_heatmaps(viz_data)
        if have_complete_norm_abs_sidecars:
            viz_norm_abs = np.concatenate(sampled_norm_abs_list)
            print("\n[i] Generating |norm| histogram...")
            plot_norm_abs_log_histogram(viz_data["termination_type"], viz_norm_abs)
    else:
        print("[!] No data available to plot.")
    print("=================================================================")


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

    import argparse

    parser = argparse.ArgumentParser(description="Blueprint Diagnostic Suite")
    parser.add_argument(
        "--window_tiles_width",
        type=int,
        default=1,
        help="Number of horizontal tiles generated.",
    )
    parser.add_argument(
        "--window_tiles_height",
        type=int,
        default=1,
        help="Number of vertical tiles generated.",
    )
    parser.add_argument(
        "--window_width",
        type=float,
        default=1.0,
        help="Camera window width w in coordinate units.",
    )
    parser.add_argument(
        "--window_height",
        type=float,
        default=1.0,
        help="Camera window height h in coordinate units.",
    )
    args = parser.parse_args()

    diagnose_blueprint(
        args.window_tiles_width,
        args.window_tiles_height,
        args.window_width,
        args.window_height,
    )

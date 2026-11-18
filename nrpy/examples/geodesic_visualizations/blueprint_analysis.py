"""
Validate and visualize native geodesic blueprint binary files.

This utility checks the native tile headers and streams records for integrity,
termination-status, coordinate, and optional normalization diagnostics before
plotting sampled data. Eight-neighbor domains use every photon, including across
tile edges, and are numbered by increasing seed column, then seed row.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
from contextlib import ExitStack
from typing import BinaryIO, Dict, List, NamedTuple, Optional, Tuple, cast

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
from matplotlib.backend_bases import MouseEvent
from matplotlib.colors import ListedColormap
from matplotlib.figure import Figure
from matplotlib.widgets import TextBox

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
    cfg.STOP_CONDITION_INSIDE_APPARENT_HORIZON: (
        "STOP_CONDITION_INSIDE_APPARENT_HORIZON",
        "Ray stopped inside a saved apparent horizon",
    ),
    cfg.FAILURE_RKF45_REJECTION_LIMIT: (
        "FAILURE_RKF45_REJECTION_LIMIT",
        "RKF45 rejected too many consecutive steps",
    ),
    cfg.FAILURE_RKF45_MIN_STEP_LIMIT: (
        "FAILURE_RKF45_MIN_STEP_LIMIT",
        "RKF45 reached the consecutive minimum-step limit",
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
    cfg.FAILURE_SPATIAL_INTERPOLATION: (
        "FAILURE_SPATIAL_INTERPOLATION",
        "Spatial interpolation failed for the ray",
    ),
    cfg.FAILURE_TEMPORAL_INTERPOLATION: (
        "FAILURE_TEMPORAL_INTERPOLATION",
        "Temporal interpolation failed for the ray",
    ),
    cfg.ACTIVE: (
        "ACTIVE",
        "Ray is still being processed (should not appear in final blueprints)",
    ),
    cfg.REJECTED: (
        "REJECTED",
        "Ray is in a rejected RKF45 stage (not a final status)",
    ),
}


class BlueprintDomain(NamedTuple):
    """One status-connected domain; seeds and inclusive bounds use column/row indices."""

    domain_id: int
    status: int
    photon_count: int
    seed: Tuple[int, int]
    bounds: Tuple[int, int, int, int]
    boundary_count: int
    edges: Tuple[str, ...]


def label_domains(
    statuses: npt.NDArray[np.int32],
) -> Tuple[npt.NDArray[np.intp], List[BlueprintDomain]]:
    """
    Label equal-status samples joined horizontally, vertically, or diagonally.

    Rows increase with b and columns with a. Seeds are visited column first,
    then row, so IDs follow the smallest remaining (a,b). An iterative flood
    fill visits every sample once, using O(N) time and storage. A boundary
    sample has an eight-neighbor sample of another status or touches a scan
    edge. Bounds are inclusive (first column, last column, first row, last row).

    :param statuses: Complete status grid, including all tiles.
    :return: Per-sample positive domain IDs and summaries in ID order.
    :raises ValueError: If the status grid is empty or not two-dimensional.

    Doctests:
    >>> diagonal = np.array([[0, 1], [1, 0]], dtype=np.int32)
    >>> labels, domains = label_domains(diagonal)
    >>> labels.tolist()
    [[1, 2], [2, 1]]
    >>> [(d.status, d.photon_count, d.seed) for d in domains]
    [(0, 2, (0, 0)), (1, 2, (0, 1))]
    >>> separated = np.array([[0, 1, 0], [0, 1, 0]], dtype=np.int32)
    >>> labels, domains = label_domains(separated)
    >>> labels.tolist()
    [[1, 2, 3], [1, 2, 3]]
    >>> sum(d.photon_count for d in domains) == separated.size
    True
    >>> label_domains(np.array([[0, 0, 2], [1, 1, 2]], dtype=np.int32))[0].tolist()
    [[1, 1, 3], [2, 2, 3]]
    >>> ring = np.zeros((5, 5), dtype=np.int32)
    >>> ring[1:4, 1:4] = 1
    >>> ring[2, 2] = 0
    >>> labels, domains = label_domains(ring)
    >>> labels.tolist()
    [[1, 1, 1, 1, 1], [1, 2, 2, 2, 1], [1, 2, 3, 2, 1], [1, 2, 2, 2, 1], [1, 1, 1, 1, 1]]
    >>> [(d.photon_count, d.boundary_count, d.edges) for d in domains]
    [(16, 16, ('left', 'right', 'bottom', 'top')), (8, 8, ()), (1, 1, ())]
    >>> domains[1].bounds, domains[2].seed
    ((1, 3, 1, 3), (2, 2))
    >>> label_domains(np.array([[4, 4, 9]], dtype=np.int32))[0].tolist()
    [[1, 1, 2]]
    >>> label_domains(np.array([[6]], dtype=np.int32))[1][0].boundary_count
    1
    >>> label_domains(np.zeros((5, 5), dtype=np.int32))[1][0].boundary_count
    16
    >>> label_domains(np.zeros((0, 2), dtype=np.int32))
    Traceback (most recent call last):
    ...
    ValueError: Domain labeling requires a nonempty two-dimensional status grid
    """
    if statuses.ndim != 2 or statuses.size == 0:
        raise ValueError(
            "Domain labeling requires a nonempty two-dimensional status grid"
        )
    height, width = statuses.shape
    labels = np.zeros(statuses.shape, dtype=np.intp)
    domains: List[BlueprintDomain] = []
    # Step 1: Start each flood fill at the first remaining sample in (a,b) order.
    for seed_col in range(width):
        for seed_row in range(height):
            if labels[seed_row, seed_col]:
                continue
            domain_id = len(domains) + 1
            status = int(statuses[seed_row, seed_col])
            labels[seed_row, seed_col] = domain_id
            pending = [(seed_row, seed_col)]
            first_col = last_col = seed_col
            first_row = last_row = seed_row
            count = boundary_count = 0
            # Step 2: Assign before queuing so each sample is processed once.
            while pending:
                row, col = pending.pop()
                count += 1
                first_col, last_col = min(first_col, col), max(last_col, col)
                first_row, last_row = min(first_row, row), max(last_row, row)
                boundary = row in (0, height - 1) or col in (0, width - 1)
                for near_row in range(max(0, row - 1), min(height, row + 2)):
                    for near_col in range(max(0, col - 1), min(width, col + 2)):
                        if statuses[near_row, near_col] != status:
                            boundary = True
                        elif labels[near_row, near_col] == 0:
                            labels[near_row, near_col] = domain_id
                            pending.append((near_row, near_col))
                boundary_count += int(boundary)
            edges = tuple(
                name
                for name, touched in (
                    ("left", first_col == 0),
                    ("right", last_col == width - 1),
                    ("bottom", first_row == 0),
                    ("top", last_row == height - 1),
                )
                if touched
            )
            domains.append(
                BlueprintDomain(
                    domain_id,
                    status,
                    count,
                    (seed_col, seed_row),
                    (first_col, last_col, first_row, last_row),
                    boundary_count,
                    edges,
                )
            )
    return labels, domains


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


def _calculate_subsample_rate(total_expected_records: int, max_viz_points: int) -> int:
    """
    Calculate a visualization subsampling rate.

    :param total_expected_records: Number of records across all tiles.
    :param max_viz_points: Maximum number of points to retain.
    :return: Subsampling rate, always at least one.
    """
    return max(1, (total_expected_records + max_viz_points - 1) // max_viz_points)


def plot_domains(
    labels: npt.NDArray[np.intp],
    domains: List[BlueprintDomain],
    a_values: npt.NDArray[np.float64],
    b_values: npt.NDArray[np.float64],
    scan_bounds: Tuple[float, float, float, float],
) -> Figure:
    """
    Display domains at native sample resolution with click and ID selection.

    The standard Matplotlib toolbar provides pan, zoom, and Home. Clicking a
    sample or entering its domain ID zooms to that domain's enclosing rectangle.
    Colors cycle after twenty domains; selection highlights the chosen domain
    in black and reports its exact ID.
    Nearest-neighbor drawing preserves cells and holes when zooming.

    :param labels: Per-sample domain IDs with rows increasing with b.
    :param domains: Domain summaries in ID order.
    :param a_values: Column sample centers in the original camera coordinates.
    :param b_values: Row sample centers in the original camera coordinates.
    :param scan_bounds: Exact (a_min, a_max, b_min, b_max) from the header.
    :return: Domain figure, with callbacks retained for interactive use.

    Doctests:
    A simulated click on the upper-right sample selects the correct domain
    and zooms in the original cropped camera coordinates, without a GUI.
    >>> from matplotlib.backends.backend_agg import FigureCanvasAgg
    >>> from unittest.mock import patch
    >>> labels, domains = label_domains(np.array([[0, 0], [0, 1]], dtype=np.int32))
    >>> figure = Figure()
    >>> canvas = FigureCanvasAgg(figure)
    >>> axes = figure.subplots()
    >>> with patch.object(plt, "subplots", return_value=(figure, axes)), patch.object(plt, "show"):
    ...     figure = plot_domains(labels, domains, np.array([0.375, 0.625]),
    ...         np.array([0.6875, 0.8125]), (0.25, 0.75, 0.625, 0.875))
    >>> canvas.draw()
    >>> x, y = axes.transData.transform((0.625, 0.8125))
    >>> canvas.callbacks.process("button_press_event", MouseEvent("button_press_event", canvas, x, y, button=1))
    >>> figure.texts[0].get_text().splitlines()[0]
    'Domain 2: STOP_CONDITION_TERMINAL_PLANE (1)'
    >>> np.allclose(axes.get_xlim(), (0.375, 0.875)) and np.allclose(axes.get_ylim(), (0.6875, 0.9375))
    True
    >>> axes.images[0].origin, axes.images[0].get_interpolation()
    ('lower', 'nearest')
    >>> np.allclose(axes.images[0].cmap(axes.images[0].norm(2)), (0, 0, 0, 1))
    True
    >>> _ = axes.set_xlim(0.25, 0.75)
    >>> _ = axes.set_ylim(0.625, 0.875)
    >>> canvas.draw()
    >>> x, y = axes.transData.transform((0.625, 0.8125))
    >>> canvas.callbacks.process("button_press_event", MouseEvent("button_press_event", canvas, x, y, button=1))
    >>> np.allclose(axes.get_xlim(), (0.375, 0.875))
    True
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    fig.subplots_adjust(bottom=0.17, top=0.82)
    colors = plt.get_cmap("tab20")(np.arange(len(domains)) % 20)
    image = ax.imshow(
        labels,
        origin="lower",
        extent=scan_bounds,
        interpolation="nearest",
        cmap=ListedColormap(colors),
        vmin=0.5,
        vmax=len(domains) + 0.5,
        aspect="equal",
    )
    ax.set_xlabel("a (dimensionless)")
    ax.set_ylabel("b (dimensionless)")
    ax.set_title(f"{len(domains)} domains; eight-neighbor connectivity")
    information = fig.text(
        0.05,
        0.96,
        "Click a domain or enter its ID to zoom. Toolbar: pan, zoom, Home.\n"
        "Colors cycle after 20 domains; the selected domain is black.",
        va="top",
    )
    (marker,) = ax.plot([], [], "wx", markersize=10)
    domain_box = TextBox(fig.add_axes((0.25, 0.05, 0.25, 0.05)), "Domain ID: ")
    da = (scan_bounds[1] - scan_bounds[0]) / labels.shape[1]
    db = (scan_bounds[3] - scan_bounds[2]) / labels.shape[0]

    def select_domain(text: str) -> None:
        """
        Show the requested domain's status and zoom to its sample cells.

        :param text: Domain ID entered in the text box.
        """
        try:
            domain_id = int(text)
        except ValueError:
            information.set_text(
                f"Enter an integer domain ID from 1 to {len(domains)}."
            )
            fig.canvas.draw_idle()
            return
        if not 1 <= domain_id <= len(domains):
            information.set_text(f"Enter a domain ID from 1 to {len(domains)}.")
            fig.canvas.draw_idle()
            return
        domain = domains[domain_id - 1]
        selected_colors = colors.copy()
        selected_colors[domain_id - 1] = (0, 0, 0, 1)
        image.set_cmap(ListedColormap(selected_colors))
        status_name, _ = _termination_type_info(domain.status)
        information.set_text(
            f"Domain {domain_id}: {status_name} ({domain.status})\n"
            f"{domain.photon_count:,} photons ({100 * domain.photon_count / labels.size:.6g}%); "
            f"{domain.boundary_count:,} boundary samples; "
            f"scan edges: {', '.join(domain.edges) or 'none'}"
        )
        col, row = domain.seed
        marker.set_data([a_values[col]], [b_values[row]])
        first_col, last_col, first_row, last_row = domain.bounds
        toolbar = fig.canvas.toolbar
        if toolbar is not None:
            toolbar.push_current()
        ax.set_xlim(a_values[first_col] - da, a_values[last_col] + da)
        ax.set_ylim(b_values[first_row] - db, b_values[last_row] + db)
        if toolbar is not None:
            toolbar.push_current()
        fig.canvas.draw_idle()

    def on_click(event: MouseEvent) -> None:
        """
        Select an image sample while leaving toolbar gestures available.

        :param event: Matplotlib mouse-button event.
        """
        if event.inaxes is ax and event.button == 1 and ax.get_navigate_mode() is None:
            domain_id = image.get_cursor_data(event)
            if domain_id is not None:
                text = str(int(domain_id))
                if domain_box.text == text:
                    select_domain(text)
                else:
                    domain_box.set_val(text)

    domain_box.on_submit(select_domain)
    # The canvas retains on_click, whose closure keeps the TextBox alive.
    fig.canvas.mpl_connect("button_press_event", on_click)
    plt.show()
    return fig


def plot_heatmaps(
    data: npt.NDArray[np.void],
    nonterminal_crossing_mask: Optional[npt.NDArray[np.bool_]] = None,
) -> None:
    """
    Generate diagnostic heatmaps for the three recorded event classes.

    Nonterminal coordinates are plotted only for records that crossed the
    nonterminal plane. Terminal coordinates are plotted only for records whose
    final status is ``STOP_CONDITION_TERMINAL_PLANE``. These plots never place pixels;
    image placement belongs to the renderer and uses normalized image-sample
    fractions.

    :param data: Structured blueprint records to plot.
    :param nonterminal_crossing_mask: Optional mask of records with a
        nonterminal-plane crossing.
    :raises ValueError: If the crossing mask and blueprint records have
        different shapes.
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    # A nonzero finite lambda marks a recorded nonterminal crossing. Its
    # presence does not depend on the ray's final termination status.
    nonterminal_mask = (
        cast(
            npt.NDArray[np.bool_],
            (
                np.isfinite(data["non_terminal_plane_lambda"])
                & (data["non_terminal_plane_lambda"] != 0.0)
                & np.isfinite(data["y_nt"])
                & np.isfinite(data["z_nt"])
            ),
        )
        if nonterminal_crossing_mask is None
        else nonterminal_crossing_mask
    )
    if nonterminal_mask.shape != data.shape:
        raise ValueError("Nonterminal crossing mask does not match blueprint records")
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
    plt.figure(figsize=(10, 6))
    termination_colors = plt.get_cmap("tab10")
    for enum_value in np.unique(termination_types):
        values = norm_abs_values[termination_types == enum_value]
        positive_values = values[values > 0.0]
        if len(positive_values) > 0:
            enum_int = int(enum_value)
            termination_name, _ = _termination_type_info(enum_int)
            plt.hist(
                np.log10(positive_values),
                bins=100,
                alpha=0.5,
                color=termination_colors(enum_int % 10),
                label=f"Type {enum_int}: {termination_name}",
            )
    plt.title("Log-Scale Normalization Magnitude by Termination Type")
    plt.xlabel("$\\log_{10}(|g_{{\\mu\\nu}}p^\\mu p^\\nu|)$")
    plt.yscale("log", base=10)
    plt.ylabel("$\\log_{10}(N_{\\mathrm{photons}})$")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.show()


def diagnose_blueprint(show_plots: bool = True) -> None:
    """
    Read every metadata-described tile and print integrity/physics diagnostics.

    Write blueprint_domains.txt next to the input tiles. Native records are
    ordered by column within each row, and every tile has the same dimensions.
    Domain labeling retains a status grid and a label grid for every photon;
    only the existing physical diagnostic plots are subsampled.

    :param show_plots: Display domain, plane, sphere, and normalization figures.
    :raises FileNotFoundError: If an expected tile is missing.
    :raises ValueError: If a blueprint binary file or normalization sidecar is malformed.

    Doctests:
    A diagonal pair across both tile seams stays one domain even when chunks
    split sample rows. The report accounts for every photon.
    >>> import contextlib
    >>> import csv
    >>> import io
    >>> import struct
    >>> import tempfile
    >>> from pathlib import Path
    >>> from unittest.mock import patch
    >>> grid = np.zeros((4, 4), dtype=np.int32)
    >>> grid[1, 1] = grid[2, 2] = cfg.STOP_CONDITION_TERMINAL_PLANE
    >>> with tempfile.TemporaryDirectory() as directory:
    ...     for tx in range(2):
    ...         for ty in range(2):
    ...             records = np.zeros(4, dtype=cfg.BLUEPRINT_DTYPE)
    ...             records["termination_type"] = grid[2*ty:2*ty+2, 2*tx:2*tx+2].ravel()
    ...             records["image_width_fraction"] = np.tile((2*tx + np.arange(2) + 0.5) / 4, 2)
    ...             records["image_height_fraction"] = np.repeat((2*ty + np.arange(2) + 0.5) / 4, 2)
    ...             header = struct.pack(cfg.BLUEPRINT_HEADER_FORMAT,
    ...                 cfg.BLUEPRINT_MAGIC, cfg.BLUEPRINT_SCHEMA_VERSION,
    ...                 cfg.BLUEPRINT_HEADER_SIZE, cfg.BLUEPRINT_RECORD_SIZE,
    ...                 tx, ty, 2, 2, 4, 1.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    ...             path = Path(directory) / f"light_blueprint_{tx:02d}_{ty:02d}.bin"
    ...             _ = path.write_bytes(header + records.tobytes())
    ...     with patch.dict(diagnose_blueprint.__globals__, script_dir=directory):
    ...         with patch.object(cfg, "CHUNK_SIZE", 3), contextlib.redirect_stdout(io.StringIO()):
    ...             diagnose_blueprint(show_plots=False)
    ...     lines = (Path(directory) / "blueprint_domains.txt").read_text().splitlines()
    ...     rows = list(csv.DictReader((line for line in lines if not line.startswith("#")), delimiter=chr(9)))
    >>> [(int(r["domain_id"]), int(r["status"]), int(r["photon_count"])) for r in rows]
    [(1, 0, 14), (2, 1, 2)]
    >>> [float(r["percentage"]) for r in rows]
    [87.5, 12.5]
    >>> [float(rows[1][key]) for key in ("seed_a", "seed_b", "a_min_sample", "a_max_sample", "b_min_sample", "b_max_sample")]
    [0.375, 0.375, 0.375, 0.625, 0.375, 0.625]
    >>> rows[1]["status_name"], rows[1]["scan_edges"], int(rows[1]["boundary_count"])
    ('STOP_CONDITION_TERMINAL_PLANE', 'none', 2)
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
    print(f" a range: [{first_header.a_min!r}, {first_header.a_max!r}]")
    print(f" b range: [{first_header.b_min!r}, {first_header.b_max!r}]")
    print("=================================================================")

    blueprint_files: List[str] = []
    headers: List[blueprint_io.BlueprintHeader] = []
    crossing_indices_by_file: Dict[str, npt.NDArray[np.uint64]] = {}
    missing_sidecars: List[str] = []
    total_expected_records = 0
    # Step 1: Recover native tile dimensions from the repeated first-row b value.
    _, _, first_chunk = next(
        blueprint_io.iter_blueprint_chunks(
            first_tile,
            1,
            allow_active_termination=True,
        )
    )
    grid_records = np.memmap(
        first_tile,
        dtype=cfg.BLUEPRINT_DTYPE,
        mode="r",
        offset=first_header.header_size,
        shape=(first_header.record_count,),
    )
    tile_width = int(
        np.count_nonzero(
            grid_records["image_height_fraction"]
            == first_chunk["image_height_fraction"][0]
        )
    )
    tile_height, remainder = divmod(first_header.record_count, tile_width)
    if remainder:
        raise ValueError("Blueprint record count does not form complete sample rows")
    del grid_records
    statuses = np.empty(
        (tile_height * tiles_height, tile_width * tiles_width), dtype=np.int32
    )
    a_values = np.empty(tile_width * tiles_width, dtype=np.float64)
    b_values = np.empty(tile_height * tiles_height, dtype=np.float64)
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
            if header.record_count != first_header.record_count:
                raise ValueError(
                    f"Blueprint '{filename}' has an unexpected tile record count"
                )
            grid_records = np.memmap(
                filepath,
                dtype=cfg.BLUEPRINT_DTYPE,
                mode="r",
                offset=header.header_size,
                shape=(tile_height, tile_width),
            )
            a_values[tile_x * tile_width : (tile_x + 1) * tile_width] = grid_records[
                "image_width_fraction"
            ][0, :]
            b_values[tile_y * tile_height : (tile_y + 1) * tile_height] = grid_records[
                "image_height_fraction"
            ][:, 0]
            del grid_records
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

            crossing_name = cfg.NON_TERMINAL_CROSSINGS_FILENAME_TEMPLATE.format(
                tile_x=tile_x, tile_y=tile_y
            )
            crossing_path = os.path.join(script_dir, crossing_name)
            if os.path.isfile(crossing_path):
                complete_records = (
                    os.path.getsize(crossing_path) // cfg.PLANE_CROSSING_RECORD_SIZE
                )
                crossing_records = np.fromfile(
                    crossing_path,
                    dtype=cfg.PLANE_CROSSING_DTYPE,
                    count=complete_records,
                )
                indices = crossing_records["photon_index"]
                crossing_indices_by_file[filepath] = np.unique(
                    indices[indices < header.record_count]
                )
            else:
                crossing_indices_by_file[filepath] = np.empty(0, dtype=np.uint64)

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
    horizon_counts = np.zeros(4, dtype=np.int64)
    first_records: List[npt.NDArray[np.void]] = []
    sampled_records: List[npt.NDArray[np.void]] = []
    sampled_norm_abs: List[npt.NDArray[np.float64]] = []
    sampled_nonterminal_crossings: List[npt.NDArray[np.bool_]] = []
    for filepath, header in zip(blueprint_files, headers):
        tile_statuses = statuses[
            header.tile_y * tile_height : (header.tile_y + 1) * tile_height,
            header.tile_x * tile_width : (header.tile_x + 1) * tile_width,
        ]
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
            crossing_indices = crossing_indices_by_file[filepath]
            for _, record_start, chunk_data in blueprint_io.iter_blueprint_chunks(
                filepath,
                cfg.CHUNK_SIZE,
                allow_active_termination=True,
            ):
                tile_statuses.flat[record_start : record_start + len(chunk_data)] = (
                    chunk_data["termination_type"]
                )
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
                horizon_stops = (
                    chunk_data["termination_type"]
                    == cfg.STOP_CONDITION_INSIDE_APPARENT_HORIZON
                )
                if np.any(horizon_stops):
                    horizon_ids = chunk_data["y_t"][horizon_stops].astype(np.intp)
                    horizon_counts += np.bincount(
                        horizon_ids, minlength=len(horizon_counts)
                    )
                total_rays += len(chunk_data)
                crossing_start = np.searchsorted(crossing_indices, record_start)
                crossing_end = np.searchsorted(
                    crossing_indices, record_start + len(chunk_data)
                )
                local_crossings = crossing_indices[crossing_start:crossing_end]
                nonterminal_mask = np.zeros(len(chunk_data), dtype=bool)
                if local_crossings.size:
                    nonterminal_mask[local_crossings.astype(np.intp) - record_start] = (
                        True
                    )
                nonterminal_mask &= (
                    np.isfinite(chunk_data["y_nt"])
                    & np.isfinite(chunk_data["z_nt"])
                    & np.isfinite(chunk_data["non_terminal_plane_t"])
                )
                in_view += int(np.sum(nonterminal_mask))
                if len(first_records) < 9:
                    first_records.extend(chunk_data[: 9 - len(first_records)])
                sampled_records.append(chunk_data[::subsample_rate])
                sampled_nonterminal_crossings.append(nonterminal_mask[::subsample_rate])
                if norm_abs_chunk is not None:
                    sampled_norm_abs.append(norm_abs_chunk[::subsample_rate])
            if sidecar_file is not None and sidecar_file.read(1):
                raise ValueError(
                    f"Normalization sidecar for '{filepath}' has extra data"
                )

    print("--- Raw Termination Enums in Binary ---")
    for enum_value in sorted(enum_counts):
        count = enum_counts[enum_value]
        name, _ = _termination_type_info(enum_value)
        percentage = count / total_rays * 100.0 if total_rays else 0.0
        print(
            f"  Raw Enum {enum_value:2d} ({name}): "
            f"{count:12,} rays ({percentage:6.2f}%)"
        )

    if np.any(horizon_counts):
        print("\n--- Containing Apparent Horizons ---")
        for horizon_id in (cfg.PUNCTURE_1, cfg.PUNCTURE_2, cfg.MERGED_PUNCTURE):
            label = (
                "common horizon"
                if horizon_id == cfg.MERGED_PUNCTURE
                else f"horizon {horizon_id}"
            )
            print(f"  {label}: {horizon_counts[horizon_id]:,} rays")

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

    # Step 2: Label all samples before making any subsampled diagnostic plots.
    print("\nLabeling eight-neighbor domains using every photon...")
    labels, domains = label_domains(statuses)
    scan_bounds = (
        first_header.a_min,
        first_header.a_max,
        first_header.b_min,
        first_header.b_max,
    )
    report_path = os.path.join(script_dir, "blueprint_domains.txt")
    with open(report_path, "w", encoding="utf-8") as report:
        report.write(
            "# Photon domains: equal termination status, eight-neighbor connectivity\n"
            f"# Scan bounds (a_min a_max b_min b_max): {' '.join(repr(v) for v in scan_bounds)}\n"
            f"# Grid: {len(a_values)} columns x {len(b_values)} rows; {labels.size} photons; {len(domains)} domains\n"
            "# IDs follow increasing seed a, then b. Each photon belongs to one domain.\n"
            "# Bounds below enclose sample centers; they are not domain outlines.\n"
            "# Boundary samples have a different-status eight-neighbor or touch a scan edge.\n"
            "# ACTIVE domains, if present, contain unfinished rays.\n"
            "domain_id\tstatus\tstatus_name\tphoton_count\tpercentage\tseed_a\tseed_b\t"
            "a_min_sample\ta_max_sample\tb_min_sample\tb_max_sample\tboundary_count\tscan_edges\n"
        )
        for domain in domains:
            col, row = domain.seed
            first_col, last_col, first_row, last_row = domain.bounds
            status_name, _ = _termination_type_info(domain.status)
            coordinates = (
                a_values[col],
                b_values[row],
                a_values[first_col],
                a_values[last_col],
                b_values[first_row],
                b_values[last_row],
            )
            report.write(
                f"{domain.domain_id}\t{domain.status}\t{status_name}\t{domain.photon_count}\t"
                f"{100 * domain.photon_count / labels.size:.9g}\t"
                + "\t".join(repr(float(value)) for value in coordinates)
                + f"\t{domain.boundary_count}\t{','.join(domain.edges) or 'none'}\n"
            )
    print(f"{len(domains):,} domains; report written to {report_path}")

    if show_plots:
        plot_domains(labels, domains, a_values, b_values, scan_bounds)
    if show_plots and sampled_records:
        viz_data = np.concatenate(sampled_records)
        plot_heatmaps(viz_data, np.concatenate(sampled_nonterminal_crossings))
        if have_sidecars:
            plot_norm_abs_log_histogram(
                viz_data["termination_type"], np.concatenate(sampled_norm_abs)
            )


if __name__ == "__main__":
    import argparse
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed:
        sys.exit(1)
    print(f"Doctest passed: All {results.attempted} test(s) passed")

    parser = argparse.ArgumentParser(
        description=(
            "Diagnose blueprint binary files using tile counts and fields of view "
            "stored in their headers."
        )
    )
    parser.add_argument(
        "--no-plots",
        action="store_true",
        help="Write the domain report without opening figures.",
    )
    args = parser.parse_args()

    first_blueprint_tile = os.path.join(script_dir, "light_blueprint_00_00.bin")
    if not os.path.exists(first_blueprint_tile):
        print(
            f"No native blueprint binary files found in {script_dir}; nothing to diagnose."
        )
    else:
        diagnose_blueprint(show_plots=not args.no_plots)

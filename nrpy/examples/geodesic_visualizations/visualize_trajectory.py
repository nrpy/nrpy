# nrpy/examples/geodesic_visualizations/visualize_trajectory.py
"""
Defines the single-ray trajectory visualization suite.

This script parses the trajectory text file to ensure numerical stability and validates
the structural integrity of the output. It creates a 3D plot of the integrated
geodesic relative to a schematic black-hole reference surface.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import logging
import os
import sys
from typing import Optional, cast

import numpy as np
import numpy.typing as npt

# Mute matplotlib's verbose font debugging output at the module level
# to prevent console spam during standard execution.
logging.getLogger("matplotlib.font_manager").setLevel(logging.WARNING)
logging.getLogger("matplotlib").setLevel(logging.WARNING)
logging.getLogger("PIL").setLevel(logging.WARNING)


def plot_trajectory(
    data: "npt.NDArray[np.float64]",
    r_horizon: float = 2.0,
    particle_type: str = "Test Particle",
    plot_norm_error: bool = False,
) -> None:
    """
    Create a 3D visualization of the particle trajectory and reference sphere.

    This function plots the spatial coordinates extracted from the integration state
    vector and superimposes a spherical reference surface at the requested radius.

    :param data: The parsed 2D NumPy array containing the trajectory metrics.
    :param r_horizon: The radial coordinate used for the reference sphere.
    :param particle_type: String descriptor of the particle for plot labels.
    :param plot_norm_error: Whether to color the path by logarithmic norm error.
    :raises ValueError: If norm-error plotting receives missing or non-finite values.
    """
    # pylint: disable=import-outside-toplevel, import-error, no-name-in-module
    import matplotlib.pyplot as plt  # type: ignore[import-not-found, import-untyped, unused-ignore]
    from matplotlib.colors import (
        Normalize,  # type: ignore[import-not-found, import-untyped, unused-ignore]
    )
    from mpl_toolkits.mplot3d.art3d import (  # type: ignore[import-not-found, import-untyped, unused-ignore]
        Line3DCollection,
    )
    from mpl_toolkits.mplot3d.axes3d import (  # type: ignore[import-not-found, import-untyped, unused-ignore]
        Axes3D,
    )

    # Step 1: Descriptive Physical Variable Mapping.
    # Extract the spatial coordinates (x, y, z) from the dataset.
    # Based on the C output format: col 2 is x, col 3 is y, col 4 is z.
    x_pts = data[:, 2]
    y_pts = data[:, 3]
    z_pts = data[:, 4]

    if plot_norm_error:
        if data.shape[1] < 11:
            raise ValueError(
                "norm-error plotting requires trajectory data with at least "
                "11 columns."
            )
        norm_error = data[:, 10]
        if not np.all(np.isfinite(norm_error)):
            raise ValueError("norm-error column must contain only finite values.")

        norm_error_magnitude = np.abs(norm_error)
        positive_magnitudes = norm_error_magnitude[norm_error_magnitude > 0.0]
        error_floor = (
            np.min(positive_magnitudes)
            if positive_magnitudes.size > 0
            else sys.float_info.epsilon
        )
        log_norm_error = np.log10(np.maximum(norm_error_magnitude, error_floor))
        color_min = float(np.min(log_norm_error))
        color_max = float(np.max(log_norm_error))
        if color_min == color_max:
            color_min -= 0.5
            color_max += 0.5
        color_norm = Normalize(vmin=color_min, vmax=color_max)
        points = np.column_stack((x_pts, y_pts, z_pts))

    fig = plt.figure(figsize=(8, 8))
    ax = cast(Axes3D, fig.add_subplot(111, projection="3d"))

    # Step 2: Plot the geodesic path.
    if plot_norm_error:
        if len(points) > 1:
            segments = np.stack((points[:-1], points[1:]), axis=1)
            path_collection = Line3DCollection(
                segments,
                cmap="viridis",
                norm=color_norm,
                linewidths=1.5,
            )
            # Segment i-1 to i represents norm error measured at accepted step i.
            path_collection.set_array(log_norm_error[1:])
            path_collection.set_label(f"{particle_type.capitalize()} Trajectory")
            ax.add_collection3d(path_collection)
            color_mapper = path_collection
        else:
            color_mapper = plt.cm.ScalarMappable(
                norm=color_norm,
                cmap="viridis",
            )
            color_mapper.set_array(log_norm_error)
            ax.scatter(
                x_pts,
                y_pts,
                z_pts,
                c=log_norm_error,
                cmap="viridis",
                norm=color_norm,
                s=20,
                depthshade=False,
                label=f"{particle_type.capitalize()} Trajectory",
            )
        ax.auto_scale_xyz(x_pts, y_pts, z_pts)
        colorbar = fig.colorbar(color_mapper, ax=ax, pad=0.1)
        colorbar.set_label(r"$\log_{10}(|\mathrm{norm\ error}|)$")
    else:
        ax.plot(
            x_pts,
            y_pts,
            z_pts,
            label=f"{particle_type.capitalize()} Trajectory",
            color="blue",
            linewidth=1.5,
        )

    # Step 3: Mark Integration Boundaries.
    # Drop markers for the simulation start and endpoints to easily verify
    # initial conditions and termination locations.
    ax.scatter(
        x_pts[0], y_pts[0], z_pts[0], color="green", marker="o", s=50, label="Start"
    )
    ax.scatter(
        x_pts[-1], y_pts[-1], z_pts[-1], color="red", marker="x", s=50, label="End"
    )

    # Step 4: Construct the reference sphere surface.
    # Parameterize angles for creating the spherical reference surface.
    # u is the azimuthal angle [0, 2pi], v is the polar angle [0, pi].
    u_val = np.linspace(0, 2 * np.pi, 20)
    v_val = np.linspace(0, np.pi, 10)

    # Create 2D meshes mapping the spherical coordinates.
    u, v = np.meshgrid(u_val, v_val, indexing="ij")

    # Convert to Cartesian coordinates to plot the reference sphere.
    xh = r_horizon * np.cos(u) * np.sin(v)
    yh = r_horizon * np.sin(u) * np.sin(v)
    zh = r_horizon * np.cos(v)
    ax.plot_surface(
        xh, yh, zh, color="black", alpha=0.3, label="Reference sphere (r=2M)"
    )

    # Step 5: Format and Display.
    ax.set_xlabel("x (M)")
    ax.set_ylabel("y (M)")
    ax.set_zlabel("z (M)")
    ax.set_aspect("equal")
    title = f"{particle_type.capitalize()} Geodesic in Kerr-Schild Cartesian Spacetime"
    if plot_norm_error:
        title += " (colored by norm error)"
    ax.set_title(title)
    ax.legend()

    # Visual disclaimer about the reference surface approximation.
    fig.text(
        0.5,
        0.02,
        "* Reference sphere rendered at r=2M; not the exact Kerr horizon for nonzero spin.",
        ha="center",
        fontsize=9,
        color="gray",
    )

    # Location to save the rendered matplotlib figure.
    plot_filename = (
        "trajectory_norm_plot.png" if plot_norm_error else "trajectory_plot.png"
    )
    plot_path = os.path.abspath(plot_filename)
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"\n[i] Visualization successfully saved to:\n    {plot_path}")

    print("[i] Displaying 3D plot. Close the window to continue...")
    plt.show()


def visualize_trajectory(
    traj_path: Optional[str] = None,
    particle_type: str = "Test Particle",
    plot_norm_error: bool = False,
) -> None:
    """
    Read the trajectory data file and orchestrate the diagnostic visualization.

    This function validates the existence and integrity of the output file
    before passing the parsed data to the plotting routine.

    :param traj_path: Path to the .txt trajectory file. Defaults to 'trajectory.txt'.
    :param particle_type: String representing the particle type for dynamic labeling.
    :param plot_norm_error: Whether to color the path by logarithmic norm error.
    """
    if traj_path is None:
        # Default expectation: the script is executed in the output directory.
        traj_path = "trajectory.txt"

    print("=================================================================")
    print(" TRAJECTORY DIAGNOSTICS & VISUALIZATION")
    print(f" File: {os.path.abspath(traj_path)}")
    print(f" Particle Type: {particle_type.capitalize()}")
    print("=================================================================")

    if not os.path.exists(traj_path):
        print(
            f"[!] ERROR: Trajectory file not found at:\n    {os.path.abspath(traj_path)}"
        )
        print("\n    Did you compile and run the C executable first?")
        print("    1. Type `make` in this directory.")
        print("    2. Run the executable generated for this project.")
        return

    # Step 1: Load and Validate Data.
    try:
        # Parse the trajectory metrics into a 2D NumPy array, ignoring comments.
        data = np.atleast_2d(np.loadtxt(traj_path, comments="#"))

        if data.size == 0:
            print(f"[!] ERROR: '{traj_path}' is empty.")
            print(
                "    Integration may have failed immediately. Check initial conditions."
            )
            return

        if data.shape[1] < 5:
            print(
                f"[!] ERROR: '{traj_path}' has {data.shape[1]} columns; "
                "at least 5 are required for trajectory visualization."
            )
            return

        if plot_norm_error and data.shape[1] < 11:
            print(
                f"[!] ERROR: '{traj_path}' has {data.shape[1]} columns; "
                "norm-error plotting requires at least 11 columns."
            )
            return

        if plot_norm_error and not np.all(np.isfinite(data[:, 10])):
            print(f"[!] ERROR: '{traj_path}' has non-finite norm-error values.")
            return

        total_steps = len(data)
        print(f"Total integration steps loaded: {total_steps:,}\n")

    except (ValueError, OSError) as e:
        print(f"[!] ERROR: Failed to parse trajectory data: {e}")
        return

    # Step 2: Check Trajectory Bounds (Diagnostics).
    # Direct sanity check of the initial and final states.
    print("Boundary State Diagnostics.")
    print(
        f"  Initial Pos (x, y, z): ({data[0, 2]:>8.4f}, {data[0, 3]:>8.4f}, "
        f"{data[0, 4]:>8.4f})"
    )
    print(
        f"  Final Pos   (x, y, z): ({data[-1, 2]:>8.4f}, {data[-1, 3]:>8.4f}, "
        f"{data[-1, 4]:>8.4f})"
    )
    print(f"  Total Path Parameter: {data[-1, 0]:>8.4f}")

    # Step 3: Launch Visualization.
    print("\n[i] Creating 3D representation...")
    plot_trajectory(
        data,
        r_horizon=2.0,
        particle_type=particle_type,
        plot_norm_error=plot_norm_error,
    )
    print("=================================================================")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Standalone Geodesic Trajectory Visualizer"
    )
    parser.add_argument(
        "--traj_path",
        type=str,
        default="trajectory.txt",
        help="Path to the numerical trajectory output file.",
    )
    parser.add_argument(
        "--particle_type",
        type=str,
        default="Test Particle",
        help="Type of particle integrated. Affects plot labels.",
    )
    parser.add_argument(
        "--plot_norm_error",
        action="store_true",
        help="Color the 3D trajectory by log10 absolute norm error.",
    )

    args = parser.parse_args()

    visualize_trajectory(
        traj_path=args.traj_path,
        particle_type=args.particle_type,
        plot_norm_error=args.plot_norm_error,
    )

# nrpy/examples/geodesic_visualizations/visualize_trajectory.py
"""
Defines the single-ray trajectory visualization suite.

This script parses the trajectory text file to ensure numerical stability and validates
the structural integrity of the output. It creates a 3D plot of the particle's path
interacting with the event horizon to visually confirm gravitational lensing effects.

Author: Dalton J. Moone
"""

import argparse
import logging
import os
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
) -> None:
    """
    Create a 3D visualization of the particle trajectory and the black hole horizon.

    This function plots the spatial coordinates extracted from the integration state
    vector and superimposes a spherical representation of the event horizon.

    :param data: The parsed 2D NumPy array containing the trajectory metrics.
    :param r_horizon: The radial coordinate representing the event horizon.
    :param particle_type: String descriptor of the particle for plot labels.
    """
    # pylint: disable=import-outside-toplevel, import-error, no-name-in-module
    import matplotlib.pyplot as plt  # type: ignore[import-not-found, import-untyped, unused-ignore]
    from mpl_toolkits.mplot3d.axes3d import (  # type: ignore[import-not-found, import-untyped, unused-ignore]
        Axes3D,
    )

    # Step 1: Descriptive Physical Variable Mapping.
    # Extract the spatial coordinates (x, y, z) from the dataset.
    # Based on the C output format: col 2 is x, col 3 is y, col 4 is z.
    x_pts = data[:, 2]
    y_pts = data[:, 3]
    z_pts = data[:, 4]

    fig = plt.figure(figsize=(8, 8))
    ax = cast(Axes3D, fig.add_subplot(111, projection="3d"))

    # Step 2: Plot the geodesic path.
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
    # initial conditions and termination locations (e.g., horizon intercepts).
    ax.scatter(
        x_pts[0], y_pts[0], z_pts[0], color="green", marker="o", s=50, label="Start"
    )
    ax.scatter(
        x_pts[-1], y_pts[-1], z_pts[-1], color="red", marker="x", s=50, label="End"
    )

    # Step 4: Construct the Event Horizon Surface.
    # Parameterize angles for creating the spherical horizon surface.
    # u is the azimuthal angle [0, 2pi], v is the polar angle [0, pi].
    u_val = np.linspace(0, 2 * np.pi, 20)
    v_val = np.linspace(0, np.pi, 10)

    # Create 2D meshes mapping the spherical coordinates.
    u, v = np.meshgrid(u_val, v_val, indexing="ij")

    # Convert to Cartesian coordinates to plot the black hole.
    xh = r_horizon * np.cos(u) * np.sin(v)
    yh = r_horizon * np.sin(u) * np.sin(v)
    zh = r_horizon * np.cos(v)
    ax.plot_surface(xh, yh, zh, color="black", alpha=0.3, label="Horizon (r=2M)")

    # Step 5: Format and Display.
    ax.set_xlabel("x (M)")
    ax.set_ylabel("y (M)")
    ax.set_zlabel("z (M)")
    ax.set_aspect("equal")
    ax.set_title(
        f"{particle_type.capitalize()} Geodesic in Kerr-Schild Cartesian Spacetime"
    )
    ax.legend()

    # Visual disclaimer about the horizon approximation.
    fig.text(
        0.5,
        0.02,
        "* Horizon rendered as a sphere at r=2M. Visually approximate for a_spin = 0.",
        ha="center",
        fontsize=9,
        color="gray",
    )

    # Location to save the rendered matplotlib figure.
    plot_path = os.path.abspath("trajectory_plot.png")
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    print(f"\n[i] Visualization successfully saved to:\n    {plot_path}")

    print("[i] Displaying 3D plot. Close the window to continue...")
    plt.show()


def visualize_trajectory(
    traj_path: Optional[str] = None, particle_type: str = "Test Particle"
) -> None:
    """
    Read the trajectory data file and orchestrate the diagnostic visualization.

    This function validates the existence and integrity of the output file
    before passing the parsed data to the plotting routine.

    :param traj_path: Path to the .txt trajectory file. Defaults to 'trajectory.txt'.
    :param particle_type: String representing the particle type for dynamic labeling.
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
        print("    2. Run the executable (e.g., `./photon_geodesic_integrator`).")
        return

    # Step 1: Load and Validate Data.
    try:
        # Parse the trajectory metrics into a 2D NumPy array, ignoring comments.
        data = np.loadtxt(traj_path, comments="#")

        if data.size == 0:
            print(f"[!] ERROR: '{traj_path}' is empty.")
            print(
                "    Integration may have failed immediately. Check initial conditions."
            )
            return

        total_steps = len(data)
        print(f"Total integration steps loaded: {total_steps:,}\n")

    except (ValueError, OSError) as e:
        print(f"[!] ERROR: Failed to parse trajectory data: {e}")
        return

    # Step 2: Check Trajectory Bounds (Diagnostics).
    # Direct sanity check of the initial and final states.
    print("Step 1: Boundary State Diagnostics.")
    print(
        f"  Initial Pos (x, y, z): ({data[0, 2]:>8.4f}, {data[0, 3]:>8.4f}, "
        f"{data[0, 4]:>8.4f})"
    )
    print(
        f"  Final Pos   (x, y, z): ({data[-1, 2]:>8.4f}, {data[-1, 3]:>8.4f}, "
        f"{data[-1, 4]:>8.4f})"
    )
    print(f"  Total Affine/Proper Parameter: {data[-1, 0]:>8.4f}")

    # Step 3: Launch Visualization.
    print("\n[i] Creating 3D representation...")
    plot_trajectory(data, r_horizon=2.0, particle_type=particle_type)
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

    args = parser.parse_args()

    visualize_trajectory(traj_path=args.traj_path, particle_type=args.particle_type)

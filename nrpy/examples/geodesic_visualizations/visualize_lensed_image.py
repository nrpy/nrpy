# nrpy/examples/geodesic_visualizations/visualize_lensed_image.py
"""
Defines the visualization script for the Photon Geodesic Integrator.

This script reads binary blueprints to render the final lensed image. It maps physical
photon trajectories onto a static background texture and terminal-plane accretion-disk
geometry, using r_min and r_max to define disk-texture bounds.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import argparse
import math
import os
import sys
import urllib.request
from typing import Union

import numpy as np
import numpy.typing as npt

try:
    import blueprint_config_and_schema as cfg  # type: ignore
    import blueprint_io  # type: ignore
    import render_lensed_image as rli  # type: ignore
except ImportError:
    from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
    from nrpy.examples.geodesic_visualizations import blueprint_io
    from nrpy.examples.geodesic_visualizations import render_lensed_image as rli


def main() -> None:
    """
    Parse parameters and orchestrate lensed image creation.

    :raises FileNotFoundError: If a tile is missing from an otherwise present set.
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # The default output path is strictly local.
    default_output_path = os.path.join(script_dir, "lensed_output.png")

    # The argument parser manages the visualization configuration and physical bounds.
    parser = argparse.ArgumentParser(
        description="Create a lensed image from tiled geodesic light blueprints."
    )

    # The output argument defaults to the local file.
    parser.add_argument(
        "--output",
        type=str,
        default=default_output_path,
        help="Path where the output PNG will be saved.",
    )

    # Physical parameters govern terminal-plane accretion-disk texture geometry.
    # Image placement comes only from normalized image-sample fractions stored in
    # records; plane coordinates remain diagnostics and terminal-texture samples.
    parser.add_argument(
        "--terminal-plane-radius",
        nargs=2,
        type=float,
        default=(0.0, 1.0),
        metavar=("MIN_RADIUS", "MAX_RADIUS"),
        help="Inner and outer physical radii of the terminal-plane disk texture.",
    )
    parser.add_argument(
        "--pixel-width",
        type=int,
        default=750,
        help="Pixel width of the final output image; does not change ray-grid resolution.",
    )
    parser.add_argument(
        "--sphere-image",
        type=str,
        default=None,
        help="Optional path to a custom celestial sphere texture image.",
    )
    parser.add_argument(
        "--terminal-plane-image",
        type=str,
        default=None,
        help="Optional path to a custom terminal-plane texture image.",
    )
    parser.add_argument(
        "--debug",
        nargs="?",
        const="",
        default=None,
        metavar="key",
        help=(
            "Color evolution-measure, RKF45, maximum-time, slot, and generic failures. "
            "Pass 'key' to also add the color key to the image."
        ),
    )

    # The parsed argument namespace contains all runtime configurations.
    args = parser.parse_args()
    if args.debug not in (None, "", "key"):
        parser.error("--debug accepts no argument or the optional argument 'key'.")
    terminal_plane_min_radius, terminal_plane_max_radius = args.terminal_plane_radius
    if terminal_plane_min_radius < 0.0:
        parser.error("--terminal-plane-radius MIN_RADIUS must be nonnegative.")
    if terminal_plane_max_radius <= terminal_plane_min_radius:
        parser.error(
            "--terminal-plane-radius MAX_RADIUS must be greater than MIN_RADIUS."
        )
    if args.pixel_width < 1:
        parser.error("--pixel-width must be positive.")

    # Absolute script directory ensures paths resolve independently of execution context.
    script_dir = os.path.dirname(os.path.abspath(__file__))

    first_tile = os.path.join(script_dir, "light_blueprint_00_00.bin")
    if not os.path.exists(first_tile):
        print(
            f"No native blueprint artifacts found in {script_dir}; nothing to render."
        )
        return

    first_header = blueprint_io.read_blueprint_header(first_tile, expected_tile=(0, 0))

    # Discover tile grid and fields of view from headers. Ray samples carry
    # normalized image coordinates, so no serialized pixel-grid dimensions are
    # needed for placement.
    blueprint_files = []
    for i in range(first_header.tiles_width):
        for j in range(first_header.tiles_height):
            blueprint_filename = f"light_blueprint_{i:02d}_{j:02d}.bin"
            blueprint_filepath = os.path.join(script_dir, blueprint_filename)

            if not os.path.exists(blueprint_filepath):
                raise FileNotFoundError(
                    f"Expected tile missing at '{blueprint_filepath}'."
                )
            blueprint_files.append(blueprint_filepath)

    print(f"Loading {len(blueprint_files)} blueprint tiles...")
    print(
        "Blueprint angular aspect ratio: "
        f"{first_header.alpha_w:.6g} x {first_header.alpha_h:.6g} radians"
    )

    # The custom sphere image overrides the default downloaded texture when provided.
    if args.sphere_image is None:
        sphere_image: str = os.path.join(script_dir, cfg.SPHERE_TEXTURE_FILE)

        # Download the celestial sphere texture if it doesn't already exist locally.
        if not os.path.exists(sphere_image):
            print(f"Downloading {cfg.SPHERE_TEXTURE_FILE}...")
            starmap_url = "https://raw.githubusercontent.com/Moone02/nrpy-visual-assets/96a39ba8510e401ea8ec836154fca5db3b13f4d3/noirlab2430b.tif"
            try:
                urllib.request.urlretrieve(starmap_url, sphere_image)
                print("Download complete.")
            except Exception as e:  # pylint: disable=broad-exception-caught
                print(
                    f"FATAL: Failed to download {cfg.SPHERE_TEXTURE_FILE} from {starmap_url}: {e}"
                )
                sys.exit(1)
    else:
        sphere_image = args.sphere_image
        print(f"Using custom celestial sphere texture: {sphere_image}")

    # Physical span encompasses full mathematical diameter of accretion disk geometry.
    source_physical_width = 2.0 * terminal_plane_max_radius

    # Preserve blueprint image aspect ratio while optionally resampling width.
    pixel_height = max(
        1, math.ceil(args.pixel_width * first_header.alpha_h / first_header.alpha_w)
    )

    terminal_plane_texture: Union[str, npt.NDArray[np.float64]]
    if args.terminal_plane_image is None:
        print("Creating terminal-plane disk array...")

        # Texture array represents the equatorial source disk using defined radii.
        disk_texture = rli.generate_source_disk_array(
            disk_physical_width=source_physical_width,
            disk_inner_radius=terminal_plane_min_radius,
            disk_outer_radius=terminal_plane_max_radius,
            colormap=cfg.COLORMAP,
        )

        # Cast the uint8 array to float64 to satisfy mypy type constraints.
        terminal_plane_texture = disk_texture.astype(np.float64)
    else:
        terminal_plane_texture = args.terminal_plane_image
        print(f"Using custom terminal-plane texture: {terminal_plane_texture}")

    print(f"Rendering image to: {args.output}...")
    if args.debug is not None:
        print("Debug failure-color key:")
        for term_type, label, color in rli.DEBUG_FAILURE_INFO:
            color_hex = f"#{color[0]:02X}{color[1]:02X}{color[2]:02X}"
            print(f"  {term_type}: {label} = {color_hex}")

    # Static image generator merges geodesic blueprint with texture maps for final image.
    rli.generate_static_lensed_image(
        output_filename=args.output,
        output_pixel_width=args.pixel_width,
        output_pixel_height=pixel_height,
        source_image_width=source_physical_width,
        sphere_image=sphere_image,
        source_image=terminal_plane_texture,
        blueprint_filenames=blueprint_files,
        custom_sphere_image=args.sphere_image is not None,
        custom_source_image=args.terminal_plane_image is not None,
        enable_debug=args.debug is not None,
        include_debug_key=args.debug == "key",
    )

    print("Visualization complete!")


if __name__ == "__main__":
    main()

# nrpy/examples/geodesic_visualizations/visualize_lensed_image.py
"""
Defines the visualization script for the Photon Geodesic Integrator.

This script reads binary blueprints to render the final lensed image. It maps physical
photon trajectories onto a static background texture and an accretion disk geometry,
using r_min and r_max to define source bounds.

Author: Dalton J. Moone
"""

import argparse
import math
import os
import sys
import urllib.request

import numpy as np

try:
    import blueprint_config_and_schema as cfg  # type: ignore
    import render_lensed_image as rli  # type: ignore
except ImportError:
    from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
    from nrpy.examples.geodesic_visualizations import render_lensed_image as rli


def main() -> None:
    """Parse parameters and orchestrate lensed image creation."""
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

    parser.add_argument(
        "--window_tiles_width",
        type=int,
        default=1,
        help="Number of horizontal tiles.",
    )

    parser.add_argument(
        "--window_tiles_height",
        type=int,
        default=1,
        help="Number of vertical tiles.",
    )

    # Physical parameters govern the source accretion disk geometry and window bounds.
    parser.add_argument(
        "--source_r_min",
        type=float,
        default=6.0,
        help="Inner physical radius r_min of the source disk.",
    )
    parser.add_argument(
        "--source_r_max",
        type=float,
        default=20.0,
        help="Outer physical radius r_max of the source disk.",
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
    parser.add_argument(
        "--pixel_width",
        type=int,
        default=750,
        help="Pixel width of the final static output image.",
    )

    # The parsed arguments struct contains all runtime configurations.
    args = parser.parse_args()

    # Absolute script directory ensures paths resolve independently of execution context.
    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Deterministically build the expected file paths based on tile dimensions.
    blueprint_zips = []
    for i in range(args.window_tiles_width):
        for j in range(args.window_tiles_height):
            zip_filename = f"light_blueprint_{i:02d}_{j:02d}.zip"
            zip_filepath = os.path.join(script_dir, zip_filename)

            # Existence check prevents reading missing or uncompiled blueprints.
            if not os.path.exists(zip_filepath):
                print(f"WARNING: Expected tile missing at '{zip_filepath}'.")
            else:
                blueprint_zips.append(zip_filepath)

    if not blueprint_zips:
        print("ERROR: No blueprint zip files found to process.")
        return

    print(f"Loading {len(blueprint_zips)} blueprint tiles...")

    # The texture path locates the high-resolution background celestial sphere image.
    starmap_path = os.path.join(script_dir, cfg.SPHERE_TEXTURE_FILE)

    # Download the celestial sphere texture if it doesn't already exist locally.
    if not os.path.exists(starmap_path):
        print(f"Downloading {cfg.SPHERE_TEXTURE_FILE}...")
        starmap_url = "https://raw.githubusercontent.com/Moone02/nrpy-visual-assets/96a39ba8510e401ea8ec836154fca5db3b13f4d3/noirlab2430b.tif"
        try:
            urllib.request.urlretrieve(starmap_url, starmap_path)
            print("Download complete.")
        except Exception as e:  # pylint: disable=broad-exception-caught
            print(
                f"FATAL: Failed to download {cfg.SPHERE_TEXTURE_FILE} from {starmap_url}: {e}"
            )
            sys.exit(1)

    # Physical span encompasses full mathematical diameter of accretion disk geometry.
    source_physical_width = 2.0 * args.source_r_max

    # Calculate the vertical pixel dimension to maintain the physical aspect ratio.
    pixel_height = math.ceil(
        (args.window_height / args.window_width) * args.pixel_width
    )

    print("Creating source disk array...")

    # Texture array represents the equatorial source disk using defined radii.
    disk_texture = rli.generate_source_disk_array(
        disk_physical_width=source_physical_width,
        disk_inner_radius=args.source_r_min,
        disk_outer_radius=args.source_r_max,
        colormap=cfg.COLORMAP,
    )

    # Cast the uint8 array to float64 to satisfy mypy type constraints.
    disk_texture_float = disk_texture.astype(np.float64)

    print(f"Rendering image to: {args.output}...")

    # Static image generator merges geodesic blueprint with texture maps for final image.
    rli.generate_static_lensed_image(
        output_filename=args.output,
        output_pixel_width=args.pixel_width,
        output_pixel_height=pixel_height,
        source_image_width=source_physical_width,
        sphere_image=starmap_path,
        source_image=disk_texture_float,
        blueprint_filenames=blueprint_zips,
        window_width=args.window_width,
        window_height=args.window_height,
    )

    print("Visualization complete!")


if __name__ == "__main__":
    main()

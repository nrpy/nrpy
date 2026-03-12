"""STATIC IMAGE RENDERING ENGINE."""

import os
from typing import Union

import numpy as np

from nrpy.helpers.geodesic_visualizations import config_and_types as cfg

# STATIC IMAGE RENDERING ENGINE
# This script translates raw geodesic data (blueprints) into visual images. It performs
# texture mapping for the celestial sphere (background stars) and the accretion disk.
# The core logic utilizes an accumulator-based approach to map ray endpoints from
# curved spacetime onto a 2D pixel grid via the camera's local window coordinates.


def _load_texture(image_input: Union[str, np.ndarray]) -> np.ndarray:
    """
    Load and normalize textures into float64 arrays.

    Converts 0-255 RGB integers to 0.0-1.0 floats for blending math.

    :param image_input: A file path (str) or a NumPy array to be loaded.
    :return: A float64 NumPy array normalized to the [0.0, 1.0] range.
    :raises FileNotFoundError: Raised if the provided file path does not exist.
    :raises TypeError: Raised if the input is not a string or NumPy array.
    """
    # pylint: disable=import-outside-toplevel
    from PIL import Image

    if isinstance(image_input, str):
        if not os.path.exists(image_input):
            raise FileNotFoundError(f"Texture file not found: {image_input}")
        with Image.open(image_input) as img:
            return np.array(img.convert("RGB")).astype(np.float64) / 255.0
    if isinstance(image_input, np.ndarray):
        texture_array = image_input.astype(np.float64)
        if np.max(texture_array) > 1.0:
            texture_array /= 255.0  # Normalize if input is 0-255
        return texture_array
    raise TypeError("Image input must be a file path (str) or a NumPy array.")


def generate_source_disk_array(
    pixel_width: int = 1024,
    disk_physical_width: float = 40.0,
    disk_inner_radius: float = 6.0,
    disk_outer_radius: float = 20.0,
    disk_temp_power_law: float = -0.75,
    colormap: str = "hot",
) -> np.ndarray:
    """
    Generate a synthetic top-down texture of an accretion disk.

    Calculates temperature based on radius $r$ and applies a colormap
    with smooth falloff at the inner and outer boundaries.

    :param pixel_width: The width/height of the generated square texture in pixels.
    :param disk_physical_width: The total physical width represented by the texture.
    :param disk_inner_radius: The inner radius of the disk where emission starts.
    :param disk_outer_radius: The outer radius of the disk where emission ends.
    :param disk_temp_power_law: The exponent for the temperature radial power law.
    :param colormap: The Matplotlib colormap name used for temperature mapping.
    :return: A uint8 RGB NumPy array of the generated accretion disk texture.
    """
    # pylint: disable=import-outside-toplevel
    import matplotlib.pyplot as plt

    half_width = disk_physical_width / 2.0
    # Create coordinate grid representing the flat source plane
    y_coords = np.linspace(-half_width, half_width, pixel_width)
    z_coords = np.linspace(half_width, -half_width, pixel_width)

    # 'xy' indexing creates the 2D spatial context for radius calculation
    yy, zz = np.meshgrid(y_coords, z_coords)
    radii = np.sqrt(yy**2 + zz**2)

    # Temperature model: power law decay from the center
    temperature = (radii / (disk_inner_radius + 1e-12)) ** disk_temp_power_law
    temperature[np.isnan(temperature)] = 0

    # Anti-aliasing / Smooth edges for the disk boundaries
    transition_width = 2.0 * (disk_physical_width / pixel_width)
    inner_falloff = np.clip(
        (radii - (disk_inner_radius - transition_width)) / transition_width, 0, 1
    )
    outer_falloff = 1.0 - np.clip((radii - disk_outer_radius) / transition_width, 0, 1)

    temperature *= inner_falloff * outer_falloff

    # Map normalized temperature to RGB color space
    max_temp = np.max(temperature)
    norm_temperature = temperature / max_temp if max_temp > 0 else temperature

    colormap_func = plt.colormaps[colormap]
    colors = colormap_func(norm_temperature)

    return (colors[:, :, :3] * 255).astype(np.uint8)


def generate_static_lensed_image(
    output_filename: str,
    output_pixel_width: int,
    source_image_width: float,
    sphere_image: Union[str, np.ndarray],
    source_image: Union[str, np.ndarray],
    blueprint_filename: str,
    window_width: float,
    chunk_size: int = cfg.CHUNK_SIZE,
    display_image: bool = True,
) -> None:
    """
    Generate a static lensed image from geodesic blueprint data.

    The function reads the binary blueprint, identifies where each ray terminated,
    and maps that endpoint to a pixel on the final image. Uses 'pixel_accumulator'
    to handle potential super-sampling if multiple rays hit the same pixel area.

    :param output_filename: Path where the final rendered image will be saved.
    :param output_pixel_width: Desired width of the output image in pixels.
    :param source_image_width: The physical width of the accretion disk source plane.
    :param sphere_image: The background celestial sphere texture (path or array).
    :param source_image: The accretion disk texture (path or array).
    :param blueprint_filename: Path to the binary blueprint file containing ray data.
    :param window_width: The physical width of the camera's local window.
    :param chunk_size: Number of records to read from the binary file at once.
    :param display_image: If True, opens the resulting image using the default viewer.
    :raises FileNotFoundError: Raised if the specified blueprint file is not found.
    """
    # pylint: disable=import-outside-toplevel
    from PIL import Image

    print(f"--- Generating Static Lensed Image: '{output_filename}' ---")

    if not os.path.exists(blueprint_filename):
        raise FileNotFoundError(f"Blueprint file not found: {blueprint_filename}")

    # Establish the FOV (Field of View) bounds for the camera window
    half_w = window_width / 2.0
    y_w_min, y_w_max = -half_w, half_w
    z_w_min, z_w_max = -half_w, half_w

    window_y_range = y_w_max - y_w_min
    window_z_range = z_w_max - z_w_min
    output_pixel_height = int(output_pixel_width * (window_z_range / window_y_range))

    # Load textures for mapping
    source_texture = _load_texture(source_image)
    sphere_texture = _load_texture(sphere_image)

    source_h, source_w, _ = source_texture.shape
    sphere_h, sphere_w, _ = sphere_texture.shape

    # Accumulators allow for averaging colors if multiple rays map to one pixel
    pixel_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width, 3), dtype=np.float64
    )
    count_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width), dtype=np.int32
    )

    # Determine how many rays are in the binary file
    file_size_bytes = os.path.getsize(blueprint_filename)
    record_size_bytes = cfg.BLUEPRINT_DTYPE.itemsize
    total_records = file_size_bytes // record_size_bytes

    with open(blueprint_filename, "rb") as f:
        records_processed = 0
        while records_processed < total_records:
            # Load a chunk of rays into memory to keep RAM usage stable
            chunk_data = np.fromfile(f, dtype=cfg.BLUEPRINT_DTYPE, count=chunk_size)
            if chunk_data.size == 0:
                break

            records_processed += chunk_data.size

            # Filter rays that actually hit within the defined camera window
            mask_in_view = (
                (chunk_data["y_w"] >= y_w_min)
                & (chunk_data["y_w"] < y_w_max)
                & (chunk_data["z_w"] >= z_w_min)
                & (chunk_data["z_w"] < z_w_max)
            )
            rays = chunk_data[mask_in_view]
            if rays.size == 0:
                continue

            # Convert local window coordinates $(y_w, z_w)$ to discrete pixel indices $(px, py)$
            px = (
                (rays["y_w"] - y_w_min) / window_y_range * (output_pixel_width - 1)
            ).astype(np.int32)
            py = (
                (z_w_max - rays["z_w"]) / window_z_range * (output_pixel_height - 1)
            ).astype(np.int32)

            # --- CASE 1: Source Plane Hits (Accretion Disk) ---
            is_source = rays["termination_type"] == cfg.TERM_SOURCE_PLANE
            if np.any(is_source):
                s_hits = rays[is_source]
                # Normalize $(y_s, z_s)$ coordinates to 0-1 for UV texture mapping
                norm_y = (
                    s_hits["y_s"] + (source_image_width / 2.0)
                ) / source_image_width
                norm_z = (
                    s_hits["z_s"] + (source_image_width / 2.0)
                ) / source_image_width

                # Map normalized coordinates to disk texture pixels
                px_s = np.clip(norm_y * (source_w - 1), 0, source_w - 1).astype(
                    np.int32
                )
                py_s = np.clip((1.0 - norm_z) * (source_h - 1), 0, source_h - 1).astype(
                    np.int32
                )

                # Add the sampled disk color to the accumulator at the camera pixel location
                np.add.at(
                    pixel_accumulator,
                    (py[is_source], px[is_source]),
                    source_texture[py_s, px_s],
                )

            # --- CASE 2: Celestial Sphere Hits (Background Stars) ---
            is_sphere = rays["termination_type"] == cfg.TERM_SPHERE
            if np.any(is_sphere):
                sph_hits = rays[is_sphere]
                # Map angular $\phi$ and $\theta$ to 0-1 range for spherical texture mapping
                norm_phi = (sph_hits["final_phi"] / (2 * np.pi)) % 1.0
                norm_theta = sph_hits["final_theta"] / np.pi

                # Map to celestial sphere texture pixels
                px_sph = np.clip(norm_phi * (sphere_w - 1), 0, sphere_w - 1).astype(
                    np.int32
                )
                py_sph = np.clip(norm_theta * (sphere_h - 1), 0, sphere_h - 1).astype(
                    np.int32
                )

                # Add the sampled star color to the accumulator
                np.add.at(
                    pixel_accumulator,
                    (py[is_sphere], px[is_sphere]),
                    sphere_texture[py_sph, px_sph],
                )

            # Track how many rays contributed to each pixel for final averaging
            np.add.at(count_accumulator, (py, px), 1)

    # Finalize the image: Divide accumulated color by hit count and clip to valid range
    hit_mask = count_accumulator > 0
    final_image_float = np.zeros_like(pixel_accumulator)
    final_image_float[hit_mask] = (
        pixel_accumulator[hit_mask] / count_accumulator[hit_mask, np.newaxis]
    )

    # Convert float64 to uint8 (0-255) and save
    img = Image.fromarray(
        (np.clip(final_image_float, 0, 1) * 255).astype(np.uint8), "RGB"
    )
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    img.save(output_filename)

    if display_image:
        img.show()


if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()
    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

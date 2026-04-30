# nrpy/examples/geodesic_visualizations/render_lensed_image.py
"""
Defines the static image rendering engine.

Translates raw geodesic data (blueprints) into visual images. Performs texture mapping
for the celestial sphere (background stars) and the accretion disk. Uses an
accumulator-based approach to map ray endpoints from curved spacetime onto a 2D
pixel grid via the camera's local window coordinates.

Author: Dalton J. Moone
"""

import os
import zipfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Callable, List, Optional, Tuple, TypeVar, Union

try:
    import numba as nb  # type: ignore # pylint: disable=import-error
except ImportError:

    _F = TypeVar("_F", bound=Callable[..., object])

    class _DummyNumba:
        """Dummy."""

        def njit(self, *_args: object, **_kwargs: object) -> Callable[[_F], _F]:
            """
            Return a decorator that returns the unmodified function.

            :param _args: Positional arguments for the njit decorator.
            :param _kwargs: Keyword arguments for the njit decorator.
            :return: A pass-through decorator.
            """

            def decorator(func: _F) -> _F:
                return func

            return decorator

    nb = _DummyNumba()
import numpy as np
import numpy.typing as npt

import nrpy.examples.geodesic_visualizations.blueprint_config_and_schema as cfg

# Global variables for zero-copy memory sharing across process workers
_WORKER_SOURCE_TEX: Optional[npt.NDArray[np.float64]] = None
_WORKER_SPHERE_TEX: Optional[npt.NDArray[np.float64]] = None


def _init_worker(
    source_tex: npt.NDArray[np.float64], sphere_tex: npt.NDArray[np.float64]
) -> None:
    # pylint: disable=global-statement
    # ProcessPoolExecutor grants workers instant, zero-copy read access to these global
    # arrays without any IPC serialization or pickling overhead.
    global _WORKER_SOURCE_TEX, _WORKER_SPHERE_TEX
    _WORKER_SOURCE_TEX = source_tex
    _WORKER_SPHERE_TEX = sphere_tex


def _load_texture(
    image_input: Union[str, npt.NDArray[np.float64]],
) -> npt.NDArray[np.float64]:
    """
    Load and normalize textures into float64 arrays.

    Converts 0-255 RGB integers to 0.0-1.0 floats for blending math.

    >>> import numpy as np
    >>> _load_texture(np.array([[[255.0, 0.0, 127.5]]])).tolist()
    [[[1.0, 0.0, 0.5]]]
    >>> _load_texture(np.array([[[1.0, 0.0, 0.5]]])).tolist()
    [[[1.0, 0.0, 0.5]]]
    >>> _load_texture(123)
    Traceback (most recent call last):
        ...
    TypeError: Image input must be a file path (str) or a NumPy array.

    :param image_input: A file path (str) or a NumPy array to be loaded.
    :return: A float64 NumPy array normalized to the [0.0, 1.0] range.
    :raises FileNotFoundError: Raised if the provided file path does not exist.
    :raises TypeError: Raised if the input is not a string or NumPy array.
    """
    if isinstance(image_input, str):
        if not os.path.exists(image_input):
            raise FileNotFoundError(f"Texture file not found: {image_input}")
        # pylint: disable=import-outside-toplevel, import-error
        from PIL import Image  # type: ignore

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
) -> npt.NDArray[np.uint8]:
    """
    Generate a synthetic top-down texture of an accretion disk.

    Calculates temperature based on radius r and applies a colormap
    with smooth falloff at the inner and outer boundaries.

    :param pixel_width: The width/height of the generated square texture in pixels.
    :param disk_physical_width: The total physical width represented by the texture.
    :param disk_inner_radius: The inner radius of the disk where emission starts.
    :param disk_outer_radius: The outer radius of the disk where emission ends.
    :param disk_temp_power_law: The exponent for the temperature radial power law.
    :param colormap: The Matplotlib colormap name used for temperature mapping.
    :return: A uint8 RGB NumPy array of the generated accretion disk texture.
    """
    # pylint: disable=import-outside-toplevel, import-error
    import matplotlib.pyplot as plt  # type: ignore

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

    # Explicitly type the result
    final_colors: npt.NDArray[np.uint8] = (colors[:, :, :3] * 255).astype(np.uint8)
    return final_colors


@nb.njit(nogil=True, fastmath=True)  # type: ignore
def _accumulate_ray_hits_jit(
    rays_y_w: npt.NDArray[np.float64],
    rays_z_w: npt.NDArray[np.float64],
    rays_term_type: npt.NDArray[np.int32],
    rays_y_s: npt.NDArray[np.float64],
    rays_z_s: npt.NDArray[np.float64],
    rays_phi: npt.NDArray[np.float64],
    rays_theta: npt.NDArray[np.float64],
    y_w_min: float,
    y_w_max: float,
    z_w_min: float,
    z_w_max: float,
    window_y_range: float,
    window_z_range: float,
    output_pixel_width: int,
    output_pixel_height: int,
    source_image_width: float,
    source_tex: npt.NDArray[np.float64],
    sphere_tex: npt.NDArray[np.float64],
    local_pixel_acc: npt.NDArray[np.float64],
    local_count_acc: npt.NDArray[np.int32],
    term_source: int,
    term_sphere: int,
) -> None:
    # JIT-compiled core math loop for mapping rays to pixels.
    # Bypasses the GIL and executes in optimized machine code.
    # Maintains float64 precision internally to avoid catastrophic cancellation.
    source_h = source_tex.shape[0]
    source_w = source_tex.shape[1]
    sphere_h = sphere_tex.shape[0]
    sphere_w = sphere_tex.shape[1]

    for i, y_w in enumerate(rays_y_w):
        z_w = rays_z_w[i]

        # Filter out rays that don't hit the camera window bounds
        if not (y_w_min <= y_w < y_w_max and z_w_min <= z_w < z_w_max):
            continue

        # Convert local window coordinates to discrete pixel indices
        px = int((y_w - y_w_min) / window_y_range * (output_pixel_width - 1))
        py = int((z_w_max - z_w) / window_z_range * (output_pixel_height - 1))

        term = rays_term_type[i]

        # --- CASE 1: Source Plane Hits (Accretion Disk) ---
        if term == term_source:
            y_s = rays_y_s[i]
            z_s = rays_z_s[i]
            norm_y = (y_s + (source_image_width / 2.0)) / source_image_width
            norm_z = (z_s + (source_image_width / 2.0)) / source_image_width

            px_s = int(norm_y * (source_w - 1))
            py_s = int((1.0 - norm_z) * (source_h - 1))

            px_s = max(0, min(px_s, source_w - 1))
            py_s = max(0, min(py_s, source_h - 1))

            local_pixel_acc[py, px, 0] += source_tex[py_s, px_s, 0]
            local_pixel_acc[py, px, 1] += source_tex[py_s, px_s, 1]
            local_pixel_acc[py, px, 2] += source_tex[py_s, px_s, 2]
            local_count_acc[py, px] += 1

        # --- CASE 2: Celestial Sphere Hits (Background Stars) ---
        elif term == term_sphere:
            phi = rays_phi[i]
            theta = rays_theta[i]

            norm_phi = (phi / (2 * np.pi)) % 1.0
            norm_theta = theta / np.pi

            px_sph = int(norm_phi * (sphere_w - 1))
            py_sph = int(norm_theta * (sphere_h - 1))

            px_sph = max(0, min(px_sph, sphere_w - 1))
            py_sph = max(0, min(py_sph, sphere_h - 1))

            local_pixel_acc[py, px, 0] += sphere_tex[py_sph, px_sph, 0]
            local_pixel_acc[py, px, 1] += sphere_tex[py_sph, px_sph, 1]
            local_pixel_acc[py, px, 2] += sphere_tex[py_sph, px_sph, 2]
            local_count_acc[py, px] += 1


def _process_blueprint_tile(
    zip_filename: str,
    chunk_bytes: int,
    y_w_min: float,
    y_w_max: float,
    z_w_min: float,
    z_w_max: float,
    window_y_range: float,
    window_z_range: float,
    output_pixel_width: int,
    output_pixel_height: int,
    source_image_width: float,
) -> Tuple[
    npt.NDArray[np.int64],
    npt.NDArray[np.int64],
    npt.NDArray[np.float32],
    npt.NDArray[np.int32],
]:
    # Worker function to process a blueprint tile archive and return sparse arrays.
    # Reads binary chunks, passes arrays to the JIT-compiled math function,
    # and compresses the result to avoid IPC memory overhead.
    if _WORKER_SOURCE_TEX is None or _WORKER_SPHERE_TEX is None:
        raise RuntimeError(
            "Worker globals not initialized. Ensure initializer is set in ProcessPoolExecutor."
        )

    # High-precision float64 arrays strictly for internal math
    local_pixel_acc = np.zeros(
        (output_pixel_height, output_pixel_width, 3), dtype=np.float64
    )
    local_count_acc = np.zeros(
        (output_pixel_height, output_pixel_width), dtype=np.int32
    )

    with zipfile.ZipFile(zip_filename, "r") as zf:
        bin_files = [name for name in zf.namelist() if name.endswith(".bin")]
        if not bin_files:
            return np.array([]), np.array([]), np.array([]), np.array([])

        with zf.open(bin_files[0], "r") as f:
            while True:
                raw_bytes = f.read(chunk_bytes)
                if not raw_bytes:
                    break

                chunk_data = np.frombuffer(raw_bytes, dtype=cfg.BLUEPRINT_DTYPE)

                # Pass 1D slices of the structured array to the JIT compiler
                _accumulate_ray_hits_jit(
                    chunk_data["y_w"],
                    chunk_data["z_w"],
                    chunk_data["termination_type"],
                    chunk_data["y_s"],
                    chunk_data["z_s"],
                    chunk_data["final_phi"],
                    chunk_data["final_theta"],
                    y_w_min,
                    y_w_max,
                    z_w_min,
                    z_w_max,
                    window_y_range,
                    window_z_range,
                    output_pixel_width,
                    output_pixel_height,
                    source_image_width,
                    _WORKER_SOURCE_TEX,
                    _WORKER_SPHERE_TEX,
                    local_pixel_acc,
                    local_count_acc,
                    cfg.TERM_SOURCE_PLANE,
                    cfg.TERM_SPHERE,
                )

    # --- SPARSE MAP-REDUCE COMPRESSION ---
    # Downcast the output payload to float32 only AFTER the float64 math is done
    hit_mask = local_count_acc > 0
    flat_y, flat_x = np.nonzero(hit_mask)
    flat_colors = local_pixel_acc[hit_mask].astype(np.float32)
    flat_counts = local_count_acc[hit_mask]

    return flat_y, flat_x, flat_colors, flat_counts


def _save_image_to_file(img: "Image.Image", output_filename: str) -> None:  # type: ignore[name-defined]
    """
    Save a PIL Image to a file, creating the parent directory if necessary.

    >>> class DummyImage:
    ...     def save(self, path: str) -> None:
    ...         with open(path, "w", encoding="utf-8") as f:
    ...             pass
    >>> dummy_img = DummyImage()
    >>> _save_image_to_file(dummy_img, "test_out_bare.png")
    >>> os.path.exists("test_out_bare.png")
    True
    >>> _save_image_to_file(dummy_img, "test_out_dir/test_nested.png")
    >>> os.path.exists("test_out_dir/test_nested.png")
    True
    >>> os.remove("test_out_bare.png")
    >>> os.remove("test_out_dir/test_nested.png")
    >>> os.rmdir("test_out_dir")

    :param img: The PIL Image object to save.
    :param output_filename: The target file path.
    """
    output_dir = os.path.dirname(output_filename)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    img.save(output_filename)


def generate_static_lensed_image(
    output_filename: str,
    output_pixel_width: int,
    output_pixel_height: int,
    source_image_width: float,
    sphere_image: Union[str, npt.NDArray[np.float64]],
    source_image: Union[str, npt.NDArray[np.float64]],
    blueprint_filenames: List[str],
    window_width: float,
    window_height: float,
    chunk_size: int = cfg.CHUNK_SIZE,
    display_image: bool = True,
) -> None:
    """
    Generate a static lensed image from geodesic blueprint data.

    Reads binary blueprints streaming from compressed archives, identifies where each
    ray terminated, and maps that endpoint to a pixel on the final image. Uses a
    process pool executor to process individual blueprint tiles in parallel to improve
    loading speeds.

    :param output_filename: Path where the final rendered image will be saved.
    :param output_pixel_width: Desired width of the output image in pixels.
    :param output_pixel_height: Desired height of the output image in pixels.
    :param source_image_width: The physical width of the accretion disk source plane.
    :param sphere_image: The background celestial sphere texture (path or array).
    :param source_image: The accretion disk texture (path or array).
    :param blueprint_filenames: List of paths to the zipped blueprint files containing ray data.
    :param window_height: The physical height of the camera's local window.
    :param window_width: The physical width of the camera's local window.
    :param chunk_size: Number of records to read from the binary stream at once.
    :param display_image: If True, opens the resulting image using the default viewer.
    """
    # pylint: disable=import-outside-toplevel, import-error
    from PIL import Image

    print(f"--- Generating Static Lensed Image: '{output_filename}' ---")

    # Establish the FOV (Field of View) bounds for the camera window
    half_w = window_width / 2.0
    half_h = window_height / 2.0
    y_w_min, y_w_max = -half_w, half_w
    z_w_min, z_w_max = -half_h, half_h

    window_y_range = y_w_max - y_w_min
    window_z_range = z_w_max - z_w_min

    # Load textures for mapping
    source_texture = _load_texture(source_image)
    sphere_texture = _load_texture(sphere_image)

    # Accumulators allow for averaging colors if multiple rays map to one pixel
    pixel_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width, 3), dtype=np.float64
    )
    count_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width), dtype=np.int32
    )

    record_size_bytes = cfg.BLUEPRINT_DTYPE.itemsize
    chunk_bytes = chunk_size * record_size_bytes

    # I/O Throttling
    # Cap workers to prevent disk thrashing.
    max_io_workers = min(8, (os.cpu_count() or 4))

    print(
        f"  -> Dispatching {len(blueprint_filenames)} tiles for parallel processing ({max_io_workers} workers)..."
    )

    with ProcessPoolExecutor(
        max_workers=max_io_workers,
        initializer=_init_worker,
        initargs=(source_texture, sphere_texture),
    ) as executor:
        futures = []
        for zip_filename in blueprint_filenames:
            futures.append(
                executor.submit(
                    _process_blueprint_tile,
                    zip_filename,
                    chunk_bytes,
                    y_w_min,
                    y_w_max,
                    z_w_min,
                    z_w_max,
                    window_y_range,
                    window_z_range,
                    output_pixel_width,
                    output_pixel_height,
                    source_image_width,
                )
            )

        # Merge the sparse 1D local accumulators back into the main global accumulator
        for i, future in enumerate(as_completed(futures)):
            flat_y, flat_x, flat_colors, flat_counts = future.result()

            if len(flat_counts) > 0:
                np.add.at(pixel_accumulator, (flat_y, flat_x), flat_colors)
                np.add.at(count_accumulator, (flat_y, flat_x), flat_counts)

            print(f"  -> Completed tile {i + 1}/{len(blueprint_filenames)}.")

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
    _save_image_to_file(img, output_filename)

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

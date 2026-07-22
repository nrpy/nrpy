# nrpy/examples/geodesic_visualizations/render_lensed_image.py
"""
Defines the static image rendering engine.

Translates raw geodesic data (blueprints) into visual images. Performs texture mapping
for the celestial sphere (background stars) and the accretion disk. Uses an
accumulator-based approach to map ray endpoints from curved spacetime onto a 2D
pixel grid via the camera's local window coordinates.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import os
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

try:
    import blueprint_config_and_schema as cfg  # type: ignore
    import blueprint_io  # type: ignore
except ImportError:
    from nrpy.examples.geodesic_visualizations import blueprint_config_and_schema as cfg
    from nrpy.examples.geodesic_visualizations import blueprint_io

# Global variables for zero-copy memory sharing across process workers
_WORKER_SOURCE_TEX: Optional[npt.NDArray[np.float64]] = None
_WORKER_SPHERE_TEX: Optional[npt.NDArray[np.float64]] = None

# Fixed high-contrast colors for opt-in termination diagnostics.
DEBUG_FAILURE_INFO: Tuple[Tuple[int, str, Tuple[int, int, int]], ...] = (
    (
        cfg.TERM_EVOLUTION_MEASURE_EXCEEDED,
        "TERM_EVOLUTION_MEASURE_EXCEEDED",
        (255, 0, 0),
    ),
    (cfg.TERM_RKF45_REJECTION_LIMIT, "TERM_RKF45_REJECTION_LIMIT", (255, 0, 168)),
    (cfg.TERM_T_MAX_EXCEEDED, "TERM_T_MAX_EXCEEDED", (255, 122, 0)),
    (cfg.TERM_SLOT_MANAGER_ERROR, "TERM_SLOT_MANAGER_ERROR", (182, 255, 0)),
    (cfg.TERM_FAILURE, "TERM_FAILURE", (0, 229, 255)),
)


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
    warn_if_not_png: bool = False,
    center_crop_to_square: bool = False,
    warn_if_non_equirectangular: bool = False,
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
    :param warn_if_not_png: Whether to print a warning if the provided path is not a PNG.
    :param center_crop_to_square: Whether to center-crop a loaded file image to a square.
    :param warn_if_non_equirectangular: Whether to warn when a loaded file image is not
        approximately 2:1 in aspect ratio.
    :return: A float64 NumPy array normalized to the [0.0, 1.0] range.
    :raises FileNotFoundError: Raised if the provided file path does not exist.
    :raises TypeError: Raised if the input is not a string or NumPy array.
    """
    if isinstance(image_input, str):
        if not os.path.exists(image_input):
            raise FileNotFoundError(f"Texture file not found: {image_input}")
        if warn_if_not_png and not image_input.lower().endswith(".png"):
            print(
                f"WARNING: Custom texture path '{image_input}' does not end in '.png'. Attempting to load it anyway."
            )
        # pylint: disable=import-outside-toplevel, import-error
        from PIL import Image

        with Image.open(image_input) as img:
            texture_array = np.array(img.convert("RGB")).astype(np.float64)

        if center_crop_to_square:
            texture_array = _center_crop_texture_to_square(texture_array, image_input)
        if warn_if_non_equirectangular:
            _warn_if_texture_not_equirectangular(texture_array, image_input)
        return texture_array / 255.0
    if isinstance(image_input, np.ndarray):
        texture_array = image_input.astype(np.float64)
        if np.max(texture_array) > 1.0:
            texture_array /= 255.0  # Normalize if input is 0-255
        return texture_array
    raise TypeError("Image input must be a file path (str) or a NumPy array.")


def _center_crop_texture_to_square(
    texture_array: npt.NDArray[np.float64], texture_path: str
) -> npt.NDArray[np.float64]:
    """
    Center-crop a texture to a square using the smaller image dimension.

    :param texture_array: The RGB image data to crop.
    :param texture_path: Path used only for warning context.
    :return: The original array if already square, otherwise a centered square crop.
    """
    texture_height = texture_array.shape[0]
    texture_width = texture_array.shape[1]
    if texture_height == texture_width:
        return texture_array

    print(
        f"WARNING: Custom source image '{texture_path}' is {texture_width}x{texture_height}; center-cropping to a square."
    )
    cropped_size = min(texture_height, texture_width)

    if texture_width > texture_height:
        left = (texture_width - cropped_size) // 2
        right = left + cropped_size
        return texture_array[:, left:right, :]

    top = (texture_height - cropped_size) // 2
    bottom = top + cropped_size
    return texture_array[top:bottom, :, :]


def _warn_if_texture_not_equirectangular(
    texture_array: npt.NDArray[np.float64], texture_path: str
) -> None:
    """
    Warn if a texture does not appear approximately equirectangular.

    :param texture_array: The RGB image data to inspect.
    :param texture_path: Path used only for warning context.
    """
    texture_height = texture_array.shape[0]
    texture_width = texture_array.shape[1]
    aspect_ratio = texture_width / texture_height
    if not np.isclose(aspect_ratio, 2.0, rtol=0.05, atol=0.05):
        print(
            f"WARNING: Custom celestial sphere image '{texture_path}' is {texture_width}x{texture_height}; expected an approximately 2:1 texture."
        )


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

    # Explicitly type the result
    final_colors: npt.NDArray[np.uint8] = (colors[:, :, :3] * 255).astype(np.uint8)
    return final_colors


@nb.njit(inline="always", nogil=True, fastmath=True)  # type: ignore
def pixel_indices_from_ordinal(
    ordinal: int,
    tile_x: int,
    tile_y: int,
    scan_density: int,
    tiles_width: int,
    tiles_height: int,
    output_pixel_width: int,
    output_pixel_height: int,
) -> Tuple[int, int]:
    """
    Map one dense launch-ray ordinal to a center-preserving output pixel.

    :param ordinal: Tile-local row-major ray ordinal.
    :param tile_x: Horizontal tile index.
    :param tile_y: Vertical tile index.
    :param scan_density: Rays per tile side.
    :param tiles_width: Number of tiles across the complete window.
    :param tiles_height: Number of tiles down the complete window.
    :param output_pixel_width: Width of the output image.
    :param output_pixel_height: Height of the output image.
    :return: ``(pixel_x, pixel_y)`` in output-image coordinates.
    """
    local_row = ordinal // scan_density
    local_col = ordinal % scan_density
    global_x = tile_x * scan_density + local_col
    global_y = tile_y * scan_density + local_row
    total_width = tiles_width * scan_density
    total_height = tiles_height * scan_density
    pixel_x = ((2 * global_x + 1) * output_pixel_width) // (2 * total_width)
    pixel_y = (
        output_pixel_height
        - 1
        - (((2 * global_y + 1) * output_pixel_height) // (2 * total_height))
    )
    return pixel_x, pixel_y


@nb.njit(nogil=True, fastmath=True)  # type: ignore
def _accumulate_ray_hits_jit(
    rays_term_type: npt.NDArray[np.int32],
    rays_y_s: npt.NDArray[np.float64],
    rays_z_s: npt.NDArray[np.float64],
    rays_phi: npt.NDArray[np.float64],
    rays_theta: npt.NDArray[np.float64],
    record_start: int,
    tile_x: int,
    tile_y: int,
    scan_density: int,
    tiles_width: int,
    tiles_height: int,
    output_pixel_width: int,
    output_pixel_height: int,
    source_image_width: float,
    source_tex: npt.NDArray[np.float64],
    sphere_tex: npt.NDArray[np.float64],
    local_pixel_acc: npt.NDArray[np.float64],
    local_count_acc: npt.NDArray[np.int32],
    local_debug_term_types: npt.NDArray[np.int8],
    term_source: int,
    term_sphere: int,
    term_evolution_measure: int,
    term_fail_generic: int,
    enable_debug: bool,
) -> None:
    # JIT-compiled core math loop for mapping rays to pixels.
    # Bypasses the GIL and executes in optimized machine code.
    # Maintains float64 precision internally to avoid catastrophic cancellation.
    source_h = source_tex.shape[0]
    source_w = source_tex.shape[1]
    sphere_h = sphere_tex.shape[0]
    sphere_w = sphere_tex.shape[1]

    for i, term in enumerate(rays_term_type):
        px, py = pixel_indices_from_ordinal(
            record_start + i,
            tile_x,
            tile_y,
            scan_density,
            tiles_width,
            tiles_height,
            output_pixel_width,
            output_pixel_height,
        )

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

        # Preserve requested integration failures for the opt-in debug overlay.
        if (
            enable_debug
            and term_evolution_measure <= term <= term_fail_generic
            and term < local_debug_term_types[py, px]
        ):
            local_debug_term_types[py, px] = term


def _process_blueprint_tile(
    blueprint_filename: str,
    chunk_records: int,
    output_pixel_width: int,
    output_pixel_height: int,
    source_image_width: float,
    enable_debug: bool,
) -> Tuple[
    npt.NDArray[np.int64],
    npt.NDArray[np.int64],
    npt.NDArray[np.float32],
    npt.NDArray[np.int32],
    npt.NDArray[np.int64],
    npt.NDArray[np.int64],
    npt.NDArray[np.int8],
]:
    # Worker function to process one native blueprint tile and return sparse arrays.
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
    if enable_debug:
        local_debug_term_types = np.full(
            (output_pixel_height, output_pixel_width),
            cfg.TERM_ACTIVE,
            dtype=np.int8,
        )
    else:
        local_debug_term_types = np.empty((0, 0), dtype=np.int8)

    for header, record_start, chunk_data in blueprint_io.iter_blueprint_chunks(
        blueprint_filename, chunk_records
    ):
        _accumulate_ray_hits_jit(
            chunk_data["termination_type"],
            chunk_data["y_s"],
            chunk_data["z_s"],
            chunk_data["final_phi"],
            chunk_data["final_theta"],
            record_start,
            header.tile_x,
            header.tile_y,
            header.scan_density,
            header.tiles_width,
            header.tiles_height,
            output_pixel_width,
            output_pixel_height,
            source_image_width,
            _WORKER_SOURCE_TEX,
            _WORKER_SPHERE_TEX,
            local_pixel_acc,
            local_count_acc,
            local_debug_term_types,
            cfg.TERM_SOURCE_PLANE,
            cfg.TERM_COORD_RADIUS_EXCEEDED,
            cfg.TERM_EVOLUTION_MEASURE_EXCEEDED,
            cfg.TERM_FAILURE,
            enable_debug,
        )

    # --- SPARSE MAP-REDUCE COMPRESSION ---
    # Downcast the output payload to float32 only AFTER the float64 math is done
    hit_mask = local_count_acc > 0
    flat_y, flat_x = np.nonzero(hit_mask)
    flat_colors = local_pixel_acc[hit_mask].astype(np.float32)
    flat_counts = local_count_acc[hit_mask]

    if enable_debug:
        debug_mask = local_debug_term_types < cfg.TERM_ACTIVE
        debug_flat_y, debug_flat_x = np.nonzero(debug_mask)
        debug_flat_term_types = local_debug_term_types[debug_mask]
    else:
        debug_flat_y = np.array([], dtype=np.int64)
        debug_flat_x = np.array([], dtype=np.int64)
        debug_flat_term_types = np.array([], dtype=np.int8)

    return (
        flat_y,
        flat_x,
        flat_colors,
        flat_counts,
        debug_flat_y,
        debug_flat_x,
        debug_flat_term_types,
    )


def _add_debug_legend(img: "Image.Image") -> None:  # type: ignore[name-defined]
    """
    Add the fixed failure-color key in the image's bottom-right corner.

    :param img: Image receiving the failure-color legend.
    """
    # pylint: disable=import-outside-toplevel, import-error
    from PIL import ImageDraw, ImageFont

    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    padding = 4
    swatch_size = 8
    line_height = 11
    labels = [failure_info[1] for failure_info in DEBUG_FAILURE_INFO]
    text_width = max(int(draw.textlength(label, font=font)) for label in labels)
    legend_width = 3 * padding + swatch_size + text_width
    legend_height = 2 * padding + line_height * len(DEBUG_FAILURE_INFO)
    left = max(0, img.width - legend_width - padding)
    top = max(0, img.height - legend_height - padding)

    # Draw last so the key is visually distinct from ray-traced pixels.
    draw.rectangle((left, top, img.width - 1, img.height - 1), fill=(0, 0, 0))
    for index, (_, label, color) in enumerate(DEBUG_FAILURE_INFO):
        y_coord = top + padding + index * line_height
        draw.rectangle(
            (
                left + padding,
                y_coord,
                left + padding + swatch_size,
                y_coord + swatch_size,
            ),
            fill=color,
        )
        draw.text(
            (left + 2 * padding + swatch_size, y_coord - 2),
            label,
            fill=(255, 255, 255),
            font=font,
        )


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
    chunk_size: int = cfg.CHUNK_SIZE,
    display_image: bool = True,
    custom_sphere_image: bool = False,
    custom_source_image: bool = False,
    enable_debug: bool = False,
    include_debug_key: bool = False,
) -> None:
    """
    Generate a static lensed image from geodesic blueprint data.

    Reads native blueprint records in chunks, identifies where each ray terminated,
    and maps each dense launch ordinal to a pixel on the final image. Uses a process
    pool executor to process individual blueprint tiles in parallel.

    :param output_filename: Path where the final rendered image will be saved.
    :param output_pixel_width: Desired width of the output image in pixels.
    :param output_pixel_height: Desired height of the output image in pixels.
    :param source_image_width: The physical width of the accretion disk source plane.
    :param sphere_image: The background celestial sphere texture (path or array).
    :param source_image: The accretion disk texture (path or array).
    :param blueprint_filenames: List of paths to native blueprint files.
    :param chunk_size: Number of records to read from the binary stream at once.
    :param display_image: If True, opens the resulting image using the default viewer.
    :param custom_sphere_image: Whether the sphere texture came from a custom CLI override.
    :param custom_source_image: Whether the source texture came from a custom CLI override.
    :param enable_debug: Whether to overlay requested integration failures with fixed colors.
    :param include_debug_key: Whether to add the debug failure-color key to the image.
    :raises ValueError: If the blueprint set is incomplete or has inconsistent metadata.
    """
    # pylint: disable=import-outside-toplevel, import-error
    from PIL import Image

    print(f"--- Generating Static Lensed Image: '{output_filename}' ---")

    # Load textures for mapping
    source_texture = _load_texture(
        source_image,
        warn_if_not_png=custom_source_image,
        center_crop_to_square=custom_source_image,
    )
    sphere_texture = _load_texture(
        sphere_image,
        warn_if_not_png=custom_sphere_image,
        warn_if_non_equirectangular=custom_sphere_image,
    )

    # Accumulators allow for averaging colors if multiple rays map to one pixel
    pixel_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width, 3), dtype=np.float64
    )
    count_accumulator = np.zeros(
        (output_pixel_height, output_pixel_width), dtype=np.int32
    )
    if enable_debug:
        debug_term_types = np.full(
            (output_pixel_height, output_pixel_width), cfg.TERM_ACTIVE, dtype=np.int8
        )
    else:
        debug_term_types = np.empty((0, 0), dtype=np.int8)

    if not blueprint_filenames:
        raise ValueError("At least one blueprint artifact is required")
    headers = [
        blueprint_io.read_blueprint_header(filename) for filename in blueprint_filenames
    ]
    tiles_width = headers[0].tiles_width
    tiles_height = headers[0].tiles_height
    expected_tiles = {
        (tile_x, tile_y)
        for tile_x in range(tiles_width)
        for tile_y in range(tiles_height)
    }
    actual_tiles = {(header.tile_x, header.tile_y) for header in headers}
    if len(headers) != len(expected_tiles) or actual_tiles != expected_tiles:
        raise ValueError("Blueprint artifacts do not form the complete tile set")
    if any(
        header.tiles_width != tiles_width or header.tiles_height != tiles_height
        for header in headers
    ):
        raise ValueError("Blueprint headers disagree about the tile grid")

    chunk_records = chunk_size

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
        for blueprint_filename in blueprint_filenames:
            futures.append(
                executor.submit(
                    _process_blueprint_tile,
                    blueprint_filename,
                    chunk_records,
                    output_pixel_width,
                    output_pixel_height,
                    source_image_width,
                    enable_debug,
                )
            )

        # Merge the sparse 1D local accumulators back into the main global accumulator
        for i, future in enumerate(as_completed(futures)):
            (
                flat_y,
                flat_x,
                flat_colors,
                flat_counts,
                debug_flat_y,
                debug_flat_x,
                debug_flat_term_types,
            ) = future.result()

            if len(flat_counts) > 0:
                np.add.at(pixel_accumulator, (flat_y, flat_x), flat_colors)
                np.add.at(count_accumulator, (flat_y, flat_x), flat_counts)
            if enable_debug and len(debug_flat_term_types) > 0:
                np.minimum.at(
                    debug_term_types,
                    (debug_flat_y, debug_flat_x),
                    debug_flat_term_types,
                )

            print(f"  -> Completed tile {i + 1}/{len(blueprint_filenames)}.")

    # Finalize the image: Divide accumulated color by hit count and clip to valid range
    hit_mask = count_accumulator > 0
    final_image_float = np.zeros_like(pixel_accumulator)
    final_image_float[hit_mask] = (
        pixel_accumulator[hit_mask] / count_accumulator[hit_mask, np.newaxis]
    )
    if enable_debug:
        for term_type, _, color in DEBUG_FAILURE_INFO:
            final_image_float[debug_term_types == term_type] = np.array(color) / 255.0

    # Convert float64 to uint8 (0-255) and save
    img = Image.fromarray(
        (np.clip(final_image_float, 0, 1) * 255).astype(np.uint8), "RGB"
    )
    if enable_debug and include_debug_key:
        _add_debug_legend(img)
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

import numpy as np
import os
from PIL import Image
from typing import Union
import config_and_types as cfg

def _load_texture(image_input: Union[str, np.ndarray]) -> np.ndarray:
    if isinstance(image_input, str):
        if not os.path.exists(image_input):
            raise FileNotFoundError(f"Texture file not found: {image_input}")
        with Image.open(image_input) as img:
            return np.array(img.convert("RGB")).astype(np.float64) / 255.0
    elif isinstance(image_input, np.ndarray):
        texture_array = image_input.astype(np.float64)
        if np.max(texture_array) > 1.0:
            texture_array /= 255.0
        return texture_array
    else:
        raise TypeError("Image input must be a file path (str) or a NumPy array.")

def generate_source_disk_array(
    pixel_width: int = 1024,
    disk_physical_width: float = 40.0,
    disk_inner_radius: float = 6.0,
    disk_outer_radius: float = 20.0,
    disk_temp_power_law: float = -0.75,
    colormap: str = 'hot'
) -> np.ndarray:
    import matplotlib.pyplot as plt
    
    half_width = disk_physical_width / 2.0
    y_coords = np.linspace(-half_width, half_width, pixel_width)
    # Mapping z_coords such that the top of the array is the max value
    z_coords = np.linspace(half_width, -half_width, pixel_width)
    
    # 'xy' indexing: first dimension is y (columns), second is z (rows)
    yy, zz = np.meshgrid(y_coords, z_coords)
    radii = np.sqrt(yy**2 + zz**2)

    temperature = (radii / (disk_inner_radius + 1e-12))**disk_temp_power_law
    temperature[np.isnan(temperature)] = 0
    
    transition_width = 2.0 * (disk_physical_width / pixel_width)
    inner_falloff = np.clip((radii - (disk_inner_radius - transition_width)) / transition_width, 0, 1)
    outer_falloff = 1.0 - np.clip((radii - disk_outer_radius) / transition_width, 0, 1)
    
    temperature *= (inner_falloff * outer_falloff)
    
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
    display_image: bool = True 
) -> None:
    print(f"--- Generating Static Lensed Image: '{output_filename}' ---")
    
    if not os.path.exists(blueprint_filename):
        raise FileNotFoundError(f"Blueprint file not found: {blueprint_filename}")

    half_w = window_width / 2.0
    y_w_min, y_w_max = -half_w, half_w
    z_w_min, z_w_max = -half_w, half_w
    
    window_y_range = y_w_max - y_w_min
    window_z_range = z_w_max - z_w_min
    output_pixel_height = int(output_pixel_width * (window_z_range / window_y_range))
    
    source_texture = _load_texture(source_image)
    sphere_texture = _load_texture(sphere_image)
    
    source_h, source_w, _ = source_texture.shape
    sphere_h, sphere_w, _ = sphere_texture.shape

    pixel_accumulator = np.zeros((output_pixel_height, output_pixel_width, 3), dtype=np.float64)
    count_accumulator = np.zeros((output_pixel_height, output_pixel_width), dtype=np.int32)

    file_size_bytes = os.path.getsize(blueprint_filename)
    record_size_bytes = cfg.BLUEPRINT_DTYPE.itemsize
    total_records = file_size_bytes // record_size_bytes

    with open(blueprint_filename, 'rb') as f:
        records_processed = 0
        while records_processed < total_records:
            chunk_data = np.fromfile(f, dtype=cfg.BLUEPRINT_DTYPE, count=chunk_size)
            if chunk_data.size == 0: break
                
            records_processed += chunk_data.size

            mask_in_view = (
                (chunk_data['y_w'] >= y_w_min) & (chunk_data['y_w'] < y_w_max) &
                (chunk_data['z_w'] >= z_w_min) & (chunk_data['z_w'] < z_w_max)
            )
            rays = chunk_data[mask_in_view]
            if rays.size == 0: continue

            # y_w_min -> px=0 (Left). z_w_max -> py=0 (Top, so z_min is Bottom).
            px = ((rays['y_w'] - y_w_min) / window_y_range * (output_pixel_width - 1)).astype(np.int32)
            py = ((z_w_max - rays['z_w']) / window_z_range * (output_pixel_height - 1)).astype(np.int32)

            # 1. Source Plane Hits
            is_source = (rays['termination_type'] == cfg.TERM_SOURCE_PLANE)
            if np.any(is_source):
                s_hits = rays[is_source]
                # y_s_min -> left. z_s_min -> bottom.
                norm_y = (s_hits['y_s'] + (source_image_width / 2.0)) / source_image_width
                norm_z = (s_hits['z_s'] + (source_image_width / 2.0)) / source_image_width
                
                px_s = np.clip(norm_y * (source_w - 1), 0, source_w - 1).astype(np.int32)
                py_s = np.clip((1.0 - norm_z) * (source_h - 1), 0, source_h - 1).astype(np.int32)
                
                np.add.at(pixel_accumulator, (py[is_source], px[is_source]), source_texture[py_s, px_s])

            # 2. Celestial Sphere Hits
            is_sphere = rays['termination_type'] == cfg.TERM_SPHERE
            if np.any(is_sphere):
                sph_hits = rays[is_sphere]
                norm_phi = (sph_hits['final_phi'] / (2 * np.pi)) % 1.0 
                norm_theta = sph_hits['final_theta'] / np.pi
                
                px_sph = np.clip(norm_phi * (sphere_w - 1), 0, sphere_w - 1).astype(np.int32)
                py_sph = np.clip(norm_theta * (sphere_h - 1), 0, sphere_h - 1).astype(np.int32)
                
                np.add.at(pixel_accumulator, (py[is_sphere], px[is_sphere]), sphere_texture[py_sph, px_sph])

            np.add.at(count_accumulator, (py, px), 1)

    hit_mask = count_accumulator > 0
    final_image_float = np.zeros_like(pixel_accumulator)
    final_image_float[hit_mask] = pixel_accumulator[hit_mask] / count_accumulator[hit_mask, np.newaxis]
    
    img = Image.fromarray((np.clip(final_image_float, 0, 1) * 255).astype(np.uint8), 'RGB')
    os.makedirs(os.path.dirname(output_filename), exist_ok=True)
    img.save(output_filename)
    
    if display_image:
        img.show()
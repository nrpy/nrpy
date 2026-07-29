# nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py
"""
Central schema and configuration for geodesic-visualization post-processing.

This module defines Python-side binary layout matching C `blueprint_data_t`,
shared termination-type constants expected in serialized blueprint files, and
common visualization defaults used by renderer and diagnostic scripts.

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import struct

import numpy as np

# Native, same-build artifact contract. Cross-endian persistence is intentionally
# outside this contract. Schema v6 stores the two fields of view once per tile
# header; per-ray records retain event diagnostics, final-state metadata, and
# normalized image-sample coordinates.
# Nonterminal coordinates are crossing diagnostics; terminal coordinates are
# terminal-plane texture samples.
BLUEPRINT_MAGIC = b"NRPYBP01"
BLUEPRINT_SCHEMA_VERSION = 6
# Header stores tile identity/counts and FOV metadata. Runtime tile/full pixel
# dimensions are commondata used during initialization and are not serialized.
BLUEPRINT_HEADER_FORMAT = "=8sIIIIIIIQdd"
BLUEPRINT_HEADER_SIZE = 60
BLUEPRINT_RECORD_SIZE = 100
if struct.calcsize(BLUEPRINT_HEADER_FORMAT) != BLUEPRINT_HEADER_SIZE:
    raise RuntimeError("BLUEPRINT_HEADER_FORMAT does not match header size")

# Step 1: Core data structures.
# This dtype MUST match the 'blueprint_data_t' struct in the C code.
# It defines how individual ray results (endpoints, times, and types) are stored in
# binary format.
BLUEPRINT_DTYPE = np.dtype(
    [
        (
            "termination_type",
            "=i4",
        ),  # Enum indicating escape, terminal-plane hit, or integration failure
        ("y_nt", "=f8"),  # Local horizontal nonterminal-plane diagnostic
        ("z_nt", "=f8"),  # Local vertical nonterminal-plane diagnostic
        ("y_t", "=f8"),  # Local horizontal terminal-plane coordinate
        ("z_t", "=f8"),  # Local vertical terminal-plane coordinate
        ("final_theta", "=f8"),  # Final polar angle on the celestial sphere
        ("final_phi", "=f8"),  # Final azimuthal angle on the celestial sphere
        (
            "non_terminal_plane_lambda",
            "=f8",
        ),  # Affine parameter at nonterminal-plane intersection
        (
            "non_terminal_plane_t",
            "=f8",
        ),  # Coordinate time at nonterminal-plane intersection
        ("L_f", "=f8"),  # Physical affine parameter when the photon terminated
        ("t_f", "=f8"),  # Coordinate time when the photon terminated
        ("image_width_fraction", "=f8"),  # Normalized width sample coordinate
        ("image_height_fraction", "=f8"),  # Normalized height sample coordinate
    ],
    align=False,
)
if BLUEPRINT_DTYPE.itemsize != BLUEPRINT_RECORD_SIZE:
    raise RuntimeError("BLUEPRINT_DTYPE does not match blueprint_data_t size")
BLUEPRINT_FIELDS = BLUEPRINT_DTYPE.fields
assert BLUEPRINT_FIELDS is not None
if BLUEPRINT_FIELDS["termination_type"][1] != 0:
    raise RuntimeError("termination_type offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["y_nt"][1] != 4:
    raise RuntimeError("y_nt offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["t_f"][1] != 76:
    raise RuntimeError("t_f offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["image_width_fraction"][1] != 84:
    raise RuntimeError("image_width_fraction offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["image_height_fraction"][1] != 92:
    raise RuntimeError("image_height_fraction offset changed in BLUEPRINT_DTYPE")
BLUEPRINT_NORM_ABS_DTYPE = np.dtype("=f8")
BLUEPRINT_NORM_ABS_FILENAME_TEMPLATE = (
    "light_blueprint_norm_abs_{tile_x:02d}_{tile_y:02d}.bin"
)

# Step 2: Termination enums.
# These integers identify the fate of a photon ray.
# They must remain synchronized with 'termination_type_t' in the C-header files.
STOP_CONDITION_COORD_RADIUS_EXCEEDED = 0  # Coordinate-radius stop condition
STOP_CONDITION_TERMINAL_PLANE = 1  # Terminal-plane stop condition
STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED = 2  # Evolution-measure stop condition
FAILURE_RKF45_REJECTION_LIMIT = 3  # RKF45 rejected too many consecutive steps
STOP_CONDITION_T_MAX_EXCEEDED = 4  # Maximum coordinate-time stop condition
FAILURE_SLOT_MANAGER_ERROR = 5  # Slot manager failed to handle the ray
FAILURE_GENERIC = 6  # Unspecified integration failure
ACTIVE = 7  # Ray is still being processed (should not appear in final blueprints)
REJECTED = 8  # Ray is in a rejected RKF45 stage (not a final status)

# Step 3: Physics and scene parameters.
MASS_OF_BLACK_HOLE = 1.0  # Normalized mass ($M$)

# Step 4: Texture and disk generation parameters.
SPHERE_TEXTURE_FILE = "noirlab2430b.tif"  # Background image for escaped rays
DISK_INNER_RADIUS = 6.0  # Inner edge of the disk (usually near ISCO)
DISK_OUTER_RADIUS = 25.0  # Outer edge of the disk
COLORMAP = "afmhot"  # Matplotlib colormap for disk temperature
DISK_TEMP_POWER_LAW = -1.5  # Radial temperature decay: $T \propto r^{power}$
SOURCE_PHYSICAL_WIDTH = 2 * DISK_OUTER_RADIUS  # Diameter of terminal-plane disk texture

# Step 5: Rendering parameters.
STATIC_IMAGE_PIXEL_WIDTH = 700  # Resolution of the final lensed image
CHUNK_SIZE = 10_000_000  # Number of rays to process in memory at once

if __name__ == "__main__":
    import doctest
    import sys

    results = doctest.testmod()

    if results.failed > 0:
        print(f"Doctest failed: {results.failed} of {results.attempted} test(s)")
        sys.exit(1)
    else:
        print(f"Doctest passed: All {results.attempted} test(s) passed")

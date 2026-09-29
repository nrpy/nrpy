# nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py
"""
Field definitions and configuration for geodesic-visualization post-processing.

This module defines Python-side binary layout matching C `blueprint_data_t`,
shared termination-type constants expected in serialized blueprint files, and
common visualization defaults used by renderer and diagnostic scripts.

Final blueprint records contain physical stop conditions and terminal failures,
never the internal active or rejected RKF45 states.

>>> FINAL_TERMINATION_TYPES
(0, 1, 2, 3, 4, 5, 6, 9, 10, 11)
>>> ACTIVE in FINAL_TERMINATION_TYPES or REJECTED in FINAL_TERMINATION_TYPES
False

Author: Dalton J. Moone
        daltonmoone **at** gmail **dot** com
"""

import struct

import numpy as np

# Native same-build binary layout. Cross-endian persistence is intentionally
# unsupported. Version 7 appends axial angular momentum, coordinate time, and
# signed plane distance at the first accepted state after a nonterminal crossing.
# Version 6 contains the preceding fields and remains readable.
# Nonterminal coordinates are crossing diagnostics; terminal coordinates are
# terminal-plane texture samples.
BLUEPRINT_MAGIC = b"NRPYBP01"
BLUEPRINT_SCHEMA_VERSION = 7
# Header stores tile identity/counts and FOV metadata. Runtime tile/full pixel
# dimensions are commondata used during initialization and are not serialized.
BLUEPRINT_HEADER_FORMAT = "=8sIIIIIIIQdd"
BLUEPRINT_HEADER_SIZE = 60
BLUEPRINT_RECORD_SIZE_V6 = 100
BLUEPRINT_RECORD_SIZE = 124
if struct.calcsize(BLUEPRINT_HEADER_FORMAT) != BLUEPRINT_HEADER_SIZE:
    raise RuntimeError("BLUEPRINT_HEADER_FORMAT does not match header size")

# Step 1: Core data structures.
# This dtype MUST match the 'blueprint_data_t' struct in the C code.
# It defines how individual ray results (endpoints, times, and types) are stored in
# binary format.
BLUEPRINT_DTYPE_V6 = np.dtype(
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
BLUEPRINT_DTYPE = np.dtype(
    BLUEPRINT_DTYPE_V6.descr
    + [
        ("non_terminal_post_step_Lz", "=f8"),
        ("non_terminal_post_step_t", "=f8"),
        ("non_terminal_post_step_distance", "=f8"),
    ],
    align=False,
)
BLUEPRINT_DTYPES_BY_VERSION = {6: BLUEPRINT_DTYPE_V6, 7: BLUEPRINT_DTYPE}
if BLUEPRINT_DTYPE_V6.itemsize != BLUEPRINT_RECORD_SIZE_V6:
    raise RuntimeError("BLUEPRINT_DTYPE_V6 does not match version-6 record size")
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
if BLUEPRINT_FIELDS["non_terminal_post_step_Lz"][1] != 100:
    raise RuntimeError("non_terminal_post_step_Lz offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["non_terminal_post_step_t"][1] != 108:
    raise RuntimeError("non_terminal_post_step_t offset changed in BLUEPRINT_DTYPE")
if BLUEPRINT_FIELDS["non_terminal_post_step_distance"][1] != 116:
    raise RuntimeError(
        "non_terminal_post_step_distance offset changed in BLUEPRINT_DTYPE"
    )
BLUEPRINT_NORM_ABS_DTYPE = np.dtype("=f8")
BLUEPRINT_NORM_ABS_FILENAME_TEMPLATE = (
    "light_blueprint_norm_abs_{tile_x:02d}_{tile_y:02d}.bin"
)
PLANE_CROSSING_RECORD_SIZE = 108
PLANE_CROSSING_DTYPE = np.dtype(
    [
        ("photon_index", "=u8"),
        ("interpolation_degree", "=u4"),
        ("integration_param", "=f8"),
        ("y_local", "=f8"),
        ("z_local", "=f8"),
        ("state", "=f8", (9,)),
    ],
    align=False,
)
if PLANE_CROSSING_DTYPE.itemsize != PLANE_CROSSING_RECORD_SIZE:
    raise RuntimeError(
        "PLANE_CROSSING_DTYPE does not match plane_crossing_record_t size"
    )
NON_TERMINAL_CROSSINGS_FILENAME_TEMPLATE = (
    "light_blueprint_non_terminal_crossings_{tile_x:02d}_{tile_y:02d}.bin"
)
TERMINAL_CROSSINGS_FILENAME_TEMPLATE = (
    "light_blueprint_terminal_crossings_{tile_x:02d}_{tile_y:02d}.bin"
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
FAILURE_SPATIAL_INTERPOLATION = 9  # Spatial interpolation failed for this ray
FAILURE_TEMPORAL_INTERPOLATION = 10  # Temporal interpolation failed for this ray
FAILURE_PLANE_INTERPOLATION_HISTORY = (
    11  # Too few distinct accepted states for a quadratic crossing
)

# Only completed physical stops and failures may be serialized. ACTIVE and
# REJECTED are internal RKF45 states and therefore deliberately absent.
FINAL_TERMINATION_TYPES = (
    STOP_CONDITION_COORD_RADIUS_EXCEEDED,
    STOP_CONDITION_TERMINAL_PLANE,
    STOP_CONDITION_EVOLUTION_MEASURE_EXCEEDED,
    FAILURE_RKF45_REJECTION_LIMIT,
    STOP_CONDITION_T_MAX_EXCEEDED,
    FAILURE_SLOT_MANAGER_ERROR,
    FAILURE_GENERIC,
    FAILURE_SPATIAL_INTERPOLATION,
    FAILURE_TEMPORAL_INTERPOLATION,
    FAILURE_PLANE_INTERPOLATION_HISTORY,
)

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

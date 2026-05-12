# nrpy/examples/geodesic_visualizations/blueprint_config_and_schema.py
"""
Central schema and configuration for geodesic-visualization post-processing.

This module defines Python-side binary layout matching C `blueprint_data_t`,
shared termination-type constants expected in serialized blueprint files, and
common visualization defaults used by renderer and diagnostic scripts.

Author: Dalton J. Moone
"""

import numpy as np

# Step 1: Core data structures.
# This dtype MUST match the 'blueprint_data_t' struct in the C code.
# It defines how individual ray results (endpoints, times, and types) are stored in
# binary format.
BLUEPRINT_DTYPE = np.dtype(
    [
        (
            "termination_type",
            np.int32,
        ),  # Enum indicating where the ray ended (e.g., Sphere, Disk)
        ("y_w", "f8"),  # Local horizontal coordinate on the camera window
        ("z_w", "f8"),  # Local vertical coordinate on the camera window
        ("y_s", "f8"),  # Horizontal coordinate on the source (accretion disk) plane
        ("z_s", "f8"),  # Vertical coordinate on the source (accretion disk) plane
        ("final_theta", "f8"),  # Final polar angle on the celestial sphere
        ("final_phi", "f8"),  # Final azimuthal angle on the celestial sphere
        ("L_w", "f8"),  # Affine parameter at window intersection
        ("t_w", "f8"),  # Coordinate time at window intersection
        ("L_s", "f8"),  # Affine parameter at source plane intersection
        ("t_s", "f8"),  # Coordinate time at source plane intersection
    ],
    align=False,
)

# Step 2: Termination enums.
# These integers identify the fate of a photon ray.
# They must remain synchronized with 'termination_type_t' in the C-header files.
TERM_SPHERE = 0  # Ray escaped to the far-field (Celestial Sphere)
TERM_SOURCE_PLANE = 1  # Ray hit the accretion disk / source plane
TERM_FAIL_PT_BIG = 2  # Numerical failure: momentum component $p_t$ exceeded bounds
TERM_FAIL_RKF45 = 3  # Numerical failure: RKF45 integrator could not converge
TERM_FAIL_T_MAX = 4  # Ray exceeded the maximum allowed integration time
TERM_FAIL_SLOT = 5  # Slot manager error (temporal binning failure)
TERM_FAIL_GENERIC = 6  # Unspecified integration failure
TERM_ACTIVE = 7  # Ray is still being processed (should not appear in final blueprints)

# Step 3: Physics and scene parameters.
MASS_OF_BLACK_HOLE = 1.0  # Normalized mass ($M$)
WINDOW_WIDTH = 1.0  # Physical width of the camera's projection window
WINDOW_HEIGHT = 1.0  # Physical height of the camera's projection window

# Step 4: Texture and disk generation parameters.
SPHERE_TEXTURE_FILE = "noirlab2430b.tif"  # Background image for escaped rays
DISK_INNER_RADIUS = 6.0  # Inner edge of the disk (usually near ISCO)
DISK_OUTER_RADIUS = 25.0  # Outer edge of the disk
COLORMAP = "afmhot"  # Matplotlib colormap for disk temperature
DISK_TEMP_POWER_LAW = -1.5  # Radial temperature decay: $T \propto r^{power}$
SOURCE_PHYSICAL_WIDTH = 2 * DISK_OUTER_RADIUS  # Total diameter of the source image

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

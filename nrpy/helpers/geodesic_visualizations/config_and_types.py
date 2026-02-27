"""
PIPELINE CONFIGURATION AND DATA STRUCTURE DEFINITIONS
This script serves as the "source of truth" for both the Python rendering scripts and 
the C-based numerical integrator. It defines the structured NumPy data type (DTYPE) 
required to parse binary ray-tracing results and establishes global constants for 
physics, file paths, and texture parameters.
"""

import numpy as np
import os

# --- 1. Core Data Structures ---
# This dtype MUST match the 'blueprint_data_t' struct in the C code.
# It defines how individual ray results (endpoints, times, and types) are stored in binary.
BLUEPRINT_DTYPE = np.dtype([
    ('termination_type', np.int32), # Enum indicating where the ray ended (e.g., Sphere, Disk)
    ('y_w', 'f8'),                 # Local horizontal coordinate on the camera window
    ('z_w', 'f8'),                 # Local vertical coordinate on the camera window
    ('y_s', 'f8'),                 # Horizontal coordinate on the source (accretion disk) plane
    ('z_s', 'f8'),                 # Vertical coordinate on the source (accretion disk) plane
    ('final_theta', 'f8'),         # Final polar angle on the celestial sphere
    ('final_phi', 'f8'),           # Final azimuthal angle on the celestial sphere
    ('L_w', 'f8'),                 # Affine parameter at window intersection
    ('t_w', 'f8'),                 # Coordinate time at window intersection
    ('L_s', 'f8'),                 # Affine parameter at source plane intersection
    ('t_s', 'f8'),                 # Coordinate time at source plane intersection
], align=False)

# --- 2. Termination Enums ---
# These integers identify the fate of a photon ray.
# They must remain synchronized with 'termination_type_t' in the C-header files.
TERM_SPHERE = 0           # Ray escaped to the far-field (Celestial Sphere)
TERM_SOURCE_PLANE = 1     # Ray hit the accretion disk / source plane
TERM_FAIL_PT_BIG = 2      # Numerical failure: momentum component $p_t$ exceeded bounds
TERM_FAIL_RKF45 = 3       # Numerical failure: RKF45 integrator could not converge
TERM_FAIL_T_MAX = 4       # Ray exceeded the maximum allowed integration time
TERM_FAIL_SLOT = 5        # Slot manager error (temporal binning failure)
TERM_FAIL_GENERIC = 6     # Unspecified integration failure
TERM_ACTIVE = 7           # Ray is still being processed (should not appear in final blueprints)

# --- 3. Global Directories ---
# Standardizes paths for project assets and generated output.
HOME_DIR = os.path.expanduser('~')
BASE_PROJECT_DIR = os.path.join(HOME_DIR, "Desktop", "Test_PR", "nrpy", "project")
LIGHT_INTEGRATOR_DIR = os.path.join(BASE_PROJECT_DIR, "photon_geodesic_integrator")
OUTPUT_BASEDIR = os.path.join(HOME_DIR, "Desktop", "Test_PR", "nrpy", "Generated_nrpy_images")

# --- 4. Physics & Scene Parameters ---
MASS_OF_BLACK_HOLE = 1.0 # Normalized mass ($M$)
WINDOW_WIDTH = 1.0       # Physical width of the camera's projection window

# --- 5. Texture & Disk Generation Parameters ---
SPHERE_TEXTURE_FILE = "starmap_2020.png" # Background image for escaped rays
DISK_INNER_RADIUS = 6.0                  # Inner edge of the disk (usually near ISCO)
DISK_OUTER_RADIUS = 25.0                 # Outer edge of the disk
COLORMAP = 'afmhot'                      # Matplotlib colormap for disk temperature
DISK_TEMP_POWER_LAW = -1.5               # Radial temperature decay: $T \propto r^{power}$
SOURCE_PHYSICAL_WIDTH = 2 * DISK_OUTER_RADIUS # Total diameter of the source image

# --- 6. Rendering Parameters ---
STATIC_IMAGE_PIXEL_WIDTH = 700           # Resolution of the final lensed image
CHUNK_SIZE = 10_000_000                  # Number of rays to process in memory at once
"""
Common parameters for NRPyElliptic.

Author: Thiago Assumpção
        assumpcaothiago **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed NRPy+ core modules:
import nrpy.params as par  # NRPy+: Parameter interface

# Parameters common to/needed by all NRPyElliptic Python modules

# Step P1: Define time parameter
time = par.register_CodeParameter("REAL", __name__, "time", 0.0)

# Step P2.a: Define the damping parameter used in the hyperbolic relaxation method
eta_damping = par.register_CodeParameter(
    "REAL", __name__, "eta_damping", 10.0, commondata=True
)

# Step P2.b: Define the parameter that controls minimum value of wavespeed in the numerical domain
MINIMUM_GLOBAL_WAVESPEED = par.register_CodeParameter(
    "REAL", __name__, "MINIMUM_GLOBAL_WAVESPEED", 0.7, commondata=True
)

# Step P3: Define puncture parameters

# Step P3.a: bare masses
bare_mass_0 = par.register_CodeParameter(
    "REAL", __name__, "bare_mass_0", 0.0, commondata=True
)
bare_mass_1 = par.register_CodeParameter(
    "REAL", __name__, "bare_mass_1", 0.0, commondata=True
)

# Step P3.b: position of the punctures in the z axis
zPunc = par.register_CodeParameter("REAL", __name__, "zPunc", 5.0, commondata=True)

# Step P3.c.1: linear momentum 0
P0U = par.register_CodeParameters(
    "REAL", __name__, ["P0_x", "P0_y", "P0_z"], 0.0, commondata=True
)
# Step P3.c.2: linear momentum 1
P1U = par.register_CodeParameters(
    "REAL", __name__, ["P1_x", "P1_y", "P1_z"], 0.0, commondata=True
)

# Step P3.d.1: angular momentum 0
S0U = par.register_CodeParameters(
    "REAL", __name__, ["S0_x", "S0_y", "S0_z"], 0.0, commondata=True
)
# Step P3.d.2: angular momentum 1
S1U = par.register_CodeParameters(
    "REAL", __name__, ["S1_x", "S1_y", "S1_z"], 0.0, commondata=True
)

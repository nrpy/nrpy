"""
Common parameters for scalar wave evolutions

Author: Zachariah B. Etienne
        zachetie **at** gmail **dot* com

License: BSD 2-Clause
"""

# Step P1: Import needed NRPy+ core modules:
import nrpy.params as par  # NRPy+: Parameter interface

thismodule = __name__

# Parameters common to/needed by all WaveEquation Python modules

# Step P2: Define the C parameter wavespeed. The `wavespeed`
#          variable is a proper SymPy variable, so it can be
#          used in below expressions. In the C code, it acts
#          just like a usual parameter, whose value is
#          specified in the parameter file.
wavespeed = par.register_CodeParameter(
    "REAL", thismodule, "wavespeed", 1.0, commondata=True
)

"""
Generate input sets for sebob consistency check.

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
"""

import numpy as np

NUM_SETS = 10000
NUM_INPUTS = 6
Q_RANGE = np.linspace(1.0, 10.0, NUM_SETS)
np.random.shuffle(Q_RANGE)
CHI1_RANGE = np.linspace(-0.99, 0.99, NUM_SETS)
np.random.shuffle(CHI1_RANGE)
CHI2_RANGE = np.linspace(-0.99, 0.99, NUM_SETS)
np.random.shuffle(CHI2_RANGE)
OMEGA_RANGE = np.linspace(0.01, 0.015, NUM_SETS)
np.random.shuffle(OMEGA_RANGE)
M = np.tile(50, NUM_SETS)
DT = np.tile(2.4627455127717882e-05, NUM_SETS)
OUTPUT_FILE = "input_sets.csv"
print(f"Generating {NUM_SETS} input sets...")
input_data = np.column_stack((Q_RANGE, CHI1_RANGE, CHI2_RANGE, OMEGA_RANGE, M, DT))
np.savetxt(OUTPUT_FILE, input_data, delimiter=",", fmt="%.18e")
print(f"Saved input sets to '{OUTPUT_FILE}'")

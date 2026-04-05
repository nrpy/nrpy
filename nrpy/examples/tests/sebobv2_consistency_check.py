"""
Set up on-the-fly accuracy comparison for sebobv2.
usage: sebobv2_consistency_check.py [-h] --current-exec CURRENT_EXEC --trusted-exec TRUSTED_EXEC

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
"""

import os
import subprocess
import sys
from io import StringIO
from pathlib import Path
from typing import Any, Tuple, Union

import numpy as np
from numpy.typing import NDArray

# --- Configuration ---
PERTURBATION_MAGNITUDE = 1e-14


# --- Helper Functions ---
def run_sebobv2(
    executable_path: str, inputs: NDArray[np.float64]
) -> NDArray[np.float64]:
    """
    Run sebob executable with a single set of inputs.

    :param executable_path: Path to the directory containing the executable;
                                the executable is expected to be named the same as this directory.
    :param inputs: List of inputs to the sebob executable.
    :return: Output of the sebob executable.
    """
    executable_dir = Path(executable_path)
    approximant = executable_dir.name
    executable_file = str(executable_dir / approximant)
    parameters_file = executable_file + ".par"
    inputs_str = [f"{elt:.15f}" for elt in inputs]
    result = subprocess.run(
        [executable_file, parameters_file] + inputs_str,
        capture_output=True,
        text=True,
        check=True,
    )

    return np.loadtxt(StringIO(result.stdout))


def calculate_rmse(
    data1: NDArray[np.float64], data2: NDArray[np.float64]
) -> Union[float, Any]:
    """
    Calculate the Root Mean Square Error (RMSE) between two datasets.

    :param data1: First dataset.
    :param data2: Second dataset.
    :return: RMSE value.
    """
    h22_1 = data1[:, 1] + 1j * data1[:, 2]
    amp1, phase1 = np.abs(h22_1), np.unwrap(np.angle(h22_1))
    h22_2 = data2[:, 1] + 1j * data2[:, 2]
    amp2, phase2 = np.abs(h22_2), np.unwrap(np.angle(h22_2))
    t_min = max(data1[0, 0], data2[0, 0])
    t_max = min(data1[-1, 0], data2[-1, 0])
    # sample at 10% of the total number of points
    sampled_times = np.linspace(t_min, t_max, data1.shape[0] // 10)
    sampled_amp_data1 = np.interp(sampled_times, data1[:, 0], amp1)
    sampled_amp_data2 = np.interp(sampled_times, data2[:, 0], amp2)
    sampled_phase_data1 = np.interp(sampled_times, data1[:, 0], phase1)
    sampled_phase_data2 = np.interp(sampled_times, data2[:, 0], phase2)
    rmse_amp = np.sqrt(
        np.mean(((sampled_amp_data1 - sampled_amp_data2) / sampled_amp_data1) ** 2)
    )
    rmse_phase = np.sqrt(
        np.mean(
            ((sampled_phase_data1 - sampled_phase_data2) / sampled_phase_data1) ** 2
        )
    )
    return rmse_amp + rmse_phase


def process_input_set(
    nominal_args: Tuple[NDArray[np.float64], str, str],
) -> Tuple[Union[float, Any], Union[float, Any]]:
    """
    Process a single input set to get both baseline and test error.

    :param nominal_args: Tuple containing the nominal input paramters, path to trusted executable, and path to current executable.
    :return: Tuple containing the baseline error and test error.
    """
    nominal_inputs, nominal_trusted_exec, nominal_current_exec = nominal_args

    # 1. Run trusted code with nominal inputs
    trusted_output = run_sebobv2(nominal_trusted_exec, nominal_inputs)
    # 2. Run current code with nominal inputs
    current_output = run_sebobv2(nominal_current_exec, nominal_inputs)

    # 3. Create perturbed inputs only for mass ratio and spins and run trusted code again
    rng_perturbation = np.random.default_rng(0)
    perturbation = (
        rng_perturbation.choice([-1, 1], size=3, replace=True)
        * rng_perturbation.uniform(1, 3, size=3)
        * PERTURBATION_MAGNITUDE
    )
    perturbed_inputs = nominal_inputs.copy()
    perturbed_inputs[0] = nominal_inputs[0] * (1 + perturbation[0])
    perturbed_inputs[1] = nominal_inputs[1] * (1 + perturbation[1])
    perturbed_inputs[2] = nominal_inputs[2] * (1 + perturbation[2])
    perturbed_output = run_sebobv2(nominal_trusted_exec, perturbed_inputs)

    # Calculate errors
    baseline_error = calculate_rmse(trusted_output, perturbed_output)
    test_error = calculate_rmse(trusted_output, current_output)
    return baseline_error, test_error


# --- Main Logic ---
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Run sebob consistency check.")
    parser.add_argument(
        "--trusted-exec",
        type=str,
        required=True,
        help="Path to trusted sebob executable.",
    )
    parser.add_argument(
        "--current-exec",
        type=str,
        required=True,
        help="Path to current sebob executable.",
    )
    args = parser.parse_args()
    num_sets = 10
    rng = np.random.default_rng(0)
    q = np.linspace(1.01, 10, num_sets)
    rng.shuffle(q)
    chi_1 = np.linspace(-0.9, 0.9, num_sets)
    rng.shuffle(chi_1)
    chi_2 = np.linspace(-0.9, 0.9, num_sets)
    rng.shuffle(chi_2)
    M = 50
    omega_0 = 0.011
    dt = 2.4627455127717882e-05

    cdir = os.getcwd()
    # go to the directory where the trusted sebob executable is located
    os.chdir(args.trusted_exec)
    subprocess.run(
        ["make", "clean"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    subprocess.run(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    os.chdir(cdir)
    # go to the directory where the current sebob executable is located
    os.chdir(args.current_exec)
    subprocess.run(
        ["make", "clean"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    subprocess.run(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    os.chdir(cdir)
    trusted_exec = args.trusted_exec
    current_exec = args.current_exec
    print(f"Starting accuracy comparison for {num_sets} input sets...")
    baseline_errors = []
    test_errors = []
    for i in range(num_sets):
        print(f"Processing input set {i+1}/{num_sets}...")
        inputs_set = np.array([q[i], chi_1[i], chi_2[i], omega_0, M, dt])
        task = (inputs_set, trusted_exec, current_exec)
        baseline_err, test_err = process_input_set(task)
        baseline_errors.append(baseline_err)
        test_errors.append(test_err)

    baseline_median = np.median(baseline_errors)
    test_median = np.median(test_errors)

    print("\n--- Test Results ---")
    print(f"Test Error Median:      {test_median:.6e}")
    print(f"Baseline Error Median:  {baseline_median:.6e}")
    if test_median <= baseline_median:
        print("\nPASSED: Median error is within roundoff baseline.\n")
        sys.exit(0)
    else:
        print(
            f"\nFAILED: Median error ({test_median:.6e}) exceeds the baseline ({baseline_median:.6e}).\n"
        )
        sys.exit(1)

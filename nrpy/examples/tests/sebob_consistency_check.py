"""
Set up on-the-fly accuracy comparison for sebob.
usage: sebob_consistency_check.py [-h] --current-exec CURRENT_EXEC --trusted-exec TRUSTED_EXEC --inputs INPUTS

Authors: Siddharth Mahesh
        sm0193 **at** mix **dot** wvu **dot** edu
"""

import os
import subprocess
import sys
from io import StringIO
from typing import Any, Tuple, Union

import numpy as np
from numpy.typing import NDArray

# --- Configuration ---
PERTURBATION_MAGNITUDE = 1e-14


# --- Helper Functions ---
def run_sebob(executable_path: str, inputs: NDArray[np.float64]) -> NDArray[np.float64]:
    """
    Run sebob executable with a single set of inputs.

    :param executable_path: Path to the sebob executable.
    :param inputs: List of inputs to the sebob executable.
    :return: Output of the sebob executable.
    """
    parameters_file = executable_path + ".par"
    inputs_str = [f"{elt:.15f}" for elt in inputs]
    result = subprocess.run(
        [executable_path, parameters_file] + inputs_str,
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
    t_max = min(data1[0, -1], data2[0, -1])
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
    trusted_output = run_sebob(nominal_trusted_exec, nominal_inputs)
    # 2. Run current code with nominal inputs
    current_output = run_sebob(nominal_current_exec, nominal_inputs)

    # 3. Create perturbed inputs only for mass ratio and spins and run trusted code again
    perturbation = (
        np.random.choice([-1, 1], size=3, replace=True)
        * np.random.uniform(1, 3, size=3)
        * PERTURBATION_MAGNITUDE
    )
    perturbed_inputs = nominal_inputs.copy()
    perturbed_inputs[0] = nominal_inputs[0] * (1 + perturbation[0])
    perturbed_inputs[1] = nominal_inputs[1] * (1 + perturbation[1])
    perturbed_inputs[2] = nominal_inputs[2] * (1 + perturbation[2])
    perturbed_output = run_sebob(nominal_trusted_exec, perturbed_inputs)

    # Calculate errors
    baseline_error = calculate_rmse(trusted_output, perturbed_output)
    test_error = calculate_rmse(trusted_output, current_output)
    return baseline_error, test_error


# --- Main Logic ---
if __name__ == "__main__":
    # equal mass non-spinning binary with ~20M separation and ~.1M spacing
    inputs_set = np.array([[1.0, 0.0, 0.0, 0.0118, 50.0, 2.4627455127717882e-05]])
    num_sets = len(inputs_set)
    cdir = os.getcwd()
    # go to the directory where the trusted sebob executable is located
    os.chdir("./project/seobnrv5_aligned_spin_inspiral")
    subprocess.run(
        ["make", "clean"],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    subprocess.run(["make"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
    # go to the directory where the current sebob executable is located
    # for now, using the same trusted and current executables
    # the next PR will have the trusted executable in a different directory
    os.chdir(cdir)
    trusted_exec = (
        "./project/seobnrv5_aligned_spin_inspiral/seobnrv5_aligned_spin_inspiral"
    )
    current_exec = (
        "./project/seobnrv5_aligned_spin_inspiral/seobnrv5_aligned_spin_inspiral"
    )
    print(f"Starting accuracy comparison for {num_sets} input sets...")
    baseline_errors = []
    test_errors = []
    for i in range(num_sets):
        print(f"Processing input set {i+1}/{num_sets}...")
        task = (inputs_set[i], trusted_exec, current_exec)
        baseline_err, test_err = process_input_set(task)
        baseline_errors.append(baseline_err)
        test_errors.append(test_err)

    # Final comparison
    baseline_median = np.median(baseline_errors)
    test_median = np.median(test_errors)

    print("\n--- Test Results ---")
    print(f"Test Error Median:      {test_median:.6e}")
    print(f"Baseline Error Median:  {baseline_median:.6e}")
    if test_median <= baseline_median:
        print("\nPASSED: Median error is within roundoff baseline.")
        sys.exit(0)
    else:
        print(
            f"\nFAILED: Median error ({test_median:.6e}) exceeds the baseline ({baseline_median:.6e})."
        )
        sys.exit(1)

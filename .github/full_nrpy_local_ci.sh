#!/bin/bash

# Exit immediately if a command exits with a non-zero status, except in specific cases
set -e

# Function to handle each test step
run_test_step() {
  local step_name="$1"
  local command="$2"
  local failure_message="$3"

  echo "-={ Step $step_counter: $step_name }=-"
  eval "$command"
  if [ $? -ne 0 ]; then
    failed_tests+=("$failure_message in $python_file")
  fi
  ((step_counter++))
}

# Clear cache
echo "Clearing cache..."
rm -rf ~/.cache/nrpy
echo "Cache cleared."

# Display Python version
python_version=$(python --version 2>&1)
echo "Python version: $python_version"

# Install dependencies
echo "Installing dependencies..."
python -m pip install --upgrade pip
pip install -U -r requirements.txt

if [[ "$python_version" != *"3.6.7"* && "$python_version" != *"3.7.13"* ]]; then
  pip install -U -r requirements-dev.txt
else
  # Install packages from requirements-dev.txt except those with a fixed version.
  grep -v '==' requirements-dev.txt | xargs pip install -U
  pip install -U clang-format
fi
echo "Dependencies installed."

# Install NRPyLaTeX and additional dependencies
echo "Installing NRPyLaTeX and additional dependencies..."
pip install -U nrpylatex
pip install -U ipython setuptools
echo "NRPyLaTeX and additional dependencies installed."

# Display SymPy version
echo "Running CI tests with SymPy version = $(isympy --version)"

# Initialize array to store failed tests
failed_tests=()

# Find all relevant Python files recursively and store them in a temporary file
echo "Locating Python files..."
find . -type f -name '*.py' \
  ! -name '__init__.py' \
  ! -path './project/*' \
  ! -path './build/*' \
  ! -path '*/tests/*' \
  ! -path './nrpy/examples/visualization_scripts/*' > python_files.txt

# Count total number of Python files
total_files=$(wc -l < python_files.txt)
echo "Found $total_files Python file(s) to analyze."

# Initialize file counter
file_counter=1

# Loop through each Python file
while IFS= read -r python_file; do
  echo ""

  # Initialize step counter for each file
  step_counter=1

  echo "-={ $python_file : file $file_counter of $total_files for static analysis }=-"

  # Run Doctests if not in ./ (root directory) or nrpy/examples/
  if [[ ! $python_file =~ ^./[^/]*\.py && ! $python_file =~ nrpy/examples/.* ]]; then
    run_test_step "Doctests" "PYTHONPATH=.:$PYTHONPATH python \"$python_file\"" "doctests"
  fi


  # Conditional static analysis for specific Python versions
  if [[ "$python_version" != *"3.6.7"* && "$python_version" != *"3.7.13"* && "$python_version" != *"3.8.18"* ]]; then
    # Step: black
    run_test_step "black" "black --check \"$python_file\"" "black"

    # Step: isort
    run_test_step "isort" "isort --check-only \"$python_file\"" "isort"

    # Step: mypy
    run_test_step "mypy" "PYTHONPATH=.:$PYTHONPATH mypy --strict --pretty --allow-untyped-calls \"$python_file\"" "mypy"
  fi

  # Step: pylint
  echo "-={ Step $step_counter: pylint }=-"
  if [[ "$python_version" == *"3.6.7"* || "$python_version" == *"3.7.13"* ]]; then
    pylint_rcfile=".pylintrc_python36"
  else
    pylint_rcfile=".pylintrc"
  fi

  # Run pylint and capture the output, prevent 'set -e' from exiting on failure
  pylint_output=$(PYTHONPATH=.:$PYTHONPATH pylint --rcfile="$pylint_rcfile" "$python_file" 2>&1 || true)

  # Extract the Pylint score by searching for "rated at" and extracting the number before "/10"
  pylint_score=$(echo "$pylint_output" | grep "rated at" | grep -Eo '[0-9]+\.[0-9]+' | head -1 || echo "0")

  echo "Pylint score is $pylint_score"

  if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
    echo "$pylint_output"
    echo "Pylint score is below 9.91, failing..."
    failed_tests+=("pylint in $python_file")
  fi
  ((step_counter++))

  # Step: pydocstyle
  run_test_step "pydocstyle" "pydocstyle \"$python_file\"" "pydocstyle"

  # Step: darglint
  run_test_step "darglint" "darglint -v 2 \"$python_file\"" "darglint"

  # Increment file counter
  ((file_counter++))

done < python_files.txt

# Remove temporary file
rm python_files.txt

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo ""
  echo "The following tests failed: ${failed_tests[*]}"
  exit 1
fi

echo ""
echo "All static analysis tests passed successfully."

# Run example scripts and compile related projects
echo ""
echo "Running example scripts and compiling related projects..."

example_scripts=(
  "nrpy/examples/wave_equation_cartesian.py project/wavetoy"
  "nrpy/examples/wave_equation_curvilinear.py project/curviwavetoy"
  "nrpy/examples/nrpyelliptic_conformally_flat.py project/nrpyelliptic_conformally_flat"
  "nrpy/examples/two_blackholes_collide.py project/two_blackholes_collide"
  "nrpy/examples/blackhole_spectroscopy.py project/blackhole_spectroscopy"
  "nrpy/examples/spinning_blackhole.py project/spinning_blackhole"
  "nrpy/examples/nrpypn_quasicircular_momenta.py project/nrpypn_quasicircular_momenta"
  "nrpy/examples/wave_equation_multicoord_wavetoy.py project/multicoords_curviwavetoy"
  "nrpy/examples/carpet_baikal_thorns.py project/carpet_baikal_thorns"
  "nrpy/examples/carpet_wavetoy_thorns.py project/carpet_wavetoy_thorns"
  "nrpy/examples/carpetx_baikal_thorns.py project/carpetx_baikal_thorns"
  "nrpy/examples/carpetx_wavetoy_thorns.py project/carpetx_wavetoy_thorns"
  "nrpy/examples/seobnrv5_aligned_spin_inspiral.py project/seobnrv5_aligned_spin_inspiral"
  "nrpy/examples/superB_two_blackholes_collide.py project/superB_two_blackholes_collide"
  "nrpy/examples/superB_blackhole_spectroscopy.py project/superB_blackhole_spectroscopy"
  "nrpy/examples/tovola_neutron_star.py project/tovola_neutron_star"
  "nrpy/examples/hydro_without_hydro.py project/hydro_without_hydro"
  "nrpy/examples/manga_bhah_lib.py project/bhah_lib"
)

for script in "${example_scripts[@]}"; do
  IFS=' ' read -r script_path project_path <<< "$script"
  echo ""
  echo "Running Python script: $script_path"
  PYTHONPATH=.:$PYTHONPATH python "$script_path"
  if [[ $? -ne 0 ]]; then
    echo "Error: Python script $script_path failed."
    exit 1
  fi

  if [[ $script_path != *"superB"* && $script_path != *"carpet"* ]]; then
    echo "Compiling project in $project_path..."
    (cd "$project_path" && make && make clean)
    if [[ $? -ne 0 ]]; then
      echo "Error: Compilation in $project_path failed."
      exit 1
    fi
    echo "Compilation in $project_path succeeded."
  fi
done

echo ""
echo "All example scripts executed and related projects compiled successfully."

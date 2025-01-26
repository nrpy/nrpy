#!/bin/bash

# Exit immediately if a command exits with a non-zero status
set -e

# Ensure the script is called with a single file argument
if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <python_file>"
  exit 1
fi

python_file=$1

# Check if the file exists
if [ ! -f "$python_file" ]; then
  echo "File '$python_file' does not exist."
  exit 1
fi

# Change to the project's root directory if necessary
# cd /path/to/project/root

# Display Python version
echo "Python version: $(python3 --version 2>&1)"

# Initialize array to store failed tests
failed_tests=()

# Initialize step counter
step_counter=1

echo "-={ $python_file }=-"

# Function to handle each test step
run_test_step() {
  local step_name=$1
  local command=$2
  local failure_message=$3

  echo "-={ Step $step_counter: $step_name }=-"
  if ! eval "$command"; then
    failed_tests+=("$failure_message in $python_file")
  fi
  ((step_counter++))
}

# Static analysis steps: black, isort, mypy, pylint, pydocstyle, darglint, doctests

# Step 1: black
run_test_step "black" "black --check \"$python_file\"" "black"

# Step 2: isort
run_test_step "isort" "isort --check-only \"$python_file\"" "isort"

# Step 3: mypy
run_test_step "mypy" "mypy --strict --pretty --allow-untyped-calls \"$python_file\"" "mypy"

# Step 4: pylint
echo "-={ Step $step_counter: pylint }=-"
pylint_score=$(pylint --rcfile=.pylintrc "$python_file" | tail -2 | grep -Eo '[0-9\.]+' | head -1 || echo "0")
echo "Pylint score: $pylint_score"
if (( $(echo "$pylint_score < 9.91" | bc -l) )); then
  pylint --rcfile=.pylintrc "$python_file" || true
  echo "Pylint score below 9.91, failing..."
  failed_tests+=("pylint in $python_file")
fi
((step_counter++))

# Step 5: pydocstyle
run_test_step "pydocstyle" "pydocstyle \"$python_file\"" "pydocstyle"

# Step 6: darglint
run_test_step "darglint" "darglint -v 2 \"$python_file\"" "darglint"

# Step 7: doctests
run_test_step "doctests" "python3 \"$python_file\"" "doctests"

# Exit with failure if any tests failed
if [ ${#failed_tests[@]} -ne 0 ]; then
  echo "Failed tests: ${failed_tests[*]}"
  exit 1
fi

echo "All tests passed."
exit 0
